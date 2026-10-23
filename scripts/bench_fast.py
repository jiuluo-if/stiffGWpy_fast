# -*- coding: utf-8 -*-
"""Compare one LSODA run with fast-solver timings on selected cases.

The benchmark records one LSODA runtime, the first fast call (including JIT
compilation when needed), warm fast-call samples, and LSODA-to-fast runtime
ratios. Timed fast calls use the fast preset's grid and quadrature, with its
tolerance and frequency resolution unless overridden by CLI options. By
default they measure the unsplit grid; ``--kink-split`` adds the exact
reheating breakpoint used by the formal fast preset. Safety-wrapper diagnostics
use the formal fast preset's kink split and record that configuration
separately. LSODA is a runtime/regression baseline here, not an accuracy
reference. The reported 95th percentile uses NumPy's linear interpolation;
increase ``--reps`` beyond its default of 3 for a more informative tail sample.

中文说明：脚本在选定参数点上测量一次 LSODA 调用、fast 首次调用（包括可能发生的 JIT 编译）和多次预热
fast 调用，并报告两者的运行时间比。计时调用使用 fast 档位的网格与求积方法，容差和频率分辨率可由命令行覆盖。
默认测量不含转折点切分的网格；``--kink-split`` 会加入精确再加热转折点。安全包装器诊断使用正式 fast 档位的
转折点切分，并在结果中单独记录该配置。LSODA 在此只作运行时间/回归基线，不作为精度参照。95 分位数按 NumPy
默认线性插值计算；默认仅重复 3 次，增加 ``--reps`` 才能更可靠地观察尾部耗时。

Usage:
    python scripts/bench_fast.py [--reps 5] [--cases 0 3] [--json out.jsonl]
"""
import argparse
import json
import os
import statistics
import subprocess
import sys
import time

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy_fast import fast_sgwb  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    "A baseline r=1e-2 cr=1 T_re=2e3 k10=1e-2": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "B low r=1e-3":                           dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    "C r=3.6e-2":                             dict(r=3.6e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "D high r=0.1":                           dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    "E low T_re=10 GeV":                      dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    "F high T_re=1e4 GeV":                    dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "G low kappa10=1e-3":                     dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-3),
    "H high kappa10=1":                       dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    "I cr=0 kappa10=1 T_re=1e3":              dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    "J cr=0 baseline":                        dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    "K extreme r=0.1 T_re=10 k10=1":          dict(r=1e-1, cr=1, T_re=1e1, kappa10=1.0),
    "L cr=0 T_re=1e4 k10=0.1":                dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
}


def env_meta():
    meta = {
        'python': sys.version.split()[0],
        'platform': sys.platform,
        'p95_method': 'numpy.percentile default linear',
        'resources': telemetry(
            threads=int(os.environ.get('FAST_THREADS', '2'))),
        'threads': fast_sgwb._THREADS,
        'col_step': fast_sgwb._COL_STEP,
        'numpy': np.__version__,
    }
    for mod in ('numba', 'scipy'):
        try:
            meta[mod] = __import__(mod).__version__
        except Exception:
            pass
    try:
        out = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                             text=True, timeout=10)
        if out.returncode == 0:
            meta['commit'] = out.stdout.strip()
    except Exception:
        pass
    return meta


def p95(values):
    """Return the linearly interpolated 95th percentile.

    中文：按 NumPy 默认线性插值方法计算 95 分位数。
    """
    return float(np.percentile(values, 95))


def run_case(name, kw, reps, kink_split=False, tol=None, freq_res=None):
    # Benchmark the documented fast/plain-grid preset explicitly.  The module
    # settings are process-global, so reset them before every case and honor
    # FAST_THREADS after applying the preset's own default thread value.
    # 显式基准测试文档中的 fast/plain-grid 档位。模块设置是进程全局状态，
    # 因此每个案例开始时都要重置，并在应用预设默认线程数后读取 FAST_THREADS。
    preset = fast_sgwb.apply_accuracy_mode('fast')
    tol = preset['tol'] if tol is None else tol
    freq_res = preset['freq_res'] if freq_res is None else freq_res
    if os.environ.get('FAST_THREADS'):
        fast_sgwb.set_threads(int(os.environ['FAST_THREADS']))
    rec = {
        'case': name,
        'kw': {k: repr(v) for k, v in kw.items()},
        'tol': tol,
        'freq_res': freq_res,
        'freq_grid': preset['freq_grid'],
        'frequency_quadrature': 'pchip',
        'kink_split': bool(kink_split),
    }
    m = LCDM_SG(**kw)
    if m.derived_param['N_inf'] is None:
        rec['skip'] = 'N_inf None (invalid combo)'
        return rec
    t0 = time.perf_counter()
    m.SGWB_iter(engine='lsoda')
    rec['t_orig'] = time.perf_counter() - t0

    mf = LCDM_SG(**kw)
    t0 = time.perf_counter()
    fast_sgwb.SGWB_iter_fast(
        mf, tol=tol, freq_res=freq_res, kink_split=kink_split,
        freq_grid=preset['freq_grid'], frequency_quadrature='pchip')
    rec['t_fast_cold_s'] = time.perf_counter() - t0

    warm = []
    last = None
    for _ in range(reps):
        mf = LCDM_SG(**kw)
        t0 = time.perf_counter()
        fast_sgwb.SGWB_iter_fast(
            mf, tol=tol, freq_res=freq_res, kink_split=kink_split,
            freq_grid=preset['freq_grid'], frequency_quadrature='pchip')
        warm.append(time.perf_counter() - t0)
        last = mf
    rec['t_fast_warm_ms'] = [t * 1e3 for t in warm]
    rec['t_fast_min_ms'] = min(warm) * 1e3
    rec['t_fast_med_ms'] = statistics.median(warm) * 1e3
    rec['t_fast_p95_ms'] = p95(warm) * 1e3
    rec['speedup_min'] = rec['t_orig'] / min(warm)
    rec['speedup_med'] = rec['t_orig'] / statistics.median(warm)
    rec['speedup_p95'] = rec['t_orig'] / p95(warm)
    # Exercise the public safety wrapper once so benchmark records include
    # numerical-failure observability (direct ``SGWB_iter_fast`` calls above
    # intentionally measure the fast kernel without fallback overhead). This
    # diagnostic always uses the formal fast preset's exact kink split.
    # 另调用一次公开安全包装器，记录数值失败和回退情况；上面的直接调用用于测速，
    # 刻意不计入 fallback 开销。此诊断始终使用正式 fast 档位的精确转折点切分。
    mf_fb = LCDM_SG(**kw)
    mf_fb.SGWB_iter(engine='fast', fallback=True, tol=tol, freq_res=freq_res)
    rec['diagnostic_config'] = {
        'accuracy_mode': 'fast',
        'tol': tol,
        'freq_res': freq_res,
        'freq_grid': preset['freq_grid'],
        'frequency_quadrature': 'pchip',
        'kink_split': True,
        'fallback': True,
    }
    rec['fast_evals'] = int(getattr(mf_fb, 'fast_evals', 0))
    rec['fast_failures'] = int(getattr(mf_fb, 'fast_failures', 0))
    rec['fast_guard_rejections'] = int(getattr(mf_fb, 'fast_guard_rejections', 0))
    rec['lsoda_fallbacks'] = int(getattr(mf_fb, 'lsoda_fallbacks', 0))
    rec['fallback_fraction'] = (float(mf_fb.lsoda_fallbacks / mf_fb.fast_evals)
                                if mf_fb.fast_evals else 0.0)
    rel_dn = 0.0
    if m.SGWB_converge and last.SGWB_converge and m.DN_gw[-1] != 0:
        rel_dn = abs(m.DN_gw[-1] - last.DN_gw[-1]) / abs(m.DN_gw[-1])
    rec['DN_gw_last_rel'] = rel_dn
    rec['conv_orig'] = bool(m.SGWB_converge)
    rec['conv_fast'] = bool(last.SGWB_converge)
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='LSODA versus fast runtime benchmark / LSODA 与 fast 运行时间对比')
    ap.add_argument('--reps', type=int, default=3, help='warm fast repeats per case')
    ap.add_argument('--cases', nargs='+', type=int, default=None,
                    help='case indices (default: all 12)')
    ap.add_argument('--json', default=None, help='JSONL output path (optional)')
    ap.add_argument('--kink-split', action='store_true',
                    help='include the exact reheating breakpoint / 加入精确再加热转折点')
    ap.add_argument('--tol', type=float, default=None,
                    help='fast outer convergence tolerance (default: fast preset) / fast 外层收敛容差')
    ap.add_argument('--freq-res', type=float, default=None,
                    help='fast frequency resolution (default: fast preset) / fast 频率分辨率')
    args = ap.parse_args(argv)
    if args.reps < 1:
        ap.error('--reps must be at least 1')
    names = list(CASES)
    which = args.cases if args.cases is not None else list(range(len(names)))
    rows = []
    print('%-34s %11s %12s %12s %9s %9s %9s %11s %9s' %
          ('case', 'orig (s)', 'cold (s)', 'min (ms)', 'med (ms)',
           'spd min', 'spd med', 'p95 (ms)', 'fallback'))
    print('-' * 130)
    for i in which:
        rec = run_case(names[i], CASES[names[i]], args.reps, args.kink_split,
                       args.tol, args.freq_res)
        rows.append(rec)
        if rec.get('skip'):
            print('%-34s %s' % (names[i], 'invalid combo'))
            continue
        print('%-34s %11.3f %12.3f %12.3f %9.3f %8.0fx %8.0fx %11.3f %9.3g' %
              (names[i], rec['t_orig'], rec['t_fast_cold_s'], rec['t_fast_min_ms'],
               rec['t_fast_med_ms'], rec['speedup_min'], rec['speedup_med'],
               rec['t_fast_p95_ms'], rec['fallback_fraction']))
    if args.json:
        payload = {'meta': env_meta(), 'reps': args.reps, 'results': rows}
        with open(args.json, 'w', encoding='utf-8') as fh:
            fh.write(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
        print('wrote %s' % args.json)


if __name__ == '__main__':
    main()
