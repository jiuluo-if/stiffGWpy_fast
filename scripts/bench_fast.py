# -*- coding: utf-8 -*-
"""
bench_fast.py -- reproducible performance comparison for stiffgwpy.

Runs the original LSODA-based ``SGWB_iter()`` and the accelerated
``fast_sgwb.SGWB_iter_fast()`` on the same parameter grid.  Per case it
reports the cold (JIT/import) fast time, the warm distribution (min / median /
p95 over ``--reps`` repeats) and the resulting speedups, plus environment and
git-commit metadata.  ``--json`` emits a machine-readable record.

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

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy import fast_sgwb  # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG  # noqa: E402

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
        'cpu_count': os.cpu_count(),
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
    return sorted(values)[int(0.95 * (len(values) - 1))]


def run_case(name, kw, reps):
    # Benchmark the documented fast/plain-grid preset explicitly.  The module
    # settings are process-global, so reset them before every case and honor
    # FAST_THREADS after applying the preset's own default thread value.
    fast_sgwb.apply_accuracy_mode('fast')
    if os.environ.get('FAST_THREADS'):
        fast_sgwb.set_threads(int(os.environ['FAST_THREADS']))
    rec = {'case': name, 'kw': {k: repr(v) for k, v in kw.items()}}
    m = LCDM_SG(**kw)
    if m.derived_param['N_inf'] is None:
        rec['skip'] = 'N_inf None (invalid combo)'
        return rec
    t0 = time.perf_counter()
    m.SGWB_iter()
    rec['t_orig'] = time.perf_counter() - t0

    mf = LCDM_SG(**kw)
    t0 = time.perf_counter()
    fast_sgwb.SGWB_iter_fast(mf)
    rec['t_fast_cold_s'] = time.perf_counter() - t0

    warm = []
    last = None
    for _ in range(reps):
        mf = LCDM_SG(**kw)
        t0 = time.perf_counter()
        fast_sgwb.SGWB_iter_fast(mf)
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
    # intentionally measure the fast kernel without fallback overhead).
    mf_fb = LCDM_SG(**kw)
    mf_fb.SGWB_iter(engine='fast', fallback=True)
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
    ap = argparse.ArgumentParser(description='reproducible stiffgwpy benchmark')
    ap.add_argument('--reps', type=int, default=3, help='warm fast repeats per case')
    ap.add_argument('--cases', nargs='+', type=int, default=None,
                    help='case indices (default: all 12)')
    ap.add_argument('--json', default=None, help='JSONL output path (optional)')
    args = ap.parse_args(argv)
    names = list(CASES)
    which = args.cases if args.cases is not None else list(range(len(names)))
    rows = []
    print('%-34s %11s %12s %12s %9s %9s %9s %11s %9s' %
          ('case', 'orig (s)', 'cold (s)', 'min (ms)', 'med (ms)',
           'spd min', 'spd med', 'p95 (ms)', 'fallback'))
    print('-' * 130)
    for i in which:
        rec = run_case(names[i], CASES[names[i]], args.reps)
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
