"""真实嵌套频率求值与 fast 局部积分误差估计器的覆盖率校验。

第九原则要求：不能在同一个 interpolant 上换多个求积算法冒充独立估计。
本脚本按

    base native grid -> 选 top-N 局部误差区间 -> 在每个区间新增真实
    frequency solve（把区间中点并入积分 support grid）-> 重新求积

得到 ``E_nested = |DN_refined - DN_base|``，再与同网格独立 reference 的
``actual_error`` 比较，统计 coverage / false-safe / monotonicity。

注意 ``eval_freqs`` 只把节点加入 solve 网格、不改变积分 support grid，
因此本脚本用受限的 monkeypatch 直接替换 ``goal_oriented_freqs`` 的返回值，
让新增中点真正进入积分网格；替换在 finally 中恢复，库默认行为不变。

``E_nested < actual_error``（coverage > 1）说明局部估计器 false-safe
（低估真实频率积分残差），必须显式记录，不能用来宣称精度达标。
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

# 六个正式 benchmark 代表点，参数与 benchmark_head_matrix.py 保持一致。
CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def _dn(model):
    """与同网格 reference artifact 对齐的 DN_gw：``m.DN_gw`` 的末元素。"""
    return float(np.asarray(model.DN_gw).reshape(-1)[-1])


def _head_commit():
    """记录生成证据时的 HEAD，便于与 release manifest 的 commit 字段对账。"""
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                              text=True, timeout=10).stdout.strip()
    except Exception:
        return None


def _solve(case, support_grid=None, quadrature='pchip'):
    """独立求解一次；``support_grid`` 非空时替换 goal 网格（实验用）。"""
    from stiffgwpy_fast import freq_adaptive as FA
    original = FA.goal_oriented_freqs
    if support_grid is not None:
        frozen = np.asarray(support_grid, dtype=float)
        FA.goal_oriented_freqs = lambda *a, **k: frozen
    model = LCDM_SG(**case)
    start = time.perf_counter()
    try:
        result = FS.SGWB_iter_fast(model, kink_split=True,
                                   frequency_quadrature=quadrature)
        runtime = time.perf_counter() - start
    finally:
        FA.goal_oriented_freqs = original
    return model, runtime, result is not None


def _reference_dn(name, ref_dir):
    """读取同网格独立 reference 的 DN（缺失时返回 None）。"""
    path = os.path.join(ref_dir, 'frequency_same_grid_reference_%s.json' % name)
    if not os.path.isfile(path):
        return None
    with open(path, 'r', encoding='utf-8') as handle:
        payload = json.load(handle)
    value = payload.get('reference_pchip_dn_same_grid')
    return None if value is None else float(value)


def run_point(name, case, n_values, reference_dn, reps, quadrature):
    """对一个代表点跑 base 与每档 N 的真实嵌套加密重解。"""
    record = {'point': name, 'kw': case, 'n_values': list(n_values)}
    record['quadrature'] = quadrature
    base_model, base_runtime, base_ok = _solve(case, quadrature=quadrature)
    if not base_ok:
        record['status'] = 'base_failed'
        return record
    f_base = np.asarray(base_model.f, dtype=float)
    dn_base = _dn(base_model)
    local = np.asarray(base_model.quadrature_error_local_by_interval, dtype=float)
    record.update({
        'status': 'ok',
        'n_freq_base': int(f_base.size),
        'dn_base': dn_base,
        'dn_eff_minus_orig': float(base_model.cosmo_param['DN_eff'] - base_model.DN_eff_orig),
        'base_runtime_s': base_runtime,
        'predicted_dn_error': float(base_model.estimated_DN_quadrature_error),
        'predicted_dn_error_rel': float(base_model.estimated_DN_quadrature_error_rel),
        'n_outer_iterations': int(len(getattr(base_model, 'DN_gw', []))),
        'outer_full_reuse_used': bool(getattr(base_model, 'outer_full_reuse_used', False)),
    })
    if reference_dn is not None:
        record['reference_dn'] = reference_dn
        record['actual_error'] = abs(dn_base - reference_dn)
        record['actual_error_rel'] = abs(dn_base - reference_dn) / abs(reference_dn)
    # 局部误差数组按 support grid 的递减区间排列：区间 i 是 (f[i+1], f[i])。
    order = np.argsort(local)[::-1]
    refinements = []
    for n_sel in n_values:
        n_sel = int(min(n_sel, local.size))
        chosen = np.sort(order[:n_sel])
        mids = 0.5 * (f_base[chosen] + f_base[chosen + 1])
        refined = np.unique(np.concatenate((f_base, mids)))
        refined = np.sort(refined)[::-1]
        best = None
        for _ in range(max(1, int(reps))):
            model, runtime, ok = _solve(case, support_grid=refined,
                                        quadrature=quadrature)
            if not ok:
                continue
            dn_refined = _dn(model)
            e_nested = abs(dn_refined - dn_base)
            candidate = {
                'n_intervals': int(n_sel),
                'n_midpoints': int(mids.size),
                'n_freq_refined': int(refined.size),
                'dn_refined': dn_refined,
                'e_nested': e_nested,
                'e_nested_rel': (e_nested / abs(dn_base)) if dn_base else None,
                'refined_runtime_s': runtime,
            }
            if best is None or candidate['refined_runtime_s'] < best['refined_runtime_s']:
                best = candidate
        if best is None:
            best = {'n_intervals': int(n_sel), 'status': 'refined_failed'}
        if 'actual_error' in record and best.get('e_nested'):
            best['coverage'] = record['actual_error'] / best['e_nested']
            best['false_safe'] = bool(best['coverage'] > 1.0)
        refinements.append(best)
    record['refinements'] = refinements
    # 单调性：N 增大时 |DN_refined - DN_base| 应单调不减（积分网格单调加密）。
    seq = [r.get('e_nested') for r in refinements]
    record['monotone_e_nested'] = all(
        a is None or b is None or b >= a - 1e-15 for a, b in zip(seq, seq[1:]))
    record['e_nested_sequence'] = seq
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description='nested native-frequency quadrature 覆盖率校验')
    parser.add_argument('--cases', default='default')
    parser.add_argument('--n-values', default='4,8,12')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--reps', type=int, default=1,
                        help='每个 N 重复求解次数，取最快一次作为 runtime 代表值')
    parser.add_argument('--reference-dir', default='docs')
    parser.add_argument('--quadrature', default='pchip',
                        help='pchip 是正式默认；simpson 作为网格加密灵敏度的阳性对照')
    parser.add_argument('--json', default=None)
    args = parser.parse_args(argv)

    # 固定到唯一正式 fast 预设：goal 频率网格 + exact kink split。
    FS.apply_accuracy_mode('fast')
    if args.threads:
        FS.set_threads(args.threads)
    n_values = [int(x) for x in args.n_values.split(',') if x.strip()]
    names = sorted(CASES) if args.cases == 'all' else args.cases.split(',')

    out = {
        'commit': _head_commit(),
        'threads': args.threads,
        'n_values': n_values,
        'reps': args.reps,
        'quadrature': args.quadrature,
        'resources': {
            'logical_cpu_count': os.cpu_count(),
            'numba_threads': int(os.environ.get('NUMBA_NUM_THREADS', '0') or 0),
            'fast_threads': int(os.environ.get('FAST_THREADS', '0') or 0),
            'numpy': np.__version__,
        },
        'points': {},
    }
    for name in names:
        reference_dn = _reference_dn(name, args.reference_dir)
        record = run_point(name, CASES[name], n_values, reference_dn, args.reps,
                           args.quadrature)
        out['points'][name] = record
        if record.get('status') == 'ok':
            print('%-12s dn_base=%.10e actual_rel=%s predicted_rel=%.3e' % (
                name, record['dn_base'],
                ('%.3e' % record['actual_error_rel']) if 'actual_error_rel' in record else 'n/a',
                record['predicted_dn_error_rel']))
            for item in record['refinements']:
                if 'e_nested' not in item:
                    print('    N=%-3d refined_failed' % item['n_intervals'])
                    continue
                print('    N=%-3d nfreq=%d e_nested=%.3e rel=%.3e coverage=%s false_safe=%s' % (
                    item['n_intervals'], item['n_freq_refined'], item['e_nested'],
                    item['e_nested_rel'],
                    ('%.3f' % item['coverage']) if 'coverage' in item else 'n/a',
                    item.get('false_safe')))
        else:
            print('%-12s %s' % (name, record.get('status')))
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as handle:
            json.dump(out, handle, indent=2, ensure_ascii=False)
            handle.write('\n')
        print('wrote %s' % args.json)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
