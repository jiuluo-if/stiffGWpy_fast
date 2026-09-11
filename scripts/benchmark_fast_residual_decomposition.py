"""Per-frequency decomposition of the fast DN_gw error versus Oracle C.

``fast`` and the frozen-amplitude reference share the same ``z_tail`` tail
convention, so a fast-vs-reference comparison cancels the shared tail defect.
Oracle C removes that defect analytically
(``1 + sin(2 theta_f)/omega_f``).  This script decomposes the remaining true
fast DN_gw error into candidate sources:

1. frequency quadrature: the same fast per-node values integrated with the fast
   composite-Simpson rule versus the reference PCHIP rule;
2. tail amplitude: the per-node weighted difference against the WKB-corrected
   anchor, split into analytic-tail and non-tail nodes;
3. the shape of that per-frequency contribution (which nodes dominate).

Reference-only diagnostic: no formal kernel is changed.
"""

import argparse
import json
import os
import sys

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from scipy import integrate, interpolate  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _simpson(freqs, integrand):
    """升序 native 网格上的复合 Simpson 积分。"""
    return float(integrate.simpson(np.asarray(integrand, dtype=np.float64),
                                   x=np.asarray(freqs, dtype=np.float64)))


def _pchip(freqs, integrand):
    """PCHIP 解析积分，对应 reference integrate_spectrum 的插值部分。"""
    x = np.asarray(freqs, dtype=np.float64)
    y = np.asarray(integrand, dtype=np.float64)
    return float(interpolate.PchipInterpolator(x, y).integrate(x[0], x[-1]))


def _simpson_weights(freqs):
    """复用 fast 的 build_Wmat 得到升序网格上的精确 Simpson 权重。"""
    x = np.asarray(freqs, dtype=np.float64)
    weights = np.zeros((x.size, x.size))
    FS.build_Wmat(x.size, x, np.diff(x), weights)
    return weights[-1]


def _row_integral(rows, key):
    """从 anchor 行提取 Ogw - Oj 的 today 积分核。"""
    ogw = np.asarray([row[key]['Ogw_today'] for row in rows], dtype=np.float64)
    oj = np.asarray([row[key]['Oj_today'] for row in rows], dtype=np.float64)
    return ogw - oj


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _pearson(a, b):
    """返回 Pearson 相关系数；样本不足或零方差时返回 None。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 3 or a.std() == 0.0 or b.std() == 0.0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _top_contributors(freqs, contrib, count=5):
    """按 |贡献| 排序返回头部频率、占比与覆盖 50%/90% 所需节点数。"""
    order = np.argsort(-np.abs(contrib))
    total = float(np.abs(contrib).sum())
    ranked = np.abs(contrib)[order]
    cumulative = np.cumsum(ranked)

    def _cover(fraction):
        if total <= 0.0:
            return None
        hit = np.searchsorted(cumulative, fraction * total) + 1
        return int(min(hit, max(ranked.size, 1)))

    top = []
    for index in order[:count]:
        top.append({
            'frequency': float(freqs[index]),
            'contribution': float(contrib[index]),
            'share': (float(abs(contrib[index]) / total) if total > 0.0 else 0.0),
        })
    return {
        'top': top,
        'n_for_50pct': _cover(0.5),
        'n_for_90pct': _cover(0.9),
        'n_freq': int(freqs.size),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff'])
    parser.add_argument('--anchor-dir', default='docs')
    parser.add_argument('--out',
                        default='docs/fast_residual_decomposition.json')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    results = []
    for point in args.points:
        anchor_path = os.path.join(args.anchor_dir,
                                   f'oracle_c_wkb_{point}.json')
        with open(anchor_path, encoding='utf-8') as handle:
            anchor = json.load(handle)
        anchor_rows = anchor['rows']
        model = LCDM_SG(**CASES[point])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature='simpson')
        fast_freqs = np.asarray(model.f, dtype=np.float64)
        fast_ig = (np.asarray(model.Ogw_today, dtype=np.float64)
                   - np.asarray(model.Oj_today, dtype=np.float64))
        anchor_freqs = np.asarray([row['frequency'] for row in anchor_rows],
                                  dtype=np.float64)
        # fast 与 anchor 必须落在同一 native 网格，否则逐节点差没有意义
        if (fast_freqs.size != anchor_freqs.size
                or not np.allclose(np.sort(fast_freqs), np.sort(anchor_freqs),
                                   rtol=0.0, atol=1e-9)):
            raise RuntimeError(f'fast and anchor grids differ for {point}')

        order_f = np.argsort(fast_freqs)
        order_a = np.argsort(anchor_freqs)
        freqs = anchor_freqs[order_a]
        ig_fast = fast_ig[order_f]
        ig_frozen = _row_integral(anchor_rows, 'frozen_today_shallow')[order_a]
        ig_deep = _row_integral(anchor_rows, 'deep_today')[order_a]
        # 未入尾的频率没有 WKB 修正列，按定义回落 frozen 锚点
        raw_wkb = np.asarray([
            (row['wkb_today']['Ogw_today'] - row['wkb_today']['Oj_today'])
            if row['wkb_today'] is not None
            else (row['frozen_today_shallow']['Ogw_today']
                  - row['frozen_today_shallow']['Oj_today'])
            for row in anchor_rows], dtype=np.float64)
        ig_wkb = raw_wkb[order_a]
        correction = np.asarray(
            [row['wkb_correction'] if row.get('wkb_correction') else 1.0
             for row in anchor_rows], dtype=np.float64)[order_a]
        used_tail = np.asarray(
            [bool(row['used_tail_shallow']) for row in anchor_rows])[order_a]

        omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
        scale = gp.Neff0 * gp.ln10 / omega_nu
        dn = {}
        for key, values in (('fast', ig_fast), ('frozen', ig_frozen),
                            ('wkb', ig_wkb), ('deep', ig_deep)):
            dn[key] = {
                'simpson': scale * _simpson(freqs, values),
                'pchip': scale * _pchip(freqs, values),
            }

        weights = _simpson_weights(freqs)
        contrib = weights * (ig_fast - ig_wkb) * gp.ln10 * scale
        tail_mask = used_tail
        abs_tail = float(np.abs(contrib[tail_mask]).sum())
        abs_nontail = float(np.abs(contrib[~tail_mask]).sum())
        entry = {
            'point': point,
            'n_freq': int(freqs.size),
            'used_tail_count': int(used_tail.sum()),
            'dn': dn,
            'fast_vs_wkb_simpson_rel': _rel(dn['fast']['simpson'],
                                            dn['wkb']['simpson']),
            'fast_vs_wkb_pchip_rel': _rel(dn['fast']['pchip'],
                                          dn['wkb']['pchip']),
            'fast_simpson_vs_pchip_rel': _rel(dn['fast']['simpson'],
                                              dn['fast']['pchip']),
            'fast_vs_deep_simpson_rel': _rel(dn['fast']['simpson'],
                                             dn['deep']['simpson']),
            'frozen_vs_wkb_simpson_rel': _rel(dn['frozen']['simpson'],
                                              dn['wkb']['simpson']),
            'contribution_total': float(contrib.sum()),
            'contribution_tail_signed': float(contrib[tail_mask].sum()),
            'contribution_nontail_signed': float(contrib[~tail_mask].sum()),
            'contribution_tail_abs_share': abs_tail / max(abs_tail + abs_nontail,
                                                          1e-300),
            'top_contributors': _top_contributors(freqs, contrib),
            'tail_correction_pearson': _pearson(contrib[tail_mask],
                                                correction[tail_mask] - 1.0),
        }
        results.append(entry)

    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=2),
        'rows': results,
        'semantics': ('fast DN_gw error decomposition: fast-vs-WKB, '
                      'Simpson-vs-PCHIP quadrature, tail/non-tail per-node '
                      'contribution on the shared native grid'),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    for row in results:
        print('%-8s n=%2d tail=%2d | fast-vs-wkb sim=%.3e pch=%.3e | '
              'simp-vs-pch=%.3e | tail|abs|share=%.3f'
              % (row['point'], row['n_freq'], row['used_tail_count'],
                 row['fast_vs_wkb_simpson_rel'], row['fast_vs_wkb_pchip_rel'],
                 row['fast_simpson_vs_pchip_rel'],
                 row['contribution_tail_abs_share']))
        top = row['top_contributors']
        print('         top1 f=%.4f share=%.3f | n(50%%)=%s n(90%%)=%s | '
              'tail corr=%s'
              % (top['top'][0]['frequency'], top['top'][0]['share'],
                 top['n_for_50pct'], top['n_for_90pct'],
                 row['tail_correction_pearson']))
    print('wrote', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
