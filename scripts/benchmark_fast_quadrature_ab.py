"""A/B the fast frequency quadrature: default Simpson versus PCHIP.

The residual decomposition (`docs/fast_residual_decomposition_assessment.md`)
showed that fast's propagation kernel is accurate to ~1e-5 while the default
Simpson frequency quadrature contributes 4.2e-4..1.3e-2.  This script measures
what making PCHIP the default would buy in true DN accuracy and what it would
cost in warm runtime.  No solver default is changed.

Measurements alternate between the two quadratures inside each repeat (with
the order flipped on odd repeats) so machine drift cannot favour one method.
"""

import argparse
import json
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def _solve(point, method):
    """跑一次完整 SGWB_iter_fast，返回 (model, DN_gw)。"""
    model = LCDM_SG(**CASES[point])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      frequency_quadrature=method)
    dn = float(np.asarray(model.DN_gw, dtype=float).reshape(-1)[-1])
    return model, dn


def _load_anchor(anchor_dir, point):
    path = os.path.join(anchor_dir, f'oracle_c_wkb_{point}.json')
    if not os.path.exists(path):
        return None
    with open(path, encoding='utf-8') as handle:
        return json.load(handle)['summary']


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff'])
    parser.add_argument('--methods', nargs='+', default=['simpson', 'pchip'])
    parser.add_argument('--warmup', type=int, default=3)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--anchor-dir', default='docs')
    parser.add_argument('--out', default='docs/fast_quadrature_ab.json')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    rows = []
    for point in args.points:
        anchor = _load_anchor(args.anchor_dir, point)
        entry = {'point': point, 'methods': {}}
        for method in args.methods:
            model, dn = _solve(point, method)
            record = {
                'DN_gw': dn,
                'DN_eff': float(model.cosmo_param['DN_eff']),
            }
            if anchor is not None:
                record['vs_wkb_rel'] = _rel(dn, anchor['DN_wkb_corrected'])
                record['vs_deep_rel'] = _rel(dn, anchor['DN_deep'])
            record['times_ms'] = []
            entry['methods'][method] = record
        for method in args.methods:
            for _ in range(args.warmup):
                _solve(point, method)
        for repeat in range(args.repeats):
            order = args.methods if repeat % 2 == 0 else list(reversed(args.methods))
            for method in order:
                start = time.perf_counter()
                _solve(point, method)
                entry['methods'][method]['times_ms'].append(
                    (time.perf_counter() - start) * 1e3)
        for method in args.methods:
            times = np.asarray(entry['methods'][method]['times_ms'])
            entry['methods'][method]['warm_median_ms'] = float(np.median(times))
            entry['methods'][method]['warm_p95_ms'] = float(
                np.percentile(times, 95))
            entry['methods'][method]['warm_min_ms'] = float(times.min())
            entry['methods'][method]['n_repeats'] = int(times.size)
        if len(args.methods) >= 2:
            base = entry['methods'][args.methods[0]]['warm_median_ms']
            alt = entry['methods'][args.methods[1]]['warm_median_ms']
            entry['runtime_ratio'] = alt / base
            entry['runtime_ratio_label'] = '%s_over_%s' % (args.methods[1],
                                                           args.methods[0])
        rows.append(entry)

    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=2),
        'warmup': args.warmup,
        'repeats': args.repeats,
        'rows': rows,
        'semantics': ('fast SGWB_iter_fast wall time and self-consistent '
                      'DN_gw for each frequency quadrature, same native grid '
                      'and accuracy mode'),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    for row in rows:
        parts = []
        for method, record in row['methods'].items():
            parts.append('%s: med=%.2f ms p95=%.2f DN=%.6e vs_wkb=%.2e'
                         % (method, record['warm_median_ms'],
                            record['warm_p95_ms'], record['DN_gw'],
                            record.get('vs_wkb_rel', float('nan'))))
        print('%-8s %s' % (row['point'], ' | '.join(parts)))
        if 'runtime_ratio' in row:
            print('         runtime %s = %.3f' % (row['runtime_ratio_label'],
                                                  row['runtime_ratio']))
    print('wrote', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
