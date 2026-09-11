"""Quantify the true fast DN_gw error against the Oracle C WKB anchor.

``fast`` and the frozen-amplitude reference share the same z_tail frozen tail
convention, so a fast-vs-reference comparison cancels the shared tail defect.
Oracle C supplies an analytic WKB correction that removes that defect, which
lets this script measure the *true* fast error: the difference between the
self-consistent fast ``DN_gw`` and the WKB-corrected reference ``DN_gw`` at
the same ``DN_eff``.

The reference anchors are read from the committed ``docs/oracle_c_wkb_*.json``
artifacts (same native grid, same ``DN_eff``); fast is re-solved at the
current HEAD.  Reference-only diagnostic: no formal kernel is changed.
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff'])
    parser.add_argument('--anchor-dir', default='docs')
    parser.add_argument('--out', default='docs/fast_true_error.json')
    parser.add_argument('--quadrature', default='pchip',
                        help='fast frequency integral under audit')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    rows = []
    for point in args.points:
        anchor_path = os.path.join(
            args.anchor_dir, f'oracle_c_wkb_{point}.json')
        with open(anchor_path, encoding='utf-8') as handle:
            anchor = json.load(handle)
        summary = anchor['summary']
        model = LCDM_SG(**CASES[point])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature=args.quadrature)
        fast_dn = float(np.asarray(model.DN_gw, dtype=float).reshape(-1)[-1])
        fast_dn_eff = float(model.cosmo_param['DN_eff'])
        rows.append({
            'point': point,
            'fast_DN_gw': fast_dn,
            'fast_DN_eff': fast_dn_eff,
            'reference_DN_eff': float(anchor['DN_eff']),
            'DN_frozen_z5': summary['DN_frozen'],
            'DN_wkb_z5': summary['DN_wkb_corrected'],
            'DN_deep_z10': summary['DN_deep'],
            'anchor_n_freq': summary['n_freq'],
            'fast_vs_frozen_z5_rel': _rel(fast_dn, summary['DN_frozen']),
            'fast_vs_wkb_z5_rel': _rel(fast_dn, summary['DN_wkb_corrected']),
            'fast_vs_deep_z10_rel': _rel(fast_dn, summary['DN_deep']),
            'frozen_to_wkb_rel': summary['frozen_minus_deep_rel'],
        })
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=2),
        'rows': rows,
        'semantics': ('fast self-consistent DN_gw versus the frozen and '
                      'Oracle C WKB-corrected reference DN_gw on the same '
                      'native grid and DN_eff'),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    for row in rows:
        print('%-8s fast=%.10e | vs frozen=%.3e  vs wkb=%.3e  vs deep=%.3e'
              % (row['point'], row['fast_DN_gw'],
                 row['fast_vs_frozen_z5_rel'], row['fast_vs_wkb_z5_rel'],
                 row['fast_vs_deep_z10_rel']))
    print('wrote', args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
