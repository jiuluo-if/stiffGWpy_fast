"""Measure independent-reference sensitivity to the frozen-tail hand-off.

The fast spectrum and its ``DN_eff`` are held fixed.  Only the continuous-
sigma DOP853 reference ``z_tail`` changes, so the result is an oracle-floor
diagnostic rather than a fast-solver benchmark.
"""

import argparse
import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.benchmark_same_grid_reference import CASES
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

Z_TAILS = (5.0, 6.0, 7.0, 8.0, 10.0)


def run_point(point, rtol=1e-10, workers=4):
    """Run one fixed-background tail sweep and return an auditable record."""
    kw = CASES[point]
    FS.apply_accuracy_mode('fast')
    FS.set_z_tail(5.0)
    fast = LCDM_SG(**kw)
    t0 = time.perf_counter()
    FS.SGWB_iter_fast(fast, kink_split=True, freq_grid='goal',
                      frequency_quadrature='simpson')
    fast_runtime = time.perf_counter() - t0
    freqs = np.asarray(fast.f, dtype=float)
    dn_eff = float(fast.cosmo_param['DN_eff'])

    rows = []
    for z_tail in Z_TAILS:
        ref_model = LCDM_SG(**kw)
        t0 = time.perf_counter()
        result = REF.run_reference(
            ref_model, dn_eff=dn_eff, freq_res=1.0, z_tail=z_tail,
            rtol=rtol, freq_subset=freqs, self_consistent=False,
            workers=workers)
        runtime = time.perf_counter() - t0
        rows.append({
            'z_tail': z_tail,
            'DN_gw': float(result['DN_gw']),
            'DN_eff': float(result['DN_eff']),
            'runtime_s': runtime,
            'used_tail_fraction': float(np.mean(result['used_tail'])),
            'quadrature_error': float(result['quadrature_error']),
            'interpolation_error': float(result['interpolation_error']),
            'n_freq': int(result['n_freq']),
        })

    summary = REF.summarize_tail_convergence(rows)
    summary.update({
        'point': point,
        'parameters': kw,
        'fast_DN_eff': dn_eff,
        'fast_DN_gw': float(fast.DN_gw[-1]),
        'fast_runtime_s': fast_runtime,
        'rtol': rtol,
        'n_freq': int(freqs.size),
        'frequency_grid': 'same native fast goal grid',
        'oracle_semantics': 'fixed DN_eff; only reference z_tail varies',
        'acceptance_note': (
            'The observed deepest-tail difference is the reported systematic '
            'bound; fitted decay is descriptive and not used to shrink it.'),
    })
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', default=['default'],
                        choices=sorted(CASES))
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--workers', type=int, default=4,
                        help='reserved for compatibility; reference uses its default')
    parser.add_argument('--out', default='docs/oracle_tail_convergence_head.json')
    args = parser.parse_args(argv)

    records = [run_point(point, rtol=args.rtol, workers=args.workers)
               for point in args.points]
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'z_tails': list(Z_TAILS),
        'points': records,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    print('wrote %s' % args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
