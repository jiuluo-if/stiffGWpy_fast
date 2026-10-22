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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'lowT', 'highT', 'stiff'])
    parser.add_argument('--anchor-dir', default='docs')
    parser.add_argument('--out', default='docs/fast_true_error.json')
    parser.add_argument('--quadrature', default='pchip',
                        help='fast frequency integral under audit')
    parser.add_argument('--h', type=float, default=None,
                        help='override the fast background-grid step')
    parser.add_argument('--z-tail', type=float, default=None,
                        help='override the numerical-to-analytic handoff depth')
    parser.add_argument('--tol', type=float, default=None,
                        help='override the fast preset outer convergence tolerance')
    parser.add_argument('--freq-res', type=float, default=None,
                        help='override the fast preset frequency-grid resolution')
    parser.add_argument('--transition-refine', action='store_true',
                        help='use the validation transition-refinement path')
    parser.add_argument('--transition-refine-radius', type=float, default=None,
                        help='refinement radius in e-folds around reheating')
    parser.add_argument('--transition-refine-factor', type=int, default=None,
                        help='number of substeps per background step near reheating')
    parser.add_argument('--warmup', type=int, default=0,
                        help='successful solver calls to discard per point')
    parser.add_argument('--repeats', type=int, default=1,
                        help='measured solver calls per point')
    args = parser.parse_args(argv)
    if args.warmup < 0 or args.repeats < 1:
        parser.error('--warmup must be >= 0 and --repeats must be >= 1')

    preset = FS.apply_accuracy_mode('fast')
    config = FS.resolve_config(
        'fast', h=args.h, z_tail=args.z_tail,
        kink_split=not args.transition_refine,
        transition_refine_radius=args.transition_refine_radius,
        transition_refine_factor=args.transition_refine_factor)
    tol = preset['tol'] if args.tol is None else args.tol
    freq_res = preset['freq_res'] if args.freq_res is None else args.freq_res
    rows = []
    for point in args.points:
        anchor_path = os.path.join(
            args.anchor_dir, f'oracle_c_wkb_{point}.json')
        with open(anchor_path, encoding='utf-8') as handle:
            anchor = json.load(handle)
        summary = anchor['summary']
        anchor_config = anchor.get('configuration', {})
        expected_config = {
            'path': 'transition-refine' if args.transition_refine else 'kink-split',
            'h': config.h,
            'transition_refine_radius': config.transition_refine_radius,
            'transition_refine_factor': config.transition_refine_factor,
        }
        legacy_defaults = {
            'transition_refine_radius': 2.0,
            'transition_refine_factor': 8,
        }
        for key, expected in expected_config.items():
            recorded = anchor_config.get(key, legacy_defaults.get(key))
            if recorded != expected:
                raise RuntimeError(
                    'fast/reference %s mismatch for %s: %r != %r; '
                    'regenerate the Oracle C anchor with the same solver config'
                    % (key, point, expected, recorded))
        results = []
        runtimes = []
        for call_index in range(args.warmup + args.repeats):
            model = LCDM_SG(**CASES[point])
            started = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=config.kink_split,
                              transition_refine=args.transition_refine,
                              config=config, tol=tol, freq_res=freq_res,
                              freq_grid='goal',
                              frequency_quadrature=args.quadrature)
            elapsed = time.perf_counter() - started
            if not model.SGWB_converge or not hasattr(model, 'DN_gw'):
                reason = getattr(model, 'fast_failure_reason', None) or 'not converged'
                raise RuntimeError('fast solver failed for %s: %s' % (point, reason))
            results.append((
                float(np.asarray(model.DN_gw, dtype=float).reshape(-1)[-1]),
                float(model.cosmo_param['DN_eff']),
            ))
            if call_index >= args.warmup:
                runtimes.append(elapsed)
        fast_dn, fast_dn_eff = results[-1]
        if any(result != results[-1] for result in results[:-1]):
            raise RuntimeError('fast solver result changed across repeats for %s' % point)
        anchor_dn_eff = float(anchor['DN_eff'])
        if fast_dn_eff != anchor_dn_eff:
            raise RuntimeError(
                'fast/reference DN_eff mismatch for %s: %.17g != %.17g; '
                'regenerate the Oracle C anchor with the same solver config'
                % (point, fast_dn_eff, anchor_dn_eff))
        anchor_freqs = np.asarray(anchor['frequencies'], dtype=float)
        fast_freqs = np.sort(np.asarray(model.f, dtype=float))
        if not np.array_equal(fast_freqs, anchor_freqs):
            raise RuntimeError(
                'fast/reference frequency grid mismatch for %s; '
                'regenerate the Oracle C anchor with the same solver config'
                % point)
        rows.append({
            'point': point,
            'fast_DN_gw': fast_dn,
            'fast_DN_eff': fast_dn_eff,
            'reference_DN_eff': anchor_dn_eff,
            'DN_frozen_z5': summary['DN_frozen'],
            'DN_wkb_z5': summary['DN_wkb_corrected'],
            'DN_deep_z10': summary['DN_deep'],
            'anchor_n_freq': summary['n_freq'],
            'fast_n_freq': int(fast_freqs.size),
            'anchor_DN_eff_match': True,
            'frequency_grid_exact_match': True,
            'fast_runtime_s_median': float(np.median(runtimes)),
            'fast_runtime_s_samples': runtimes,
            'fast_vs_frozen_z5_rel': _rel(fast_dn, summary['DN_frozen']),
            'fast_vs_wkb_z5_rel': _rel(fast_dn, summary['DN_wkb_corrected']),
            'fast_vs_deep_z10_rel': _rel(fast_dn, summary['DN_deep']),
            'frozen_to_wkb_rel': summary['frozen_minus_deep_rel'],
        })
    payload = {
        'schema_version': 1,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(workers=1, threads=2),
        'configuration': {
            'profile': 'fast',
            'path': ('transition-refine' if args.transition_refine else 'kink-split'),
            'h': config.h,
            'col_step': config.col_step,
            'z_tail': config.z_tail,
            'transition_refine_radius': config.transition_refine_radius,
            'transition_refine_factor': config.transition_refine_factor,
            'phase_max': config.phase_max,
            'threads': config.threads,
            'freq_grid': 'goal',
            'freq_res': freq_res,
            'tol': tol,
            'frequency_quadrature': args.quadrature,
            'warmup': args.warmup,
            'repeats': args.repeats,
        },
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
