"""Full native-grid Prüfer amplitude-phase certification against Cartesian.

Extends the sparse-frequency screening in ``benchmark_prufer_oracle.py`` to
the complete formal ``freq_grid='goal'`` native frequency set, including an
outer self-consistency replay and (for the default point) a deterministic
replay of the full spectrum.  Reference-only: no formal kernel is changed.
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

from scripts.benchmark_prufer_oracle import (  # noqa: E402
    CASES,
    compare_outer,
    spectrum_prufer,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _statistics(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {'median': None, 'p95': None, 'max': None, 'n': 0}
    return {
        'median': float(np.median(finite)),
        'p95': float(np.percentile(finite, 95)),
        'max': float(np.max(finite)),
        'n': int(finite.size),
    }


def fullgrid_spectrum_compare(model, freqs, dn_eff, z_tail, rtol=1e-10):
    start = time.perf_counter()
    prufer = spectrum_prufer(model, freqs, dn_eff, z_tail, rtol=rtol)
    ogw, oj, opgw, used_tail = REF.spectrum_reference(
        model, freqs, dn_eff, z_tail=z_tail, rtol=rtol, workers=1)
    spectrum_s = time.perf_counter() - start
    omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
    g2, _, _ = REF.integrate_spectrum(freqs, ogw, oj)
    cartesian_dn = float(gp.Neff0 * g2 / omega_nu)
    prufer_summary = {
        'DN_gw': float(prufer['DN_gw']),
        'quadrature_error': float(prufer['quadrature_error']),
        'interpolation_error': (
            None if prufer['interpolation_error'] is None
            else float(prufer['interpolation_error'])),
        'used_tail_fraction': float(prufer['used_tail_fraction']),
        'rows': prufer['rows'],
    }
    return {
        'z_tail': float(z_tail),
        'prufer_DN_gw': float(prufer['DN_gw']),
        'cartesian_DN_gw': cartesian_dn,
        'DN_relative_error': abs(prufer['DN_gw'] - cartesian_dn)
        / max(abs(cartesian_dn), 1e-300),
        'Ogw_relative_error': _statistics(
            np.abs(prufer['Ogw'] - ogw) / np.maximum(np.abs(ogw), 1e-300)),
        'Oj_relative_error': _statistics(
            np.abs(prufer['Oj'] - oj) / np.maximum(np.abs(oj), 1e-300)),
        'Opgw_relative_error': _statistics(
            np.abs(prufer['Opgw'] - opgw) / np.maximum(np.abs(opgw), 1e-300)),
        'used_tail_fraction_prufer': float(prufer['used_tail_fraction']),
        'used_tail_fraction_cartesian': float(np.mean(used_tail)),
        'runtime_s': spectrum_s,
        'spectrum_prufer': prufer_summary,
    }


def replay_spectrum(model, freqs, dn_eff, z_tail, rtol=1e-10):
    first = spectrum_prufer(model, freqs, dn_eff, z_tail, rtol=rtol)
    second = spectrum_prufer(model, freqs, dn_eff, z_tail, rtol=rtol)
    ogw_a, oj_a, opgw_a, used_a = REF.spectrum_reference(
        model, freqs, dn_eff, z_tail=z_tail, rtol=rtol, workers=1)
    ogw_b, oj_b, opgw_b, used_b = REF.spectrum_reference(
        model, freqs, dn_eff, z_tail=z_tail, rtol=rtol, workers=1)
    return {
        'prufer_DN_bitwise': float(first['DN_gw']) == float(second['DN_gw']),
        'prufer_Ogw_bitwise': bool(np.array_equal(first['Ogw'], second['Ogw'])),
        'prufer_Oj_bitwise': bool(np.array_equal(first['Oj'], second['Oj'])),
        'prufer_Opgw_bitwise': bool(np.array_equal(first['Opgw'], second['Opgw'])),
        'cartesian_bitwise': bool(np.array_equal(ogw_a, ogw_b)
                                  and np.array_equal(oj_a, oj_b)
                                  and np.array_equal(opgw_a, opgw_b)
                                  and np.array_equal(used_a, used_b)),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default'])
    parser.add_argument('--z-tails', nargs='+', type=float, default=[5.0])
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--replay', action='store_true',
                        help='deterministic full-spectrum replay for default')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    for point in args.points:
        model = LCDM_SG(**CASES[point])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature='simpson')
        freqs = np.sort(np.asarray(model.f, dtype=float))
        dn_eff = float(model.cosmo_param['DN_eff'])
        spectra = [fullgrid_spectrum_compare(
            model, freqs, dn_eff, z_tail, rtol=args.rtol)
            for z_tail in args.z_tails]
        outer = [compare_outer(model, freqs, z_tail, rtol=args.rtol)
                 for z_tail in args.z_tails]
        payload = {
            'schema_version': 1,
            'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
            'point': point,
            'n_freq': int(freqs.size),
            'frequencies': freqs.tolist(),
            'z_tails': [float(value) for value in args.z_tails],
            'DN_eff': dn_eff,
            'resources': telemetry(workers=1, threads=2),
            'spectra': spectra,
            'outer': outer,
            'replay': replay_spectrum(model, freqs, dn_eff, args.z_tails[0],
                                      rtol=args.rtol) if args.replay else None,
            'semantics': 'full native-grid Prüfer versus Cartesian DOP853',
        }
        out = f'docs/oracle_prufer_fullgrid_{point}.json'
        with open(out, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
