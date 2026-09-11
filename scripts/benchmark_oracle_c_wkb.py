"""Oracle C: higher-order WKB tail correction versus the frozen-amplitude tail.

The reference tail uses a frozen-amplitude analytic continuation
(``Th = coeff * exp(-z + N - N_inf)``).  That continuation already carries the
leading-order WKB amplitude evolution ``d ln A / dN = 1.5*sigma - 2`` (the
physical amplitude ``h = A exp(-z + N)`` is frozen at leading order).  The
residual defect comes from the fast oscillating term in the exact Prüfer
amplitude equation, ``d ln h / dN = cos(2*theta)``, whose stationary-phase
boundary term integrates to ``sin(2*theta_f) / omega_f``.

This script tests the resulting fully analytic correction

    transfer^2_wkb = transfer^2_frozen * (1 + sin(2*theta_f) / omega_f)

against a deep numerical handoff (``z=10``, where the local adiabaticity
``eps ~ 1e-4``).  Three today estimates are compared on the same native grid:

1. ``frozen(z5)``: current pipeline, handoff at z=5 then frozen amplitude;
2. ``wkb(z5)``: same handoff with the analytic correction above;
3. ``deep(z10)``: handoff at z=10, frozen thereafter (numerical reference).

Reference-only diagnostic: no formal kernel is changed.
"""

import argparse
import json
import math
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

from scripts.benchmark_prufer_oracle import CASES, solve_prufer  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def wkb_transfer_correction(theta_handoff, z_handoff):
    """Analytic O(1/omega) boundary correction to the frozen transfer-squared.

    ``d ln h / dN = cos(2*theta)`` with ``dtheta/dN = -omega`` integrates to
    ``sin(2*theta_f) / omega_f`` at leading order, so the transfer-squared
    (which carries ``h^2``) is corrected by ``1 + sin(2*theta_f)/omega_f``.
    """
    omega = math.exp(z_handoff)
    return 1.0 + math.sin(2.0 * theta_handoff) / omega


def _today_from_transfer(bg, freq, transfer_squared, tensor_power):
    x_squared = transfer_squared * 10.0 ** (
        2.0 * (freq - REF.f_hor_abs(bg, freq, bg.N_inf)))
    oj = -transfer_squared * tensor_power / 3.0
    opgw = x_squared * tensor_power / 36.0
    return {
        'Ogw_today': float(3.0 * opgw + oj),
        'Oj_today': float(oj),
        'Opgw_today': float(opgw),
    }


def wkb_corrected_today(model, bg, freq, row):
    """Analytic WKB-corrected today values from a handoff row (used_tail only)."""
    tensor_power = model.derived_param['A_t'] \
        * (10.0 ** freq / gp.f_piv) ** model.derived_param['nt']
    n_inf = bg.N_inf
    # solve_prufer 存的是 amplitude_handoff = A/sqrt(2)，其平方正好是
    # frozen tail 使用的相平均功率 A^2/2。
    frozen_ts = row['amplitude_handoff'] ** 2 \
        * math.exp(-2.0 * row['z_handoff'] + 2.0 * row['event_N']
                   - 2.0 * n_inf)
    correction = wkb_transfer_correction(row['phase_handoff'],
                                         row['z_handoff'])
    return _today_from_transfer(bg, freq, frozen_ts * correction, tensor_power)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default'])
    parser.add_argument('--freq-count', type=int, default=8)
    parser.add_argument('--rtol', type=float, default=1e-10)
    parser.add_argument('--out-dir', default='docs')
    args = parser.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    for point in args.points:
        model = LCDM_SG(**CASES[point])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature='simpson')
        freqs = np.sort(np.asarray(model.f, dtype=float))
        picks = np.linspace(0, freqs.size - 1, args.freq_count, dtype=int)
        freqs = freqs[picks]
        dn_eff = float(model.cosmo_param['DN_eff'])
        bg = REF._Background(model, dn_eff)
        omega_nu = gp.Omega_nh2 / model.derived_param['h'] ** 2
        rows = []
        start = time.perf_counter()
        for freq in freqs:
            shallow = solve_prufer(model, float(freq), dn_eff, 5.0,
                                   rtol=args.rtol)
            deep = solve_prufer(model, float(freq), dn_eff, 10.0,
                                rtol=args.rtol)
            row = {
                'frequency': float(freq),
                'used_tail_shallow': bool(shallow['used_tail']),
                'used_tail_deep': bool(deep['used_tail']),
                'frozen_today_shallow': {
                    k: shallow[k] for k in
                    ('Ogw_today', 'Oj_today', 'Opgw_today')},
                'deep_today': {
                    k: deep[k] for k in
                    ('Ogw_today', 'Oj_today', 'Opgw_today')},
            }
            if shallow['used_tail']:
                wkb = wkb_corrected_today(model, bg, float(freq), shallow,
                                          )
                row['wkb_today'] = wkb
                row['omega_handoff'] = math.exp(shallow['z_handoff'])
                row['theta_handoff'] = float(shallow['phase_handoff'])
                row['wkb_correction'] = wkb_transfer_correction(
                    shallow['phase_handoff'], shallow['z_handoff'])
                row['eps_handoff'] = abs(
                    1.5 * REF._H2_and_sigma(bg, shallow['event_N'])[1] - 1.0
                ) / math.exp(shallow['z_handoff'])
            else:
                row['wkb_today'] = None
            rows.append(row)
        elapsed = time.perf_counter() - start

        def dn_of(values, key):
            ogw = np.asarray([v[key]['Ogw_today'] for v in values])
            oj = np.asarray([v[key]['Oj_today'] for v in values])
            g2, _, _ = REF.integrate_spectrum(freqs, ogw, oj)
            return float(gp.Neff0 * g2 / omega_nu)

        used = [r for r in rows if r['wkb_today'] is not None]
        dn_frozen = dn_of(rows, 'frozen_today_shallow')
        dn_deep = dn_of(rows, 'deep_today')
        ogw_w = np.asarray([
            r['wkb_today']['Ogw_today'] if r['wkb_today'] is not None
            else r['frozen_today_shallow']['Ogw_today'] for r in rows])
        oj_w = np.asarray([
            r['wkb_today']['Oj_today'] if r['wkb_today'] is not None
            else r['frozen_today_shallow']['Oj_today'] for r in rows])
        g2w, _, _ = REF.integrate_spectrum(freqs, ogw_w, oj_w)
        dn_wkb = float(gp.Neff0 * g2w / omega_nu)
        compared = [r for r in used if r['used_tail_deep']]

        def residual_stats(key):
            vals = []
            for r in compared:
                a = r[key]['Ogw_today']
                b = r['deep_today']['Ogw_today']
                vals.append(abs(a - b) / max(abs(b), 1e-300))
            vals = np.asarray(vals)
            return {
                'n': int(vals.size),
                'median': float(np.median(vals)) if vals.size else None,
                'max': float(np.max(vals)) if vals.size else None,
            }

        eps_vals = np.asarray([r['eps_handoff'] for r in used])
        payload = {
            'schema_version': 1,
            'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
            'point': point,
            'frequencies': freqs.tolist(),
            'DN_eff': dn_eff,
            'resources': telemetry(workers=1, threads=2),
            'rows': rows,
            'summary': {
                'n_freq': int(freqs.size),
                'used_tail_count': len(used),
                'compared_count': len(compared),
                'DN_frozen': dn_frozen,
                'DN_wkb_corrected': dn_wkb,
                'DN_deep': dn_deep,
                'frozen_minus_deep_rel': abs(dn_frozen - dn_deep)
                / max(abs(dn_deep), 1e-300),
                'wkb_minus_deep_rel': abs(dn_wkb - dn_deep)
                / max(abs(dn_deep), 1e-300),
                'per_mode_frozen_vs_deep': residual_stats(
                    'frozen_today_shallow'),
                'per_mode_wkb_vs_deep': residual_stats('wkb_today'),
                'eps_handoff_median': float(np.median(eps_vals)),
                'eps_handoff_max': float(np.max(eps_vals)),
                'runtime_s': elapsed,
            },
            'semantics': ('Oracle C: frozen-amplitude (z=5) vs analytic '
                          'O(1/omega) WKB-corrected (z=5) vs deep (z=10) '
                          'numerical handoff'),
        }
        out = os.path.join(args.out_dir,
                           f'oracle_c_wkb_{point}.json')
        with open(out, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('wrote', out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
