# -*- coding: utf-8 -*-
"""Compare fast frequency quadratures against the independent reference.

The reference is evaluated on the exact native fast grid.  This isolates the
frequency-integration error from the separate question of how the grid itself
was selected.
"""

import argparse
import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')

import numpy as np
from scipy.integrate import simpson

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                     T_re=1e3, kappa10=1e-3),
    'positive_tilt': dict(r=1e-2, cr=1, n_t=0.2,
                          T_re=2e3, kappa10=1e-2),
    'negative_tilt': dict(r=1e-2, cr=1, n_t=-0.4,
                          T_re=2e3, kappa10=1e-2),
    'sobol_000': dict(r=3.4845505440304265e-06, n_t=-0.12252988945692778,
                      cr=1.0, T_re=34989.571710260374,
                      DN_re=18.787082182243466, kappa10=5.28549585835966e-06),
    'sobol_002': dict(r=0.008513729433124471, n_t=-0.3924837075173855,
                      cr=1.0, T_re=118660.64687585819,
                      DN_re=11.174977160990238,
                      kappa10=0.0017982905763510582),
    'sobol_006': dict(r=0.002986405620277968, n_t=-0.16131426580250263,
                      cr=0.0, T_re=2832.147926998191,
                      DN_re=3.932188209146261, kappa10=0.01617670894806547),
    'sobol_010': dict(r=0.08346659409929896, n_t=-0.024900889955461025,
                      cr=0.0, T_re=4977.754805490329,
                      DN_re=26.626055845990777,
                      kappa10=0.15852849886128825),
    'sobol_015': dict(r=1.5975020128853773e-06, n_t=0.18561450485140085,
                      cr=0.0, T_re=114872.20639509047,
                      DN_re=5.934050939977169,
                      kappa10=2.79610302157714e-05),
}


def _rel(a, b):
    scale = np.maximum(np.abs(b), 1e-300)
    return np.abs(a - b) / scale


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--point', default='default', choices=sorted(CASES))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--rtol', type=float, default=1e-9)
    ap.add_argument('--z-tail', type=float, default=8.0)
    ap.add_argument('--seed-n', type=int, default=64,
                    help='goal-grid seed override for candidate A/B runs')
    ap.add_argument('--out', default='docs/frequency_same_grid_reference.json')
    args = ap.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    FS.set_z_tail(5.0)
    kw = CASES[args.point]

    from stiffgwpy_fast import freq_adaptive as FA
    original_goal = FA.goal_oriented_freqs

    def candidate_goal(*call_args, **kwargs):
        kwargs['seed_n'] = args.seed_n
        return original_goal(*call_args, **kwargs)

    FA.goal_oriented_freqs = candidate_goal

    try:
        fast = LCDM_SG(**kw)
        t0 = time.perf_counter()
        FS.SGWB_iter_fast(fast, kink_split=True, freq_grid='goal', frequency_quadrature='simpson')
        fast_time = time.perf_counter() - t0
        freqs = np.asarray(fast.f, dtype=float)
        dn_eff = float(fast.cosmo_param['DN_eff'])
        fast_omega = np.asarray(fast.Ogw_today, dtype=float)
        fast_j = np.asarray(fast.Oj_today, dtype=float)
        fast_integrand = fast_omega - fast_j

        ref_model = LCDM_SG(**kw)
        t0 = time.perf_counter()
        ref_omega, ref_j, _, used_tail = REF.spectrum_reference(
            ref_model, freqs, dn_eff, z_tail=args.z_tail, rtol=args.rtol,
            workers=args.workers)
        ref_time = time.perf_counter() - t0
        ref_integrand = ref_omega - ref_j
    finally:
        FA.goal_oriented_freqs = original_goal

    order = np.argsort(freqs)
    omega_nu = gp.Omega_nh2 / fast.derived_param['h']**2
    fast_simpson_g2 = float(simpson(fast_integrand[order], x=freqs[order]) * FS.ln10)
    fast_pchip_g2 = float(FS.integrate_frequency_pchip(freqs, fast_integrand) * FS.ln10)
    ref_g2, ref_quad_err, ref_interp_err = REF.integrate_spectrum(
        freqs, ref_omega, ref_j)
    fast_simpson_dn = gp.Neff0 * fast_simpson_g2 / omega_nu
    fast_pchip_dn = gp.Neff0 * fast_pchip_g2 / omega_nu
    ref_pchip_dn = gp.Neff0 * ref_g2 / omega_nu

    records = {
        'point': args.point,
        'kw': kw,
        'n_freq': int(freqs.size),
        'seed_n': args.seed_n,
        'dn_eff': dn_eff,
        'used_tail_fraction': float(np.mean(used_tail)),
        'fast_runtime_s': fast_time,
        'reference_runtime_s': ref_time,
        'fast_simpson_dn': fast_simpson_dn,
        'fast_pchip_dn': fast_pchip_dn,
        'reference_pchip_dn_same_grid': float(ref_pchip_dn),
        'reference_quadrature_error': ref_quad_err,
        'reference_interpolation_error': ref_interp_err,
        'dn_relative_error': {
            'simpson': float(abs(fast_simpson_dn - ref_pchip_dn) / abs(ref_pchip_dn)),
            'pchip': float(abs(fast_pchip_dn - ref_pchip_dn) / abs(ref_pchip_dn)),
        },
        'spectrum_relative_error': {
            'simpson_max': float(_rel(fast_integrand, ref_integrand).max()),
            'simpson_p95': float(np.percentile(_rel(fast_integrand, ref_integrand), 95)),
        },
    }
    print(json.dumps(records, ensure_ascii=False, indent=2), flush=True)
    with open(args.out, 'w', encoding='utf-8') as fh:
        json.dump(records, fh, ensure_ascii=False, indent=2)
        fh.write('\n')
    print('wrote %s' % args.out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
