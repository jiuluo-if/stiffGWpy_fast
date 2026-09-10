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
    'highT': dict(r=1e-2, cr=1, T_re=1e7, kappa10=1e-2),
    'stiff': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2, w=0.6),
    'low_r': dict(r=1e-4, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e1),
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
    ap.add_argument('--out', default='docs/frequency_same_grid_reference.json')
    args = ap.parse_args(argv)

    FS.apply_accuracy_mode('fast')
    FS.set_z_tail(5.0)
    kw = CASES[args.point]

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
