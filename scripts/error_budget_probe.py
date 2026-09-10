# -*- coding: utf-8 -*-
"""对 fast HEAD 做低成本的 DN_gw 组件隔离探针。"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=20)
    ap.add_argument('--affinity-count', type=int, default=20)
    ap.add_argument('--case', choices=('default', 'lowT', 'highT', 'stiff'), default='default')
    ap.add_argument('--out', default='docs/error_budget_probe_head.json')
    args = ap.parse_args()
    os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
    os.environ['FAST_THREADS'] = str(args.threads)
    import psutil
    process = psutil.Process()
    process.cpu_affinity(process.cpu_affinity()[:args.affinity_count])
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    from stiffgwpy_fast import exact_background as EB
    from stiffgwpy_fast import fast_sgwb as FS
    from stiffgwpy_fast import global_param as gp
    from stiffgwpy_fast import reference as REF
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    kw = {
        'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
        'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
        'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
        'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    }[args.case]
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    original_primitive = EB.fast_phi_s2_split
    original_reuse = FS._OUTER_FULL_REUSE_ENABLED
    baseline_spectrum = None

    def solve(label, primitive='fast', reuse=True, frequency_quadrature='simpson'):
        nonlocal baseline_spectrum
        EB.fast_phi_s2_split = (original_primitive if primitive == 'fast'
                                else EB.exact_phi_s2_split)
        FS._OUTER_FULL_REUSE_ENABLED = reuse
        model = LCDM_SG(**kw)
        t0 = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                          frequency_quadrature=frequency_quadrature)
        runtime_s = time.perf_counter() - t0
        g2_simpson = float(model.DN_gw[-1])
        g2_pchip, qerr, ierr = REF.integrate_spectrum(
            model.f, model.Ogw_today, model.Oj_today)
        omega_nu = gp.Omega_nh2 / model.derived_param['h']**2
        g2_pchip_dn = float(gp.Neff0 * g2_pchip / omega_nu)
        spectrum = np.asarray(model.Ogw_today - model.Oj_today, dtype=float)
        if label == 'baseline':
            baseline_spectrum = spectrum.copy()
        spectrum_rel = None
        if baseline_spectrum is not None:
            spectrum_rel = float(np.max(
                np.abs(spectrum - baseline_spectrum)
                / np.maximum(np.abs(baseline_spectrum), 1e-300)))
        return {
            'label': label,
            'primitive': primitive,
            'reuse': reuse,
            'frequency_quadrature': frequency_quadrature,
            'runtime_s': runtime_s,
            'n_freq': int(len(model.f)),
            'DN_gw_solver': g2_simpson,
            'DN_gw_pchip_same_spectrum': g2_pchip_dn,
            'same_spectrum_quadrature_delta_rel': abs(g2_pchip_dn - g2_simpson) / abs(g2_simpson),
            'DN_rel_vs_baseline': None,
            'spectrum_max_rel_vs_baseline': spectrum_rel,
            'reference_integral_abs_error': qerr,
            'reference_interpolation_error': ierr,
            'outer_full_reuse_used': bool(getattr(model, 'outer_full_reuse_used', False)),
            'fast_failure_reason': getattr(model, 'fast_failure_reason', None),
        }

    try:
        rows = [
            solve('baseline', 'fast', True),
            solve('pchip_quadrature', 'fast', True, 'pchip'),
            solve('exact_background_primitive', 'exact', True),
            solve('outer_reuse_disabled', 'fast', False),
        ]
        baseline_dn = rows[0]['DN_gw_solver']
        for row in rows:
            row['DN_rel_vs_baseline'] = abs(
                row['DN_gw_solver'] - baseline_dn) / abs(baseline_dn)
    finally:
        EB.fast_phi_s2_split = original_primitive
        FS._OUTER_FULL_REUSE_ENABLED = original_reuse
    payload = {'commit': os.popen('git rev-parse HEAD').read().strip(),
               'case': args.case,
               'threads': args.threads, 'affinity': process.cpu_affinity(),
               'rows': rows}
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
