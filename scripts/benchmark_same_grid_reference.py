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

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from scipy.interpolate import PchipInterpolator  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                     T_re=1e3, kappa10=1e-3),
    'positive_tilt': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
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


def _panel_error_distribution(baseline, reference):
    """Distribute two-interval Simpson errors without cancellation."""
    baseline = np.asarray(baseline, dtype=float)
    reference = np.asarray(reference, dtype=float)
    actual = np.zeros_like(baseline)
    for start in range(0, baseline.size - 1, 2):
        stop = min(start + 2, baseline.size)
        panel_error = abs(np.sum(baseline[start:stop]) -
                          np.sum(reference[start:stop]))
        weights = np.abs(baseline[start:stop])
        weight_sum = float(np.sum(weights))
        if weight_sum > 0.0:
            actual[start:stop] = panel_error * weights / weight_sum
        else:
            actual[start:stop] = panel_error / (stop - start)
    return actual


def ensemble_local_error(error_vectors):
    """Return the pointwise maximum of compatible local error vectors."""
    vectors = [np.asarray(vector, dtype=float) for vector in error_vectors]
    if not vectors:
        raise ValueError('at least one local error vector is required')
    shape = vectors[0].shape
    if any(vector.shape != shape for vector in vectors):
        raise ValueError('local error vectors must have the same shape')
    return np.max(np.stack(vectors, axis=0), axis=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--point', default='default', choices=sorted(CASES))
    ap.add_argument('--workers', type=int, default=1)
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

    omega_nu = gp.Omega_nh2 / fast.derived_param['h']**2
    quadrature_methods = ('simpson', 'pchip', 'log_pchip', 'gauss2',
                          'gauss3', 'gauss5', 'natural_cubic', 'chebyshev')
    fast_g2_by_method = {
        method: float(FS.integrate_frequency_quadrature(
            freqs, fast_integrand, method) * FS.ln10)
        for method in quadrature_methods
    }
    ref_g2, ref_quad_err, ref_interp_err = REF.integrate_spectrum(
        freqs, ref_omega, ref_j)
    ref_pchip_dn = gp.Neff0 * ref_g2 / omega_nu
    fast_dn_by_method = {
        method: gp.Neff0 * value / omega_nu
        for method, value in fast_g2_by_method.items()
    }
    estimator_rel_by_method = {}
    estimator_max_rel_by_method = {}
    for method in quadrature_methods:
        errors, _, _ = FS.estimate_frequency_quadrature_local(
            freqs, fast_integrand, method)
        estimator_rel_by_method[method] = float(
            np.sum(errors) / max(abs(fast_g2_by_method['simpson'] / FS.ln10),
                                 1e-300))
        estimator_max_rel_by_method[method] = float(
            np.max(errors) / max(abs(fast_g2_by_method['simpson'] / FS.ln10),
                                 1e-300))
    order = np.argsort(freqs)
    native_x = freqs[order]
    fast_y = fast_integrand[order]
    ref_y = ref_integrand[order]
    local_errors = [FS.estimate_frequency_quadrature_local(
        native_x, fast_y, method, allocation='panel_envelope')[0]
                    for method in quadrature_methods]
    fast_local_error, fast_local_candidate, fast_local_baseline = (
        FS.estimate_frequency_quadrature_local(
            native_x, fast_y, 'pchip', allocation='panel_envelope'))
    fast_local_error = ensemble_local_error(local_errors)
    ref_spline = PchipInterpolator(native_x, ref_y)
    ref_local = np.asarray([
        ref_spline.integrate(left, right)
        for left, right in zip(native_x[:-1], native_x[1:])])
    # Q2 actual: error of the formal Simpson panel against the independent
    # reference, distributed over its two native intervals.  The previous
    # version compared PCHIP against the reference, which measured spectrum
    # interpolation/ODE differences rather than the PCHIP-Simpson estimator.
    local_actual_error = _panel_error_distribution(
        fast_local_baseline, ref_local)
    local_scale = max(abs(ref_g2), 1e-300)
    if np.std(fast_local_error) > 0.0 and np.std(local_actual_error) > 0.0:
        local_corr = float(np.corrcoef(
            fast_local_error, local_actual_error)[0, 1])
    else:
        local_corr = None
    local_estimator = {
        'method': 'ensemble_minus_local_simpson_panel',
        'actual_sum_rel': float(np.sum(local_actual_error) / local_scale),
        'predicted_sum_rel': float(np.sum(fast_local_error) / local_scale),
        'actual_max_rel': float(np.max(local_actual_error) / local_scale),
        'predicted_max_rel': float(np.max(fast_local_error) / local_scale),
        'pearson_r': local_corr,
        'coverage': float(np.mean(fast_local_error >= local_actual_error)),
        'false_safe_rate': float(np.mean(fast_local_error < local_actual_error)),
    }
    local_ratio = local_actual_error / np.maximum(fast_local_error, 1e-300)
    local_estimator.update({
        'actual_over_predicted_p95': float(np.percentile(local_ratio, 95)),
        'actual_over_predicted_p99': float(np.percentile(local_ratio, 99)),
        'actual_over_predicted_max': float(np.max(local_ratio)),
        'required_safety_factor_95': float(np.percentile(local_ratio, 95)),
        'required_safety_factor_99': float(np.percentile(local_ratio, 99)),
        'mixed_pchip_reference_actual_sum_rel': float(
            np.sum(np.abs(fast_local_candidate - ref_local)) / local_scale),
    })

    records = {
        'point': args.point,
        'resources': telemetry(workers=args.workers, threads=2),
        'kw': kw,
        'n_freq': int(freqs.size),
        'seed_n': args.seed_n,
        'dn_eff': dn_eff,
        'used_tail_fraction': float(np.mean(used_tail)),
        'fast_runtime_s': fast_time,
        'reference_runtime_s': ref_time,
        'fast_simpson_dn': fast_dn_by_method['simpson'],
        'fast_pchip_dn': fast_dn_by_method['pchip'],
        'fast_dn_by_method': fast_dn_by_method,
        'estimator_local_sum_rel_by_method': estimator_rel_by_method,
        'estimator_local_max_rel_by_method': estimator_max_rel_by_method,
        'local_estimator': local_estimator,
        'reference_pchip_dn_same_grid': float(ref_pchip_dn),
        'reference_quadrature_error': ref_quad_err,
        'reference_interpolation_error': ref_interp_err,
        'dn_relative_error': {
            method: float(abs(value - ref_pchip_dn) / abs(ref_pchip_dn))
            for method, value in fast_dn_by_method.items()
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
