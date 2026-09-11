# -*- coding: utf-8 -*-
"""Replay local DN quadrature estimator coverage on existing reference data."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from scipy import integrate, interpolate  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                     T_re=1e3, kappa10=1e-3),
    'positive_tilt': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
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
}

ESTIMATOR_METHODS = ('pchip', 'gauss2', 'gauss3', 'gauss5')


def _reference_path(name):
    return ROOT / 'docs' / ('frequency_same_grid_reference_%s.json' % name)


def fixed_spectrum_dense_reference(freqs, integrand, subdivisions=64):
    """Integrate a fixed native spectrum on each interval by dense Simpson.

    The dense curve is the PCHIP interpolant of the *same* native samples.
    Consequently this is a quadrature-only reference: no ODE solve or new
    spectrum samples enter the comparison.  ``subdivisions`` is forced even
    so SciPy's composite Simpson rule applies on every interval.
    """
    x = np.asarray(freqs, dtype=np.float64)
    y = np.asarray(integrand, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size or x.size < 2:
        raise ValueError('freqs and integrand must be equal 1-D arrays with at least 2 points')
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('freqs and integrand must be finite')
    order = np.argsort(x)
    xa = x[order]
    ya = y[order]
    if np.any(np.diff(xa) <= 0.0):
        raise ValueError('freqs must contain unique nodes')
    n_sub = int(subdivisions)
    if n_sub < 2:
        raise ValueError('subdivisions must be at least 2')
    if n_sub % 2:
        n_sub += 1
    spline = interpolate.PchipInterpolator(xa, ya)
    local = np.empty(xa.size - 1, dtype=np.float64)
    for i, (left, right) in enumerate(zip(xa[:-1], xa[1:])):
        dense_x = np.linspace(left, right, n_sub + 1)
        local[i] = integrate.simpson(spline(dense_x), x=dense_x)
    return local


def _panel_error_distribution(baseline, reference):
    """Distribute two-interval Simpson panel errors without cancellation."""
    baseline = np.asarray(baseline, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if baseline.shape != reference.shape:
        raise ValueError('baseline and reference must have the same shape')
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


def summarize_estimator_coverage(rows, estimator):
    """Summarize conservative coverage for one estimator column."""
    valid = [row for row in rows if row['actual_rel'] is not None]
    covered = [row for row in valid if row['predicted_rel'] >= row['actual_rel']]
    safe = [row for row in valid if row['predicted_rel'] < row['actual_rel']]
    ratios = [row['predicted_rel'] / max(row['actual_rel'], 1e-300)
              for row in valid]
    inverse = [row['actual_rel'] / max(row['predicted_rel'], 1e-300)
               for row in valid]
    return {
        'estimator': estimator,
        'n': len(valid),
        'covered': len(covered),
        'coverage': len(covered) / len(valid) if valid else None,
        'false_safe': len(safe),
        'max_actual_rel': max((row['actual_rel'] for row in valid),
                              default=None),
        'median_prediction_to_actual': float(np.median(ratios)) if ratios else None,
        'max_prediction_to_actual': max(ratios, default=None),
        'actual_over_prediction_p95': float(np.percentile(inverse, 95)) if inverse else None,
        'actual_over_prediction_p99': float(np.percentile(inverse, 99)) if inverse else None,
    }


def main():
    FS.apply_accuracy_mode('fast')
    rows = []
    for name, kw in CASES.items():
        reference_file = _reference_path(name)
        if not reference_file.exists():
            continue
        with reference_file.open(encoding='utf-8') as handle:
            reference = json.load(handle)
        model = LCDM_SG(**kw)
        result = FS.SGWB_iter_fast(
            model, kink_split=True, freq_grid='goal',
            frequency_quadrature='simpson')
        if result is None:
            rows.append({'case': name, 'failure': model.fast_failure_reason})
            continue
        dense_local = fixed_spectrum_dense_reference(
            model.f, model.Ogw_today - model.Oj_today)
        scale = max(abs(FS.integrate_frequency_quadrature(
            model.f, model.Ogw_today - model.Oj_today, 'simpson')), 1e-300)
        pure_actual = None
        method_rows = {}
        for method in ESTIMATOR_METHODS:
            errors, _, baseline = FS.estimate_frequency_quadrature_local(
                model.f, model.Ogw_today - model.Oj_today, method,
                allocation='panel_envelope')
            if pure_actual is None:
                pure_actual = _panel_error_distribution(baseline, dense_local)
            method_rows[method] = {
                'predicted_sum_rel': float(np.sum(errors) / scale),
                'predicted_max_rel': float(np.max(errors) / scale),
            }
        pchip = method_rows['pchip']
        simpson_dn = float(reference['fast_simpson_dn'])
        actual = reference.get('dn_relative_error', {}).get('pchip')
        rows.append({
            'case': name,
            'n_freq': int(model.f.size),
            'actual_rel': float(actual) if actual is not None else None,
            'predicted_sum_rel': pchip['predicted_sum_rel'],
            'predicted_max_rel': pchip['predicted_max_rel'],
            'pure_quadrature_actual_sum_rel': float(np.sum(pure_actual) / scale),
            'pure_quadrature_predicted_sum_rel': pchip['predicted_sum_rel'],
            'pure_quadrature_actual_max_rel': float(np.max(pure_actual) / scale),
            'pure_quadrature_predicted_max_rel': pchip['predicted_max_rel'],
            'pure_quadrature_covered': bool(
                pchip['predicted_max_rel'] >= float(np.max(pure_actual) / scale)),
            'pure_quadrature_by_method': method_rows,
            'reference_commit': reference.get('commit'),
            'simpson_dn_reference': simpson_dn,
            'failure': getattr(model, 'fast_failure_reason', None),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(threads=2),
        'estimator_method': 'pchip_minus_local_simpson_panel',
        'rows': rows,
        'summary': {
            key: summarize_estimator_coverage(
                [{**row, 'predicted_rel': row[key], 'actual_rel': row.get('actual_rel')}
                 for row in rows if row.get('actual_rel') is not None], key)
            for key in ('predicted_sum_rel', 'predicted_max_rel')
        },
        'pure_quadrature_summary': {
            'sum': summarize_estimator_coverage([
                {'actual_rel': row['pure_quadrature_actual_sum_rel'],
                 'predicted_rel': row['pure_quadrature_predicted_sum_rel']}
                for row in rows], 'pchip_minus_simpson_sum'),
            'max': summarize_estimator_coverage([
                {'actual_rel': row['pure_quadrature_actual_max_rel'],
                 'predicted_rel': row['pure_quadrature_predicted_max_rel']}
                for row in rows], 'pchip_minus_simpson_max'),
        },
        'pure_quadrature_method_summaries': {
            method: {
                'sum': summarize_estimator_coverage([
                    {'actual_rel': row['pure_quadrature_actual_sum_rel'],
                     'predicted_rel': row['pure_quadrature_by_method'][method][
                         'predicted_sum_rel']}
                    for row in rows], method + '_minus_simpson_sum'),
                'max': summarize_estimator_coverage([
                    {'actual_rel': row['pure_quadrature_actual_max_rel'],
                     'predicted_rel': row['pure_quadrature_by_method'][method][
                         'predicted_max_rel']}
                    for row in rows], method + '_minus_simpson_max'),
            }
            for method in ESTIMATOR_METHODS
        },
    }
    output = ROOT / 'docs' / 'quadrature_estimator_coverage_head.json'
    with output.open('w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
