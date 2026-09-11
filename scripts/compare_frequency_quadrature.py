# -*- coding: utf-8 -*-
"""比较同一 native fast spectrum 的频率积分表示，不改变 solver 默认行为。"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

try:
    from scripts._resource_budget import apply_environment, limit_affinity, telemetry
except ImportError:
    from _resource_budget import apply_environment, limit_affinity, telemetry

apply_environment()
import numpy as np  # noqa: E402
from scipy import interpolate  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def integrals(freqs, omega):
    methods = ('simpson', 'pchip', 'log_pchip', 'gauss2', 'gauss3',
               'gauss5', 'natural_cubic', 'chebyshev')
    return {
        method: FS.integrate_frequency_quadrature(freqs, omega, method)
        for method in methods
    }


def interpolation_diagnostics(freqs, omega):
    x, y = np.asarray(freqs)[np.argsort(freqs)], np.asarray(omega)[np.argsort(freqs)]
    dense = np.linspace(x[0], x[-1], max(256, 4 * x.size))
    diagnostics = {}
    for method in ('pchip', 'log_pchip', 'gauss2', 'gauss3', 'gauss5',
                   'natural_cubic', 'chebyshev'):
        if method == 'natural_cubic':
            spline = interpolate.CubicSpline(x, y)
            values = spline(dense)
        elif method == 'log_pchip':
            if np.any(y <= 0.0):
                diagnostics[method] = {'monotone': None, 'overshoot_rel': None}
                continue
            spline = interpolate.PchipInterpolator(x, np.log(y))
            values = np.exp(spline(dense))
        else:
            spline = interpolate.PchipInterpolator(x, y)
            values = spline(dense)
        interval_points = np.linspace(0.0, 1.0, 9)
        shape_preserving = True
        for left, right, y0, y1 in zip(x[:-1], x[1:], y[:-1], y[1:]):
            segment = spline(left + (right - left) * interval_points)
            tolerance = 1e-12 * max(1.0, abs(y0), abs(y1))
            if y1 >= y0:
                shape_preserving &= bool(np.all(np.diff(segment) >= -tolerance))
            else:
                shape_preserving &= bool(np.all(np.diff(segment) <= tolerance))
        scale = max(float(np.max(np.abs(y))), np.finfo(float).tiny)
        overshoot = max(0.0, float(np.max(values) - np.max(y)),
                        float(np.min(y) - np.min(values))) / scale
        diagnostics[method] = {
            'shape_preserving': bool(shape_preserving),
            'overshoot_rel': overshoot,
        }
    diagnostics['simpson'] = {'shape_preserving': None, 'overshoot_rel': 0.0}
    return diagnostics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--affinity-count', type=int, default=2)
    ap.add_argument('--reps', type=int, default=50)
    ap.add_argument('--out', default='docs/frequency_quadrature_methods_head.json')
    args = ap.parse_args()
    os.environ['FAST_THREADS'] = str(args.threads)
    import psutil
    process = psutil.Process()
    available = process.cpu_affinity()
    if len(available) < args.affinity_count:
        raise SystemExit('当前进程可用 CPU 少于 affinity-count')
    limit_affinity(process, args.affinity_count)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    rows = []
    for name, kw in CASES.items():
        model = LCDM_SG(**kw)
        start = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
        runtime_s = time.perf_counter() - start
        omega = model.Ogw_today - model.Oj_today
        values = integrals(model.f, omega)
        runtime_by_method = {}
        for method in values:
            samples = []
            for _ in range(args.reps):
                start_method = time.perf_counter()
                FS.integrate_frequency_quadrature(model.f, omega, method)
                samples.append(time.perf_counter() - start_method)
            runtime_by_method[method] = {
                'median_s': float(np.median(samples)),
                'p95_s': float(np.percentile(samples, 95)),
            }
        omega_nu = gp.Omega_nh2 / model.derived_param['h']**2
        dn = {key: float(gp.Neff0 * FS.ln10 * value / omega_nu)
              for key, value in values.items()}
        simpson = dn['simpson']
        rows.append({
            'case': name,
            'runtime_s': runtime_s,
            'n_freq': int(len(model.f)),
            'DN_gw_solver': float(model.DN_gw[-1]),
            'DN_by_method': dn,
            'relative_to_simpson': {key: abs(value - simpson) / abs(simpson)
                                    for key, value in dn.items()},
            'runtime_by_method': runtime_by_method,
            'interpolation_diagnostics': interpolation_diagnostics(model.f, omega),
            'reuse': bool(getattr(model, 'outer_full_reuse_used', False)),
            'failure': getattr(model, 'fast_failure_reason', None),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=args.threads),
        'threading_layer': __import__('numba').threading_layer(),
        'rows': rows,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
