# -*- coding: utf-8 -*-
"""比较同一 native fast spectrum 的频率积分表示，不改变 solver 默认行为。"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
import numpy as np
from scipy import integrate, interpolate

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


def _ascending(freqs, values):
    order = np.argsort(np.asarray(freqs))
    return np.asarray(freqs, dtype=float)[order], np.asarray(values, dtype=float)[order]


def _pchip(x, y):
    return interpolate.PchipInterpolator(x, y)


def _gauss_integral(x, y):
    spl = _pchip(x, y)
    nodes, weights = np.polynomial.legendre.leggauss(8)
    total = 0.0
    for left, right in zip(x[:-1], x[1:]):
        mid = 0.5 * (left + right)
        half = 0.5 * (right - left)
        total += half * np.dot(weights, spl(mid + half * nodes))
    return float(total)


def _loglog_pchip_integral(x, y):
    positive = np.maximum(y, np.finfo(float).tiny)
    spl = _pchip(x, np.log10(positive))
    nodes, weights = np.polynomial.legendre.leggauss(8)
    total = 0.0
    for left, right in zip(x[:-1], x[1:]):
        mid = 0.5 * (left + right)
        half = 0.5 * (right - left)
        total += half * np.dot(weights, np.power(10.0, spl(mid + half * nodes)))
    return float(total)


def _chebyshev_integral(x, y):
    lo, hi = float(x[0]), float(x[-1])
    scaled = 2.0 * (x - lo) / (hi - lo) - 1.0
    degree = min(20, x.size - 1)
    coeff = np.polynomial.chebyshev.chebfit(scaled, y, degree)
    anti = np.polynomial.chebyshev.chebint(coeff)
    return float((hi - lo) * 0.5 *
                 (np.polynomial.chebyshev.chebval(1.0, anti) -
                  np.polynomial.chebyshev.chebval(-1.0, anti)))


def integrals(freqs, omega):
    x, y = _ascending(freqs, omega)
    methods = {
        'simpson': float(integrate.simpson(y, x=x)),
        'pchip': float(_pchip(x, y).integrate(x[0], x[-1])),
        'natural_cubic': float(interpolate.CubicSpline(x, y).integrate(x[0], x[-1])),
        'loglog_pchip': _loglog_pchip_integral(x, y),
        'gauss_legendre_pchip': _gauss_integral(x, y),
        'chebyshev': _chebyshev_integral(x, y),
    }
    return methods


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=20)
    ap.add_argument('--affinity-count', type=int, default=20)
    ap.add_argument('--out', default='docs/frequency_quadrature_methods_head.json')
    args = ap.parse_args()
    os.environ['FAST_THREADS'] = str(args.threads)
    import psutil
    process = psutil.Process()
    available = process.cpu_affinity()
    if len(available) < args.affinity_count:
        raise SystemExit('当前进程可用 CPU 少于 affinity-count')
    process.cpu_affinity(available[:args.affinity_count])
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
            'reuse': bool(getattr(model, 'outer_full_reuse_used', False)),
            'failure': getattr(model, 'fast_failure_reason', None),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'threads': args.threads,
        'affinity': process.cpu_affinity(),
        'threading_layer': __import__('numba').threading_layer(),
        'rows': rows,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
