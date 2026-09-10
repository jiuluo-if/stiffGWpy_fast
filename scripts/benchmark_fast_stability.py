# -*- coding: utf-8 -*-
"""Run a current-HEAD fast-only stability screen over named and Sobol points."""
from __future__ import annotations

import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
import numpy as np
import psutil
from scipy.stats import qmc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SINGLE = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0, T_re=1e3, kappa10=1e-3),
    'negative_tilt': dict(r=1e-2, cr=1, n_t=-0.4, T_re=2e3, kappa10=1e-2),
}


def sample_points(n=16):
    sampler = qmc.Sobol(d=6, scramble=True, seed=20260910)
    rows = []
    for i, u in enumerate(sampler.random(n)):
        rows.append((f'sobol_{i:03d}', {
            'r': 1e-6 * (1e5 ** float(u[0])),
            'n_t': -0.5 + float(u[1]),
            'cr': 1.0 if u[2] >= 0.5 else 0.0,
            'T_re': 1e1 * (1e5 ** float(u[3])),
            'DN_re': 30.0 * float(u[4]),
            'kappa10': 1e-6 * (1e6 ** float(u[5])),
        }))
    return rows


def main():
    process = psutil.Process()
    available = process.cpu_affinity()
    if len(available) < 20:
        raise SystemExit('当前进程可用 CPU 少于 20')
    process.cpu_affinity(available[:20])
    FS.apply_accuracy_mode('fast')
    FS.set_threads(20)
    config = FS.FastSolverConfig(h=0.005, col_step=8, z_tail=5.0,
                                 phase_max=0.25, freq_grid='goal', threads=20)
    points = list(SINGLE.items()) + sample_points()
    rows = []
    for label, kw in points:
        model = LCDM_SG(**kw)
        start = time.perf_counter()
        try:
            result = FS.SGWB_iter_fast(model, config=config, kink_split=True)
            error = None
        except Exception as exc:  # pragma: no cover - diagnostic guard
            result = None
            error = f'{type(exc).__name__}: {exc}'
        elapsed = time.perf_counter() - start
        reason = getattr(model, 'fast_failure_reason', None)
        status = 'ok' if result is model and model.SGWB_converge else (
            'guard' if reason == 'shared_Neff_guard' else 'failed')
        rows.append({
            'label': label,
            'status': status,
            'error': error,
            'failure_reason': reason,
            'runtime_s': elapsed,
            'n_freq': int(getattr(model, 'f', np.empty(0)).size),
            'DN_gw': float(model.DN_gw[-1]) if hasattr(model, 'DN_gw') else None,
            'estimated_DN_quadrature_error': float(
                getattr(model, 'estimated_DN_quadrature_error', np.nan)),
            'estimated_DN_quadrature_error_rel': float(
                getattr(model, 'estimated_DN_quadrature_error_rel', np.nan)),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'threads': 20,
        'affinity': process.cpu_affinity(),
        'threading_layer': __import__('numba').threading_layer(),
        'n_points': len(rows),
        'rows': rows,
        'guard_count': sum(row['status'] == 'guard' for row in rows),
        'failure_count': sum(row['status'] == 'failed' for row in rows),
    }
    with open('docs/benchmark_fast_stability_head.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
