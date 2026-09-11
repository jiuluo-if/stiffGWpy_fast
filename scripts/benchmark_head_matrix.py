# -*- coding: utf-8 -*-
"""固定环境下冻结六个代表点的 warm runtime 与 Simpson/PCHIP 差异。"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time

try:
    from scripts._resource_budget import apply_environment, limit_affinity, telemetry
except ImportError:
    from _resource_budget import apply_environment, limit_affinity, telemetry

apply_environment()
import numpy as np  # noqa: E402
import psutil  # noqa: E402

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


def pchip_dn(model):
    from scipy import interpolate

    order = np.argsort(model.f)
    x = np.asarray(model.f)[order]
    y = np.asarray(model.Ogw_today - model.Oj_today)[order]
    integral = interpolate.PchipInterpolator(x, y).integrate(x[0], x[-1])
    omega_nu = gp.Omega_nh2 / model.derived_param['h']**2
    return float(gp.Neff0 * FS.ln10 * integral / omega_nu)


def main(threads=2):
    quadrature = os.environ.get('FAST_QUADRATURE', 'pchip')
    process = psutil.Process()
    limit_affinity(process, threads)
    os.environ['FAST_THREADS'] = str(threads)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    rows = []
    for name, kw in CASES.items():
        timings = []
        last = None
        for _ in range(25):
            model = LCDM_SG(**kw)
            start = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                              frequency_quadrature=quadrature)
            timings.append((time.perf_counter() - start) * 1000.0)
            last = model
        warm = timings[1:]
        ordered = sorted(warm)
        rows.append({
            'case': name,
            'cold_ms': timings[0],
            'warm_min_ms': min(warm),
            'warm_median_ms': statistics.median(warm),
            'warm_p95_ms': ordered[int(0.95 * (len(ordered) - 1))],
            'n_freq': int(len(last.f)),
            'quadrature': getattr(last, 'frequency_quadrature_used', None),
            'DN_gw_default': float(last.DN_gw[-1]),
            'DN_gw_scipy_pchip_same_spectrum': pchip_dn(last),
            'scipy_pchip_delta_rel': abs(pchip_dn(last) - last.DN_gw[-1]) / abs(last.DN_gw[-1]),
            'reuse': bool(getattr(last, 'outer_full_reuse_used', False)),
            'failure': getattr(last, 'fast_failure_reason', None),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=threads),
        'repeats': 25,
        'rows': rows,
    }
    with open('docs/benchmark_head_matrix.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=2)
    main(parser.parse_args().threads)
