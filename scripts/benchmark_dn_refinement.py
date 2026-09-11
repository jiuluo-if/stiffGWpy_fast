# -*- coding: utf-8 -*-
"""原型：依据局部 DN 积分误差挑选 midpoint，不改正式 goal builder。"""
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
import psutil  # noqa: E402
from scipy import interpolate  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def local_priority(freqs, omega, count):
    order = np.argsort(freqs)
    x = np.asarray(freqs)[order]
    y = np.asarray(omega)[order]
    spl = interpolate.PchipInterpolator(x, y)
    score = np.empty(x.size - 1)
    for i in range(score.size):
        pchip = float(spl.integrate(x[i], x[i + 1]))
        trap = 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1])
        score[i] = abs(pchip - trap)
    chosen = np.argsort(score)[::-1][:count]
    return 0.5 * (x[chosen] + x[chosen + 1])


def solve(eval_freqs=None):
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    start = time.perf_counter()
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      eval_freqs=eval_freqs)
    runtime_s = time.perf_counter() - start
    omega = model.Ogw_today - model.Oj_today
    x = np.sort(model.f)
    y = omega[np.argsort(model.f)]
    pchip_dn = gp.Neff0 * FS.ln10 * float(
        interpolate.PchipInterpolator(x, y).integrate(x[0], x[-1])) / (
            gp.Omega_nh2 / model.derived_param['h']**2)
    return model, runtime_s, pchip_dn


def main(threads=2):
    process = psutil.Process()
    limit_affinity(process, threads)
    os.environ['FAST_THREADS'] = str(threads)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    base, runtime_s, pchip_dn = solve()
    additions = local_priority(base.f, base.Ogw_today - base.Oj_today, 20)
    rows = [{
        'add_count': 0,
        'n_freq': int(len(base.f)),
        'runtime_s': runtime_s,
        'DN_gw_simpson': float(base.DN_gw[-1]),
        'DN_gw_pchip': pchip_dn,
        'selected_nodes': [],
    }]
    for count in (4, 10, 20):
        model, elapsed, pchip = solve(additions[:count])
        rows.append({
            'add_count': count,
            'n_freq': int(len(model.f)),
            'runtime_s': elapsed,
            'DN_gw_simpson': float(model.DN_gw[-1]),
            'DN_gw_pchip': pchip,
            'selected_nodes': additions[:count].tolist(),
        })
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=threads),
        'threading_layer': __import__('numba').threading_layer(),
        'rows': rows,
    }
    with open('docs/benchmark_dn_refinement_head.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=2)
    main(parser.parse_args().threads)
