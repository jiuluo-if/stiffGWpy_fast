# -*- coding: utf-8 -*-
"""测试不同 goal-grid seed 的 DN 收敛，不改变正式默认 seed。"""
from __future__ import annotations

import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
import psutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def main():
    process = psutil.Process()
    process.cpu_affinity(process.cpu_affinity()[:20])
    os.environ['FAST_THREADS'] = '20'
    FS.apply_accuracy_mode('fast')
    FS.set_threads(20)
    original = FA.goal_oriented_freqs
    kw = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    rows = []
    # These seeds target the requested 76/80/90/110 neighborhood; the actual
    # unique node count is recorded because reserves and native nodes vary.
    for seed in (64, 68, 78, 98):
        def builder(*args, _seed=seed, **kwargs):
            kwargs['seed_n'] = _seed
            return original(*args, **kwargs)
        FA.goal_oriented_freqs = builder
        for quadrature in ('simpson', 'pchip'):
            model = LCDM_SG(**kw)
            start = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                              frequency_quadrature=quadrature)
            runtime_s = time.perf_counter() - start
            rows.append({
                'seed_n': seed,
                'quadrature': quadrature,
                'n_freq': int(len(model.f)),
                'runtime_s': runtime_s,
                'DN_gw': float(model.DN_gw[-1]),
                'freq_grid_error': float(getattr(model, 'freq_grid_error', 0.0)),
                'reuse': bool(getattr(model, 'outer_full_reuse_used', False)),
                'failure': getattr(model, 'fast_failure_reason', None),
            })
    FA.goal_oriented_freqs = original
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'threads': 20,
        'affinity': process.cpu_affinity(),
        'threading_layer': __import__('numba').threading_layer(),
        'rows': rows,
    }
    with open('docs/benchmark_goal_seeds_head.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
