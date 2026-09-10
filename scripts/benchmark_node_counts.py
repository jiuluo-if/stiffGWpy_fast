# -*- coding: utf-8 -*-
"""Diagnostic 76/80/90/110-node convergence sweep on the current fast HEAD.

The seed values are chosen only to hit exact node counts; this is a measurement
baseline, not a production refinement rule.
"""
from __future__ import annotations

import json
import os
import sys
import time

os.environ.setdefault('NUMBA_THREADING_LAYER', 'workqueue')
import numpy as np
import psutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

TARGET_SEEDS = ((76, 64), (80, 68), (90, 79), (110, 101))
KW = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
REFERENCE_DN = 0.002262832966946746


def main():
    process = psutil.Process()
    available = process.cpu_affinity()
    if len(available) < 20:
        raise SystemExit('当前进程可用 CPU 少于 20')
    process.cpu_affinity(available[:20])
    FS.apply_accuracy_mode('fast')
    FS.set_threads(20)
    original = FA.goal_oriented_freqs
    rows = []
    try:
        for target, seed in TARGET_SEEDS:
            def builder(*args, _seed=seed, **kwargs):
                kwargs['seed_n'] = _seed
                return original(*args, **kwargs)

            FA.goal_oriented_freqs = builder
            model = LCDM_SG(**KW)
            start = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
            elapsed = time.perf_counter() - start
            dn = float(model.DN_gw[-1])
            rows.append({
                'target_n_freq': target,
                'seed_n': seed,
                'n_freq': int(model.f.size),
                'runtime_s': elapsed,
                'DN_gw': dn,
                'DN_rel_to_dense_reference': abs(dn - REFERENCE_DN) / REFERENCE_DN,
                'reuse': bool(getattr(model, 'outer_full_reuse_used', False)),
                'failure': getattr(model, 'fast_failure_reason', None),
                'omega_max': float(np.max(model.Ogw_today - model.Oj_today)),
                'omega_nonzero': int(np.count_nonzero(model.Ogw_today - model.Oj_today)),
                'f_min': float(np.min(model.f)),
                'f_max': float(np.max(model.f)),
            })
    finally:
        FA.goal_oriented_freqs = original
    dns = [row['DN_gw'] for row in rows]
    errors = [row['DN_rel_to_dense_reference'] for row in rows]
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'threads': 20,
        'affinity': process.cpu_affinity(),
        'threading_layer': __import__('numba').threading_layer(),
        'reference_DN_gw': REFERENCE_DN,
        'rows': rows,
        'DN_monotonic': all(dns[i] >= dns[i + 1] for i in range(len(dns) - 1)),
        'absolute_error_monotonic': all(errors[i] >= errors[i + 1] for i in range(len(errors) - 1)),
    }
    with open('docs/benchmark_node_counts_head.json', 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
