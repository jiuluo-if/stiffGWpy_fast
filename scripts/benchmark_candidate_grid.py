# -*- coding: utf-8 -*-
"""A/B test a candidate goal-grid density without changing production defaults."""
from __future__ import annotations

import argparse
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

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'low_r': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=78)
    ap.add_argument('--threads', type=int, default=20)
    ap.add_argument('--out', default='docs/benchmark_candidate_grid_head.json')
    args = ap.parse_args()

    process = psutil.Process()
    process.cpu_affinity(process.cpu_affinity()[:args.threads])
    os.environ['FAST_THREADS'] = str(args.threads)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    original = FA.goal_oriented_freqs
    rows = []

    def builder(*call_args, **kwargs):
        kwargs['seed_n'] = args.seed
        return original(*call_args, **kwargs)

    FA.goal_oriented_freqs = builder
    try:
        for name, kw in CASES.items():
            for quadrature in ('simpson', 'pchip'):
                model = LCDM_SG(**kw)
                start = time.perf_counter()
                result = FS.SGWB_iter_fast(
                    model, kink_split=True, freq_grid='goal',
                    frequency_quadrature=quadrature)
                rows.append({
                    'case': name,
                    'seed_n': args.seed,
                    'quadrature': quadrature,
                    'runtime_s': time.perf_counter() - start,
                    'n_freq': int(len(model.f)) if result is not None else None,
                    'DN_gw': (float(model.DN_gw[-1])
                              if result is not None else None),
                    'estimated_DN_quadrature_error_rel': (
                        float(getattr(model, 'estimated_DN_quadrature_error_rel', 0.0))
                        if result is not None else None),
                    'reuse': (bool(getattr(model, 'outer_full_reuse_used', False))
                              if result is not None else None),
                    'failure': getattr(model, 'fast_failure_reason', None),
                })
    finally:
        FA.goal_oriented_freqs = original

    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'seed_n': args.seed,
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
