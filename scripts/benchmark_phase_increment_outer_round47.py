"""Alternating full end-to-end baseline/candidate benchmark for Round 47."""
from __future__ import annotations

import json
import os
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, limit_affinity, telemetry  # noqa: E402

apply_environment()

import psutil  # noqa: E402

from scripts.benchmark_phase_increment_kernel_twin import (  # noqa: E402
    CASES,
    solve_kernel_phase_increment,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def timed(case, candidate):
    original = FS.solve_kernel
    if candidate:
        FS.solve_kernel = lambda *args, **kw: solve_kernel_phase_increment(
            *args, 32)
    try:
        model = LCDM_SG(**CASES[case])
        start = time.perf_counter()
        result = FS.SGWB_iter_fast(
            model, kink_split=True, freq_grid='goal',
            frequency_quadrature='pchip')
        return time.perf_counter() - start, result, model
    finally:
        FS.solve_kernel = original


def main():
    parser = __import__('argparse').ArgumentParser()
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', default='docs/phase_increment_outer_round47_20260923.json')
    args = parser.parse_args()
    process = psutil.Process()
    selected = limit_affinity(process, args.threads)
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    rows = []
    for case in CASES:
        # Warm both dispatch paths before alternating timed measurements.
        timed(case, False)
        timed(case, True)
        base = []
        cand = []
        for i in range(args.repeats):
            first_candidate = (i % 2) == 1
            first = timed(case, first_candidate)[0]
            second = timed(case, not first_candidate)[0]
            if first_candidate:
                cand.append(first)
                base.append(second)
            else:
                base.append(first)
                cand.append(second)
        bm = statistics.median(base)
        cm = statistics.median(cand)
        rows.append({
            'case': case,
            'threads': args.threads,
            'repeats': args.repeats,
            'baseline_median_ms': bm * 1e3,
            'candidate_median_ms': cm * 1e3,
            'baseline_p95_ms': float(__import__('numpy').percentile(base, 95) * 1e3),
            'candidate_p95_ms': float(__import__('numpy').percentile(cand, 95) * 1e3),
            'candidate_over_baseline': cm / bm,
            'baseline_samples_ms': [v * 1e3 for v in base],
            'candidate_samples_ms': [v * 1e3 for v in cand],
        })
    payload = {
        'experiment': 'round47_phase_increment_full_outer_alternating',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': telemetry(process, threads=args.threads),
        'affinity_selected': selected,
        'settings': {'threads': args.threads, 'repeats': args.repeats,
                     'kink_split': True, 'freq_grid': 'goal',
                     'frequency_quadrature': 'pchip', 'reanchor': 32},
        'records': rows,
    }
    out = ROOT / args.out
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                   encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
