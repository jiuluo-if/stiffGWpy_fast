"""Formal total-runtime A/B for the standalone phase recurrence prototype."""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import solve_kernel_recurrence  # noqa: E402
from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _run(point, candidate, threads):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    model = LCDM_SG(**CASES[point])
    original = FS.solve_kernel
    try:
        if candidate:
            FS.solve_kernel = solve_kernel_recurrence
        start = time.perf_counter()
        result = FS.SGWB_iter_fast(
            model, kink_split=True, freq_grid='goal',
            frequency_quadrature='pchip')
        elapsed = time.perf_counter() - start
    finally:
        FS.solve_kernel = original
    status = 'ok' if result is not None else getattr(
        model, 'fast_failure_reason', 'failure')
    if status != 'ok':
        return elapsed, {'status': status}
    return elapsed, {
        'status': status,
        'n_freq': int(len(model.f)),
        'spectrum': model.log10OmegaGW.copy(),
        'DN_gw': float(model.DN_gw[-1]),
    }


def run_case(point, repeats, threads):
    for candidate in (False, True):
        _run(point, candidate, threads)
    base_times = []
    candidate_times = []
    base_last = candidate_last = None
    for _ in range(repeats):
        elapsed, base_last = _run(point, False, threads)
        base_times.append(elapsed)
        elapsed, candidate_last = _run(point, True, threads)
        candidate_times.append(elapsed)
    if base_last['status'] != candidate_last['status']:
        return {
            'point': point,
            'status': 'status_mismatch',
            'base_status': base_last['status'],
            'candidate_status': candidate_last['status'],
        }
    row = {
        'point': point,
        'status': base_last['status'],
        'repeats': repeats,
        'baseline_median_ms': statistics.median(base_times) * 1e3,
        'candidate_median_ms': statistics.median(candidate_times) * 1e3,
        'baseline_p95_ms': float(np.percentile(base_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
        'candidate_over_baseline': statistics.median(candidate_times)
        / statistics.median(base_times),
    }
    if base_last['status'] == 'ok':
        delta = np.abs(candidate_last['spectrum'] - base_last['spectrum'])
        row.update({
            'n_freq': base_last['n_freq'],
            'spectrum_delta_dex_p95': float(np.percentile(delta, 95)),
            'spectrum_delta_dex_max': float(np.max(delta)),
            'DN_rel': abs(candidate_last['DN_gw'] - base_last['DN_gw'])
            / max(abs(base_last['DN_gw']), 1e-300),
        })
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', nargs='+', choices=sorted(CASES),
                        default=['default', 'highT', 'stiff', 'high_kappa', 'lowT'])
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--out', default='docs/phase_recurrence_total_runtime_20260917.json')
    args = parser.parse_args(argv)
    rows = [run_case(point, args.repeats, args.threads) for point in args.points]
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_phase_exp_recurrence_total_runtime',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'points': args.points, 'repeats': args.repeats,
                     'threads': args.threads, 'kink_split': True,
                     'frequency_quadrature': 'pchip'},
        'records': rows,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
