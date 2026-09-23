"""Round 52 cheap closed-form defect-bound screen.

The Round 51 embedded matrix defect was too expensive.  This diagnostic uses a
conservative endpoint-only Magnus defect bound for the same log-mean two-step
map.  It is a mathematical guard, not a learned predictor, and is never wired
into the production solver.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from scripts.benchmark_uniform_asymptotic_defect_round51 import (  # noqa: E402
    _uniform_block_defect,
    _uniform_block_map_numba,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402

BOUND_SAFETY = 16.0


def _defect_bound(z0: float, z1: float, h: float) -> float:
    """Conservative endpoint-only bound for a two-interval block.

    ``h`` is the native half-block step, so the total block width is ``2*h``.
    For ``A(w)=A0+w*A1``, ``||[A0,A1]||_2=2`` and a linear endpoint path gives
    the first omitted Magnus term bounded by ``H^2*|dw|/6``.  The exponential
    stability factor and a fixed safety factor cover higher terms in this
    diagnostic.  The bound is deliberately conservative; it is not calibrated
    to pass a target number of blocks.
    """
    if z0 == z1:
        return 0.0
    w0 = math.exp(z0)
    w1 = math.exp(z1)
    h_total = 2.0 * h
    w_max = max(w0, w1)
    matrix_norm = math.sqrt(2.0 + 2.0 * w_max * w_max)
    first_omitted = h_total * h_total * abs(w1 - w0) / 6.0
    return BOUND_SAFETY * math.exp(h_total * matrix_norm) * first_omitted


@njit(inline="always")
def _defect_bound_numba(z0, z1, h):
    if z0 == z1:
        return 0.0
    w0 = math.exp(z0)
    w1 = math.exp(z1)
    h_total = 2.0 * h
    w_max = max(w0, w1)
    matrix_norm = math.sqrt(2.0 + 2.0 * w_max * w_max)
    first_omitted = h_total * h_total * abs(w1 - w0) / 6.0
    return BOUND_SAFETY * math.exp(h_total * matrix_norm) * first_omitted


@njit
def _transfer_timing_bound(z_nodes, h, repeats, threshold):
    xh, yh = 0.125, -0.875
    accepted = 0
    fallback = 0
    for _ in range(repeats):
        xh, yh = 0.125, -0.875
        for k in range(0, len(z_nodes) - 2, 2):
            bound = _defect_bound_numba(z_nodes[k], z_nodes[k + 2], h)
            if bound <= threshold:
                xh, yh = _uniform_block_map_numba(
                    xh, yh, z_nodes[k], z_nodes[k + 2], 2.0 * h)
                accepted += 1
            else:
                xh, yh = FS._phase_segment(
                    xh, yh, z_nodes[k], z_nodes[k + 1], h, 0.25)
                xh, yh = FS._phase_segment(
                    xh, yh, z_nodes[k + 1], z_nodes[k + 2], h, 0.25)
                fallback += 1
    return xh, yh, accepted, fallback


@njit
def _transfer_timing_baseline(z_nodes, h, repeats):
    xh, yh = 0.125, -0.875
    for _ in range(repeats):
        xh, yh = 0.125, -0.875
        for k in range(0, len(z_nodes) - 2, 2):
            xh, yh = FS._phase_segment(
                xh, yh, z_nodes[k], z_nodes[k + 1], h, 0.25)
            xh, yh = FS._phase_segment(
                xh, yh, z_nodes[k + 1], z_nodes[k + 2], h, 0.25)
    return xh, yh


def _path_summary(case_name, threshold, mode_limit=16, step_limit=256):
    _, common = _prepared(case_name, 2)
    _, phi_grid, _, _, _, j0s, z0s, *_ = common
    h = float(common[14])
    kink_index = common[19]
    accepted = 0
    fallback = 0
    max_bound_to_actual = 0.0
    max_bound = 0.0
    modes = np.unique(np.linspace(0, len(j0s) - 1, mode_limit, dtype=int))
    for mode in modes:
        j0 = int(j0s[mode])
        phi0 = float(phi_grid[j0])
        k = j0
        steps = 0
        while k + 1 < len(phi_grid) and steps < step_limit:
            z0 = float(z0s[mode] + phi_grid[k] - phi0)
            if z0 >= 5.0:
                break
            if (k + 2 < len(phi_grid) and k != kink_index and k + 1 != kink_index):
                z1 = float(z0s[mode] + phi_grid[k + 2] - phi0)
                bound = _defect_bound(z0, z1, h)
                actual = _uniform_block_defect(z0, z1, h)
                max_bound = max(max_bound, bound)
                max_bound_to_actual = max(
                    max_bound_to_actual, actual / max(bound, 1.0e-300))
                if bound <= threshold and 2.0 * h * max(
                        _log_mean_endpoint(z0, z1), 1.0) <= 0.5 and z1 < 5.0:
                    accepted += 1
                    k += 2
                    steps += 2
                    continue
            fallback += 1
            k += 1
            steps += 1
    return {
        "case": case_name,
        "threshold": float(threshold),
        "accepted_blocks": int(accepted),
        "fallback_steps": int(fallback),
        "max_bound": float(max_bound),
        "max_actual_to_bound": float(max_bound_to_actual),
    }


def _log_mean_endpoint(z0, z1):
    mid = 0.5 * (z0 + z1)
    half = 0.5 * (z1 - z0)
    if half == 0.0:
        return math.exp(mid)
    if abs(half) < 1.0e-4:
        h2 = half * half
        return math.exp(mid) * (1.0 + h2 / 6.0 + h2 * h2 / 120.0)
    return math.exp(mid) * math.sinh(half) / half


def _timing(repeats=50, threshold=1e-8):
    z_nodes = np.linspace(0.0, 4.0, 4097, dtype=np.float64)
    _transfer_timing_bound(z_nodes, 0.005, 1, threshold)
    base_times = []
    bound_times = []
    for _ in range(repeats):
        started = time.perf_counter()
        _transfer_timing_baseline(z_nodes, 0.005, 1)
        base_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        state = _transfer_timing_bound(z_nodes, 0.005, 1, threshold)
        bound_times.append(time.perf_counter() - started)
    return {
        "repeats": int(repeats),
        "baseline_median_ms": float(np.median(base_times) * 1e3),
        "baseline_p95_ms": float(np.percentile(base_times, 95) * 1e3),
        "bound_guard_median_ms": float(np.median(bound_times) * 1e3),
        "bound_guard_p95_ms": float(np.percentile(bound_times, 95) * 1e3),
        "bound_guard_over_full_chain": float(np.median(bound_times) / max(np.median(base_times), 1e-300)),
        "accepted_blocks": int(state[2]),
        "fallback_steps": int(state[3]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--threshold", type=float, nargs="+", default=[1e-8, 1e-6])
    args = parser.parse_args()
    rows = [_path_summary(case, threshold)
            for threshold in args.threshold for case in CASES]
    payload = {
        "experiment": "round52_uniform_endpoint_defect_bound_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "bound_safety_factor": BOUND_SAFETY,
        "rows": rows,
        "timing": {str(threshold): _timing(threshold=threshold)
                   for threshold in args.threshold},
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
