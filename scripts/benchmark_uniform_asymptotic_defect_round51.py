"""Round 51 standalone uniform-integral two-interval defect screen.

For a linear-in-grid coordinate ``z`` the coefficient ``w=exp(z)`` has an
analytic interval average (the logarithmic mean).  This diagnostic combines
two native intervals into one constant-coefficient transfer using that
average.  An embedded operator defect compares the one-block map with the
product of the two interval-average maps.  The candidate is never installed
in production; it is only eligible for later work reduction if the defect and
the comparison with a finer deterministic transfer both pass.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(inline="always")
def _log_mean_exp_numba(z0, z1):
    mid = 0.5 * (z0 + z1)
    half = 0.5 * (z1 - z0)
    if half == 0.0:
        return math.exp(mid)
    ah = abs(half)
    if ah < 1.0e-4:
        h2 = half * half
        correction = 1.0 + h2 / 6.0 + (h2 * h2) / 120.0
    else:
        correction = math.sinh(half) / half
    return math.exp(mid) * correction


@njit(inline="always")
def _uniform_block_map_numba(xh, yh, z0, z1, h):
    w = _log_mean_exp_numba(z0, z1)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


@njit(inline="always")
def _map_entries_numba(w, h):
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return c - si, -w * si, w * si, c + si


@njit(inline="always")
def _uniform_block_defect_numba(z0, z1, h):
    if z0 == z1:
        return 0.0
    zm = 0.5 * (z0 + z1)
    one = _map_entries_numba(_log_mean_exp_numba(z0, z1), 2.0 * h)
    left = _map_entries_numba(_log_mean_exp_numba(z0, zm), h)
    right = _map_entries_numba(_log_mean_exp_numba(zm, z1), h)
    ar, br, cr, dr = right
    al, bl, cl, dl = left
    two = (
        ar * al + br * cl,
        ar * bl + br * dl,
        cr * al + dr * cl,
        cr * bl + dr * dl,
    )
    numerator = 0.0
    denominator = 0.0
    for a, b in zip(one, two):
        numerator += (a - b) * (a - b)
        denominator += b * b
    return math.sqrt(numerator) / max(math.sqrt(denominator), 1.0e-300)


@njit
def _transfer_timing(z_nodes, h, repeats, candidate):
    xh, yh = 0.125, -0.875
    for _ in range(repeats):
        xh, yh = 0.125, -0.875
        for k in range(0, len(z_nodes) - 2, 2):
            if candidate:
                xh, yh = _uniform_block_map_numba(
                    xh, yh, z_nodes[k], z_nodes[k + 2], 2.0 * h)
            else:
                xh, yh = FS._phase_segment(
                    xh, yh, z_nodes[k], z_nodes[k + 1], h, 0.25)
                xh, yh = FS._phase_segment(
                    xh, yh, z_nodes[k + 1], z_nodes[k + 2], h, 0.25)
    return xh, yh


@njit
def _transfer_timing_guarded(z_nodes, h, repeats, threshold):
    xh, yh = 0.125, -0.875
    accepted = 0
    fallback = 0
    for _ in range(repeats):
        xh, yh = 0.125, -0.875
        for k in range(0, len(z_nodes) - 2, 2):
            defect = _uniform_block_defect_numba(
                z_nodes[k], z_nodes[k + 2], h)
            if defect <= threshold:
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


def _log_mean_exp(z0: float, z1: float) -> float:
    """Return ``(exp(z1)-exp(z0))/(z1-z0)`` without cancellation."""
    mid = 0.5 * (z0 + z1)
    half = 0.5 * (z1 - z0)
    if half == 0.0:
        return math.exp(mid)
    ah = abs(half)
    if ah < 1.0e-4:
        h2 = half * half
        correction = 1.0 + h2 / 6.0 + (h2 * h2) / 120.0
    else:
        correction = math.sinh(half) / half
    return math.exp(mid) * correction


def _map_entries(w: float, h: float):
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return c - si, -w * si, w * si, c + si


def _uniform_block_map(xh, yh, z0, z1, h):
    wbar = _log_mean_exp(z0, z1)
    a, b, c, d = _map_entries(wbar, h)
    return a * xh + b * yh, c * xh + d * yh


def _compose(right, left):
    ar, br, cr, dr = right
    al, bl, cl, dl = left
    return (
        ar * al + br * cl,
        ar * bl + br * dl,
        cr * al + dr * cl,
        cr * bl + dr * dl,
    )


def _uniform_block_defect(z0, z1, h):
    """Relative Frobenius defect of one block versus two half blocks."""
    if z0 == z1:
        return 0.0
    zm = 0.5 * (z0 + z1)
    one = _map_entries(_log_mean_exp(z0, z1), 2.0 * h)
    left = _map_entries(_log_mean_exp(z0, zm), h)
    right = _map_entries(_log_mean_exp(zm, z1), h)
    two = _compose(right, left)
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(one, two)))
    denominator = max(math.sqrt(sum(a * a for a in two)), 1.0e-300)
    return numerator / denominator


def _fine_segment(xh, yh, z0, z1, h, subdivisions=8):
    dz = (z1 - z0) / subdivisions
    hh = h / subdivisions
    for i in range(subdivisions):
        left = z0 + i * dz
        right = z0 + (i + 1) * dz
        xh, yh = FS._phase_segment(xh, yh, left, right, hh, 0.0)
    return xh, yh


def _run_case(case_name, threshold, phase_cap, mode_limit=16, step_limit=256):
    model, common = _prepared(case_name, 2)
    _, phi_grid, _, s2, s2inv, j0s, z0s, *_ = common
    h_step = float(common[14])
    kink_index = common[19]
    rows = []
    for mode in np.unique(np.linspace(0, len(j0s) - 1, mode_limit, dtype=int)):
        j0 = int(j0s[mode])
        phi0 = float(phi_grid[j0])
        base_x, base_y = 0.0, math.exp(float(z0s[mode])) * float(s2inv[j0])
        cand_x, cand_y = base_x, base_y
        fine_x, fine_y = base_x, base_y
        accepted = 0
        fallback = 0
        max_defect = 0.0
        max_power_rel = 0.0
        max_fine_rel = 0.0
        k = j0
        steps = 0
        while k < len(phi_grid) - 1 and steps < step_limit:
            z_left = float(z0s[mode] + phi_grid[k] - phi0)
            z_right = float(z0s[mode] + phi_grid[k + 1] - phi0)
            if z_left >= 5.0:
                break
            if (k + 2 < len(phi_grid) and k != kink_index and k + 1 != kink_index):
                z_two = float(z0s[mode] + phi_grid[k + 2] - phi0)
                defect = _uniform_block_defect(z_left, z_two, h_step)
                phase = (h_step + h_step) * max(_log_mean_exp(z_left, z_two), 1.0)
                max_defect = max(max_defect, defect)
                use_block = defect <= threshold and phase <= phase_cap and z_two < 5.0
            else:
                use_block = False
                defect = 0.0
                z_two = z_right
            if use_block:
                z_next = z_two
                z_mid = z_right
                base_x, base_y = FS._phase_segment(
                    base_x, base_y, z_left, z_mid, h_step, 0.25)
                base_x, base_y = FS._phase_segment(
                    base_x, base_y, z_mid, z_next, h_step, 0.25)
                fine_x, fine_y = _fine_segment(
                    fine_x, fine_y, z_left, z_mid, h_step)
                fine_x, fine_y = _fine_segment(
                    fine_x, fine_y, z_mid, z_next, h_step)
                cand_x, cand_y = _uniform_block_map(
                    cand_x, cand_y, z_left, z_two, h_step + h_step)
                next_k = k + 2
                accepted += 1
            else:
                base_x, base_y = FS._phase_segment(
                    base_x, base_y, z_left, z_right, h_step, 0.25)
                fine_x, fine_y = _fine_segment(
                    fine_x, fine_y, z_left, z_right, h_step)
                cand_x, cand_y = FS._phase_segment(
                    cand_x, cand_y, z_left, z_right, h_step, 0.25)
                next_k = k + 1
                fallback += 1
            base_power = base_x * base_x + base_y * base_y
            cand_power = cand_x * cand_x + cand_y * cand_y
            fine_power = fine_x * fine_x + fine_y * fine_y
            max_power_rel = max(
                max_power_rel,
                abs(cand_power - base_power) / max(abs(base_power), 1.0e-300),
            )
            max_fine_rel = max(
                max_fine_rel,
                abs(cand_power - fine_power) / max(abs(fine_power), 1.0e-300),
            )
            k = next_k
            steps += 1 if not use_block else 2
        rows.append({
            "mode": int(mode),
            "j0": j0,
            "accepted_blocks": accepted,
            "fallback_steps": fallback,
            "native_intervals": steps,
            "max_embedded_defect": max_defect,
            "max_power_relative_to_production": max_power_rel,
            "max_power_relative_to_fine_midpoint": max_fine_rel,
        })
    return {
        "case": case_name,
        "threshold": float(threshold),
        "phase_cap": float(phase_cap),
        "mode_limit": int(mode_limit),
        "step_limit": int(step_limit),
        "rows": rows,
        "accepted_blocks": int(sum(row["accepted_blocks"] for row in rows)),
        "fallback_steps": int(sum(row["fallback_steps"] for row in rows)),
        "max_embedded_defect": max(row["max_embedded_defect"] for row in rows),
        "max_power_relative_to_production": max(
            row["max_power_relative_to_production"] for row in rows),
        "max_power_relative_to_fine_midpoint": max(
            row["max_power_relative_to_fine_midpoint"] for row in rows),
    }


def _run_transfer_timing(repeats=50):
    import time

    z_nodes = np.linspace(0.0, 4.0, 4097, dtype=np.float64)
    _transfer_timing(z_nodes, 0.005, 1, False)
    _transfer_timing(z_nodes, 0.005, 1, True)
    baseline_times = []
    candidate_times = []
    guarded_times = []
    guarded_state = None
    guarded_accepted = 0
    guarded_fallback = 0
    for _ in range(repeats):
        started = time.perf_counter()
        base_state = _transfer_timing(z_nodes, 0.005, 1, False)
        baseline_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        candidate_state = _transfer_timing(z_nodes, 0.005, 1, True)
        candidate_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        guarded_state = _transfer_timing_guarded(z_nodes, 0.005, 1, 1e-8)
        guarded_times.append(time.perf_counter() - started)
        guarded_accepted = guarded_state[2]
        guarded_fallback = guarded_state[3]
    return {
        "repeats": int(repeats),
        "baseline_median_ms": float(np.median(baseline_times) * 1.0e3),
        "candidate_median_ms": float(np.median(candidate_times) * 1.0e3),
        "baseline_p95_ms": float(np.percentile(baseline_times, 95) * 1.0e3),
        "candidate_p95_ms": float(np.percentile(candidate_times, 95) * 1.0e3),
        "guarded_median_ms": float(np.median(guarded_times) * 1.0e3),
        "guarded_p95_ms": float(np.percentile(guarded_times, 95) * 1.0e3),
        "candidate_over_baseline": float(
            np.median(candidate_times) / np.median(baseline_times)),
        "guarded_over_baseline": float(
            np.median(guarded_times) / np.median(baseline_times)),
        "baseline_state": [float(v) for v in base_state],
        "candidate_state": [float(v) for v in candidate_state],
        "guarded_state": [float(v) for v in guarded_state[:2]],
        "guarded_accepted_blocks": int(guarded_accepted),
        "guarded_fallback_blocks": int(guarded_fallback),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--threshold", type=float, nargs="+", default=[1e-8, 1e-6, 1e-4])
    parser.add_argument("--phase-cap", type=float, default=0.5)
    args = parser.parse_args()
    rows = [
        _run_case(case, threshold, args.phase_cap)
        for threshold in args.threshold
        for case in CASES
    ]
    payload = {
        "experiment": "round51_uniform_integral_two_interval_defect_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "method": "log-mean exp coefficient with embedded one-block/two-block matrix defect",
        "settings": vars(args),
        "transfer_timing": _run_transfer_timing(),
        "rows": rows,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
