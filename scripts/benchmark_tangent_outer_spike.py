"""Standalone tangent-propagation cost screen for the outer loop.

This is deliberately a feasibility prototype, not an outer-loop replacement.
It derives the analytic Jacobian of one constant-z transfer with respect to z,
then carries a tangent state alongside the production state.  The screen asks
whether that extra response state could plausibly cost less than the second
hard-case propagation before any final-spectrum semantics are attempted.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _make_args, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(inline="always", cache=True)
def _tangent_step(xh, yh, tx, ty, z_mid, h, dz):
    w = math.exp(z_mid)
    w2 = w * w
    if w2 > 1.0 + 1.0e-12:
        omega = math.sqrt(w2 - 1.0)
        q = omega * h
        c = math.cos(q)
        sine = math.sin(q)
        si = sine / omega
        domega = w2 / omega
        dc = -sine * h * domega
        dsi = (c * h * domega / omega
               - sine * domega / (omega * omega))
    elif w2 < 1.0 - 1.0e-12:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
        domega = -w2 / omega
        dx = h * domega
        dc = x * dx
        dsi = h * x * dx / 3.0
    else:
        c = 1.0
        si = h
        dc = -h * h
        dsi = -h * h * h / 3.0
    ax = c - si
    bx = -w * si
    ay = w * si
    by = c + si
    d_bx = -(w * si + w * dsi)
    d_ay = w * si + w * dsi
    d_ax = dc - dsi
    d_by = dc + dsi
    next_x = ax * xh + bx * yh
    next_y = ay * xh + by * yh
    tangent_x = ax * tx + bx * ty + dz * (d_ax * xh + d_bx * yh)
    tangent_y = ay * tx + by * ty + dz * (d_ay * xh + d_by * yh)
    return next_x, next_y, tangent_x, tangent_y


def _local_jacobian_error():
    z = 1.2
    h = 0.005
    eps = 1.0e-6
    tangent = _tangent_step(0.2, -0.1, 0.0, 0.0, z, h, 1.0)
    plus = FS.scaled_step(0.2, -0.1, z + eps, h)
    minus = FS.scaled_step(0.2, -0.1, z - eps, h)
    fd = ((plus[0] - minus[0]) / (2.0 * eps),
          (plus[1] - minus[1]) / (2.0 * eps))
    return float(max(abs(tangent[2] - fd[0]), abs(tangent[3] - fd[1])))


@njit(inline="always", cache=True)
def _tangent_segment(xh, yh, tx, ty, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * math.exp(z_mid) / phase_max))
        if n_sub < 1:
            n_sub = 1
    h_sub = h_step / n_sub
    if n_sub == 1:
        return _tangent_step(xh, yh, tx, ty, z_mid, h_step, 1.0e-6)
    dz = (z_end - z_start) / n_sub
    for sub in range(n_sub):
        z_sub_mid = z_start + (sub + 0.5) * dz
        xh, yh, tx, ty = _tangent_step(
            xh, yh, tx, ty, z_sub_mid, h_sub, 1.0e-6)
    return xh, yh, tx, ty


@njit(parallel=True, cache=True)
def propagate_tangent(Nv, Phi_grid, S2inv, j0s, z0s, h_arr, h,
                      z_tail, phase_max, kink_index, kink_fraction, phi_re,
                      output):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0) * S2inv[j0]
        tx, ty = 0.0, 0.0
        k = j0
        zz = z0
        while k < len(Nv) - 1 and zz < z_tail:
            h_step = h_arr[k] if h_arr is not None else h
            z_node = z0 + Phi_grid[k] - Phi0
            z_end = z0 + Phi_grid[k + 1] - Phi0
            if k == kink_index and 0.0 < kink_fraction < 1.0:
                z_break = z0 + phi_re - Phi0
                h_left = h_step * kink_fraction
                xh, yh, tx, ty = _tangent_segment(
                    xh, yh, tx, ty, z_node, z_break, h_left, phase_max)
                xh, yh, tx, ty = _tangent_segment(
                    xh, yh, tx, ty, z_break, z_end,
                    h_step - h_left, phase_max)
            else:
                xh, yh, tx, ty = _tangent_segment(
                    xh, yh, tx, ty, z_node, z_end, h_step, phase_max)
            k += 1
            zz = z0 + Phi_grid[k] - Phi0
        output[mode, 0] = xh
        output[mode, 1] = yh
        output[mode, 2] = tx
        output[mode, 3] = ty


def _run_cost(case_name, threads=2, repeats=15):
    model, common = _prepared(case_name, threads)
    args = _make_args(common)
    Nv, Phi_grid, _, _, S2inv, j0s, z0s = common[:7]
    h_arr, phase_max, kink_index, kink_fraction, phi_re = (
        common[17], common[18], common[19], common[20], common[21])
    tangent_output = np.zeros((len(j0s), 4), dtype=np.float64)
    FS.solve_kernel(*args)
    propagate_tangent(
        Nv, Phi_grid, S2inv, j0s, z0s, h_arr, common[14], common[15],
        phase_max, kink_index, kink_fraction, phi_re, tangent_output)
    baseline_times = []
    tangent_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        FS.solve_kernel(*_make_args(common))
        baseline_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        propagate_tangent(
            Nv, Phi_grid, S2inv, j0s, z0s, h_arr, common[14], common[15],
            phase_max, kink_index, kink_fraction, phi_re, tangent_output)
        tangent_times.append(time.perf_counter() - start)
    return {
        "case": case_name,
        "threads": int(threads),
        "repeats": int(repeats),
        "baseline_full_kernel_median_ms": statistics.median(baseline_times) * 1e3,
        "tangent_no_assembly_median_ms": statistics.median(tangent_times) * 1e3,
        "tangent_over_full_baseline": statistics.median(tangent_times) / statistics.median(baseline_times),
        "baseline_full_kernel_p95_ms": float(np.percentile(baseline_times, 95) * 1e3),
        "tangent_no_assembly_p95_ms": float(np.percentile(tangent_times, 95) * 1e3),
        "local_jacobian_error": _local_jacobian_error(),
        "tangent_finite": bool(np.isfinite(tangent_output).all()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=15)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "tangent_outer_response_cost_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "cases": [_run_cost(name, args.threads, args.repeats)
                  for name in (args.case or list(CASES))],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
