"""Standalone scalar-adjoint cost screen for the outer loop.

An adjoint can estimate one scalar observable without carrying a tangent for
every output component, but the user-facing API still returns the spectrum.
This prototype measures the propagation cost and verifies the transpose map;
it is not a production outer-loop implementation.
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
def _adjoint_step(lx, ly, z_mid, h):
    w = math.exp(z_mid)
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
    return (c - si) * lx + w * si * ly, -w * si * lx + (c + si) * ly


def _local_adjoint_error():
    z = 1.2
    h = 0.005
    lx, ly = 0.2, -0.1
    actual = _adjoint_step(lx, ly, z, h)
    first = FS.scaled_step(1.0, 0.0, z, h)
    second = FS.scaled_step(0.0, 1.0, z, h)
    expected = (first[0] * lx - second[0] * ly,
                -first[1] * lx + second[1] * ly)
    return float(max(abs(actual[0] - expected[0]), abs(actual[1] - expected[1])))


@njit(inline="always", cache=True)
def _adjoint_segment(lx, ly, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * math.exp(z_mid) / phase_max))
        if n_sub < 1:
            n_sub = 1
    h_sub = h_step / n_sub
    if n_sub == 1:
        return _adjoint_step(lx, ly, z_mid, h_step)
    dz = (z_end - z_start) / n_sub
    for sub in range(n_sub):
        z_sub_mid = z_start + (sub + 0.5) * dz
        lx, ly = _adjoint_step(lx, ly, z_sub_mid, h_sub)
    return lx, ly


@njit(parallel=True, cache=True)
def propagate_adjoint(Nv, Phi_grid, j0s, z0s, h_arr, h, z_tail,
                      phase_max, kink_index, kink_fraction, phi_re, output):
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Phi0 = Phi_grid[j0]
        lx, ly = 1.0, 1.0
        k = j0
        zz = z0
        while k < len(Nv) - 1 and zz < z_tail:
            h_step = h_arr[k] if h_arr is not None else h
            z_node = z0 + Phi_grid[k] - Phi0
            z_end = z0 + Phi_grid[k + 1] - Phi0
            if k == kink_index and 0.0 < kink_fraction < 1.0:
                z_break = z0 + phi_re - Phi0
                h_left = h_step * kink_fraction
                lx, ly = _adjoint_segment(
                    lx, ly, z_node, z_break, h_left, phase_max)
                lx, ly = _adjoint_segment(
                    lx, ly, z_break, z_end, h_step - h_left, phase_max)
            else:
                lx, ly = _adjoint_segment(
                    lx, ly, z_node, z_end, h_step, phase_max)
            k += 1
            zz = z0 + Phi_grid[k] - Phi0
        output[mode, 0] = lx
        output[mode, 1] = ly


def _run_cost(case_name, threads=2, repeats=15):
    model, common = _prepared(case_name, threads)
    args = _make_args(common)
    Nv, Phi_grid, _, _, _, j0s, z0s = common[:7]
    h_arr, phase_max, kink_index, kink_fraction, phi_re = (
        common[17], common[18], common[19], common[20], common[21])
    output = np.zeros((len(j0s), 2), dtype=np.float64)
    FS.solve_kernel(*args)
    propagate_adjoint(
        Nv, Phi_grid, j0s, z0s, h_arr, common[14], common[15], phase_max,
        kink_index, kink_fraction, phi_re, output)
    baseline_times = []
    adjoint_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        FS.solve_kernel(*_make_args(common))
        baseline_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        propagate_adjoint(
            Nv, Phi_grid, j0s, z0s, h_arr, common[14], common[15], phase_max,
            kink_index, kink_fraction, phi_re, output)
        adjoint_times.append(time.perf_counter() - start)
    return {
        "case": case_name,
        "threads": int(threads),
        "repeats": int(repeats),
        "baseline_full_kernel_median_ms": statistics.median(baseline_times) * 1e3,
        "adjoint_no_assembly_median_ms": statistics.median(adjoint_times) * 1e3,
        "adjoint_over_full_baseline": statistics.median(adjoint_times) / statistics.median(baseline_times),
        "baseline_full_kernel_p95_ms": float(np.percentile(baseline_times, 95) * 1e3),
        "adjoint_no_assembly_p95_ms": float(np.percentile(adjoint_times, 95) * 1e3),
        "local_adjoint_error": _local_adjoint_error(),
        "adjoint_finite": bool(np.isfinite(output).all()),
        "full_spectrum_reconstructable": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=15)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "observable_adjoint_outer_cost_screen",
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
