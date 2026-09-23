"""Standalone exact-map exp recurrence screen for phase substeps.

The candidate keeps the production subdivision and every sin/cos/sqrt transfer
unchanged.  Only the positive coefficient ``w=exp(z)`` is generated once at
the first substep and advanced by a geometric recurrence.  This is distinct
from the rejected phase-angle recurrence: no phase angle is approximated or
re-anchored.
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
from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(inline="always", cache=False)
def _scaled_step_with_w(xh, yh, w, h):
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


@njit(cache=False)
def _segment_exp_recurrence(xh, yh, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    n_sub = FS._phase_substeps(h_step, z_mid, phase_max)
    if n_sub == 1:
        return FS.scaled_step(xh, yh, z_mid, h_step)
    h_sub = h_step / n_sub
    dz_sub = (z_end - z_start) / n_sub
    w = math.exp(z_start + 0.5 * dz_sub)
    w_ratio = math.exp(dz_sub)
    for _ in range(n_sub):
        xh, yh = _scaled_step_with_w(xh, yh, w, h_sub)
        w *= w_ratio
    return xh, yh


@njit(cache=False)
def _path_baseline(z_start, z_end, h_steps, phase_max):
    xh, yh = 0.0, 1.0
    for k in range(len(h_steps)):
        xh, yh = FS._phase_segment(
            xh, yh, z_start[k], z_end[k], h_steps[k], phase_max)
    return xh, yh


@njit(cache=False)
def _path_candidate(z_start, z_end, h_steps, phase_max):
    xh, yh = 0.0, 1.0
    for k in range(len(h_steps)):
        xh, yh = _segment_exp_recurrence(
            xh, yh, z_start[k], z_end[k], h_steps[k], phase_max)
    return xh, yh


def _actual_path(case, mode=20):
    model, common = _prepared(case, 2)
    _, phi, phi_mid, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    starts, ends, steps = [], [], []
    total_substeps = 0
    for k in range(j0, len(model.Nv) - 1):
        left = z0 + float(phi[k]) - float(phi[j0])
        right = z0 + float(phi[k + 1]) - float(phi[j0])
        if left >= 5.0:
            break
        h_step = float(model.Nv[k + 1] - model.Nv[k])
        n_sub = int(FS._phase_substeps(h_step, float(
            z0 + phi_mid[k] - phi[j0]), 0.25))
        starts.append(left)
        ends.append(right)
        steps.append(h_step)
        total_substeps += n_sub
    return (np.asarray(starts, dtype=np.float64),
            np.asarray(ends, dtype=np.float64),
            np.asarray(steps, dtype=np.float64), total_substeps)


def _probe(case, repeats=25):
    z_start, z_end, h_steps, total_substeps = _actual_path(case)
    baseline = _path_baseline(z_start, z_end, h_steps, 0.25)
    candidate = _path_candidate(z_start, z_end, h_steps, 0.25)
    base_power = baseline[0] ** 2 + baseline[1] ** 2
    cand_power = candidate[0] ** 2 + candidate[1] ** 2
    base_times, cand_times = [], []
    for _ in range(repeats):
        start = time.perf_counter()
        _path_baseline(z_start, z_end, h_steps, 0.25)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _path_candidate(z_start, z_end, h_steps, 0.25)
        cand_times.append(time.perf_counter() - start)
    bmed = statistics.median(base_times)
    cmed = statistics.median(cand_times)
    return {
        "case": case,
        "native_segments": int(len(h_steps)),
        "phase_substeps": int(total_substeps),
        "finite": bool(np.isfinite(candidate).all()),
        "baseline_median_ms": bmed * 1.0e3,
        "candidate_median_ms": cmed * 1.0e3,
        "candidate_over_baseline": cmed / bmed,
        "baseline_p95_ms": float(np.percentile(base_times, 95) * 1.0e3),
        "candidate_p95_ms": float(np.percentile(cand_times, 95) * 1.0e3),
        "power_relative_error": abs(cand_power - base_power)
        / max(abs(base_power), 1.0e-300),
        "component_relative_error": max(
            abs(candidate[0] - baseline[0]) / max(abs(baseline[0]), 1.0e-300),
            abs(candidate[1] - baseline[1]) / max(abs(baseline[1]), 1.0e-300)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=25)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    rows = [_probe(case, args.repeats) for case in CASES]
    payload = {
        "experiment": "round61_exp_coefficient_recurrence",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "settings": vars(args),
        "phase_max": 0.25,
        "rows": rows,
    }
    (ROOT / args.out).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
