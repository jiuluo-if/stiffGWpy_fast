"""Round 66 standalone tangent correction for one outer background update.

The production outer map changes the background through the scalar ``DN_eff``
update.  This prototype propagates a first variation of the Cartesian transfer
state alongside the old state and predicts the updated endpoint without a
second transfer propagation.  It is deliberately fail-closed when the
discrete horizon start, tail handoff, kink split, or phase subdivision changes.
It is a feasibility screen only; production code is untouched.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _pin_current_process(requested=2):
    """Pin with the standard library; use optional psutil only on Windows."""
    try:
        available = sorted(os.sched_getaffinity(0))
        selected = set(available[:min(int(requested), len(available))])
        os.sched_setaffinity(0, selected)
        return sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        try:
            import psutil
            process = psutil.Process()
            available = list(process.cpu_affinity())
            selected = available[:min(int(requested), len(available))]
            process.cpu_affinity(selected)
            return list(process.cpu_affinity())
        except (ImportError, AttributeError, OSError):
            return []


@njit(cache=False, inline="always")
def _step_with_variation(z_value, delta_z, h_step, x_value, y_value,
                         dx_value, dy_value):
    """Apply midpoint transfer and its first variation with respect to z."""
    w = math.exp(z_value)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        angle = omega * h_step
        sine = math.sin(angle)
        cosine = math.cos(angle)
        c = cosine
        si = sine / omega
        domega = w2 / omega
        dc = -sine * h_step * domega
        dsi = domega * (omega * h_step * cosine - sine) / (omega * omega)
    else:
        x_argument = math.sqrt(1.0 - w2) * h_step
        c = 1.0 + 0.5 * x_argument * x_argument
        si = h_step * (1.0 + x_argument * x_argument / 6.0)
        dc = -w2 * h_step * h_step
        dsi = -w2 * h_step * h_step * h_step / 3.0
    a00 = c - si
    a01 = -w * si
    a10 = w * si
    a11 = c + si
    da00 = dc - dsi
    da01 = -w * (si + dsi)
    da10 = w * (si + dsi)
    da11 = dc + dsi
    next_x = a00 * x_value + a01 * y_value
    next_y = a10 * x_value + a11 * y_value
    next_dx = (a00 * dx_value + a01 * dy_value
               + delta_z * (da00 * x_value + da01 * y_value))
    next_dy = (a10 * dx_value + a11 * dy_value
               + delta_z * (da10 * x_value + da11 * y_value))
    return next_x, next_y, next_dx, next_dy


@njit(cache=False)
def _propagate(old_z, new_z, h_values, x_value, y_value,
               dx_value, dy_value, tangent):
    for index in range(old_z.size):
        if tangent:
            x_value, y_value, dx_value, dy_value = _step_with_variation(
                old_z[index], new_z[index] - old_z[index], h_values[index],
                x_value, y_value, dx_value, dy_value)
        else:
            x_value, y_value, _, _ = _step_with_variation(
                new_z[index], 0.0, h_values[index],
                x_value, y_value, 0.0, 0.0)
    return x_value, y_value, dx_value, dy_value


def _copy_snapshot(snapshot):
    return tuple(value.copy() if isinstance(value, np.ndarray) else value
                 for value in snapshot)


def _capture(case_name):
    saved = (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
             FS._OUTER_FULL_REUSE_SIGMA_TOL,
             FS._OUTER_FULL_REUSE_FHOR_TOL)
    snapshots = []
    try:
        def solve(*args, **kwargs):
            result = saved[2](*args, **kwargs)
            snapshots.append(_copy_snapshot(args))
            return result

        FS.solve_kernel = solve
        _pin_current_process(2)
        FS.apply_accuracy_mode("fast")
        FS.set_threads(2)
        FS._OUTER_FULL_REUSE_SIGMA_TOL = 0.0
        FS._OUTER_FULL_REUSE_FHOR_TOL = 0.0
        model = LCDM_SG(**CASES[case_name])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid="goal")
    finally:
        (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
         FS._OUTER_FULL_REUSE_SIGMA_TOL,
         FS._OUTER_FULL_REUSE_FHOR_TOL) = saved
    if len(snapshots) < 2:
        raise RuntimeError("outer run did not produce two propagation snapshots")
    return snapshots[0], snapshots[1]


def _append_segment(old_z0, old_z1, new_z0, new_z1, h_step, n_sub,
                    old_z, new_z, h_values):
    for sub in range(n_sub):
        fraction = (2.0 * sub + 1.0) / (2.0 * n_sub)
        old_z.append(old_z0 + (old_z1 - old_z0) * fraction)
        new_z.append(new_z0 + (new_z1 - new_z0) * fraction)
        h_values.append(h_step / n_sub)


def _path_pair(first, second, mode):
    old_phi = np.asarray(first[1], dtype=np.float64)
    new_phi = np.asarray(second[1], dtype=np.float64)
    old_mid = np.asarray(first[2], dtype=np.float64)
    new_mid = np.asarray(second[2], dtype=np.float64)
    old_j0 = int(first[5][mode])
    new_j0 = int(second[5][mode])
    if old_j0 != new_j0:
        raise ValueError("frequency start index changed")
    if (int(first[23]) != int(second[23])
            or abs(float(first[24]) - float(second[24])) > 1.0e-13
            or abs(float(first[25]) - float(second[25])) > 1.0e-13):
        raise ValueError("kink representation changed")
    old_z0 = float(first[6][mode])
    new_z0 = float(second[6][mode])
    z_tail = float(first[15])
    old_phi0 = float(old_phi[old_j0])
    new_phi0 = float(new_phi[new_j0])
    old_z = []
    new_z = []
    h_values = []
    k = old_j0
    while k < old_phi.size - 1:
        old_start = old_z0 + float(old_phi[k]) - old_phi0
        new_start = new_z0 + float(new_phi[k]) - new_phi0
        if old_start >= z_tail:
            break
        h_step = (float(first[19][k]) if first[19] is not None
                  else float(first[14]))
        old_mid_value = old_z0 + float(old_mid[k]) - old_phi0
        new_mid_value = new_z0 + float(new_mid[k]) - new_phi0
        n_sub_old = FS._phase_substeps(h_step, old_mid_value,
                                       float(first[21]))
        n_sub_new = FS._phase_substeps(h_step, new_mid_value,
                                       float(second[21]))
        if n_sub_old != n_sub_new:
            raise ValueError("phase subdivision changed")
        if k == int(first[23]) and 0.0 < float(first[24]) < 1.0:
            old_break = old_z0 + float(first[25]) - old_phi0
            new_break = new_z0 + float(second[25]) - new_phi0
            old_end = old_z0 + float(old_phi[k + 1]) - old_phi0
            new_end = new_z0 + float(new_phi[k + 1]) - new_phi0
            left_h = h_step * float(first[24])
            _append_segment(old_start, old_break, new_start, new_break,
                            left_h, FS._phase_substeps(
                                left_h, 0.5 * (old_start + old_break),
                                float(first[21])), old_z, new_z, h_values)
            right_h = h_step - left_h
            _append_segment(old_break, old_end, new_break, new_end,
                            right_h, FS._phase_substeps(
                                right_h, 0.5 * (old_break + old_end),
                                float(first[21])), old_z, new_z, h_values)
        else:
            old_end = old_z0 + float(old_phi[k + 1]) - old_phi0
            new_end = new_z0 + float(new_phi[k + 1]) - new_phi0
            _append_segment(old_start, old_end, new_start, new_end,
                            h_step, n_sub_old, old_z, new_z, h_values)
        k += 1
        if old_z0 + float(old_phi[k]) - old_phi0 >= z_tail:
            break
    if k >= old_phi.size - 1:
        kend = old_phi.size - 1
    else:
        kend = max(old_j0, k - 1)
    new_k = k
    if new_k >= new_phi.size - 1:
        new_kend = new_phi.size - 1
    else:
        new_kend = max(new_j0, new_k - 1)
    if kend != new_kend:
        raise ValueError("tail handoff index changed")
    return (np.asarray(old_z), np.asarray(new_z), np.asarray(h_values),
            old_j0, kend)


def _probe(case_name, mode=20, repeats=20):
    first, second = _capture(case_name)
    old_z, new_z, h_values, j0, kend = _path_pair(first, second, mode)
    old_y = math.exp(float(first[6][mode])) * float(first[4][j0])
    new_y = math.exp(float(second[6][mode])) * float(second[4][j0])
    predicted = _propagate(old_z, new_z, h_values, 0.0, old_y,
                           0.0, new_y - old_y, True)
    actual = _propagate(new_z, new_z, h_values, 0.0, new_y,
                        0.0, 0.0, False)
    predicted_state = np.asarray(predicted[:2]) + np.asarray(predicted[2:])
    actual_state = np.asarray(actual[:2])
    state_error = float(np.linalg.norm(predicted_state - actual_state)
                        / max(float(np.linalg.norm(actual_state)), 1.0e-300))
    base_times = []
    tangent_times = []
    for _ in range(repeats):
        start = time.perf_counter()
        _propagate(old_z, new_z, h_values, 0.0, old_y,
                   0.0, new_y - old_y, True)
        tangent_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _propagate(new_z, new_z, h_values, 0.0, new_y,
                   0.0, 0.0, False)
        base_times.append(time.perf_counter() - start)
    return {
        "case": case_name,
        "mode": int(mode),
        "eligible": True,
        "native_intervals": int(kend - j0),
        "substeps": int(h_values.size),
        "finite": bool(np.isfinite(predicted_state).all()),
        "state_relative_error": state_error,
        "tangent_over_direct": statistics.median(tangent_times) / statistics.median(base_times),
        "tangent_median_us": statistics.median(tangent_times) * 1.0e6,
        "direct_median_us": statistics.median(base_times) * 1.0e6,
        "z0_delta": float(second[6][mode] - first[6][mode]),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    affinity = _pin_current_process(2)
    rows = []
    for case_name in args.case or ["default", "lowT", "highT", "stiff", "high_kappa"]:
        try:
            rows.append(_probe(case_name))
        except (RuntimeError, ValueError, FloatingPointError) as exc:
            rows.append({"case": case_name, "eligible": False,
                         "reason": str(exc)})
    payload = {
        "schema_version": 1,
        "candidate": "event_aware_outer_tangent_endpoint_round66",
        "generated_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "resources": telemetry(workers=1, threads=2),
        "production_unchanged": True,
        "rows": rows,
        "semantics": (
            "First-variation endpoint correction over the production midpoint "
            "transfer; discrete event changes fail closed."),
    }
    payload["resources"]["affinity"] = affinity
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
