"""Round 53 scaled Riccati/log-amplitude feasibility screen.

This is a deterministic mathematical prototype.  It stores either ``x/y`` or
``y/x`` plus the logarithm and sign of the denominator, switching charts when
the current ratio exceeds one.  It is deliberately standalone and never
changes the production Cartesian transfer kernel.
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
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def _state_from_cartesian(x, y):
    if abs(y) >= abs(x):
        return (0, x / y if y != 0.0 else 0.0,
                math.log(max(abs(y), 1e-300)), 1.0 if y >= 0.0 else -1.0)
    return (1, y / x if x != 0.0 else 0.0,
            math.log(max(abs(x), 1e-300)), 1.0 if x >= 0.0 else -1.0)


def _cartesian_from_state(state):
    chart, ratio, logamp, sign = state
    amp = math.exp(logamp)
    if chart == 0:
        return sign * amp * ratio, sign * amp
    return sign * amp, sign * amp * ratio


def _riccati_step(state, w, h):
    chart, ratio, logamp, sign = state
    if chart == 0:
        f0 = -2.0 * ratio - w * (1.0 + ratio * ratio)
        ratio_mid = ratio + 0.5 * h * f0
        ratio_new = ratio + h * (
            -2.0 * ratio_mid - w * (1.0 + ratio_mid * ratio_mid))
        logamp_new = logamp + h * (1.0 + w * ratio_mid)
        if abs(ratio_new) <= 1.0:
            return chart, ratio_new, logamp_new, sign
        sign_new = sign * (1.0 if ratio_new >= 0.0 else -1.0)
        return (1, 1.0 / ratio_new,
                logamp_new + math.log(abs(ratio_new)), sign_new)
    f0 = 2.0 * ratio + w * (1.0 + ratio * ratio)
    ratio_mid = ratio + 0.5 * h * f0
    ratio_new = ratio + h * (2.0 * ratio_mid + w * (1.0 + ratio_mid * ratio_mid))
    logamp_new = logamp + h * (-1.0 - w * ratio_mid)
    if abs(ratio_new) <= 1.0:
        return chart, ratio_new, logamp_new, sign
    sign_new = sign * (1.0 if ratio_new >= 0.0 else -1.0)
    return (0, 1.0 / ratio_new,
            logamp_new + math.log(abs(ratio_new)), sign_new)


def _riccati_segment(state, z0, z1, h, phase_max):
    zmid = 0.5 * (z0 + z1)
    if phase_max <= 0.0 or zmid <= 0.0:
        n_sub = 1
    else:
        n_sub = max(1, int(math.ceil(h * math.exp(zmid) / phase_max)))
    hsub = h / n_sub
    dz = (z1 - z0) / n_sub
    for sub in range(n_sub):
        zs = z0 + dz * (sub + 0.5)
        state = _riccati_step(state, math.exp(zs), hsub)
    return state


@njit(inline="always")
def _riccati_step_numba(chart, ratio, logamp, sign, w, h):
    if chart == 0:
        ratio_mid = ratio + 0.5 * h * (-2.0 * ratio - w * (1.0 + ratio * ratio))
        ratio_new = ratio + h * (
            -2.0 * ratio_mid - w * (1.0 + ratio_mid * ratio_mid))
        logamp_new = logamp + h * (1.0 + w * ratio_mid)
    else:
        ratio_mid = ratio + 0.5 * h * (2.0 * ratio + w * (1.0 + ratio * ratio))
        ratio_new = ratio + h * (
            2.0 * ratio_mid + w * (1.0 + ratio_mid * ratio_mid))
        logamp_new = logamp + h * (-1.0 - w * ratio_mid)
    if abs(ratio_new) <= 1.0:
        return chart, ratio_new, logamp_new, sign
    sign_new = sign * (1.0 if ratio_new >= 0.0 else -1.0)
    if chart == 0:
        return 1, 1.0 / ratio_new, logamp_new + math.log(abs(ratio_new)), sign_new
    return 0, 1.0 / ratio_new, logamp_new + math.log(abs(ratio_new)), sign_new


@njit
def _riccati_chain(z_nodes, h, repeats):
    chart, ratio, logamp, sign = 0, -0.14285714285714285, math.log(0.875), -1.0
    for _ in range(repeats):
        chart, ratio, logamp, sign = 0, -0.14285714285714285, math.log(0.875), -1.0
        for k in range(len(z_nodes) - 1):
            z0 = z_nodes[k]
            z1 = z_nodes[k + 1]
            zmid = 0.5 * (z0 + z1)
            n_sub = 1
            if zmid > 0.0:
                n_sub = max(1, int(math.ceil(h * math.exp(zmid) / 0.25)))
            hsub = h / n_sub
            dz = (z1 - z0) / n_sub
            for sub in range(n_sub):
                w = math.exp(z0 + dz * (sub + 0.5))
                chart, ratio, logamp, sign = _riccati_step_numba(
                    chart, ratio, logamp, sign, w, hsub)
    return chart, ratio, logamp, sign


@njit
def _cartesian_chain(z_nodes, h, repeats):
    xh, yh = 0.125, -0.875
    for _ in range(repeats):
        xh, yh = 0.125, -0.875
        for k in range(len(z_nodes) - 1):
            xh, yh = FS._phase_segment(
                xh, yh, z_nodes[k], z_nodes[k + 1], h, 0.0)
    return xh, yh


def _timing(repeats=50):
    z_nodes = np.linspace(0.0, 5.0, 4097, dtype=np.float64)
    _riccati_chain(z_nodes, 0.005, 1)
    _cartesian_chain(z_nodes, 0.005, 1)
    base_times = []
    riccati_times = []
    for _ in range(repeats):
        started = time.perf_counter()
        base_state = _cartesian_chain(z_nodes, 0.005, 1)
        base_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        riccati_state = _riccati_chain(z_nodes, 0.005, 1)
        riccati_times.append(time.perf_counter() - started)
    return {
        "repeats": int(repeats),
        "baseline_median_ms": float(np.median(base_times) * 1e3),
        "candidate_median_ms": float(np.median(riccati_times) * 1e3),
        "baseline_p95_ms": float(np.percentile(base_times, 95) * 1e3),
        "candidate_p95_ms": float(np.percentile(riccati_times, 95) * 1e3),
        "candidate_over_baseline": float(
            np.median(riccati_times) / np.median(base_times)),
        "baseline_state": [float(v) for v in base_state],
        "candidate_state": [int(riccati_state[0]), float(riccati_state[1]),
                            float(riccati_state[2]), float(riccati_state[3])],
    }


def _path_screen(case_name, mode_limit=16, step_limit=256):
    _, common = _prepared(case_name, 2)
    _, phi_grid, _, s2, s2inv, j0s, z0s, *_ = common
    h = float(common[14])
    rows = []
    for mode in np.unique(np.linspace(0, len(j0s) - 1, mode_limit, dtype=int)):
        j0 = int(j0s[mode])
        phi0 = float(phi_grid[j0])
        xh, yh = 0.0, math.exp(float(z0s[mode])) * float(s2inv[j0])
        state = _state_from_cartesian(xh, yh)
        max_power_rel = 0.0
        finite = True
        steps = 0
        for k in range(j0, min(len(phi_grid) - 1, j0 + step_limit)):
            z0 = float(z0s[mode] + phi_grid[k] - phi0)
            z1 = float(z0s[mode] + phi_grid[k + 1] - phi0)
            xh, yh = FS._phase_segment(xh, yh, z0, z1, h, 0.25)
            state = _riccati_segment(state, z0, z1, h, 0.25)
            xc, yc = _cartesian_from_state(state)
            power = xh * xh + yh * yh
            candidate_power = xc * xc + yc * yc
            max_power_rel = max(
                max_power_rel,
                abs(candidate_power - power) / max(abs(power), 1e-300))
            finite = finite and all(math.isfinite(v) for v in (xc, yc, power))
            steps += 1
        rows.append({"mode": int(mode), "steps": steps, "finite": finite,
                     "max_power_relative": max_power_rel})
    return {
        "case": case_name,
        "finite_all": all(row["finite"] for row in rows),
        "max_power_relative": max(row["max_power_relative"] for row in rows),
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "experiment": "round53_scaled_riccati_log_amplitude_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "method": "ratio Riccati plus log-amplitude with deterministic chart switching",
        "timing": _timing(),
        "rows": [_path_screen(case) for case in CASES],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
