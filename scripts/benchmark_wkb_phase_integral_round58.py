"""Standalone full-native-interval WKB phase-integral transfer screen.

This is different from the historical boundary-WKB handoff: it applies a
local WKB fundamental matrix on every native interval, using the analytic
phase integral of the constant-q transformed equation.  It is diagnostic
only and never changes the production kernel.
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

import numpy as np
from scipy.integrate import solve_ivp

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

from numba import njit  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(cache=False, inline="always")
def _primitive(w, alpha):
    root = math.sqrt(w * w - alpha)
    return root - math.sqrt(alpha) * math.atan(root / math.sqrt(alpha))


@njit(cache=False, inline="always")
def _basis(z, q, n_value, theta):
    alpha = q + 1.0 + 0.25 * q * q
    w = math.exp(z)
    omega2 = w * w - alpha
    if omega2 <= 0.0:
        return math.nan, math.nan, math.nan, math.nan
    omega = math.sqrt(omega2)
    amp = 1.0 / math.sqrt(omega)
    amp_prime = -0.5 * q * w * w / (omega * omega) * amp
    c = math.cos(theta)
    s = math.sin(theta)
    u1 = amp * c
    u2 = amp * s
    up1 = amp_prime * c - amp * omega * s
    up2 = amp_prime * s + amp * omega * c
    scale = math.exp(0.5 * q * n_value)
    x1 = scale * u1
    x2 = scale * u2
    xp1 = scale * (0.5 * q * u1 + up1)
    xp2 = scale * (0.5 * q * u2 + up2)
    return x1, x2, -(xp1 + x1) / w, -(xp2 + x2) / w


@njit(cache=False, inline="always")
def wkb_transfer(x_value, y_value, z_start, z_end, h_step):
    q = (z_end - z_start) / h_step
    alpha = q + 1.0 + 0.25 * q * q
    w_start = math.exp(z_start)
    w_end = math.exp(z_end)
    if abs(q) < 1.0e-10 or w_start * w_start <= alpha or w_end * w_end <= alpha:
        return math.nan, math.nan
    phase = (_primitive(w_end, alpha) - _primitive(w_start, alpha)) / q
    a11, a12, a21, a22 = _basis(z_start, q, 0.0, 0.0)
    b11, b12, b21, b22 = _basis(z_end, q, h_step, phase)
    det = a11 * a22 - a12 * a21
    c1 = (a22 * x_value - a12 * y_value) / det
    c2 = (-a21 * x_value + a11 * y_value) / det
    return b11 * c1 + b12 * c2, b21 * c1 + b22 * c2


@njit(cache=False)
def _wkb_loop(z_start, z_end, h_steps):
    x_value, y_value = 0.0, 1.0
    for k in range(len(h_steps)):
        x_value, y_value = wkb_transfer(
            x_value, y_value, z_start[k], z_end[k], h_steps[k])
    return x_value, y_value


@njit(cache=False)
def _cart_loop(z_start, z_end, h_steps):
    x_value, y_value = 0.0, 1.0
    for k in range(len(h_steps)):
        x_value, y_value = FS.scaled_step(
            x_value, y_value, 0.5 * (z_start[k] + z_end[k]), h_steps[k])
    return x_value, y_value


def _oracle(x_value, y_value, z_start, z_end, h_step):
    q = (z_end - z_start) / h_step

    def rhs(n_value, state):
        w = math.exp(z_start + q * n_value)
        return (-state[0] - w * state[1], w * state[0] + state[1])

    sol = solve_ivp(rhs, (0.0, h_step), (x_value, y_value), method="DOP853",
                    rtol=2.0e-12, atol=2.0e-14, max_step=h_step / 8.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(sol.y[0, -1]), float(sol.y[1, -1])


def _intervals(case, mode=20, limit=32):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    z_start, z_end, h_steps = [], [], []
    for k in range(j0, min(len(model.Nv) - 1, j0 + limit)):
        left = z0 + float(phi[k]) - float(phi[j0])
        right = z0 + float(phi[k + 1]) - float(phi[j0])
        if right >= 5.0:
            break
        z_start.append(left)
        z_end.append(right)
        h_steps.append(float(model.Nv[k + 1] - model.Nv[k]))
    return (np.asarray(z_start), np.asarray(z_end), np.asarray(h_steps))


def _probe(case, repeats):
    z_start, z_end, h_steps = _intervals(case)
    baseline = _cart_loop(z_start, z_end, h_steps)
    candidate = _wkb_loop(z_start, z_end, h_steps)
    finite = bool(np.isfinite(candidate).all())
    if not finite:
        return {
            "case": case,
            "intervals": int(len(h_steps)),
            "finite": False,
            "rejection": "turning_point_or_forbidden_interval",
            "baseline_power_rel_vs_dop853": None,
            "candidate_power_rel_vs_dop853": None,
            "baseline_median_ms": None,
            "candidate_median_ms": None,
            "candidate_over_baseline": None,
            "baseline_p95_ms": None,
            "candidate_p95_ms": None,
        }
    oracle = (0.0, 1.0)
    for left, right, h_step in zip(z_start, z_end, h_steps):
        oracle = _oracle(oracle[0], oracle[1], float(left), float(right), float(h_step))
    base_power = baseline[0] ** 2 + baseline[1] ** 2
    cand_power = candidate[0] ** 2 + candidate[1] ** 2
    oracle_power = oracle[0] ** 2 + oracle[1] ** 2
    base_times, cand_times = [], []
    for _ in range(repeats):
        started = time.perf_counter()
        _cart_loop(z_start, z_end, h_steps)
        base_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        _wkb_loop(z_start, z_end, h_steps)
        cand_times.append(time.perf_counter() - started)
    bmed = statistics.median(base_times)
    cmed = statistics.median(cand_times)
    return {
        "case": case,
        "intervals": int(len(h_steps)),
        "finite": finite,
        "baseline_median_ms": bmed * 1.0e3,
        "candidate_median_ms": cmed * 1.0e3,
        "candidate_over_baseline": cmed / bmed,
        "baseline_p95_ms": float(np.percentile(base_times, 95) * 1.0e3),
        "candidate_p95_ms": float(np.percentile(cand_times, 95) * 1.0e3),
        "candidate_power_rel_vs_dop853": abs(cand_power - oracle_power)
        / max(abs(oracle_power), 1.0e-300),
        "baseline_power_rel_vs_dop853": abs(base_power - oracle_power)
        / max(abs(oracle_power), 1.0e-300),
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
        "experiment": "round58_full_native_wkb_phase_integral",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "settings": vars(args),
        "rows": rows,
    }
    (ROOT / args.out).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
