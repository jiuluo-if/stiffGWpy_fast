"""Standalone exact variable-coefficient Bessel transfer screen.

For one native interval with constant ``q = dz/dN`` and ``w=exp(z)``, the
scaled Cartesian system has a closed Bessel fundamental matrix.  This
prototype applies that matrix over a native interval, so it removes all
phase substeps mathematically.  It is diagnostic only: production remains
on the midpoint transfer map.
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
from scipy.special import jv, jvp, yv, yvp

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def _fundamental(z_value: float, q: float, n_value: float) -> np.ndarray:
    """Return the 2x2 Bessel fundamental matrix at ``(z, N)``."""
    if abs(q) < 1.0e-8:
        raise ValueError("constant-q Bessel chart is singular at q=0")
    alpha = q + 1.0 + 0.25 * q * q
    if alpha <= 0.0:
        raise ValueError("Bessel order is non-real for this q")
    abs_q = abs(q)
    order = math.sqrt(alpha) / abs_q
    argument = math.exp(z_value) / abs_q
    scale = math.exp(0.5 * q * n_value)
    j = float(jv(order, argument))
    y = float(yv(order, argument))
    dj = float(jvp(order, argument, 1))
    dy = float(yvp(order, argument, 1))
    # x = exp(qN/2) u(w), y = -(x_N + x) / w, w = exp(z).
    xj = scale * j
    xy = scale * y
    xnj = scale * (0.5 * q * j + q * argument * dj)
    xny = scale * (0.5 * q * y + q * argument * dy)
    w = math.exp(z_value)
    return np.array(((xj, xy), (-(xnj + xj) / w, -(xny + xy) / w)))


def bessel_transfer(x_value: float, y_value: float, z_start: float,
                    z_end: float, h_step: float) -> tuple[float, float]:
    """Apply the exact constant-q Bessel transfer for one native interval."""
    q = (z_end - z_start) / h_step
    left = _fundamental(z_start, q, 0.0)
    right = _fundamental(z_end, q, h_step)
    result = right @ np.linalg.solve(left, np.array((x_value, y_value)))
    return float(result[0]), float(result[1])


def _ode_transfer(x_value: float, y_value: float, z_start: float,
                  z_end: float, h_step: float) -> tuple[float, float]:
    q = (z_end - z_start) / h_step

    def rhs(n_value, state):
        w = math.exp(z_start + q * n_value)
        return (-state[0] - w * state[1], w * state[0] + state[1])

    sol = solve_ivp(rhs, (0.0, h_step), (x_value, y_value),
                    method="DOP853", rtol=2.0e-12, atol=2.0e-14,
                    max_step=h_step / 8.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(sol.y[0, -1]), float(sol.y[1, -1])


def _actual_intervals(case: str, mode: int = 20, limit: int = 32):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    rows = []
    for k in range(j0, min(len(model.Nv) - 1, j0 + limit)):
        z_start = z0 + float(phi[k]) - float(phi[j0])
        z_end = z0 + float(phi[k + 1]) - float(phi[j0])
        h_step = float(model.Nv[k + 1] - model.Nv[k])
        if z_end >= 5.0:
            break
        rows.append((z_start, z_end, h_step))
    return rows


def _probe(case: str, mode: int, repeats: int):
    intervals = _actual_intervals(case, mode)
    x0, y0 = 0.0, 1.0
    baseline = (x0, y0)
    candidate = (x0, y0)
    oracle = (x0, y0)
    for z_start, z_end, h_step in intervals:
        z_mid = 0.5 * (z_start + z_end)
        baseline = FS.scaled_step(baseline[0], baseline[1], z_mid, h_step)
        candidate = bessel_transfer(candidate[0], candidate[1], z_start, z_end, h_step)
        oracle = _ode_transfer(oracle[0], oracle[1], z_start, z_end, h_step)
    base_power = baseline[0] ** 2 + baseline[1] ** 2
    cand_power = candidate[0] ** 2 + candidate[1] ** 2
    oracle_power = oracle[0] ** 2 + oracle[1] ** 2
    baseline_times, candidate_times = [], []
    for _ in range(repeats):
        state = (x0, y0)
        started = time.perf_counter()
        for z_start, z_end, h_step in intervals:
            state = FS.scaled_step(state[0], state[1],
                                   0.5 * (z_start + z_end), h_step)
        baseline_times.append(time.perf_counter() - started)
        state = (x0, y0)
        started = time.perf_counter()
        for z_start, z_end, h_step in intervals:
            state = bessel_transfer(state[0], state[1], z_start, z_end, h_step)
        candidate_times.append(time.perf_counter() - started)
    bmed = statistics.median(baseline_times)
    cmed = statistics.median(candidate_times)
    return {
        "case": case,
        "mode": mode,
        "intervals": len(intervals),
        "baseline_median_ms": bmed * 1.0e3,
        "candidate_median_ms": cmed * 1.0e3,
        "candidate_over_baseline": cmed / bmed,
        "baseline_p95_ms": float(np.percentile(baseline_times, 95) * 1.0e3),
        "candidate_p95_ms": float(np.percentile(candidate_times, 95) * 1.0e3),
        "candidate_power_rel_vs_baseline": abs(cand_power - base_power)
        / max(abs(base_power), 1.0e-300),
        "candidate_power_rel_vs_dop853": abs(cand_power - oracle_power)
        / max(abs(oracle_power), 1.0e-300),
        "baseline_power_rel_vs_dop853": abs(base_power - oracle_power)
        / max(abs(oracle_power), 1.0e-300),
        "finite": bool(np.isfinite(candidate).all()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=25)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    rows = [_probe(case, 20, args.repeats) for case in CASES]
    payload = {
        "experiment": "round57_bessel_variable_coefficient_transfer",
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
