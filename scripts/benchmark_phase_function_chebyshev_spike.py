"""Standalone Kummer phase-function/Chebyshev feasibility screen.

This is intentionally not a production solver.  It solves the positive
Kummer amplitude equation on smooth post-horizon background windows, fits the
amplitude and accumulated phase with Chebyshev polynomials, and evaluates the
Kummer residual from the fitted derivatives.  Raw ``y'/y`` Riccati
coordinates, turning-point windows, and the production path are not used.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from scipy.integrate import cumulative_trapezoid, solve_ivp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402


def _window(z, nv, kink_index):
    for z_lo, z_hi in ((2.0, 4.0), (1.5, 3.0), (2.0, 3.0)):
        left = np.flatnonzero(z >= z_lo)
        right = np.flatnonzero(z >= z_hi)
        if not left.size or not right.size:
            continue
        start = int(left[0])
        end = int(right[0])
        if end <= start + 8:
            continue
        if kink_index < 0 or not (start <= kink_index < end):
            return start, end
    return None


def _fit_mode(Nv, sigma, z, start, end, degree):
    nodes = np.asarray(Nv[start:end + 1], dtype=np.float64)
    z_nodes = np.asarray(z[start:end + 1], dtype=np.float64)
    sigma_nodes = np.asarray(sigma[start:end + 1], dtype=np.float64)
    z_prime = 1.5 * sigma_nodes - 1.0
    z_second = np.gradient(z_prime, nodes)
    omega2 = np.exp(2.0 * z_nodes) - 1.0 - z_prime + 0.5 * z_second - 0.25 * z_prime * z_prime
    if not np.isfinite(omega2).all() or float(np.min(omega2)) <= 0.0:
        return {"status": "turning_or_nonfinite"}
    def rhs(n_value, state):
        q = float(np.interp(n_value, nodes, omega2))
        return [state[1], -q * state[0], state[3], -q * state[2]]

    solution = solve_ivp(
        rhs, (float(nodes[0]), float(nodes[-1])), [1.0, 0.0, 0.0, 1.0],
        method="DOP853", rtol=1e-9, atol=1e-11, dense_output=True,
        max_step=float(np.max(np.diff(nodes)))
    )
    if not solution.success:
        return {"status": "ode_failure", "message": solution.message}
    fundamental = np.asarray(solution.sol(nodes), dtype=np.float64)
    rho = np.sqrt(fundamental[0] ** 2 + fundamental[2] ** 2)
    if not np.isfinite(rho).all() or float(np.min(rho)) <= 0.0:
        return {"status": "nonpositive_solution"}
    phase = cumulative_trapezoid(1.0 / (rho * rho), nodes, initial=0.0)
    scale = 2.0 / (nodes[-1] - nodes[0])
    x = (nodes - nodes[0]) * scale - 1.0
    rho_coeff = np.polynomial.chebyshev.chebfit(x, rho, degree)
    phase_coeff = np.polynomial.chebyshev.chebfit(x, phase, degree)
    rho_fit = np.polynomial.chebyshev.chebval(x, rho_coeff)
    rho_prime_fit = np.polynomial.chebyshev.chebval(
        x, np.polynomial.chebyshev.chebder(rho_coeff)) * scale
    rho_second_fit = np.polynomial.chebyshev.chebval(
        x, np.polynomial.chebyshev.chebder(rho_coeff, 2)) * scale * scale
    residual = rho_second_fit + omega2 * rho_fit - 1.0 / (rho_fit ** 3)
    phase_fit = np.polynomial.chebyshev.chebval(x, phase_coeff)
    return {
        "status": "ok",
        "nodes": int(nodes.size),
        "degree": int(degree),
        "rho_min": float(np.min(rho)),
        "rho_fit_min": float(np.min(rho_fit)),
        "rho_fit_relative_max": float(np.max(np.abs(rho_fit - rho) / rho)),
        "phase_fit_absolute_max": float(np.max(np.abs(phase_fit - phase))),
        "chebyshev_residual_max": float(np.max(
            np.abs(residual) / np.maximum(np.abs(omega2 * rho), 1e-12))),
    }


def _run_case(name, mode_count=4, degree=12):
    model, common = _prepared(name, 2)
    Nv = np.asarray(model.Nv, dtype=np.float64)
    Phi = np.asarray(common[1], dtype=np.float64)
    j0s = np.asarray(common[5], dtype=np.int64)
    z0s = np.asarray(common[6], dtype=np.float64)
    sigma = np.asarray(model.sigma, dtype=np.float64)
    kink_index = int(common[20])
    picks = np.unique(np.linspace(0, len(j0s) - 1, mode_count, dtype=int))
    rows = []
    for mode in picks:
        z = z0s[mode] + Phi - Phi[j0s[mode]]
        window = _window(z, Nv, kink_index)
        if window is None:
            rows.append({"mode": int(mode), "status": "no_smooth_window"})
            continue
        row = _fit_mode(Nv, sigma, z, window[0], window[1], degree)
        row.update({"mode": int(mode), "window": [int(window[0]), int(window[1])]})
        rows.append(row)
    usable = [row for row in rows if row["status"] == "ok"]
    return {
        "case": name,
        "mode_count_requested": int(mode_count),
        "usable_modes": len(usable),
        "positive_modes": sum(row.get("rho_fit_min", 0.0) > 0.0 for row in usable),
        "chebyshev_residual_max": max(
            (row["chebyshev_residual_max"] for row in usable), default=float("inf")
        ),
        "rho_fit_relative_max": max(
            (row["rho_fit_relative_max"] for row in usable), default=float("inf")
        ),
        "phase_fit_absolute_max": max(
            (row["phase_fit_absolute_max"] for row in usable), default=float("inf")
        ),
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--mode-count", type=int, default=4)
    parser.add_argument("--degree", type=int, default=12)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = args.case or list(CASES)
    payload = {
        "candidate": "kummer_phase_function_chebyshev_feasibility",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "cases": [_run_case(name, args.mode_count, args.degree) for name in cases],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
