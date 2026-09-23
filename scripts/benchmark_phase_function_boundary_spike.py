"""Standalone boundary-conditioned Kummer phase-function feasibility spike.

The previous phase screen used an arbitrary fundamental basis.  This prototype
instead selects a complex scalar solution with an outgoing-WKB amplitude and
derivative at the smooth right boundary, integrates it backward, and tests
whether the resulting Kummer amplitude is compressible on a useful interval.
It is diagnostic only; production propagation is untouched.
"""

from __future__ import annotations

import argparse
import json
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


def _state_matrix(rho, rho_prime, phase, z, z_prime):
    theta_prime = 1.0 / (rho * rho)
    scale = np.exp(0.5 * z)
    c = np.cos(phase)
    s = np.sin(phase)
    u = np.array([[rho * c, rho * s],
                  [rho_prime * c - rho * theta_prime * s,
                   rho_prime * s + rho * theta_prime * c]])
    x = scale * u[0]
    x_prime = scale * (u[1] + 0.5 * z_prime * u[0])
    w = np.exp(z)
    y = -(x_prime + x) / w
    return np.vstack((x, y))


def _cartesian_transfer(nodes, z_nodes):
    def rhs(n_value, state):
        z = float(np.interp(n_value, nodes, z_nodes))
        w = np.exp(z)
        matrix = np.array([[-1.0, -w], [w, 1.0]])
        return (matrix @ state.reshape(2, 2)).reshape(-1)

    initial = np.eye(2).reshape(-1)
    solution = solve_ivp(
        rhs, (float(nodes[0]), float(nodes[-1])), initial,
        method="DOP853", rtol=2e-11, atol=2e-13,
        max_step=float(np.max(np.diff(nodes))),
    )
    if not solution.success:
        return None
    return solution.y[:, -1].reshape(2, 2)


def _fit_boundary_conditioned(Nv, sigma, z, start, end, degree):
    nodes = np.asarray(Nv[start:end + 1], dtype=np.float64)
    z_nodes = np.asarray(z[start:end + 1], dtype=np.float64)
    sigma_nodes = np.asarray(sigma[start:end + 1], dtype=np.float64)
    z_prime = 1.5 * sigma_nodes - 1.0
    z_second = np.gradient(z_prime, nodes)
    omega2 = (np.exp(2.0 * z_nodes) - 1.0 - z_prime
              + 0.5 * z_second - 0.25 * z_prime * z_prime)
    if not np.isfinite(omega2).all() or float(np.min(omega2)) <= 0.0:
        return {"status": "turning_or_nonfinite"}

    def rhs(n_value, state):
        q = float(np.interp(n_value, nodes, omega2))
        return [state[1], -q * state[0], state[3], -q * state[2]]

    omega = np.sqrt(omega2)
    omega_prime = np.gradient(omega, nodes)
    amp = float(omega[-1] ** -0.5)
    amp_prime = float(-0.5 * omega_prime[-1] / omega[-1] * amp)
    solution = solve_ivp(
        rhs, (float(nodes[-1]), float(nodes[0])),
        [amp, 0.0, amp_prime, omega[-1] * amp],
        method="DOP853", rtol=1e-10, atol=1e-12,
        dense_output=True, max_step=float(np.max(np.diff(nodes))),
    )
    if not solution.success:
        return {"status": "ode_failure", "message": solution.message}

    state = np.asarray(solution.sol(nodes), dtype=np.float64)
    rho = np.hypot(state[0], state[2])
    phase = cumulative_trapezoid(1.0 / (rho * rho), nodes, initial=0.0)
    if not np.isfinite(rho).all() or float(np.min(rho)) <= 0.0:
        return {"status": "nonpositive_solution"}

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
    z_prime_fit = np.gradient(z_nodes, nodes)
    start_matrix = _state_matrix(
        float(rho_fit[0]), float(rho_prime_fit[0]), float(phase_fit[0]),
        float(z_nodes[0]), float(z_prime_fit[0]))
    end_matrix = _state_matrix(
        float(rho_fit[-1]), float(rho_prime_fit[-1]), float(phase_fit[-1]),
        float(z_nodes[-1]), float(z_prime_fit[-1]))
    phase_transfer = end_matrix @ np.linalg.inv(start_matrix)
    cartesian_transfer = _cartesian_transfer(nodes, z_nodes)
    if cartesian_transfer is None:
        return {"status": "cartesian_ode_failure"}
    transfer_scale = max(float(np.linalg.norm(cartesian_transfer)), 1e-300)
    transfer_error = float(np.linalg.norm(
        phase_transfer - cartesian_transfer) / transfer_scale)
    return {
        "status": "ok",
        "nodes": int(nodes.size),
        "degree": int(degree),
        "rho_min": float(np.min(rho)),
        "rho_fit_min": float(np.min(rho_fit)),
        "rho_fit_relative_max": float(np.max(np.abs(rho_fit - rho) / rho)),
        "phase_fit_absolute_max": float(np.max(np.abs(phase_fit - phase))),
        "kummer_residual_max": float(np.max(
            np.abs(residual) / np.maximum(np.abs(omega2 * rho), 1e-12))),
        "cartesian_transfer_relative_error": transfer_error,
        "window_z_length": float(z_nodes[-1] - z_nodes[0]),
        "window_N_length": float(nodes[-1] - nodes[0]),
    }


def _run_case(case_name, mode_count=4, degree=12):
    model, common = _prepared(case_name, 2)
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
        row = _fit_boundary_conditioned(
            Nv, sigma, z, window[0], window[1], degree)
        row.update({"mode": int(mode), "window": [int(window[0]), int(window[1])]})
        rows.append(row)
    usable = [row for row in rows if row["status"] == "ok"]
    return {
        "candidate": "boundary_conditioned_kummer_phase",
        "case": case_name,
        "production_unchanged": True,
        "mode_count_requested": int(mode_count),
        "usable_modes": len(usable),
        "positive_modes": sum(row.get("rho_fit_min", 0.0) > 0.0 for row in usable),
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--mode-count", type=int, default=4)
    parser.add_argument("--degree", type=int, default=12)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "boundary_conditioned_kummer_phase",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "cases": [_run_case(name, args.mode_count, args.degree)
                  for name in (args.case or list(CASES))],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
