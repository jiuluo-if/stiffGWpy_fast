"""Standalone direct Ermakov--Kummer phase-function screen.

The candidate integrates the positive amplitude equation
``rho'' + q*rho = rho**-3`` and reconstructs the Cartesian transfer matrix
over the complete post-horizon oscillatory chain.  It is diagnostic only:
production subdivision, guards, assembly, and outer semantics are untouched.
Unlike the rejected Chebyshev and boundary-WKB screens, no polynomial fit or
right-boundary WKB handoff is used.
"""

from __future__ import annotations

import argparse
import json
import math
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
from scipy.integrate import solve_ivp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402


def _oscillatory_window(case, mode=20, z_lo=2.0, z_hi=5.0):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    z = z0 + np.asarray(phi, dtype=np.float64) - float(phi[j0])
    indices = np.flatnonzero((z >= z_lo) & (z <= z_hi))
    if indices.size < 16:
        raise ValueError(f"insufficient oscillatory window for {case}")
    start, end = int(indices[0]), int(indices[-1])
    nodes = np.asarray(model.Nv[start:end + 1], dtype=np.float64)
    sigma = np.asarray(model.sigma[start:end + 1], dtype=np.float64)
    z_nodes = z[start:end + 1]
    z_prime = 1.5 * sigma - 1.0
    z_second = np.gradient(z_prime, nodes)
    q = np.exp(2.0 * z_nodes) - 1.0 - z_prime + 0.5 * z_second - 0.25 * z_prime * z_prime
    if not np.isfinite(q).all() or float(np.min(q)) <= 0.0:
        raise ValueError(f"non-positive q in oscillatory window for {case}")
    return nodes, z_nodes, z_prime, q


def _phase_matrix(nodes, z, z_prime, q, max_step):
    q_nodes = np.asarray(q, dtype=np.float64)

    def rhs(n_value, state):
        q_value = float(np.interp(n_value, nodes, q_nodes))
        rho, rho_prime, theta = state
        return rho_prime, rho ** -3 - q_value * rho, rho ** -2

    h0 = float(nodes[1] - nodes[0])
    q_prime0 = float((q[1] - q[0]) / h0)
    rho0 = float(q[0] ** -0.25)
    state0 = [rho0, -0.25 * q_prime0 / q[0] * rho0, 0.0]
    solution = solve_ivp(
        rhs,
        (float(nodes[0]), float(nodes[-1])),
        state0,
        method="DOP853",
        rtol=2e-11,
        atol=2e-13,
        dense_output=True,
        max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    states = np.asarray(solution.sol([nodes[0], nodes[-1]]), dtype=np.float64)

    def state_matrix(index):
        rho, rho_prime, theta = states[:, index]
        theta_prime = rho ** -2
        scale = math.exp(0.5 * float(z[index]))
        c, s = math.cos(theta), math.sin(theta)
        u = np.array([[rho * c, rho * s],
                      [rho_prime * c - rho * theta_prime * s,
                       rho_prime * s + rho * theta_prime * c]])
        x = scale * u[0]
        x_prime = scale * (u[1] + 0.5 * float(z_prime[index]) * u[0])
        w = math.exp(float(z[index]))
        y = -(x_prime + x) / w
        return np.vstack((x, y))

    start_matrix = state_matrix(0)
    end_matrix = state_matrix(1)
    transfer = end_matrix @ np.linalg.inv(start_matrix)
    return transfer, solution, states


def _cartesian_matrix(nodes, z):
    def rhs(n_value, flat):
        w = math.exp(float(np.interp(n_value, nodes, z)))
        matrix = np.array([[-1.0, -w], [w, 1.0]])
        return (matrix @ flat.reshape(2, 2)).reshape(-1)

    solution = solve_ivp(
        rhs,
        (float(nodes[0]), float(nodes[-1])),
        np.eye(2).reshape(-1),
        method="DOP853",
        rtol=2e-11,
        atol=2e-13,
        max_step=float(np.max(np.diff(nodes))),
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[:, -1].reshape(2, 2), solution


def _probe(case):
    nodes, z, z_prime, q = _oscillatory_window(case)
    start = time.perf_counter()
    candidate, phase_solution, states = _phase_matrix(
        nodes, z, z_prime, q, float(np.max(np.diff(nodes))))
    candidate_s = time.perf_counter() - start
    start = time.perf_counter()
    reference, reference_solution = _cartesian_matrix(nodes, z)
    reference_s = time.perf_counter() - start
    stable, _, stable_states = _phase_matrix(
        nodes, z, z_prime, q, float(np.max(np.diff(nodes)) / 2.0))
    scale = max(float(np.linalg.norm(reference)), 1e-300)
    state_scale = max(float(np.linalg.norm(states[:, -1])), 1.0)
    return {
        "case": case,
        "native_intervals": int(nodes.size - 1),
        "finite": bool(np.isfinite(candidate).all()),
        "rho_min": float(np.min(states[0])),
        "candidate_vs_cartesian_relative_error": float(
            np.linalg.norm(candidate - reference) / scale),
        "max_step_stability_relative_error": float(
            np.linalg.norm(stable - candidate) / scale),
        "phase_state_stability_relative_error": float(
            np.linalg.norm(stable_states[:, -1] - states[:, -1]) / state_scale),
        "candidate_seconds": candidate_s,
        "reference_seconds": reference_s,
        "candidate_over_reference": candidate_s / reference_s,
        "candidate_nfev": int(phase_solution.nfev),
        "reference_nfev": int(reference_solution.nfev),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = {
        "candidate": "direct_ermakov_kummer_phase_function_round62",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "cases": [_probe(case) for case in (args.case or list(CASES))],
    }
    (ROOT / args.output).write_text(json.dumps(payload, indent=2) + "\n",
                                    encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
