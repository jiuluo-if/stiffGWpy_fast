"""Standalone Levin interaction-picture envelope feasibility screen.

For the transformed scalar equation ``u'' + omega(N)^2 u = 0``, this screen
factors out the rapid WKB carrier and integrates the resulting complex envelope
system.  It is a mathematical prototype only; production Cartesian transfer
is unchanged.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

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
        if end > start + 8 and not (start <= kink_index < end):
            return start, end
    return None


def _basis(z, z_prime, omega, omega_prime, phase):
    root = np.sqrt(omega)
    inv_root = 1.0 / root
    e_plus = np.exp(1j * phase)
    e_minus = np.exp(-1j * phase)
    da = -0.5 * omega_prime / omega * inv_root
    db = 0.5 * omega_prime / omega * root
    B = np.array([[inv_root * e_plus, inv_root * e_minus],
                  [1j * root * e_plus, -1j * root * e_minus]], dtype=complex)
    B_prime = np.array([
        [(da + 1j * omega * inv_root) * e_plus,
         (da - 1j * omega * inv_root) * e_minus],
        [(1j * db - omega * root) * e_plus,
         (-1j * db - omega * root) * e_minus],
    ], dtype=complex)
    return B, B_prime


def _screen_mode(Nv, sigma, z, start, end, degree):
    nodes = np.asarray(Nv[start:end + 1], dtype=np.float64)
    z_nodes = np.asarray(z[start:end + 1], dtype=np.float64)
    sigma_nodes = np.asarray(sigma[start:end + 1], dtype=np.float64)
    z_prime = 1.5 * sigma_nodes - 1.0
    z_second = np.gradient(z_prime, nodes)
    omega2 = (np.exp(2.0 * z_nodes) - 1.0 - z_prime
              + 0.5 * z_second - 0.25 * z_prime * z_prime)
    if not np.isfinite(omega2).all() or np.min(omega2) <= 0.0:
        return {"status": "turning_or_nonfinite"}
    omega = np.sqrt(omega2)
    omega_prime = np.gradient(omega, nodes)
    phase = cumulative_trapezoid(omega, nodes, initial=0.0)

    def rhs(n_value, state):
        q = float(np.interp(n_value, nodes, omega2))
        om = float(np.interp(n_value, nodes, omega))
        omp = float(np.interp(n_value, nodes, omega_prime))
        th = float(np.interp(n_value, nodes, phase))
        zz = float(np.interp(n_value, nodes, z_nodes))
        B, B_prime = _basis(zz, 0.0, om, omp, th)
        A = np.array([[0.0, 1.0], [-q, 0.0]], dtype=complex)
        C = np.linalg.solve(B, A @ B - B_prime)
        return (C @ state.reshape(2, 2)).reshape(-1)

    solution = solve_ivp(
        rhs, (float(nodes[0]), float(nodes[-1])), np.eye(2, dtype=complex).reshape(-1),
        method="DOP853", rtol=2e-10, atol=2e-12, dense_output=True,
        max_step=float(np.max(np.diff(nodes))),
    )
    if not solution.success:
        return {"status": "envelope_ode_failure", "message": solution.message}
    envelope = solution.sol(nodes).reshape(2, 2, -1)
    scale = 2.0 / (nodes[-1] - nodes[0])
    x = (nodes - nodes[0]) * scale - 1.0
    fitted = np.empty_like(envelope)
    fit_error = 0.0
    for i in range(2):
        for j in range(2):
            real_coeff = np.polynomial.chebyshev.chebfit(x, envelope[i, j].real, degree)
            imag_coeff = np.polynomial.chebyshev.chebfit(x, envelope[i, j].imag, degree)
            fitted[i, j] = (
                np.polynomial.chebyshev.chebval(x, real_coeff)
                + 1j * np.polynomial.chebyshev.chebval(x, imag_coeff))
            fit_error = max(
                fit_error,
                float(np.max(np.abs(fitted[i, j] - envelope[i, j])
                             / np.maximum(np.abs(envelope[i, j]), 1e-12))),
            )

    B_start, _ = _basis(float(z_nodes[0]), float(z_prime[0]),
                         float(omega[0]), float(omega_prime[0]), float(phase[0]))
    B_end, _ = _basis(float(z_nodes[-1]), float(z_prime[-1]),
                       float(omega[-1]), float(omega_prime[-1]), float(phase[-1]))
    approximate = B_end @ fitted[:, :, -1] @ np.linalg.inv(B_start)

    def direct_rhs(n_value, state):
        q = float(np.interp(n_value, nodes, omega2))
        return (np.array([[0.0, 1.0], [-q, 0.0]]) @ state.reshape(2, 2)).reshape(-1)

    direct = solve_ivp(
        direct_rhs, (float(nodes[0]), float(nodes[-1])), np.eye(2).reshape(-1),
        method="DOP853", rtol=2e-11, atol=2e-13,
        max_step=float(np.max(np.diff(nodes))),
    )
    if not direct.success:
        return {"status": "direct_ode_failure"}
    direct_transfer = direct.y[:, -1].reshape(2, 2)
    transfer_error = float(np.linalg.norm(approximate - direct_transfer)
                           / max(np.linalg.norm(direct_transfer), 1e-300))
    return {
        "status": "ok",
        "degree": int(degree),
        "nodes": int(nodes.size),
        "envelope_fit_relative_max": fit_error,
        "transfer_relative_error": transfer_error,
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
        row = _screen_mode(Nv, sigma, z, window[0], window[1], degree)
        row.update({"mode": int(mode), "window": [int(window[0]), int(window[1])]})
        rows.append(row)
    return {
        "candidate": "levin_interaction_envelope",
        "case": case_name,
        "production_unchanged": True,
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
        "candidate": "levin_interaction_envelope",
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
