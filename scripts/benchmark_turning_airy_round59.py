"""Standalone simple-turning-point Airy connection screen.

The Round 58 WKB map failed because the real path crosses ``Omega^2=0``.
This prototype replaces only one crossing native interval by the Airy
fundamental matrix obtained from a linearized transformed potential.  It is
an accuracy feasibility screen; production propagation is untouched.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import airy

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402


def _basis(z, q, n_value, slope, n_turn):
    sign = 1.0 if slope >= 0.0 else -1.0
    beta = abs(slope) ** (1.0 / 3.0)
    arg = -sign * beta * (n_value - n_turn)
    ai, aip, bi, bip = airy(arg)
    scale_arg = -sign * beta
    u = (float(ai), float(bi))
    up = (scale_arg * float(aip), scale_arg * float(bip))
    w = math.exp(z)
    pref = math.exp(0.5 * q * n_value)
    rows = []
    for value, derivative in zip(u, up):
        x = pref * value
        xp = pref * (0.5 * q * value + derivative)
        rows.append((x, -(xp + x) / w))
    return np.array(((rows[0][0], rows[1][0]),
                     (rows[0][1], rows[1][1])))


def airy_transfer(x_value, y_value, z_start, z_end, h_step):
    q = (z_end - z_start) / h_step
    alpha = q + 1.0 + 0.25 * q * q
    q0 = math.exp(2.0 * z_start) - alpha
    q1 = math.exp(2.0 * z_end) - alpha
    slope = (q1 - q0) / h_step
    if q0 * q1 > 0.0 or abs(slope) < 1.0e-12:
        raise ValueError("interval is not a simple turning interval")
    n_turn = -q0 / slope
    left = _basis(z_start, q, 0.0, slope, n_turn)
    right = _basis(z_end, q, h_step, slope, n_turn)
    return tuple(right @ np.linalg.solve(left, np.array((x_value, y_value))))


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


def _turning_interval(case, mode=20):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    for k in range(j0, len(model.Nv) - 1):
        left = z0 + float(phi[k]) - float(phi[j0])
        right = z0 + float(phi[k + 1]) - float(phi[j0])
        h_step = float(model.Nv[k + 1] - model.Nv[k])
        q = (right - left) / h_step
        alpha = q + 1.0 + 0.25 * q * q
        q0 = math.exp(2.0 * left) - alpha
        q1 = math.exp(2.0 * right) - alpha
        if q0 * q1 <= 0.0:
            return left, right, h_step, q0, q1
    return None


def _probe(case):
    interval = _turning_interval(case)
    if interval is None:
        return {"case": case, "status": "no_simple_turning_interval"}
    left, right, h_step, q0, q1 = interval
    baseline = (0.0, 1.0)
    midpoint = 0.5 * (left + right)
    from stiffgwpy_fast import fast_sgwb as FS
    baseline = FS.scaled_step(baseline[0], baseline[1], midpoint, h_step)
    candidate = airy_transfer(0.0, 1.0, left, right, h_step)
    oracle = _oracle(0.0, 1.0, left, right, h_step)
    cp = candidate[0] ** 2 + candidate[1] ** 2
    bp = baseline[0] ** 2 + baseline[1] ** 2
    op = oracle[0] ** 2 + oracle[1] ** 2
    return {
        "case": case,
        "status": "ok",
        "z_start": left,
        "z_end": right,
        "h_step": h_step,
        "potential_start": q0,
        "potential_end": q1,
        "candidate_power_rel_vs_dop853": abs(cp - op) / max(abs(op), 1e-300),
        "baseline_power_rel_vs_dop853": abs(bp - op) / max(abs(op), 1e-300),
        "finite": bool(np.isfinite(candidate).all()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    payload = {
        "experiment": "round59_simple_turning_airy_connection",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "settings": vars(args),
        "rows": [_probe(case) for case in CASES],
    }
    (ROOT / args.out).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
