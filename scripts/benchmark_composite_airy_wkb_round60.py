"""Standalone composite forbidden/Airy/WKB transfer screen.

The composite uses the production midpoint map before the simple turning
interval, the Round 59 Airy connection at that interval, and the Round 58
phase-integral WKB map afterwards.  It is diagnostic only; no production
dispatch or failure semantics are changed.
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

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _prepared  # noqa: E402
from scripts.benchmark_turning_airy_round59 import airy_transfer  # noqa: E402
from scripts.benchmark_wkb_phase_integral_round58 import wkb_transfer  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


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


def _intervals(case, mode=20, limit=4000):
    model, common = _prepared(case, 2)
    _, phi, _, _, _, j0s, z0s, *_ = common
    j0 = int(j0s[mode])
    z0 = float(z0s[mode])
    result = []
    for k in range(j0, min(len(model.Nv) - 1, j0 + limit)):
        left = z0 + float(phi[k]) - float(phi[j0])
        right = z0 + float(phi[k + 1]) - float(phi[j0])
        h_step = float(model.Nv[k + 1] - model.Nv[k])
        result.append((left, right, h_step))
    return result


def _probe(case):
    intervals = _intervals(case)
    baseline = (0.0, 1.0)
    candidate = (0.0, 1.0)
    oracle = (0.0, 1.0)
    counts = {"forbidden": 0, "turning_airy": 0, "oscillatory_wkb": 0}
    try:
        for left, right, h_step in intervals:
            midpoint = 0.5 * (left + right)
            baseline = FS.scaled_step(baseline[0], baseline[1], midpoint, h_step)
            q = (right - left) / h_step
            alpha = q + 1.0 + 0.25 * q * q
            q0 = math.exp(2.0 * left) - alpha
            q1 = math.exp(2.0 * right) - alpha
            if q0 * q1 <= 0.0:
                candidate = airy_transfer(candidate[0], candidate[1], left, right, h_step)
                counts["turning_airy"] += 1
            elif q1 < 0.0:
                candidate = FS.scaled_step(candidate[0], candidate[1], midpoint, h_step)
                counts["forbidden"] += 1
            else:
                candidate = wkb_transfer(candidate[0], candidate[1], left, right, h_step)
                counts["oscillatory_wkb"] += 1
            oracle = _oracle(oracle[0], oracle[1], left, right, h_step)
    except (ValueError, FloatingPointError):
        return {"case": case, "status": "candidate_failure", "counts": counts}
    bp = baseline[0] ** 2 + baseline[1] ** 2
    cp = candidate[0] ** 2 + candidate[1] ** 2
    op = oracle[0] ** 2 + oracle[1] ** 2
    return {
        "case": case,
        "status": "ok" if np.isfinite(candidate).all() else "nonfinite",
        "intervals": len(intervals),
        "counts": counts,
        "candidate_power_rel_vs_dop853": abs(cp - op) / max(abs(op), 1e-300),
        "baseline_power_rel_vs_dop853": abs(bp - op) / max(abs(op), 1e-300),
        "candidate_finite": bool(np.isfinite(candidate).all()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    payload = {
        "experiment": "round60_composite_airy_wkb_screen",
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
