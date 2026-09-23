"""Round 67: conservative transfer-perturbation certificate.

This is a standalone feasibility screen for the existing goal-grid outer
reuse path.  It does not change production code.  For an unchanged event
sequence it bounds the endpoint state difference by propagating a norm bound
through the same midpoint transfer matrices used by ``solve_kernel``::

    d[i+1] <= ||M_new|| d[i] + ||M_new - M_old|| ||x_old[i]||.

The Frobenius norm is intentionally conservative.  A bound is useful only if
it is small for every frequency contributing to the spectrum; a single
ineligible event or a vacuous bound rejects the candidate fail-closed.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
from pathlib import Path

import numpy as np

try:
    from scripts._resource_budget import telemetry
    from scripts.benchmark_outer_sensitivity_round66 import (
        _capture,
        _path_pair,
        _pin_current_process,
    )
except ImportError:  # pragma: no cover
    from _resource_budget import telemetry
    from benchmark_outer_sensitivity_round66 import (
        _capture,
        _path_pair,
        _pin_current_process,
    )

ROOT = Path(__file__).resolve().parents[1]

STRICT_POWER_RELATIVE_TARGET = 1.0e-12
FAST_POWER_RELATIVE_TARGET = 1.0e-6


def transfer_matrix(z_value: float, h_step: float) -> np.ndarray:
    """Return the production midpoint transfer matrix in Python arithmetic."""
    w = math.exp(float(z_value))
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        angle = omega * h_step
        c = math.cos(angle)
        si = math.sin(angle) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h_step
        c = 1.0 + 0.5 * x * x
        si = h_step * (1.0 + x * x / 6.0)
    return np.array(((c - si, -w * si), (w * si, c + si)), dtype=np.float64)


def bound_transfer_perturbation(old_z, new_z, h_values, old_state,
                                new_state):
    """Return an a-posteriori conservative endpoint bound and direct error."""
    old_state = np.asarray(old_state, dtype=np.float64).copy()
    new_state = np.asarray(new_state, dtype=np.float64).copy()
    distance = float(np.linalg.norm(new_state - old_state))
    max_matrix_norm = 0.0
    accumulated = 0.0
    for old_value, new_value, h_step in zip(old_z, new_z, h_values):
        old_matrix = transfer_matrix(old_value, h_step)
        new_matrix = transfer_matrix(new_value, h_step)
        new_norm = float(np.linalg.norm(new_matrix, ord="fro"))
        delta_norm = float(np.linalg.norm(new_matrix - old_matrix, ord="fro"))
        max_matrix_norm = max(max_matrix_norm, new_norm)
        accumulated += delta_norm * float(np.linalg.norm(old_state))
        distance = new_norm * distance + delta_norm * float(np.linalg.norm(old_state))
        old_state = old_matrix @ old_state
        new_state = new_matrix @ new_state
    actual = float(np.linalg.norm(new_state - old_state))
    scale = max(float(np.linalg.norm(new_state)), 1.0e-300)
    relative_bound = distance / scale
    relative_actual = actual / scale
    power_bound = 2.0 * relative_bound + relative_bound * relative_bound
    return {
        "finite": bool(np.isfinite((distance, actual, power_bound)).all()),
        "endpoint_bound": float(distance),
        "endpoint_actual": actual,
        "endpoint_relative_bound": float(relative_bound),
        "endpoint_relative_actual": relative_actual,
        "power_relative_bound": float(power_bound),
        "max_new_matrix_frobenius": max_matrix_norm,
        "summed_local_perturbation": float(accumulated),
    }


def _mode_certificate(first, second, mode: int):
    old_z, new_z, h_values, j0, _ = _path_pair(first, second, mode)
    old_y = math.exp(float(first[6][mode])) * float(first[4][j0])
    new_y = math.exp(float(second[6][mode])) * float(second[4][j0])
    old_initial = np.array((0.0, old_y), dtype=np.float64)
    new_initial = np.array((0.0, new_y), dtype=np.float64)
    row = bound_transfer_perturbation(
        old_z, new_z, h_values, old_initial, new_initial)
    row.update({
        "mode": int(mode),
        "native_intervals": int(old_z.size),
        "substeps": int(h_values.size),
        "z0_delta": float(second[6][mode] - first[6][mode]),
        "path_eligible": True,
        "bound_finite": bool(row["finite"]),
    })
    return row


def _probe(case_name: str):
    first, second = _capture(case_name)
    rows = []
    for mode in range(len(first[5])):
        try:
            row = _mode_certificate(first, second, mode)
            row["bound_eligible"] = bool(
                row["finite"] and
                row["power_relative_bound"] <= FAST_POWER_RELATIVE_TARGET)
        except (ValueError, FloatingPointError, RuntimeError) as exc:
            row = {"mode": int(mode), "path_eligible": False,
                   "bound_finite": False, "bound_eligible": False,
                   "finite": False, "reason": str(exc)}
        rows.append(row)
    bounds = [row["power_relative_bound"] for row in rows
              if "power_relative_bound" in row and math.isfinite(
                  row["power_relative_bound"])]
    return {
        "case": case_name,
        "mode_count": len(rows),
        "event_sequence_eligible": all(
            row.get("path_eligible", False) for row in rows),
        "all_bounds_finite": all(
            row.get("bound_finite", False) for row in rows),
        "fast_target_eligible": all(
            row.get("bound_eligible", False) for row in rows),
        "strict_target_eligible": bool(rows) and all(
            row.get("bound_finite", False)
            and row.get("power_relative_bound", math.inf)
            <= STRICT_POWER_RELATIVE_TARGET for row in rows),
        "max_power_relative_bound": max(bounds) if bounds else math.inf,
        "median_power_relative_bound": statistics.median(bounds) if bounds else math.inf,
        "max_power_relative_actual": max(
            (row["endpoint_relative_actual"] * 2.0
             + row["endpoint_relative_actual"] ** 2 for row in rows
             if "endpoint_relative_actual" in row), default=math.inf),
        "rows": rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append",
                        choices=("default", "lowT", "highT", "stiff", "high_kappa"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    affinity = _pin_current_process(2)
    cases = args.case or ["default", "lowT", "highT", "stiff", "high_kappa"]
    rows = []
    for case_name in cases:
        try:
            rows.append(_probe(case_name))
        except (RuntimeError, ValueError, FloatingPointError) as exc:
            rows.append({"case": case_name, "path_eligible": False,
                         "bound_eligible": False,
                         "reason": str(exc)})
    payload = {
        "schema_version": 1,
        "candidate": "outer_transfer_perturbation_certificate_round67",
        "generated_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "resources": telemetry(workers=1, threads=2),
        "resources_affinity": affinity,
        "production_unchanged": True,
        "targets": {
            "strict_power_relative": STRICT_POWER_RELATIVE_TARGET,
            "fast_screen_power_relative": FAST_POWER_RELATIVE_TARGET,
        },
        "rows": rows,
        "decision": (
            "REJECTED_BY_VACUOUS_BOUND_OR_EVENT_CHANGE"
            if not all(row.get("strict_target_eligible", False) for row in rows)
            else "CONTINUE_TO_FULL_OUTER_GATE"),
        "semantics": (
            "Frobenius-norm variation bound for the production midpoint "
            "transfer; event-sequence changes fail closed. Standalone only."),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
