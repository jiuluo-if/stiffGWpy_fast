"""Replay full-grid Riccati failures independently at each collocation order.

中文：独立检查 gamma=0 全网格残差失败点是否存在较高阶配点可通过的情形。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import prototype_nonosc_phase_res4_mode as phase  # noqa: E402
from scripts._resource_budget import telemetry  # noqa: E402

FULLGRID = ROOT / "docs" / "nonosc_phase_res4_fullgrid_gamma0_20261104.json"
REFERENCE = ROOT / "docs" / "reference_prufer_cartesian_paired_fullgrid_20261104.json"
OUTPUT = ROOT / "docs" / "nonosc_phase_res4_residual_order_probe_20261104.json"


def probe_mode(row, reference_log, params):
    freq = float(row["freq_log10_hz"])
    phase_start = float(row["attempts"][-1]["phase_start"])
    phase.PARAMS.update({key: value for key, value in params.items() if key != "DN_eff"})
    phase.configure_case(params["DN_eff"], freq, phase_start)
    t0 = time.perf_counter()
    low = phase.solve_ivp(
        phase.low_rhs,
        (phase.NSTART, phase.bg.N_inf),
        (phase.Z0, 1.0, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=[1e-12, 1e-13, 1e-13],
        events=[phase.low_event],
        dense_output=False,
    )
    item = {
        "index": int(row["index"]),
        "freq_log10_hz": freq,
        "phase_start": phase_start,
        "attempt_count_from_fullgrid": len(row["attempts"]),
        "orders": [],
    }
    if not low.success or not low.t_events[0].size:
        item.update(status="handoff_not_reached", message=low.message)
        return item

    n_phase = float(low.t_events[0][0])
    low_state = low.y_events[0][0].tolist()
    item["low_rhs_nfev"] = int(low.nfev)
    item["low_seconds"] = time.perf_counter() - t0
    for order in (16, 24, 32):
        order_started = time.perf_counter()
        try:
            state, bounds, segments, _, tail_used, zend = phase.candidate(
                order, low_state, n_phase, phase_start
            )
            end_n = bounds[-1][0]
            spectrum = (
                phase.spectrum_from_tail(state, end_n, zend)
                if tail_used
                else phase.spectrum_today(state, zend)
            )
            value = float(spectrum["OmegaGW"])
            if not math.isfinite(value) or value <= 0.0:
                raise FloatingPointError("phase spectrum is not finite and positive")
            log_value = math.log10(value)
            item["orders"].append(
                {
                    "order": order,
                    "status": "completed",
                    "tail_used": tail_used,
                    "log10_omega": log_value,
                    "abs_dex_error_vs_cartesian": abs(log_value - reference_log),
                    "max_collocation_residual_scaled": max(
                        segment["collocation_residual_scaled"] for segment in segments
                    ),
                    "max_dense_residual_scaled": max(
                        segment["dense_residual_scaled"] for segment in segments
                    ),
                    "seconds": time.perf_counter() - order_started,
                }
            )
        except Exception as exc:  # failures are retained with their actual reason
            item["orders"].append(
                {
                    "order": order,
                    "status": "failed",
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                    "seconds": time.perf_counter() - order_started,
                }
            )
    item["status"] = "order_probe_complete"
    return item


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    fullgrid = json.loads(FULLGRID.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    phase_source = Path(phase.__file__)
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    if hashlib.sha256(reference_source.read_bytes()).hexdigest() != reference["reference_module_sha256"]:
        raise RuntimeError("Cartesian archive and current Reference module hashes differ")
    if fullgrid["reference_archive_sha256"] != hashlib.sha256(REFERENCE.read_bytes()).hexdigest():
        raise RuntimeError("gamma-zero screen and Cartesian archive digests differ")
    if fullgrid["phase_script_sha256"] != hashlib.sha256(phase_source.read_bytes()).hexdigest():
        raise RuntimeError("gamma-zero screen was produced by a different phase script")

    rows_by_index = {int(row["index"]): row for row in reference["rows"]}
    rows = [
        row
        for row in fullgrid["rows"]
        if row["status"] != "completed" and row["attempts"][-1]["status"] == "residual_rejected"
    ]
    phase.PHASE_GAMMA_MIN = float(fullgrid["solver_settings"]["minimum_phase_interval_gamma"])
    params = dict(fullgrid["parameters"])
    started = time.perf_counter()
    results = [
        probe_mode(row, rows_by_index[int(row["index"])]["cartesian_log10_omega"], params)
        for row in rows
    ]
    order_counts = {
        str(order): sum(
            any(item["order"] == order and item["status"] == "completed" for item in row["orders"])
            for row in results
        )
        for order in (16, 24, 32)
    }
    payload = {
        "schema_version": 1,
        "experiment": "independent_collocation_order_probe_on_gamma_zero_residual_failures",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(phase_source.read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(reference_source.read_bytes()).hexdigest(),
        "fullgrid_artifact_sha256": hashlib.sha256(FULLGRID.read_bytes()).hexdigest(),
        "reference_archive_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        "parameters": params,
        "minimum_gamma": phase.PHASE_GAMMA_MIN,
        "mode_count": len(results),
        "order_success_counts": order_counts,
        "order32_successes": order_counts["32"],
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "rows": results,
        "limitations": [
            "Each order is replayed independently at the final phase start recorded by the full-grid screen.",
            "This probes residual-gate failures only; it is not a complete parameter-domain result.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
