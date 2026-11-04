"""Screen split Chebyshev segments on the remaining sample2 residual failures.

中文：在最小二乘兜底后，对剩余残差模式试验2段和4段传播；不改变验收门槛。
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import analyze_sample2_root_phase_order_failures as analyzer  # noqa: E402
from scripts import prototype_nonosc_phase_res4_mode as phase  # noqa: E402
from scripts._resource_budget import telemetry  # noqa: E402

ROOT_GRID = ROOT / "docs" / "sample2_phase_reference_root_res8_tail5_gamma0_20261105.json"
RESIDUAL_PROBE = ROOT / "docs" / "sample2_phase_root_residual_order_probe_20261105.json"
DUAL_REFERENCE = ROOT / "docs" / "sample2_phase_residual_orders_dual_reference_20261105.json"
START_OFFSETS = ROOT / "docs" / "sample2_phase_start_offset_probe_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_phase_residual_segment_split_probe_20261105.json"
PARTITIONS = (2, 4)
ORDERS = (16, 24, 32)


def _split_segment(original, parts):
    def run(a, b, regime, order, state, z_left):
        current = state.copy()
        z = z_left
        rows = []
        edges = [a + (b - a) * index / parts for index in range(parts + 1)]
        for left, right in zip(edges[:-1], edges[1:]):
            current, z, row = original(left, right, regime, order, current, z)
            rows.append(row)
        return current, z, {
            "N": [a, b],
            "z": [z_left, z],
            "regime": regime,
            "order": order,
            "collocation_residual_scaled": max(row["collocation_residual_scaled"] for row in rows),
            "dense_residual_scaled": max(row["dense_residual_scaled"] for row in rows),
            "subsegments": rows,
        }

    return run


def main():
    root_grid = json.loads(ROOT_GRID.read_text(encoding="utf-8"))
    residual_probe = json.loads(RESIDUAL_PROBE.read_text(encoding="utf-8"))
    dual_reference = json.loads(DUAL_REFERENCE.read_text(encoding="utf-8"))
    start_probe = json.loads(START_OFFSETS.read_text(encoding="utf-8"))
    phase_sha = hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest()
    if root_grid["phase_script_sha256"] != phase_sha or residual_probe["phase_script_sha256"] != phase_sha:
        raise RuntimeError("phase artifacts and current prototype hash differ")

    params = dict(root_grid["params"])
    dn_eff = float(params.pop("DN_eff"))
    phase.PARAMS.update(params)
    phase.ZTAIL = float(root_grid["z_tail"])
    phase.PHASE_GAMMA_MIN = 0.0
    accepted = {
        int(row["index"]) for row in dual_reference["rows"] if row["validated"]
    } | {
        int(row["index"]) for row in start_probe["rows"] if row["validated"]
    }
    selected = [
        row
        for row in root_grid["rows"]
        if row["status"] != "completed"
        and row["attempts"][-1]["status"] == "residual_rejected"
        and int(row["index"]) not in accepted
    ]
    rows = []
    started = time.perf_counter()
    original = phase.phase_segment
    for root_row in selected:
        phase_start = float(root_row["attempts"][-1]["phase_start"])
        phase.configure_case(dn_eff, root_row["freq_log10_hz"], phase_start)
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
            "index": int(root_row["index"]),
            "freq_log10_hz": float(root_row["freq_log10_hz"]),
            "phase_start": phase_start,
            "reference_log10_omega": float(root_row["reference_log10_omega"]),
            "reference_used_tail": bool(root_row["reference_used_tail"]),
            "partitions": [],
        }
        if not low.success or not low.t_events[0].size:
            item.update(status="handoff_not_reached", message=low.message)
            rows.append(item)
            continue
        item["low_rhs_nfev"] = int(low.nfev)
        for parts in PARTITIONS:
            phase.phase_segment = _split_segment(original, parts)
            trials = [
                analyzer._run_order(
                    root_row["freq_log10_hz"],
                    phase_start,
                    low.y_events[0][0].tolist(),
                    float(low.t_events[0][0]),
                    order,
                    item["reference_log10_omega"],
                )
                for order in ORDERS
            ]
            accepted_orders = [
                trial
                for trial in trials
                if trial["status"] == "completed"
                and trial["abs_dex_error_vs_cartesian"] <= analyzer.SPECTRUM_TOL_DEX
                and bool(trial["tail_used"]) == item["reference_used_tail"]
            ]
            item["partitions"].append(
                {
                    "segments_per_original_branch_interval": parts,
                    "completed_order_count": len(accepted_orders),
                    "orders": trials,
                }
            )
        phase.phase_segment = original
        item["status"] = "probe_complete"
        rows.append(item)
    phase.phase_segment = original
    any_screened = {
        row["index"]
        for row in rows
        if any(part["completed_order_count"] for part in row.get("partitions", []))
    }
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    payload = {
        "schema_version": 1,
        "experiment": "sample2_remaining_residual_modes_segment_split_screen",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "analyzer_script_sha256": hashlib.sha256(Path(analyzer.__file__).read_bytes()).hexdigest(),
        "phase_script_sha256": phase_sha,
        "reference_module_sha256": hashlib.sha256(reference_source.read_bytes()).hexdigest(),
        "root_grid_artifact_sha256": hashlib.sha256(ROOT_GRID.read_bytes()).hexdigest(),
        "residual_probe_artifact_sha256": hashlib.sha256(RESIDUAL_PROBE.read_bytes()).hexdigest(),
        "dual_reference_artifact_sha256": hashlib.sha256(DUAL_REFERENCE.read_bytes()).hexdigest(),
        "phase_start_offset_artifact_sha256": hashlib.sha256(START_OFFSETS.read_bytes()).hexdigest(),
        "params": root_grid["params"],
        "z_tail": root_grid["z_tail"],
        "partitions": list(PARTITIONS),
        "orders": list(ORDERS),
        "acceptance_gates": {
            "collocation_scaled_max": 1e-10,
            "dense_scaled_max": 1e-9,
            "max_cartesian_spectrum_error_dex": analyzer.SPECTRUM_TOL_DEX,
            "tail_semantics_must_match": True,
        },
        "mode_count_probed": len(rows),
        "mode_count_with_any_screened_order": len(any_screened),
        "screened_mode_indices": sorted(any_screened),
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "This combines interval splitting with the current least-squares Newton fallback on the remaining residual failures.",
            "The screened Cartesian comparison is not a fresh dual-reference validation.",
            "This does not complete the phase integral map or root.",
        ],
        "rows": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
