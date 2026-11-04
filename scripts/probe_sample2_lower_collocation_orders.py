"""Probe lower Riccati collocation orders on sample2 residual failures.

中文：对88个既有残差失败频率运行8/12/16阶，保留失败并检验相邻阶收敛。
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
OUTPUT = ROOT / "docs" / "sample2_phase_root_lower_order_probe_20261105.json"
ORDERS = (8, 12, 16)


def main():
    root_grid = json.loads(ROOT_GRID.read_text(encoding="utf-8"))
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    analyzer_source = ROOT / "scripts" / "analyze_sample2_root_phase_order_failures.py"
    phase.PARAMS.update({key: value for key, value in root_grid["params"].items() if key != "DN_eff"})
    phase.ZTAIL = float(root_grid["z_tail"])
    phase.PHASE_GAMMA_MIN = 0.0
    rows = [
        row
        for row in root_grid["rows"]
        if row["status"] != "completed" and row["attempts"][-1]["status"] == "residual_rejected"
    ]

    results = []
    started = time.perf_counter()
    for row in rows:
        phase_start = float(row["attempts"][-1]["phase_start"])
        phase.configure_case(root_grid["params"]["DN_eff"], row["freq_log10_hz"], phase_start)
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
            "freq_log10_hz": float(row["freq_log10_hz"]),
            "phase_start": phase_start,
            "reference_log10_omega": float(row["reference_log10_omega"]),
            "reference_used_tail": bool(row["reference_used_tail"]),
            "orders": [],
        }
        if not low.success or not low.t_events[0].size:
            item.update(status="handoff_not_reached", message=low.message)
            results.append(item)
            continue
        n_phase = float(low.t_events[0][0])
        low_state = low.y_events[0][0].tolist()
        item["orders"] = [
            analyzer._run_order(
                row["freq_log10_hz"],
                phase_start,
                low_state,
                n_phase,
                order,
                item["reference_log10_omega"],
            )
            for order in ORDERS
        ]
        pairs = []
        for first, second in ((0, 1), (1, 2)):
            a, b = item["orders"][first], item["orders"][second]
            if a["status"] != "completed" or b["status"] != "completed":
                continue
            difference = abs(a["log10_omega"] - b["log10_omega"])
            reference_error = max(a["abs_dex_error_vs_cartesian"], b["abs_dex_error_vs_cartesian"])
            tail_match = (
                bool(a["tail_used"]) == item["reference_used_tail"]
                and bool(b["tail_used"]) == item["reference_used_tail"]
            )
            pairs.append(
                {
                    "orders": [a["order"], b["order"]],
                    "abs_dex_difference": difference,
                    "max_abs_dex_error_vs_cartesian": reference_error,
                    "tail_semantics_match": tail_match,
                    "accepted": (
                        difference <= analyzer.PAIR_TOL_DEX
                        and reference_error <= analyzer.SPECTRUM_TOL_DEX
                        and tail_match
                    ),
                }
            )
        item["adjacent_order_pairs"] = pairs
        item["status"] = "probe_complete"
        results.append(item)

    accepted = [
        (row, pair)
        for row in results
        for pair in row.get("adjacent_order_pairs", [])
        if pair["accepted"]
    ]
    payload = {
        "schema_version": 1,
        "experiment": "sample2_root_residual_failure_lower_order_probe",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "analyzer_script_sha256": hashlib.sha256(analyzer_source.read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(reference_source.read_bytes()).hexdigest(),
        "root_grid_artifact_sha256": hashlib.sha256(ROOT_GRID.read_bytes()).hexdigest(),
        "params": root_grid["params"],
        "z_tail": root_grid["z_tail"],
        "minimum_gamma": 0.0,
        "orders": list(ORDERS),
        "pair_acceptance": {
            "max_adjacent_order_difference_dex": analyzer.PAIR_TOL_DEX,
            "max_independent_cartesian_error_dex": analyzer.SPECTRUM_TOL_DEX,
            "requires_each_order_collocation_and_dense_residual_gates": True,
        },
        "mode_count_probed": len(results),
        "order_success_counts": {
            str(order): sum(
                any(x["order"] == order and x["status"] == "completed" for x in row["orders"])
                for row in results
            )
            for order in ORDERS
        },
        "accepted_pair_count": len(accepted),
        "accepted_pair_mode_indices": sorted({row["index"] for row, _ in accepted}),
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "This is a lower-order residual-failure screen at one fixed sample2 root background.",
            "Adjacent-order checks still require independent-reference validation before being counted as expanded coverage.",
            "This does not form a complete phase integral map or root.",
        ],
        "rows": results,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
