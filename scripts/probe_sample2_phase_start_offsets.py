"""Refine phase handoff starts for sample2 modes still failing residual gates.

中文：对未验证的剩余频率细化内部相位起点，并用双参考筛选通过的候选。
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import analyze_sample2_root_phase_order_failures as analyzer  # noqa: E402
from scripts import prototype_nonosc_phase_res4_mode as phase  # noqa: E402
from scripts import validate_sample2_residual_orders as validation  # noqa: E402
from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

ROOT_GRID = ROOT / "docs" / "sample2_phase_reference_root_res8_tail5_gamma0_20261105.json"
RESIDUAL_PROBE = ROOT / "docs" / "sample2_phase_root_residual_order_probe_20261105.json"
DUAL_REFERENCE = ROOT / "docs" / "sample2_phase_residual_orders_dual_reference_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_phase_start_offset_probe_20261105.json"
START_OFFSETS = tuple(index / 128.0 for index in range(-64, 64))
ORDERS = (16, 24, 32)


def select_best_candidate(trials, reference_used_tail):
    """Select the lowest-residual candidate that already passes screening gates."""
    eligible = []
    for trial in trials:
        candidate = trial["candidate"]
        if candidate.get("status") != "completed":
            continue
        try:
            archive_error = float(candidate["abs_dex_error_vs_cartesian"])
            candidate_log = float(candidate["log10_omega"])
        except (KeyError, TypeError, ValueError, OverflowError):
            continue
        if not math.isfinite(archive_error) or archive_error > analyzer.SPECTRUM_TOL_DEX:
            continue
        if not math.isfinite(candidate_log):
            continue
        if bool(candidate.get("tail_used")) != bool(reference_used_tail):
            continue
        if not validation._residual_within_gate(
            candidate, "max_collocation_residual_scaled", validation.COLLOCATION_TOL
        ):
            continue
        if not validation._residual_within_gate(
            candidate, "max_dense_residual_scaled", validation.DENSE_TOL
        ):
            continue
        eligible.append(
            {
                "phase_start": float(trial["phase_start"]),
                "candidate": candidate,
                "residual_merit": max(
                    candidate["max_collocation_residual_scaled"],
                    candidate["max_dense_residual_scaled"],
                ),
            }
        )
    return min(eligible, key=lambda row: row["residual_merit"]) if eligible else None


def _run_low(freq, phase_start, dn_eff):
    phase.configure_case(dn_eff, freq, phase_start)
    started = time.perf_counter()
    result = phase.solve_ivp(
        phase.low_rhs,
        (phase.NSTART, phase.bg.N_inf),
        (phase.Z0, 1.0, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=[1e-12, 1e-13, 1e-13],
        events=[phase.low_event],
        dense_output=False,
    )
    return result, time.perf_counter() - started


def main(argv=None):
    del argv
    root_grid = json.loads(ROOT_GRID.read_text(encoding="utf-8"))
    residual_probe = json.loads(RESIDUAL_PROBE.read_text(encoding="utf-8"))
    dual_reference = json.loads(DUAL_REFERENCE.read_text(encoding="utf-8"))
    if root_grid["phase_script_sha256"] != hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest():
        raise RuntimeError("root-grid artifact and current phase script hashes differ")
    if residual_probe["phase_script_sha256"] != hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest():
        raise RuntimeError("residual probe and current phase script hashes differ")

    params = dict(root_grid["params"])
    dn_eff = float(params.pop("DN_eff"))
    phase.PARAMS.update(params)
    phase.ZTAIL = float(root_grid["z_tail"])
    phase.PHASE_GAMMA_MIN = 0.0
    model = LCDM_SG(**params, DN_eff=dn_eff)
    accepted_before = {
        int(row["index"]) for row in dual_reference["rows"] if row["validated"]
    }
    probe_by_index = {int(row["index"]): row for row in residual_probe["rows"]}
    selected = [
        row
        for row in root_grid["rows"]
        if row["status"] != "completed"
        and row["attempts"][-1]["status"] == "residual_rejected"
        and int(row["index"]) not in accepted_before
    ]

    started = time.perf_counter()
    results = []
    for row in selected:
        index = int(row["index"])
        base_start = float(probe_by_index[index]["phase_start"])
        item = {
            "index": index,
            "freq_log10_hz": float(row["freq_log10_hz"]),
            "base_phase_start": base_start,
            "reference_log10_omega": float(row["reference_log10_omega"]),
            "reference_used_tail": bool(row["reference_used_tail"]),
            "phase_start_trials": [],
        }
        candidates = []
        for offset in START_OFFSETS:
            phase_start = base_start + offset
            if not 0.0 <= phase_start <= 5.0:
                item["phase_start_trials"].append(
                    {"phase_start": phase_start, "status": "outside_scan_interval"}
                )
                continue
            try:
                low, low_seconds = _run_low(row["freq_log10_hz"], phase_start, dn_eff)
            except Exception as exc:
                item["phase_start_trials"].append(
                    {
                        "phase_start": phase_start,
                        "status": "low_solver_exception",
                        "exception_type": type(exc).__name__,
                        "message": str(exc),
                    }
                )
                continue
            trial = {
                "phase_start": phase_start,
                "low_seconds": low_seconds,
                "low_rhs_nfev": int(low.nfev),
                "orders": [],
            }
            if not low.success or not low.t_events[0].size:
                trial.update(status="handoff_not_reached", message=low.message)
                item["phase_start_trials"].append(trial)
                continue
            trial["status"] = "handoff_reached"
            for order in ORDERS:
                candidate = analyzer._run_order(
                    row["freq_log10_hz"],
                    phase_start,
                    low.y_events[0][0].tolist(),
                    float(low.t_events[0][0]),
                    order,
                    item["reference_log10_omega"],
                )
                trial["orders"].append(candidate)
                candidates.append({"phase_start": phase_start, "candidate": candidate})
            item["phase_start_trials"].append(trial)

        best = select_best_candidate(candidates, item["reference_used_tail"])
        item["selected_candidate"] = best
        if best is None:
            item["validated"] = False
            item["validation"] = {"accepted": False, "reasons": ["no_screened_candidate"]}
            results.append(item)
            continue

        reference_started = time.perf_counter()
        try:
            cartesian = validation._reference_value(
                model,
                item["freq_log10_hz"],
                dn_eff,
                float(root_grid["z_tail"]),
                float(root_grid["z_start"]),
                1e-11,
                "cartesian",
            )
            prufer = validation._reference_value(
                model,
                item["freq_log10_hz"],
                dn_eff,
                float(root_grid["z_tail"]),
                float(root_grid["z_start"]),
                1e-9,
                "prufer",
            )
            item["cartesian_tight"] = cartesian
            item["prufer"] = prufer
            item["validation"] = validation.classify_candidate(best["candidate"], cartesian, prufer)
        except Exception as exc:
            item["validation"] = {
                "accepted": False,
                "reasons": ["reference_exception"],
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
        item["reference_seconds"] = time.perf_counter() - reference_started
        item["validated"] = bool(item["validation"]["accepted"])
        results.append(item)

    accepted = [row for row in results if row["validated"]]
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    payload = {
        "schema_version": 1,
        "experiment": "sample2_residual_failure_phase_start_offset_dual_reference",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest(),
        "analyzer_script_sha256": hashlib.sha256(Path(analyzer.__file__).read_bytes()).hexdigest(),
        "validation_script_sha256": hashlib.sha256(Path(validation.__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(reference_source.read_bytes()).hexdigest(),
        "root_grid_artifact_sha256": hashlib.sha256(ROOT_GRID.read_bytes()).hexdigest(),
        "residual_probe_artifact_sha256": hashlib.sha256(RESIDUAL_PROBE.read_bytes()).hexdigest(),
        "dual_reference_artifact_sha256": hashlib.sha256(DUAL_REFERENCE.read_bytes()).hexdigest(),
        "params": root_grid["params"],
        "z_tail": root_grid["z_tail"],
        "z_start": root_grid["z_start"],
        "phase_start_offsets": list(START_OFFSETS),
        "orders": list(ORDERS),
        "candidate_residual_gates": {
            "collocation_scaled_max": validation.COLLOCATION_TOL,
            "dense_scaled_max": validation.DENSE_TOL,
        },
        "independent_reference_gates": {
            "candidate_cartesian_prufer_error_dex_max": validation.SPECTRUM_TOL_DEX,
            "cartesian_prufer_agreement_dex_max": validation.REFERENCE_AGREEMENT_TOL_DEX,
            "tail_semantics_must_match": True,
        },
        "mode_count_scanned": len(results),
        "mode_count_with_screened_candidate": sum(row["selected_candidate"] is not None for row in results),
        "validated_mode_count": len(accepted),
        "validated_mode_indices": sorted(row["index"] for row in accepted),
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "This refines only phase-start handoff on previously residual-rejected frequencies at one fixed sample2 root background.",
            "The 800-bin phase map remains incomplete and no integral map or new root is claimed.",
            "Fast production solver inputs and outputs are not changed by this prototype scan.",
        ],
        "rows": results,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
