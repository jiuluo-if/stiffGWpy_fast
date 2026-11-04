"""Independently validate order-16 phase modes against Cartesian and Prüfer.

中文：对原残差门槛通过的单阶相位结果执行双参考复核；失败项保留，不替换数值。
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

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PROBE = ROOT / "docs" / "sample2_phase_root_residual_order_probe_20261105.json"
ROOT_GRID = ROOT / "docs" / "sample2_phase_reference_root_res8_tail5_gamma0_20261105.json"
OUTPUT = ROOT / "docs" / "sample2_phase_order16_dual_reference_20261105.json"
COLLOCATION_TOL = 1e-10
DENSE_TOL = 1e-9
SPECTRUM_TOL_DEX = 1e-3
REFERENCE_AGREEMENT_TOL_DEX = 1e-4


def classify_candidate(candidate, cartesian, prufer):
    """Return explicit pass/failure reasons for one order and two references."""
    reasons = []
    if candidate.get("status") != "completed":
        reasons.append("candidate_not_completed")
    if candidate.get("max_collocation_residual_scaled") is None or candidate.get(
        "max_collocation_residual_scaled", math.inf
    ) > COLLOCATION_TOL:
        reasons.append("collocation_residual_gate")
    if candidate.get("max_dense_residual_scaled") is None or candidate.get(
        "max_dense_residual_scaled", math.inf
    ) > DENSE_TOL:
        reasons.append("dense_residual_gate")

    references = (cartesian, prufer)
    if any(item is None for item in references):
        reasons.append("reference_missing")
    else:
        try:
            candidate_log = float(candidate["log10_omega"])
            cartesian_log = float(cartesian["log10_omega"])
            prufer_log = float(prufer["log10_omega"])
            reference_difference = abs(cartesian_log - prufer_log)
            cartesian_error = abs(candidate_log - cartesian_log)
            prufer_error = abs(candidate_log - prufer_log)
            if not all(math.isfinite(value) for value in (candidate_log, cartesian_log, prufer_log)):
                reasons.append("nonfinite_log_spectrum")
            if reference_difference > REFERENCE_AGREEMENT_TOL_DEX:
                reasons.append("references_disagree")
            if cartesian_error > SPECTRUM_TOL_DEX or prufer_error > SPECTRUM_TOL_DEX:
                reasons.append("candidate_reference_error")
            if bool(candidate["tail_used"]) != bool(cartesian["used_tail"]):
                reasons.append("cartesian_tail_semantics_mismatch")
            if bool(candidate["tail_used"]) != bool(prufer["used_tail"]):
                reasons.append("prufer_tail_semantics_mismatch")
        except (KeyError, TypeError, ValueError, OverflowError):
            reasons.append("malformed_spectrum_record")
            reference_difference = cartesian_error = prufer_error = None

    return {
        "accepted": not reasons,
        "reasons": reasons,
        "cartesian_prufer_difference_dex": reference_difference if not reasons else None,
        "candidate_cartesian_error_dex": cartesian_error if not reasons else None,
        "candidate_prufer_error_dex": prufer_error if not reasons else None,
    }


def _reference_value(model, freq, dn_eff, z_tail, z_start, rtol, backend):
    started = time.perf_counter()
    result = R.solve_reference_mode(
        model,
        freq,
        dn_eff,
        z_tail=z_tail,
        z_start=z_start,
        rtol=rtol,
        backend=backend,
    )
    power = float(result["Ogw_today"] - result["Oj_today"])
    if not math.isfinite(power) or power <= 0.0:
        raise FloatingPointError("reference returned non-positive or non-finite OmegaGW")
    return {
        "log10_omega": math.log10(power),
        "used_tail": bool(result["used_tail"]),
        "n_steps": int(result["n_steps"]),
        "seconds": time.perf_counter() - started,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    probe = json.loads(PROBE.read_text(encoding="utf-8"))
    root_grid = json.loads(ROOT_GRID.read_text(encoding="utf-8"))
    params = dict(probe["params"])
    dn_eff = float(params.pop("DN_eff"))
    model = LCDM_SG(**params, DN_eff=dn_eff)
    reference_path = ROOT / "stiffgwpy_fast" / "reference.py"
    probe_path = ROOT / "scripts" / "analyze_sample2_root_phase_order_failures.py"
    phase_path = ROOT / "scripts" / "prototype_nonosc_phase_res4_mode.py"
    selected = []
    for row in probe["rows"]:
        candidate = next(
            (item for item in row["orders"] if item["order"] == 16 and item["status"] == "completed"),
            None,
        )
        if candidate is not None:
            selected.append((row, candidate))

    started = time.perf_counter()
    results = []
    for row, candidate in selected:
        item = {
            "index": int(row["index"]),
            "freq_log10_hz": float(row["freq_log10_hz"]),
            "phase_start": float(row["phase_start"]),
            "candidate": candidate,
        }
        for name, rtol, backend in (
            ("cartesian_tight", 1e-11, "cartesian"),
            ("prufer", 1e-9, "prufer"),
        ):
            try:
                item[name] = _reference_value(
                    model,
                    item["freq_log10_hz"],
                    dn_eff,
                    float(probe["z_tail"]),
                    float(root_grid["z_start"]),
                    rtol,
                    backend,
                )
            except Exception as exc:
                item[name] = {
                    "status": "reference_failed",
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                }
        cartesian = item["cartesian_tight"] if "log10_omega" in item["cartesian_tight"] else None
        prufer = item["prufer"] if "log10_omega" in item["prufer"] else None
        item["validation"] = classify_candidate(candidate, cartesian, prufer)
        results.append(item)

    accepted = [row for row in results if row["validation"]["accepted"]]
    payload = {
        "schema_version": 1,
        "experiment": "sample2_residual_order16_dual_reference_validation",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "probe_script_sha256": hashlib.sha256(probe_path.read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(phase_path.read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        "probe_artifact_sha256": hashlib.sha256(PROBE.read_bytes()).hexdigest(),
        "root_grid_artifact_sha256": hashlib.sha256(ROOT_GRID.read_bytes()).hexdigest(),
        "params": probe["params"],
        "z_tail": probe["z_tail"],
        "z_start": root_grid["z_start"],
        "candidate_order": 16,
        "candidate_residual_gates": {
            "collocation_scaled_max": COLLOCATION_TOL,
            "dense_scaled_max": DENSE_TOL,
        },
        "acceptance_gates": {
            "each_candidate_reference_error_dex_max": SPECTRUM_TOL_DEX,
            "cartesian_prufer_agreement_dex_max": REFERENCE_AGREEMENT_TOL_DEX,
            "tail_semantics_must_match": True,
            "requires_candidate_residual_gates": True,
        },
        "candidate_count": len(selected),
        "validated_count": len(accepted),
        "failed_count": len(results) - len(accepted),
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "Validation is limited to order-16 residual-failure modes at this fixed sample2 root background.",
            "This does not complete the 800-frequency phase map, integral, or self-consistency root.",
        ],
        "rows": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
