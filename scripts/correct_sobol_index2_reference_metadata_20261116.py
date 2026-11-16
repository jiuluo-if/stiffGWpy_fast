"""Correct only index labels in the completed Sobol-2 Reference artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index2_reference_20261116.json"
PREFLIGHT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_holdout_20261116.json"
CALCULATION_SCRIPT = ROOT / "scripts" / "verify_cr_positive_sobol_index2_reference_calculation_20261116.py"
CORRECTED_RUNNER = ROOT / "scripts" / "verify_cr_positive_sobol_index2_reference_20261116.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    sample = next(row for row in preflight["samples"] if row["sobol_index"] == 2)
    if result["sobol_index"] != 2:
        raise RuntimeError("raw calculation is not the expected Sobol index 2")
    if result["status"] != "complete_sobol_index1_direct_reference_holdout":
        raise RuntimeError("unexpected raw status; refuse to rewrite metadata")
    if result["source_sha256"]["script"] != sha(CALCULATION_SCRIPT):
        raise RuntimeError("calculation source snapshot does not match the recorded source SHA")
    if result["source_sha256"]["preflight_input"] != sha(PREFLIGHT):
        raise RuntimeError("Sobol preflight input SHA changed")
    if result["parameters"] != sample["parameters"]:
        raise RuntimeError("Sobol index 2 parameters do not match the preflight record")
    if result["current_fast"]["frequency_sha256_f64le"] != sample["fast"]["frequency_sha256_f64le"]:
        raise RuntimeError("direct Reference grid does not match the Sobol-2 Fast grid")
    if result["current_fast"]["root_DN_eff"] != sample["fast"]["root_DN_eff"]:
        raise RuntimeError("direct Reference root does not match the Sobol-2 Fast root")
    if result["reference"]["source_sha256"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("Reference source SHA differs from current Reference")

    old_status = result["status"]
    old_scope = result["scope"]
    result["status"] = "complete_sobol_index2_direct_reference_holdout"
    result["scope"] = (
        "Follow-up Sobol index 2 from the fixed 8-point cr>0 cutoff design. "
        "Current Fast root/grid were verified against preflight, then Reference z8/z10/z12 "
        "were directly solved at that fixed root and exact frequency grid. Reference source remains unchanged."
    )
    result["metadata_label_correction"] = {
        "scientific_values_changed": False,
        "corrected_fields": ["status", "scope"],
        "previous_status": old_status,
        "previous_scope": old_scope,
        "reason": "The run used a copied index-1 runner and its status/scope strings were not fully relabeled; Fast/Reference computations, parameters, root and grid were independently checked against Sobol index 2 before this metadata-only correction.",
        "calculation_script_snapshot": str(CALCULATION_SCRIPT.relative_to(ROOT)).replace("\\", "/"),
        "calculation_script_sha256": sha(CALCULATION_SCRIPT),
        "corrected_reusable_runner": str(CORRECTED_RUNNER.relative_to(ROOT)).replace("\\", "/"),
        "corrected_reusable_runner_sha256": sha(CORRECTED_RUNNER),
        "correction_script": str(Path(__file__).relative_to(ROOT)).replace("\\", "/"),
        "correction_script_sha256": sha(Path(__file__)),
    }
    result["source_sha256"]["calculation_script_snapshot"] = sha(CALCULATION_SCRIPT)
    result["source_sha256"]["corrected_runner"] = sha(CORRECTED_RUNNER)
    result["source_sha256"]["metadata_correction_script"] = sha(Path(__file__))
    RESULT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "result": str(RESULT),
        "status": result["status"],
        "sobol_index": result["sobol_index"],
        "scientific_values_changed": result["metadata_label_correction"]["scientific_values_changed"],
        "root_error_vs_reference_z12": result["comparison"]["fast_root_relative_error_vs_reference_z12_map"],
        "spectrum_max_abs_dex": result["comparison"]["fast_z5_vs_reference_z12"]["max_abs_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
