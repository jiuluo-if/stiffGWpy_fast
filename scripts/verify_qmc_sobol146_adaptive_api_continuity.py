"""Verify that exposing adaptive options preserved the validated Fast map."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FAST_PATH = ROOT / "docs" / "qmc_sobol146_fast_adaptive_grid_20261105.json"
REFERENCE_PATH = ROOT / "docs" / "qmc_sobol146_fast_adaptive_reference_20261105.json"
OUTPUT = ROOT / "docs" / "qmc_sobol146_adaptive_api_continuity_20261105.json"
CASE = "adaptive_sigma_exact_tail10_1e-06"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> None:
    fast_record = json.loads(FAST_PATH.read_text(encoding="utf-8"))
    ref_record = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    fast = fast_record["results"][CASE]
    if fast["frequency_sha256"] != ref_record["frequency_sha256"]:
        raise AssertionError("high-level Fast wrapper changed the validated frequency array")
    dn_delta = abs(float(fast["DN_gw"]) - float(ref_record["fast_DN_gw"]))
    if dn_delta != 0.0:
        raise AssertionError("high-level Fast wrapper changed the validated DN_gw")
    if fast_record["fast_module_sha256"] != ref_record["fast_module_sha256"]:
        raise AssertionError("Fast numerical kernel source changed after the Reference check")
    if not fast["local_error_budget"]["uncertified_terms"]:
        raise AssertionError("expected the one-point result to retain its scoped calibration status")

    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_adaptive_high_level_api_continuity",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "fast_grid_artifact": FAST_PATH.name,
        "fast_grid_artifact_sha256": _sha(FAST_PATH),
        "reference_validation_artifact": REFERENCE_PATH.name,
        "reference_validation_artifact_sha256": _sha(REFERENCE_PATH),
        "reference_validation_original_fast_snapshot_sha256": ref_record["fast_artifact_sha256"],
        "fast_case": CASE,
        "checks": {
            "frequency_sha256_equal": True,
            "frequency_count_equal": (
                int(fast["frequency_count"]) == int(ref_record["frequency_count"])
            ),
            "DN_gw_bitwise_equal": dn_delta == 0.0,
            "fast_core_sha256_equal": True,
            "absolute_DN_gw_delta": dn_delta,
            "fast_core_sha256": fast_record["fast_module_sha256"],
            "adaptive_high_level_api_sha256": fast_record["model_api_module_sha256"],
        },
        "runtime_configuration": {
            "freq_grid": fast["freq_grid"],
            "freq_grid_target": fast["frequency_grid_target"],
            "freq_grid_max_points": fast["frequency_grid_max_points"],
            "sigma_exact": fast["sigma_exact"],
            "z_tail": fast["z_tail"],
            "frequency_count": fast["frequency_count"],
            "DN_gw": fast["DN_gw"],
        },
        "script_sha256": _sha(Path(__file__)),
        "resources": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "This links the original independent Reference comparison to the public high-level invocation through exact core-source, frequency-grid, and DN_gw equality checks.",
            "It does not rerun the 321-mode Reference or reconstruct the prior Fast spectrum array; the low-level Fast core and all numerical settings are unchanged.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"checks": record["checks"], "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()
