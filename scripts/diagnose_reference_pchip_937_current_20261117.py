"""Diagnose current-source 937-node Reference PCHIP underflow and stable arithmetic."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.diagnose_reference_pchip_multigrid_20261116 import diagnose_case  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_current_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "reference_pchip_937_current_20261117.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source.get("status") != "complete_fixed_background_direct_nested_grid_screen":
        raise RuntimeError("Current-source direct nested-grid screen is incomplete")
    for rel, expected in source["source_sha256"].items():
        if sha256(ROOT / rel) != expected:
            raise RuntimeError(f"Input source hash differs for {rel}")

    count = "937"
    frequencies = np.asarray(source["frequency_grids"][count]["log10_hz"], dtype=np.float64)
    spectrum = np.asarray(
        source["reference_fixed_background"][count]["log10_OmegaGW"],
        dtype=np.float64,
    )
    if frequencies.size != 937 or spectrum.size != 937:
        raise RuntimeError("Expected all 937 directly solved Reference nodes")
    diagnostic = diagnose_case(
        "current_source_reference_937",
        frequencies,
        spectrum,
        parameters=source["parameters"],
        frozen_dn=float(source["background"]["common_DN_eff"]),
        z_tail=float(source["background"]["reference_tail"]),
        rtol=float(source["background"]["reference_rtol"]),
        source_sha256={
            **source["source_sha256"],
            "input_artifact": sha256(INPUT),
            "diagnostic_runner": sha256(Path(__file__)),
            "stable_log_pchip": sha256(ROOT / "scripts" / "stable_log_pchip.py"),
            "diagnostic_library": sha256(ROOT / "scripts" / "diagnose_reference_pchip_multigrid_20261116.py"),
        },
    )
    output = {
        "schema_version": 1,
        "status": "complete_current_source_reference_pchip_937_diagnostic",
        "scope": "Use the exact 937 Reference nodes directly solved on the current-source nested grid; compare SciPy linear PCHIP arithmetic against log arithmetic for the same positive linear PCHIP polynomial. Reference production code is not modified.",
        "input_artifact": INPUT.relative_to(ROOT).as_posix(),
        "source_sha256": {
            **source["source_sha256"],
            "input_artifact": sha256(INPUT),
            "diagnostic_runner": sha256(Path(__file__)),
            "stable_log_pchip": sha256(ROOT / "scripts" / "stable_log_pchip.py"),
            "diagnostic_library": sha256(ROOT / "scripts" / "diagnose_reference_pchip_multigrid_20261116.py"),
        },
        "diagnostic": diagnostic,
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "normalized_zero_nodes": diagnostic["normalized_zero_nodes"],
        "normalized_positive_subnormal_nodes": diagnostic["normalized_positive_subnormal_nodes"],
        "original_warnings": diagnostic["original_linear_pchip"]["warnings"],
        "relative_integral_difference": diagnostic["stable_log_arithmetic_same_linear_pchip"]["relative_integral_difference"],
        "midpoint_max_abs_dex": diagnostic["midpoint_interpolation"]["finite_linear_vs_stable_max_abs_dex"],
        "midpoint_p95_abs_dex": diagnostic["midpoint_interpolation"]["finite_linear_vs_stable_p95_abs_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
