"""Separate frequency-grid mapping error and root shifts across a cr=0 amplitude ladder."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import underflow_amplitude_ladder_20261115 as ladder  # noqa: E402

AMPLITUDES = (1e-316, 1e-318, 5e-324)
COUNTS = (235, 469, 937)
BASE_RESULT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_amplitude_grid_convergence_20261115.json"


def nested_midpoints(freqs_desc: np.ndarray) -> np.ndarray:
    ascending = np.asarray(freqs_desc, dtype=np.float64)[::-1]
    refined = np.empty(2 * ascending.size - 1, dtype=np.float64)
    refined[::2] = ascending
    refined[1::2] = 0.5 * (ascending[:-1] + ascending[1:])
    return refined[::-1]


def grid_hash(freqs: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()


def compact(row: dict[str, object], grid_name: str) -> dict[str, object]:
    return {
        "grid": grid_name,
        "frequency_count": row["frequency_count"],
        "frequency_sha256_f64le": grid_hash(np.asarray(row.get("_freqs"), dtype=np.float64)),
        "solve_mode": row["solve_mode"],
        "r": row["r"],
        "evaluation_background_DN_eff": row["evaluation_background_DN_eff"],
        "root_DN_eff": row["root_DN_eff"],
        "mapped_DN_gw": row["mapped_DN_gw"],
        "root_map_residual_relative": row["root_map_residual_relative"],
        "map_residual_vs_evaluation_background_relative": row["map_residual_vs_evaluation_background_relative"],
        "runtime_s": row["runtime_s"],
        "fast_failure_reason": row["fast_failure_reason"],
    }


def main() -> int:
    base = json.loads(BASE_RESULT.read_text(encoding="utf-8"))
    grids = {235: np.asarray(base["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)}
    grids[469] = nested_midpoints(grids[235])
    grids[937] = nested_midpoints(grids[469])
    if tuple(len(grids[n]) for n in COUNTS) != COUNTS:
        raise AssertionError("nested direct grids have unexpected lengths")

    cases = []
    for r in AMPLITUDES:
        root235_result = ladder.run_case(r, grids[235])
        root235 = float(root235_result["root_DN_eff"])
        per_grid = []
        for count in COUNTS:
            freqs = grids[count]
            root_row = ladder.run_case(r, freqs)
            fixed_row = ladder.run_case(r, freqs, frozen_dn=root235)
            root_row["_freqs"] = freqs
            fixed_row["_freqs"] = freqs
            per_grid.append({
                "count": count,
                "self_consistent": compact(root_row, str(count)),
                "fixed_at_235_root": compact(fixed_row, str(count)),
            })
            root_row.pop("_freqs")
            fixed_row.pop("_freqs")
        reference_root = float(per_grid[0]["self_consistent"]["root_DN_eff"])
        for row in per_grid:
            root = float(row["self_consistent"]["root_DN_eff"])
            row["self_consistent"]["relative_root_change_vs_235"] = root / reference_root - 1.0
            mapped = float(row["fixed_at_235_root"]["mapped_DN_gw"])
            base_mapped = float(per_grid[0]["fixed_at_235_root"]["mapped_DN_gw"])
            row["fixed_at_235_root"]["relative_mapping_change_vs_235"] = (
                mapped / base_mapped - 1.0
            )
        reference = ladder.reference_fixed_background(grids[235], root235, amplitudes=(r,))
        ref_row = next(item for item in reference["cases"] if item["r"] == r)
        fast235 = next(row for row in per_grid if row["count"] == 235)["self_consistent"]
        fast_spectrum = np.asarray(root235_result["log10_OmegaGW"], dtype=np.float64)
        ref_spectrum = np.asarray(ref_row["log10_OmegaGW"], dtype=np.float64)
        reference_screen = {
            "frequency_count": 235,
            "frozen_DN_eff": root235,
            "reference_z_tail": reference["z_tail"],
            "reference_rtol": reference["rtol"],
            "runtime_s": reference["runtime_s"],
            "max_abs_spectrum_error_dex": float(np.max(np.abs(fast_spectrum - ref_spectrum))),
            "map_relative_error": float(fast235["mapped_DN_gw"] / ref_row["mapped_DN_gw"] - 1.0),
            "reference_map_relative_residual": float(ref_row["mapped_DN_gw"] / root235 - 1.0),
            "runtime_warnings": reference["runtime_warnings"],
            "normalized_nodes_underflowed_to_zero": ref_row.get("normalized_nodes_underflowed_to_zero"),
            "quadrature_error_relative": ref_row.get("quadrature_error_relative"),
            "interpolation_error_relative": ref_row.get("interpolation_error_relative"),
        }
        cases.append({"r": r, "root_235": root235, "grids": per_grid, "reference_at_235_root": reference_screen})

    output = {
        "schema_version": 1,
        "status": "complete_direct_fast_amplitude_grid_screen",
        "scope": "Three discrete cr=0 amplitudes, direct Fast on nested 235/469/937 native frequency grids; each grid is both solved self-consistently and mapped at that amplitude's 235-point root.",
        "parameters_except_r": {key: value for key, value in ladder.BASE_PARAMS.items() if key != "r"},
        "amplitudes": list(AMPLITUDES),
        "grids": {
            str(count): {
                "count": int(grids[count].size),
                "sha256_f64le": grid_hash(grids[count]),
            }
            for count in COUNTS
        },
        "cases": cases,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "amplitude_ladder_script": "scripts/underflow_amplitude_ladder_20261115.py",
                "diagnostic_script": "scripts/underflow_amplitude_grid_convergence_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "r": case["r"],
            "roots": [(row["count"], row["self_consistent"]["root_DN_eff"])
                      for row in case["grids"]],
            "fixed_map_changes": [(row["count"], row["fixed_at_235_root"]["relative_mapping_change_vs_235"])
                                  for row in case["grids"]],
            "reference_235": case["reference_at_235_root"],
        } for case in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
