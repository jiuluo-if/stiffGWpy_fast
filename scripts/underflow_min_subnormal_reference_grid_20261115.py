"""Direct 469/937 Reference comparison at the minimum-subnormal amplitude root."""

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

R = 5e-324
COUNTS = (469, 937)
BASE_RESULT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
GRID_RESULT = ROOT / "docs" / "research" / "results" / "underflow_amplitude_grid_convergence_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"


def nested_midpoints(freqs_desc: np.ndarray) -> np.ndarray:
    ascending = np.asarray(freqs_desc, dtype=np.float64)[::-1]
    refined = np.empty(2 * ascending.size - 1, dtype=np.float64)
    refined[::2] = ascending
    refined[1::2] = 0.5 * (ascending[:-1] + ascending[1:])
    return refined[::-1]


def main() -> int:
    base = json.loads(BASE_RESULT.read_text(encoding="utf-8"))
    convergence = json.loads(GRID_RESULT.read_text(encoding="utf-8"))
    grids = {235: np.asarray(base["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)}
    grids[469] = nested_midpoints(grids[235])
    grids[937] = nested_midpoints(grids[469])
    sample = next(case for case in convergence["cases"] if case["r"] == R)
    root235 = float(sample["root_235"])
    cases = []
    for count in COUNTS:
        freqs = grids[count]
        fast = ladder.run_case(R, freqs, frozen_dn=root235)
        reference = ladder.reference_fixed_background(freqs, root235, amplitudes=(R,))
        ref = reference["cases"][0]
        fast_log = np.asarray(fast["log10_OmegaGW"], dtype=np.float64)
        ref_log = np.asarray(ref["log10_OmegaGW"], dtype=np.float64)
        signed = fast_log - ref_log
        absolute = np.abs(signed)
        cases.append({
            "frequency_count": int(freqs.size),
            "frequency_log10_hz": freqs.tolist(),
            "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
            "frozen_DN_eff": root235,
            "fast": {
                "mapped_DN_gw": fast["mapped_DN_gw"],
                "root_map_residual_relative": fast["map_residual_vs_evaluation_background_relative"],
                "runtime_s": fast["runtime_s"],
                "log10_OmegaGW": fast_log.tolist(),
            },
            "reference": {
                "mapped_DN_gw": ref["mapped_DN_gw"],
                "map_relative_residual": ref["mapped_DN_gw"] / root235 - 1.0,
                "runtime_s": reference["runtime_s"],
                "z_tail": reference["z_tail"],
                "rtol": reference["rtol"],
                "used_tail_count": reference["used_tail_count"],
                "runtime_warnings": reference["runtime_warnings"],
                "normalized_nodes_underflowed_to_zero": ref["normalized_nodes_underflowed_to_zero"],
                "quadrature_error_relative": ref["quadrature_error_relative"],
                "interpolation_error_relative": ref["interpolation_error_relative"],
                "log10_OmegaGW": ref["log10_OmegaGW"],
            },
            "spectrum_error": {
                "max_abs_dex": float(np.max(absolute)),
                "p95_abs_dex": float(np.quantile(absolute, 0.95)),
                "mean_abs_dex": float(np.mean(absolute)),
                "max_error_index": int(np.argmax(absolute)),
                "max_error_log10_hz": float(freqs[np.argmax(absolute)]),
                "max_error_signed_dex": float(signed[np.argmax(absolute)]),
                "signed_error_dex": signed.tolist(),
            },
        })

    output = {
        "schema_version": 1,
        "status": "complete_min_subnormal_direct_reference_grid_screen",
        "scope": "Direct Fast and Cartesian Reference at 469/937 nodes, frozen to the minimum-subnormal sample's 235-node self-consistent root.",
        "parameters": {**ladder.BASE_PARAMS, "r": R},
        "frequency_grids": {str(row["frequency_count"]): {
            "count": row["frequency_count"],
            "sha256_f64le": row["frequency_sha256_f64le"],
        } for row in cases},
        "reference_tail_and_tolerance": {"z_tail": 8.0, "rtol": 1e-9},
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
                "diagnostic_script": "scripts/underflow_min_subnormal_reference_grid_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "count": row["frequency_count"],
            "fast_map_residual": row["fast"]["root_map_residual_relative"],
            "reference_map_residual": row["reference"]["map_relative_residual"],
            "spectrum_error_max_p95_dex": [row["spectrum_error"]["max_abs_dex"], row["spectrum_error"]["p95_abs_dex"]],
            "reference_warnings": len(row["reference"]["runtime_warnings"]),
            "reference_runtime_s": row["reference"]["runtime_s"],
        } for row in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
