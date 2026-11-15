"""Compare Fast z_tail hand-offs on one fixed underflow background and grid."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import underflow_grid_convergence_20261115 as grid_run  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_z5_z8_compare_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    baseline = json.loads(INPUT.read_text(encoding="utf-8"))
    frozen_dn = float(baseline["background"]["common_DN_eff"])
    original_tail = grid_run.Z_TAIL_FAST
    grids = {}
    try:
        grid_run.Z_TAIL_FAST = 8.0
        for count in (235, 469, 937):
            row = baseline["fast_fixed_background"][str(count)]
            reference = baseline["reference_fixed_background"][str(count)]
            freqs = np.asarray(
                baseline["frequency_grids"][str(count)]["log10_hz"],
                dtype=np.float64,
            )
            expected_sha = baseline["frequency_grids"][str(count)][
                "sha256_f64le"
            ]
            if grid_run._grid_sha256(freqs) != expected_sha:
                raise AssertionError(f"{count}-node frequency grid hash mismatch")
            started = time.perf_counter()
            fast_z8 = grid_run._fast_run(freqs, frozen_dn=frozen_dn)
            elapsed = time.perf_counter() - started

            fast_z5_log = np.asarray(row["log10_OmegaGW"], dtype=np.float64)
            fast_z8_log = np.asarray(fast_z8["log10_OmegaGW"], dtype=np.float64)
            ref_log = np.asarray(reference["log10_OmegaGW"], dtype=np.float64)
            tail_delta = fast_z8_log - fast_z5_log
            errors_z5 = np.abs(fast_z5_log - ref_log)
            errors_z8 = np.abs(fast_z8_log - ref_log)
            grids[str(count)] = {
                "frequency_count": int(freqs.size),
                "frequency_sha256_f64le": expected_sha,
                "z_tail_5": {
                    "runtime_s": row["runtime_s"],
                    "map_DN_gw": row["mapped_DN_gw"],
                    "spectrum_max_abs_error_dex": float(np.max(errors_z5)),
                    "spectrum_p95_abs_error_dex": float(np.quantile(errors_z5, 0.95)),
                    "log10_OmegaGW": fast_z5_log.tolist(),
                },
                "z_tail_8": {
                    "runtime_s": elapsed,
                    "map_DN_gw": fast_z8["mapped_DN_gw"],
                    "map_residual_at_frozen_background": fast_z8[
                        "map_residual_at_frozen_background"
                    ],
                    "spectrum_max_abs_error_dex": float(np.max(errors_z8)),
                    "spectrum_p95_abs_error_dex": float(np.quantile(errors_z8, 0.95)),
                    "log10_OmegaGW": fast_z8_log.tolist(),
                },
                "z8_minus_z5": {
                    "spectrum_max_abs_delta_dex": float(np.max(np.abs(tail_delta))),
                    "spectrum_p95_abs_delta_dex": float(np.quantile(np.abs(tail_delta), 0.95)),
                    "log10_OmegaGW_delta": tail_delta.tolist(),
                },
            }
    finally:
        grid_run.Z_TAIL_FAST = original_tail
    output = {
        "schema_version": 1,
        "status": "complete_nested_grid_fixed_background_tail_comparison",
        "parameters": baseline["parameters"],
        "background": {
            "DN_eff": frozen_dn,
            "Fast_tail_z5": 5.0,
            "Fast_tail_z8": 8.0,
            "Reference_tail": 8.0,
            "Reference_rtol": baseline["background"]["reference_rtol"],
        },
        "grids": grids,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": baseline["environment"]["threads"],
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "grid_runner.py": _sha(ROOT / "scripts" / "underflow_grid_convergence_20261115.py"),
            "input_json": _sha(INPUT),
            "tail_compare_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "grid_summaries": {
            count: {
                "z5_max_p95_dex": [
                    row["z_tail_5"]["spectrum_max_abs_error_dex"],
                    row["z_tail_5"]["spectrum_p95_abs_error_dex"],
                ],
                "z8_max_p95_dex": [
                    row["z_tail_8"]["spectrum_max_abs_error_dex"],
                    row["z_tail_8"]["spectrum_p95_abs_error_dex"],
                ],
                "z8_map_DN_gw": row["z_tail_8"]["map_DN_gw"],
                "z8_runtime_s": row["z_tail_8"]["runtime_s"],
            }
            for count, row in grids.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
