"""Check the Reference tail hand-off convergence on the 235-node holdout."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import underflow_grid_convergence_20261115 as grid_run  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_reference_tail_z8_z10_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    freqs = np.asarray(data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    common_dn = float(data["background"]["common_DN_eff"])
    log_amplitude = grid_run._log_tensor_amplitude(
        grid_run.PARAMS,
        {"A_t": 0.0},
    )
    reference_model = grid_run._reference_unit_model()
    original_tail = grid_run.Z_TAIL_REFERENCE
    try:
        grid_run.Z_TAIL_REFERENCE = 10.0
        ref_z10 = grid_run._reference_map(
            reference_model, freqs, common_dn, log_amplitude
        )
    finally:
        grid_run.Z_TAIL_REFERENCE = original_tail
    ref_z8 = data["reference_fixed_background"]["235"]
    log8 = np.asarray(ref_z8["log10_OmegaGW"], dtype=np.float64)
    log10 = np.asarray(ref_z10["log10_OmegaGW"], dtype=np.float64)
    delta = log10 - log8
    output = {
        "schema_version": 1,
        "status": "complete_reference_tail_convergence_check",
        "parameters": data["parameters"],
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": data["frequency_grids"]["235"]["sha256_f64le"],
        "frozen_DN_eff": common_dn,
        "reference_rtol": data["background"]["reference_rtol"],
        "z_tail_8": {
            "runtime_s": ref_z8["runtime_s"],
            "mapped_DN_gw": ref_z8["mapped_DN_gw"],
            "log10_OmegaGW": log8.tolist(),
            "warnings": ref_z8["runtime_warnings"],
        },
        "z_tail_10": {
            "runtime_s": ref_z10["runtime_s"],
            "mapped_DN_gw": ref_z10["mapped_DN_gw"],
            "log10_OmegaGW": log10.tolist(),
            "warnings": ref_z10["runtime_warnings"],
        },
        "z10_minus_z8": {
            "relative_map_delta": float(ref_z10["mapped_DN_gw"] / ref_z8["mapped_DN_gw"] - 1.0),
            "max_abs_spectrum_delta_dex": float(np.max(np.abs(delta))),
            "p95_abs_spectrum_delta_dex": float(np.quantile(np.abs(delta), 0.95)),
            "spectrum_delta_dex": delta.tolist(),
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "exact_background.py": _sha(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "grid_runner.py": _sha(ROOT / "scripts" / "underflow_grid_convergence_20261115.py"),
            "input_json": _sha(INPUT),
            "tail_check_script": _sha(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "z8_map": ref_z8["mapped_DN_gw"],
        "z10_map": ref_z10["mapped_DN_gw"],
        "map_relative_delta": output["z10_minus_z8"]["relative_map_delta"],
        "spectrum_max_p95_dex": [
            output["z10_minus_z8"]["max_abs_spectrum_delta_dex"],
            output["z10_minus_z8"]["p95_abs_spectrum_delta_dex"],
        ],
        "warnings": ref_z10["runtime_warnings"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
