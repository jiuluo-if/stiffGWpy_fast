"""Compare Reference's original PCHIP integral with stable log-node arithmetic."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.stable_log_pchip import LogPchip  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_amplitude_ladder_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_reference_pchip_warning_diagnostic_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    freqs = np.asarray(data["frequency_log10_hz"], dtype=np.float64)
    ref = data["reference_fixed_background"]
    h = float(data["parameters_except_r"]["H0"]) / 100.0
    map_scale = gp.Neff0 / (gp.Omega_nh2 / (h * h))
    cases = []
    for row in ref["cases"]:
        r = float(row["r"])
        if r == 0.0:
            continue
        log_y = np.asarray(row["log10_OmegaGW"], dtype=np.float64) * math.log(10.0)
        shift = float(np.max(log_y))
        normalized = np.exp(log_y - shift)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            g2_scaled, qerr, ierr = integrate_spectrum(
                freqs, normalized, np.zeros_like(normalized)
            )
        original_log_g2 = shift + math.log(g2_scaled)

        stable = LogPchip(freqs[::-1], log_y[::-1])
        stable_log_g2 = math.log(math.log(10.0)) + stable.integrate_log()
        original_map = map_scale * math.exp(original_log_g2)
        stable_map = map_scale * math.exp(stable_log_g2)
        cases.append({
            "r": r,
            "frequency_count": int(freqs.size),
            "input_positive_finite": bool(np.isfinite(log_y).all()),
            "linear_normalized_zero_nodes": int(np.count_nonzero(normalized == 0.0)),
            "linear_normalized_nonzero_subnormal_nodes": int(np.count_nonzero(
                (normalized > 0.0) & (normalized < np.finfo(float).tiny)
            )),
            "original_pchip_runtime_warnings": [str(item.message) for item in caught],
            "stable_log_pchip_runtime_warnings": [],
            "original_map_DN_gw": original_map,
            "stable_same_linear_pchip_map_DN_gw": stable_map,
            "stable_vs_original_map_relative_difference": stable_map / original_map - 1.0,
            "reference_map_DN_gw": row["mapped_DN_gw"],
            "original_vs_saved_reference_relative_difference": original_map / float(row["mapped_DN_gw"]) - 1.0,
            "quadrature_error_relative": None if qerr is None else float(qerr / g2_scaled),
            "interpolation_error_relative": None if ierr is None else float(ierr / g2_scaled),
        })

    output = {
        "schema_version": 1,
        "status": "complete_stable_linear_pchip_warning_diagnostic_not_adopted",
        "scope": "Use the direct unit-amplitude Reference samples already solved at the shared frozen root; compare original normalized linear PCHIP quadrature to positive PCHIP integrated from log nodes.",
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "frozen_DN_eff": data["reference_fixed_background"]["frozen_DN_eff"],
        "reference_z_tail": data["reference_fixed_background"]["z_tail"],
        "reference_rtol": data["reference_fixed_background"]["rtol"],
        "cases": cases,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "stable_log_pchip.py": sha(ROOT / "scripts" / "stable_log_pchip.py"),
            "input_json": sha(INPUT),
            "diagnostic_script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "r": row["r"],
            "zero_nodes": row["linear_normalized_zero_nodes"],
            "warnings": len(row["original_pchip_runtime_warnings"]),
            "stable_map_relative_difference": row["stable_vs_original_map_relative_difference"],
            "original_vs_saved_reference_relative_difference": row["original_vs_saved_reference_relative_difference"],
        } for row in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
