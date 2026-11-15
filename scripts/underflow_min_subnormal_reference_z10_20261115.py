"""Check Reference tail-depth convergence at the 937-node minimum-subnormal root."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
FAST_Z8_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_tail_depth_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    fast_z8_data = json.loads(FAST_Z8_INPUT.read_text(encoding="utf-8"))
    row8 = next(row for row in data["cases"] if row["frequency_count"] == 937)
    fast8 = next(row for row in fast_z8_data["cases"] if row["frequency_count"] == 937)
    freqs = np.asarray(row8["frequency_log10_hz"], dtype=np.float64)
    root = float(row8["frozen_DN_eff"])
    params = {**data["parameters"], "A_s": 1.0, "r": 1.0, "n_t": 0.0, "DN_eff": 0.0}
    model = LCDM_SG(**params)
    model.gen_expansion()
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            model, freqs, root, z_tail=10.0, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        if np.any(~np.isfinite(base)) or np.any(base <= 0.0):
            raise FloatingPointError("Reference z10 unit-transfer spectrum is not positive finite")
        log_amp = math.log(data["parameters"]["r"]) + math.log(data["parameters"]["A_s"])
        log_power = log_amp + data["parameters"]["n_t"] * (
            math.log(10.0) * freqs - math.log(gp.f_piv)
        )
        log_integrand = np.log(base) + log_power
        shift = float(np.max(log_integrand))
        normalized = np.exp(log_integrand - shift)
        g2_scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    runtime = time.perf_counter() - started
    log_g2 = shift + math.log(g2_scaled)
    h = float(model.derived_param["h"])
    omega_nu = gp.Omega_nh2 / (h * h)
    mapped = gp.Neff0 * math.exp(log_g2) / omega_nu
    ref8_log = np.asarray(row8["reference"]["log10_OmegaGW"], dtype=np.float64)
    fast8_log = np.asarray(fast8["fast_z8"]["log10_OmegaGW"], dtype=np.float64)
    ref10_log = log_integrand / math.log(10.0)
    d_ref = ref10_log - ref8_log
    d_fast = fast8_log - ref10_log
    output = {
        "schema_version": 1,
        "status": "complete_min_subnormal_937_reference_z8_z10_tail_check",
        "scope": "937-node direct Reference at minimum-subnormal sample's frozen 235-root; compare z10 against already computed z8 and Fast z8.",
        "parameters": data["parameters"],
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": row8["frequency_sha256_f64le"],
        "frozen_DN_eff": root,
        "reference": {
            "z_tail": 10.0,
            "rtol": 1e-9,
            "runtime_s": runtime,
            "used_tail_count": int(np.count_nonzero(used_tail)),
            "mapped_DN_gw": float(mapped),
            "map_relative_residual": float(mapped / root - 1.0),
            "runtime_warnings": [str(item.message) for item in caught],
            "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(normalized == 0.0)),
            "quadrature_error_relative": None if qerr is None else float(qerr / g2_scaled),
            "interpolation_error_relative": None if ierr is None else float(ierr / g2_scaled),
            "log10_OmegaGW": ref10_log.tolist(),
        },
        "tail_comparison": {
            "reference_z8_to_z10_max_abs_dex": float(np.max(np.abs(d_ref))),
            "reference_z8_to_z10_p95_abs_dex": float(np.quantile(np.abs(d_ref), 0.95)),
            "fast_z8_to_reference_z10_max_abs_dex": float(np.max(np.abs(d_fast))),
            "fast_z8_to_reference_z10_p95_abs_dex": float(np.quantile(np.abs(d_fast), 0.95)),
            "fast_z8_vs_reference_z10_max_error_index": int(np.argmax(np.abs(d_fast))),
            "fast_z8_vs_reference_z10_max_error_log10_hz": float(freqs[np.argmax(np.abs(d_fast))]),
        },
        "fast_z8": fast8["fast_z8"],
        "reference_z8": row8["reference"],
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "reference.py": "stiffgwpy_fast/reference.py",
                "input_z8_json": "docs/research/results/underflow_min_subnormal_reference_grid_20261115.json",
                "fast_z8_json": "docs/research/results/underflow_min_subnormal_tail_depth_20261115.json",
                "diagnostic_script": "scripts/underflow_min_subnormal_reference_z10_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "tail_comparison": output["tail_comparison"],
                      "reference_z10_map_residual": output["reference"]["map_relative_residual"],
                      "warnings": len(output["reference"]["runtime_warnings"]),
                      "runtime_s": runtime}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
