"""Direct same-grid Reference checks for the S2-overflow rectangle."""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
fast_sgwb = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
SGWB_iter_fast = fast_sgwb.SGWB_iter_fast
resolve_config = fast_sgwb.resolve_config
OUTPUT = Path(__file__).with_name("validate_h2_identity_region_z8.json")
CASES = [
    (
        "training_low_density_low_temperature_corner",
        {"Omega_bh2": 1.0e5, "T_re": 1000.0},
    ),
    (
        "training_high_density_high_temperature_corner",
        {"Omega_bh2": 1.0e6, "T_re": 4000.0},
    ),
    (
        "independent_interior_holdout",
        {"Omega_bh2": 5.0e5, "T_re": 2500.0},
    ),
    (
        "previously_successful_region_control",
        {"Omega_bh2": 3.0e5, "T_re": 4000.0},
    ),
]
BASE_PARAMS = {
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    config = resolve_config("fast", threads=2)
    for case_id, overrides in CASES:
        params = dict(BASE_PARAMS)
        params.update(overrides)
        model = LCDM_SG(**params)
        start = time.perf_counter()
        result = SGWB_iter_fast(
            model, config=config, tol=1e-6, physical_guard=True
        )
        fast_seconds = time.perf_counter() - start
        row = {
            "case_id": case_id,
            "params": params,
            "fast_success": result is model,
            "fast_failure_reason": model.fast_failure_reason,
            "fast_seconds": fast_seconds,
        }
        if result is model:
            freqs = np.asarray(model.f, dtype=np.float64)
            fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
            total_dn_eff = float(model.cosmo_param["DN_eff"])
            fast_root = float(model.DN_gw[-1])
            reference_start = time.perf_counter()
            reference = run_reference(
                model,
                dn_eff=total_dn_eff,
                freq_subset=freqs,
                z_tail=8.0,
                rtol=1e-9,
                self_consistent=False,
                workers=2,
            )
            reference_seconds = time.perf_counter() - reference_start
            ref_log = np.asarray(reference["log10OmegaGW"], dtype=np.float64)
            ref_linear = np.asarray(reference["Ogw"]) - np.asarray(reference["Oj"])
            valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
            abs_error = np.abs(fast_log[valid] - ref_log[valid])
            reference_map = float(reference["DN_gw"])
            row.update({
                "fast_root": fast_root,
                "total_DN_eff": total_dn_eff,
                "reference_map_at_fast_root": reference_map,
                "fixed_point_residual": reference_map - total_dn_eff,
                "root_relative_error_proxy": abs(reference_map - fast_root) / max(abs(reference_map), 1e-300),
                "reference_seconds": reference_seconds,
                "reference_z_tail": 8.0,
                "reference_rtol": 1e-9,
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": float(np.max(abs_error)),
                "spectrum_p95_abs_dex": float(np.quantile(abs_error, 0.95)),
                "raw": {
                    "freqs_log10_hz": freqs.tolist(),
                    "fast_log10_OmegaGW": fast_log.tolist(),
                    "reference_log10_OmegaGW": ref_log.tolist(),
                    "reference_Ogw": np.asarray(reference["Ogw"]).tolist(),
                    "reference_Oj": np.asarray(reference["Oj"]).tolist(),
                    "valid_positive_reference_mask": valid.tolist(),
                },
            })
        rows.append(row)
        artifact = {
            "status": "in_progress",
            "reference_code_modified": False,
            "fast_config": {
                "accuracy_mode": "fast",
                "h": float(config.h),
                "z_tail": float(config.z_tail),
                "phase_max": float(config.phase_max),
                "threads": 2,
                "frequency_grid": "goal",
            },
            "source_sha256": {
                name: _sha256(ROOT / "stiffgwpy_fast" / name)
                for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
            },
            "cases": rows,
        }
        OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
        print(case_id, row.get("fast_root"), row.get("spectrum_max_abs_dex"))

    artifact["status"] = "complete_selected_region_reference_checks"
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
