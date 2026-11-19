"""Same-grid Reference check for an old-Fast resource-timeout boundary point."""

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
OUTPUT = Path(__file__).with_name("validate_h2_identity_resource_boundary_z8.json")
PARAMS = {
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 4000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    config = fast_sgwb.resolve_config("fast", threads=2)
    start = time.perf_counter()
    result = fast_sgwb.SGWB_iter_fast(
        model, config=config, tol=1e-6, physical_guard=True
    )
    fast_seconds = time.perf_counter() - start
    if result is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    total_dn_eff = float(model.cosmo_param["DN_eff"])
    start = time.perf_counter()
    reference = run_reference(
        model, dn_eff=total_dn_eff, freq_subset=freqs,
        z_tail=8.0, rtol=1e-9, self_consistent=False, workers=2,
    )
    reference_seconds = time.perf_counter() - start
    ref_log = np.asarray(reference["log10OmegaGW"], dtype=np.float64)
    ref_linear = np.asarray(reference["Ogw"]) - np.asarray(reference["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
    error = np.abs(fast_log[valid] - ref_log[valid])
    artifact = {
        "status": "complete_same_grid_reference_check",
        "baseline_fast": {
            "base_sha": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
            "classification": "resource_timeout",
            "timeout_seconds": 15.0,
        },
        "params": PARAMS,
        "fast_config": {
            "h": float(config.h), "z_tail": float(config.z_tail),
            "phase_max": float(config.phase_max), "threads": 2,
            "frequency_grid": "goal",
        },
        "fast": {
            "root": float(model.DN_gw[-1]), "total_DN_eff": total_dn_eff,
            "seconds": fast_seconds,
            "h2_identity_used": bool(getattr(model, "_fast_h2_identity_used", False)),
            "nonadiabatic_handoff_selected": int(
                getattr(model, "_fast_nonadiabatic_handoff_selected", 0)
            ),
        },
        "reference_fixed_root": {
            "DN_gw_map": float(reference["DN_gw"]),
            "fixed_point_residual": float(reference["DN_gw"] - total_dn_eff),
            "root_relative_error_proxy": abs(float(reference["DN_gw"]) - float(model.DN_gw[-1]))
            / max(abs(float(reference["DN_gw"])), 1e-300),
            "seconds": reference_seconds, "z_tail": 8.0, "rtol": 1e-9,
        },
        "spectrum_error_dex": {
            "valid_count": int(valid.sum()),
            "max_abs": float(np.max(error)),
            "p95_abs": float(np.quantile(error, 0.95)),
        },
        "raw": {
            "freqs_log10_hz": freqs.tolist(),
            "fast_log10_OmegaGW": fast_log.tolist(),
            "reference_log10_OmegaGW": ref_log.tolist(),
            "reference_Ogw": np.asarray(reference["Ogw"]).tolist(),
            "reference_Oj": np.asarray(reference["Oj"]).tolist(),
            "valid_positive_reference_mask": valid.tolist(),
        },
        "source_sha256": {
            name: _sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fast_root": artifact["fast"]["root"],
        "reference_map": artifact["reference_fixed_root"]["DN_gw_map"],
        "root_error_proxy": artifact["reference_fixed_root"]["root_relative_error_proxy"],
        "spectrum_max_abs_dex": artifact["spectrum_error_dex"]["max_abs"],
        "spectrum_p95_abs_dex": artifact["spectrum_error_dex"]["p95_abs"],
        "reference_seconds": reference_seconds,
        "artifact": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
