"""One z10 convergence check at the cr=1 corner that missed the root gate."""

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
OUTPUT = Path(__file__).with_name("validate_h2_identity_cr1_corner_z10.json")
PARAMS = {
    "Omega_bh2": 1.0e6,
    "Omega_ch2": 1.0e6,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    model = LCDM_SG(**PARAMS)
    result = fast_sgwb.SGWB_iter_fast(
        model, config=fast_sgwb.resolve_config("fast", threads=2),
        tol=1e-6, physical_guard=True,
    )
    if result is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    total = float(model.cosmo_param["DN_eff"])
    start = time.perf_counter()
    reference = run_reference(
        model, dn_eff=total, freq_subset=freqs, z_tail=10.0,
        rtol=1e-9, self_consistent=False, workers=2,
    )
    reference_seconds = time.perf_counter() - start
    ref_log = np.asarray(reference["log10OmegaGW"], dtype=np.float64)
    ref_linear = np.asarray(reference["Ogw"]) - np.asarray(reference["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
    errors = np.abs(fast_log[valid] - ref_log[valid])
    ref_map = float(reference["DN_gw"])
    payload = {
        "status": "complete",
        "reference_code_modified": False,
        "z_tail": 10.0,
        "same_frequency_grid_as_fast": True,
        "params": PARAMS,
        "consistency_tilt": float(model.derived_param["nt"]),
        "derived_N_inf": float(model.derived_param["N_inf"]),
        "fast_root": float(model.DN_gw[-1]),
        "reference_map_at_fast_root": ref_map,
        "root_relative_error_proxy": abs(ref_map - float(model.DN_gw[-1]))
        / max(abs(ref_map), 1e-300),
        "reference_seconds": reference_seconds,
        "nfreq": int(freqs.size),
        "spectrum_valid_count": int(valid.sum()),
        "spectrum_max_abs_dex": float(np.max(errors)),
        "spectrum_p95_abs_dex": float(np.quantile(errors, 0.95)),
        "spectrum_errors_dex": errors.tolist(),
        "source_sha256": {
            name: _sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in payload.items() if k not in ("params", "source_sha256", "spectrum_errors_dex")}))


if __name__ == "__main__":
    main()
