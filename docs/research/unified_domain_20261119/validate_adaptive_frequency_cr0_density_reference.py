"""Independent Reference map on the nested adaptive CR0 density-edge grid."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

HERE = Path(__file__).resolve().parent
PARAMS = {
    "Omega_bh2": 13.99609375e6,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": -0.0302734375,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
OUTPUT = HERE / "validate_adaptive_frequency_cr0_density_reference.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    started = time.perf_counter()
    solved = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1.0e-6,
        physical_guard=True,
        freq_grid_target=1.0e-6,
        freq_grid_max_points=1500,
    )
    fast_seconds = time.perf_counter() - started
    row = {
        "params": PARAMS,
        "fast_success": solved is model,
        "fast_failure_reason": model.fast_failure_reason,
        "fast_seconds": fast_seconds,
    }
    if solved is model:
        row.update(
            {
                "fast_root": float(model.DN_gw[-1]),
                "fast_total_DN_eff": float(model.cosmo_param["DN_eff"]),
                "freq_grid_used": model.freq_grid_used,
                "freq_grid_error": float(model.freq_grid_error),
                "frequency_count": int(model.f.size),
                "frequencies": np.asarray(model.f, dtype=float).tolist(),
            }
        )
        started = time.perf_counter()
        reference = run_reference(
            model,
            dn_eff=float(model.cosmo_param["DN_eff"]),
            freq_subset=np.asarray(model.f, dtype=np.float64),
            z_tail=8.0,
            rtol=1.0e-9,
            self_consistent=False,
            workers=2,
        )
        reference_seconds = time.perf_counter() - started
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
        reference_log = np.asarray(reference["log10OmegaGW"], dtype=np.float64)
        valid = np.isfinite(fast_log) & np.isfinite(reference_log)
        spectrum_error = np.abs(fast_log[valid] - reference_log[valid])
        mapped_root = float(reference["DN_gw"])
        row.update(
            {
                "reference_status": "success",
                "reference_seconds": reference_seconds,
                "reference_settings": {
                    "z_tail": 8.0,
                    "rtol": 1.0e-9,
                    "self_consistent": False,
                    "workers": 2,
                },
                "reference_mapped_root": mapped_root,
                "root_relative_error": abs(mapped_root - row["fast_root"])
                / max(abs(mapped_root), 1.0e-300),
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": (
                    float(np.max(spectrum_error)) if spectrum_error.size else None
                ),
            }
        )

    payload = {
        "status": "complete",
        "reference_code_modified": False,
        "source_sha256": {
            name: hashlib.sha256(
                (ROOT / "stiffgwpy_fast" / name).read_bytes()
            ).hexdigest()
            for name in (
                "fast_sgwb.py",
                "reference.py",
                "exact_background.py",
                "_wkb_transfer.py",
            )
        },
        "case": row,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
