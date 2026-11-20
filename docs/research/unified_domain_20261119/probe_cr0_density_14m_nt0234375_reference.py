"""Probe Reference feasibility at the midpoint toward the negative-tilt edge."""
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
OUTPUT = HERE / "probe_cr0_density_14m_nt0234375_reference.json"
BASE_PARAMS = {
    "Omega_ch2": 0.1201075, "H0": 67.32117, "DN_eff": 0.0,
    "A_s": 2.100549e-9, "r": 0.01, "cr": 0.0,
    "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}
CASES = (
    ("lower_density_edge", 13.875e6),
    ("interior", 13.9375e6),
    ("upper_density_edge", 14.0e6),
)
TILT = -0.0234375


def main() -> None:
    rows = []
    for case_id, density in CASES:
        params = dict(BASE_PARAMS, Omega_bh2=density, n_t=TILT)
        model = LCDM_SG(**params)
        solved = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1.0e-6, physical_guard=True)
        row = {"case_id": case_id, "params": params, "fast_success": solved is model,
               "fast_failure_reason": model.fast_failure_reason}
        if solved is model:
            started = time.perf_counter()
            try:
                ref = run_reference(
                    model, dn_eff=float(model.cosmo_param["DN_eff"]),
                    freq_subset=np.asarray(model.f, dtype=np.float64), z_tail=8.0,
                    rtol=1.0e-9, self_consistent=False, workers=2,
                )
                elapsed = time.perf_counter() - started
                mapped = float(ref["DN_gw"])
                fast_root = float(model.DN_gw[-1])
                fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
                ref_log = np.asarray(ref["log10OmegaGW"], dtype=np.float64)
                valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (np.asarray(ref["Ogw"]) - np.asarray(ref["Oj"]) > 0.0)
                errors = np.abs(fast_log[valid] - ref_log[valid])
                row.update({
                    "reference_status": "success", "reference_seconds": elapsed,
                    "fast_root": fast_root, "reference_mapped_root": mapped,
                    "root_relative_error": abs(mapped-fast_root)/max(abs(mapped), 1.0e-300),
                    "spectrum_valid_count": int(valid.sum()),
                    "spectrum_max_abs_dex": float(np.max(errors)),
                })
            except Exception as exc:
                row.update({
                    "reference_status": "integration_failure",
                    "reference_seconds": time.perf_counter() - started,
                    "exception_type": type(exc).__name__, "message": str(exc),
                })
        rows.append(row)
        print(case_id, row.get("reference_status", "fast_failure"), flush=True)
    payload = {
        "status": "complete", "tilt": TILT,
        "reference_settings": {"z_tail": 8.0, "rtol": 1.0e-9, "self_consistent": False, "workers": 2},
        "reference_code_modified": False,
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py", "_wkb_transfer.py")
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
