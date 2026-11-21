"""Reference bisection probe at 13.9990234375M for the negative-tilt edge."""
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
    "Omega_bh2": 13.9990234375e6, "Omega_ch2": 0.1201075, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 0.01,
    "n_t": -0.0302734375, "cr": 0.0, "T_re": 2000.0,
    "DN_re": 10.0, "kappa10": 0.01,
}
OUTPUT = HERE / "probe_cr0_density_13p9990234375m_nt0302734375_reference.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    solved = SGWB_iter_fast(model, config=resolve_config("fast", threads=2), tol=1.0e-6, physical_guard=True)
    row = {"params": PARAMS, "fast_success": solved is model, "fast_failure_reason": model.fast_failure_reason}
    if solved is model:
        started = time.perf_counter()
        try:
            ref = run_reference(
                model, dn_eff=float(model.cosmo_param["DN_eff"]),
                freq_subset=np.asarray(model.f, dtype=np.float64), z_tail=8.0,
                rtol=1.0e-9, self_consistent=False, workers=2,
            )
            mapped = float(ref["DN_gw"])
            root = float(model.DN_gw[-1])
            fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
            ref_log = np.asarray(ref["log10OmegaGW"], dtype=np.float64)
            valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (np.asarray(ref["Ogw"]) - np.asarray(ref["Oj"]) > 0.0)
            errors = np.abs(fast_log[valid] - ref_log[valid])
            row.update({
                "reference_status": "success", "reference_seconds": time.perf_counter()-started,
                "fast_root": root, "reference_mapped_root": mapped,
                "root_relative_error": abs(mapped-root)/max(abs(mapped), 1.0e-300),
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": float(np.max(errors)),
            })
        except Exception as exc:
            row.update({
                "reference_status": "integration_failure", "reference_seconds": time.perf_counter()-started,
                "exception_type": type(exc).__name__, "message": str(exc),
            })
    payload = {
        "status": "complete", "reference_settings": {"z_tail": 8.0, "rtol": 1.0e-9, "self_consistent": False, "workers": 2},
        "reference_code_modified": False,
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "reference.py", "exact_background.py", "_wkb_transfer.py")
        },
        "case": row,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
