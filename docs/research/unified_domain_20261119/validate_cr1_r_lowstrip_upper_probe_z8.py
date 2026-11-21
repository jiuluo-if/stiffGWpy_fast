"""z8 Reference probe at the candidate lower-r rectangle's tight upper edge."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
fast = __import__("stiffgwpy_fast.fast_sgwb", fromlist=["resolve_config"])
reference = __import__("stiffgwpy_fast.reference", fromlist=["run_reference"]).run_reference
LCDM_SG = __import__("stiffgwpy_fast.stiff_SGWB", fromlist=["LCDM_SG"]).LCDM_SG
SCREEN = json.loads((HERE / "screen_cr1_r_lowstrip_lattice.json").read_text(encoding="utf-8"))
PARAMS = next(row["input"]["params"] for row in SCREEN["cases"]
              if row["input"]["case_id"] == "rectangle_02")
OUTPUT = HERE / "validate_cr1_r_lowstrip_upper_probe_z8.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    config = fast.resolve_config("fast", threads=2)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    frequencies = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=frequencies, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter() - started
    ref_log = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    valid = (np.isfinite(fast_log) & np.isfinite(ref_log)
             & (np.asarray(result["Ogw"]) - np.asarray(result["Oj"]) > 0.0))
    ref_map = float(result["DN_gw"])
    row = {
        "case_id": "low_T_upper_r_corner", "params": PARAMS,
        "derived_nt": float(model.derived_param["nt"]),
        "derived_N_inf": float(model.derived_param["N_inf"]),
        "fast_root": float(model.DN_gw[-1]), "reference_map_at_fast_root": ref_map,
        "root_relative_error_proxy": abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
        "reference_seconds": elapsed, "spectrum_valid_count": int(valid.sum()),
        "spectrum_max_abs_dex": float(np.max(np.abs(fast_log[valid]-ref_log[valid]))),
        "raw": {
            "freqs_log10_hz": frequencies.tolist(),
            "fast_log10_OmegaGW": fast_log.tolist(),
            "reference_log10_OmegaGW": ref_log.tolist(),
            "reference_Ogw": np.asarray(result["Ogw"]).tolist(),
            "reference_Oj": np.asarray(result["Oj"]).tolist(),
            "valid_positive_reference_mask": valid.tolist(),
        },
    }
    payload = {
        "status": "complete_cr1_r_lowstrip_upper_probe",
        "reference_code_modified": False,
        "region_probe": {"T_re": 1750.0, "r": 0.005, "kappa10": 0.01},
        "fast_config": {"threads": 2, "frequency_grid": "goal", "reference_z_tail": 8.0, "reference_rtol": 1e-9},
        "source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")
        },
        "case": row,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: row[k] for k in (
        "derived_nt", "derived_N_inf", "fast_root", "reference_map_at_fast_root",
        "root_relative_error_proxy", "spectrum_max_abs_dex")}, indent=2))


if __name__ == "__main__":
    main()
