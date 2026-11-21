"""z8 Reference probe at the upper-r, low-temperature CR1 edge."""
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
SCREEN = json.loads((HERE / "screen_cr1_r_upperband_lattice.json").read_text(encoding="utf-8"))
PARAMS = next(row["input"]["params"] for row in SCREEN["cases"]
              if row["input"]["case_id"] == "rectangle_02")
OUTPUT = HERE / "validate_cr1_r_upperband_upper_probe_z8.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    config = fast.resolve_config("fast", threads=2)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter()-started
    signal = np.asarray(result["Ogw"]) - np.asarray(result["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(signal) & (signal > 0.0)
    direct_log = np.full(signal.shape, np.nan)
    direct_log[valid] = np.log10(signal[valid])
    ref_map = float(result["DN_gw"])
    row = {
        "case_id":"low_T_high_r_corner", "params":PARAMS,
        "derived_nt":float(model.derived_param["nt"]),
        "derived_N_inf":float(model.derived_param["N_inf"]),
        "fast_root":float(model.DN_gw[-1]), "reference_map_at_fast_root":ref_map,
        "root_relative_error_proxy":abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
        "reference_seconds":elapsed,"spectrum_valid_count":int(valid.sum()),
        "spectrum_max_abs_dex_from_raw_linear":float(np.max(np.abs(fast_log[valid]-direct_log[valid]))),
        "raw":{
            "freqs_log10_hz":freqs.tolist(),
            "fast_log10_OmegaGW":fast_log.tolist(),
            "reference_log10_OmegaGW_flooring_api":np.asarray(result["log10OmegaGW"]).tolist(),
            "reference_raw_Ogw":np.asarray(result["Ogw"]).tolist(),
            "reference_raw_Oj":np.asarray(result["Oj"]).tolist(),
            "reference_direct_log10_Ogw_minus_Oj":direct_log.tolist(),
            "valid_positive_reference_mask":valid.tolist(),
        },
    }
    payload = {
        "status":"complete_cr1_r_upperband_upper_probe",
        "reference_code_modified":False,
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{
            name:hashlib.sha256((ROOT/"stiffgwpy_fast"/name).read_bytes()).hexdigest()
            for name in ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")
        },
        "case":row,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:row[k] for k in (
        "derived_nt","fast_root","reference_map_at_fast_root",
        "root_relative_error_proxy","spectrum_max_abs_dex_from_raw_linear")},indent=2))


if __name__ == "__main__":
    main()
