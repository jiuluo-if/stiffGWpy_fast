"""z8 Reference probe for a mathematically valid very-low-r CR1 point."""
from __future__ import annotations

import hashlib
import json
import math
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
H2_vec = __import__("stiffgwpy_fast.exact_background", fromlist=["H2_vec"]).H2_vec
PARAMS = {
    "Omega_bh2": 7.0e5, "Omega_ch2": 7.0e5, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "r": 1e-40, "n_t": 0.0,
    "cr": 1.0, "T_re": 1875.0, "DN_re": 10.0, "kappa10": 0.01,
}
OUTPUT = HERE / "validate_cr1_r_log_tail_probe_z8.json"


def main() -> None:
    model = LCDM_SG(**PARAMS)
    derived = model.derived_param
    config = fast.resolve_config("fast", threads=2)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    h2 = H2_vec(model.Nv, model, 0.0)
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter()-started
    ref_log = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    ref_signal = np.asarray(result["Ogw"]) - np.asarray(result["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(ref_signal) & (ref_signal > 0.0)
    direct_ref_log = np.full(ref_signal.shape, np.nan)
    direct_ref_log[valid] = np.log10(ref_signal[valid])
    ref_map = float(result["DN_gw"])
    delta_n = float(derived["N_inf"]-derived["N_re"])
    payload = {
        "status":"complete_cr1_low_r_reference_probe",
        "reference_code_modified":False,
        "math_validity":{
            "derived_nt":float(derived["nt"]), "N_inf":float(derived["N_inf"]),
            "N_re":float(derived["N_re"]), "Delta_N":delta_n,
            "derived_r_min_at_fixed_other_inputs":PARAMS["r"]*math.exp(-3*delta_n),
            "H2_positive_finite":bool(np.isfinite(h2).all() and np.min(h2)>0.0),
            "H2_min":float(np.min(h2)),
        },
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{
            name:hashlib.sha256((ROOT/"stiffgwpy_fast"/name).read_bytes()).hexdigest()
            for name in ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")
        },
        "case":{
            "params":PARAMS,"fast_root":float(model.DN_gw[-1]),
            "reference_map_at_fast_root":ref_map,
            "root_relative_error_proxy":abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
            "reference_seconds":elapsed,"spectrum_valid_count":int(valid.sum()),
            "spectrum_max_abs_dex":float(np.max(np.abs(fast_log[valid]-ref_log[valid]))),
            "spectrum_max_abs_dex_from_raw_linear":float(np.max(np.abs(fast_log[valid]-direct_ref_log[valid]))),
            "reference_log_floor_bins":int(np.sum(ref_log == -40.0)),
            "raw":{
                "freqs_log10_hz":freqs.tolist(),
                "fast_log10_OmegaGW":fast_log.tolist(),
                "reference_log10_OmegaGW":ref_log.tolist(),
                "reference_Ogw":np.asarray(result["Ogw"]).tolist(),
                "reference_Oj":np.asarray(result["Oj"]).tolist(),
                "reference_direct_log10_Ogw_minus_Oj":direct_ref_log.tolist(),
                "valid_positive_reference_mask":valid.tolist(),
            },
        },
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({**payload["math_validity"], **{k:payload["case"][k] for k in (
        "fast_root","reference_map_at_fast_root","root_relative_error_proxy",
        "spectrum_max_abs_dex_from_raw_linear","reference_log_floor_bins")}},indent=2))


if __name__ == "__main__":
    main()
