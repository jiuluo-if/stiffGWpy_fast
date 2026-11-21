"""Current-source Reference probes for a temperature-limited high-r band."""
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
PARAMS_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in SCREEN["cases"]}
CASES = [("T1875_high_r", PARAMS_BY_ID["rectangle_05"]),
         ("T2000_high_r", PARAMS_BY_ID["rectangle_08"])]
OUTPUT = HERE / "validate_cr1_r_upperband_edge_probes_z8.json"


def main() -> None:
    config = fast.resolve_config("fast", threads=2)
    rows = []
    for case_id, params in CASES:
        model = LCDM_SG(**params)
        if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
            raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
        freqs = np.asarray(model.f, dtype=np.float64)
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
        started = time.perf_counter()
        result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                           freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                           self_consistent=False, workers=2)
        elapsed = time.perf_counter()-started
        signal = np.asarray(result["Ogw"])-np.asarray(result["Oj"])
        valid = np.isfinite(fast_log) & np.isfinite(signal) & (signal > 0.0)
        direct = np.full(signal.shape, np.nan)
        direct[valid] = np.log10(signal[valid])
        ref_map = float(result["DN_gw"])
        rows.append({
            "case_id":case_id,"params":params,
            "derived_nt":float(model.derived_param["nt"]),
            "derived_N_inf":float(model.derived_param["N_inf"]),
            "fast_root":float(model.DN_gw[-1]),
            "reference_map_at_fast_root":ref_map,
            "root_relative_error_proxy":abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
            "reference_seconds":elapsed,"spectrum_valid_count":int(valid.sum()),
            "spectrum_max_abs_dex_from_raw_linear":float(np.max(np.abs(fast_log[valid]-direct[valid]))),
            "raw":{
                "freqs_log10_hz":freqs.tolist(),"fast_log10_OmegaGW":fast_log.tolist(),
                "reference_log10_OmegaGW_flooring_api":np.asarray(result["log10OmegaGW"]).tolist(),
                "reference_raw_Ogw":np.asarray(result["Ogw"]).tolist(),
                "reference_raw_Oj":np.asarray(result["Oj"]).tolist(),
                "reference_direct_log10_Ogw_minus_Oj":direct.tolist(),
                "valid_positive_reference_mask":valid.tolist(),
            },
        })
        OUTPUT.write_text(json.dumps({"status":"in_progress","cases":rows},indent=2)+"\n",encoding="utf-8")
        print(case_id,rows[-1]["root_relative_error_proxy"],rows[-1]["spectrum_max_abs_dex_from_raw_linear"],flush=True)
    payload={
        "status":"complete_cr1_upperband_temperature_edge_probes",
        "reference_code_modified":False,
        "region_probe":{"T_re":[1875.0,2000.0],"r":0.01,"kappa10":0.01},
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{name:hashlib.sha256((ROOT/"stiffgwpy_fast"/name).read_bytes()).hexdigest() for name in
                         ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")},
        "cases":rows,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
