"""Targeted current-source z8 Reference checks near the CR1 low-T edge."""
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
LATTICE = json.loads((HERE / "screen_cr1_nearby_reheat_kappa_lattice.json").read_text(encoding="utf-8"))["cases"]
SOBOL = json.loads((HERE / "screen_cr1_nearby_reheat_kappa_sobol.json").read_text(encoding="utf-8"))["cases"]
CASES = [
    ("low_T_low_kappa", LATTICE[2]["input"]["params"]),
    ("low_T_high_kappa", LATTICE[4]["input"]["params"]),
    ("independent_sobol_holdout", SOBOL[2]["input"]["params"]),
]
OUTPUT = HERE / "validate_cr1_nearby_reheat_kappa_z8.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    config = fast.resolve_config("fast", threads=2)
    for case_id, params in CASES:
        model = LCDM_SG(**params)
        if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
            raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
        freqs = np.asarray(model.f, dtype=np.float64)
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
        started = time.perf_counter()
        ref = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                        freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                        self_consistent=False, workers=2)
        elapsed = time.perf_counter() - started
        ref_log = np.asarray(ref["log10OmegaGW"], dtype=np.float64)
        valid = (np.isfinite(fast_log) & np.isfinite(ref_log)
                 & (np.asarray(ref["Ogw"]) - np.asarray(ref["Oj"]) > 0.0))
        ref_map = float(ref["DN_gw"])
        row = {
            "case_id": case_id, "params": params,
            "derived_nt": float(model.derived_param["nt"]),
            "derived_N_inf": float(model.derived_param["N_inf"]),
            "fast_root": float(model.DN_gw[-1]),
            "reference_map_at_fast_root": ref_map,
            "root_relative_error_proxy": abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
            "reference_seconds": elapsed,
            "spectrum_valid_count": int(valid.sum()),
            "spectrum_max_abs_dex": float(np.max(np.abs(fast_log[valid]-ref_log[valid]))),
            "raw": {
                "freqs_log10_hz": freqs.tolist(),
                "fast_log10_OmegaGW": fast_log.tolist(),
                "reference_log10_OmegaGW": ref_log.tolist(),
                "reference_Ogw": np.asarray(ref["Ogw"]).tolist(),
                "reference_Oj": np.asarray(ref["Oj"]).tolist(),
                "valid_positive_reference_mask": valid.tolist(),
            },
        }
        rows.append(row)
        OUTPUT.write_text(json.dumps({"status":"in_progress","cases":rows},indent=2)+"\n",encoding="utf-8")
        print(case_id, row["root_relative_error_proxy"], row["spectrum_max_abs_dex"], flush=True)
    payload = {
        "status":"complete_cr1_nearby_reheat_kappa_reference_checks",
        "reference_code_modified":False,
        "validity_contract":"cr=1 model-derived tilt; finite N_inf; positive finite H2; physical_guard acceptance",
        "region":{"T_re":[1687.5,1750.0],"kappa10":[0.0125,0.015]},
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{name:sha256(ROOT/"stiffgwpy_fast"/name) for name in
                         ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")},
        "cases":rows,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
