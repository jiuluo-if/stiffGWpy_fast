"""z8 Reference checks for the CR1 lower-r adjacent rectangle."""
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
LATTICE = json.loads((HERE / "screen_cr1_r_lowstrip_lattice.json").read_text(encoding="utf-8"))["cases"]
SOBOL = json.loads((HERE / "screen_cr1_r_lowstrip_sobol.json").read_text(encoding="utf-8"))["cases"]
LATTICE_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in LATTICE}
SOBOL_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in SOBOL}
PROBE = json.loads((HERE / "validate_cr1_r_lowstrip_upper_probe_z8.json").read_text(encoding="utf-8"))
CASES = [
    ("low_T_low_r_corner", LATTICE_BY_ID["rectangle_00"]),
    ("high_T_low_r_corner", LATTICE_BY_ID["rectangle_06"]),
    ("high_T_high_r_corner", LATTICE_BY_ID["rectangle_08"]),
    ("independent_sobol_holdout", SOBOL_BY_ID["holdout_00"]),
]
OUTPUT = HERE / "validate_cr1_r_lowstrip_z8.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_reference_case(case_id: str, params: dict[str, float], config) -> dict[str, object]:
    model = LCDM_SG(**params)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter() - started
    ref_log = np.asarray(result["log10OmegaGW"], dtype=np.float64)
    valid = (np.isfinite(fast_log) & np.isfinite(ref_log)
             & (np.asarray(result["Ogw"]) - np.asarray(result["Oj"]) > 0.0))
    ref_map = float(result["DN_gw"])
    return {
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
            "reference_Ogw": np.asarray(result["Ogw"]).tolist(),
            "reference_Oj": np.asarray(result["Oj"]).tolist(),
            "valid_positive_reference_mask": valid.tolist(),
        },
    }


def main() -> None:
    config = fast.resolve_config("fast", threads=2)
    probe_params = LATTICE_BY_ID["rectangle_02"]
    probe_model = LCDM_SG(**probe_params)
    if fast.SGWB_iter_fast(probe_model, config=config, tol=1e-6, physical_guard=True) is not probe_model:
        raise RuntimeError(f"Fast failed at saved upper probe: {probe_model.fast_failure_reason}")
    probe_row = PROBE["case"]
    if (float(probe_model.DN_gw[-1]) != probe_row["fast_root"]
            or not np.array_equal(np.asarray(probe_model.f), np.asarray(probe_row["raw"]["freqs_log10_hz"]))
            or not np.array_equal(np.asarray(probe_model.log10OmegaGW), np.asarray(probe_row["raw"]["fast_log10_OmegaGW"]))):
        raise RuntimeError("saved upper-probe Fast inputs or outputs no longer match current source")
    rows = [probe_row]
    for case_id, params in CASES:
        row = run_reference_case(case_id, params, config)
        rows.append(row)
        OUTPUT.write_text(json.dumps({"status":"in_progress","cases":rows},indent=2)+"\n",encoding="utf-8")
        print(case_id, row["root_relative_error_proxy"], row["spectrum_max_abs_dex"], flush=True)
    payload = {
        "status":"complete_cr1_r_lowstrip_reference_checks",
        "reference_code_modified":False,
        "validity_contract":"cr=1; derived nt=-r/8; finite N_inf; positive finite H2; physical_guard acceptance",
        "region":{"T_re":[1750.0,2000.0],"r":[0.001,0.005],"kappa10":0.01},
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{name:sha256(ROOT/"stiffgwpy_fast"/name) for name in
                         ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")},
        "cases":rows,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
