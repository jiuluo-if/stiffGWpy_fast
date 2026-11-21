"""Reference checks for the adjacent CR1 high-r rectangle."""
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
LATTICE = json.loads((HERE / "screen_cr1_r_upperband_candidate_lattice.json").read_text(encoding="utf-8"))["cases"]
SOBOL = json.loads((HERE / "screen_cr1_r_upperband_candidate_sobol.json").read_text(encoding="utf-8"))["cases"]
LATTICE_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in LATTICE}
SOBOL_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in SOBOL}
OUTPUT = HERE / "validate_cr1_r_upperband_candidate_z8.json"
SOURCE_NAMES = ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")


def source_hashes() -> dict[str, str]:
    return {name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in SOURCE_NAMES}


def direct_log(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(signal) & (signal > 0.0)
    result = np.full(signal.shape, np.nan)
    result[valid] = np.log10(signal[valid])
    return result, valid


def make_row(case_id: str, params: dict[str, float], config) -> dict[str, object]:
    model = LCDM_SG(**params)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=np.float64)
    spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=freqs, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter()-started
    signal = np.asarray(result["Ogw"])-np.asarray(result["Oj"])
    ref_log, valid = direct_log(signal)
    ref_map = float(result["DN_gw"])
    derived = model.derived_param
    return {
        "case_id":case_id,"params":params,
        "derived_nt":float(derived["nt"]),"derived_N_inf":float(derived["N_inf"]),
        "fast_root":float(model.DN_gw[-1]),"reference_map_at_fast_root":ref_map,
        "root_relative_error_proxy":abs(ref_map-float(model.DN_gw[-1]))/abs(ref_map),
        "reference_seconds":elapsed,"spectrum_valid_count":int(valid.sum()),
        "spectrum_max_abs_dex_from_raw_linear":float(np.max(np.abs(spectrum[valid]-ref_log[valid]))),
        "reference_log_floor_bins":int(np.sum(np.asarray(result["log10OmegaGW"])==-40.0)),
        "raw":{
            "freqs_log10_hz":freqs.tolist(),"fast_log10_OmegaGW":spectrum.tolist(),
            "reference_log10_OmegaGW_flooring_api":np.asarray(result["log10OmegaGW"]).tolist(),
            "reference_raw_Ogw":np.asarray(result["Ogw"]).tolist(),
            "reference_raw_Oj":np.asarray(result["Oj"]).tolist(),
            "reference_direct_log10_Ogw_minus_Oj":ref_log.tolist(),
            "valid_positive_reference_mask":valid.tolist(),
        },
    }


def reused_row(artifact_name: str, selector, case_id: str) -> dict[str, object]:
    data = json.loads((HERE / artifact_name).read_text(encoding="utf-8"))
    if data.get("source_sha256") != source_hashes():
        raise RuntimeError(f"Reference source hash mismatch in {artifact_name}")
    row = selector(data)
    result = dict(row)
    result["case_id"] = case_id
    return result


def main() -> None:
    edge = json.loads((HERE / "validate_cr1_r_upperband_edge_probes_z8.json").read_text(encoding="utf-8"))
    if edge.get("source_sha256") != source_hashes():
        raise RuntimeError("upper-r edge probe source hash mismatch")
    lowstrip = json.loads((HERE / "validate_cr1_r_lowstrip_z8.json").read_text(encoding="utf-8"))
    if lowstrip.get("source_sha256") != source_hashes():
        raise RuntimeError("lower-r Reference source hash mismatch")
    old_audit = json.loads((HERE / "audit_cr1_r_upperband_old_success_z8.json").read_text(encoding="utf-8"))
    if old_audit.get("reference_sha256") != source_hashes()["reference.py"]:
        raise RuntimeError("old-success Reference source hash mismatch")
    old_case = dict(old_audit["current"])
    old_case["case_id"] = "old_fast_success_holdout"
    old_case["params"] = old_audit["params"]
    model = LCDM_SG(**old_audit["params"])
    old_case["derived_nt"] = float(model.derived_param["nt"])
    old_case["derived_N_inf"] = float(model.derived_param["N_inf"])
    old_case["spectrum_valid_count"] = old_case.get("spectrum_valid_count", 76)
    old_case["spectrum_max_abs_dex_from_raw_linear"] = old_case["spectrum_max_abs_dex_from_raw_linear"]

    rows = [
        dict(next(row for row in edge["cases"] if row["case_id"] == "T1875_high_r"), case_id="low_T_high_r_corner"),
        dict(next(row for row in edge["cases"] if row["case_id"] == "T2000_high_r"), case_id="high_T_high_r_corner"),
        dict(next(row for row in lowstrip["cases"] if row["case_id"] == "high_T_high_r_corner"), case_id="high_T_low_r_corner"),
        old_case,
    ]
    config = fast.resolve_config("fast", threads=2)
    new_cases = [
        ("low_T_low_r_corner", LATTICE_BY_ID["rectangle_00"]),
        ("independent_sobol_holdout", SOBOL_BY_ID["holdout_02"]),
    ]
    for case_id, params in new_cases:
        row = make_row(case_id, params, config)
        rows.append(row)
        OUTPUT.write_text(json.dumps({"status":"in_progress","cases":rows},indent=2)+"\n",encoding="utf-8")
        print(case_id,row["root_relative_error_proxy"],row["spectrum_max_abs_dex_from_raw_linear"],flush=True)
    payload={
        "status":"complete_cr1_upperband_reference_checks",
        "reference_code_modified":False,
        "validity_contract":"cr=1; nt=-r/8; finite N_inf; positive finite H2; physical_guard acceptance",
        "region":{"T_re":[1875.0,2000.0],"r":[0.005,0.01],"kappa10":0.01},
        "reference_log_interpretation":"For entries below the Reference -40 log floor, use log10(Ogw-Oj) from preserved positive raw linear outputs.",
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":source_hashes(),"cases":rows,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
