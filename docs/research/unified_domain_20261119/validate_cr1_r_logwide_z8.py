"""Reference checks spanning the CR1 cutoff-to-observable-r rectangle."""
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
LATTICE = json.loads((HERE / "screen_cr1_r_logwide_lattice.json").read_text(encoding="utf-8"))["cases"]
SOBOL = json.loads((HERE / "screen_cr1_r_logwide_sobol.json").read_text(encoding="utf-8"))["cases"]
LATTICE_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in LATTICE}
SOBOL_BY_ID = {row["input"]["case_id"]: row["input"]["params"] for row in SOBOL}
NEW_CASES = [
    ("low_r_cutoff_low_T_corner", LATTICE_BY_ID["rectangle_00"]),
    ("old_fast_success_mid_log_r", LATTICE_BY_ID["rectangle_07"]),
    ("independent_sobol_low_r_holdout", SOBOL_BY_ID["holdout_02"]),
]
OUTPUT = HERE / "validate_cr1_r_logwide_z8.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_row(case_id: str, params: dict[str, float], config) -> dict[str, object]:
    model = LCDM_SG(**params)
    if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
        raise RuntimeError(f"Fast failed at {case_id}: {model.fast_failure_reason}")
    frequencies = np.asarray(model.f, dtype=np.float64)
    fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
    started = time.perf_counter()
    result = reference(model, dn_eff=float(model.cosmo_param["DN_eff"]),
                       freq_subset=frequencies, z_tail=8.0, rtol=1e-9,
                       self_consistent=False, workers=2)
    elapsed = time.perf_counter()-started
    integrand = np.asarray(result["Ogw"]) - np.asarray(result["Oj"])
    valid = np.isfinite(fast_log) & np.isfinite(integrand) & (integrand > 0.0)
    direct_log = np.full(integrand.shape, np.nan)
    direct_log[valid] = np.log10(integrand[valid])
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
        "spectrum_max_abs_dex_from_raw_linear": float(np.max(np.abs(fast_log[valid]-direct_log[valid]))),
        "reference_log_floor_bins": int(np.sum(np.asarray(result["log10OmegaGW"]) == -40.0)),
        "raw": {
            "freqs_log10_hz": frequencies.tolist(),
            "fast_log10_OmegaGW": fast_log.tolist(),
            "reference_log10_OmegaGW_flooring_api": np.asarray(result["log10OmegaGW"]).tolist(),
            "reference_raw_Ogw": np.asarray(result["Ogw"]).tolist(),
            "reference_raw_Oj": np.asarray(result["Oj"]).tolist(),
            "reference_direct_log10_Ogw_minus_Oj": direct_log.tolist(),
            "valid_positive_reference_mask": valid.tolist(),
        },
    }


def reuse_case(artifact_name: str, case_id: str | None = None) -> dict[str, object]:
    data = json.loads((HERE / artifact_name).read_text(encoding="utf-8"))
    current_hashes = {
        name: sha256(ROOT / "stiffgwpy_fast" / name)
        for name in ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")
    }
    if data.get("source_sha256") != current_hashes:
        raise RuntimeError(f"source hash mismatch in reusable Reference artifact {artifact_name}")
    if "case" in data:
        row = data["case"]
        if "spectrum_max_abs_dex_from_raw_linear" not in row:
            raw = row["raw"]
            valid = np.asarray(raw["valid_positive_reference_mask"], dtype=bool)
            signal = np.asarray(raw["reference_Ogw"]) - np.asarray(raw["reference_Oj"])
            direct = np.log10(signal[valid])
            row["spectrum_max_abs_dex_from_raw_linear"] = float(np.max(np.abs(np.asarray(raw["fast_log10_OmegaGW"])[valid]-direct)))
            row["reference_log_floor_bins"] = int(np.sum(np.asarray(raw["reference_log10_OmegaGW"]) == -40.0))
        return row
    row = next(item for item in data["cases"] if item["case_id"] == case_id)
    raw = row["raw"]
    valid = np.asarray(raw["valid_positive_reference_mask"], dtype=bool)
    signal = np.asarray(raw["reference_Ogw"]) - np.asarray(raw["reference_Oj"])
    direct = np.log10(signal[valid])
    row["spectrum_max_abs_dex_from_raw_linear"] = float(np.max(np.abs(np.asarray(raw["fast_log10_OmegaGW"])[valid]-direct)))
    row["reference_log_floor_bins"] = int(np.sum(np.asarray(raw["reference_log10_OmegaGW"] ) == -40.0))
    return row


def main() -> None:
    config = fast.resolve_config("fast", threads=2)
    rows = [
        reuse_case("validate_cr1_r_cutoff_probe_z8.json"),
        reuse_case("validate_cr1_r_lowstrip_z8.json", "low_T_upper_r_corner"),
        reuse_case("validate_cr1_r_lowstrip_z8.json", "high_T_high_r_corner"),
        reuse_case("validate_cr1_r_log_tail_probe_z8.json"),
    ]
    for case_id, params in NEW_CASES:
        row = make_row(case_id, params, config)
        rows.append(row)
        OUTPUT.write_text(json.dumps({"status":"in_progress","cases":rows},indent=2)+"\n",encoding="utf-8")
        print(case_id, row["root_relative_error_proxy"], row["spectrum_max_abs_dex_from_raw_linear"], flush=True)
    payload = {
        "status":"complete_cr1_logwide_reference_checks",
        "reference_code_modified":False,
        "reference_log_interpretation":"Reference log output has a -40 floor; below it, spectrum errors use log10(Ogw-Oj) from preserved positive raw linear outputs.",
        "validity_contract":"cr=1; nt=-r/8; finite N_inf with N_inf>=N_re; positive finite H2; physical_guard acceptance",
        "region":{"T_re":[1750.0,2000.0],"r":[4e-43,0.005],"kappa10":0.01,"r_sampling":"log10 uniform"},
        "fast_config":{"threads":2,"frequency_grid":"goal","reference_z_tail":8.0,"reference_rtol":1e-9},
        "source_sha256":{name:sha256(ROOT/"stiffgwpy_fast"/name) for name in
                         ("exact_background.py","fast_sgwb.py","freq_adaptive.py","reference.py")},
        "cases":rows,
    }
    OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
