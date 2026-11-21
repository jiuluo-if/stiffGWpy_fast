"""Representative z8 Reference checks for the CR0 density/reheat tilt region."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
fast = __import__("stiffgwpy_fast.fast_sgwb", fromlist=["resolve_config"])
reference = __import__("stiffgwpy_fast.reference", fromlist=["run_reference"]).run_reference
LCDM_SG = __import__("stiffgwpy_fast.stiff_SGWB", fromlist=["LCDM_SG"]).LCDM_SG
OUTPUT = HERE / "validate_cr0_density_reheat_negative_tilt_z8.json"
SOURCES = ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")


def hashes() -> dict[str, str]:
    return {
        name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
        for name in SOURCES
    }


def main() -> None:
    inputs = [
        ("low_T_low_tilt_corner", "screen_cr0_density_reheat_negative_tilt_lattice.json", "rectangle_00"),
        ("low_T_high_tilt_corner", "screen_cr0_density_reheat_negative_tilt_lattice.json", "rectangle_02"),
        ("high_T_low_tilt_corner", "screen_cr0_density_reheat_negative_tilt_lattice.json", "rectangle_06"),
        ("high_T_high_tilt_corner", "screen_cr0_density_reheat_negative_tilt_lattice.json", "rectangle_08"),
        ("old_success_lattice_holdout", "screen_cr0_density_reheat_negative_tilt_lattice.json", "holdout_01"),
        ("old_success_sobol_holdout", "screen_cr0_density_reheat_negative_tilt_sobol.json", "holdout_06"),
        ("independent_sobol_holdout", "screen_cr0_density_reheat_negative_tilt_sobol.json", "holdout_02"),
    ]
    selected = []
    for label, filename, case_id in inputs:
        data = json.loads((HERE / filename).read_text(encoding="utf-8"))
        case = next(row for row in data["cases"] if row["input"]["case_id"] == case_id)
        selected.append((label, case["input"]["params"]))

    config = fast.resolve_config("fast", threads=2)
    rows: list[dict[str, object]] = []
    for label, params in selected:
        model = LCDM_SG(**params)
        if fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True) is not model:
            raise RuntimeError(f"Fast failed at {label}: {model.fast_failure_reason}")
        freqs = np.asarray(model.f, dtype=np.float64)
        started = time.perf_counter()
        result = reference(
            model,
            dn_eff=float(model.cosmo_param["DN_eff"]),
            freq_subset=freqs,
            z_tail=8.0,
            rtol=1e-9,
            self_consistent=False,
            workers=2,
        )
        elapsed = time.perf_counter() - started
        signal = np.asarray(result["Ogw"]) - np.asarray(result["Oj"])
        valid = np.isfinite(signal) & (signal > 0.0)
        ref_log = np.full(signal.shape, np.nan)
        ref_log[valid] = np.log10(signal[valid])
        fast_log = np.asarray(model.log10OmegaGW)
        row = {
            "case_id": label,
            "params": params,
            "derived_nt": float(model.derived_param["nt"]),
            "derived_N_inf": float(model.derived_param["N_inf"]),
            "fast_root": float(model.DN_gw[-1]),
            "reference_map_at_fast_root": float(result["DN_gw"]),
            "root_relative_error_proxy": abs(float(result["DN_gw"]) - float(model.DN_gw[-1])) / abs(float(result["DN_gw"])),
            "reference_seconds": elapsed,
            "spectrum_valid_count": int(valid.sum()),
            "spectrum_max_abs_dex_from_raw_linear": float(np.max(np.abs(fast_log[valid] - ref_log[valid]))),
            "reference_log_floor_bins": int(np.sum(np.asarray(result["log10OmegaGW"]) == -40.0)),
            "raw": {
                "freqs_log10_hz": freqs.tolist(),
                "fast_log10_OmegaGW": fast_log.tolist(),
                "reference_raw_Ogw": np.asarray(result["Ogw"]).tolist(),
                "reference_raw_Oj": np.asarray(result["Oj"]).tolist(),
                "reference_direct_log10_Ogw_minus_Oj": ref_log.tolist(),
                "valid_positive_reference_mask": valid.tolist(),
            },
        }
        rows.append(row)
        OUTPUT.write_text(json.dumps({"status": "in_progress", "cases": rows}, indent=2) + "\n", encoding="utf-8")
        print(label, row["root_relative_error_proxy"], row["spectrum_max_abs_dex_from_raw_linear"], flush=True)

    payload = {
        "status": "complete_cr0_density_reheat_negative_tilt_reference_checks",
        "reference_code_modified": False,
        "validity_contract": "cr=0; n_t independent; positive finite H2; finite N_inf; physical_guard accepted; n_t<=0",
        "region": {
            "Omega_bh2": 1.0e6, "Omega_ch2": 0.1201075,
            "T_re_GeV": [4.0e3, 1.0e7], "n_t": [-0.5, -0.25], "kappa10": 0.01,
        },
        "fast_config": {"threads": 2, "frequency_grid": "goal", "reference_z_tail": 8.0, "reference_rtol": 1e-9},
        "source_sha256": hashes(),
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
