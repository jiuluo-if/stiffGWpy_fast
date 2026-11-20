"""Selected direct z8 Reference checks for a CR0 high-density rectangle."""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "validate_cr0_high_density_extension_z8.json"
BASE_PARAMS = {
    "Omega_bh2": 2.0e6,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = (
    ("upper_density_zero_tilt_corner", 2.0e6, 0.0),
    ("upper_density_negative_tilt_corner", 2.0e6, -0.5),
    ("lower_density_zero_tilt_corner", 1.0e6, 0.0),
    ("independent_sobol_holdout", 1.5929892305284739e6, -0.1713485741056502),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _map(model: object, dn_eff: float, freqs: np.ndarray) -> dict[str, object]:
    started = time.perf_counter()
    result = run_reference(
        model,
        dn_eff=dn_eff,
        freq_subset=freqs,
        z_tail=8.0,
        rtol=1.0e-9,
        self_consistent=False,
        workers=2,
    )
    return {"seconds": time.perf_counter() - started, "result": result}


def main() -> None:
    rows: list[dict[str, object]] = []
    config = fast.resolve_config("fast", threads=2)
    for case_id, density, tilt in CASES:
        params = dict(BASE_PARAMS, Omega_bh2=density, n_t=tilt)
        model = LCDM_SG(**params)
        started = time.perf_counter()
        solved = fast.SGWB_iter_fast(
            model, config=config, tol=1.0e-6, physical_guard=True
        )
        row: dict[str, object] = {
            "case_id": case_id,
            "params": params,
            "derived_nt": float(model.derived_param["nt"]),
            "derived_N_inf": float(model.derived_param["N_inf"]),
            "fast_success": solved is model,
            "fast_failure_reason": model.fast_failure_reason,
            "fast_seconds": time.perf_counter() - started,
            "reference_z_tail": 8.0,
            "reference_rtol": 1.0e-9,
        }
        if solved is not model:
            rows.append(row)
            break

        freqs = np.asarray(model.f, dtype=np.float64)
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
        fast_root = float(model.DN_gw[-1])
        dn = float(model.cosmo_param["DN_eff"])
        first = _map(model, dn, freqs)
        first_result = first["result"]
        ref_log = np.asarray(first_result["log10OmegaGW"], dtype=np.float64)
        ref_linear = np.asarray(first_result["Ogw"]) - np.asarray(first_result["Oj"])
        valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
        abs_error = np.abs(fast_log[valid] - ref_log[valid])
        mapped = float(first_result["DN_gw"])
        map_history = [{"input_root": dn, "mapped_root": mapped, "seconds": first["seconds"]}]
        residual = math.inf
        for _ in range(4):
            closure = _map(model, mapped, freqs)
            new_mapped = float(closure["result"]["DN_gw"])
            map_history.append(
                {"input_root": mapped, "mapped_root": new_mapped, "seconds": closure["seconds"]}
            )
            residual = abs(new_mapped - mapped) / max(abs(new_mapped), 1.0e-300)
            mapped = new_mapped
            if residual < 2.0e-5:
                break
        row.update(
            {
                "fast_root": fast_root,
                "total_DN_eff": dn,
                "closed_reference_root": mapped,
                "root_relative_error": abs(mapped - fast_root) / max(abs(mapped), 1.0e-300),
                "reference_map_history": map_history,
                "reference_closure_residual": residual,
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": float(np.max(abs_error)),
                "spectrum_p95_abs_dex": float(np.quantile(abs_error, 0.95)),
                "raw": {
                    "freqs_log10_hz": freqs.tolist(),
                    "fast_log10_OmegaGW": fast_log.tolist(),
                    "reference_log10_OmegaGW_at_fast_root": ref_log.tolist(),
                    "reference_Ogw_at_fast_root": np.asarray(first_result["Ogw"]).tolist(),
                    "reference_Oj_at_fast_root": np.asarray(first_result["Oj"]).tolist(),
                    "valid_positive_reference_mask": valid.tolist(),
                },
            }
        )
        rows.append(row)
        payload = {
            "status": "in_progress",
            "reference_code_modified": False,
            "baseline_fast_screen": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
            "screen_candidate_sha": "702d54347e3c489d1632142dc3151cbd54dd6657",
            "fast_config": {
                "accuracy_mode": "fast",
                "h": float(config.h),
                "z_tail": float(config.z_tail),
                "phase_max": float(config.phase_max),
                "threads": 2,
                "frequency_grid": "goal",
                "tol": 1.0e-6,
            },
            "source_sha256": {
                name: _sha256(ROOT / "stiffgwpy_fast" / name)
                for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
            },
            "cases": rows,
        }
        OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(case_id, row.get("root_relative_error"), row.get("spectrum_max_abs_dex"), flush=True)
        if row["root_relative_error"] > 2.0e-4 or row["spectrum_max_abs_dex"] > 1.0e-3:
            break

    payload["status"] = "complete_selected_region_reference_checks"
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
