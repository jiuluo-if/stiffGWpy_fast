"""Current-source z8 Reference checks at narrowed CR1 reheating points."""
from __future__ import annotations

import hashlib
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
OUTPUT = HERE / "validate_cr1_b175_reheat_kappa_narrow_z8.json"
LATTICE = json.loads(
    (HERE / "screen_cr1_b175_reheat_kappa_lower_lattice.json").read_text(
        encoding="utf-8"
    )
)["cases"]
SOBOL = json.loads(
    (HERE / "screen_cr1_b175_reheat_kappa_narrow_sobol.json").read_text(
        encoding="utf-8"
    )
)["cases"]


def _case_at(t_re: float, kappa: float) -> dict[str, object]:
    row = next(
        item for item in LATTICE
        if item["input"]["params"]["T_re"] == t_re
        and item["input"]["params"]["kappa10"] == kappa
    )
    return row["input"]["params"]


CASES = [
    ("narrow_low_reheat_low_kappa", _case_at(1750.0, 0.01)),
    ("narrow_low_reheat_high_kappa", _case_at(1750.0, 0.015)),
    ("narrow_sobol_holdout", SOBOL[4]["input"]["params"]),
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    config = fast.resolve_config("fast", threads=2)
    for case_id, params in CASES:
        model = LCDM_SG(**params)
        derived = model.derived_param
        row = {
            "case_id": case_id,
            "params": params,
            "consistency_tilt": float(derived["nt"]),
            "derived_N_inf": None if derived["N_inf"] is None else float(derived["N_inf"]),
        }
        if derived["N_inf"] is None:
            row.update({"fast_success": False, "fast_failure_reason": "invalid_cutoff"})
            rows.append(row)
            continue
        started = time.perf_counter()
        solved = fast.SGWB_iter_fast(model, config=config, tol=1e-6, physical_guard=True)
        row.update({
            "fast_success": solved is model,
            "fast_failure_reason": model.fast_failure_reason,
            "fast_seconds": time.perf_counter() - started,
        })
        if solved is model:
            freqs = np.asarray(model.f, dtype=np.float64)
            fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
            total_dn_eff = float(model.cosmo_param["DN_eff"])
            started = time.perf_counter()
            ref = reference(
                model, dn_eff=total_dn_eff, freq_subset=freqs, z_tail=8.0,
                rtol=1e-9, self_consistent=False, workers=2,
            )
            elapsed = time.perf_counter() - started
            ref_log = np.asarray(ref["log10OmegaGW"], dtype=np.float64)
            positive = np.asarray(ref["Ogw"]) - np.asarray(ref["Oj"]) > 0.0
            valid = np.isfinite(fast_log) & np.isfinite(ref_log) & positive
            errors = np.abs(fast_log[valid] - ref_log[valid])
            ref_map = float(ref["DN_gw"])
            row.update({
                "fast_root": float(model.DN_gw[-1]),
                "reference_map_at_fast_root": ref_map,
                "fixed_point_residual": ref_map - total_dn_eff,
                "root_relative_error_proxy": abs(ref_map - float(model.DN_gw[-1]))
                / max(abs(ref_map), 1e-300),
                "reference_seconds": elapsed,
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": float(np.max(errors)),
                "spectrum_p95_abs_dex": float(np.quantile(errors, 0.95)),
                "raw": {
                    "freqs_log10_hz": freqs.tolist(),
                    "fast_log10_OmegaGW": fast_log.tolist(),
                    "reference_log10_OmegaGW": ref_log.tolist(),
                    "reference_Ogw": np.asarray(ref["Ogw"]).tolist(),
                    "reference_Oj": np.asarray(ref["Oj"]).tolist(),
                    "valid_positive_reference_mask": valid.tolist(),
                },
            })
        rows.append(row)
        payload = {
            "status": "in_progress",
            "reference_code_modified": False,
            "validity_contract": "cr=1 consistency-derived tilt; finite derived N_inf; positive finite H2; physical guard",
            "region": {"T_re": [1750.0, 2000.0], "kappa10": [0.01, 0.015]},
            "fast_config": {
                "h": float(config.h), "z_tail": float(config.z_tail),
                "phase_max": float(config.phase_max), "threads": 2,
                "frequency_grid": "goal", "reference_z_tail": 8.0,
                "reference_rtol": 1e-9,
            },
            "source_sha256": {
                name: _sha256(ROOT / "stiffgwpy_fast" / name)
                for name in ("exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "reference.py")
            },
            "cases": rows,
        }
        OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(case_id, row.get("root_relative_error_proxy"),
              row.get("spectrum_max_abs_dex"), flush=True)
    payload["status"] = "complete_cr1_narrow_reheat_kappa_reference_checks"
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
