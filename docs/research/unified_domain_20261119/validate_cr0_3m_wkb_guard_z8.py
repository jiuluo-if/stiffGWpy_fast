"""Selected z8 Reference checks after the singular-basis fallback."""

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
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "validate_cr0_3m_wkb_guard_z8_multi.json"
BASE_PARAMS = {
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = (
    ("upper_density_zero_tilt", 3.0e6, 0.0),
    ("upper_density_negative_tilt", 3.0e6, -0.5),
    ("independent_sobol_upper_density", 2.989492731168866e6, -0.01380485063418746),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        root_origin = float(model.cosmo_param["DN_eff"])
        dn = float(model.cosmo_param["DN_eff"])
        ref_started = time.perf_counter()
        fast_root_result = run_reference(
            model,
            dn_eff=dn,
            freq_subset=freqs,
            z_tail=8.0,
            rtol=1.0e-9,
            self_consistent=False,
            workers=2,
        )
        mapped = float(fast_root_result["DN_gw"])
        history: list[dict[str, float]] = [
            {
                "input_dn_eff": dn,
                "mapped_dn_gw": mapped,
                "seconds": time.perf_counter() - ref_started,
            }
        ]
        residual = abs(mapped - (dn - root_origin)) / max(abs(mapped), 1e-300)
        dn = root_origin + mapped
        for _ in range(5):
            if residual < 2.0e-5:
                break
            ref_started = time.perf_counter()
            result = run_reference(
                model,
                dn_eff=dn,
                freq_subset=freqs,
                z_tail=8.0,
                rtol=1.0e-9,
                self_consistent=False,
                workers=2,
            )
            mapped = float(result["DN_gw"])
            history.append(
                {
                    "input_dn_eff": dn,
                    "mapped_dn_gw": mapped,
                    "seconds": time.perf_counter() - ref_started,
                }
            )
            residual = abs(mapped - (dn - root_origin)) / max(abs(mapped), 1e-300)
            dn = root_origin + mapped
            if residual < 2.0e-5:
                break

        ref_log = np.asarray(fast_root_result["log10OmegaGW"], dtype=np.float64)
        ref_ogw = np.asarray(fast_root_result["Ogw"], dtype=np.float64)
        ref_oj = np.asarray(fast_root_result["Oj"], dtype=np.float64)
        ref_linear = ref_ogw - ref_oj
        valid = np.isfinite(fast_log) & np.isfinite(ref_log) & (ref_linear > 0.0)
        error = np.abs(fast_log[valid] - ref_log[valid])
        row.update(
            {
                "fast_root": fast_root,
                "closed_reference_root": dn - root_origin,
                "root_relative_error": abs(dn - (root_origin + fast_root)) / max(abs(dn), 1e-300),
                "reference_closure_residual": residual,
                "reference_map_history": history,
                "spectrum_valid_count": int(valid.sum()),
                "spectrum_max_abs_dex": float(np.max(error)),
                "spectrum_p95_abs_dex": float(np.quantile(error, 0.95)),
                "raw": {
                    "freqs_log10_hz": freqs.tolist(),
                    "fast_log10_OmegaGW": fast_log.tolist(),
                    "reference_log10_OmegaGW_at_fast_root": ref_log.tolist(),
                    "reference_Ogw_at_fast_root": ref_ogw.tolist(),
                    "reference_Oj_at_fast_root": ref_oj.tolist(),
                    "valid_positive_reference_mask": valid.tolist(),
                },
            }
        )
        rows.append(row)
        OUTPUT.write_text(
            json.dumps(
                {
                    "status": "in_progress",
                    "reference_code_modified": False,
                    "cases": rows,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        print(
            case_id,
            row["root_relative_error"],
            row["spectrum_max_abs_dex"],
            flush=True,
        )
        if row["root_relative_error"] > 2.0e-4 or row["spectrum_max_abs_dex"] > 1.0e-3:
            break

    payload = {
        "status": "complete_selected_region_reference_checks",
        "reference_code_modified": False,
        "baseline_fast_sha": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "source_sha256": {
            name: _sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "_wkb_transfer.py", "reference.py")
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"raw artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
