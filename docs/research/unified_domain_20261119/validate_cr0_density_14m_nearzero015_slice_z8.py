"""Selected z8 Reference checks for the 14M near-zero-tilt CR0 box."""
from __future__ import annotations

import hashlib
import importlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
fast = importlib.import_module("stiffgwpy_fast.fast_sgwb")
run_reference = importlib.import_module("stiffgwpy_fast.reference").run_reference
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
closure_spec = importlib.util.spec_from_file_location(
    "cr0_reference_origin_helper", HERE / "close_cr0_3m_reference_origin_corrected.py"
)
if closure_spec is None or closure_spec.loader is None:
    raise RuntimeError("could not load the tested Reference origin helper")
closure_module = importlib.util.module_from_spec(closure_spec)
closure_spec.loader.exec_module(closure_module)
reference_total_from_map = closure_module._reference_total_from_map
OUTPUT = HERE / "validate_cr0_density_14m_nearzero015_slice_z8.json"
BASE_PARAMS = {
    "Omega_ch2": 0.1201075, "H0": 67.32117, "DN_eff": 0.0,
    "A_s": 2.100549e-9, "r": 0.01, "cr": 0.0,
    "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}
CASES = (
    ("lower_density_low_tilt", 13.875e6, -0.015625),
    ("upper_density_low_tilt", 14.0e6, -0.015625),
    ("interior", 13.9375e6, -0.0078125),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_map(model: object, dn_eff: float, freqs: np.ndarray) -> tuple[dict[str, object], float]:
    started = time.perf_counter()
    result = run_reference(
        model, dn_eff=dn_eff, freq_subset=freqs, z_tail=8.0,
        rtol=1.0e-9, self_consistent=False, workers=2,
    )
    return result, time.perf_counter() - started


def main() -> None:
    rows: list[dict[str, object]] = []
    config = fast.resolve_config("fast", threads=2)
    for case_id, density, tilt in CASES:
        params = dict(BASE_PARAMS, Omega_bh2=density, n_t=tilt)
        model = LCDM_SG(**params)
        started = time.perf_counter()
        solved = fast.SGWB_iter_fast(model, config=config, tol=1.0e-6, physical_guard=True)
        row: dict[str, object] = {
            "case_id": case_id, "params": params,
            "n_inf": float(model.derived_param["N_inf"]),
            "fast_success": solved is model, "fast_failure_reason": model.fast_failure_reason,
            "fast_seconds": time.perf_counter() - started,
            "reference_z_tail": 8.0, "reference_rtol": 1.0e-9,
        }
        if solved is not model:
            rows.append(row)
            break
        freqs = np.asarray(model.f, dtype=np.float64)
        fast_log = np.asarray(model.log10OmegaGW, dtype=np.float64)
        root = float(model.DN_gw[-1])
        origin = float(getattr(model, "DN_eff_orig", params["DN_eff"]))
        total = float(model.cosmo_param["DN_eff"])
        first, elapsed = run_map(model, total, freqs)
        mapped = float(first["DN_gw"])
        history = [{"input_dn_eff": total, "mapped_dn_gw": mapped, "seconds": elapsed}]
        residual = abs(mapped - (total - origin)) / max(abs(mapped), 1.0e-300)
        total = reference_total_from_map(origin, mapped)
        closure_iteration_failure = None
        for _ in range(5):
            if residual < 2.0e-5:
                break
            try:
                result, elapsed = run_map(model, total, freqs)
            except Exception as exc:
                closure_iteration_failure = f"{type(exc).__name__}: {exc}"
                break
            mapped = float(result["DN_gw"])
            history.append({"input_dn_eff": total, "mapped_dn_gw": mapped, "seconds": elapsed})
            residual = abs(mapped - (total - origin)) / max(abs(mapped), 1.0e-300)
            total = reference_total_from_map(origin, mapped)
        ref_log = np.asarray(first["log10OmegaGW"], dtype=np.float64)
        positive = np.asarray(first["Ogw"]) - np.asarray(first["Oj"]) > 0.0
        valid = np.isfinite(fast_log) & np.isfinite(ref_log) & positive
        errors = np.abs(fast_log[valid] - ref_log[valid])
        root_error = abs(total - (origin + root)) / max(abs(total), 1.0e-300)
        row.update({
            "origin_DN_eff": origin, "fast_total_DN_eff": float(model.cosmo_param["DN_eff"]),
            "fast_root": root, "closed_reference_total_DN_eff": total,
            "closed_reference_root": total-origin, "root_relative_error": root_error,
            "map_closure_residual": residual, "reference_map_history": history,
            "closure_iteration_failure_after_converged_map": closure_iteration_failure,
            "spectrum_valid_count": int(valid.sum()),
            "spectrum_max_abs_dex": float(np.max(errors)),
            "spectrum_p95_abs_dex": float(np.quantile(errors, 0.95)),
            "raw": {
                "freqs_log10_hz": freqs.tolist(), "fast_log10_OmegaGW": fast_log.tolist(),
                "reference_log10_OmegaGW_at_fast_root": ref_log.tolist(),
                "reference_Ogw_at_fast_root": np.asarray(first["Ogw"]).tolist(),
                "reference_Oj_at_fast_root": np.asarray(first["Oj"]).tolist(),
                "valid_positive_reference_mask": valid.tolist(),
            },
        })
        rows.append(row)
        print(case_id, f"root={root_error:.6g}", f"spectrum={row['spectrum_max_abs_dex']:.6g}", flush=True)
        if root_error > 2.0e-4 or row["spectrum_max_abs_dex"] > 1.0e-3:
            break
    payload = {
        "status": "complete_selected_current_source_reference_checks" if rows and all(r.get("fast_success") for r in rows) else "partial_or_fast_failure",
        "closure_stop_rule": "stop after a reference map when relative fixed-point residual <2e-5",
        "reference_code_modified": False,
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "_wkb_transfer.py", "reference.py")
        },
        "source_base_commit": "d4ee1b9c8f9fad401c5e22c2316094a608074727",
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("saved", OUTPUT)


if __name__ == "__main__":
    main()
