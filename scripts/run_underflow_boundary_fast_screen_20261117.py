"""English description: Run guarded Fast solves for the fixed-background amplitude strata. 中文说明：对固定背景振幅分层执行带物理守卫的 Fast 预筛。"""
from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_preflight_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_fast_screen_20261117.json"
SOURCE_FILES = (
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/reference.py",
    "stiffgwpy_fast/stiff_SGWB.py",
    "stiffgwpy_fast/LCDM_stiff_Neff.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def run_case(row: dict[str, object]) -> dict[str, object]:
    params = {key: value for key, value in row["parameters"].items()}
    model = LCDM_SG(**params)
    started = time.perf_counter()
    error = None
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = model.SGWB_iter(
                engine="fast",
                fallback=False,
                tol=1e-6,
                z_tail=8.0,
                threads=1,
                accuracy_mode="fast",
                freq_grid="goal",
                physical_guard=True,
            )
        except Exception as exc:
            result = None
            error = f"{type(exc).__name__}: {exc}"
    elapsed = time.perf_counter() - started
    logs = np.asarray(getattr(model, "log10OmegaGW", []), dtype=np.float64)
    mapped_values = np.asarray(getattr(model, "DN_gw", []), dtype=np.float64)
    mapped = float(mapped_values[-1]) if mapped_values.size else None
    root = float(model.cosmo_param["DN_eff"])
    input_dn = float(params["DN_eff"])
    result_success = result is model
    residual = None
    if result_success and mapped is not None:
        final_increment = root - input_dn
        residual = (
            0.0 if final_increment == 0.0 and mapped == 0.0
            else (mapped - final_increment) / max(abs(final_increment), np.nextafter(0.0, 1.0))
        )
    actual_f = np.asarray(getattr(model, "f", []), dtype=np.float64)
    return {
        "case_id": row["case_id"],
        "parameters": params,
        "amplitude_label": row["amplitude_label"],
        "A_t_stored_before_solve": row["A_t_stored"],
        "preflight_frequency_grid_used": row["frequency_grid_used"],
        "preflight_frequency_count": row["frequency_count"],
        "preflight_frequency_sha256_f64le": row["frequency_sha256_f64le"],
        "fast_returned_model": result_success,
        "fast_failure_reason": getattr(model, "fast_failure_reason", None),
        "fast_failure_detail": getattr(model, "fast_failure_detail", None),
        "exception": error,
        "runtime_s": elapsed,
        "DN_eff_final": root,
        "mapped_DN_gw": mapped,
        "root_map_residual_relative": residual,
        "frequency_grid_used": getattr(model, "freq_grid_used", None),
        "frequency_count": int(actual_f.size),
        "frequency_sha256_f64le": digest(actual_f) if actual_f.size else None,
        "frequency_log10_hz": actual_f.tolist(),
        "frequency_matches_preflight": bool(
            actual_f.size == row["frequency_count"]
            and digest(actual_f) == row["frequency_sha256_f64le"]
        ),
        "log_spectrum_finite_count": int(np.count_nonzero(np.isfinite(logs))),
        "log_spectrum_negative_infinity_count": int(np.count_nonzero(np.isneginf(logs))),
        "log_spectrum_nan_count": int(np.count_nonzero(np.isnan(logs))),
        "log_spectrum_max": float(np.nanmax(logs)) if logs.size and np.isfinite(logs).any() else None,
        "log10_OmegaGW": [
            float(value) if math.isfinite(float(value)) else None for value in logs
        ],
        "runtime_warnings": [str(item.message) for item in caught],
    }


def main() -> int:
    preflight = json.loads(INPUT.read_text(encoding="utf-8"))
    if preflight.get("status") != "complete":
        raise RuntimeError("Boundary preflight is not complete")
    for rel, expected in preflight["source_sha256"].items():
        if sha256(ROOT / rel) != expected:
            raise RuntimeError(f"Source changed after preflight: {rel}")
    selected_labels = {
        "normal_edge",
        "subnormal_edge",
        "subnormal_interior",
        "rounded_zero_product",
        "min_subnormal_r",
        "exact_zero",
    }
    test_tilts = {1.6, 17.38}
    selected = [
        row for row in preflight["cases"]
        if row["actual_n_t"] in test_tilts
        and row["parameters"]["T_re"] == 2000.0
        and row["parameters"]["kappa10"] == 0.01
        and row["amplitude_label"] in selected_labels
    ]
    if len(selected) != len(test_tilts) * len(selected_labels):
        raise RuntimeError("Expected one preflight row per tilt and amplitude stratum")
    results = [run_case(row) for row in selected]
    output = {
        "schema_version": 1,
        "status": "complete_fast_guarded_amplitude_strata_screen",
        "scope": "Six cr=0 amplitude strata at n_t=1.6 and 17.38, T_re=2000,kappa10=.01; one thread, current fast profile, physical guard enabled. This is a screen, not Reference certification or a domain claim.",
        "preflight_artifact": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "preflight_sha256": sha256(INPUT),
        "source_sha256": {rel: sha256(ROOT / rel) for rel in SOURCE_FILES},
        "runner_sha256": sha256(Path(__file__)),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "cases": results,
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": output["status"],
        "cases": [
            {
                key: row[key]
                for key in (
                    "case_id", "fast_returned_model", "fast_failure_reason", "exception",
                    "DN_eff_final", "mapped_DN_gw", "root_map_residual_relative",
                    "frequency_grid_used", "frequency_count", "frequency_matches_preflight",
                    "log_spectrum_finite_count", "log_spectrum_negative_infinity_count",
                    "runtime_warnings", "runtime_s",
                )
            }
            for row in results
        ],
        "output": str(OUTPUT),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
