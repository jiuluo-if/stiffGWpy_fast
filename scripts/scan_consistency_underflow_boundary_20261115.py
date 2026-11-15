"""Scan cr>0 underflow cutoffs and a separate live S2-overflow control."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import _has_exact_zero_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 1e-316,
    "n_t": 0.0,
    "cr": 1.0,
    "DN_re": 20.0,
    "DN_eff": 0.0,
}
OUTPUT = ROOT / "docs" / "research" / "results" / "consistency_underflow_boundary_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_case(label: str, params: dict[str, float]) -> dict[str, object]:
    model = LCDM_SG(**params)
    started = time.perf_counter()
    derived_error = None
    try:
        derived = model.derived_param
    except Exception as exc:
        derived = {}
        derived_error = f"{type(exc).__name__}: {exc}"
    derived_seconds = time.perf_counter() - started
    log_product = None
    if params["r"] > 0.0 and params["A_s"] > 0.0:
        log_product = math.log(params["r"]) + math.log(params["A_s"])
    start_fast = time.perf_counter()
    error = None
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=0.99, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=False,
        )
    except Exception as exc:
        result = None
        error = f"{type(exc).__name__}: {exc}"
    fast_seconds = time.perf_counter() - start_fast
    if derived_error:
        stage = "derived_background"
    elif model.fast_failure_reason == "nonfinite_background_primitive":
        stage = f"primitive:{getattr(model, '_fast_background_failure_detail', 'unknown')}"
    elif model.fast_failure_reason in ("invalid_r", "invalid_cutoff"):
        stage = "pre_primitive_physical_guard"
    elif error:
        stage = "fast_solver_exception"
    elif result is None:
        stage = "other_fast_rejection"
    else:
        stage = "fast_completed"
    d = derived
    delta_n = None
    if log_product is not None and params["cr"] > 0 and d.get("N_re") is not None:
        delta_n = (
            (gp.N_10 - float(d["N_re"])) * 4.0 / 3.0
            + math.log(gp.M_Pl) * 4.0 / 3.0
            + (math.log(45.0 / 2.0) + log_product
               - math.log(model.rhorad_re + model.rhostiff_re)) / 3.0
        )
    return {
        "label": label,
        "parameters": params,
        "exact_zero_tensor_limit": bool(_has_exact_zero_tensor_amplitude(model.cosmo_param)),
        "A_t_stored": float(d.get("A_t", params["r"] * params["A_s"])),
        "log_r_As": log_product,
        "V_inf": d.get("V_inf"),
        "Delta_N_from_log_product": delta_n,
        "N_inf": d.get("N_inf"),
        "derived_error": derived_error,
        "derived_runtime_s": derived_seconds,
        "fast_returned_model": result is not None,
        "fast_failure_reason": model.fast_failure_reason,
        "fast_error": error,
        "failure_stage": stage,
        "background_failure_detail": getattr(model, "_fast_background_failure_detail", None),
        "background_failure_index": getattr(model, "_fast_background_failure_index", None),
        "background_failure_log_scale": getattr(model, "_fast_background_failure_log_scale", None),
        "fast_runtime_s": fast_seconds,
    }


def main() -> int:
    cases: list[dict[str, object]] = []
    for t_re in (10.0, 2e3, 1e4):
        for kappa in (1e-3, 1e-2, 1.0):
            cases.append(run_case(
                f"underflow_T{t_re:g}_k{kappa:g}",
                {**BASE, "T_re": t_re, "kappa10": kappa},
            ))
    # This known background-only stress point has positive N_inf and must reach
    # primitive construction; it is separate from the invalid-cutoff samples.
    cases.append(run_case(
        "normal_amplitude_S2_overflow_control",
        {
            **BASE,
            "r": 1e-2,
            "n_t": 0.5,
            "T_re": 2e3,
            "kappa10": 1e-2,
            "Omega_bh2": 1e6,
        },
    ))
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_underflow_boundary_and_s2_control",
        "frequency_grid": "goal",
        "fast_options": {"tol": 0.99, "z_tail": 5.0, "threads": 1},
        "cases": cases,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "diagnostic_script": "scripts/scan_consistency_underflow_boundary_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "label": row["label"],
            "Delta_N_from_log_product": row["Delta_N_from_log_product"],
            "N_inf": row["N_inf"],
            "failure_stage": row["failure_stage"],
            "background_failure_detail": row["background_failure_detail"],
        } for row in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
