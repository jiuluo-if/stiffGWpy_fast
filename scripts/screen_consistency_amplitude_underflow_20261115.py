"""Classify amplitude-zero and underflow behavior across cr branches."""

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
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import _has_exact_zero_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "T_re": 2000.0,
    "DN_re": 20.0,
    "kappa10": 0.01,
    "DN_eff": 0.0,
}
CASES = (
    ("cr0_nonzero_product_underflow", 0.0, 1e-316),
    ("cr0_exact_zero", 0.0, 0.0),
    ("cr1_nonzero_product_underflow", 1.0, 1e-316),
    ("cr1_exact_zero", 1.0, 0.0),
)
OUTPUT = ROOT / "docs" / "research" / "results" / "consistency_amplitude_underflow_20261115.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    records = []
    for name, cr, r in CASES:
        params = {**BASE_PARAMS, "cr": cr, "r": r}
        model = LCDM_SG(**params)
        started = time.perf_counter()
        gen_error = None
        try:
            model.gen_expansion()
        except Exception as exc:
            gen_error = f"{type(exc).__name__}: {exc}"
        gen_runtime = time.perf_counter() - started
        derived_error = None
        try:
            derived = model.derived_param
        except Exception as exc:
            derived = {}
            derived_error = f"{type(exc).__name__}: {exc}"
        at = float(derived.get("A_t", params["r"] * params["A_s"]))
        log_at = float(_log_tensor_amplitude(model.cosmo_param, {"A_t": at}))
        physical_log_product = (
            math.log(params["r"]) + math.log(params["A_s"])
            if params["r"] > 0.0 and params["A_s"] > 0.0
            else None
        )
        consistency_delta_n = None
        if (
            cr > 0.0
            and params["r"] > 0.0
            and params["A_s"] > 0.0
            and derived_error is None
        ):
            consistency_delta_n = (
                (gp.N_10 - float(derived["N_re"])) * 4.0 / 3.0
                + math.log(gp.M_Pl) * 4.0 / 3.0
                + (
                    math.log(45.0 / 2.0)
                    + physical_log_product
                    - math.log(model.rhorad_re + model.rhostiff_re)
                ) / 3.0
            )
        is_exact_zero = bool(_has_exact_zero_tensor_amplitude(model.cosmo_param))
        root_started = time.perf_counter()
        solver_error = None
        try:
            result = model.SGWB_iter(
                engine="fast",
                fallback=False,
                tol=0.99,
                z_tail=5.0,
                threads=1,
                accuracy_mode="fast",
                freq_grid="goal",
                physical_guard=False,
            )
        except Exception as exc:
            result = None
            solver_error = f"{type(exc).__name__}: {exc}"
        solver_runtime = time.perf_counter() - root_started
        spectrum = getattr(model, "log10OmegaGW", None)
        if derived_error is not None:
            failure_stage = "derived_background"
        elif solver_error is not None:
            failure_stage = "fast_solver_exception"
        elif result is None:
            failure_stage = "fast_guard_or_solver_rejection"
        else:
            failure_stage = "fast_completed"
        if model.fast_failure_reason in ("invalid_r", "invalid_cutoff"):
            s2_stage = "not_reached_pre_primitive_guard"
        elif result is not None:
            s2_stage = "primitives_and_propagation_completed"
        else:
            s2_stage = "not_identified"
        records.append({
            "case": name,
            "parameters": params,
            "exact_zero_tensor_limit": is_exact_zero,
            "A_t_stored": at,
            "V_inf": (
                None if derived.get("V_inf") is None
                else float(derived["V_inf"])
            ),
            "physical_log_r_As": physical_log_product,
            "consistency_Delta_N_from_physical_log_product": consistency_delta_n,
            "log_A_t": log_at if math.isfinite(log_at) else None,
            "log_A_t_is_finite": math.isfinite(log_at),
            "N_inf": (
                None if derived.get("N_inf") is None
                else float(derived["N_inf"])
            ),
            "gen_expansion_error": gen_error,
            "derived_background_error": derived_error,
            "gen_expansion_runtime_s": gen_runtime,
            "fast_returned_model": result is not None,
            "fast_failure_reason": model.fast_failure_reason,
            "fast_solver_error": solver_error,
            "failure_stage": failure_stage,
            "s2_stage": s2_stage,
            "fast_runtime_s": solver_runtime,
            "frequency_count": int(getattr(model, "f", np.empty(0)).size),
            "log_spectrum_finite_count": (
                None if spectrum is None
                else int(np.count_nonzero(np.isfinite(spectrum)))
            ),
        })
    output = {
        "schema_version": 1,
        "status": "complete_amplitude_underflow_classification_screen",
        "cases": records,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": _sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "LCDM_stiff_Neff.py": _sha(ROOT / "stiffgwpy_fast" / "LCDM_stiff_Neff.py"),
            "stiff_SGWB.py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [
            {
                "case": row["case"],
                "A_t_stored": row["A_t_stored"],
                "log_A_t_is_finite": row["log_A_t_is_finite"],
                "N_inf": row["N_inf"],
                "returned": row["fast_returned_model"],
                "failure": row["fast_failure_reason"],
                "failure_stage": row["failure_stage"],
                "derived_background_error": row["derived_background_error"],
            }
            for row in records
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
