"""Screen Fast at the normal/subnormal/exact-zero tensor-amplitude boundaries."""

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

from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
R_BOUNDARY = float(np.finfo(np.float64).tiny / BASE["A_s"])
AMPLITUDES = (
    R_BOUNDARY * (1.0 + 1e-4),
    R_BOUNDARY,
    R_BOUNDARY * (1.0 - 1e-4),
    1e-300,
    1e-307,
    1e-315,
    5e-324,
    0.0,
)
OUTPUT = ROOT / "docs" / "research" / "results" / "amplitude_representation_boundary_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    cases = []
    for r in AMPLITUDES:
        params = {**BASE, "r": float(r)}
        model = LCDM_SG(**params)
        started = time.perf_counter()
        error = None
        try:
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                threads=1, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
        except Exception as exc:
            result = None
            error = f"{type(exc).__name__}: {exc}"
        runtime = time.perf_counter() - started
        derived = model.derived_param
        at = float(derived["A_t"])
        log_at = float(_log_tensor_amplitude(model.cosmo_param, derived))
        spectrum = getattr(model, "log10OmegaGW", None)
        finite_count = None if spectrum is None else int(np.count_nonzero(np.isfinite(spectrum)))
        linear = getattr(model, "Ogw_today", None)
        if linear is not None:
            linear_signal = np.asarray(linear) - np.asarray(model.Oj_today)
            linear_zero_count = int(np.count_nonzero(linear_signal == 0.0))
            linear_positive_subnormal_count = int(np.count_nonzero(
                (linear_signal > 0.0) & (linear_signal < np.finfo(float).tiny)
            ))
        else:
            linear_zero_count = None
            linear_positive_subnormal_count = None
        cases.append({
            "r": float(r),
            "parameters": params,
            "stored_A_t": at,
            "A_t_class": (
                "exact_zero_input" if r == 0.0
                else "rounded_zero_nonzero_product" if at == 0.0
                else "normal" if at >= np.finfo(float).tiny
                else "subnormal"
            ),
            "log_A_t": log_at if math.isfinite(log_at) else None,
            "exact_zero_tensor_limit": bool(r == 0.0),
            "fast_returned_model": result is not None,
            "failure_reason": model.fast_failure_reason,
            "error": error,
            "root_DN_eff": float(model.cosmo_param["DN_eff"]) if result is not None else None,
            "mapped_DN_gw": float(model.DN_gw[-1]) if result is not None else None,
            "root_converged_delta": getattr(model, "dn_converged_delta", None),
            "frequency_count": int(np.asarray(getattr(model, "f", [])).size),
            "frequency_sha256_f64le": hashlib.sha256(np.asarray(getattr(model, "f", []), dtype="<f8").tobytes()).hexdigest(),
            "finite_log_spectrum_count": finite_count,
            "linear_exact_zero_count": linear_zero_count,
            "linear_positive_subnormal_count": linear_positive_subnormal_count,
            "runtime_s": runtime,
            "log10_OmegaGW": (
                None if spectrum is None else [
                    float(value) if np.isfinite(value) else None
                    for value in np.asarray(spectrum, dtype=np.float64)
                ]
            ),
        })
    base_row = next(row for row in cases if row["r"] == R_BOUNDARY * (1.0 + 1e-4))
    root_per_r = None
    if base_row["root_DN_eff"] not in (None, 0.0):
        root_per_r = float(base_row["root_DN_eff"]) / float(base_row["r"])
    for row in cases:
        row["root_per_r_relative_to_normal_side"] = (
            None if root_per_r is None or row["root_DN_eff"] is None or row["r"] == 0.0
            else float(row["root_DN_eff"]) / float(row["r"]) / root_per_r - 1.0
        )
    output = {
        "schema_version": 1,
        "status": "complete_cr0_normal_subnormal_underflow_exact_zero_boundary_screen",
        "scope": "Direct Fast self-consistent roots at n_t=0 across the stored A_t normal/subnormal boundary, deeper subnormal products, rounded-to-zero nonzero product and exact r=0.",
        "parameters_except_r": BASE,
        "r_normal_subnormal_boundary": R_BOUNDARY,
        "amplitudes": list(AMPLITUDES),
        "fast_options": {"tol": 1e-7, "z_tail": 5.0, "threads": 1, "physical_guard": True},
        "cases": cases,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "diagnostic_script": "scripts/screen_amplitude_representation_boundary_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "boundary_r": R_BOUNDARY,
        "cases": [{
            "r": row["r"], "class": row["A_t_class"],
            "A_t": row["stored_A_t"], "root": row["root_DN_eff"],
            "root_per_r_relative": row["root_per_r_relative_to_normal_side"],
            "failure": row["failure_reason"],
            "finite_log_nodes": row["finite_log_spectrum_count"],
        } for row in cases],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
