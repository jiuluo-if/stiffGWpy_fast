"""Compare log-scaled subnormal maps on direct Fast and Reference grids."""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
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
R = 1e-315
OUTPUT = ROOT / "docs" / "research" / "results" / "subnormal_map_grid_convergence_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log_reference_map(model, freqs, log_amp, dn_eff, z_tail):
    model.cosmo_param["DN_eff"] = float(dn_eff)
    model.DN_eff_orig = None
    model.gen_expansion()
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            model, freqs, float(dn_eff), z_tail=z_tail, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        log_integrand = np.log(base) + log_amp
        shift = float(np.max(log_integrand))
        normalized = np.exp(log_integrand - shift)
        g2, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    h = float(model.derived_param["h"])
    log_dn = math.log(gp.Neff0 / (gp.Omega_nh2 / h**2)) + shift + math.log(g2)
    return {
        "z_tail": z_tail,
        "rtol": 1e-9,
        "runtime_s": time.perf_counter() - started,
        "log_DN_gw": log_dn,
        "DN_gw_float64": math.exp(log_dn),
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "quadrature_error_relative": None if qerr is None else float(qerr / g2),
        "interpolation_error_relative": None if ierr is None else float(ierr / g2),
        "warnings": [str(item.message) for item in caught],
    }


def main() -> int:
    physical = LCDM_SG(**{**BASE, "r": R})
    physical.gen_expansion()
    log_amp = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    rows = []
    try:
        for freq_res in (1.0, 2.0):
            grid_model = LCDM_SG(**{**BASE, "r": R})
            grid_model.gen_expansion()
            freqs = np.asarray(
                construct_fast_frequency_grid(grid_model, freq_res=freq_res, log_amplitude=log_amp),
                dtype=np.float64,
            )
            FA.goal_oriented_freqs = lambda *args, _f=freqs, **kwargs: _f.copy()
            FA.construct_fast_frequency_grid = lambda *args, _f=freqs, **kwargs: _f.copy()
            fast_model = LCDM_SG(**{**BASE, "r": R})
            started = time.perf_counter()
            fast_result = fast_model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                threads=1, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
            fast = {
                "returned": fast_result is not None,
                "root_DN_eff": None if fast_result is None else float(fast_model.cosmo_param["DN_eff"]),
                "mapped_DN_gw": None if fast_result is None else float(fast_model.DN_gw[-1]),
                "root_converged_delta": getattr(fast_model, "dn_converged_delta", None),
                "failure_reason": fast_model.fast_failure_reason,
                "failure_detail": getattr(fast_model, "fast_failure_detail", None),
                "runtime_s": time.perf_counter() - started,
            }
            refs = {}
            for z_tail in (8.0, 10.0):
                refs[f"z{int(z_tail)}"] = log_reference_map(
                    LCDM_SG(**{**BASE, "r": 1.0, "A_s": 1.0}),
                    freqs, log_amp, 0.0, z_tail,
                )
            map_change = math.expm1(refs["z10"]["log_DN_gw"] - refs["z8"]["log_DN_gw"])
            fast_ref_error = None
            if fast["mapped_DN_gw"] is not None and fast["mapped_DN_gw"] > 0.0:
                fast_ref_error = math.expm1(
                    math.log(fast["mapped_DN_gw"]) - refs["z10"]["log_DN_gw"]
                )
            rows.append({
                "freq_res": freq_res,
                "frequency_count": int(freqs.size),
                "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
                "fast": fast,
                "reference": refs,
                "reference_z8_to_z10_relative_map_change": map_change,
                "fast_z5_vs_reference_z10_map_relative_error": fast_ref_error,
            })
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    output = {
        "schema_version": 1,
        "status": "complete_cr0_nt0_subnormal_map_frequency_and_tail_screen",
        "scope": "Direct Fast self-consistent and direct unit-amplitude Cartesian Reference maps for r=1e-315 on freq_res 1/2 grids, z8/z10, with log-domain amplitude restoration.",
        "parameters": {**BASE, "r": R},
        "fast_options": {"tol": 1e-7, "z_tail": 5.0, "threads": 1},
        "reference_options": {"backend": "cartesian", "rtol": 1e-9, "workers": 1},
        "cases": rows,
        "environment": {"python": sys.version, "platform": sys.platform, "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "script": "scripts/screen_subnormal_map_grid_convergence_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "freq_res": row["freq_res"], "n": row["frequency_count"],
            "fast_map": row["fast"]["mapped_DN_gw"],
            "reference_z8": row["reference"]["z8"]["DN_gw_float64"],
            "reference_z10": row["reference"]["z10"]["DN_gw_float64"],
            "z8_to_z10": row["reference_z8_to_z10_relative_map_change"],
            "fast_vs_z10": row["fast_z5_vs_reference_z10_map_relative_error"],
        } for row in rows],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
