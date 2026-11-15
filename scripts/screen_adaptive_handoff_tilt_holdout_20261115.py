"""Screen the z_today<5 handoff mask on a held-out blue tilt."""

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

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, construct_fast_frequency_grid  # noqa: E402
from stiffgwpy_fast.reference import _Background, f_hor_abs, integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_tilt_holdout_20261115.json"
THRESHOLDS = (5.0, 6.0, 7.0, 8.0)
PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 5e-324,
    "n_t": 17.30,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}


def run_fast(freqs: np.ndarray, z_tail: float, frozen_dn: float | None = None) -> dict[str, object]:
    model = LCDM_SG(**PARAMS)
    if frozen_dn is not None:
        model.cosmo_param["DN_eff"] = frozen_dn
        model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = frozen_dn
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    if frozen_dn is not None:
        model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        try:
            result = model.SGWB_iter(
                engine="fast", fallback=False,
                tol=1e-7 if frozen_dn is None else 0.99,
                z_tail=z_tail, threads=1, accuracy_mode="fast",
                freq_grid="goal", physical_guard=False,
            )
        finally:
            FA.goal_oriented_freqs = original_goal
            FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast z_tail={z_tail} failed: {model.fast_failure_reason}")
    return {
        "root": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "Ogw": np.asarray(model.Ogw_today, dtype=np.float64),
        "Oj": np.asarray(model.Oj_today, dtype=np.float64),
        "log_spectrum": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "runtime_s": elapsed,
        "runtime_warnings": [str(item.message) for item in caught],
    }


def main() -> int:
    physical = LCDM_SG(**PARAMS)
    physical.gen_expansion()
    log_amp = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    freqs = np.asarray(construct_fast_frequency_grid(physical, log_amplitude=log_amp), dtype=np.float64)
    root_solution = run_fast(freqs, 5.0)
    root = float(root_solution["root"])
    fast5 = run_fast(freqs, 5.0, frozen_dn=root)
    fast8 = run_fast(freqs, 8.0, frozen_dn=root)

    unit_params = {**PARAMS, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0}
    reference_model = LCDM_SG(**unit_params)
    reference_model.gen_expansion()
    with warnings.catch_warnings(record=True) as caught_ref:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            reference_model, freqs, root, z_tail=8.0, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        unit_base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        if np.any(~np.isfinite(unit_base)) or np.any(unit_base <= 0.0):
            raise FloatingPointError("Reference holdout base spectrum is not positive finite")
        log_power = log_amp + PARAMS["n_t"] * (
            math.log(10.0) * freqs - math.log(gp.f_piv)
        )
        log_reference = (np.log(unit_base) + log_power) / gp.ln10
        shift = float(np.max(np.log(unit_base) + log_power))
        normalized = np.exp(np.log(unit_base) + log_power - shift)
        g2_scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    log_map_scale = math.log(gp.Neff0 / (gp.Omega_nh2 / float(reference_model.derived_param["h"]) ** 2))
    # ``integrate_spectrum`` already returns the d(ln f) factor ln(10).
    reference_map = math.exp(shift + math.log(g2_scaled) + log_map_scale)

    bg = _Background(reference_model, root)
    f_hor_today = f_hor_abs(bg, bg.N_inf)
    z_today = (freqs - f_hor_today) * gp.ln10
    adaptive_candidates = []
    for threshold in THRESHOLDS:
        use_z8 = z_today < threshold
        hybrid_log = np.where(use_z8, fast8["log_spectrum"], fast5["log_spectrum"])
        error = hybrid_log - log_reference
        hybrid_integrand = np.where(use_z8, fast8["Ogw"] - fast8["Oj"], fast5["Ogw"] - fast5["Oj"])
        hybrid_map = gp.Neff0 * gp.ln10 * float(np.sum(FS._pchip_integrals_vectorized(freqs, hybrid_integrand))) / (
            gp.Omega_nh2 / float(reference_model.derived_param["h"]) ** 2
        )
        adaptive_candidates.append({
            "z_today_threshold": threshold,
            "z8_selected_mode_count": int(np.count_nonzero(use_z8)),
            "z8_selected_mode_fraction": float(np.mean(use_z8)),
            "max_abs_spectrum_error_vs_reference_z8_dex": float(np.max(np.abs(error))),
            "p95_abs_spectrum_error_vs_reference_z8_dex": float(np.quantile(np.abs(error), 0.95)),
            "max_error_log10_hz": float(freqs[np.argmax(np.abs(error))]),
            "mapped_DN_gw": float(hybrid_map),
            "relative_map_error_vs_reference": float(hybrid_map / reference_map - 1.0),
            "relative_fixed_root_residual": float(hybrid_map / root - 1.0),
        })
    adaptive_candidate = adaptive_candidates[0]
    output = {
        "schema_version": 1,
        "status": "complete_adaptive_handoff_tilt_holdout_screen_not_production",
        "scope": "Independent n_t holdout at minimum positive subnormal r; directly solve Fast z5 self-consistent root, fixed-root z5/z8, and Reference z8 on one native Fast goal grid.",
        "parameters": PARAMS,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "frozen_DN_eff": root,
        "frequency_log10_hz": freqs.tolist(),
        "z_today_range": [float(np.min(z_today)), float(np.max(z_today))],
        "z8_selected_mode_count": adaptive_candidate["z8_selected_mode_count"],
        "z8_selected_mode_fraction": adaptive_candidate["z8_selected_mode_fraction"],
        "root_solution": {
            "root_DN_eff": root,
            "mapped_DN_gw": root_solution["mapped_DN_gw"],
            "outer_map_relative_residual": float(root_solution["mapped_DN_gw"] / root - 1.0),
            "runtime_s": root_solution["runtime_s"],
            "runtime_warnings": root_solution["runtime_warnings"],
        },
        "fixed_background_runs": {
            "z5_runtime_s_not_paired": fast5["runtime_s"],
            "z8_runtime_s_not_paired": fast8["runtime_s"],
            "z5_warnings": fast5["runtime_warnings"],
            "z8_warnings": fast8["runtime_warnings"],
        },
        "reference": {
            "z_tail": 8.0,
            "rtol": 1e-9,
            "mapped_DN_gw": reference_map,
            "map_relative_residual": float(reference_map / root - 1.0),
            "used_tail_count": int(np.count_nonzero(used_tail)),
            "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(normalized == 0.0)),
            "quadrature_error_relative": None if qerr is None else float(qerr / g2_scaled),
            "interpolation_error_relative": None if ierr is None else float(ierr / g2_scaled),
            "runtime_warnings": [str(item.message) for item in caught_ref],
        },
        "adaptive_candidate": {
            **adaptive_candidate,
            "criterion": "use z8 row if z_today is below threshold; otherwise use z5 row",
            "runtime_status": "not measured; z5/z8 rows came from separate full direct runs",
        },
        "adaptive_candidates": adaptive_candidates,
        "interpretation": "This is an independent tilt holdout for an error-selection criterion, not an implementation or performance result.",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "diagnostic_script": "scripts/screen_adaptive_handoff_tilt_holdout_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": output["frequency_count"],
        "root": root,
        "z8_fraction": output["z8_selected_mode_fraction"],
        "thresholds": [{"threshold": row["z_today_threshold"], "z8_fraction": row["z8_selected_mode_fraction"], "max_error": row["max_abs_spectrum_error_vs_reference_z8_dex"], "map_error": row["relative_map_error_vs_reference"]} for row in adaptive_candidates],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
