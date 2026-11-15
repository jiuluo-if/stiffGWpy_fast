"""Distinguish finite log amplitude from an unrepresentable self-consistency integral."""

from __future__ import annotations

import hashlib
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

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
R_MIN_SUBNORMAL = float(np.nextafter(0.0, 1.0))
OUTPUT = ROOT / "docs" / "research" / "results" / "dn_gw_integral_underflow_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fast_case(r: float, freqs: np.ndarray) -> dict[str, object]:
    model = LCDM_SG(**{**BASE, "r": r})
    log_amp = _log_tensor_amplitude(model.cosmo_param, model.derived_param)
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    return {
        "r": r,
        "stored_A_t": float(model.derived_param["A_t"]),
        "log_A_t": float(log_amp) if math.isfinite(log_amp) else None,
        "fast_returned_model": result is not None,
        "fast_failure_reason": model.fast_failure_reason,
        "fast_failure_detail": getattr(model, "fast_failure_detail", None),
        "root_DN_eff": None if result is None else float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": None if result is None else float(model.DN_gw[-1]),
        "dn_converged_delta": getattr(model, "dn_converged_delta", None),
        "runtime_s": time.perf_counter() - started,
        "frequency_count": int(np.asarray(getattr(model, "f", [])).size),
    }


def main() -> int:
    min_params = {**BASE, "r": R_MIN_SUBNORMAL}
    physical = LCDM_SG(**min_params)
    physical.gen_expansion()
    log_amp = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    freqs = np.asarray(construct_fast_frequency_grid(physical, log_amplitude=log_amp), dtype=np.float64)

    below_root = fast_case(1e-315, freqs)
    two_subnormal_units = fast_case(2.0 * R_MIN_SUBNORMAL, freqs)
    three_subnormal_units = fast_case(3.0 * R_MIN_SUBNORMAL, freqs)
    underflow_root = fast_case(R_MIN_SUBNORMAL, freqs)
    exact_zero = fast_case(0.0, freqs)

    unit_params = {**BASE, "r": 1.0, "A_s": 1.0, "n_t": 0.0}
    unit_model = LCDM_SG(**unit_params)
    unit_model.gen_expansion()
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            unit_model, freqs, 0.0, z_tail=8.0, rtol=1e-9,
            workers=1, backend="cartesian",
        )
        base_spectrum = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        log_integrand = np.log(base_spectrum) + log_amp
        shift = float(np.max(log_integrand))
        scaled = np.exp(log_integrand - shift)
        g2_scaled, qerr, ierr = integrate_spectrum(freqs, scaled, np.zeros_like(scaled))
    reference_runtime = time.perf_counter() - started
    h = float(unit_model.derived_param["h"])
    log_map = math.log(gp.Neff0 / (gp.Omega_nh2 / (h * h))) + shift + math.log(g2_scaled)
    min_subnormal_log = math.log(R_MIN_SUBNORMAL)
    half_min_subnormal_log = min_subnormal_log - math.log(2.0)
    root_fraction_min_subnormal = math.exp(log_map - min_subnormal_log)
    for case in (below_root, two_subnormal_units, three_subnormal_units, underflow_root):
        case_log_map = log_map + case["log_A_t"] - log_amp
        case["reference_scaled_log_DN_gw"] = case_log_map
        case["reference_scaled_DN_gw_float64"] = math.exp(case_log_map)
        case["fast_vs_reference_scaled_relative_error"] = (
            None if case["mapped_DN_gw"] is None or case["mapped_DN_gw"] == 0.0
            else math.expm1(math.log(case["mapped_DN_gw"]) - case_log_map)
        )
        if case["mapped_DN_gw"] is not None and case["mapped_DN_gw"] > 0.0:
            half_ulp_ratio = 0.5 * (R_MIN_SUBNORMAL / case["mapped_DN_gw"])
            case["float64_root_rounding_relative_bound"] = (
                half_ulp_ratio / (1.0 - half_ulp_ratio)
            )
        else:
            case["float64_root_rounding_relative_bound"] = None
    output = {
        "schema_version": 1,
        "status": "complete_log_amplitude_vs_integral_underflow_classification",
        "scope": "cr=0, n_t=0. Direct Fast around the subnormal-output threshold, compared with a unit-transfer Reference log integral on the identical native 235-node grid; exact r=0 remains separate.",
        "parameters_except_r": BASE,
        "r_min_subnormal": R_MIN_SUBNORMAL,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "diagnostic_script": "scripts/analyze_dn_gw_integral_underflow_20261115.py",
            }.items()
        },
        "fast_cases": {
            "nearby_representable_root_r_1e-315": below_root,
            "r_two_input_subnormal_units": two_subnormal_units,
            "r_three_input_subnormal_units": three_subnormal_units,
            "minimum_positive_r": underflow_root,
            "exact_zero_r": exact_zero,
        },
        "reference_unit_transfer_at_minimum_positive_r": {
            "input_log_A_t_finite": bool(math.isfinite(log_amp)),
            "log10_A_t": log_amp / math.log(10.0),
            "reference_z_tail": 8.0,
            "reference_rtol": 1e-9,
            "runtime_s": reference_runtime,
            "mapped_log_DN_gw": log_map,
            "mapped_log10_DN_gw": log_map / math.log(10.0),
            "minimum_subnormal_log": min_subnormal_log,
            "half_minimum_subnormal_log": half_min_subnormal_log,
            "log_root_below_min_subnormal_by": min_subnormal_log - log_map,
            "root_fraction_of_minimum_subnormal": root_fraction_min_subnormal,
            "rounds_to_zero_under_nearest_ties_even": bool(log_map < half_min_subnormal_log),
            "mapped_DN_gw_float64": math.exp(log_map),
            "used_tail_count": int(np.count_nonzero(used_tail)),
            "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(scaled == 0.0)),
            "quadrature_error_relative": None if qerr is None else float(qerr / g2_scaled),
            "interpolation_error_relative": None if ierr is None else float(ierr / g2_scaled),
            "warnings": [str(item.message) for item in caught],
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "r1e-315": below_root,
        "r2sub": two_subnormal_units,
        "r3sub": three_subnormal_units,
        "rmin": underflow_root,
        "r0": exact_zero,
        "reference_log10_root": output["reference_unit_transfer_at_minimum_positive_r"]["mapped_log10_DN_gw"],
        "root_below_min_subnormal_by_ln": output["reference_unit_transfer_at_minimum_positive_r"]["log_root_below_min_subnormal_by"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
