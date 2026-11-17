"""English description: Compare physical-power and unit-transfer log scaling on one fixed free-tilt background. 中文说明：同背景、同频网格下对照普通功率传播与 log 振幅缩放，并区分精确零。"""

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

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude, _requires_log_tensor_power  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_current_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "log_amplitude_scale_equivalence_20261117.json"
FROZEN_DN = 0.027412228089632867
N_T = 1.6
Z_TAIL = 8.0
THREADS = 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def error_stats(a: np.ndarray, b: np.ndarray) -> dict[str, float | int]:
    delta = np.abs(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64))
    finite = np.isfinite(delta)
    if not finite.any():
        raise FloatingPointError("No finite spectrum nodes are comparable")
    values = delta[finite]
    return {
        "finite_count": int(values.size),
        "nonfinite_count": int(delta.size - values.size),
        "max_abs_dex": float(np.max(values)),
        "p95_abs_dex": float(np.quantile(values, 0.95)),
    }


def array_hash(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def run_case(
    base_params: dict[str, float], freqs: np.ndarray, r: float, *, force_log_mode: bool
) -> dict[str, object]:
    params = dict(base_params)
    params["r"] = float(r)
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = FROZEN_DN
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    original_requires = FA._requires_log_tensor_power

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = FROZEN_DN
        return original_gen(*args, **kwargs)

    if force_log_mode:
        def force_log_requirement(*args, **kwargs):
            return True
        FA._requires_log_tensor_power = force_log_requirement
    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    model.gen_expansion = fixed_background_gen
    start = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=Z_TAIL,
            threads=THREADS,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
        FA._requires_log_tensor_power = original_requires
    runtime = time.perf_counter() - start

    derived = model.derived_param
    log_amplitude = _log_tensor_amplitude(model.cosmo_param, derived)
    requires_log = (
        math.isfinite(log_amplitude)
        and _requires_log_tensor_power(log_amplitude, derived["nt"], freqs, gp.f_piv)
    )
    linear = np.asarray(model.Ogw_today, dtype=np.float64) - np.asarray(model.Oj_today, dtype=np.float64)
    logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    mapped = float(model.DN_gw[-1]) if hasattr(model, "DN_gw") else float("nan")
    return {
        "r": float(r),
        "A_t_stored": float(derived["A_t"]),
        "log_amplitude": log_amplitude if math.isfinite(log_amplitude) else None,
        "log_mode_expected_without_forcing": bool(
            derived["A_t"] == 0.0 or not math.isfinite(float(derived["A_t"])) or requires_log
        ),
        "force_log_mode": force_log_mode,
        "requires_log_power_at_fixed_grid": bool(requires_log),
        "result_success": result is model,
        "failure_reason": model.fast_failure_reason,
        "runtime_s": runtime,
        "fixed_DN_eff": FROZEN_DN,
        "mapped_DN_gw": mapped if math.isfinite(mapped) else None,
        "map_residual_relative": mapped / FROZEN_DN - 1.0 if math.isfinite(mapped) else None,
        "frequency_count": int(model.f.size) if hasattr(model, "f") else 0,
        "frequency_sha256_f64le": digest(np.asarray(model.f, dtype=np.float64)) if hasattr(model, "f") else None,
        "final_background_grid_sha256": {
            "Nv": array_hash(np.asarray(model.Nv, dtype=np.float64)),
            "f_hor": array_hash(np.asarray(model.f_hor, dtype=np.float64)),
        },
        "linear_spectrum_zero_count": int(np.count_nonzero(linear == 0.0)),
        "linear_spectrum_positive_count": int(np.count_nonzero(linear > 0.0)),
        "log_spectrum_finite_count": int(np.count_nonzero(np.isfinite(logs))),
        "log_spectrum_negative_infinity_count": int(np.count_nonzero(np.isneginf(logs))),
        "log10_OmegaGW": [float(value) if math.isfinite(float(value)) else None for value in logs],
    }


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source.get("status") != "complete_fixed_background_direct_nested_grid_screen":
        raise RuntimeError("Current-source nested-grid artifact is incomplete")
    for rel, expected in source["source_sha256"].items():
        if sha256(ROOT / rel) != expected:
            raise RuntimeError(f"Input source hash differs for {rel}")
    frequencies = np.asarray(source["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    if frequencies.size != 235 or digest(frequencies) != source["frequency_grids"]["235"]["sha256_f64le"]:
        raise RuntimeError("Exact 235-node grid hash mismatch")

    base_params = {key: float(value) for key, value in source["parameters"].items()}
    base_params["n_t"] = N_T
    cases = {
        "representable_standard": run_case(base_params, frequencies, 1e-300, force_log_mode=False),
        "representable_forced_log": run_case(base_params, frequencies, 1e-300, force_log_mode=True),
        "nonzero_underflow_log": run_case(base_params, frequencies, 1e-316, force_log_mode=False),
        "exact_zero": run_case(base_params, frequencies, 0.0, force_log_mode=False),
    }
    standard = cases["representable_standard"]
    forced_log = cases["representable_forced_log"]
    underflow = cases["nonzero_underflow_log"]
    exact_zero = cases["exact_zero"]
    if not all(cases[key]["result_success"] for key in cases):
        raise RuntimeError(f"A Fast scaling branch failed: {[(k,v['failure_reason']) for k,v in cases.items()]}")
    if any(
        cases[key]["final_background_grid_sha256"] != standard["final_background_grid_sha256"]
        for key in cases
    ):
        raise RuntimeError("Amplitude-scale cases used different final background grids")
    if not all(cases[key]["frequency_sha256_f64le"] == digest(frequencies) for key in cases):
        raise RuntimeError("Scaling comparison did not use the exact common grid")
    if standard["force_log_mode"] or standard["requires_log_power_at_fixed_grid"]:
        raise RuntimeError("Control case unexpectedly selected log-domain power")
    if not forced_log["force_log_mode"] or not underflow["log_mode_expected_without_forcing"]:
        raise RuntimeError("Expected log-amplitude case did not select the intended branch")
    if exact_zero["log_spectrum_finite_count"] != 0 or exact_zero["log_spectrum_negative_infinity_count"] != frequencies.size:
        raise RuntimeError("Exact zero tensor limit was conflated with nonzero underflow")

    standard_logs = np.asarray(standard["log10_OmegaGW"], dtype=np.float64)
    forced_logs = np.asarray(forced_log["log10_OmegaGW"], dtype=np.float64)
    underflow_logs = np.asarray(underflow["log10_OmegaGW"], dtype=np.float64)
    expected_log_ratio = (float(underflow["log_amplitude"]) - float(standard["log_amplitude"])) / math.log(10.0)
    output = {
        "schema_version": 1,
        "status": "complete_fixed_background_log_amplitude_scale_equivalence",
        "scope": "cr=0 only: compare standard and forced unit-transfer log propagation at identical nonzero parameters, background, frequency grid and threads; then test separate nonzero product underflow and exact-zero limits. No cr>0 background dependency is bypassed.",
        "input_artifact": INPUT.relative_to(ROOT).as_posix(),
        "n_t": N_T,
        "frozen_DN_eff": FROZEN_DN,
        "z_tail": Z_TAIL,
        "threads": THREADS,
        "frequency_count": int(frequencies.size),
        "frequency_sha256_f64le": digest(frequencies),
        "cases": cases,
        "same_parameter_standard_vs_forced_log_spectrum_dex": error_stats(standard_logs, forced_logs),
        "same_parameter_standard_vs_forced_log_map_relative_difference": (
            float(standard["mapped_DN_gw"]) / float(forced_log["mapped_DN_gw"]) - 1.0
        ),
        "nonzero_underflow_vs_standard_log_ratio_expected_dex": expected_log_ratio,
        "nonzero_underflow_vs_standard_spectrum_shape_error_dex": error_stats(
            underflow_logs, standard_logs + expected_log_ratio
        ),
        "nonzero_underflow_vs_standard_map_ratio_relative_error": (
            float(underflow["mapped_DN_gw"])
            / (float(standard["mapped_DN_gw"]) * math.pow(10.0, expected_log_ratio))
            - 1.0
        ),
        "source_sha256": {
            **source["source_sha256"],
            "input_artifact": sha256(INPUT),
            "comparison_runner": sha256(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    }
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "standard_log_A_t": standard["log_amplitude"],
        "standard_at_subnormal": standard["A_t_stored"],
        "underflow_A_t_stored": underflow["A_t_stored"],
        "same_parameter_spectrum_delta": output["same_parameter_standard_vs_forced_log_spectrum_dex"],
        "same_parameter_map_delta": output["same_parameter_standard_vs_forced_log_map_relative_difference"],
        "underflow_shape_error": output["nonzero_underflow_vs_standard_spectrum_shape_error_dex"],
        "underflow_map_ratio_error": output["nonzero_underflow_vs_standard_map_ratio_relative_error"],
        "linear_zero_but_log_finite": {
            key: [cases[key]["linear_spectrum_zero_count"], cases[key]["log_spectrum_finite_count"]]
            for key in ("representable_standard", "representable_forced_log", "nonzero_underflow_log", "exact_zero")
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
