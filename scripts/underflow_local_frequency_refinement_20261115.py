"""Directly refine Reference/Fast frequencies around the 937-node outlier."""

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
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
TAIL_INPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_z5_z8_compare_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_local_frequency_refinement_20261115.json"
WINDOW = (886, 901)
FAST_TAIL = 8.0
REFERENCE_TAIL = 10.0
REFERENCE_RTOL = 1e-9
STEP_VARIANTS = ((0.005, False), (0.005, True), (0.0025, True), (0.00125, True), (0.000625, True))
TAIL_VARIANTS = (5.0, 8.0)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    tail_data = json.loads(TAIL_INPUT.read_text(encoding="utf-8"))
    params = dict(data["parameters"])
    dn_eff = float(data["background"]["common_DN_eff"])
    freqs = np.asarray(data["frequency_grids"]["937"]["log10_hz"], dtype=np.float64)
    lo, hi = WINDOW
    if not 0 <= lo < hi < freqs.size:
        raise ValueError("local refinement window lies outside the 937-node grid")
    candidates = 0.5 * (freqs[lo:hi] + freqs[lo + 1 : hi + 1])
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = dn_eff
    model.gen_expansion()
    original_gen = model.gen_expansion

    def locked_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = dn_eff
        return original_gen(*args, **kwargs)

    model.gen_expansion = locked_background_gen
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=FAST_TAIL,
            threads=1,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
            eval_freqs=candidates,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    fast_runtime = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    fast_freqs = np.asarray(model.f, dtype=np.float64)
    fast_logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    fast_indices = np.searchsorted(-fast_freqs, -candidates)
    if not np.array_equal(fast_freqs[fast_indices], candidates):
        raise AssertionError("Fast did not solve every local refinement frequency directly")
    fast_direct = fast_logs[fast_indices]

    unit_params = dict(params)
    unit_params.update({"r": 1.0, "A_s": 1.0, "n_t": 0.0})
    reference_model = LCDM_SG(**unit_params)
    reference_model.gen_expansion()
    log_amplitude = _log_tensor_amplitude(params, {"A_t": 0.0})
    reference_logs = np.empty(candidates.size, dtype=np.float64)
    nfev = np.empty(candidates.size, dtype=np.int64)
    reference_started = time.perf_counter()
    for index, freq in enumerate(candidates):
        mode = REF.solve_reference_mode(
            reference_model,
            float(freq),
            dn_eff,
            z_tail=REFERENCE_TAIL,
            rtol=REFERENCE_RTOL,
            backend="cartesian",
        )
        base = 3.0 * mode["Opgw_today"] if mode["used_tail"] else (
            mode["Ogw_today"] - mode["Oj_today"]
        )
        log_power = log_amplitude + float(params["n_t"]) * (
            math.log(10.0) * float(freq) - math.log(gp.f_piv)
        )
        reference_logs[index] = (math.log(base) + log_power) / math.log(10.0)
        nfev[index] = int(mode["n_steps"])
    reference_runtime = time.perf_counter() - reference_started

    tail_worst_reference_checks = []
    tail_reference_started = time.perf_counter()
    for count in (235, 469, 937):
        tail_row = tail_data["grids"][str(count)]
        grid_freqs = np.asarray(
            data["frequency_grids"][str(count)]["log10_hz"], dtype=np.float64
        )
        z8_fast_logs = np.asarray(
            tail_row["z_tail_8"]["log10_OmegaGW"], dtype=np.float64
        )
        z8_reference_logs = np.asarray(
            data["reference_fixed_background"][str(count)]["log10_OmegaGW"],
            dtype=np.float64,
        )
        local_errors = np.abs(z8_fast_logs - z8_reference_logs)
        index = int(np.argmax(local_errors))
        freq = float(grid_freqs[index])
        mode = REF.solve_reference_mode(
            reference_model,
            freq,
            dn_eff,
            z_tail=REFERENCE_TAIL,
            rtol=REFERENCE_RTOL,
            backend="cartesian",
        )
        base = 3.0 * mode["Opgw_today"] if mode["used_tail"] else (
            mode["Ogw_today"] - mode["Oj_today"]
        )
        log_power = _log_tensor_amplitude(params, {"A_t": 0.0}) + float(
            params["n_t"]
        ) * (math.log(10.0) * freq - math.log(gp.f_piv))
        reference_z10_log = (math.log(base) + log_power) / math.log(10.0)
        tail_worst_reference_checks.append({
            "grid_count": count,
            "frequency_index": index,
            "log10_hz": freq,
            "fast_z8_log10_OmegaGW": float(z8_fast_logs[index]),
            "reference_z8_log10_OmegaGW": float(z8_reference_logs[index]),
            "reference_z10_log10_OmegaGW": reference_z10_log,
            "fast_z8_abs_error_vs_reference_z10_dex": abs(
                float(z8_fast_logs[index]) - reference_z10_log
            ),
            "reference_z8_to_z10_abs_delta_dex": abs(
                float(z8_reference_logs[index]) - reference_z10_log
            ),
            "reference_nfev": int(mode["n_steps"]),
            "reference_handoff_eps": mode["eps_handoff"],
        })
    tail_reference_runtime = time.perf_counter() - tail_reference_started

    grid_logs = np.asarray(data["fast_fixed_background"]["937"]["log10_OmegaGW"])
    ref_grid_logs = np.asarray(data["reference_fixed_background"]["937"]["log10_OmegaGW"])
    outlier_index = 893
    outlier_frequency = float(freqs[outlier_index])
    reference_outlier = float(ref_grid_logs[outlier_index])
    step_refinements = []
    for tail_z in TAIL_VARIANTS:
        for h, sigma_exact in STEP_VARIANTS:
            step_model = LCDM_SG(**params)
            step_model.cosmo_param["DN_eff"] = dn_eff
            step_model.gen_expansion()
            original_step_gen = step_model.gen_expansion

            def locked_step_background_gen(*args, **kwargs):
                step_model.cosmo_param["DN_eff"] = dn_eff
                return original_step_gen(*args, **kwargs)

            step_model.gen_expansion = locked_step_background_gen
            FA.goal_oriented_freqs = fixed_grid
            FA.construct_fast_frequency_grid = fixed_grid
            step_started = time.perf_counter()
            try:
                step_result = step_model.SGWB_iter(
                    engine="fast",
                    fallback=False,
                    tol=0.99,
                    z_tail=tail_z,
                    threads=1,
                    accuracy_mode="fast",
                    freq_grid="goal",
                    physical_guard=False,
                    h=h,
                    sigma_exact=sigma_exact,
                )
            finally:
                FA.goal_oriented_freqs = original_goal
                FA.construct_fast_frequency_grid = original_construct
            step_runtime = time.perf_counter() - step_started
            if step_result is None:
                raise RuntimeError(
                    f"Fast step refinement failed at h={h}: {step_model.fast_failure_reason}"
                )
            step_index = int(np.flatnonzero(step_model.f == outlier_frequency)[0])
            step_log = float(step_model.log10OmegaGW[step_index])
            step_refinements.append({
                "z_tail": tail_z,
                "h": h,
                "sigma_exact": sigma_exact,
                "runtime_s": step_runtime,
                "mapped_DN_gw": float(step_model.DN_gw[-1]),
                "log10_OmegaGW": step_log,
                "error_vs_direct_reference_dex": abs(step_log - reference_outlier),
            })

    fast_interp = PchipInterpolator(freqs[::-1], grid_logs[::-1])(candidates)
    reference_interp = PchipInterpolator(freqs[::-1], ref_grid_logs[::-1])(candidates)
    error = fast_direct - reference_logs
    output = {
        "schema_version": 1,
        "status": "complete_direct_local_refinement_1873_equivalent_nodes",
        "parameters": params,
        "frozen_DN_eff": dn_eff,
        "base_grid_count": int(freqs.size),
        "base_grid_sha256_f64le": data["frequency_grids"]["937"]["sha256_f64le"],
        "window_937_indices_inclusive": [lo, hi],
        "candidate_count": int(candidates.size),
        "candidate_log10_hz": candidates.tolist(),
        "fast_tail_z": FAST_TAIL,
        "reference_tail_z": REFERENCE_TAIL,
        "reference_rtol": REFERENCE_RTOL,
        "fast_runtime_s": fast_runtime,
        "reference_runtime_s": reference_runtime,
        "reference_z10_tail_check_runtime_s": tail_reference_runtime,
        "z10_reference_at_each_grid_z8_worst_node": tail_worst_reference_checks,
        "reference_nfev": nfev.tolist(),
        "direct_fast_log10_OmegaGW": fast_direct.tolist(),
        "direct_reference_log10_OmegaGW": reference_logs.tolist(),
        "fast_minus_reference_dex": error.tolist(),
        "abs_error_summary_dex": {
            "max": float(np.max(np.abs(error))),
            "p95": float(np.quantile(np.abs(error), 0.95)),
        },
        "pchip_prediction_error_against_direct_reference_dex": {
            "max": float(np.max(np.abs(reference_interp - reference_logs))),
            "p95": float(np.quantile(np.abs(reference_interp - reference_logs), 0.95)),
        },
        "pchip_prediction_error_against_direct_fast_dex": {
            "max": float(np.max(np.abs(fast_interp - fast_direct))),
            "p95": float(np.quantile(np.abs(fast_interp - fast_direct), 0.95)),
        },
        "937_node_outlier": {
            "index": outlier_index,
            "log10_hz": outlier_frequency,
            "reference_log10_OmegaGW": reference_outlier,
            "step_refinements": step_refinements,
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "grid_input_json": _sha(INPUT),
            "tail_input_json": _sha(TAIL_INPUT),
            "diagnostic_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "candidate_count": int(candidates.size),
        "direct_max_p95_error_dex": [output["abs_error_summary_dex"]["max"], output["abs_error_summary_dex"]["p95"]],
        "reference_interpolation_max_p95_dex": [output["pchip_prediction_error_against_direct_reference_dex"]["max"], output["pchip_prediction_error_against_direct_reference_dex"]["p95"]],
        "fast_runtime_s": fast_runtime,
        "reference_runtime_s": reference_runtime,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
