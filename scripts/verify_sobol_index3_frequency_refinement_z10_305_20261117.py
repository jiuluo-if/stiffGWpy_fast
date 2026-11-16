"""Directly extend the Sobol index 3 z10 grid from 153 to 305 nodes."""

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
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index3_frequency_refinement_z10_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index3_frequency_refinement_z10_305_20261117.json"
Z_TAIL = 10.0
REFERENCE_RTOL = 1e-9
ROOT_RESIDUAL_TOL = 2e-4
MAX_ROOT_ITER = 5


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def raw_log_spectrum(result: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    signal = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
    valid = np.isfinite(signal) & (signal > 0.0)
    logs = np.full(signal.shape, np.nan, dtype=np.float64)
    logs[valid] = np.log10(signal[valid])
    return logs, valid


def solve_reference(model: LCDM_SG, root: float, freqs: np.ndarray) -> tuple[dict[str, object], float, list[str]]:
    caught: list[warnings.WarningMessage]
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model, dn_eff=root, freq_subset=freqs, z_tail=Z_TAIL,
            rtol=REFERENCE_RTOL, self_consistent=False, physical_guard=True,
        )
    runtime = time.perf_counter() - start
    return result, runtime, [f"{item.category.__name__}: {item.message}" for item in caught]


def main() -> int:
    previous = json.loads(INPUT.read_text(encoding="utf-8"))
    if previous["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("153-node Reference artifact source SHA differs from current Reference")
    if previous["source_sha256"]["fast_sgwb.py"] != sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"):
        raise RuntimeError("153-node Reference artifact Fast source SHA differs from current Fast")
    params = {key: float(value) for key, value in previous["parameters"].items()}
    root = float(previous["fixed_root"])
    model = LCDM_SG(**params)
    fast_start = time.perf_counter()
    fast_result = model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    fast_runtime = time.perf_counter() - fast_start
    if fast_result is not model:
        raise RuntimeError(f"Sobol index 3 current Fast failed: {model.fast_failure_reason}")
    freqs77 = np.asarray(model.f, dtype=np.float64)
    if digest(freqs77) != previous["frequency_grids"]["base_sha256_f64le"]:
        raise RuntimeError("Replayed Fast 77-node grid differs from z10 refinement artifact")
    if not math.isclose(float(model.cosmo_param["DN_eff"]), root, rel_tol=2e-12, abs_tol=0.0):
        raise RuntimeError("Replayed Fast root differs from z10 refinement artifact")
    mid77 = (freqs77[:-1] + freqs77[1:]) / 2.0
    freqs153 = np.empty(2 * freqs77.size - 1, dtype=np.float64)
    freqs153[::2], freqs153[1::2] = freqs77, mid77
    if digest(freqs153) != previous["frequency_grids"]["refined_sha256_f64le"]:
        raise RuntimeError("Reconstructed 153-node grid differs from direct z10 artifact")
    mid153 = (freqs153[:-1] + freqs153[1:]) / 2.0
    freqs305 = np.empty(2 * freqs153.size - 1, dtype=np.float64)
    freqs305[::2], freqs305[1::2] = freqs153, mid153
    direction = 1.0 if freqs153[-1] > freqs153[0] else -1.0
    if not np.all(direction * np.diff(freqs305) > 0.0):
        raise RuntimeError("305-node nested grid lost native ordering")

    base_logs = np.asarray(previous["spectra"]["refined_z10_log10_OmegaGW_at_fast_root"], dtype=np.float64)
    if base_logs.size != freqs153.size or not np.all(np.isfinite(base_logs)):
        raise RuntimeError("Direct 153-node z10 spectrum is incomplete")
    base_signal = np.power(10.0, base_logs)
    if not np.all(np.isfinite(base_signal) & (base_signal > 0.0)):
        raise RuntimeError("Direct 153-node raw spectrum is not positive and finite")
    pchip = PchipInterpolator(freqs153[::-1], base_signal[::-1], extrapolate=False)
    predicted = np.asarray(pchip(mid153), dtype=np.float64)
    prediction_valid = np.isfinite(predicted) & (predicted > 0.0)

    partial: dict[str, object] = {
        "schema_version": 1,
        "status": "running_sobol_index3_z10_305_frequency_refinement",
        "scope": "Continue the directly solved Sobol index 3 z10 nested grid from 153 to 305 nodes by inserting every adjacent midpoint. Directly solve the added nodes and continue Reference root mapping until a root-relative residual closes.",
        "input_artifact": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "fixed_root": root,
        "frequency_grids": {
            "base_77_count": int(freqs77.size),
            "base_77_sha256_f64le": digest(freqs77),
            "base_153_count": int(freqs153.size),
            "base_153_sha256_f64le": digest(freqs153),
            "refined_305_count": int(freqs305.size),
            "refined_305_sha256_f64le": digest(freqs305),
            "construction": "77 native nodes, 76 direct midpoint nodes, then 152 direct midpoints of the 153-node mesh",
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_artifact": sha(INPUT),
            "script": sha(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
            "fast_replay_runtime_s": fast_runtime,
        },
        "root_iteration": {
            "root_residual_tolerance": ROOT_RESIDUAL_TOL,
            "max_iterations": MAX_ROOT_ITER,
            "iterations": [],
        },
    }
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    direct_model = LCDM_SG(**params)
    direct, runtime, caught = solve_reference(direct_model, root, freqs305)
    direct_logs, direct_valid = raw_log_spectrum(direct)
    if not np.all(direct_valid):
        raise RuntimeError("Direct 305-node z10 spectrum contains invalid raw nodes")
    mapped_initial = float(direct["DN_gw"])
    common_delta = direct_logs[::2] - base_logs
    new_node_logs = direct_logs[1::2]
    midpoint_delta = np.log10(predicted[prediction_valid]) - new_node_logs[prediction_valid]
    fixed_record = {
        "input_DN_eff": root,
        "mapped_DN_gw": mapped_initial,
        "relative_map_residual_vs_fast_root": mapped_initial / root - 1.0,
        "runtime_s": runtime,
        "quadrature_error": float(direct["quadrature_error"]),
        "interpolation_error": float(direct["interpolation_error"]),
        "raw_signal_positive_finite_nodes": int(np.count_nonzero(direct_valid)),
        "runtime_warnings": caught,
        "map_relative_change_vs_153_nodes": mapped_initial / float(previous["fixed_root_reference"]["refined_153_nodes"]["mapped_DN_gw"]) - 1.0,
        "map_relative_change_vs_77_nodes": mapped_initial / float(previous["fixed_root_reference"]["base_77_nodes"]["mapped_DN_gw"]) - 1.0,
        "common_153_node_spectrum_delta_max_abs_dex": float(np.max(np.abs(common_delta))),
        "common_153_node_spectrum_delta_p95_abs_dex": float(np.percentile(np.abs(common_delta), 95.0)),
        "new_midpoint_direct_vs_153_node_linear_pchip_dex": {
            "new_node_count": int(mid153.size),
            "compared_nodes": int(midpoint_delta.size),
            "uncompared_nodes": int(mid153.size - midpoint_delta.size),
            "max_abs_dex": float(np.max(np.abs(midpoint_delta))) if midpoint_delta.size else None,
            "p95_abs_dex": float(np.percentile(np.abs(midpoint_delta), 95.0)) if midpoint_delta.size else None,
        },
    }
    partial["status"] = "running_sobol_index3_z10_305_frequency_refinement"
    partial["fixed_root_reference"] = fixed_record
    partial["spectra"] = {"direct_z10_log10_OmegaGW_at_fast_root": direct_logs.tolist()}
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    x = mapped_initial
    for iteration in range(1, MAX_ROOT_ITER + 1):
        iteration_model = LCDM_SG(**params)
        mapped_result, iter_runtime, iter_warnings = solve_reference(iteration_model, x, freqs305)
        mapped_logs, mapped_valid = raw_log_spectrum(mapped_result)
        if not np.all(mapped_valid):
            raise RuntimeError("305-node root-iteration spectrum contains invalid raw nodes")
        mapped = float(mapped_result["DN_gw"])
        residual = mapped / x - 1.0
        spectrum_delta = mapped_logs - direct_logs
        record = {
            "iteration": iteration, "input_DN_eff": x, "mapped_DN_gw": mapped,
            "relative_map_residual": residual, "runtime_s": iter_runtime,
            "quadrature_error": float(mapped_result["quadrature_error"]),
            "interpolation_error": float(mapped_result["interpolation_error"]),
            "raw_signal_positive_finite_nodes": int(np.count_nonzero(mapped_valid)),
            "runtime_warnings": iter_warnings,
            "spectrum_change_vs_305_map_at_fast_root": {
                "compared_nodes": int(spectrum_delta.size),
                "max_abs_dex": float(np.max(np.abs(spectrum_delta))),
                "p95_abs_dex": float(np.percentile(np.abs(spectrum_delta), 95.0)),
            },
        }
        partial["root_iteration"]["iterations"].append(record)
        partial["spectra"]["direct_z10_log10_OmegaGW_at_updated_root"] = mapped_logs.tolist()
        partial["status"] = (
            "complete_sobol_index3_z10_305_frequency_root_closure"
            if abs(residual) <= ROOT_RESIDUAL_TOL
            else "running_sobol_index3_z10_305_frequency_refinement"
        )
        OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if abs(residual) <= ROOT_RESIDUAL_TOL:
            break
        x = mapped
    else:
        partial["status"] = "reference_305_grid_root_iteration_limit_without_closure"
        OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        raise RuntimeError("305-node Reference root iteration did not close")

    print(json.dumps({
        "output": str(OUTPUT), "nodes": int(freqs305.size),
        "map_change_vs_153": fixed_record["map_relative_change_vs_153_nodes"],
        "map_change_vs_77": fixed_record["map_relative_change_vs_77_nodes"],
        "new_node_pchip_max_abs_dex": fixed_record["new_midpoint_direct_vs_153_node_linear_pchip_dex"]["max_abs_dex"],
        "root_residual": partial["root_iteration"]["iterations"][-1]["relative_map_residual"],
        "root": partial["root_iteration"]["iterations"][-1]["mapped_DN_gw"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
