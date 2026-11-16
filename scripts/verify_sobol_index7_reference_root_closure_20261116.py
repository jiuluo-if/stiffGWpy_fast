"""Continue the direct Reference fixed-point iteration for Sobol index 7."""

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

from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index7_reference_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index7_reference_root_closure_20261116.json"
REFERENCE_RTOL = 1e-9
ROOT_RESIDUAL_TOL = 2e-4
MAX_ITER = 5


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


def main() -> int:
    direct = json.loads(INPUT.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in direct["parameters"].items()}
    fast_root = float(direct["current_fast"]["root_DN_eff"])
    start_fast = time.perf_counter()
    model = LCDM_SG(**params)
    result = model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    if result is not model:
        raise RuntimeError(f"Sobol index 7 current Fast failed: {model.fast_failure_reason}")
    fast_runtime = time.perf_counter() - start_fast
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_root_now = float(model.cosmo_param["DN_eff"])
    if digest(freqs) != direct["current_fast"]["frequency_sha256_f64le"]:
        raise RuntimeError("Sobol index 7 Fast frequency grid differs from the direct Reference artifact")
    if not math.isclose(fast_root_now, fast_root, rel_tol=2e-12, abs_tol=0.0):
        raise RuntimeError("Sobol index 7 Fast root differs from the direct Reference artifact")

    initial_map = direct["reference"]["z_tails"]["12.0"]
    start_root = float(initial_map["DN_gw"])
    initial_logs = np.asarray([np.nan if v is None else float(v) for v in initial_map["raw_log10_OmegaGW"]])
    partial = {
        "schema_version": 1,
        "status": "running_sobol_index7_reference_root_iteration",
        "scope": "Continue fixed-point iteration from the directly solved z12 Reference map at the current Fast root. Every additional map is directly solved at the previous Reference-mapped root on the exact same native frequency grid; stop only when the relative map residual is below the root screen.",
        "input_artifact": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "fixed_frequency_grid": {"count": int(freqs.size), "sha256_f64le": digest(freqs)},
        "current_fast": {
            "root_DN_eff": fast_root_now,
            "runtime_s": fast_runtime,
            "fixed_point_residual_relative": float(model.DN_gw[-1] / fast_root_now - 1.0),
        },
        "initial_z12_reference_at_fast_root": {
            "input_DN_eff": fast_root,
            "mapped_DN_gw": start_root,
            "relative_map_residual": start_root / fast_root - 1.0,
            "runtime_s": float(initial_map["runtime_s"]),
        },
        "reference": {
            "z_tail": 12.0,
            "rtol": REFERENCE_RTOL,
            "root_residual_tolerance": ROOT_RESIDUAL_TOL,
            "max_iterations": MAX_ITER,
            "iterations": [],
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_artifact": sha(INPUT),
            "script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    x = start_root
    for iteration in range(1, MAX_ITER + 1):
        ref_model = LCDM_SG(**params)
        caught: list[warnings.WarningMessage]
        start = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ref = run_reference(
                ref_model, dn_eff=x, freq_subset=freqs, z_tail=12.0,
                rtol=REFERENCE_RTOL, self_consistent=False, physical_guard=True,
            )
        mapped = float(ref["DN_gw"])
        residual = mapped / x - 1.0
        logs, valid = raw_log_spectrum(ref)
        spectrum_valid = valid & np.isfinite(initial_logs)
        delta = logs[spectrum_valid] - initial_logs[spectrum_valid]
        record = {
            "iteration": iteration,
            "input_DN_eff": x,
            "mapped_DN_gw": mapped,
            "relative_map_residual": residual,
            "runtime_s": time.perf_counter() - start,
            "quadrature_error": float(ref["quadrature_error"]),
            "interpolation_error": float(ref["interpolation_error"]),
            "raw_signal_positive_finite_nodes": int(np.count_nonzero(valid)),
            "raw_signal_invalid_nodes": int(valid.size - np.count_nonzero(valid)),
            "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught],
            "spectrum_change_vs_z12_at_fast_root": {
                "compared_nodes": int(delta.size),
                "max_abs_dex": float(np.max(np.abs(delta))) if delta.size else None,
                "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)) if delta.size else None,
            },
        }
        partial["reference"]["iterations"].append(record)
        partial["status"] = (
            "complete_sobol_index7_reference_root_iteration"
            if abs(residual) <= ROOT_RESIDUAL_TOL
            else "running_sobol_index7_reference_root_iteration"
        )
        OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if abs(residual) <= ROOT_RESIDUAL_TOL:
            break
        x = mapped
    else:
        partial["status"] = "reference_root_iteration_limit_without_closure"
        OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        raise RuntimeError(f"z12 Reference root residual exceeds {ROOT_RESIDUAL_TOL} after {MAX_ITER} iterations")

    print(json.dumps({
        "output": str(OUTPUT),
        "iterations": len(partial["reference"]["iterations"]),
        "final_root_residual": partial["reference"]["iterations"][-1]["relative_map_residual"],
        "final_root": partial["reference"]["iterations"][-1]["mapped_DN_gw"],
        "final_spectrum_change_max_abs_dex": partial["reference"]["iterations"][-1]["spectrum_change_vs_z12_at_fast_root"]["max_abs_dex"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
