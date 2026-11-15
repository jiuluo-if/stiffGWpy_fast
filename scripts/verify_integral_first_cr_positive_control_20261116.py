"""Verify the integral-first change leaves the physical cr>0 path untouched."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_consistency_branch_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "integral_first_cr_positive_control_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def main() -> int:
    prior = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    params = dict(prior["parameters"])
    if params["cr"] <= 0.0 or params["A_s"] <= 0.0 or params["r"] <= 0.0:
        raise RuntimeError("the selected control must remain on a positive-amplitude cr>0 branch")
    ref = prior["reference_fixed_root_physical_amplitude"]["10"]
    root = float(prior["fast_root_DN_eff"])

    root_model = LCDM_SG(**params)
    started = time.perf_counter()
    root_result = root_model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    root_runtime = time.perf_counter() - started
    if root_result is not root_model:
        raise RuntimeError(f"current cr>0 Fast control failed: {root_model.fast_failure_reason}")
    freqs = np.asarray(root_model.f, dtype=np.float64)
    current_root = float(root_model.cosmo_param["DN_eff"])
    grid_sha = digest(freqs)
    if freqs.size != prior["frequency_count"] or grid_sha != prior["frequency_grid_sha256_f64le"]:
        raise RuntimeError("current cr>0 goal grid differs from the direct Reference artifact")

    fixed_model = LCDM_SG(**{**params, "DN_eff": root})
    fixed_model.gen_expansion()
    stored_amplitude = float(fixed_model.derived_param["A_t"])
    log_amplitude = float(_log_tensor_amplitude(fixed_model.cosmo_param, fixed_model.derived_param))
    if not (stored_amplitude > 0.0 and np.isfinite(stored_amplitude) and np.isfinite(log_amplitude)):
        raise RuntimeError("the physical cr>0 control unexpectedly entered an underflowed-A_t state")

    original_gen = fixed_model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    original_map = FS._integrate_log_scaled_dn_map
    calls = []

    def locked_gen(*args, **kwargs):
        fixed_model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    def map_spy(*args, **kwargs):
        calls.append(True)
        return original_map(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    FS._integrate_log_scaled_dn_map = map_spy
    fixed_model.gen_expansion = locked_gen
    fixed_started = time.perf_counter()
    try:
        fixed_result = fixed_model.SGWB_iter(
            engine="fast", fallback=False, tol=0.99, z_tail=8.0,
            threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
        FS._integrate_log_scaled_dn_map = original_map
    fixed_runtime = time.perf_counter() - fixed_started
    if fixed_result is not fixed_model:
        raise RuntimeError(f"current cr>0 fixed-root Fast failed: {fixed_model.fast_failure_reason}")

    fast_spectrum = np.asarray(fixed_model.log10OmegaGW, dtype=np.float64)
    ref_spectrum = np.asarray(ref["log10_OmegaGW"], dtype=np.float64)
    delta = fast_spectrum - ref_spectrum
    worst = int(np.argmax(np.abs(delta)))
    ref_map = float(ref["mapped_DN_gw"])
    fast_map = float(fixed_model.DN_gw[-1])
    output = {
        "schema_version": 1,
        "status": "complete_current_integral_first_cr_positive_physical_control",
        "scope": "Physical-amplitude cr>0 branch control at r=.01; current Fast self-consistent root and fixed-root z8 spectrum/map compared on the exact 76-node direct Reference z10 grid. The new log-scaled integral helper is monitored and must remain unused.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_sha,
        "reference_artifact": {
            "path": str(REFERENCE_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "z_tail": 10.0,
            "rtol": 1e-9,
            "frozen_DN_eff": root,
            "mapped_DN_gw": ref_map,
            "reference_source_sha256_matches_current": (
                prior["source_sha256"]["reference.py"] == sha(ROOT / "stiffgwpy_fast" / "reference.py")
            ),
            "root_spectrum_max_abs_dex": prior["spectrum_error_summaries_dex"]["fast_z8_vs_reference_z10"]["max_abs_dex"],
        },
        "current_fast": {
            "self_consistent_returned": True,
            "self_consistent_root_DN_eff": current_root,
            "self_consistent_map_DN_gw": float(root_model.DN_gw[-1]),
            "self_consistent_root_delta": float(root_model.dn_converged_delta),
            "self_consistent_runtime_s": root_runtime,
            "fixed_root_z8_mapped_DN_gw": fast_map,
            "fixed_root_z8_residual_relative": fast_map / root - 1.0,
            "fixed_root_z8_runtime_s": fixed_runtime,
            "failure_reason": fixed_model.fast_failure_reason,
            "failure_detail": getattr(fixed_model, "fast_failure_detail", None),
        },
        "physical_background": {
            "stored_A_t": stored_amplitude,
            "log_A_t": log_amplitude,
            "nonzero_physical_A_t": True,
            "integral_first_helper_call_count": len(calls),
            "cr_positive_background_preserved": fixed_model.cosmo_param["cr"] == params["cr"],
            "N_inf": float(fixed_model.derived_param["N_inf"]),
            "V_inf": float(fixed_model.derived_param["V_inf"]),
        },
        "spectrum_error_vs_reference_z10": {
            "max_abs_dex": float(np.max(np.abs(delta))),
            "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
            "worst_log10f": float(freqs[worst]),
            "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
            "fast_log10_OmegaGW": fast_spectrum.tolist(),
            "reference_log10_OmegaGW": ref_spectrum.tolist(),
        },
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
                "reference.py": "stiffgwpy_fast/reference.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "reference_control_input": "docs/research/results/underflow_tail_depth_consistency_branch_20261115.json",
                "script": "scripts/verify_integral_first_cr_positive_control_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "current_root": current_root,
        "reference_root": root,
        "fixed_residual": output["current_fast"]["fixed_root_z8_residual_relative"],
        "spectrum_max": output["spectrum_error_vs_reference_z10"]["max_abs_dex"],
        "integral_first_calls": len(calls),
        "stored_A_t": stored_amplitude,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
