"""Verify a valid cr>0 point immediately above the derived cutoff."""

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
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BOUNDARY_INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_cutoff_amplitude_boundary_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_cutoff_valid_sample_20261116.json"
RATIO_ABOVE_CUTOFF = 1.1
CASE = {"T_re": 2000.0, "kappa10": 0.01}
BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "n_t": 0.0,
    "cr": 1.0,
    "DN_re": 20.0,
    "DN_eff": 0.0,
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def main() -> int:
    boundary = json.loads(BOUNDARY_INPUT.read_text(encoding="utf-8"))
    boundary_case = next(
        row for row in boundary["cases"]
        if row["T_re"] == CASE["T_re"] and row["kappa10"] == CASE["kappa10"]
    )
    r_cutoff = float(boundary_case["cutoff_root_r"])
    params = {**BASE, **CASE, "r": r_cutoff * RATIO_ABOVE_CUTOFF}

    root_model = LCDM_SG(**params)
    initial = root_model.derived_param
    initial_delta_n = float(initial["N_inf"] - initial["N_re"])
    stored_at = float(initial["A_t"])
    if initial["N_inf"] is None or initial_delta_n <= 0.0:
        raise RuntimeError("selected point is not on the valid cutoff side")
    if stored_at <= 0.0 or not math.isfinite(stored_at):
        raise RuntimeError("selected point does not have a representable positive A_t")

    started = time.perf_counter()
    result = root_model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    fast_runtime = time.perf_counter() - started
    if result is not root_model:
        raise RuntimeError(f"current Fast self-consistent solve failed: {root_model.fast_failure_reason}")
    root = float(root_model.cosmo_param["DN_eff"])
    freqs = np.asarray(root_model.f, dtype=np.float64)
    grid_sha = digest(freqs)

    # Compute Reference directly on this newly generated Fast grid and exact Fast root.
    reference_model = LCDM_SG(**{**params, "DN_eff": root})
    warnings_seen: list[str] = []
    ref_results = {}
    ref_runtime = {}
    for z_tail in (8.0, 10.0, 12.0):
        tic = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ref_results[z_tail] = run_reference(
                reference_model, dn_eff=root, freq_subset=freqs,
                z_tail=z_tail, rtol=1e-9, self_consistent=False,
                physical_guard=True,
            )
        ref_runtime[z_tail] = time.perf_counter() - tic
        warnings_seen.extend(f"z_tail={z_tail}: {w.category.__name__}: {w.message}" for w in caught)

    # Re-evaluate Fast at the same root and exact same native grid.
    fixed_model = LCDM_SG(**{**params, "DN_eff": root})
    original_gen = fixed_model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    calls: list[bool] = []
    original_log_map = FS._integrate_log_scaled_dn_map

    def locked_gen(*args, **kwargs):
        fixed_model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    def map_spy(*args, **kwargs):
        calls.append(True)
        return original_log_map(*args, **kwargs)

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
        FS._integrate_log_scaled_dn_map = original_log_map
    fixed_runtime = time.perf_counter() - fixed_started
    if fixed_result is not fixed_model:
        raise RuntimeError(f"current Fast fixed-root solve failed: {fixed_model.fast_failure_reason}")

    ref8 = ref_results[8.0]
    ref10 = ref_results[10.0]
    ref12 = ref_results[12.0]
    fast_spectrum = np.asarray(fixed_model.log10OmegaGW, dtype=np.float64)
    raw_ref_signal = np.asarray(ref10["Ogw"] - ref10["Oj"], dtype=np.float64)
    if not np.all(np.isfinite(raw_ref_signal) & (raw_ref_signal > 0.0)):
        raise FloatingPointError("direct Reference raw spectrum contains non-positive or non-finite values")
    ref_spectrum = np.log10(raw_ref_signal)
    raw_ref8_signal = np.asarray(ref8["Ogw"] - ref8["Oj"], dtype=np.float64)
    if not np.all(np.isfinite(raw_ref8_signal) & (raw_ref8_signal > 0.0)):
        raise FloatingPointError("direct Reference z8 raw spectrum contains non-positive or non-finite values")
    ref8_spectrum = np.log10(raw_ref8_signal)
    raw_ref12_signal = np.asarray(ref12["Ogw"] - ref12["Oj"], dtype=np.float64)
    if not np.all(np.isfinite(raw_ref12_signal) & (raw_ref12_signal > 0.0)):
        raise FloatingPointError("direct Reference z12 raw spectrum contains non-positive or non-finite values")
    ref12_spectrum = np.log10(raw_ref12_signal)
    delta = fast_spectrum - ref_spectrum
    delta_z12 = fast_spectrum - ref12_spectrum
    worst = int(np.argmax(np.abs(delta)))
    worst_z12 = int(np.argmax(np.abs(delta_z12)))
    fast_map = float(fixed_model.DN_gw[-1])
    ref_map = float(ref10["DN_gw"])
    fast_fixed_background = LCDM_SG(**{**params, "DN_eff": root})
    fast_fixed_background.gen_expansion()
    final_derived = fast_fixed_background.derived_param
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_valid_cutoff_direct_reference_comparison",
        "scope": "One physical-amplitude cr>0 sample at 10% above the analytically derived cutoff for T_re=2000 and kappa10=.01. Fast self-consistency and same-root Fast-vs-direct-Reference z8/z10 comparison use the exact newly generated Fast frequency grid. This is a discrete point, not a domain certification.",
        "parameters": params,
        "cutoff": {
            "input_artifact": str(BOUNDARY_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "r_cutoff": r_cutoff,
            "ratio_above_cutoff": RATIO_ABOVE_CUTOFF,
            "r_over_r_cutoff": params["r"] / r_cutoff,
            "initial_Delta_N": initial_delta_n,
            "initial_N_inf": float(initial["N_inf"]),
            "initial_N_re": float(initial["N_re"]),
        },
        "frequency_grid": {
            "count": int(freqs.size),
            "log10f_sha256_f64le": grid_sha,
            "log10f": freqs.tolist(),
        },
        "current_fast": {
            "self_consistent_returned": True,
            "self_consistent_root_DN_eff": root,
            "self_consistent_map_DN_gw": float(root_model.DN_gw[-1]),
            "self_consistent_root_delta": float(root_model.dn_converged_delta),
            "runtime_s": fast_runtime,
            "same_root_fixed_background_DN_gw": fast_map,
            "same_root_fixed_background_residual_relative": fast_map / root - 1.0,
            "same_root_runtime_s": fixed_runtime,
            "integral_first_helper_call_count": len(calls),
        },
        "direct_reference": {
            "source_sha256": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fixed_DN_eff": root,
            "frequency_count": int(freqs.size),
            "z8": {
                "DN_gw": float(ref8["DN_gw"]),
                "quadrature_error": float(ref8["quadrature_error"]),
                "interpolation_error": float(ref8["interpolation_error"]),
                "used_tail_fraction": float(np.mean(np.asarray(ref8["used_tail"], dtype=bool))),
                "runtime_s": ref_runtime[8.0],
            },
            "z10": {
                "DN_gw": ref_map,
                "quadrature_error": float(ref10["quadrature_error"]),
                "interpolation_error": float(ref10["interpolation_error"]),
                "used_tail_fraction": float(np.mean(np.asarray(ref10["used_tail"], dtype=bool))),
                "runtime_s": ref_runtime[10.0],
            },
            "z12": {
                "DN_gw": float(ref12["DN_gw"]),
                "quadrature_error": float(ref12["quadrature_error"]),
                "interpolation_error": float(ref12["interpolation_error"]),
                "used_tail_fraction": float(np.mean(np.asarray(ref12["used_tail"], dtype=bool))),
                "runtime_s": ref_runtime[12.0],
            },
            "z8_to_z10_map_relative_change": float(ref_map / float(ref8["DN_gw"]) - 1.0),
            "z10_to_z12_map_relative_change": float(float(ref12["DN_gw"]) / ref_map - 1.0),
            "z8_to_z10_max_abs_dex": float(np.max(np.abs(ref8_spectrum - ref_spectrum))),
            "z10_to_z12_max_abs_dex": float(np.max(np.abs(ref12_spectrum - ref_spectrum))),
            "public_log_spectrum_floor": -40.0,
            "z10_public_log_spectrum_floor_nodes": int(np.count_nonzero(np.asarray(ref10["log10OmegaGW"]) == -40.0)),
            "raw_positive_log_spectrum_used_for_comparison": True,
            "z8_raw_log10_OmegaGW": ref8_spectrum.tolist(),
            "z10_raw_log10_OmegaGW": ref_spectrum.tolist(),
            "z12_raw_log10_OmegaGW": ref12_spectrum.tolist(),
            "warnings": warnings_seen,
        },
        "comparison": {
            "fast_z8_vs_reference_z10_max_abs_dex": float(np.max(np.abs(delta))),
            "fast_z8_vs_reference_z10_p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
            "fast_z8_vs_reference_z10_worst_log10f": float(freqs[worst]),
            "fast_z8_vs_reference_z10_passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
            "fast_z8_vs_reference_z12_max_abs_dex": float(np.max(np.abs(delta_z12))),
            "fast_z8_vs_reference_z12_p95_abs_dex": float(np.percentile(np.abs(delta_z12), 95.0)),
            "fast_z8_vs_reference_z12_worst_log10f": float(freqs[worst_z12]),
            "fast_z8_vs_reference_z12_passes_1e-3_dex": bool(np.max(np.abs(delta_z12)) <= 1e-3),
            "fast_fixed_root_map_relative_error_vs_reference_z10": float(fast_map / ref_map - 1.0),
            "fast_fixed_root_map_relative_error_vs_reference_z12": float(fast_map / float(ref12["DN_gw"]) - 1.0),
            "reference_z8_to_z10_map_change": float(ref_map / float(ref8["DN_gw"]) - 1.0),
            "reference_z10_to_z12_map_change": float(float(ref12["DN_gw"]) / ref_map - 1.0),
            "fast_log10_OmegaGW": fast_spectrum.tolist(),
            "reference_z10_log10_OmegaGW": ref_spectrum.tolist(),
        },
        "physical_background_at_fast_root": {
            "A_t": float(final_derived["A_t"]),
            "N_inf": None if final_derived["N_inf"] is None else float(final_derived["N_inf"]),
            "N_re": float(final_derived["N_re"]),
            "cr": float(fast_fixed_background.cosmo_param["cr"]),
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
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "boundary_input": "docs/research/results/cr_positive_cutoff_amplitude_boundary_20261116.json",
                "script": "scripts/verify_cr_positive_cutoff_valid_sample_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "r": params["r"],
        "initial_Delta_N": initial_delta_n,
        "DN_eff": root,
        "fast_fixed_root_map_rel_error": output["comparison"]["fast_fixed_root_map_relative_error_vs_reference_z10"],
        "max_spectrum_error_dex": output["comparison"]["fast_z8_vs_reference_z10_max_abs_dex"],
        "max_spectrum_error_vs_z12_dex": output["comparison"]["fast_z8_vs_reference_z12_max_abs_dex"],
        "reference_tail_map_change": output["direct_reference"]["z8_to_z10_map_relative_change"],
        "reference_z10_to_z12_map_change": output["direct_reference"]["z10_to_z12_map_relative_change"],
        "warnings": warnings_seen,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
