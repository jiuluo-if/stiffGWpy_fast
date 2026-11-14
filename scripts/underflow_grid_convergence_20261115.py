"""Direct nested-grid underflow screen for Fast and its scaled Reference oracle."""

from __future__ import annotations

import argparse
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

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.freq_adaptive import (  # noqa: E402
    _log_tensor_amplitude,
    construct_fast_frequency_grid,
)
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-316,
    "n_t": 17.38,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
Z_TAIL_FAST = 5.0
Z_TAIL_REFERENCE = 8.0
REFERENCE_RTOL = 1e-9
THREADS = 1
GRID_COUNTS = (235, 469, 937)
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _grid_sha256(freqs: np.ndarray) -> str:
    return _sha256_bytes(np.asarray(freqs, dtype="<f8").tobytes(order="C"))


def _source_hashes() -> dict[str, str]:
    paths = (
        "stiffgwpy_fast/fast_sgwb.py",
        "stiffgwpy_fast/freq_adaptive.py",
        "stiffgwpy_fast/exact_background.py",
        "stiffgwpy_fast/reference.py",
        "stiffgwpy_fast/stiff_SGWB.py",
        "scripts/underflow_grid_convergence_20261115.py",
    )
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def _nested_midpoints(freqs_desc: np.ndarray) -> np.ndarray:
    ascending = np.asarray(freqs_desc, dtype=np.float64)[::-1]
    refined = np.empty(2 * ascending.size - 1, dtype=np.float64)
    refined[::2] = ascending
    refined[1::2] = 0.5 * (ascending[:-1] + ascending[1:])
    return refined[::-1]


def _fast_run(freqs: np.ndarray, *, frozen_dn: float | None = None) -> dict:
    model = LCDM_SG(**PARAMS)
    model.gen_expansion()
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    original_gen = model.gen_expansion

    def fixed_grid(*args, **kwargs):
        return np.asarray(freqs, dtype=np.float64).copy()

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = float(frozen_dn)
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    if frozen_dn is not None:
        model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99 if frozen_dn is not None else 1e-7,
            z_tail=Z_TAIL_FAST,
            threads=THREADS,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(
            f"Fast failed on {freqs.size} nodes: {model.fast_failure_reason}"
        )
    if not np.array_equal(model.f, freqs):
        raise AssertionError("Fast did not solve the requested native frequency grid")
    mapped = float(model.DN_gw[-1])
    return {
        "frequency_count": int(model.f.size),
        "frequency_sha256_f64le": _grid_sha256(model.f),
        "root_DN_eff": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": mapped,
        "map_residual_at_frozen_background": (
            None if frozen_dn is None else mapped - float(frozen_dn)
        ),
        "runtime_s": elapsed,
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64).tolist(),
        "linear_OmegaGW": np.asarray(model.Ogw_today - model.Oj_today, dtype=np.float64).tolist(),
        "frequency_log10_hz": np.asarray(model.f, dtype=np.float64).tolist(),
        "failure_reason": model.fast_failure_reason,
        "grid_requested": model.freq_grid_requested,
        "grid_used": model.freq_grid_used,
    }


def _reference_unit_model() -> LCDM_SG:
    params = dict(PARAMS)
    params.update({"DN_eff": 0.0, "A_s": 1.0, "r": 1.0, "n_t": 0.0})
    model = LCDM_SG(**params)
    model.gen_expansion()
    return model


def _reference_map(
    model: LCDM_SG, freqs: np.ndarray, frozen_dn: float, log_amplitude: float
) -> dict:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        Ogw, Oj, Opgw, used_tail = spectrum_reference(
            model,
            freqs,
            frozen_dn,
            z_tail=Z_TAIL_REFERENCE,
            rtol=REFERENCE_RTOL,
            workers=THREADS,
            backend="cartesian",
        )
        base = np.where(used_tail, 3.0 * Opgw, Ogw - Oj)
        if np.any(~np.isfinite(base)) or np.any(base <= 0.0):
            raise FloatingPointError("Reference unit-transfer spectrum is not positive and finite")
        log_power = log_amplitude + PARAMS["n_t"] * (
            math.log(10.0) * freqs - math.log(gp.f_piv)
        )
        log_integrand = np.log(base) + log_power
        shift = float(np.max(log_integrand))
        normalized = np.exp(log_integrand - shift)
        g2_scaled, qerr_scaled, ierr_scaled = integrate_spectrum(
            freqs, normalized, np.zeros_like(normalized)
        )
    elapsed = time.perf_counter() - started
    log_g2 = shift + math.log(g2_scaled)
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    dn_gw = gp.Neff0 * math.exp(log_g2) / omega_nu
    return {
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": _grid_sha256(freqs),
        "input_DN_eff": float(frozen_dn),
        "mapped_DN_gw": float(dn_gw),
        "relative_map_residual": float(dn_gw / frozen_dn - 1.0),
        "runtime_s": elapsed,
        "log10_OmegaGW": (log_integrand / math.log(10.0)).tolist(),
        "log_integrand_ln": log_integrand.tolist(),
        "normalized_integrand": normalized.tolist(),
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "normalized_nodes_underflowed_to_zero": int(np.count_nonzero(normalized == 0.0)),
        "quadrature_error_relative": (
            None if qerr_scaled is None else float(qerr_scaled / g2_scaled)
        ),
        "interpolation_error_relative": (
            None if ierr_scaled is None else float(ierr_scaled / g2_scaled)
        ),
        "runtime_warnings": [str(item.message) for item in caught],
    }


def _interpolation_diagnostic(
    coarse_freqs: np.ndarray,
    coarse_values: np.ndarray,
    fine_freqs: np.ndarray,
    fine_values: np.ndarray,
) -> dict:
    x = coarse_freqs[::-1]
    y = coarse_values[::-1]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        interpolator = PchipInterpolator(x, y)
    coarse_set = set(np.asarray(coarse_freqs, dtype=np.float64).tolist())
    new_mask = np.array([float(f) not in coarse_set for f in fine_freqs])
    new_freqs = fine_freqs[new_mask]
    actual = fine_values[new_mask]
    predicted = np.asarray(interpolator(new_freqs), dtype=np.float64)
    valid = (actual > 0.0) & (predicted > 0.0)
    errors = np.abs(np.log10(actual[valid]) - np.log10(predicted[valid]))
    return {
        "new_node_count": int(new_freqs.size),
        "positive_comparable_count": int(np.count_nonzero(valid)),
        "zero_interpolated_count": int(np.count_nonzero(predicted == 0.0)),
        "runtime_warnings": [str(item.message) for item in caught],
        "max_abs_error_dex": None if errors.size == 0 else float(np.max(errors)),
        "p95_abs_error_dex": None if errors.size == 0 else float(np.quantile(errors, 0.95)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--grid-counts",
        nargs="+",
        type=int,
        choices=GRID_COUNTS,
        default=GRID_COUNTS,
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    counts = tuple(args.grid_counts)
    if counts != GRID_COUNTS[: len(counts)]:
        parser.error("grid counts must be a prefix of 235 469 937")

    physical = LCDM_SG(**PARAMS)
    physical.gen_expansion()
    log_amplitude = _log_tensor_amplitude(physical.cosmo_param, physical.derived_param)
    base_grid = construct_fast_frequency_grid(physical, log_amplitude=log_amplitude)
    grids = {235: np.asarray(base_grid, dtype=np.float64)}
    for previous, next_count in zip(counts, counts[1:]):
        grids[next_count] = _nested_midpoints(grids[previous])
    if tuple((count, grids[count].size) for count in counts) != tuple(zip(counts, counts)):
        raise AssertionError("nested frequency counts are not 235/469/937")

    fast_roots = {
        count: _fast_run(grid)
        for count, grid in grids.items()
    }
    common_dn = float(fast_roots[235]["root_DN_eff"])
    fast_fixed = {
        count: _fast_run(grid, frozen_dn=common_dn)
        for count, grid in grids.items()
    }
    unit_reference = _reference_unit_model()
    if not np.array_equal(physical.Nv, unit_reference.Nv):
        raise AssertionError("unit Reference surrogate changed the physical background grid")
    reference_fixed = {
        count: _reference_map(unit_reference, grid, common_dn, log_amplitude)
        for count, grid in grids.items()
    }
    fast_root_checks = {
        count: (
            fast_fixed[count]
            if math.isclose(float(fast_roots[count]["root_DN_eff"]), common_dn, rel_tol=0.0, abs_tol=0.0)
            else _fast_run(grid, frozen_dn=float(fast_roots[count]["root_DN_eff"]))
        )
        for count, grid in grids.items()
    }
    reference_root_checks = {
        count: (
            reference_fixed[count]
            if math.isclose(float(fast_roots[count]["root_DN_eff"]), common_dn, rel_tol=0.0, abs_tol=0.0)
            else _reference_map(
                unit_reference,
                grid,
                float(fast_roots[count]["root_DN_eff"]),
                log_amplitude,
            )
        )
        for count, grid in grids.items()
    }
    for count in counts:
        fast_roots[count]["fixed_map_residual_at_own_root"] = float(
            fast_root_checks[count]["mapped_DN_gw"] - fast_roots[count]["root_DN_eff"]
        )
        reference_root_checks[count]["map_residual_at_fast_root"] = float(
            reference_root_checks[count]["mapped_DN_gw"] - fast_roots[count]["root_DN_eff"]
        )

    common_log_scale = max(
        float(np.max(reference_fixed[count]["log_integrand_ln"]))
        for count in counts
    )
    normalized_reference = {}
    for count in counts:
        normalized = np.exp(
            np.asarray(reference_fixed[count]["log_integrand_ln"]) - common_log_scale
        )
        normalized_reference[count] = normalized
        reference_fixed[count]["normalized_integrand_common_scale"] = normalized.tolist()
        reference_fixed[count]["common_scale_zero_count"] = int(
            np.count_nonzero(normalized == 0.0)
        )

    interpolation = {}
    pairs = tuple(
        (coarse, fine)
        for coarse, fine in ((235, 469), (469, 937), (235, 937))
        if coarse in grids and fine in grids
    )
    for coarse_count, fine_count in pairs:
        coarse = normalized_reference[coarse_count]
        fine = normalized_reference[fine_count]
        interpolation[f"{coarse_count}_to_{fine_count}"] = _interpolation_diagnostic(
            grids[coarse_count], coarse, grids[fine_count], fine
        )

    spectrum_errors = {}
    for count in counts:
        fast_log = np.asarray(fast_fixed[count]["log10_OmegaGW"])
        ref_log = np.asarray(reference_fixed[count]["log10_OmegaGW"])
        errors = np.abs(fast_log - ref_log)
        spectrum_errors[str(count)] = {
            "max_abs_error_dex": float(np.max(errors)),
            "p95_abs_error_dex": float(np.quantile(errors, 0.95)),
            "count": int(errors.size),
        }

    output = {
        "schema_version": 1,
        "status": "complete_fixed_background_direct_nested_grid_screen",
        "experiment": "current-source Fast and scaled Cartesian Reference on nested 235/469/937 grids",
        "parameters": PARAMS,
        "background": {
            "common_DN_eff": common_dn,
            "scope": "cr=0, DN_re=10; unit Reference surrogate matches Nv and f_hor exactly",
            "fast_tail": Z_TAIL_FAST,
            "reference_tail": Z_TAIL_REFERENCE,
            "reference_rtol": REFERENCE_RTOL,
        },
        "frequency_grids": {
            str(count): {
                "count": int(grid.size),
                "sha256_f64le": _grid_sha256(grid),
                "log10_hz": grid.tolist(),
            }
            for count, grid in grids.items()
        },
        "fast_self_consistent": {str(k): v for k, v in fast_roots.items()},
        "fast_fixed_background": {str(k): v for k, v in fast_fixed.items()},
        "fast_root_residual_checks": {str(k): v for k, v in fast_root_checks.items()},
        "reference_fixed_background": {str(k): v for k, v in reference_fixed.items()},
        "reference_map_at_fast_root": {str(k): v for k, v in reference_root_checks.items()},
        "spectrum_error_fast_vs_reference_dex": spectrum_errors,
        "reference_pchip_new_node_diagnostic": interpolation,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": THREADS,
        },
        "source_sha256": _source_hashes(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "common_DN_eff": common_dn,
        "fast_roots": {k: v["root_DN_eff"] for k, v in output["fast_self_consistent"].items()},
        "reference_maps": {k: v["mapped_DN_gw"] for k, v in output["reference_fixed_background"].items()},
        "reference_zero_nodes": {
            k: v["normalized_nodes_underflowed_to_zero"]
            for k, v in output["reference_fixed_background"].items()
        },
        "warnings": {
            k: v["runtime_warnings"]
            for k, v in output["reference_fixed_background"].items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
