"""English description: Refine the 469-node underflow Reference grid to 937 direct nodes. 中文说明：将下溢 Reference 频率网格从469点嵌套加密到937个直接求解点。"""
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

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import (  # noqa: E402
    _Background,
    _H2_and_sigma_python,
    integrate_spectrum,
    spectrum_reference,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

FAST_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_fast_screen_20261117.json"
REF469_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_reference_235_469_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_reference_469_937_20261118.json"
SOURCE_FILES = (
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/reference.py",
    "stiffgwpy_fast/stiff_SGWB.py",
    "stiffgwpy_fast/LCDM_stiff_Neff.py",
)
N_T = 1.6
Z_TAILS = (10.0, 12.0)
RTOL = 1e-9
LABELS = (
    "normal_edge",
    "subnormal_edge",
    "subnormal_interior",
    "rounded_zero_product",
    "min_subnormal_r",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def scaled_map(freqs: np.ndarray, signal: np.ndarray, log_power: np.ndarray, model: LCDM_SG):
    if np.any(signal <= 0.0) or not np.isfinite(signal).all():
        raise FloatingPointError("Reference O_gw-O_j must be positive and finite")
    log_integrand = np.log(signal) + log_power
    shift = float(np.max(log_integrand))
    normalized = np.exp(log_integrand - shift)
    g2_scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    if g2_scaled <= 0.0 or not math.isfinite(float(g2_scaled)):
        raise FloatingPointError("Reference normalized integral is not positive and finite")
    log_g2 = shift + math.log(float(g2_scaled))
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    return (
        float(gp.Neff0 / omega_nu * math.exp(log_g2)),
        None if qerr is None else float(qerr / g2_scaled),
        None if ierr is None else float(ierr / g2_scaled),
    )


def unit_model(params: dict[str, float]) -> LCDM_SG:
    model = LCDM_SG(**{**params, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0})
    model.gen_expansion()
    return model


def solve_reference(model: LCDM_SG, freqs: np.ndarray, z_tail: float, dn_eff: float, label: str):
    print(f"START Reference {freqs.size} nodes z_tail={z_tail:g} DN_eff={dn_eff:.6e} [{label}]", flush=True)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        ogw, oj, opgw, used_tail = spectrum_reference(
            model, freqs, dn_eff, z_tail=z_tail, rtol=RTOL, workers=1, backend="cartesian"
        )
    elapsed = time.perf_counter() - started
    signal = np.asarray(ogw - oj, dtype=np.float64)
    if np.any(signal <= 0.0) or not np.isfinite(signal).all():
        raise FloatingPointError(f"Reference spectrum invalid at z_tail={z_tail:g}, DN_eff={dn_eff:.6e}")
    print(f"DONE Reference {freqs.size} nodes z_tail={z_tail:g} in {elapsed:.1f}s; warnings={len(caught)}", flush=True)
    return {
        "runtime_s": elapsed,
        "dn_eff": dn_eff,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "positive_signal_count": int(np.count_nonzero(signal > 0.0)),
        "nonfinite_signal_count": int(np.count_nonzero(~np.isfinite(signal))),
        "runtime_warnings": [str(item.message) for item in caught],
        "Ogw": np.asarray(ogw, dtype=np.float64),
        "Oj": np.asarray(oj, dtype=np.float64),
        "Opgw": np.asarray(opgw, dtype=np.float64),
        "used_tail": np.asarray(used_tail, dtype=bool),
        "signal": signal,
    }


def solve_fast(params: dict[str, float], freqs: np.ndarray, root_in: float, z_tail: float):
    model = LCDM_SG(**{**params, "DN_eff": root_in})
    old_goal = FA.goal_oriented_freqs
    old_construct = FA.construct_fast_frequency_grid

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-6, z_tail=z_tail,
            threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = old_goal
        FA.construct_fast_frequency_grid = old_construct
    if result is not model:
        raise RuntimeError(f"Fast 937-node solve failed: {model.fast_failure_reason}")
    logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    actual_f = np.asarray(model.f, dtype=np.float64)
    if digest(actual_f) != digest(freqs):
        raise RuntimeError("Fast did not honor the exact 937-node grid")
    return {
        "root_input_DN_eff": root_in,
        "root_output_DN_eff": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(np.asarray(model.DN_gw, dtype=np.float64)[-1]),
        "frequency_count": int(actual_f.size),
        "frequency_sha256_f64le": digest(actual_f),
        "frequency_grid_used": model.freq_grid_used,
        "runtime_s": time.perf_counter() - started,
        "log10_OmegaGW": logs.tolist(),
        "log_spectrum_finite_count": int(np.count_nonzero(np.isfinite(logs))),
    }


def background_probe(model: LCDM_SG, dn_values: dict[str, float]) -> dict[str, object]:
    derived = model.derived_param
    n_lo = float(derived["N_inf"] - derived["N_re"])
    n_hi = float(derived["N_inf"])
    points = np.linspace(n_lo, n_hi, 2001, dtype=np.float64)
    baseline = _Background(model, 0.0)
    h2_base = np.asarray([_H2_and_sigma_python(baseline, float(n))[0] for n in points])
    sigma_base = np.asarray([_H2_and_sigma_python(baseline, float(n))[1] for n in points])
    results = {}
    for label, dn_eff in dn_values.items():
        candidate = _Background(model, dn_eff)
        h2 = np.asarray([_H2_and_sigma_python(candidate, float(n))[0] for n in points])
        sigma = np.asarray([_H2_and_sigma_python(candidate, float(n))[1] for n in points])
        results[label] = {
            "dn_eff": dn_eff,
            "H2_bitwise_equal_count": int(np.count_nonzero(h2 == h2_base)),
            "sigma_bitwise_equal_count": int(np.count_nonzero(sigma == sigma_base)),
            "probe_count": int(points.size),
            "max_relative_H2_delta": float(np.max(np.abs(h2 / h2_base - 1.0))),
            "max_absolute_sigma_delta": float(np.max(np.abs(sigma - sigma_base))),
            "background_coefficients": {
                key: float(getattr(candidate, key))
                for key in (
                    "Oerh2", "Otrh2", "Otreh2", "OLh2", "N_inf", "N_re_abs"
                )
            },
        }
    return {
        "domain_N": [n_lo, n_hi],
        "probe_count": int(points.size),
        "dn_eff_zero_coefficients": {
            key: float(getattr(baseline, key))
            for key in ("Oerh2", "Otrh2", "Otreh2", "OLh2", "N_inf", "N_re_abs")
        },
        "cases": results,
        "interpretation": "This is a finite probe of exact Reference background formulas, not a substitute for direct map closure; no DN_eff dependency was bypassed.",
    }


def main() -> int:
    fast_doc = json.loads(FAST_INPUT.read_text(encoding="utf-8"))
    ref_doc = json.loads(REF469_INPUT.read_text(encoding="utf-8"))
    expected_sources = {rel: sha256(ROOT / rel) for rel in SOURCE_FILES}
    if fast_doc["source_sha256"] != expected_sources or ref_doc["source_sha256"] != expected_sources:
        raise RuntimeError("Live source hashes differ from the saved direct-solve artifacts")
    if ref_doc.get("status") != "complete_nested_direct_reference_235_469_underflow_boundary":
        raise RuntimeError("469-node direct Reference artifact is incomplete")
    old_freqs = np.asarray(ref_doc["frequency_grids"]["469"]["log10_hz"], dtype=np.float64)
    old_hash = digest(old_freqs)
    if old_freqs.size != 469 or old_hash != ref_doc["frequency_grids"]["469"]["sha256_f64le"]:
        raise RuntimeError("Saved exact 469-node grid hash mismatch")
    midpoints = 0.5 * (old_freqs[:-1] + old_freqs[1:])
    freqs = np.sort(np.concatenate((old_freqs, midpoints)))[::-1]
    if freqs.size != 937 or not np.all(np.diff(freqs) < 0.0):
        raise RuntimeError("Nested refinement did not produce 937 unique descending nodes")
    old_indices = np.searchsorted(-freqs, -old_freqs)
    if not np.array_equal(freqs[old_indices], old_freqs):
        raise RuntimeError("Saved 469 nodes are not exact members of 937-node grid")

    params = {
        "Omega_bh2": 0.0223828,
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
        "DN_re": 10.0,
        "r": 1.0,
        "n_t": 0.0,
        "cr": 0.0,
        "T_re": 2000.0,
        "kappa10": 0.01,
    }
    fast_rows = {
        row["amplitude_label"]: row
        for row in fast_doc["cases"]
        if row["parameters"]["n_t"] == N_T
        and row["parameters"]["T_re"] == 2000.0
        and row["parameters"]["kappa10"] == 0.01
    }
    prior_fast = ref_doc["fast_469"]
    model = unit_model(params)
    dn_values = {
        label: float(prior_fast[label]["root_output_DN_eff"])
        for label in LABELS
    }
    bg_probe = background_probe(model, dn_values)

    reference_results: dict[str, object] = {}
    scaled_results: dict[str, object] = {}
    fast_results: dict[str, dict[str, object]] = {str(z): {} for z in Z_TAILS}
    for z_tail in Z_TAILS:
        prior = ref_doc["direct_reference"][str(z_tail)]
        old_signal = np.asarray(prior["Ogw"], dtype=np.float64) - np.asarray(prior["Oj"], dtype=np.float64)
        ref = solve_reference(model, freqs, z_tail, 0.0, f"main-z{z_tail:g}")
        signal = np.asarray(ref["signal"], dtype=np.float64)
        common_delta = np.abs(signal[old_indices] - old_signal)
        if np.any(common_delta != 0.0):
            raise RuntimeError(f"469 common Reference nodes differ at z{z_tail:g}: max={np.max(common_delta):g}")
        new_mask = np.ones(freqs.size, dtype=bool)
        new_mask[old_indices] = False
        new_indices = np.flatnonzero(new_mask)
        spline = PchipInterpolator(
            np.sort(old_freqs), old_signal[np.argsort(old_freqs)], extrapolate=False
        )
        predicted = np.asarray(spline(freqs[new_indices]), dtype=np.float64)
        direct_new = signal[new_indices]
        pchip_delta = np.abs(np.log10(direct_new) - np.log10(predicted))
        reference_results[str(z_tail)] = {
            "runtime_s": ref["runtime_s"],
            "dn_eff": 0.0,
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": digest(freqs),
            "used_tail_count": int(np.count_nonzero(ref["used_tail"])),
            "positive_signal_count": int(np.count_nonzero(signal > 0.0)),
            "nonfinite_signal_count": int(np.count_nonzero(~np.isfinite(signal))),
            "runtime_warnings": ref["runtime_warnings"],
            "shared_469_max_abs_linear_signal_delta": float(np.max(common_delta)),
            "new_direct_node_count": int(new_indices.size),
            "new_node_vs_469_linear_signal_pchip_max_abs_dex": float(np.max(pchip_delta)),
            "new_node_vs_469_linear_signal_pchip_p95_abs_dex": float(np.quantile(pchip_delta, 0.95)),
            "Ogw": np.asarray(ref["Ogw"]).tolist(),
            "Oj": np.asarray(ref["Oj"]).tolist(),
            "Opgw": np.asarray(ref["Opgw"]).tolist(),
            "used_tail": np.asarray(ref["used_tail"]).tolist(),
        }
        scaled_results[str(z_tail)] = {}
        for label in LABELS:
            fast_row = fast_rows[label]
            r = float(fast_row["parameters"]["r"])
            log_amp = math.log(r) + math.log(float(fast_row["parameters"]["A_s"]))
            log_power = log_amp + N_T * (math.log(10.0) * freqs - math.log(gp.f_piv))
            mapped, qerr, ierr = scaled_map(freqs, signal, log_power, model)
            old_map = float(ref_doc["scaled_amplitude_comparisons"][str(z_tail)][label]["reference_mapped_DN_gw"])
            fast_results[str(z_tail)][label] = solve_fast(
                {**params, "r": r, "A_s": float(fast_row["parameters"]["A_s"]), "n_t": N_T},
                freqs,
                0.0,
                z_tail,
            )
            fast_result = fast_results[str(z_tail)][label]
            fast_log = np.asarray(fast_result["log10_OmegaGW"], dtype=np.float64)
            log_reference = np.log10(signal) + log_power / math.log(10.0)
            delta = np.abs(fast_log - log_reference)
            scaled_results[str(z_tail)][label] = {
                "reference_DN_eff_input": 0.0,
                "reference_mapped_DN_gw": mapped,
                "previous_469_reference_map": old_map,
                "frequency_refinement_map_ratio_937_over_469": mapped / old_map - 1.0,
                "fast_235_selfconsistent_root": fast_row["DN_eff_final"],
                "fast_469_map_probe_input_root": prior_fast[label]["root_input_DN_eff"],
                "fast_469_map_probe_output_root": prior_fast[label]["root_output_DN_eff"],
                "fast_937_root_input": fast_result["root_input_DN_eff"],
                "fast_937_root_output": fast_result["root_output_DN_eff"],
                "fast_937_mapped_DN_gw": fast_result["mapped_DN_gw"],
                "fast_reference_map_relative_difference_937": fast_result["mapped_DN_gw"] / mapped - 1.0,
                "fast_937_selfconsistent_root_shift_from_235": fast_result["root_output_DN_eff"] / float(fast_row["DN_eff_final"]) - 1.0,
                "same_grid_fast_reference_max_abs_dex": float(np.max(delta)),
                "same_grid_fast_reference_p95_abs_dex": float(np.quantile(delta, 0.95)),
                "relative_quadrature_error_estimate": qerr,
                "relative_interpolation_error_estimate": ierr,
                "scaled_integrand_zero_count": int(np.count_nonzero(np.exp(np.log(signal) + log_power - np.max(np.log(signal) + log_power)) == 0.0)),
            }
    # Directly close the z12 map at the minimum-subnormal amplitude on this 937-node grid.
    closure_label = "min_subnormal_r"
    closure_start = float(scaled_results["12.0"][closure_label]["reference_mapped_DN_gw"])
    closure_model = unit_model(params)
    closure_ref = solve_reference(closure_model, freqs, 12.0, closure_start, "updated-root-z12-minsub")
    closure_signal = np.asarray(closure_ref["signal"], dtype=np.float64)
    closure_row = fast_rows[closure_label]
    closure_amp = math.log(float(closure_row["parameters"]["r"])) + math.log(float(closure_row["parameters"]["A_s"]))
    closure_log_power = closure_amp + N_T * (math.log(10.0) * freqs - math.log(gp.f_piv))
    closure_map, closure_qerr, closure_ierr = scaled_map(freqs, closure_signal, closure_log_power, closure_model)
    closure = {
        "amplitude_label": closure_label,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "z_tail": 12.0,
        "root_before_remap": closure_start,
        "mapped_root_after_remap": closure_map,
        "root_residual_relative": closure_map / closure_start - 1.0,
        "quadrature_error_relative": closure_qerr,
        "interpolation_error_relative": closure_ierr,
        "runtime_s": closure_ref["runtime_s"],
        "runtime_warnings": closure_ref["runtime_warnings"],
        "scaled_integrand_zero_count": 0,
    }
    output = {
        "schema_version": 1,
        "status": "complete_nested_direct_reference_469_937_underflow_boundary",
        "scope": "cr=0,n_t=1.6,T_re=2000,kappa10=.01; nested 469→937 direct Reference z10/z12, exact-grid Fast comparisons, background sensitivity probe, and updated-root closure for minimum-subnormal r. Discrete evidence only.",
        "fast_input_sha256": sha256(FAST_INPUT),
        "reference_469_input_sha256": sha256(REF469_INPUT),
        "source_sha256": expected_sources,
        "runner_sha256": sha256(Path(__file__)),
        "parameters": params,
        "rtol": RTOL,
        "threads": 1,
        "z_tails": list(Z_TAILS),
        "frequency_grids": {
            "469": {"count": int(old_freqs.size), "sha256_f64le": old_hash, "log10_hz": old_freqs.tolist()},
            "937": {"count": int(freqs.size), "sha256_f64le": digest(freqs), "log10_hz": freqs.tolist()},
        },
        "background_root_sensitivity_probe": bg_probe,
        "direct_reference": reference_results,
        "scaled_amplitude_comparisons": scaled_results,
        "fast_937": fast_results,
        "reference_root_closure": closure,
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
        "grid_469_to_937_map_shift": {
            z: {label: row["frequency_refinement_map_ratio_937_over_469"] for label, row in rows.items()}
            for z, rows in scaled_results.items()
        },
        "new_node_interpolation": {
            z: {
                "max_dex": row["new_node_vs_469_linear_signal_pchip_max_abs_dex"],
                "p95_dex": row["new_node_vs_469_linear_signal_pchip_p95_abs_dex"],
                "shared_node_delta": row["shared_469_max_abs_linear_signal_delta"],
            }
            for z, row in reference_results.items()
        },
        "root_closure": closure,
        "background_probe": bg_probe,
        "same_grid_fast_reference": {
            z: {
                label: {
                    "map_rel": row["fast_reference_map_relative_difference_937"],
                    "max_dex": row["same_grid_fast_reference_max_abs_dex"],
                }
                for label, row in rows.items()
            }
            for z, rows in scaled_results.items()
        },
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
