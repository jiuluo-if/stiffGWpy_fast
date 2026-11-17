"""English description: Refine low-tilt underflow Reference from 235 to nested 469 direct nodes. 中文说明：将低倾斜下溢 Reference 从235点加密到嵌套469个直接求解节点。"""
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
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

FAST_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_fast_screen_20261117.json"
REF_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_reference_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_reference_235_469_20261117.json"
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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def scaled_map(freqs: np.ndarray, signal: np.ndarray, log_power: np.ndarray, model: LCDM_SG):
    if np.any(signal <= 0.0) or not np.isfinite(signal).all():
        raise FloatingPointError("Reference O_gw-O_j must stay positive and finite")
    log_integrand = np.log(signal) + log_power
    shift = float(np.max(log_integrand))
    normalized = np.exp(log_integrand - shift)
    g2scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    if g2scaled <= 0.0 or not math.isfinite(float(g2scaled)):
        raise FloatingPointError("Reference normalized integral is not positive finite")
    logg2 = shift + math.log(float(g2scaled))
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    mapped = gp.Neff0 / omega_nu * math.exp(logg2)
    return float(mapped), float(qerr / g2scaled) if qerr is not None else None, float(ierr / g2scaled) if ierr is not None else None


def reference_grid(params, freqs: np.ndarray, z_tail: float, dn_eff: float):
    unit_params = {**params, "r": 1.0, "A_s": 1.0, "n_t": 0.0, "DN_eff": 0.0}
    model = LCDM_SG(**unit_params)
    model.gen_expansion()
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        ogw, oj, opgw, used_tail = spectrum_reference(
            model, freqs, dn_eff, z_tail=z_tail, rtol=RTOL, workers=1, backend="cartesian"
        )
    return {
        "runtime_s": time.perf_counter() - started,
        "dn_eff": dn_eff,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "Ogw": np.asarray(ogw, dtype=np.float64),
        "Oj": np.asarray(oj, dtype=np.float64),
        "Opgw": np.asarray(opgw, dtype=np.float64),
        "used_tail": np.asarray(used_tail, dtype=bool),
        "runtime_warnings": [str(item.message) for item in caught],
        "model": model,
    }


def run_fast_on_grid(params: dict[str, float], freqs: np.ndarray, root_in: float):
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
            engine="fast", fallback=False, tol=1e-6, z_tail=12.0,
            threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = old_goal
        FA.construct_fast_frequency_grid = old_construct
    if result is not model:
        raise RuntimeError(f"Fast 469-node solve failed: {model.fast_failure_reason}")
    logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    return {
        "runtime_s": time.perf_counter() - started,
        "root_input_DN_eff": root_in,
        "root_output_DN_eff": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(np.asarray(model.DN_gw, dtype=np.float64)[-1]),
        "frequency_count": int(model.f.size),
        "frequency_sha256_f64le": digest(np.asarray(model.f, dtype=np.float64)),
        "frequency_grid_used": model.freq_grid_used,
        "log10_OmegaGW": logs.tolist(),
        "log_spectrum_finite_count": int(np.count_nonzero(np.isfinite(logs))),
        "runtime_warnings": [],
    }


def main() -> int:
    fast_doc = json.loads(FAST_INPUT.read_text(encoding="utf-8"))
    ref_doc = json.loads(REF_INPUT.read_text(encoding="utf-8"))
    expected_sources = {rel: sha256(ROOT / rel) for rel in SOURCE_FILES}
    if fast_doc.get("source_sha256") != expected_sources or ref_doc.get("source_sha256") != expected_sources:
        raise RuntimeError("Fast/Reference artifact source hashes are stale")
    if ref_doc.get("status") != "complete_direct_reference_scaled_underflow_boundary_screen":
        raise RuntimeError("235-node direct Reference artifact is incomplete")
    fast_rows = {
        row["amplitude_label"]: row
        for row in fast_doc["cases"]
        if row["parameters"]["n_t"] == N_T
        and row["parameters"]["T_re"] == 2000.0
        and row["parameters"]["kappa10"] == 0.01
    }
    labels = ("normal_edge", "subnormal_edge", "subnormal_interior", "rounded_zero_product", "min_subnormal_r")
    if any(label not in fast_rows for label in labels):
        raise RuntimeError("One or more amplitude strata are missing from the Fast screen")
    old_key = next(key for key, values in ref_doc["frequency_grids_log10_hz"].items() if len(values) == 235)
    old_freqs = np.asarray(ref_doc["frequency_grids_log10_hz"][old_key], dtype=np.float64)
    old_freq_hash = digest(old_freqs)
    if old_freq_hash != old_key:
        raise RuntimeError("Saved 235-node frequency hash does not match its exact nodes")
    midpoints = 0.5 * (old_freqs[:-1] + old_freqs[1:])
    freqs = np.sort(np.concatenate((old_freqs, midpoints)))[::-1]
    if freqs.size != 469 or not np.all(np.diff(freqs) < 0.0):
        raise RuntimeError("Nested refinement did not produce 469 unique descending nodes")
    old_indices = np.searchsorted(-freqs, -old_freqs)
    if not np.array_equal(freqs[old_indices], old_freqs):
        raise RuntimeError("Original 235 nodes are not exact members of the 469-node grid")

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
    reference_235_z12 = reference_grid(params, old_freqs, 12.0, 0.0)
    reference_235_z12_signal = (
        np.asarray(reference_235_z12["Ogw"]) - np.asarray(reference_235_z12["Oj"])
    )
    if np.any(reference_235_z12_signal <= 0.0) or not np.isfinite(reference_235_z12_signal).all():
        raise FloatingPointError("235-node z12 Reference signal must be positive and finite")
    raw: dict[str, dict[str, object]] = {}
    comparisons: dict[str, dict[str, object]] = {}
    fast469: dict[str, object] = {}
    for z_tail in Z_TAILS:
        ref = reference_grid(params, freqs, z_tail, 0.0)
        ogw = np.asarray(ref["Ogw"])
        oj = np.asarray(ref["Oj"])
        signal = ogw - oj
        if np.any(signal <= 0.0) or not np.isfinite(signal).all():
            raise FloatingPointError(f"Nonpositive or nonfinite Reference signal at z_tail={z_tail}")
        if z_tail == 12.0:
            old_signal = reference_235_z12_signal
        else:
            previous = ref_doc["reference_unit_transfer"][old_key][str(z_tail)]
            old_signal = np.asarray(previous["Ogw"], dtype=np.float64) - np.asarray(previous["Oj"], dtype=np.float64)
        common_error = float(np.max(np.abs(signal[old_indices] - old_signal)))
        if common_error != 0.0:
            raise RuntimeError(f"Shared direct Reference nodes changed by {common_error:g} at z_tail={z_tail}")
        midpoint_mask = np.ones(freqs.size, dtype=bool)
        midpoint_mask[old_indices] = False
        new_idx = np.flatnonzero(midpoint_mask)
        old_spline = PchipInterpolator(np.sort(old_freqs), old_signal[np.argsort(old_freqs)], extrapolate=False)
        predicted = np.asarray(old_spline(freqs[new_idx]), dtype=np.float64)
        direct = signal[new_idx]
        interpolation_delta = np.abs(np.log10(direct) - np.log10(predicted))
        grid_results: dict[str, object] = {
            "runtime_s": ref["runtime_s"],
            "dn_eff": 0.0,
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": digest(freqs),
            "used_tail_count": int(np.count_nonzero(ref["used_tail"])),
            "positive_signal_count": int(np.count_nonzero(signal > 0.0)),
            "nonfinite_signal_count": int(np.count_nonzero(~np.isfinite(signal))),
            "runtime_warnings": ref["runtime_warnings"],
            "common_235_max_abs_linear_signal_delta": common_error,
            "new_direct_node_count": int(new_idx.size),
            "new_node_vs_235_linear_signal_pchip_max_abs_dex": float(np.max(interpolation_delta)),
            "new_node_vs_235_linear_signal_pchip_p95_abs_dex": float(np.quantile(interpolation_delta, 0.95)),
            "Ogw": ogw.tolist(),
            "Oj": oj.tolist(),
            "Opgw": np.asarray(ref["Opgw"]).tolist(),
            "used_tail": np.asarray(ref["used_tail"]).tolist(),
        }
        raw[str(z_tail)] = grid_results
        comparisons[str(z_tail)] = {}
        for label in labels:
            fast_row = fast_rows[label]
            if not fast_row["fast_returned_model"]:
                raise RuntimeError(f"Selected Fast case failed its guard: {label}")
            r = float(fast_row["parameters"]["r"])
            log_amp = math.log(r) + math.log(float(fast_row["parameters"]["A_s"]))
            log_power = log_amp + N_T * (math.log(10.0) * freqs - math.log(gp.f_piv))
            mapped, qerr, ierr = scaled_map(freqs, signal, log_power, ref["model"])
            old_mapped, _, _ = scaled_map(old_freqs, old_signal, log_amp + N_T * (
                math.log(10.0) * old_freqs - math.log(gp.f_piv)
            ), ref["model"])
            comparisons[str(z_tail)][label] = {
                "reference_DN_eff_input": 0.0,
                "reference_mapped_DN_gw": mapped,
                "fast_235_root_DN_eff": fast_row["DN_eff_final"],
                "reference_root_from_same_start": mapped,
                "reference_235_mapped_DN_gw": old_mapped,
                "frequency_refinement_map_ratio_469_over_235": mapped / old_mapped - 1.0,
                "relative_quadrature_error_estimate": qerr,
                "relative_interpolation_error_estimate": ierr,
                "scaled_integrand_zero_count": 0,
            }
            # Same exact 469 grid and tail for direct Fast/Reference comparison.
            if label not in fast469:
                fast469[label] = run_fast_on_grid(
                    {**params, "r": r, "A_s": float(fast_row["parameters"]["A_s"]), "n_t": N_T},
                    freqs,
                    float(fast_row["DN_eff_final"]),
                )
            fast_result = fast469[label]
            if fast_result["frequency_sha256_f64le"] != digest(freqs):
                raise RuntimeError(f"Fast did not honor the exact 469-node grid for {label}")
            fast_logs = np.asarray(fast_result["log10_OmegaGW"], dtype=np.float64)
            log_reference = np.log10(signal) + log_power / math.log(10.0)
            delta = np.abs(fast_logs - log_reference)
            comparisons[str(z_tail)][label].update({
                "fast_469_map_DN_gw": fast_result["mapped_DN_gw"],
                "fast_reference_map_relative_difference_469": float(fast_result["mapped_DN_gw"]) / mapped - 1.0,
                "same_grid_fast_reference_max_abs_dex": float(np.max(delta)),
                "same_grid_fast_reference_p95_abs_dex": float(np.quantile(delta, 0.95)),
                "fast_469_runtime_s": fast_result["runtime_s"],
            })
        ref["model"] = None

    # Close the z12 root at the minimum-subnormal point by one fresh direct 469-node map.
    closure_label = "min_subnormal_r"
    closure_root = float(comparisons["12.0"][closure_label]["reference_mapped_DN_gw"])
    closure_ref = reference_grid(params, freqs, 12.0, closure_root)
    closure_signal = np.asarray(closure_ref["Ogw"]) - np.asarray(closure_ref["Oj"])
    closure_fast_row = fast_rows[closure_label]
    closure_log_amp = math.log(float(closure_fast_row["parameters"]["r"])) + math.log(
        float(closure_fast_row["parameters"]["A_s"])
    )
    closure_log_power = closure_log_amp + N_T * (
        math.log(10.0) * freqs - math.log(gp.f_piv)
    )
    closure_mapped, closure_qerr, closure_ierr = scaled_map(
        freqs, closure_signal, closure_log_power, closure_ref["model"]
    )
    closure = {
        "amplitude_label": closure_label,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "z_tail": 12.0,
        "root_before_remap": closure_root,
        "mapped_root_after_remap": closure_mapped,
        "root_residual_relative": closure_mapped / closure_root - 1.0,
        "quadrature_error_relative": closure_qerr,
        "interpolation_error_relative": closure_ierr,
        "runtime_s": closure_ref["runtime_s"],
        "runtime_warnings": closure_ref["runtime_warnings"],
        "scaled_integrand_zero_count": 0,
    }
    output = {
        "schema_version": 1,
        "status": "complete_nested_direct_reference_235_469_underflow_boundary",
        "scope": "cr=0,n_t=1.6,T_re=2000,kappa10=.01; nested 235→469 direct Reference grids at z10/z12, matched exact-grid Fast checks, and z12 updated-root closure for minimum-subnormal r. Not a continuous-domain claim.",
        "fast_input_sha256": sha256(FAST_INPUT),
        "reference_235_input_sha256": sha256(REF_INPUT),
        "source_sha256": expected_sources,
        "runner_sha256": sha256(Path(__file__)),
        "parameters": params,
        "rtol": RTOL,
        "threads": 1,
        "z_tails": list(Z_TAILS),
        "frequency_grids": {
            "235": {"count": int(old_freqs.size), "sha256_f64le": old_freq_hash, "log10_hz": old_freqs.tolist()},
            "469": {"count": int(freqs.size), "sha256_f64le": digest(freqs), "log10_hz": freqs.tolist()},
        },
        "direct_reference": raw,
        "direct_reference_235_z12": {
            "runtime_s": reference_235_z12["runtime_s"],
            "frequency_count": int(old_freqs.size),
            "frequency_sha256_f64le": digest(old_freqs),
            "used_tail_count": int(np.count_nonzero(reference_235_z12["used_tail"])),
            "runtime_warnings": reference_235_z12["runtime_warnings"],
            "Ogw": np.asarray(reference_235_z12["Ogw"]).tolist(),
            "Oj": np.asarray(reference_235_z12["Oj"]).tolist(),
            "Opgw": np.asarray(reference_235_z12["Opgw"]).tolist(),
            "used_tail": np.asarray(reference_235_z12["used_tail"]).tolist(),
        },
        "scaled_amplitude_comparisons": comparisons,
        "fast_469": fast469,
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
        "grid_235_to_469_map_shift": {
            z: {label: value["frequency_refinement_map_ratio_469_over_235"] for label, value in rows.items()}
            for z, rows in comparisons.items()
        },
        "new_node_interpolation": {
            z: {
                "max_dex": row["new_node_vs_235_linear_signal_pchip_max_abs_dex"],
                "p95_dex": row["new_node_vs_235_linear_signal_pchip_p95_abs_dex"],
                "common_node_delta": row["common_235_max_abs_linear_signal_delta"],
            }
            for z, row in raw.items()
        },
        "root_closure": closure,
        "same_grid_fast_reference": {
            z: {
                label: {
                    "map_rel": value["fast_reference_map_relative_difference_469"],
                    "max_dex": value["same_grid_fast_reference_max_abs_dex"],
                }
                for label, value in rows.items()
            }
            for z, rows in comparisons.items()
        },
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
