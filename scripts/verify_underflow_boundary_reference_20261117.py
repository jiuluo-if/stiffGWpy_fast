"""English description: Compare boundary-stratum Fast spectra with direct scaled Reference solves. 中文说明：将边界振幅分层 Fast 谱与直接单位传递 Reference 结果对照。"""
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

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import integrate_spectrum, spectrum_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

FAST_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_fast_screen_20261117.json"
PREFLIGHT_INPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_preflight_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_boundary_reference_20261117.json"
SOURCE_FILES = (
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/reference.py",
    "stiffgwpy_fast/stiff_SGWB.py",
    "stiffgwpy_fast/LCDM_stiff_Neff.py",
)
TARGET_N_T = 1.6
TARGET_T_RE = 2000.0
TARGET_KAPPA = 0.01
TILTS = (TARGET_N_T,)
Z_TAILS = (8.0, 10.0)
RTOL = 1e-9


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def finite_json(values: np.ndarray) -> list[float | None]:
    return [float(value) if math.isfinite(float(value)) else None for value in values]


def scaled_map(
    freqs: np.ndarray, base: np.ndarray, log_power: np.ndarray, model: LCDM_SG
) -> tuple[float, np.ndarray, float | None, float | None, int]:
    if not np.isfinite(base).all() or np.any(base <= 0.0):
        raise FloatingPointError("Direct Reference signal must be positive and finite on this grid")
    log_integrand = np.log(base) + log_power
    shift = float(np.max(log_integrand))
    normalized = np.exp(log_integrand - shift)
    g2_scaled, qerr, ierr = integrate_spectrum(freqs, normalized, np.zeros_like(normalized))
    if g2_scaled <= 0.0 or not math.isfinite(float(g2_scaled)):
        raise FloatingPointError("Scaled Reference integral is not positive and finite")
    log_g2 = shift + math.log(float(g2_scaled))
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"]) ** 2
    mapped = gp.Neff0 / omega_nu * math.exp(log_g2)
    return (
        float(mapped),
        log_integrand / math.log(10.0),
        None if qerr is None else float(qerr / g2_scaled),
        None if ierr is None else float(ierr / g2_scaled),
        int(np.count_nonzero(normalized == 0.0)),
    )


def solve_reference_grid(
    params: dict[str, float], freqs: np.ndarray, z_tail: float, dn_eff: float
) -> dict[str, object]:
    unit_params = {
        **params,
        "r": 1.0,
        "A_s": 1.0,
        "n_t": 0.0,
        "DN_eff": 0.0,
    }
    model = LCDM_SG(**unit_params)
    model.gen_expansion()
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        ogw, oj, opgw, used_tail = spectrum_reference(
            model,
            freqs,
            dn_eff,
            z_tail=z_tail,
            rtol=RTOL,
            workers=1,
            backend="cartesian",
        )
    elapsed = time.perf_counter() - started
    signal = np.asarray(ogw - oj, dtype=np.float64)
    return {
        "runtime_s": elapsed,
        "dn_eff": dn_eff,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "used_tail_count": int(np.count_nonzero(used_tail)),
        "positive_signal_count": int(np.count_nonzero(signal > 0.0)),
        "nonfinite_signal_count": int(np.count_nonzero(~np.isfinite(signal))),
        "signal_min": float(np.min(signal)),
        "signal_max": float(np.max(signal)),
        "log10_signal": finite_json(np.log10(signal)) if np.all(signal > 0.0) else None,
        "runtime_warnings": [str(item.message) for item in caught],
        "Ogw": ogw.tolist(),
        "Oj": oj.tolist(),
        "Opgw": opgw.tolist(),
        "used_tail": used_tail.tolist(),
        "model_for_scaling": model,
    }


def main() -> int:
    fast_doc = json.loads(FAST_INPUT.read_text(encoding="utf-8"))
    preflight = json.loads(PREFLIGHT_INPUT.read_text(encoding="utf-8"))
    if fast_doc.get("status") != "complete_fast_guarded_amplitude_strata_screen":
        raise RuntimeError("Fast amplitude screen is incomplete")
    if preflight.get("status") != "complete":
        raise RuntimeError("Frequency/background preflight is incomplete")
    expected_sources = {rel: sha256(ROOT / rel) for rel in SOURCE_FILES}
    if fast_doc["source_sha256"] != expected_sources or preflight["source_sha256"] != expected_sources:
        raise RuntimeError("Live source hashes differ from the Fast screen or preflight")
    fast_rows = {
        row["amplitude_label"]: row
        for row in fast_doc["cases"]
        if row["parameters"]["n_t"] == TARGET_N_T
        and row["parameters"]["T_re"] == TARGET_T_RE
        and row["parameters"]["kappa10"] == TARGET_KAPPA
    }
    if len(fast_rows) != 6:
        raise RuntimeError("Expected all six amplitude strata in Fast screen")
    amplitude_labels = tuple(fast_rows)
    grids: dict[str, np.ndarray] = {}
    for label, row in fast_rows.items():
        if not row["fast_returned_model"]:
            raise RuntimeError(f"Chosen Reference stratum was not admitted by Fast: {label}")
        freqs = np.asarray(row["frequency_log10_hz"], dtype=np.float64)
        if digest(freqs) != row["frequency_sha256_f64le"]:
            raise RuntimeError(f"Fast raw-grid hash mismatch for {label}")
        key = row["frequency_sha256_f64le"]
        grids.setdefault(key, freqs)
    if len(grids) != 2:
        raise RuntimeError(f"Expected two unique exact Fast grids, found {len(grids)}")

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
        "T_re": TARGET_T_RE,
        "kappa10": TARGET_KAPPA,
    }
    # Key retained separately so the unit Reference model remains visibly cr=0.
    grid_reference: dict[str, dict[str, object]] = {}
    scale_rows: dict[str, dict[str, dict[str, object]]] = {}
    for grid_hash, freqs in grids.items():
        grid_reference[grid_hash] = {}
        scale_rows[grid_hash] = {}
        for z_tail in Z_TAILS:
            ref = solve_reference_grid(params, freqs, z_tail, 0.0)
            signal = np.asarray(ref["Ogw"], dtype=np.float64) - np.asarray(ref["Oj"], dtype=np.float64)
            if np.any(signal <= 0.0) or not np.isfinite(signal).all():
                raise FloatingPointError(f"Reference signal invalid for {grid_hash} at z_tail={z_tail}")
            scale_rows[grid_hash][str(z_tail)] = {}
            for label in amplitude_labels:
                fast_row = fast_rows[label]
                r = float(fast_row["parameters"]["r"])
                if r == 0.0:
                    continue
                log_amplitude = math.log(r) + math.log(float(fast_row["parameters"]["A_s"]))
                log_power = log_amplitude + TARGET_N_T * (
                    math.log(10.0) * freqs - math.log(gp.f_piv)
                )
                mapped, log_spectrum, qerr_rel, ierr_rel, zeros = scaled_map(
                    freqs, signal, log_power, ref["model_for_scaling"]
                )
                same_grid = fast_row["frequency_sha256_f64le"] == grid_hash
                if same_grid:
                    fast_logs = np.asarray(fast_row["log10_OmegaGW"], dtype=np.float64)
                    delta = np.abs(fast_logs - log_spectrum)
                    max_spectrum_error = float(np.max(delta))
                    p95_spectrum_error = float(np.quantile(delta, 0.95))
                    fast_map_error = float(fast_row["DN_eff_final"]) / mapped - 1.0
                else:
                    max_spectrum_error = None
                    p95_spectrum_error = None
                    fast_map_error = None
                scale_rows[grid_hash][str(z_tail)][label] = {
                    "mapped_DN_gw": mapped,
                    "fast_DN_eff_root": fast_row["DN_eff_final"],
                    "same_frequency_grid_as_fast": same_grid,
                    "fast_reference_map_relative_difference": fast_map_error,
                    "max_spectrum_error_dex": max_spectrum_error,
                    "p95_spectrum_error_dex": p95_spectrum_error,
                    "quadrature_error_relative": qerr_rel,
                    "interpolation_error_relative": ierr_rel,
                    "scaled_integrand_zero_count": zeros,
                    "reference_log10_OmegaGW": log_spectrum.tolist(),
                }
            grid_reference[grid_hash][str(z_tail)] = {
                key: value for key, value in ref.items() if key != "model_for_scaling"
            }
            del ref

    grid_hashes = list(grids)
    coarse_hash = min(grid_hashes, key=lambda key: grids[key].size)
    fine_hash = max(grid_hashes, key=lambda key: grids[key].size)
    coarse_x = grids[coarse_hash]
    fine_x = grids[fine_hash]
    coarse_x_asc_idx = np.argsort(coarse_x)
    fine_x_asc_idx = np.argsort(fine_x)
    x_coarse = coarse_x[coarse_x_asc_idx]
    x_fine = fine_x[fine_x_asc_idx]
    overlap = (x_fine >= x_coarse[0]) & (x_fine <= x_coarse[-1])
    interpolation: dict[str, object] = {
        "coarse_grid_count": int(coarse_x.size),
        "fine_grid_count": int(fine_x.size),
        "overlap_fine_node_count": int(np.count_nonzero(overlap)),
        "coarse_grid_sha256": coarse_hash,
        "fine_grid_sha256": fine_hash,
        "interpolated_error_by_tail_and_amplitude": {},
        "fixed_root_map_difference_by_tail_and_amplitude": {},
    }
    for z_tail in Z_TAILS:
        zkey = str(z_tail)
        interpolation["interpolated_error_by_tail_and_amplitude"][zkey] = {}
        interpolation["fixed_root_map_difference_by_tail_and_amplitude"][zkey] = {}
        coarse_signal = np.asarray(grid_reference[coarse_hash][zkey]["Ogw"]) - np.asarray(
            grid_reference[coarse_hash][zkey]["Oj"]
        )
        fine_signal = np.asarray(grid_reference[fine_hash][zkey]["Ogw"]) - np.asarray(
            grid_reference[fine_hash][zkey]["Oj"]
        )
        spline = PchipInterpolator(x_coarse, coarse_signal[coarse_x_asc_idx], extrapolate=False)
        predicted = np.asarray(spline(x_fine[overlap]), dtype=np.float64)
        actual = fine_signal[fine_x_asc_idx][overlap]
        for label in amplitude_labels:
            fast_row = fast_rows[label]
            if float(fast_row["parameters"]["r"]) == 0.0:
                continue
            log_amp = math.log(float(fast_row["parameters"]["r"])) + math.log(float(fast_row["parameters"]["A_s"]))
            logpt_fine = log_amp + TARGET_N_T * (
                math.log(10.0) * x_fine[overlap] - math.log(gp.f_piv)
            )
            valid = (predicted > 0.0) & (actual > 0.0)
            actual_log = np.log10(actual[valid]) + logpt_fine[valid] / math.log(10.0)
            predicted_log = np.log10(predicted[valid]) + logpt_fine[valid] / math.log(10.0)
            d = np.abs(predicted_log - actual_log)
            coarse_map = scale_rows[coarse_hash][zkey][label]["mapped_DN_gw"]
            fine_map = scale_rows[fine_hash][zkey][label]["mapped_DN_gw"]
            interpolation["interpolated_error_by_tail_and_amplitude"][zkey][label] = {
                "valid_overlap_node_count": int(d.size),
                "max_abs_dex": float(np.max(d)) if d.size else None,
                "p95_abs_dex": float(np.quantile(d, 0.95)) if d.size else None,
            }
            interpolation["fixed_root_map_difference_by_tail_and_amplitude"][zkey][label] = (
                float(fine_map) / float(coarse_map) - 1.0
            )

    # Directly re-map the most amplitude-sensitive accepted case at its z10 updated root.
    closure_label = "min_subnormal_r"
    closure_fast = fast_rows[closure_label]
    closure_hash = closure_fast["frequency_sha256_f64le"]
    closure_freqs = grids[closure_hash]
    closure_first = scale_rows[closure_hash]["10.0"][closure_label]
    closure_root = float(closure_first["mapped_DN_gw"])
    closure_reference = solve_reference_grid(params, closure_freqs, 10.0, closure_root)
    closure_signal = np.asarray(closure_reference["Ogw"]) - np.asarray(closure_reference["Oj"])
    closure_r = float(closure_fast["parameters"]["r"])
    closure_log_amp = math.log(closure_r) + math.log(float(closure_fast["parameters"]["A_s"]))
    closure_logpt = closure_log_amp + TARGET_N_T * (
        math.log(10.0) * closure_freqs - math.log(gp.f_piv)
    )
    closure_map, _, closure_qerr, closure_ierr, closure_zero = scaled_map(
        closure_freqs, closure_signal, closure_logpt, closure_reference["model_for_scaling"]
    )
    reference_root_closure = {
        "amplitude_label": closure_label,
        "frequency_count": int(closure_freqs.size),
        "z_tail": 10.0,
        "root_before_remap": closure_root,
        "mapped_root_after_remap": closure_map,
        "root_residual_relative": closure_map / closure_root - 1.0,
        "quadrature_error_relative": closure_qerr,
        "interpolation_error_relative": closure_ierr,
        "scaled_integrand_zero_count": closure_zero,
        "second_map_runtime_s": closure_reference["runtime_s"],
        "runtime_warnings": closure_reference["runtime_warnings"],
    }

    output = {
        "schema_version": 1,
        "status": "complete_direct_reference_scaled_underflow_boundary_screen",
        "scope": "cr=0,n_t=1.6,T_re=2000,kappa10=.01; exact current Fast 75-node goal and 235-node construct grids; direct Cartesian Reference unit-transfer solves at z8/z10, with log-power scaling for each nonzero amplitude. No continuous-domain claim.",
        "fast_input_sha256": sha256(FAST_INPUT),
        "preflight_input_sha256": sha256(PREFLIGHT_INPUT),
        "parameters": params,
        "rtol": RTOL,
        "z_tails": list(Z_TAILS),
        "source_sha256": expected_sources,
        "runner_sha256": sha256(Path(__file__)),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "frequency_grids_log10_hz": {key: value.tolist() for key, value in grids.items()},
        "reference_unit_transfer": grid_reference,
        "fast_reference_comparisons": scale_rows,
        "frequency_grid_comparison": interpolation,
        "reference_root_closure": reference_root_closure,
    }
    # Internal model objects are removed from serialized reference records above.
    OUTPUT.write_text(json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    compact = {
        "output": str(OUTPUT),
        "grid_counts": [
            {"count": int(value.size), "sha256_f64le": key}
            for key, value in grids.items()
        ],
        "root_closure": reference_root_closure,
        "per_tail_amplitude": {
            z: {
                label: {
                    "map_rel": value["fast_reference_map_relative_difference"],
                    "max_dex": value["max_spectrum_error_dex"],
                    "p95_dex": value["p95_spectrum_error_dex"],
                }
                for label, value in rows.items()
            }
            for grid_hash, tails in scale_rows.items()
            for z, rows in tails.items()
        },
        "frequency_grid": interpolation,
    }
    print(json.dumps(compact, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
