"""Run a fixed Sobol preflight and direct-Reference holdouts on cr>0 cutoff points."""

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
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_holdout_20261116.json"
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
SEED = 20261116
SAMPLE_COUNT = 8
LOG10_T_BOUNDS = (1.0, 4.0)
LOG10_KAPPA_BOUNDS = (-3.0, 0.0)
CUTOFF_RATIO = 1.1
REFERENCE_INDICES = (0,)
REFERENCE_TAILS = (8.0, 10.0, 12.0)
REFERENCE_RTOL = 1e-9


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    paths = {
        "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
        "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
        "exact_background.py": "stiffgwpy_fast/exact_background.py",
        "reference.py": "stiffgwpy_fast/reference.py",
        "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
        "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
        "global_param.py": "stiffgwpy_fast/global_param.py",
        "script": "scripts/screen_cr_positive_sobol_cutoff_holdout_20261116.py",
    }
    return {name: sha(ROOT / path) for name, path in paths.items()}


def cutoff_root(params: dict[str, float]) -> tuple[float, dict[str, object]]:
    probe = LCDM_SG(**{**params, "r": float(np.nextafter(0.0, 1.0))})
    derived = probe.derived_param
    log_r_cutoff = (
        -4.0 * (gp.N_10 - float(derived["N_re"]))
        - 4.0 * math.log(gp.M_Pl)
        - math.log(45.0 / 2.0)
        - math.log(params["A_s"])
        + math.log(probe.rhorad_re + probe.rhostiff_re)
    )
    return math.exp(log_r_cutoff), derived


def reference_summary(result: dict[str, object], runtime_s: float, caught) -> dict[str, object]:
    raw_signal = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
    valid = np.isfinite(raw_signal) & (raw_signal > 0.0)
    raw_logs: list[float | None] = [
        float(math.log10(value)) if okay else None
        for value, okay in zip(raw_signal, valid)
    ]
    return {
        "DN_gw": float(result["DN_gw"]),
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "used_tail_fraction": float(np.mean(np.asarray(result["used_tail"], dtype=bool))),
        "runtime_s": runtime_s,
        "raw_signal_positive_finite_nodes": int(np.count_nonzero(valid)),
        "raw_signal_invalid_nodes": int(raw_signal.size - np.count_nonzero(valid)),
        "public_log_floor_nodes": int(np.count_nonzero(np.asarray(result["log10OmegaGW"]) == -40.0)),
        "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught],
        "raw_log10_OmegaGW": raw_logs,
    }


def compare_fast_to_reference(fast_logs: np.ndarray, ref: dict[str, object], freqs: np.ndarray) -> dict[str, object]:
    ref_values = ref["raw_log10_OmegaGW"]
    valid = np.asarray([value is not None for value in ref_values], dtype=bool) & np.isfinite(fast_logs)
    if not valid.any():
        return {
            "compared_nodes": 0,
            "uncompared_nodes": int(freqs.size),
            "max_abs_dex": None,
            "p95_abs_dex": None,
            "worst_log10f": None,
        }
    ref_logs = np.asarray([0.0 if value is None else float(value) for value in ref_values], dtype=np.float64)
    delta = fast_logs[valid] - ref_logs[valid]
    worst = int(np.argmax(np.abs(delta)))
    common = np.flatnonzero(valid)
    return {
        "compared_nodes": int(np.count_nonzero(valid)),
        "uncompared_nodes": int(freqs.size - np.count_nonzero(valid)),
        "max_abs_dex": float(np.max(np.abs(delta))),
        "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
        "worst_log10f": float(freqs[common[worst]]),
    }


def build_design() -> tuple[np.ndarray, list[dict[str, float]]]:
    sampler = qmc.Sobol(d=2, scramble=True, seed=SEED)
    unit = sampler.random_base2(m=int(math.log2(SAMPLE_COUNT)))
    points = []
    for value in unit:
        log_t = LOG10_T_BOUNDS[0] + float(value[0]) * (LOG10_T_BOUNDS[1] - LOG10_T_BOUNDS[0])
        log_kappa = LOG10_KAPPA_BOUNDS[0] + float(value[1]) * (LOG10_KAPPA_BOUNDS[1] - LOG10_KAPPA_BOUNDS[0])
        points.append({"T_re": 10.0 ** log_t, "kappa10": 10.0 ** log_kappa})
    return unit, points


def main() -> int:
    unit_points, settings = build_design()
    output: dict[str, object] = {
        "schema_version": 1,
        "status": "running_sobol_preflight",
        "scope": "Fixed scrambled Sobol design over log10(T_re) and log10(kappa10) for cr>0 valid-cutoff points. Compute the analytic cutoff using the exact consistency relation, set r=1.1*r_cut, run physical-guarded current Fast on all points, and run direct Reference z8/z10/z12 only for the predeclared first two admitted points. No invalid-cutoff sample is propagated and the Reference source is unchanged.",
        "design": {
            "sampler": "scipy.stats.qmc.Sobol",
            "scramble": True,
            "seed": SEED,
            "sample_count": SAMPLE_COUNT,
            "unit_points": unit_points.tolist(),
            "log10_T_re_bounds": list(LOG10_T_BOUNDS),
            "log10_kappa10_bounds": list(LOG10_KAPPA_BOUNDS),
            "r_over_cutoff": CUTOFF_RATIO,
            "direct_reference_indices_predeclared": list(REFERENCE_INDICES),
            "direct_reference_z_tails": list(REFERENCE_TAILS),
            "reference_rtol": REFERENCE_RTOL,
        },
        "base_parameters": BASE,
        "samples": [],
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": source_hashes(),
    }

    def checkpoint() -> None:
        OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    for index, settings_i in enumerate(settings):
        params = {**BASE, **settings_i}
        r_cut, _ = cutoff_root(params)
        params["r"] = CUTOFF_RATIO * r_cut
        model = LCDM_SG(**params)
        derived = model.derived_param
        initial_delta = None if derived["N_inf"] is None else float(derived["N_inf"] - derived["N_re"])
        at = float(derived["A_t"])
        record: dict[str, object] = {
            "sobol_index": index,
            "unit_point": unit_points[index].tolist(),
            "parameters": params,
            "tensor_tilt": {
                "supplied_n_t": params["n_t"],
                "derived_n_t": float(derived["nt"]),
                "supplied_controls_cr_positive_spectrum": False,
            },
            "cutoff": {
                "r_cutoff": r_cut,
                "r_over_r_cutoff": params["r"] / r_cut,
                "initial_Delta_N": initial_delta,
                "initial_N_inf": None if derived["N_inf"] is None else float(derived["N_inf"]),
                "initial_N_re": float(derived["N_re"]),
                "stored_A_t": at,
                "A_t_class": "rounded_zero" if at == 0.0 else "normal" if at >= np.finfo(float).tiny else "subnormal",
            },
            "fast": {},
            "direct_reference": {},
        }
        if initial_delta is None or initial_delta <= 0.0:
            record["status"] = "invalid_cutoff_preflight"
            record["fast"] = {"called": False, "failure_reason": "initial_cutoff_invalid"}
            output["samples"].append(record)
            checkpoint()
            continue
        if at <= 0.0 or not math.isfinite(at):
            record["status"] = "invalid_tensor_amplitude_preflight"
            record["fast"] = {"called": False, "failure_reason": "stored_A_t_not_positive_finite"}
            output["samples"].append(record)
            checkpoint()
            continue

        start = time.perf_counter()
        caught_fast: list[warnings.WarningMessage]
        try:
            with warnings.catch_warnings(record=True) as caught_fast:
                warnings.simplefilter("always", RuntimeWarning)
                result = model.SGWB_iter(
                    engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                    threads=1, accuracy_mode="fast", freq_grid="goal",
                    physical_guard=True,
                )
            runtime = time.perf_counter() - start
            succeeded = result is model
            freqs = np.asarray(getattr(model, "f", []), dtype=np.float64)
            record["fast"] = {
                "called": True,
                "returned_model": succeeded,
                "failure_reason": model.fast_failure_reason,
                "failure_detail": getattr(model, "fast_failure_detail", None),
                "runtime_s": runtime,
                "root_DN_eff": float(model.cosmo_param["DN_eff"]) if succeeded else None,
                "mapped_DN_gw": float(model.DN_gw[-1]) if succeeded and len(model.DN_gw) else None,
                "root_residual_relative": float(model.DN_gw[-1] / model.cosmo_param["DN_eff"] - 1.0) if succeeded else None,
                "actual_derived_n_t": float(model.derived_param["nt"]) if succeeded else None,
                "frequency_count": int(freqs.size),
                "frequency_sha256_f64le": digest(freqs) if freqs.size else None,
                "public_spectrum_finite": bool(np.isfinite(model.log10OmegaGW).all()) if succeeded else None,
                "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught_fast],
            }
            record["status"] = "fast_admitted" if succeeded else "fast_rejected"
            if succeeded and index in REFERENCE_INDICES:
                root = float(model.cosmo_param["DN_eff"])
                fixed = LCDM_SG(**{**params, "DN_eff": root})
                ref_rows: dict[str, object] = {}
                for z_tail in REFERENCE_TAILS:
                    start_ref = time.perf_counter()
                    caught_ref: list[warnings.WarningMessage]
                    with warnings.catch_warnings(record=True) as caught_ref:
                        warnings.simplefilter("always")
                        ref = run_reference(
                            fixed, dn_eff=root, freq_subset=freqs,
                            z_tail=z_tail, rtol=REFERENCE_RTOL,
                            self_consistent=False, physical_guard=True,
                        )
                    ref_rows[str(z_tail)] = reference_summary(
                        ref, time.perf_counter() - start_ref, caught_ref
                    )
                    record["direct_reference"] = {
                        "source_sha256": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
                        "fixed_DN_eff": root,
                        "frequency_sha256_f64le": digest(freqs),
                        "z_tails": ref_rows,
                    }
                    checkpoint()
                record["direct_reference"]["z8_to_z10_map_relative_change"] = (
                    float(ref_rows["10.0"]["DN_gw"]) / float(ref_rows["8.0"]["DN_gw"]) - 1.0
                )
                record["direct_reference"]["z10_to_z12_map_relative_change"] = (
                    float(ref_rows["12.0"]["DN_gw"]) / float(ref_rows["10.0"]["DN_gw"]) - 1.0
                )
                record["direct_reference"]["z10_to_z12_max_abs_dex"] = float(np.nanmax(np.abs(
                    np.asarray([np.nan if v is None else float(v) for v in ref_rows["12.0"]["raw_log10_OmegaGW"]])
                    - np.asarray([np.nan if v is None else float(v) for v in ref_rows["10.0"]["raw_log10_OmegaGW"]])
                )))
                fast_logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
                root_err = root / float(ref_rows["12.0"]["DN_gw"]) - 1.0
                record["direct_reference"]["fast_root_relative_error_vs_z12_map"] = root_err
                ref_logs = np.asarray([np.nan if v is None else float(v) for v in ref_rows["12.0"]["raw_log10_OmegaGW"]])
                valid = np.isfinite(ref_logs) & np.isfinite(fast_logs)
                delta = fast_logs[valid] - ref_logs[valid]
                record["direct_reference"]["fast_z5_spectrum_vs_z12"] = {
                    "compared_nodes": int(np.count_nonzero(valid)),
                    "uncompared_nodes": int(freqs.size - np.count_nonzero(valid)),
                    "max_abs_dex": float(np.max(np.abs(delta))) if delta.size else None,
                    "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)) if delta.size else None,
                    "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3) if delta.size else False,
                }
                record["status"] = "sobol_fast_and_direct_reference_complete"
            output["samples"].append(record)
            checkpoint()
        except Exception as exc:
            record["status"] = "sample_exception"
            record["exception"] = {"type": type(exc).__name__, "message": str(exc)}
            record["elapsed_s"] = time.perf_counter() - start
            output["samples"].append(record)
            checkpoint()

    output["status"] = "complete_sobol_cutoff_preflight_with_predeclared_reference_cases"
    checkpoint()
    print(json.dumps({
        "output": str(OUTPUT),
        "samples": [{
            "sobol_index": row["sobol_index"],
            "T_re": row["parameters"]["T_re"],
            "kappa10": row["parameters"]["kappa10"],
            "status": row["status"],
            "fast_failure_reason": row.get("fast", {}).get("failure_reason"),
            "reference_completed": bool(row.get("direct_reference")),
        } for row in output["samples"]],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
