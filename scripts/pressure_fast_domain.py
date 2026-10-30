"""Deterministic Fast parameter-domain pressure screen.

中文：生成逐参数边界、耦合边缘与 Sobol 样本；逐点记录预检、求解状态、失败类别、精度有效性和耗时。
The resource preflight is an experiment limit, never an input transformation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import psutil
import scipy
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import _has_exact_zero_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-2,
    "n_t": 0.5,
    "cr": 0.0,
    "T_re": 2e3,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
PARAMETERS = tuple(BASE)
LOG_VALUES = (0.0, 1e-300, 1e-200, 1e-100, 1e-40, 1e-20, 1e-8, 1e-2,
              1.0, 1e2, 1e20, 1e100, 1e200, 1e300)
SIGNED_VALUES = (-1e100, -1e20, -1e6, -100.0, -1.0, -1e-12, 0.0,
                 1e-12, 1.0, 100.0, 1e6, 1e20, 1e100)


def _case(label: str, **updates: float) -> dict:
    params = BASE.copy()
    params.update(updates)
    return {"label": label, "params": params}


def boundary_cases() -> list[dict]:
    """Cover zero, near-zero, logarithmic, sign and thermal-transition edges."""
    cases = []
    for cr in (0.0, 1.0):
        cases.append(_case(f"baseline_cr_{int(cr)}", cr=cr))
    for key in PARAMETERS:
        if key == "cr":
            values = (-1.0, 0.0, 1.0)
        elif key in ("r", "A_s", "kappa10"):
            values = LOG_VALUES
        elif key == "T_re":
            transitions = (float(gp.T_max), float(gp.T_i), 1e-3, 1e-2, 1e-1,
                           1.0, 1e2, 1e4, 1e6, 1e10, 1e20, 1e100)
            values = (0.0, 1e-300, *transitions)
        elif key in ("n_t", "H0", "Omega_bh2", "Omega_ch2", "DN_eff"):
            values = SIGNED_VALUES
        elif key == "DN_re":
            values = (-100.0, -1.0, -1e-12, 0.0, 1e-12, 1e-6, 1.0,
                      10.0, 50.0, 100.0, 300.0, 1e6)
        else:
            values = SIGNED_VALUES
        for cr in (0.0, 1.0):
            for index, value in enumerate(values):
                update = {key: float(value), "cr": cr}
                cases.append(_case(f"boundary_{key}_cr{int(cr)}_{index}", **update))

    for cr in (0.0, 1.0):
        for nt in (-2.0, 0.0, 1.0, 2.0, 5.0, 10.0):
            for kappa in (1e-8, 1e-2, 1.0, 1e2, 1e4):
                for tre in (1e-3, float(gp.T_max), 1e4, 1e8):
                    cases.append(_case(
                        f"coupled_cr{int(cr)}_nt{nt:g}_k{kappa:g}_T{tre:g}",
                        cr=cr, n_t=nt, kappa10=kappa, T_re=tre,
                    ))
    return cases


def coupled_cases() -> list[dict]:
    """Return the deterministic ``cr × n_t × kappa10 × T_re`` cross."""
    return [case for case in boundary_cases()
            if case["label"].startswith("coupled_")]


def qmc_cases(points: np.ndarray, prefix: str) -> list[dict]:
    """Map unit-cube samples into broad signed/log parameter ranges."""
    rows = []
    for index, u in enumerate(points):
        params = BASE.copy()
        params.update({
            "r": 10.0 ** (-200.0 + 202.0 * u[0]),
            "n_t": -10.0 + 20.0 * u[1],
            "kappa10": 10.0 ** (-100.0 + 200.0 * u[2]),
            "T_re": 10.0 ** (-50.0 + 100.0 * u[3]),
            "DN_re": -50.0 + 350.0 * u[4],
            "Omega_bh2": -10.0 + 20.0 * u[5],
            "Omega_ch2": -10.0 + 20.0 * u[6],
            "H0": -100.0 + 300.0 * u[7],
            "A_s": 10.0 ** (-100.0 + 102.0 * u[8]),
            "DN_eff": -1e6 + 2e6 * u[9],
            "cr": float(index % 2),
        })
        rows.append({"label": f"{prefix}_{index:04d}", "params": params})
    return rows


def sobol_cases(power: int, seed: int) -> list[dict]:
    """Sobol powers of two preserve sequence balance properties."""
    points = qmc.Sobol(d=10, scramble=True, seed=seed).random_base2(power)
    return qmc_cases(points, f"sobol_{seed}")


def latin_cases(power: int, seed: int) -> list[dict]:
    points = qmc.LatinHypercube(d=10, scramble=True, seed=seed).random(2**power)
    return qmc_cases(points, f"lhs_{seed}")


def _classify_reason(reason: str | None) -> str:
    if reason in ("invalid_r", "invalid_cutoff", "shared_Neff_guard"):
        return "explicit_input_or_physics_guard"
    if reason == "invalid_iterate_cutoff":
        return "outer_iterate_derived_cutoff_invalid"
    if reason == "nonpositive_h2":
        return "background_h2_nonpositive"
    if reason == "max_iter":
        return "outer_iteration_limit"
    if reason and reason.startswith("nonfinite"):
        return "nonfinite_numerical_state"
    if reason == "invalid_frequency_grid":
        return "frequency_grid_failure"
    if reason == "tensor_amplitude_underflow":
        return "floating_point_underflow"
    if reason:
        return "solver_failure"
    return "unknown_failure"


def _json_scalar(value):
    if value is None:
        return None
    scalar = float(value)
    return scalar if math.isfinite(scalar) else str(scalar)


def _finish_record(result: dict, started: float) -> dict:
    result["runtime_s"] = time.perf_counter() - started
    result["rss_mb_after"] = psutil.Process().memory_info().rss / 1024**2
    return result


def evaluate(item: dict, max_n_inf: float, threads: int) -> dict:
    params = item["params"]
    digest = hashlib.sha256(
        json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    started = time.perf_counter()
    result = {
        "label": item["label"], "params": params, "parameter_sha256": digest,
        "status": None, "classification": None, "reason": None,
        "derived": {}, "runtime_s": None, "n_freq": None,
        "DN_gw": None, "spectrum_finite_bins": None,
        "fast_accuracy_status": None, "fast_accuracy_budget": None,
    }
    try:
        model = LCDM_SG(**params)
        derived = model.derived_param
        result["derived"] = {
            key: _json_scalar(value)
            for key, value in derived.items()
            if key in ("N_inf", "nt", "A_t", "kappa_s", "Omega_sh2")
            and (value is None or np.isscalar(value))
        }
    except Exception as exc:
        result.update(status="exception", classification="parameter_derivation_exception",
                      reason=repr(exc))
        return _finish_record(result, started)

    n_inf = derived.get("N_inf")
    if n_inf is None or not math.isfinite(float(n_inf)) or float(n_inf) <= 0.0:
        result.update(status="preflight_reject", classification="invalid_derived_cutoff",
                      reason="N_inf is undefined, non-finite, or non-positive")
        return _finish_record(result, started)
    if float(n_inf) > max_n_inf:
        result.update(status="preflight_skip", classification="resource_limit_preflight",
                      reason=f"N_inf={float(n_inf):.8g} exceeds configured {max_n_inf:g}")
        return _finish_record(result, started)

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-6, z_tail=5.0,
                freq_res=1.0, threads=threads, accuracy_mode="fast",
                physical_guard=False,
            )
        result["warnings"] = [str(w.message) for w in caught]
        result["reason"] = model.fast_failure_reason
        result["background_failure_index"] = getattr(
            model, "_fast_background_failure_index", None
        )
        result["background_failure_detail"] = getattr(
            model, "_fast_background_failure_detail", None
        )
        if out is None:
            result.update(status="failed", classification=_classify_reason(model.fast_failure_reason))
        else:
            log_spectrum = np.asarray(model.log10OmegaGW, dtype=float)
            physical_spectrum = np.asarray(model.Ogw_today) - np.asarray(model.Oj_today)
            dn_gw = float(model.DN_gw[-1])
            exact_zero_tensor = (
                _has_exact_zero_tensor_amplitude(params)
                and dn_gw == 0.0 and np.all(physical_spectrum == 0.0)
                and np.isneginf(log_spectrum).all()
            )
            if exact_zero_tensor:
                output_status = "success"
                output_classification = "exact_zero_tensor_spectrum"
            elif (np.isfinite(log_spectrum).all()
                  and np.all(physical_spectrum > 0.0)):
                output_status = "success"
                output_classification = "finite_positive_spectrum"
            else:
                output_status = "returned_nonfinite_output"
                output_classification = "degenerate_or_invalid_spectrum"
            result.update(
                status=output_status,
                classification=output_classification,
                n_freq=int(log_spectrum.size),
                spectrum_finite_bins=int(np.count_nonzero(np.isfinite(log_spectrum))),
                DN_gw=dn_gw,
                f_range=[float(np.min(model.f)), float(np.max(model.f))],
                dn_gw_background_nodes=int(len(model.DN_gw)),
                outer_reuse=bool(getattr(model, "outer_full_reuse_used", False)),
            )
            accuracy_budget = getattr(model, "local_error_budget", None)
            if exact_zero_tensor:
                result["fast_accuracy_status"] = "analytic_exact_limit"
            elif accuracy_budget is not None:
                result.update(
                    fast_accuracy_status=accuracy_budget.get(
                        "certification_status", "uncertified"),
                    fast_accuracy_budget=accuracy_budget,
                )
    except Exception as exc:
        result.update(status="exception", classification="solver_exception", reason=repr(exc))
    return _finish_record(result, started)


def _coordinate_midpoint(left: float, right: float, key: str) -> float:
    if key in ("r", "A_s", "kappa10", "T_re") and left > 0.0 and right > 0.0:
        return math.exp(0.5 * (math.log(left) + math.log(right)))
    return 0.5 * (left + right)


def minimize_failure(item: dict, original: dict, max_n_inf: float,
                     threads: int, steps: int = 6) -> dict:
    """Move one coordinate at a time toward the baseline, retaining the same failure class."""
    failing = item["params"].copy()
    passing = BASE.copy()
    passing["cr"] = failing["cr"]
    trace = []
    evaluations = 0
    for key in PARAMETERS:
        if key == "cr" or failing[key] == passing[key]:
            continue
        fail_value = failing[key]
        pass_value = passing[key]
        for _ in range(steps):
            midpoint = _coordinate_midpoint(fail_value, pass_value, key)
            candidate = failing.copy()
            candidate[key] = midpoint
            trial = evaluate({"label": f"min_{item['label']}_{key}", "params": candidate},
                             max_n_inf, threads)
            evaluations += 1
            same_failure = (trial["status"] != "success"
                            and trial["classification"] == original["classification"])
            if same_failure:
                fail_value = midpoint
                failing = candidate
            elif trial["status"] == "success":
                pass_value = midpoint
            else:
                trace.append({"parameter": key, "probe": midpoint,
                              "classification": trial["classification"],
                              "status": trial["status"], "boundary_resolved": False})
                break
            trace.append({"parameter": key, "probe": midpoint,
                          "classification": trial["classification"],
                          "status": trial["status"]})
        passing[key] = pass_value
    return {"params": failing, "coordinate_evaluations": evaluations, "trace": trace}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("smoke", "coupled", "boundaries", "sobol", "latin", "all"), default="smoke")
    parser.add_argument("--power", type=int, default=5, help="Sobol sample count is 2**power")
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--max-n-inf", type=float, default=120.0)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--start-index", type=int, default=0,
                        help="first index in the selected deterministic suite")
    parser.add_argument("--max-cases", type=int, default=None,
                        help="optional bounded chunk size for resumable scans")
    parser.add_argument("--minimize-failures", action="store_true",
                        help="bisect each failed coordinate toward the default model")
    parser.add_argument("--output", type=Path, default=Path("fast_domain_pressure.jsonl"))
    args = parser.parse_args()
    if args.power < 0 or args.power > 16:
        parser.error("--power must be in [0, 16]")
    if args.max_n_inf <= 0.0 or args.threads < 1:
        parser.error("--max-n-inf and --threads must be positive")

    if args.suite == "smoke":
        cases = [_case("baseline_cr0", cr=0.0), _case("baseline_cr1", cr=1.0),
                 _case("zero_stiff", cr=0.0, kappa10=0.0),
                 _case("near_zero_r", cr=0.0, r=1e-200)]
    elif args.suite == "coupled":
        cases = coupled_cases()
    elif args.suite == "boundaries":
        cases = boundary_cases()
    elif args.suite == "sobol":
        cases = sobol_cases(args.power, args.seed)
    elif args.suite == "latin":
        cases = latin_cases(args.power, args.seed)
    else:
        cases = boundary_cases() + sobol_cases(args.power, args.seed) + latin_cases(args.power, args.seed)
    if args.start_index < 0 or (args.max_cases is not None and args.max_cases < 1):
        parser.error("--start-index must be non-negative and --max-cases must be positive")
    suite_size = len(cases)
    first_index = args.start_index
    cases = cases[first_index:first_index + args.max_cases if args.max_cases is not None else None]

    FS.apply_accuracy_mode("fast")
    process = psutil.Process()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "record_type": "manifest", "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "suite": args.suite, "seed": args.seed, "case_count": len(cases),
        "suite_total_count": suite_size, "start_index": first_index,
        "sobol_power": args.power, "max_n_inf_preflight": args.max_n_inf,
        "threads_argument": args.threads, "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "logical_cpu_count": os.cpu_count(),
        "process_cpu_affinity": (process.cpu_affinity()
                                 if hasattr(process, "cpu_affinity") else None),
        "numba_threads_ambient": int(FS.get_num_threads()),
        "blas_threads": {key: os.environ.get(key) for key in
                         ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                          "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")},
        "acceptance": {"dn_gw_relative": 2e-4, "spectrum_dex": 1e-3,
                       "requires_separate_reference": True},
        "preflight_note": "resource_limit_preflight is not a solver failure and inputs are never clipped",
    }
    counts: Counter = Counter()
    accuracy_counts: Counter = Counter()
    runtimes = []
    dn_error_estimates = []
    batch_started = time.perf_counter()
    max_rss_mb = process.memory_info().rss / 1024**2
    with output.open("w", encoding="utf-8") as stream:
        stream.write(json.dumps(metadata, ensure_ascii=False) + "\n")
        for local_index, item in enumerate(cases):
            index = first_index + local_index
            record = evaluate(item, args.max_n_inf, args.threads)
            record["sample_index"] = index
            if args.minimize_failures and record["status"] != "success":
                record["minimized_counterexample"] = minimize_failure(
                    item, record, args.max_n_inf, args.threads)
            stream.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
            stream.flush()
            counts[record["status"]] += 1
            if record.get("fast_accuracy_status") is not None:
                accuracy_counts[record["fast_accuracy_status"]] += 1
            accuracy_budget = record.get("fast_accuracy_budget") or {}
            if accuracy_budget.get("DN_gw_error") is not None:
                dn_error_estimates.append(float(accuracy_budget["DN_gw_error"]))
            if record["runtime_s"] is not None:
                runtimes.append(record["runtime_s"])
            max_rss_mb = max(max_rss_mb, process.memory_info().rss / 1024**2)
            summary_record = {key: record.get(key) for key in
                              ("sample_index", "label", "status", "classification",
                               "reason", "fast_accuracy_status", "runtime_s")}
            summary_record["DN_gw_error_estimate"] = accuracy_budget.get("DN_gw_error")
            print(json.dumps(summary_record, ensure_ascii=False), flush=True)
    batch_wall_s = time.perf_counter() - batch_started
    summary = {
        "commit": metadata["commit"], "suite": args.suite, "case_count": len(cases),
        "suite_total_count": suite_size, "start_index": first_index,
        "status_counts": dict(counts),
        "fast_accuracy_status_counts": dict(accuracy_counts),
        "fast_dn_gw_error_estimate": {
            "min": float(min(dn_error_estimates)) if dn_error_estimates else None,
            "median": float(np.percentile(dn_error_estimates, 50)) if dn_error_estimates else None,
            "p95": float(np.percentile(dn_error_estimates, 95)) if dn_error_estimates else None,
            "max": float(max(dn_error_estimates)) if dn_error_estimates else None,
        },
        "runtime_p50_s": float(np.percentile(runtimes, 50)) if runtimes else None,
        "runtime_p95_s": float(np.percentile(runtimes, 95)) if runtimes else None,
        "sum_solver_runtime_s": sum(runtimes),
        "batch_wall_s": batch_wall_s,
        "batch_throughput_cases_per_s": len(cases) / batch_wall_s if batch_wall_s else None,
        "process_rss_mb_peak_observed": max_rss_mb,
        "process_rss_mb_end": process.memory_info().rss / 1024**2,
        "jsonl": str(output),
    }
    summary_path = output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
