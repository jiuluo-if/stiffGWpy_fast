"""Alternate modulo-based and next-index tensor assembly kernels.

中文：在固定线程和 CPU 亲和性下交替测量新旧 kernel，并核对完整输出摘要。
"""

import argparse
import ast
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
from pathlib import Path

import numba
import numpy as np
import psutil
import scipy

ROOT = Path(__file__).resolve().parents[1]
psutil.Process().cpu_affinity([0, 1])
sys.path.insert(0, str(ROOT / "scripts"))
import profile_fast_breakdown as profile  # noqa: E402

fast = profile.FS
BASELINE_COMMIT = "68b17acf8e049bf17ed100422da0f2f3553f15e9"
source = subprocess.check_output(
    ["git", "show", f"{BASELINE_COMMIT}:stiffgwpy_fast/fast_sgwb.py"],
    cwd=ROOT,
    text=True,
    encoding="utf-8",
)
tree = ast.parse(source)
function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "solve_kernel")
function.name = "_baseline_solve_kernel_py"
function.decorator_list = []
ast.fix_missing_locations(function)
exec(compile(ast.Module(body=[function], type_ignores=[]), fast.__file__, "exec"), fast.__dict__)
baseline_py = fast.__dict__.pop("_baseline_solve_kernel_py")
baseline_kernel = fast.njit(parallel=True, cache=False)(baseline_py)
candidate_kernel = fast.solve_kernel


def digest(model):
    fields = ("f", "log10OmegaGW", "DN_gw", "g2", "w2")
    return {
        name: hashlib.sha256(
            np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64)).tobytes()
        ).hexdigest()
        for name in fields
    }


def timed_variant(name, case, keep_arrays=False):
    fast.solve_kernel = baseline_kernel if name == "baseline" else candidate_kernel
    try:
        _, totals, model = profile.run_once(profile.CASES[case], kink_split=True)
    finally:
        fast.solve_kernel = candidate_kernel
    arrays = (
        {field: np.asarray(getattr(model, field)).copy() for field in ("f", "log10OmegaGW", "DN_gw", "g2", "w2")}
        if keep_arrays
        else None
    )
    return {
        "arrays": arrays,
        "total_s": float(totals["total"][0]),
        "kernel_s": float(sum(totals["tensor_solve_kernel"])),
        "digest": digest(model),
        "converged": bool(getattr(model, "SGWB_converge", False)),
        "failure_reason": getattr(model, "fast_failure_reason", None),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=80)
    parser.add_argument("--cases", default="default,highT,stiff,high_kappa,lowT,low_r")
    parser.add_argument("--col-step", type=int, default=8)
    args = parser.parse_args()
    out = args.output
    repeats = args.repeats
    cases = args.cases.split(",")
    col_step = args.col_step
    if repeats < 1:
        parser.error("--repeats must be positive")
    if not 1 <= col_step <= 8:
        parser.error("--col-step must be in [1, 8]")
    unknown_cases = sorted(set(cases) - set(profile.CASES))
    if unknown_cases:
        parser.error(f"unknown case(s): {', '.join(unknown_cases)}")
    fast.apply_accuracy_mode("fast")
    fast.set_col_step(col_step)
    fast.set_threads(2)
    result = {}
    for case in cases:
        warm = {name: timed_variant(name, case, keep_arrays=True) for name in ("baseline", "candidate")}
        if warm["baseline"]["digest"] != warm["candidate"]["digest"]:
            for field in warm["baseline"]["arrays"]:
                a, b = warm["baseline"]["arrays"][field], warm["candidate"]["arrays"][field]
                same = a.shape == b.shape and np.array_equal(a, b, equal_nan=True)
                if not same:
                    delta = np.abs(a - b) if a.shape == b.shape else np.array([])
                    print(
                        f"DIFF {case} {field} shapes={a.shape}/{b.shape} equal={same} max_abs={delta.max() if delta.size else 'shape'} count={np.count_nonzero(a != b) if a.shape == b.shape else 'shape'}"
                    )
            raise RuntimeError(f"warm output mismatch: {case}")
        for row in warm.values():
            row.pop("arrays")
        samples = {"baseline": [], "candidate": []}
        for r in range(repeats):
            for name in ("baseline", "candidate") if r % 2 == 0 else ("candidate", "baseline"):
                samples[name].append(timed_variant(name, case))
        exact = {name: {json.dumps(row["digest"], sort_keys=True) for row in rows} for name, rows in samples.items()}
        if exact["baseline"] != exact["candidate"] or len(exact["candidate"]) != 1:
            raise RuntimeError(f"digest mismatch: {case}")
        if any(not x["converged"] or x["failure_reason"] for rows in samples.values() for x in rows):
            raise RuntimeError(f"convergence mismatch: {case}")
        base = statistics.median(x["total_s"] for x in samples["baseline"])
        cand = statistics.median(x["total_s"] for x in samples["candidate"])
        paired = {}
        for metric in ("total_s", "kernel_s"):
            ratios = np.asarray(
                [
                    candidate[metric] / baseline[metric] - 1.0
                    for baseline, candidate in zip(samples["baseline"], samples["candidate"])
                ]
            )
            rng = np.random.default_rng(20261021)
            boot = np.median(rng.choice(ratios, size=(10000, ratios.size), replace=True), axis=1)
            paired[metric] = {
                "median_change_pct": float(np.median(ratios) * 100.0),
                "bootstrap_ci95_pct": [float(x * 100.0) for x in np.percentile(boot, [2.5, 97.5])],
            }
        result[case] = {
            "all_exact": True,
            "paired": paired,
            "baseline_total_median_ms": base * 1e3,
            "candidate_total_median_ms": cand * 1e3,
            "total_change_pct": (cand / base - 1.0) * 100.0,
            "baseline_kernel_median_ms": statistics.median(x["kernel_s"] for x in samples["baseline"]) * 1e3,
            "candidate_kernel_median_ms": statistics.median(x["kernel_s"] for x in samples["candidate"]) * 1e3,
            "samples": samples,
        }
    paired_summary = {}
    for metric in ("total_s", "kernel_s"):
        case_ratios = {
            case: np.asarray(
                [
                    candidate[metric] / baseline[metric] - 1.0
                    for baseline, candidate in zip(
                        result[case]["samples"]["baseline"],
                        result[case]["samples"]["candidate"],
                    )
                ]
            )
            for case in cases
        }
        per_case = {}
        for case, ratios in case_ratios.items():
            rng = np.random.default_rng(20261021)
            boot = np.median(rng.choice(ratios, size=(20000, ratios.size), replace=True), axis=1)
            per_case[case] = {
                "median_change_pct": float(np.median(ratios) * 100.0),
                "bootstrap_ci95_pct": [float(x * 100.0) for x in np.percentile(boot, [2.5, 97.5])],
            }
        round_geomeans = np.exp(np.mean(np.log1p(np.vstack(list(case_ratios.values()))), axis=0)) - 1.0
        rng = np.random.default_rng(20261021)
        boot = np.median(rng.choice(round_geomeans, size=(20000, round_geomeans.size), replace=True), axis=1)
        paired_summary[metric] = {
            "per_case": per_case,
            "round_geomean_median_change_pct": float(np.median(round_geomeans) * 100.0),
            "round_geomean_bootstrap_ci95_pct": [float(x * 100.0) for x in np.percentile(boot, [2.5, 97.5])],
            "bootstrap_replicates": 20000,
            "seed": 20261021,
        }
    tracked_env = (
        "FAST_THREADS",
        "NUMBA_NUM_THREADS",
        "NUMBA_THREADING_LAYER",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )
    payload = {
        "protocol": {
            "baseline_commit": BASELINE_COMMIT,
            "candidate": "monotonic next_assembly_k instead of inner-loop modulo",
            "repeats_per_variant_per_case": repeats,
            "cases": cases,
            "col_step": col_step,
            "threads": fast.get_num_threads(),
            "threading_layer": numba.threading_layer(),
            "parallel_chunksize": numba.get_parallel_chunksize(),
            "affinity": psutil.Process().cpu_affinity(),
        },
        "metadata": {
            "run_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "baseline_source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
            "candidate_source_sha256": hashlib.sha256(Path(fast.__file__).read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numba": numba.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "environment": {name: os.environ.get(name) for name in tracked_env},
        },
        "paired_summary": paired_summary,
        "results": result,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                k: {
                    m: round(v[m], 4)
                    for m in (
                        "baseline_total_median_ms",
                        "candidate_total_median_ms",
                        "total_change_pct",
                        "baseline_kernel_median_ms",
                        "candidate_kernel_median_ms",
                    )
                }
                for k, v in result.items()
            },
            indent=2,
        )
    )
    print(out)


if __name__ == "__main__":
    main()
