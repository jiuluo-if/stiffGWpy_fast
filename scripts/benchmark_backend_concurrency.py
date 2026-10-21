"""Measure full-path throughput with serialized versus thread-safe Numba calls.

中文：在 TBB/OpenMP 后端下，比较旧串行锁与 backend-aware 并发调用的双模型吞吐。
"""

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import psutil


def parse_affinity(value):
    return [int(item) for item in value.split(",") if item]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--repeats", type=int, default=30)
parser.add_argument("--threads", type=int, default=2)
parser.add_argument("--cases", default="default,highT")
parser.add_argument("--affinity", type=parse_affinity, default=parse_affinity("0,1,2,3"))
args = parser.parse_args()
if args.repeats < 1 or args.threads < 1 or len(args.affinity) < args.threads * 2:
    parser.error("repeats/threads must be positive and affinity must expose at least 2 * threads CPUs")

process = psutil.Process()
process.cpu_affinity(args.affinity)
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numba  # noqa: E402
import numpy as np  # noqa: E402

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    "default": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "highT": dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "stiff": dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    "high_kappa": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    "lowT": dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    "low_r": dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
}
case_names = args.cases.split(",")
unknown = sorted(set(case_names) - set(CASES))
if unknown:
    parser.error(f"unknown cases: {', '.join(unknown)}")
layer = numba.threading_layer()
if layer not in ("tbb", "omp"):
    raise RuntimeError(f"requires a thread-safe Numba layer, got {layer}")
config = FS.resolve_config("fast", threads=args.threads)
original_layer_check = FS._numba_layer_is_threadsafe


def run_one(name):
    model = LCDM_SG(**CASES[name])
    result = FS.SGWB_iter_fast(model, config=config, kink_split=True)
    if result is None or not getattr(model, "SGWB_converge", False):
        raise RuntimeError(f"solver failed: {name}, {getattr(model, 'fast_failure_reason', None)}")
    fields = ("f", "log10OmegaGW", "DN_gw", "g2", "w2")
    digest = {
        field: hashlib.sha256(
            np.ascontiguousarray(np.asarray(getattr(model, field), dtype=np.float64)).tobytes()
        ).hexdigest()
        for field in fields
    }
    return {"case": name, "digest": digest, "converged": True, "failure_reason": None}


def batch(executor, serialized):
    FS._numba_layer_is_threadsafe = (lambda: False) if serialized else original_layer_check
    try:
        start = time.perf_counter()
        futures = [executor.submit(run_one, name) for name in case_names]
        rows = [future.result() for future in futures]
        return time.perf_counter() - start, rows
    finally:
        FS._numba_layer_is_threadsafe = original_layer_check


def main():
    results = {"serialized": [], "concurrent": []}
    with ThreadPoolExecutor(max_workers=len(case_names)) as executor:
        warm = {label: batch(executor, serialized=(label == "serialized")) for label in ("serialized", "concurrent")}
        if warm["serialized"][1] != warm["concurrent"][1]:
            raise RuntimeError("warm output mismatch")
        reference_digests = {row["case"]: row["digest"] for row in warm["serialized"][1]}
        for index in range(args.repeats):
            order = ("serialized", "concurrent") if index % 2 == 0 else ("concurrent", "serialized")
            for label in order:
                duration, rows = batch(executor, serialized=(label == "serialized"))
                if rows != warm[label][1]:
                    raise RuntimeError(f"output mismatch in {label} round {index}")
                results[label].append(duration)
    ratios = [new / old - 1.0 for old, new in zip(results["serialized"], results["concurrent"])]
    rng = np.random.default_rng(20261021)
    values = np.asarray(ratios)
    boot = np.median(rng.choice(values, size=(20000, values.size), replace=True), axis=1)
    payload = {
        "run_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "candidate_code_commit": "9ad880624b9cd98f5caef3bcc11173df205725cc",
        "baseline_policy_commit": "8375477e805bc2ebdb7d171189fc56fb00cd0d26",
        "solver_source_sha256": hashlib.sha256(Path(FS.__file__).read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "baseline_policy": "previous wrapper predicate; explicit thread config always takes the solve lock",
        "candidate_policy": "serialize workqueue/unknown backend; permit explicit configs on TBB/OpenMP",
        "cases": case_names,
        "repeats": args.repeats,
        "threading_layer": layer,
        "threads_per_call": args.threads,
        "affinity": process.cpu_affinity(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "numba": numba.__version__,
        "numpy": np.__version__,
        "environment": {
            key: os.environ.get(key)
            for key in (
                "NUMBA_NUM_THREADS",
                "NUMBA_THREADING_LAYER",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
            )
        },
        "all_outputs_exact": True,
        "reference_digests": reference_digests,
        "serialized_wall_s": results["serialized"],
        "concurrent_wall_s": results["concurrent"],
        "serialized_median_s": statistics.median(results["serialized"]),
        "concurrent_median_s": statistics.median(results["concurrent"]),
        "paired_median_change_pct": float(np.median(values) * 100.0),
        "bootstrap_ci95_pct": [float(value * 100.0) for value in np.percentile(boot, [2.5, 97.5])],
        "bootstrap_replicates": 20000,
        "seed": 20261021,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: payload[key]
                for key in (
                    "threading_layer",
                    "paired_median_change_pct",
                    "bootstrap_ci95_pct",
                    "serialized_median_s",
                    "concurrent_median_s",
                    "all_outputs_exact",
                )
            },
            indent=2,
        )
    )
    print(args.output)


if __name__ == "__main__":
    main()
