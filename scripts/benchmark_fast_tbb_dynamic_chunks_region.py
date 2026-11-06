"""Compare Fast's workqueue scheduling with TBB dynamic chunks on fixed cases.

This changes only Numba's parallel scheduling. The three worker processes run
the same current Fast source, thread count, model inputs, and solver settings.
"""

from __future__ import annotations

import hashlib
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
REGION_INPUT = ROOT / "docs" / "endpoint_tail_cr1_3d_region_20261105.json"
THREADS = int(os.environ.get("FAST_TBB_BENCH_THREADS", "2"))
PAIRS = int(os.environ.get("FAST_TBB_BENCH_PAIRS", "25"))
THREAD_TAG = "" if THREADS == 2 else f"_{THREADS}threads"
SMOKE = os.environ.get("FAST_TBB_BENCH_SMOKE") == "1"
CASE_FILTER = {
    name for name in os.environ.get("FAST_TBB_BENCH_CASES", "").split(",") if name
}
FILTER_TAG = "_focused" if CASE_FILTER else ""
OUTPUT = Path(os.environ.get(
    "FAST_TBB_BENCH_OUTPUT",
    ROOT / "docs" / f"fast_tbb_dynamic_config_region_{PAIRS}pairs{THREAD_TAG}{FILTER_TAG}_20261106.json",
))
VARIANTS = {
    "workqueue_static": ("workqueue", 0),
    "tbb_static": ("tbb", 0),
    "tbb_dynamic_1": ("tbb", 1),
}


def _worker_main(layer: str, default_chunk: int) -> int:
    from scripts._resource_budget import apply_environment, telemetry

    apply_environment(numba_threads=THREADS, fast_threads=THREADS)
    import numpy as np
    from numba import get_num_threads, njit, prange, threading_layer
    try:
        import psutil
    except ImportError:  # pragma: no cover - optional instrumentation only
        psutil = None

    @njit(parallel=True)
    def initialize_parallel_backend():
        values = np.zeros(2, dtype=np.int64)
        for i in prange(values.size):
            values[i] = i
        return values

    initialize_parallel_backend()

    from stiffgwpy_fast import fast_sgwb as FS
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    FS.apply_accuracy_mode("fast")
    FS.set_threads(THREADS)
    if threading_layer() != layer or get_num_threads() != THREADS:
        raise RuntimeError(
            f"worker requested layer={layer},threads={THREADS}; "
            f"got layer={threading_layer()},threads={get_num_threads()}"
        )

    process = None if psutil is None else psutil.Process()
    peak_rss_bytes = 0 if process is None else int(process.memory_info().rss)

    def solve(params: dict, chunk: int) -> dict:
        nonlocal peak_rss_bytes
        started = time.perf_counter()
        model = LCDM_SG(**params)
        configured_chunk = None if layer == "workqueue" else int(chunk)
        result = FS.SGWB_iter_fast(
            model,
            tol=1e-6,
            freq_res=1.0,
            sigma_exact=False,
            transition_refine=False,
            kink_split=True,
            freq_grid="goal",
            config=FS.resolve_config(
                "fast", threads=THREADS, parallel_chunksize=configured_chunk
            ),
            frequency_quadrature="pchip",
            physical_guard=False,
        )
        elapsed = time.perf_counter() - started
        if result is not model or not model.SGWB_converge:
            raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
        freqs = np.asarray(model.f, dtype=np.float64)
        spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
        digest = hashlib.sha256(
            freqs.astype("<f8").tobytes()
            + spectrum.astype("<f8").tobytes()
            + np.asarray([model.cosmo_param["DN_eff"]], dtype="<f8").tobytes()
        ).hexdigest()
        rss_bytes = None if process is None else int(process.memory_info().rss)
        if rss_bytes is not None:
            peak_rss_bytes = max(peak_rss_bytes, rss_bytes)
        return {
            "runtime_s": float(elapsed),
            "root": float(model.cosmo_param["DN_eff"]),
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "spectrum_sha256": hashlib.sha256(spectrum.astype("<f8").tobytes()).hexdigest(),
            "result_sha256": digest,
            "warnings": list(getattr(model, "fast_warnings", [])),
            "rss_after_mb": None if rss_bytes is None else rss_bytes / (1024.0 * 1024.0),
            "rss_peak_mb": None if process is None else peak_rss_bytes / (1024.0 * 1024.0),
        }

    print(json.dumps({
        "ready": True,
        "resource": telemetry(process=process, workers=1, threads=THREADS),
        "rss_before_mb": None if process is None else peak_rss_bytes / (1024.0 * 1024.0),
    }), flush=True)
    for line in sys.stdin:
        try:
            request = json.loads(line)
            if request.get("stop"):
                break
            answer = solve(request["parameters"], int(request["chunk_size"]))
            print(json.dumps({"ok": True, **answer}, allow_nan=False), flush=True)
        except Exception as exc:  # send worker failures to the orchestrator
            print(json.dumps({"ok": False, "error": f"{type(exc).__name__}: {exc}"}), flush=True)
    return 0


def _cases() -> dict[str, dict]:
    normal = {
        "default": {"r": 1e-2, "cr": 1, "T_re": 2e3, "kappa10": 1e-2},
        "lowT": {"r": 1e-2, "cr": 1, "T_re": 1e1, "kappa10": 1e-2},
        "highT": {"r": 1e-2, "cr": 1, "T_re": 1e4, "kappa10": 1e-2},
        "stiff": {"r": 1e-1, "cr": 1, "T_re": 2e3, "kappa10": 1e-2},
        "high_kappa": {"r": 1e-2, "cr": 1, "T_re": 2e3, "kappa10": 1.0},
    }
    region = json.loads(REGION_INPUT.read_text(encoding="utf-8"))
    region_cases = {
        name: dict(case["parameters"])
        for name, case in region["cases"].items()
    }
    cases = {**{f"trusted/{k}": v for k, v in normal.items()},
             **{f"cr1_region/{k}": v for k, v in region_cases.items()}}
    if SMOKE:
        return {
            "trusted/lowT": cases["trusted/lowT"],
            "cr1_region/corner_r1_T0_k0": cases["cr1_region/corner_r1_T0_k0"],
        }
    if CASE_FILTER:
        missing = CASE_FILTER - set(cases)
        if missing:
            raise ValueError(f"unknown benchmark case filter(s): {sorted(missing)}")
        return {name: cases[name] for name in cases if name in CASE_FILTER}
    return cases


def _spawn_worker(layer: str, chunk: int):
    environment = os.environ.copy()
    environment.update({
        "NUMBA_THREADING_LAYER": layer,
        "NUMBA_NUM_THREADS": str(THREADS),
        "FAST_THREADS": str(THREADS),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--worker", layer, str(chunk)],
        cwd=ROOT,
        env=environment,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        bufsize=1,
    )
    ready = json.loads(process.stdout.readline())
    if not ready.get("ready"):
        raise RuntimeError(f"worker failed to initialize: {ready}")
    return process, ready["resource"]


def _request(worker, params: dict, chunk: int) -> dict:
    process = worker
    process.stdin.write(json.dumps({"parameters": params, "chunk_size": chunk}) + "\n")
    process.stdin.flush()
    result = json.loads(process.stdout.readline())
    if not result.get("ok"):
        raise RuntimeError(result.get("error", "unknown Fast worker error"))
    return result


def _run_parent() -> dict:
    import numpy as np

    cases = _cases()
    workers = {}
    resources = {}
    for variant, (layer, chunk) in VARIANTS.items():
        workers[variant], resources[variant] = _spawn_worker(layer, chunk)

    try:
        warmups = {}
        for case_name, params in cases.items():
            warmups[case_name] = {}
            for variant, (_, chunk) in VARIANTS.items():
                result = _request(workers[variant], params, chunk)
                warmups[case_name][variant] = result

        measurements = {
            case_name: {variant: [] for variant in VARIANTS}
            for case_name in cases
        }
        for case_index, (case_name, params) in enumerate(cases.items()):
            variants = list(VARIANTS)
            for pair in range(PAIRS):
                offset = (pair + case_index) % len(variants)
                order = variants[offset:] + variants[:offset]
                if pair % 2:
                    order.reverse()
                for variant in order:
                    chunk = VARIANTS[variant][1]
                    result = _request(workers[variant], params, chunk)
                    measurements[case_name][variant].append(result)
    finally:
        for process in workers.values():
            try:
                process.stdin.write('{"stop":true}\n')
                process.stdin.flush()
                process.wait(timeout=10)
            except Exception:
                process.kill()

    rows = {}
    pairwise = {
        "tbb_static_over_workqueue_static": {"p50_over_1_05": [], "p95_over_1_05": []},
        "tbb_dynamic_1_over_tbb_static": {"p50_over_1_05": [], "p95_over_1_05": []},
        "tbb_dynamic_1_over_workqueue_static": {"p50_over_1_05": [], "p95_over_1_05": []},
    }
    for case_name, variants in measurements.items():
        stats = {}
        for variant, samples in variants.items():
            runtimes = [sample["runtime_s"] for sample in samples]
            stats[variant] = {
                "p50_s": float(statistics.median(runtimes)),
                "p95_s": float(np.percentile(runtimes, 95)),
            }
        reference = warmups[case_name]["workqueue_static"]["result_sha256"]
        equality = {
            variant: all(sample["result_sha256"] == reference for sample in samples)
            and warmups[case_name][variant]["result_sha256"] == reference
            for variant, samples in variants.items()
        }
        memory = {
            variant: {
                "warm_rss_after_mb": warmups[case_name][variant]["rss_after_mb"],
                "peak_rss_mb": max(
                    [sample["rss_peak_mb"] for sample in samples]
                    + [warmups[case_name][variant]["rss_peak_mb"]]
                ) if warmups[case_name][variant]["rss_peak_mb"] is not None else None,
            }
            for variant, samples in variants.items()
        }
        ratios = {
            "tbb_static_over_workqueue_static_p50": stats["tbb_static"]["p50_s"] / stats["workqueue_static"]["p50_s"],
            "tbb_static_over_workqueue_static_p95": stats["tbb_static"]["p95_s"] / stats["workqueue_static"]["p95_s"],
            "tbb_dynamic_1_over_tbb_static_p50": stats["tbb_dynamic_1"]["p50_s"] / stats["tbb_static"]["p50_s"],
            "tbb_dynamic_1_over_tbb_static_p95": stats["tbb_dynamic_1"]["p95_s"] / stats["tbb_static"]["p95_s"],
            "tbb_dynamic_1_over_workqueue_static_p50": stats["tbb_dynamic_1"]["p50_s"] / stats["workqueue_static"]["p50_s"],
            "tbb_dynamic_1_over_workqueue_static_p95": stats["tbb_dynamic_1"]["p95_s"] / stats["workqueue_static"]["p95_s"],
        }
        rows[case_name] = {
            "parameters": cases[case_name],
            "warmups": warmups[case_name],
            "timing": stats,
            "ratios": ratios,
            "memory": memory,
            "bitwise_result_match": equality,
            "candidate_same_as_workqueue": all(equality.values()),
        }
        for name, prefix in (
            ("tbb_static_over_workqueue_static", "tbb_static_over_workqueue_static"),
            ("tbb_dynamic_1_over_tbb_static", "tbb_dynamic_1_over_tbb_static"),
            ("tbb_dynamic_1_over_workqueue_static", "tbb_dynamic_1_over_workqueue_static"),
        ):
            for q in ("p50", "p95"):
                if ratios[f"{prefix}_{q}"] > 1.05:
                    pairwise[name][f"{q}_over_1_05"].append(case_name)

    return {
        "schema_version": 1,
        "status": "complete",
        "experiment": "Fast full-solve comparison of workqueue static versus TBB scheduling via FastSolverConfig",
        "decision": "PENDING_REVIEW",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "case_count": len(cases),
        "pairs_per_case_per_variant": PAIRS,
        "threads": THREADS,
        "variants": VARIANTS,
        "resources": resources,
        "source_hashes": {
            "fast_sgwb": hashlib.sha256((ROOT / "stiffgwpy_fast" / "fast_sgwb.py").read_bytes()).hexdigest(),
            "config": hashlib.sha256((ROOT / "stiffgwpy_fast" / "config.py").read_bytes()).hexdigest(),
            "exact_background": hashlib.sha256((ROOT / "stiffgwpy_fast" / "exact_background.py").read_bytes()).hexdigest(),
            "freq_adaptive": hashlib.sha256((ROOT / "stiffgwpy_fast" / "freq_adaptive.py").read_bytes()).hexdigest(),
            "region_input": hashlib.sha256(REGION_INPUT.read_bytes()).hexdigest(),
            "benchmark_script": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "cases": rows,
        "gate_summary": pairwise,
        "limitations": [
            "TBB must be selected before importing Numba; the current project resource helper defaults to workqueue.",
            "This benchmark changes scheduling only. All solver outputs are required to match bitwise.",
            "The case set is five standard profiles plus one fixed sampled cr=1 region panel, not a mathematical proof for every valid parameter.",
            "Only two Numba threads and this Windows host are covered.",
        ],
    }


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        return _worker_main(sys.argv[2], int(sys.argv[3]))

    record = _run_parent()
    record["decision"] = "REJECTED" if any(
        not row["candidate_same_as_workqueue"]
        or row["ratios"]["tbb_dynamic_1_over_workqueue_static_p50"] > 1.05
        or row["ratios"]["tbb_dynamic_1_over_workqueue_static_p95"] > 1.05
        for row in record["cases"].values()
    ) else "SCREEN_PASS"
    output = ROOT / "docs" / f"fast_tbb_dynamic_config_smoke{THREAD_TAG}_20261106.json" if SMOKE else OUTPUT
    output.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    counts = {}
    for region_name in ("trusted", "cr1_region"):
        group = [row for name, row in record["cases"].items() if name.startswith(region_name + "/")]
        counts[region_name] = {
            "case_count": len(group),
            "tbb_chunk1_p50_over_workqueue_over_1_05": sum(
                row["ratios"]["tbb_dynamic_1_over_workqueue_static_p50"] > 1.05 for row in group
            ),
            "tbb_chunk1_p95_over_workqueue_over_1_05": sum(
                row["ratios"]["tbb_dynamic_1_over_workqueue_static_p95"] > 1.05 for row in group
            ),
            "all_outputs_bitwise_equal": all(row["candidate_same_as_workqueue"] for row in group),
        }
    print(json.dumps({"decision": record["decision"], "summary": counts, "artifact": str(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
