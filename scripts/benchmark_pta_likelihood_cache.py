"""Compare cached PTA likelihoods with their pre-cache implementation.

中文：对比当前 PTA 似然缓存与缓存改动前实现，并记录逐次计时和初始化成本。
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import pickle
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import cobaya
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BASELINE_COMMIT = "1f648a50af9bb57ca929e967d85a9b888c92aba5"
REPEATS = 25
THEORY_FREQUENCIES = np.linspace(-13.0, -5.0, 241)
THEORY_SPECTRUM = -9.0 + 0.15 * np.sin(np.linspace(0.0, 4.0, THEORY_FREQUENCIES.size))
HUBBLE = 2.18e-18
PARAMETERS = {"A_BBH": -15.6, "gamma_BBH": 4.7}
DATA_DIR = ROOT / "stiffgwpy_fast" / "cobaya" / "likelihoods" / "PTA"


def _load_baseline_class(source_path: str, class_name: str):
    source = subprocess.check_output(
        ["git", "show", f"{BASELINE_COMMIT}:{source_path}"],
        cwd=ROOT,
        text=True,
    )
    namespace = {"__name__": f"_baseline_{class_name}", "__file__": source_path}
    exec(compile(source, source_path, "exec"), namespace)
    return namespace[class_name]


def _measure(adapter: str, module_path: str, source_path: str, class_name: str, attrs: dict) -> dict:
    baseline_class = _load_baseline_class(source_path, class_name)
    cached_class = getattr(importlib.import_module(module_path), class_name)
    baseline = SimpleNamespace(**attrs)
    cached = SimpleNamespace(**attrs)

    baseline_start = time.perf_counter()
    baseline_class.initialize(baseline)
    baseline_initialize_ms = 1e3 * (time.perf_counter() - baseline_start)

    cached_start = time.perf_counter()
    cached_class.initialize(cached)
    cached_initialize_ms = 1e3 * (time.perf_counter() - cached_start)

    call_args = (THEORY_FREQUENCIES, THEORY_SPECTRUM, HUBBLE)

    def baseline_call():
        return baseline_class.log_likelihood(baseline, *call_args, **PARAMETERS)

    def cached_call():
        return cached_class.log_likelihood(cached, *call_args, **PARAMETERS)

    baseline_call()
    cached_call()

    baseline_times = []
    cached_times = []
    baseline_values = []
    cached_values = []
    for repeat in range(REPEATS):
        order = ("baseline", "cached") if repeat % 2 == 0 else ("cached", "baseline")
        for mode in order:
            start = time.perf_counter()
            if mode == "baseline":
                value = float(baseline_call())
                baseline_times.append(time.perf_counter() - start)
                baseline_values.append(value)
            else:
                value = float(cached_call())
                cached_times.append(time.perf_counter() - start)
                cached_values.append(value)

    differential_cases = []
    for index, (amplitude, slope) in enumerate(
        [
            (-18.0, 0.0),
            (-18.0, 7.0),
            (-17.2, 3.3),
            (-16.5, 4.0),
            (-15.6, 4.7),
            (-15.0, 5.0),
            (-14.0, 6.0),
            (-12.5, 1.0),
            (-11.5, 6.5),
            (-11.0, 13.0 / 3.0),
        ]
    ):
        spectrum = THEORY_SPECTRUM + 0.03 * index
        parameters = {"A_BBH": amplitude, "gamma_BBH": slope}
        inputs = (THEORY_FREQUENCIES, spectrum, HUBBLE + index * 1e-20)
        baseline_value = float(baseline_class.log_likelihood(baseline, *inputs, **parameters))
        cached_value = float(cached_class.log_likelihood(cached, *inputs, **parameters))
        differential_cases.append(
            {
                "index": index,
                "parameters": parameters,
                "baseline": baseline_value,
                "cached": cached_value,
                "exact": baseline_value == cached_value,
            }
        )

    caches = getattr(cached, "_kde_by_frequency", getattr(cached, "_density_splines", ()))
    if len(caches) != attrs["Nfreqs"]:
        raise RuntimeError(f"Expected {attrs['Nfreqs']} cached frequency objects for {adapter}")
    baseline_median = statistics.median(baseline_times)
    cached_median = statistics.median(cached_times)
    return {
        "adapter": adapter,
        "repeats_per_implementation": REPEATS,
        "baseline_initialize_ms": baseline_initialize_ms,
        "cached_initialize_ms": cached_initialize_ms,
        "baseline_median_ms": 1e3 * baseline_median,
        "cached_median_ms": 1e3 * cached_median,
        "baseline_samples_ms": [1e3 * value for value in baseline_times],
        "cached_samples_ms": [1e3 * value for value in cached_times],
        "cached_over_baseline": cached_median / baseline_median,
        "cache_pickle_bytes": len(pickle.dumps(caches, protocol=pickle.HIGHEST_PROTOCOL)),
        "outputs_exact": baseline_values == cached_values,
        "max_abs_delta": max(abs(old - new) for old, new in zip(baseline_values, cached_values)),
        "baseline_output": baseline_values[0],
        "cached_output": cached_values[0],
        "baseline_timed_outputs": baseline_values,
        "cached_timed_outputs": cached_values,
        "differential_exact_matches": sum(case["exact"] for case in differential_cases),
        "differential_case_count": len(differential_cases),
        "differential_max_abs_delta": max(
            abs(case["baseline"] - case["cached"]) for case in differential_cases
        ),
        "differential_cases": differential_cases,
    }


def run() -> dict:
    started_at = datetime.now().astimezone()
    results = [
        _measure(
            "IPTA",
            "stiffgwpy_fast.cobaya.likelihoods.PTA.IPTA",
            "stiffgwpy_fast/cobaya/likelihoods/PTA/IPTA.py",
            "IPTA",
            {
                "sample_file": str(DATA_DIR / "EPTAdr2" / "EPTA_dr2new_mock.dat"),
                "freq_file": str(DATA_DIR / "EPTAdr2" / "freqs_dr2new.txt"),
                "Nfreqs": 9,
            },
        ),
        _measure(
            "NANOGrav",
            "stiffgwpy_fast.cobaya.likelihoods.PTA.NANOGrav",
            "stiffgwpy_fast/cobaya/likelihoods/PTA/NANOGrav.py",
            "NANOGrav",
            {
                "kde_file": str(DATA_DIR / "NANOGrav15yr" / "density_mock.npy"),
                "log10rhogrid_file": str(DATA_DIR / "NANOGrav15yr" / "log10rhogrid.npy"),
                "freq_file": str(DATA_DIR / "NANOGrav15yr" / "freqs.npy"),
                "Nfreqs": 14,
            },
        ),
    ]
    for result in results:
        if not result["outputs_exact"] or result["differential_exact_matches"] != result["differential_case_count"]:
            raise RuntimeError(f"Numerical mismatch in {result['adapter']} cache comparison")

    return {
        "date": started_at.date().isoformat(),
        "started_at": started_at.isoformat(timespec="seconds"),
        "benchmark_script_sha256": hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        "baseline_commit": BASELINE_COMMIT,
        "candidate_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "protocol": {
            "method": "alternating direct log_likelihood calls after per-instance initialize",
            "warmup_calls_per_implementation": 1,
            "repeats_per_implementation": REPEATS,
            "parameters": PARAMETERS,
            "theory_frequency_count": int(THEORY_FREQUENCIES.size),
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
            "cpu_count": os.cpu_count(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cobaya": cobaya.__version__,
            "thread_environment": {
                name: os.environ.get(name)
                for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
            },
        },
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    timestamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / f"benchmark_pta_likelihood_cache_{timestamp}.json",
        help="JSON output path",
    )
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
