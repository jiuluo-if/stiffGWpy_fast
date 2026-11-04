"""Measure the batch size at which Reference H2 JIT compilation amortizes."""

from __future__ import annotations

import hashlib
import json
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.19649188544593907,
    "cr": 0.0,
    "n_t": 1.0622501865308032,
    "T_re": 1761.314039059193,
    "kappa10": 0.011598483151630803,
    "DN_re": 10.0,
}
DN_EFF = 2.291411867022825
Z_TAIL = 5.0
Z_START = -12.0
RTOL = 1e-9
PAIR_COUNT = 3
BASELINE_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "reference_h2_numba_batch_crossover_20261105.json"


def _git(*args):
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _solve(model, freqs, evaluator):
    R._H2_and_sigma = evaluator
    return R.spectrum_reference(
        model,
        freqs,
        DN_EFF,
        z_tail=Z_TAIL,
        rtol=RTOL,
        workers=1,
        backend="cartesian",
        z_start=Z_START,
    )


def main():
    full_freqs = np.asarray(np.load(BASELINE_NPZ)["freqs"], dtype=float)
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    python_eval = R._H2_and_sigma_python
    candidate_eval = R._H2_and_sigma
    old_warm = R._warmup_numba_reference_background
    warm_s = [0.0]

    def timed_warm(bg):
        started = time.perf_counter()
        try:
            return old_warm(bg)
        finally:
            warm_s[0] += time.perf_counter() - started

    R._warmup_numba_reference_background = timed_warm
    cold = {}
    # A sub-threshold batch should compile only the existing spline evaluator.
    for size in (16, 32):
        indices = np.linspace(0, full_freqs.size - 1, size, dtype=int)
        freqs = full_freqs[indices]
        warm_s[0] = 0.0
        started = time.perf_counter()
        _solve(model, freqs, candidate_eval)
        cold_elapsed = time.perf_counter() - started
        cold[str(size)] = {
            "cold_wall_s": cold_elapsed,
            "reference_background_warm_s": warm_s[0],
            "h2_kernel_ready_after_call": bool(R._NUMBA_H2_SIGMA_READY),
            "output_count": int(freqs.size),
        }

    results = []
    for size in (16, 32, 64):
        indices = np.linspace(0, full_freqs.size - 1, size, dtype=int)
        freqs = full_freqs[indices]
        samples = {"python": [], "compiled": []}
        outputs = {"python": [], "compiled": []}
        for pair in range(PAIR_COUNT):
            order = ("python", "compiled") if pair % 2 == 0 else (
                "compiled", "python"
            )
            for method in order:
                evaluator = python_eval if method == "python" else candidate_eval
                started = time.perf_counter()
                result = _solve(model, freqs, evaluator)
                samples[method].append(time.perf_counter() - started)
                outputs[method].append(result)
        exact = all(
            np.array_equal(left[j], right[j])
            for left, right in zip(outputs["python"], outputs["compiled"])
            for j in (0, 1, 2, 3)
        )
        if not exact:
            raise AssertionError(f"output mismatch at batch size {size}")
        py_p50 = statistics.median(samples["python"])
        nb_p50 = statistics.median(samples["compiled"])
        cold_sample = cold.get(str(size))
        results.append({
            "batch_size": size,
            "frequency_sha256": hashlib.sha256(
                np.asarray(freqs, dtype="<f8").tobytes()
            ).hexdigest(),
            "python_runtime_s": samples["python"],
            "compiled_runtime_s": samples["compiled"],
            "python_p50_s": py_p50,
            "python_p95_s": float(np.percentile(samples["python"], 95)),
            "compiled_p50_s": nb_p50,
            "compiled_p95_s": float(np.percentile(samples["compiled"], 95)),
            "warm_candidate_over_python_p50": nb_p50 / py_p50,
            "cold_batch_runtime_s": None if cold_sample is None else cold_sample[
                "cold_wall_s"
            ],
            "cold_candidate_vs_warm_python_delta_s": (
                None if size < R._NUMBA_H2_SIGMA_MIN_BATCH_SIZE or cold_sample is None
                else cold_sample["cold_wall_s"] - py_p50
            ),
            "bitwise_outputs_match": exact,
        })

    R._warmup_numba_reference_background = old_warm
    R._H2_and_sigma = candidate_eval
    try:
        import numba
        import scipy
        numba_version, scipy_version = numba.__version__, scipy.__version__
    except ImportError:
        numba_version = scipy_version = None
    result = {
        "schema_version": 1,
        "experiment": "reference_h2_numba_batch_crossover",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "baseline_archive_sha256": _sha(BASELINE_NPZ),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "backend": "cartesian",
        "workers": 1,
        "cold_path": cold,
        "paired_runs": PAIR_COUNT,
        "batch_sizes": results,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy_version,
            "numba": numba_version,
            "workers": 1,
        },
        "limitations": [
            "Only one difficult cosmology/root and three batch sizes were checked.",
            "The crossover is a workload dispatch guide, not a precision boundary.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
