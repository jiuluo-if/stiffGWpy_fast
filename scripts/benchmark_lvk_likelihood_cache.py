"""Measure cached LVK likelihood inputs and NumPy chi-square reduction.

The driver compares the previous per-call data preparation and Python reduction
with the current adapter on deterministic synthetic spectra and the packaged
54,593-bin LVK table. It measures only ``log_likelihood`` calls, not full Cobaya
or MCMC runtime.

中文：使用随包 LVK 数据表和确定性合成谱，对比似然中的固定数据缓存及卡方归约优化。
只测量 ``log_likelihood`` 调用，不代表完整 Cobaya 或 MCMC 的加速。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402
from scipy import interpolate  # noqa: E402

from stiffgwpy_fast.cobaya.likelihoods.LIGO_SGWB.LVK_SGWB_CC import (  # noqa: E402
    LVK_SGWB_CC,
)

DATA_PATH = ROOT / "stiffgwpy_fast/cobaya/likelihoods/LIGO_SGWB/C_O1_O2_O3.dat"


def previous_log_likelihood(data, f_theory, omega_theory):
    """Reproduce the previous per-call prep and Python reduction.

    中文：重现优化前每次计算都重算数据频率并用 Python ``sum`` 归约的路径。
    """
    f_lvk = np.log10(data[:, 0])
    cf_lvk = data[:, 1]
    sigma_lvk = data[:, 2]
    omega_model = np.zeros_like(cf_lvk)
    if f_theory[-1] >= f_lvk[0]:
        support = f_theory >= -5.0
        spectrum = interpolate.interp1d(
            f_theory[support], omega_theory[support], kind="cubic")
        supported_bins = f_lvk <= f_theory[-1]
        omega_model[supported_bins] = 10.0 ** spectrum(f_lvk[supported_bins])
    chi2 = np.square(np.divide(cf_lvk - omega_model, sigma_lvk))
    return -sum(chi2) / 2.0, chi2


def _time_call(function):
    start = time.perf_counter_ns()
    value = function()
    return value, time.perf_counter_ns() - start


def _case_spectra():
    """Yield deterministic spectra with different upper frequency support."""
    for index, f_end in enumerate((1.0, 1.31, 1.8, 2.7, 3.4, 4.0)):
        f_theory = np.linspace(-5.0, f_end, 700)
        slope = 0.7 + 0.11 * index
        omega_theory = -8.0 + 0.3 * np.sin(slope * f_theory)
        yield {
            "case": index,
            "f_end_log10_hz": f_end,
            "n_theory": int(f_theory.size),
            "f_theory": f_theory,
            "omega_theory": omega_theory,
        }


def benchmark(repeats=200, warmups=2):
    """Compare the original and optimized LVK likelihood on fixed inputs.

    中文：交错测量旧实现与新实现，并记录似然差值、运行环境和逐次耗时。
    """
    data = np.loadtxt(DATA_PATH)
    likelihood = LVK_SGWB_CC.__new__(LVK_SGWB_CC)
    likelihood.CC_file = str(DATA_PATH)
    likelihood.initialize()

    truncated_f = np.linspace(1.5, 4.0, 700)
    truncated_omega = np.zeros_like(truncated_f)
    try:
        previous_log_likelihood(data, truncated_f, truncated_omega)
    except ValueError as exc:
        old_lower_support_error = type(exc).__name__
    else:
        old_lower_support_error = None
    new_truncated_loglike = likelihood.log_likelihood(
        truncated_f, truncated_omega)

    rows = []
    for case in _case_spectra():
        f_theory = case.pop("f_theory")
        omega_theory = case.pop("omega_theory")

        def baseline_call():
            return previous_log_likelihood(data, f_theory, omega_theory)[0]

        def optimized_call():
            return likelihood.log_likelihood(f_theory, omega_theory)

        old_value, old_chi2 = previous_log_likelihood(data, f_theory, omega_theory)
        new_value = optimized_call()
        fsum_value = -math.fsum(float(value) for value in old_chi2) / 2.0
        old_ns = []
        new_ns = []
        for _ in range(warmups):
            baseline_call()
            optimized_call()
        for pair in range(repeats):
            calls = ((optimized_call, new_ns), (baseline_call, old_ns)) if pair % 2 else (
                (baseline_call, old_ns), (optimized_call, new_ns))
            for function, samples in calls:
                value, elapsed_ns = _time_call(function)
                samples.append(elapsed_ns)
                if function is baseline_call:
                    old_value = value
                else:
                    new_value = value

        old_median_ns = float(np.median(old_ns))
        new_median_ns = float(np.median(new_ns))
        rows.append({
            **case,
            "old_loglike": float(previous_log_likelihood(
                data, f_theory, omega_theory)[0]),
            "optimized_loglike": float(optimized_call()),
            "abs_loglike_delta": float(abs(old_value - new_value)),
            "old_abs_error_vs_fsum": float(abs(old_value - fsum_value)),
            "optimized_abs_error_vs_fsum": float(abs(new_value - fsum_value)),
            "old_samples_ns": old_ns,
            "optimized_samples_ns": new_ns,
            "old_median_ms": old_median_ns / 1e6,
            "optimized_median_ms": new_median_ns / 1e6,
            "optimized_over_old": new_median_ns / old_median_ns,
        })

    return {
        "commit": os.popen("git rev-parse HEAD").read().strip(),
        "data_path": DATA_PATH.relative_to(ROOT).as_posix(),
        "data_sha256": hashlib.sha256(DATA_PATH.read_bytes()).hexdigest(),
        "data_rows": int(data.shape[0]),
        "resources": telemetry(workers=1, threads=None),
        "lower_support_probe": {
            "theory_range_log10_hz": [1.5, 4.0],
            "previous_error": old_lower_support_error,
            "optimized_loglike_finite": bool(np.isfinite(new_truncated_loglike)),
        },
        "protocol": {
            "warmups_per_method_case": warmups,
            "alternating_pairs_per_case": repeats,
            "timer": "time.perf_counter_ns around log_likelihood call only",
            "cases": "six deterministic sinusoidal log10(Omega) spectra with varied upper frequency support",
            "old_reduction": "Python sum over float64 chi-square array",
            "optimized_reduction": "numpy.sum(..., dtype=float64), pairwise reduction",
        },
        "results": rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Benchmark the LVK likelihood's cached inputs and reduction")
    parser.add_argument("--repeats", type=int, default=200)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--json", default=None,
                        help="output JSON path (default: timestamped file under docs/)")
    args = parser.parse_args(argv)
    if args.repeats < 1 or args.warmups < 0:
        parser.error("--repeats must be positive and --warmups non-negative")
    payload = benchmark(args.repeats, args.warmups)
    output = (Path(args.json) if args.json else ROOT / (
        "docs/benchmark_lvk_likelihood_cache_%s.json"
        % time.strftime("%Y%m%d_%H%M%S")))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    print(json.dumps({
        "output": str(output),
        "commit": payload["commit"],
        "results": [{key: row[key] for key in (
            "case", "old_median_ms", "optimized_median_ms",
            "optimized_over_old", "abs_loglike_delta")}
            for row in payload["results"]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
