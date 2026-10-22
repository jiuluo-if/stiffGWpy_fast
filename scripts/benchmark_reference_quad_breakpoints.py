"""Compare reference-spectrum quadrature breakpoint choices on fixed spectra."""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad as scipy_quad
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
reference = importlib.import_module("stiffgwpy_fast.reference")
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG


OUTPUT = ROOT / "docs" / "benchmark_reference_quad_breakpoints_20261022.json"
BASELINE_COMMIT = "7579a3e7ca60372580d668ad320be6b777ef0d4f"
PAIR_COUNT = 40
BOOTSTRAP_COUNT = 10_000


def previous_integrate_spectrum(freqs, ogw, oj):
    """Reproduce HEAD behavior: use fine-grid breakpoints for both splines."""
    freqs = np.asarray(freqs, dtype=float)
    integrand = np.asarray(ogw, dtype=float) - np.asarray(oj, dtype=float)
    order = np.argsort(freqs)
    x, y = freqs[order], integrand[order]
    spline = PchipInterpolator(x, y)

    def quadrature(func):
        return scipy_quad(
            func,
            float(x[0]),
            float(x[-1]),
            epsabs=0.0,
            epsrel=1e-10,
            limit=800,
            points=x,
        )

    value, abs_error = quadrature(lambda point: float(spline(point)))
    half = x[::2]
    half_spline = PchipInterpolator(half, y[::2])
    half_value, _ = quadrature(lambda point: float(half_spline(point)))
    ln10 = math.log(10.0)
    return (
        float(ln10 * value),
        float(ln10 * abs(abs_error)),
        float(ln10 * abs(value - half_value)),
    )


def _count_previous_calls(freqs, values):
    count = [0]
    original_quad = previous_integrate_spectrum.__globals__["scipy_quad"]

    def counted_quad(func, *args, **kwargs):
        def counted_func(point, *func_args):
            count[0] += 1
            return func(point, *func_args)

        return original_quad(counted_func, *args, **kwargs)

    previous_integrate_spectrum.__globals__["scipy_quad"] = counted_quad
    try:
        result = previous_integrate_spectrum(freqs, values, np.zeros_like(values))
    finally:
        previous_integrate_spectrum.__globals__["scipy_quad"] = original_quad
    return result, count[0]


def _count_candidate_calls(freqs, values):
    count = [0]
    original_quad = reference.quad

    def counted_quad(func, *args, **kwargs):
        def counted_func(point, *func_args):
            count[0] += 1
            return func(point, *func_args)

        return original_quad(counted_func, *args, **kwargs)

    reference.quad = counted_quad
    try:
        result = reference.integrate_spectrum(freqs, values, np.zeros_like(values))
    finally:
        reference.quad = original_quad
    return result, count[0]


def _cases():
    powerlaw_x = np.linspace(-18.0, 6.0, 401)
    yield "powerlaw-401", powerlaw_x, 10.0 ** (-14.0 - 0.3 * powerlaw_x)

    peaked_x = np.linspace(-20.0, 7.0, 247)
    peaked_y = np.exp(-0.23 * (peaked_x + 4.0)) * (
        1.0
        + 0.17 * np.exp(-0.5 * ((peaked_x + 5.7) / 0.21) ** 2)
        + 0.08 * np.exp(-0.5 * ((peaked_x - 0.6) / 0.48) ** 2)
        + 0.012 * np.sin(3.7 * peaked_x)
    )
    yield "peaked-247", peaked_x, peaked_y

    rng = np.random.default_rng(20261022)
    irregular_x = np.sort(rng.uniform(-22.0, 8.0, 240))
    irregular_y = np.exp(-0.18 * (irregular_x + 2.0)) * (
        1.0
        + 0.2 * np.exp(-0.5 * ((irregular_x + 8.0) / 0.3) ** 2)
        + 0.05 * np.sin(2.1 * irregular_x)
    )
    yield "irregular-240", irregular_x, irregular_y

    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    model.SGWB_iter(engine="fast", accuracy_mode="fast", tol=1e-7)
    yield (
        "formal-fast-default-76",
        np.asarray(model.f, dtype=float),
        np.asarray(model.Ogw_today - model.Oj_today, dtype=float),
    )


def _source_sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _head_source_sha():
    content = subprocess.check_output(
        ["git", "show", f"{BASELINE_COMMIT}:stiffgwpy_fast/reference.py"],
        cwd=ROOT,
    )
    return hashlib.sha256(content).hexdigest()


def _run_case(name, freqs, values):
    zeros = np.zeros_like(values)
    previous_result = previous_integrate_spectrum(freqs, values, zeros)
    candidate_result = reference.integrate_spectrum(freqs, values, zeros)
    _, previous_calls = _count_previous_calls(freqs, values)
    _, candidate_calls = _count_candidate_calls(freqs, values)

    # Warm both paths, then alternate order to reduce time-drift bias.
    for _ in range(3):
        previous_integrate_spectrum(freqs, values, zeros)
        reference.integrate_spectrum(freqs, values, zeros)

    previous_ms = []
    candidate_ms = []
    for pair in range(PAIR_COUNT):
        variants = ("previous", "candidate")
        if pair % 2:
            variants = variants[::-1]
        for variant in variants:
            start = time.perf_counter()
            if variant == "previous":
                previous_integrate_spectrum(freqs, values, zeros)
            else:
                reference.integrate_spectrum(freqs, values, zeros)
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            (previous_ms if variant == "previous" else candidate_ms).append(
                elapsed_ms
            )

    ratios = np.asarray(candidate_ms) / np.asarray(previous_ms)
    rng = np.random.default_rng(8401)
    samples = ratios[
        rng.integers(0, ratios.size, size=(BOOTSTRAP_COUNT, ratios.size))
    ]
    bootstrap_medians = np.median(samples, axis=1)
    return {
        "case": name,
        "point_count": int(freqs.size),
        "input_sha256": hashlib.sha256(
            np.column_stack((freqs, values)).astype("<f8").tobytes()
        ).hexdigest(),
        "previous_function_calls": previous_calls,
        "candidate_function_calls": candidate_calls,
        "function_call_reduction_fraction": 1.0 - candidate_calls / previous_calls,
        "previous_median_ms": statistics.median(previous_ms),
        "candidate_median_ms": statistics.median(candidate_ms),
        "paired_candidate_over_previous_median": float(np.median(ratios)),
        "paired_ratio_bootstrap_95pct": [
            float(value)
            for value in np.quantile(bootstrap_medians, [0.025, 0.975])
        ],
        "previous_ms_samples": previous_ms,
        "candidate_ms_samples": candidate_ms,
        "absolute_output_deltas": {
            key: abs(candidate_result[index] - previous_result[index])
            for index, key in enumerate(
                ("g2", "quadrature_error", "interpolation_error")
            )
        },
        "previous_outputs": dict(
            zip(("g2", "quadrature_error", "interpolation_error"), previous_result)
        ),
        "candidate_outputs": dict(
            zip(("g2", "quadrature_error", "interpolation_error"), candidate_result)
        ),
    }


def main():
    try:
        affinity = sorted(__import__("os").sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity = None
    result = {
        "created_at": datetime.now().astimezone().isoformat(),
        "workload": "isolated integrate_spectrum calls on synthetic spectra and one formal-fast spectrum",
        "candidate": "use half-density PCHIP knots as coarse quadrature breakpoints",
        "scope_note": "Subroutine-only timing; not an end-to-end reference or fast solver speedup.",
        "formal_fast_spectrum_generation": {
            "engine": "fast",
            "accuracy_mode": "fast",
            "tol": 1e-7,
            "parameters": {"r": 1e-2, "cr": 1, "T_re": 2e3, "kappa10": 1e-2},
        },
        "previous_source_commit": BASELINE_COMMIT,
        "previous_reference_source_sha256": _head_source_sha(),
        "candidate_reference_source_sha256": _source_sha(
            ROOT / "stiffgwpy_fast" / "reference.py"
        ),
        "model_wrapper_source_sha256": _source_sha(
            ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"
        ),
        "fast_solver_source_sha256": _source_sha(
            ROOT / "stiffgwpy_fast" / "fast_sgwb.py"
        ),
        "script_sha256": _source_sha(__file__),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "logical_cpu_count": __import__("os").cpu_count(),
            "cpu_affinity": affinity,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": importlib.import_module("numba").__version__,
        },
        "pair_count": PAIR_COUNT,
        "bootstrap_count": BOOTSTRAP_COUNT,
        "cases": [
            _run_case(name, freqs, values)
            for name, freqs, values in _cases()
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    for case in result["cases"]:
        print(
            case["case"],
            f"calls {case['previous_function_calls']} -> "
            f"{case['candidate_function_calls']}",
            f"paired ratio {case['paired_candidate_over_previous_median']:.4f}",
            f"CI {case['paired_ratio_bootstrap_95pct']}",
            f"output deltas {case['absolute_output_deltas']}",
        )


if __name__ == "__main__":
    main()
