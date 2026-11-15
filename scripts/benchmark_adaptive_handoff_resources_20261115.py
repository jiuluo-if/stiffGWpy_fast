"""Compare per-mode handoff across thread counts at the same frozen root/grid."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import psutil
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
REF10_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_resource_benchmark_20261115.json"
THREAD_COUNTS = (1, 2, 8)
VARIANTS = ("z5", "z8", "adaptive")
ORDERS = (VARIANTS, ("z8", "adaptive", "z5"), ("adaptive", "z5", "z8"))
PAIRS = 12


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_variant(
    params: dict[str, float],
    freqs: np.ndarray,
    root: float,
    threads: int,
    variant: str,
) -> tuple[dict[str, object], int]:
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = root
    if variant == "adaptive":
        model._fast_adaptive_handoff_below = 6.0
        model._fast_adaptive_handoff_to = 8.0
        z_tail = 5.0
    else:
        z_tail = float(variant[1:])
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        try:
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=0.99,
                z_tail=z_tail, threads=threads, accuracy_mode="fast",
                freq_grid="goal", physical_guard=False,
            )
        finally:
            FA.goal_oriented_freqs = original_goal
            FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None:
        raise RuntimeError(f"variant={variant}, threads={threads}: {model.fast_failure_reason}")
    return {
        "runtime_s": elapsed,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "root_residual_relative": float(model.DN_gw[-1] / root - 1.0),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "selected_mode_count": getattr(model, "_fast_adaptive_handoff_selected", None),
        "runtime_warnings": [str(item.message) for item in caught],
    }, int(psutil.Process().memory_info().rss)


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    ref_data = json.loads(REF10_INPUT.read_text(encoding="utf-8"))
    case = next(item for item in data["cases"] if item["frequency_count"] == 937)
    params = data["parameters"]
    freqs = np.asarray(case["frequency_log10_hz"], dtype=np.float64)
    root = float(case["frozen_DN_eff"])
    reference = np.asarray(ref_data["reference"]["log10_OmegaGW"], dtype=np.float64)
    process = psutil.Process()
    baseline_rss = int(process.memory_info().rss)
    resources = {}
    for threads in THREAD_COUNTS:
        first_calls = {}
        first_rss = {}
        for variant in VARIANTS:
            result, rss = run_variant(params, freqs, root, threads, variant)
            first_calls[variant] = result["runtime_s"]
            first_rss[variant] = rss
        timing = {variant: [] for variant in VARIANTS}
        rss_rows = {variant: [] for variant in VARIANTS}
        spectra = {variant: [] for variant in VARIANTS}
        maps = {variant: [] for variant in VARIANTS}
        selected = []
        for repeat in range(PAIRS):
            for variant in ORDERS[repeat % len(ORDERS)]:
                result, rss = run_variant(params, freqs, root, threads, variant)
                timing[variant].append(result["runtime_s"])
                rss_rows[variant].append(rss)
                spectra[variant].append(result["log10_OmegaGW"])
                maps[variant].append(result["mapped_DN_gw"])
                if variant == "adaptive":
                    selected.append(result["selected_mode_count"])
        summaries = {}
        for variant in VARIANTS:
            seconds = np.asarray(timing[variant], dtype=np.float64)
            rss_mib = np.asarray(rss_rows[variant], dtype=np.float64) / (1024.0 ** 2)
            spec = np.asarray(spectra[variant], dtype=np.float64)
            spread = np.max(np.abs(spec - spec[0]), axis=1)
            err = np.abs(spec[0] - reference)
            summaries[variant] = {
                "p50_runtime_s": float(np.median(seconds)),
                "p95_runtime_s": float(np.quantile(seconds, 0.95)),
                "p50_rss_mib": float(np.median(rss_mib)),
                "p95_rss_mib": float(np.quantile(rss_mib, 0.95)),
                "max_repeat_spectrum_spread_dex": float(np.max(spread)),
                "max_abs_spectrum_error_vs_reference_z10_dex": float(np.max(err)),
                "p95_abs_spectrum_error_vs_reference_z10_dex": float(np.quantile(err, 0.95)),
                "mapped_DN_gw_p50": float(np.median(maps[variant])),
                "root_residual_relative_p50": float(np.median(maps[variant]) / root - 1.0),
            }
        z5_p50 = summaries["z5"]["p50_runtime_s"]
        for variant in VARIANTS:
            summaries[variant]["p50_speedup_vs_z5"] = z5_p50 / summaries[variant]["p50_runtime_s"]
        resources[str(threads)] = {
            "first_call_runtime_s": first_calls,
            "first_call_rss_mib": {key: value / (1024.0 ** 2) for key, value in first_rss.items()},
            "adaptive_selected_mode_counts": sorted(set(selected)),
            "variants": summaries,
        }
    output = {
        "schema_version": 1,
        "status": "complete_adaptive_handoff_multi_resource_benchmark",
        "scope": "Same minimum-subnormal physical parameters, frozen 235-root, exact native 937-node grid; one per-variant warm call then 12 alternating triplets at each thread count.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "frozen_DN_eff": root,
        "threads": list(THREAD_COUNTS),
        "pairs_per_thread_count": PAIRS,
        "warm_order_cycle": [list(order) for order in ORDERS],
        "baseline_rss_mib": baseline_rss / (1024.0 ** 2),
        "adaptive_rule": {"z_today_below": 6.0, "defer_handoff_to": 8.0, "phase_resolved_active_rows": "not deferred"},
        "resources": resources,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": __import__("numba").__version__,
            "psutil": psutil.__version__,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "input_json": "docs/research/results/underflow_min_subnormal_reference_grid_20261115.json",
                "reference_z10_json": "docs/research/results/underflow_min_subnormal_reference_z10_20261115.json",
                "diagnostic_script": "scripts/benchmark_adaptive_handoff_resources_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "resources": {thread: {key: {
            "p50": value["p50_runtime_s"],
            "p95": value["p95_runtime_s"],
            "speedup_vs_z5": value["p50_speedup_vs_z5"],
            "max_error": value["max_abs_spectrum_error_vs_reference_z10_dex"],
        } for key, value in row["variants"].items()} for thread, row in resources.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
