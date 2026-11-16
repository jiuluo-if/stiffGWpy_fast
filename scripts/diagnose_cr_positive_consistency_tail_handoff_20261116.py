"""Diagnose shallow WKB handoff error on the cr>0 consistency branch."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

DEFAULT_INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_consistency_cutoff_20261116.json"
DEFAULT_OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_consistency_tail_handoff_20261116.json"
THREADS = 1
REPEATS = 12
VARIANTS = ("z5", "z8", "adaptive")
ORDERS = (VARIANTS, ("z8", "adaptive", "z5"), ("adaptive", "z5", "z8"))
DEFER_TO = 8.0
RIPPLE_TARGET_DEX = 1e-3


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def ripple_bound_dex(z: float) -> float:
    return math.log1p(math.exp(-z) + math.exp(-2.0 * z)) / math.log(10.0)


def ripple_threshold(target_dex: float) -> float:
    low, high = -10.0, 50.0
    for _ in range(100):
        mid = 0.5 * (low + high)
        if ripple_bound_dex(mid) > target_dex:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def run_variant(
    name: str,
    params: dict[str, float],
    freqs: np.ndarray,
    root: float,
    rss_probe: psutil.Process,
) -> dict[str, object]:
    model = LCDM_SG(**{**params, "DN_eff": root})
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    original_thresholds = FS._adaptive_wkb_tail_thresholds
    z_today_samples: list[np.ndarray] = []

    def locked_background(*args, **kwargs):
        model.cosmo_param["DN_eff"] = root
        return original_gen(*args, **kwargs)

    def threshold_spy(z_today, z_tail, defer_below, defer_to, phase_weight=None):
        z_today_samples.append(np.asarray(z_today, dtype=np.float64).copy())
        return original_thresholds(z_today, z_tail, defer_below, defer_to, phase_weight)

    model.gen_expansion = locked_background
    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    if name == "adaptive":
        model._fast_adaptive_handoff_below = ripple_threshold(RIPPLE_TARGET_DEX)
        model._fast_adaptive_handoff_to = DEFER_TO
        FS._adaptive_wkb_tail_thresholds = threshold_spy
        z_tail = 5.0
    else:
        model._fast_adaptive_handoff_below = -1.0e300
        model._fast_adaptive_handoff_to = float(name[1:]) + 3.0
        FS._adaptive_wkb_tail_thresholds = original_thresholds
        z_tail = float(name[1:])

    started = time.perf_counter()
    caught: list[warnings.WarningMessage]
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=0.99, z_tail=z_tail,
                threads=THREADS, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
        FS._adaptive_wkb_tail_thresholds = original_thresholds
    runtime = time.perf_counter() - started
    if result is not model:
        raise RuntimeError(f"Fast {name} failed: {model.fast_failure_reason}")

    return {
        "variant": name,
        "runtime_s": runtime,
        "rss_bytes": int(rss_probe.memory_info().rss),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "root_residual_relative": float(model.DN_gw[-1] / root - 1.0),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "handoff_eps": np.asarray(model.handoff_eps, dtype=np.float64),
        "handoff_count": int(np.count_nonzero(np.asarray(model.handoff_eps) >= 0.0)),
        "adaptive_selected_count": getattr(model, "_fast_adaptive_handoff_selected", 0),
        "adaptive_selected_fraction": getattr(model, "_fast_adaptive_handoff_fraction", 0.0),
        "z_today": z_today_samples[-1] if z_today_samples else None,
        "warnings": [str(item.message) for item in caught],
    }


def summarize(records: list[dict[str, object]], reference: np.ndarray, freqs: np.ndarray) -> dict[str, object]:
    first = records[0]
    spectrum = np.asarray(first["log10_OmegaGW"], dtype=np.float64)
    delta = spectrum - reference
    worst = int(np.argmax(np.abs(delta)))
    runtimes = np.asarray([float(record["runtime_s"]) for record in records], dtype=np.float64)
    rss_mib = np.asarray([int(record["rss_bytes"]) for record in records], dtype=np.float64) / (1024.0 ** 2)
    spreads = [
        float(np.max(np.abs(np.asarray(record["log10_OmegaGW"]) - spectrum)))
        for record in records
    ]
    z_today = first["z_today"]
    return {
        "p50_runtime_s": float(np.median(runtimes)),
        "p95_runtime_s": float(np.quantile(runtimes, 0.95)),
        "p50_rss_mib": float(np.median(rss_mib)),
        "p95_rss_mib": float(np.quantile(rss_mib, 0.95)),
        "max_repeat_spectrum_spread_dex": float(max(spreads)),
        "max_abs_spectrum_error_vs_reference_z12_dex": float(np.max(np.abs(delta))),
        "p95_abs_spectrum_error_vs_reference_z12_dex": float(np.percentile(np.abs(delta), 95.0)),
        "worst_log10f": float(freqs[worst]),
        "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
        "mapped_DN_gw": float(first["mapped_DN_gw"]),
        "root_residual_relative": float(first["root_residual_relative"]),
        "handoff_count": int(first["handoff_count"]),
        "adaptive_selected_count": int(first["adaptive_selected_count"]),
        "adaptive_selected_fraction": float(first["adaptive_selected_fraction"]),
        "z_today_min": float(np.min(z_today)) if z_today is not None else None,
        "z_today_max": float(np.max(z_today)) if z_today is not None else None,
        "warnings": sorted(set(w for record in records for w in record["warnings"])),
        "runtime_samples_s": runtimes.tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    source = json.loads(input_path.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in source["parameters"].items()}
    root = float(source["current_fast"]["root_DN_eff"])
    freqs = np.asarray(source["frequency_grid"]["log10f"], dtype=np.float64)
    reference = np.asarray(
        source["direct_reference"]["results"]["z12"]["raw_log10_OmegaGW"],
        dtype=np.float64,
    )
    if freqs.size != 77 or digest(freqs) != source["frequency_grid"]["log10f_sha256_f64le"]:
        raise RuntimeError("the saved direct Fast grid hash or expected node count changed")
    if reference.shape != freqs.shape or not np.isfinite(reference).all():
        raise RuntimeError("the direct Reference z12 spectrum is incomplete")

    rss_probe = psutil.Process()
    rss_before = int(rss_probe.memory_info().rss)
    first_calls = {}
    first_results = {}
    for variant in VARIANTS:
        record = run_variant(variant, params, freqs, root, rss_probe)
        first_calls[variant] = record
        first_results[variant] = record

    records = {variant: [] for variant in VARIANTS}
    for repeat in range(REPEATS):
        for variant in ORDERS[repeat % len(ORDERS)]:
            records[variant].append(run_variant(variant, params, freqs, root, rss_probe))

    summaries = {
        variant: summarize(records[variant], reference, freqs)
        for variant in VARIANTS
    }
    threshold = ripple_threshold(RIPPLE_TARGET_DEX)
    output = {
        "schema_version": 1,
        "status": "complete_cr_positive_consistency_fixed_root_tail_handoff_diagnostic",
        "scope": "Same physical cr>0 consistency-branch point, same frozen root, exact saved 77-node grid and one thread; compare default z5, uniform z8, and private per-mode deferral with a ripple-derived z_today threshold. The requested n_t input does not control derived nt on this branch. Reference z12 was directly computed by the hashed input run.",
        "input_artifact": str(input_path.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "frozen_DN_eff": root,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": digest(freqs),
        "reference_z12_map_DN_gw": float(source["direct_reference"]["results"]["z12"]["DN_gw"]),
        "reference_z8_to_z10_map_relative_change": source["direct_reference"]["z8_to_z10_map_relative_change"],
        "reference_z10_to_z12_map_relative_change": source["direct_reference"]["z10_to_z12_map_relative_change"],
        "reference_z10_to_z12_max_abs_dex": source["direct_reference"]["z10_to_z12_max_abs_dex_on_positive_nodes"],
        "ripple_criterion": {
            "formula": "log10(1 + exp(-z_today) + exp(-2*z_today))",
            "target_dex": RIPPLE_TARGET_DEX,
            "derived_z_today_threshold": threshold,
            "bound_at_threshold_dex": ripple_bound_dex(threshold),
            "defer_to_z_tail": DEFER_TO,
        },
        "first_solver_call_after_module_import": {
            variant: {
                key: value
                for key, value in record.items()
                if key not in ("log10_OmegaGW", "handoff_eps", "z_today")
            }
            for variant, record in first_calls.items()
        },
        "runtime_protocol": {
            "threads": THREADS,
            "warmup_calls_per_variant": 1,
            "paired_repeats_per_variant": REPEATS,
            "order_cycle": [list(order) for order in ORDERS],
            "module_import_time_excluded": True,
            "first_solver_call_may_include_jit": True,
            "rss_before_mib": rss_before / (1024.0 ** 2),
        },
        "variants": summaries,
        "adaptive_vs_uniform_z8_speedup_p50": summaries["z8"]["p50_runtime_s"] / summaries["adaptive"]["p50_runtime_s"],
        "adaptive_vs_z5_p50_runtime_ratio": summaries["adaptive"]["p50_runtime_s"] / summaries["z5"]["p50_runtime_s"],
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
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "direct_reference_input": str(input_path.relative_to(ROOT)).replace("\\", "/"),
                "diagnostic_script": "scripts/diagnose_cr_positive_consistency_tail_handoff_20261116.py",
            }.items()
        },
        "spectra": {
            "reference_z12_raw_log10_OmegaGW": reference.tolist(),
            "first_call_fast_spectra": {
                variant: np.asarray(record["log10_OmegaGW"], dtype=np.float64).tolist()
                for variant, record in first_results.items()
            },
            "steady_state_first_spectra": {
                variant: np.asarray(records[variant][0]["log10_OmegaGW"], dtype=np.float64).tolist()
                for variant in VARIANTS
            },
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "threshold": threshold,
        "adaptive_speedup_vs_z8_p50": output["adaptive_vs_uniform_z8_speedup_p50"],
        "adaptive_vs_z5_p50_ratio": output["adaptive_vs_z5_p50_runtime_ratio"],
        "variants": {
            name: {
                "p50_s": summary["p50_runtime_s"],
                "p95_s": summary["p95_runtime_s"],
                "max_dex": summary["max_abs_spectrum_error_vs_reference_z12_dex"],
                "p95_dex": summary["p95_abs_spectrum_error_vs_reference_z12_dex"],
                "rss_p50_mib": summary["p50_rss_mib"],
                "rss_p95_mib": summary["p95_rss_mib"],
            }
            for name, summary in summaries.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
