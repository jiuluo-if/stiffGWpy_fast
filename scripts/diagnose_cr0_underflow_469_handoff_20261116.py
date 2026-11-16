"""Run same-root adaptive WKB A/B on the direct 469-node underflow holdout."""

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
sys.path.insert(0, str(ROOT / "scripts"))

import diagnose_cr_positive_consistency_tail_handoff_20261116 as handoff  # noqa: E402

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr0_underflow_469_handoff_20261116.json"
THREADS = 1
REPEATS = 12


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_self_consistent(params: dict[str, float], freqs: np.ndarray) -> dict[str, object]:
    model = LCDM_SG(**params)
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid
    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    start = time.perf_counter()
    caught: list[warnings.WarningMessage]
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                threads=THREADS, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    if result is not model:
        raise RuntimeError(f"current Fast 469-node root failed: {model.fast_failure_reason}")
    return {
        "runtime_s": time.perf_counter() - start,
        "DN_eff_root": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "fixed_point_residual_relative": float(model.DN_gw[-1] / model.cosmo_param["DN_eff"] - 1.0),
        "warnings": [str(item.message) for item in caught],
    }


def summarize(records: list[dict[str, object]], ref: np.ndarray, freq: np.ndarray) -> dict[str, object]:
    spectrum = np.asarray(records[0]["log10_OmegaGW"], dtype=np.float64)
    delta = spectrum - ref
    runtimes = np.asarray([float(item["runtime_s"]) for item in records])
    rss = np.asarray([float(item["rss_bytes"]) for item in records]) / (1024.0 ** 2)
    spreads = [
        float(np.max(np.abs(np.asarray(item["log10_OmegaGW"]) - spectrum)))
        for item in records
    ]
    worst = int(np.argmax(np.abs(delta)))
    z_today_rows = [item["z_today"] for item in records if item["z_today"] is not None]
    return {
        "p50_runtime_s": float(np.median(runtimes)),
        "p95_runtime_s": float(np.quantile(runtimes, 0.95)),
        "p50_rss_mib": float(np.median(rss)),
        "p95_rss_mib": float(np.quantile(rss, 0.95)),
        "max_repeat_spectrum_spread_dex": max(spreads),
        "max_abs_spectrum_error_vs_reference_z8_dex": float(np.max(np.abs(delta))),
        "p95_abs_spectrum_error_vs_reference_z8_dex": float(np.percentile(np.abs(delta), 95.0)),
        "worst_log10_frequency": float(freq[worst]),
        "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
        "mapped_DN_gw": float(records[0]["mapped_DN_gw"]),
        "root_residual_relative_vs_frozen_root": float(records[0]["root_residual_relative"]),
        "handoff_count": int(records[0]["handoff_count"]),
        "adaptive_selected_count": int(records[0]["adaptive_selected_count"]),
        "adaptive_selected_fraction": float(records[0]["adaptive_selected_fraction"]),
        "z_today_min": float(min(np.min(row) for row in z_today_rows)) if z_today_rows else None,
        "z_today_max": float(max(np.max(row) for row in z_today_rows)) if z_today_rows else None,
        "warnings": sorted({warning for item in records for warning in item["warnings"]}),
        "runtime_samples_s": runtimes.tolist(),
    }


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    case = next(item for item in source["cases"] if item["frequency_count"] == 469)
    params = {key: float(value) for key, value in source["parameters"].items()}
    freq = np.asarray(case["frequency_log10_hz"], dtype=np.float64)
    reference = case["reference"]
    ref_spectrum = np.asarray(reference["log10_OmegaGW"], dtype=np.float64)
    grid_hash = hashlib.sha256(np.asarray(freq, dtype="<f8").tobytes()).hexdigest()
    if freq.size != 469 or grid_hash != case["frequency_sha256_f64le"]:
        raise RuntimeError("direct 469-node Reference grid hash mismatch")
    if source["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("direct Reference source SHA differs from current Reference")
    if ref_spectrum.shape != freq.shape or not np.isfinite(ref_spectrum).all():
        raise RuntimeError("direct 469-node Reference z8 spectrum is incomplete")

    frozen_root = float(case["frozen_DN_eff"])
    root_result = run_self_consistent(params, freq)
    process = psutil.Process()
    first_calls = {
        variant: handoff.run_variant(variant, params, freq, frozen_root, process)
        for variant in handoff.VARIANTS
    }
    records = {variant: [] for variant in handoff.VARIANTS}
    for repeat in range(REPEATS):
        for variant in handoff.ORDERS[repeat % len(handoff.ORDERS)]:
            records[variant].append(handoff.run_variant(variant, params, freq, frozen_root, process))
    summaries = {
        variant: summarize(records[variant], ref_spectrum, freq)
        for variant in handoff.VARIANTS
    }
    for variant in handoff.VARIANTS:
        if summaries[variant]["max_repeat_spectrum_spread_dex"] != 0.0:
            raise RuntimeError(f"{variant} spectrum changed across identical fixed-root repeats")

    output = {
        "schema_version": 1,
        "status": "complete_current_fast_469_handoff_holdout",
        "scope": "Current Fast same-root/grid/thread z5, adaptive, and z8 A/B against directly solved Reference z8 on the exact 469-node high-tilt minimum-subnormal cr=0 grid. Also run current Fast self-consistency on that same direct grid.",
        "parameters": params,
        "derived_n_t": params["n_t"],
        "frequency_count": int(freq.size),
        "frequency_sha256_f64le": grid_hash,
        "frozen_reference_DN_eff": frozen_root,
        "reference_z8": {
            "mapped_DN_gw": float(reference["mapped_DN_gw"]),
            "map_relative_residual": float(reference["map_relative_residual"]),
            "runtime_warnings": reference["runtime_warnings"],
            "normalized_nodes_underflowed_to_zero": int(reference["normalized_nodes_underflowed_to_zero"]),
            "quadrature_error_relative": reference["quadrature_error_relative"],
            "interpolation_error_relative": reference["interpolation_error_relative"],
        },
        "current_fast_self_consistent": root_result,
        "root_relative_difference_vs_reference_frozen_root": root_result["DN_eff_root"] / frozen_root - 1.0,
        "root_relative_difference_vs_reference_z8_map": (
            root_result["DN_eff_root"] / float(reference["mapped_DN_gw"]) - 1.0
        ),
        "runtime_protocol": {
            "threads": THREADS,
            "warmup_calls_per_variant": 1,
            "paired_repeats_per_variant": REPEATS,
            "order_cycle": [list(row) for row in handoff.ORDERS],
            "module_import_time_excluded": True,
            "self_consistent_root_call_may_include_jit": True,
            "variant_measurements_follow_self_consistent_root_warmup": True,
        },
        "first_variant_calls_after_root_warmup": {
            variant: {
                key: value for key, value in first_calls[variant].items()
                if key not in {"log10_OmegaGW", "handoff_eps", "z_today"}
            }
            for variant in handoff.VARIANTS
        },
        "variants": summaries,
        "adaptive_vs_z8_speedup_p50": (
            summaries["z8"]["p50_runtime_s"] / summaries["adaptive"]["p50_runtime_s"]
        ),
        "adaptive_vs_z5_p50_runtime_ratio": (
            summaries["adaptive"]["p50_runtime_s"] / summaries["z5"]["p50_runtime_s"]
        ),
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": __import__("numba").__version__,
            "psutil": psutil.__version__,
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "stiff_SGWB.py": sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "grid_and_reference_input": sha(INPUT),
            "handoff_runner": sha(ROOT / "scripts" / "diagnose_cr_positive_consistency_tail_handoff_20261116.py"),
            "diagnostic_script": sha(Path(__file__)),
        },
        "spectra": {
            "reference_z8_log10_OmegaGW": ref_spectrum.tolist(),
            **{variant: np.asarray(first_calls[variant]["log10_OmegaGW"], dtype=np.float64).tolist() for variant in handoff.VARIANTS},
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "self_consistent_root_relative_difference": output["root_relative_difference_vs_reference_frozen_root"],
        "variants": {key: {
            "max_abs_spectrum_error_dex": value["max_abs_spectrum_error_vs_reference_z8_dex"],
            "p50_runtime_s": value["p50_runtime_s"],
            "p95_runtime_s": value["p95_runtime_s"],
            "selected": value["adaptive_selected_count"],
        } for key, value in summaries.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
