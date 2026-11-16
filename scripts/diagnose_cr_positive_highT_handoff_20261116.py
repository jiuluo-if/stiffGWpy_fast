"""Test ripple-aware WKB handoff on a direct high-temperature cr>0 point."""

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

INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_cutoff_valid_highT_fixed_kappa_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_highT_handoff_20261116.json"
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
        raise RuntimeError(f"current Fast cr>0 high-T root failed: {model.fast_failure_reason}")
    return {
        "runtime_s": time.perf_counter() - start,
        "DN_eff_root": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "fixed_point_residual_relative": float(model.DN_gw[-1] / model.cosmo_param["DN_eff"] - 1.0),
        "derived_n_t": float(model.derived_param["nt"]),
        "physical_A_t": float(model.derived_param["A_t"]),
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
        "max_abs_spectrum_error_vs_reference_z12_dex": float(np.max(np.abs(delta))),
        "p95_abs_spectrum_error_vs_reference_z12_dex": float(np.percentile(np.abs(delta), 95.0)),
        "worst_log10_frequency": float(freq[worst]),
        "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
        "mapped_DN_gw": float(records[0]["mapped_DN_gw"]),
        "map_relative_error_vs_reference_z12": float(records[0]["mapped_DN_gw"] / records[0]["reference_DN_gw_z12"] - 1.0),
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
    params = {key: float(value) for key, value in source["parameters"].items()}
    freq = np.asarray(source["frequency_grid"]["log10f"], dtype=np.float64)
    direct = source["direct_reference"]
    root = float(direct["fixed_DN_eff"])
    ref_spectrum = np.asarray(direct["z12_raw_log10_OmegaGW"], dtype=np.float64)
    grid_hash = hashlib.sha256(np.asarray(freq, dtype="<f8").tobytes()).hexdigest()
    if freq.size != 77 or grid_hash != source["frequency_grid"]["log10f_sha256_f64le"]:
        raise RuntimeError("direct cr>0 77-node grid hash mismatch")
    if source["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("direct Reference source SHA differs from current Reference")
    if ref_spectrum.shape != freq.shape or not np.isfinite(ref_spectrum).all():
        raise RuntimeError("direct Reference z12 raw spectrum is incomplete")

    root_result = run_self_consistent(params, freq)
    current_model = LCDM_SG(**params)
    actual_nt = float(current_model.derived_param["nt"])
    physical_at_initial = float(current_model.derived_param["A_t"])
    process = psutil.Process()
    first_calls = {
        variant: handoff.run_variant(variant, params, freq, root, process)
        for variant in handoff.VARIANTS
    }
    for variant in handoff.VARIANTS:
        first_calls[variant]["reference_DN_gw_z12"] = float(direct["z12"]["DN_gw"])
    records = {variant: [] for variant in handoff.VARIANTS}
    for repeat in range(REPEATS):
        for variant in handoff.ORDERS[repeat % len(handoff.ORDERS)]:
            item = handoff.run_variant(variant, params, freq, root, process)
            item["reference_DN_gw_z12"] = float(direct["z12"]["DN_gw"])
            records[variant].append(item)
    summaries = {
        variant: summarize(records[variant], ref_spectrum, freq)
        for variant in handoff.VARIANTS
    }
    for variant in handoff.VARIANTS:
        if summaries[variant]["max_repeat_spectrum_spread_dex"] != 0.0:
            raise RuntimeError(f"{variant} spectrum changed across identical fixed-root repeats")

    output = {
        "schema_version": 1,
        "status": "complete_current_fast_cr_positive_highT_handoff_diagnostic",
        "scope": "Current Fast same-root/grid/thread z5, adaptive and z8 A/B versus directly solved Reference z12 on a cutoff-valid cr>0 high-temperature point. Reference remains independent and unchanged.",
        "parameters": params,
        "supplied_n_t": params["n_t"],
        "derived_n_t": actual_nt,
        "physical_A_t_at_initial_background": physical_at_initial,
        "frequency_count": int(freq.size),
        "frequency_sha256_f64le": grid_hash,
        "frozen_reference_DN_eff": root,
        "reference_convergence": {
            "z8_to_z10_map_relative_change": direct["z8_to_z10_map_relative_change"],
            "z10_to_z12_map_relative_change": direct["z10_to_z12_map_relative_change"],
            "z10_to_z12_max_abs_dex": direct["z10_to_z12_max_abs_dex"],
            "z12_quadrature_error": direct["z12"]["quadrature_error"],
            "z12_interpolation_error": direct["z12"]["interpolation_error"],
            "public_floor_nodes": direct["z10_public_log_spectrum_floor_nodes"],
            "raw_positive_log_spectrum_used": direct["raw_positive_log_spectrum_used_for_comparison"],
            "runtime_warnings": direct["warnings"],
        },
        "current_fast_self_consistent": root_result,
        "root_relative_difference_vs_reference_frozen_root": root_result["DN_eff_root"] / root - 1.0,
        "root_relative_difference_vs_reference_z12_map": root_result["DN_eff_root"] / float(direct["z12"]["DN_gw"]) - 1.0,
        "runtime_protocol": {
            "threads": THREADS,
            "warmup_calls_per_variant": 1,
            "paired_repeats_per_variant": REPEATS,
            "order_cycle": [list(row) for row in handoff.ORDERS],
            "module_import_time_excluded": True,
            "self_consistent_root_call_may_include_jit": True,
            "variant_measurements_follow_root_warmup": True,
        },
        "first_variant_calls_after_root_warmup": {
            variant: {
                key: value for key, value in first_calls[variant].items()
                if key not in {"log10_OmegaGW", "handoff_eps", "z_today"}
            }
            for variant in handoff.VARIANTS
        },
        "variants": summaries,
        "adaptive_vs_z8_speedup_p50": summaries["z8"]["p50_runtime_s"] / summaries["adaptive"]["p50_runtime_s"],
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
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "stiff_SGWB.py": sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "input_artifact": sha(INPUT),
            "handoff_runner": sha(ROOT / "scripts" / "diagnose_cr_positive_consistency_tail_handoff_20261116.py"),
            "diagnostic_script": sha(Path(__file__)),
        },
        "spectra": {
            "reference_z12_raw_log10_OmegaGW": ref_spectrum.tolist(),
            **{variant: np.asarray(first_calls[variant]["log10_OmegaGW"], dtype=np.float64).tolist() for variant in handoff.VARIANTS},
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "supplied_n_t": params["n_t"],
        "derived_n_t": actual_nt,
        "root_relative_error_vs_reference_z12": output["root_relative_difference_vs_reference_z12_map"],
        "variants": {key: {
            "max_spectrum_error_dex": value["max_abs_spectrum_error_vs_reference_z12_dex"],
            "p50_runtime_s": value["p50_runtime_s"],
            "p95_runtime_s": value["p95_runtime_s"],
            "selected": value["adaptive_selected_count"],
        } for key, value in summaries.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
