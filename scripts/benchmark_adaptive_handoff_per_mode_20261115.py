"""Paired warm A/B timing for uniform and z_today-gated tail handoff."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import psutil
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import underflow_adaptive_handoff_mask_screen_20261115 as screen  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_grid_20261115.json"
TAIL_Z8 = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_tail_depth_20261115.json"
REFERENCE_Z10 = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_per_mode_benchmark_20261115.json"
REPEATS = 12
THREADS = 1
VARIANTS = ("z5", "z8", "adaptive")
ORDERS = (VARIANTS, ("z8", "adaptive", "z5"), ("adaptive", "z5", "z8"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_variant(name: str, params: dict, freqs: np.ndarray, root: float) -> tuple[dict, int]:
    if name == "adaptive":
        result, _ = screen.run_fast(params, freqs, root, 5.0, adaptive_handoff=True)
    else:
        result, _ = screen.run_fast(params, freqs, root, float(name[1:]))
    return result, int(psutil.Process().memory_info().rss)


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    ref_data = json.loads(REFERENCE_Z10.read_text(encoding="utf-8"))
    case = next(row for row in source["cases"] if row["frequency_count"] == 937)
    ref_case = ref_data["reference"]
    freqs = np.asarray(case["frequency_log10_hz"], dtype=np.float64)
    params = source["parameters"]
    root = float(case["frozen_DN_eff"])
    reference_spectrum = np.asarray(ref_case["log10_OmegaGW"], dtype=np.float64)
    process = psutil.Process()
    rss_before = int(process.memory_info().rss)

    first_call_s = {}
    rss_after_first = {}
    for variant in VARIANTS:
        started = time.perf_counter()
        result, rss = run_variant(variant, params, freqs, root)
        first_call_s[variant] = time.perf_counter() - started
        rss_after_first[variant] = rss

    timings = {variant: [] for variant in VARIANTS}
    rss_samples = {variant: [] for variant in VARIANTS}
    spectra = {variant: [] for variant in VARIANTS}
    maps = {variant: [] for variant in VARIANTS}
    selected_counts = []
    for repeat in range(REPEATS):
        for variant in ORDERS[repeat % len(ORDERS)]:
            started = time.perf_counter()
            result, rss = run_variant(variant, params, freqs, root)
            timings[variant].append(time.perf_counter() - started)
            rss_samples[variant].append(rss)
            spectra[variant].append(result["log10_OmegaGW"])
            maps[variant].append(result["mapped_DN_gw"])
            if variant == "adaptive":
                selected_counts.append(result["adaptive_handoff_selected_count"])

    summaries = {}
    for variant in VARIANTS:
        spec = np.asarray(spectra[variant], dtype=np.float64)
        spread = np.max(np.abs(spec - spec[0]), axis=1)
        errors = np.abs(spec[0] - reference_spectrum)
        times = np.asarray(timings[variant], dtype=np.float64)
        rss_mib = np.asarray(rss_samples[variant], dtype=np.float64) / (1024.0 ** 2)
        summaries[variant] = {
            "p50_runtime_s": float(np.median(times)),
            "p95_runtime_s": float(np.quantile(times, 0.95)),
            "p50_rss_mib": float(np.median(rss_mib)),
            "p95_rss_mib": float(np.quantile(rss_mib, 0.95)),
            "max_repeat_spectrum_spread_dex": float(np.max(spread)),
            "max_abs_spectrum_error_vs_reference_z10_dex": float(np.max(errors)),
            "p95_abs_spectrum_error_vs_reference_z10_dex": float(np.quantile(errors, 0.95)),
            "mapped_DN_gw_p50": float(np.median(maps[variant])),
            "mapped_DN_gw_relative_to_frozen_root": float(np.median(maps[variant]) / root - 1.0),
            "paired_p50_speedup_vs_z5": None,
        }
    for variant in VARIANTS:
        summaries[variant]["paired_p50_speedup_vs_z5"] = (
            summaries["z5"]["p50_runtime_s"] / summaries[variant]["p50_runtime_s"]
        )
    output = {
        "schema_version": 1,
        "status": "complete_paired_warm_per_mode_handoff_benchmark",
        "scope": "Same minimum-subnormal sample, frozen root, exact 937-node grid and one-thread Fast profile; alternate 12 paired z5/z8/adaptive rounds after one per-variant warmup call.",
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "frozen_DN_eff": root,
        "adaptive_thresholds": {"z_today_below": 6.0, "defer_handoff_to": 8.0},
        "threads": THREADS,
        "pair_count_per_variant": REPEATS,
        "order_cycle": [list(order) for order in ORDERS],
        "first_call_after_process_start_s": first_call_s,
        "rss_after_first_call_mib": {k: v / (1024.0 ** 2) for k, v in rss_after_first.items()},
        "rss_before_mib": rss_before / (1024.0 ** 2),
        "adaptive_selected_mode_counts": sorted(set(selected_counts)),
        "variants": summaries,
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
                "screen_helper.py": "scripts/underflow_adaptive_handoff_mask_screen_20261115.py",
                "reference_z10_json": "docs/research/results/underflow_min_subnormal_reference_z10_20261115.json",
                "diagnostic_script": "scripts/benchmark_adaptive_handoff_per_mode_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "first_call_s": first_call_s,
        "variants": {k: {
            "p50_s": v["p50_runtime_s"], "p95_s": v["p95_runtime_s"],
            "speedup_vs_z5": v["paired_p50_speedup_vs_z5"],
            "max_spectrum_error": v["max_abs_spectrum_error_vs_reference_z10_dex"],
            "p95_spectrum_error": v["p95_abs_spectrum_error_vs_reference_z10_dex"],
        } for k, v in summaries.items()},
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
