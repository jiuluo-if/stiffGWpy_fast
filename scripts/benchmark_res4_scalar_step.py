"""Paired res4 accuracy/performance check for aligned scalar-step dispatch.

中文：配对比较默认 h=.005 路径与 h=.02 对齐转折细化路径，并记录同频参考误差、耗时、RSS 和吞吐。
"""

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
# 中文：先加入仓库根路径，再导入当前工作树内的求解器。
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0, DN_eff=0.0)
VARIANTS = {
    "uniform_h005": dict(h=0.005, transition=False, radius=0.0),
    "refined_h02_radius0": dict(h=0.02, transition=True, radius=0.0),
    "refined_h02_array_steps": dict(h=0.02, transition=True, radius=0.0,
                                    force_array=True),
}


def batch_case_rate(pair_count, variant_count, elapsed_seconds):
    """Count every executed variant, not merely paired baselines."""
    return pair_count * variant_count / elapsed_seconds


def run_variant(name, frequencies, reference_log, threads):
    option = VARIANTS[name]
    config = FS.resolve_config(
        "fast", freq_grid="grid_independent", h=option["h"],
        z_tail=8.0, threads=threads,
        transition_refine_radius=option["radius"], transition_refine_factor=32,
    )
    model = LCDM_SG(**PARAMS)
    start = time.perf_counter_ns()
    original_step_detector = FS._uniform_step_terminal_width
    if option.get("force_array", False):
        FS._uniform_step_terminal_width = lambda _Nv, _h: 0.0
    try:
        result = FS.SGWB_iter_fast(
            model, config=config, freq_res=4.0, tol=1e-6,
            transition_refine=option["transition"], kink_split=True,
        )
    finally:
        FS._uniform_step_terminal_width = original_step_detector
    elapsed_ms = (time.perf_counter_ns() - start) / 1e6
    if result is not model:
        raise RuntimeError(f"{name} failed: {model.fast_failure_reason}")
    f = np.asarray(model.f, dtype=float)
    y = np.asarray(model.log10OmegaGW, dtype=float)
    order = np.argsort(f)
    indexes = np.searchsorted(f[order], frequencies)
    if np.any(indexes >= f.size) or not np.array_equal(f[order][indexes], frequencies):
        raise RuntimeError("Fast and Reference native frequency grids differ")
    error = np.abs(y[order][indexes] - reference_log)
    return {
        "elapsed_ms": elapsed_ms,
        "rss_after_call_mib": psutil.Process().memory_info().rss / (1024.0 ** 2),
        "dn_eff": float(model.cosmo_param["DN_eff"]),
        "dn_gw_relative_residual": float(model.dn_converged_delta),
        "background_nodes": int(model.Nv.size),
        "spectrum_error_median_p95_max_dex": np.percentile(error, [50, 95, 100]).tolist(),
        "modes_over_5e-4_dex": int(np.count_nonzero(error > 5e-4)),
        "spectrum_sha256": hashlib.sha256(np.ascontiguousarray(y).tobytes()).hexdigest(),
        "log_spectrum": y[order][indexes].copy(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--pairs", type=int, default=15)
    parser.add_argument("--affinity", default=None,
                        help="optional comma-separated CPU IDs, for example 0,1,2,3,4,5,6,8")
    parser.add_argument(
        "--reference", type=Path,
        default=ROOT / "docs" / "res4_reference_secant_root_tail10_20261104.npz",
    )
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "docs" / "res4_scalar_step_benchmark_20261104.json",
    )
    args = parser.parse_args()
    if args.threads < 1 or args.pairs < 1:
        parser.error("threads and pairs must be positive")

    process = psutil.Process()
    available = list(process.cpu_affinity())
    selected = ([int(value) for value in args.affinity.split(",")]
                if args.affinity else available[:args.threads])
    if not selected or any(cpu not in available for cpu in selected):
        parser.error(f"affinity must be a nonempty subset of {available}")
    process.cpu_affinity(selected)
    ref = np.load(args.reference)
    frequencies = np.asarray(ref["freqs"], dtype=float)
    reference_log = np.asarray(ref["reference_log"], dtype=float)
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)

    first_raw = {name: run_variant(name, frequencies, reference_log, args.threads)
                 for name in VARIANTS}
    first_call = {
        name: {key: value for key, value in row.items() if key != "log_spectrum"}
        for name, row in first_raw.items()
    }
    samples = {name: [] for name in VARIANTS}
    batch_start = time.perf_counter()
    for pair in range(args.pairs):
        names = list(VARIANTS)
        if pair % 2:
            names.reverse()
        for name in names:
            samples[name].append(run_variant(name, frequencies, reference_log, args.threads))
    batch_seconds = time.perf_counter() - batch_start

    record = {
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "reference_artifact": str(args.reference.resolve()),
        "threads": args.threads,
        "affinity": selected,
        "pairs": args.pairs,
        "first_call": first_call,
        "batch_seconds": batch_seconds,
        "batch_cases_per_second": batch_case_rate(
            args.pairs, len(VARIANTS), batch_seconds
        ),
        "variants": {},
    }
    for name, rows in samples.items():
        record["variants"][name] = {
            "warm_p50_ms": statistics.median(row["elapsed_ms"] for row in rows),
            "warm_p95_ms": float(np.percentile([row["elapsed_ms"] for row in rows], 95)),
            "rss_after_call_max_mib": max(row["rss_after_call_mib"] for row in rows),
            **{key: value for key, value in rows[0].items()
               if key not in ("elapsed_ms", "rss_after_call_mib", "spectrum_sha256", "log_spectrum")},
            "deterministic": len({row["spectrum_sha256"] for row in rows}) == 1,
            "raw_elapsed_ms": [row["elapsed_ms"] for row in rows],
        }
    baseline = record["variants"]["uniform_h005"]
    candidate = record["variants"]["refined_h02_radius0"]
    record["candidate_over_default_p50_ratio"] = candidate["warm_p50_ms"] / baseline["warm_p50_ms"]
    record["candidate_over_default_p95_ratio"] = candidate["warm_p95_ms"] / baseline["warm_p95_ms"]
    scalar = first_raw["refined_h02_radius0"]
    array = first_raw["refined_h02_array_steps"]
    delta = np.abs(scalar["log_spectrum"] - array["log_spectrum"])
    record["scalar_vs_array_step_equivalence"] = {
        "dn_eff_abs_difference": abs(scalar["dn_eff"] - array["dn_eff"]),
        "spectrum_delta_median_p95_max_dex": np.percentile(delta, [50, 95, 100]).tolist(),
    }
    array_stats = record["variants"]["refined_h02_array_steps"]
    record["scalar_over_array_p50_ratio"] = candidate["warm_p50_ms"] / array_stats["warm_p50_ms"]
    record["scalar_over_array_p95_ratio"] = candidate["warm_p95_ms"] / array_stats["warm_p95_ms"]
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
