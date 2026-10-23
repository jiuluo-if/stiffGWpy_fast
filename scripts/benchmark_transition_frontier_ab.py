"""Compare formal fast with a matched transition-refinement candidate.

中文：在固定线程和 CPU 亲和性下配对测量正式 fast 与候选配置，并记录原始耗时和输出摘要。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, limit_affinity, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT_FIELDS = ("f", "log10OmegaGW", "DN_gw", "g2", "w2")
DEFAULT_CASES = ("edge_dnre_lo", "edge_dnre_hi")


def _digest(model: LCDM_SG) -> dict[str, str]:
    """Hash stable numerical outputs for one completed solver call.

    中文：为一次完成的求解计算稳定的数值输出摘要。
    """
    return {
        name: hashlib.sha256(
            np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64)).tobytes()
        ).hexdigest()
        for name in OUTPUT_FIELDS
    }


def _run_once(case: str, config, preset: dict, transition_refine: bool) -> dict:
    """Run one fresh model and return elapsed time and output identity.

    中文：对新模型求解一次，并返回耗时和输出摘要。
    """
    model = LCDM_SG(**CASES[case])
    start_ns = time.perf_counter_ns()
    FS.SGWB_iter_fast(
        model,
        kink_split=config.kink_split,
        transition_refine=transition_refine,
        config=config,
        tol=preset["tol"],
        freq_res=preset["freq_res"],
        freq_grid="goal",
        frequency_quadrature="pchip",
    )
    elapsed_ns = time.perf_counter_ns() - start_ns
    if not model.SGWB_converge or not hasattr(model, "DN_gw"):
        reason = getattr(model, "fast_failure_reason", None) or "not converged"
        raise RuntimeError(f"{case} solve failed: {reason}")
    return {
        "elapsed_ns": elapsed_ns,
        "digest": _digest(model),
        "DN_eff": float(model.cosmo_param["DN_eff"]),
        "DN_gw": float(np.asarray(model.DN_gw, dtype=np.float64).reshape(-1)[-1]),
        "n_freq": int(np.asarray(model.f).size),
    }


def _summarize_pairs(baseline: list[dict], candidate: list[dict], seed: int) -> dict:
    """Summarize paired timing ratios and their bootstrap interval.

    中文：汇总配对耗时比及其中位数自助法区间。
    """
    ratios = np.asarray(
        [cand["elapsed_ns"] / base["elapsed_ns"] for base, cand in zip(baseline, candidate)],
        dtype=np.float64,
    )
    rng = np.random.default_rng(seed)
    indexes = rng.integers(0, ratios.size, size=(20_000, ratios.size))
    bootstrap_medians = np.median(ratios[indexes], axis=1)
    return {
        "baseline_median_ms": statistics.median(row["elapsed_ns"] for row in baseline) / 1e6,
        "candidate_median_ms": statistics.median(row["elapsed_ns"] for row in candidate) / 1e6,
        "paired_median_ratio_candidate_over_baseline": float(np.median(ratios)),
        "paired_ratio_bootstrap_95_ci": [float(value) for value in np.percentile(bootstrap_medians, [2.5, 97.5])],
        "candidate_faster_pairs": int(np.count_nonzero(ratios < 1.0)),
        "pairs": int(ratios.size),
        "paired_ratios": ratios.tolist(),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", nargs="+", choices=sorted(CASES), default=list(DEFAULT_CASES))
    parser.add_argument("--pairs", type=int, default=30, help="measured alternating pairs per case")
    parser.add_argument("--warmups", type=int, default=2, help="calls per variant before measurement")
    parser.add_argument("--affinity-count", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h", type=float, default=0.0296875)
    parser.add_argument("--radius", type=float, default=0.015625)
    parser.add_argument("--factor", type=int, default=95)
    args = parser.parse_args(argv)
    if args.pairs < 1 or args.warmups < 0 or args.affinity_count < 1:
        parser.error("--pairs and --affinity-count must be positive; --warmups must be nonnegative")

    process = None
    try:
        import psutil

        process = psutil.Process()
    except ImportError:
        pass
    affinity = limit_affinity(process, args.affinity_count) if process is not None else []
    if not affinity:
        raise RuntimeError("psutil CPU affinity control is required for this fixed-resource benchmark")
    preset = FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    baseline_config = FS.resolve_config("fast")
    candidate_config = FS.resolve_config(
        "fast",
        h=args.h,
        kink_split=False,
        transition_refine_radius=args.radius,
        transition_refine_factor=args.factor,
    )
    variants = (
        ("formal_fast", baseline_config, False),
        ("transition_refine", candidate_config, True),
    )

    case_results = {}
    for case in args.cases:
        warmups = {}
        for name, config, refined in variants:
            rows = [_run_once(case, config, preset, refined) for _ in range(args.warmups)]
            digests = {json.dumps(row["digest"], sort_keys=True) for row in rows}
            if len(digests) > 1:
                raise RuntimeError(f"{case} {name} output changed during warmup")
            warmups[name] = rows

        samples = {name: [] for name, _, _ in variants}
        for pair_index in range(args.pairs):
            order = variants if pair_index % 2 == 0 else tuple(reversed(variants))
            pair = {}
            for name, config, refined in order:
                pair[name] = _run_once(case, config, preset, refined)
                samples[name].append(pair[name])

        stable_digests = {
            name: {json.dumps(row["digest"], sort_keys=True) for row in rows} for name, rows in samples.items()
        }
        if any(len(digests) != 1 for digests in stable_digests.values()):
            raise RuntimeError(f"{case} output changed across measured calls")
        if any(rows and rows[0]["digest"] != samples[name][0]["digest"] for name, rows in warmups.items()):
            raise RuntimeError(f"{case} output changed between warmup and measurement")
        if any(any(row["DN_eff"] != rows[0]["DN_eff"] for row in rows[1:]) for rows in samples.values()):
            raise RuntimeError(f"{case} DN_eff changed across measured calls")
        summary = _summarize_pairs(samples["formal_fast"], samples["transition_refine"], args.seed)
        case_results[case] = {
            "summary": summary,
            "output_digests": {name: rows[0]["digest"] for name, rows in samples.items()},
            "warmup_digests": {name: rows[0]["digest"] if rows else None for name, rows in warmups.items()},
            "samples": samples,
        }

    payload = {
        "schema_version": 1,
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "resources": telemetry(process, workers=1, threads=2),
        "affinity_selected": affinity,
        "warmups_per_variant": args.warmups,
        "pairs_per_case": args.pairs,
        "bootstrap": {"resamples": 20_000, "seed": args.seed, "statistic": "median of paired ratios"},
        "baseline": {"path": "kink-split", "configuration": baseline_config.__dict__},
        "candidate": {
            "path": "transition-refine",
            "h": candidate_config.h,
            "radius": candidate_config.transition_refine_radius,
            "factor": candidate_config.transition_refine_factor,
            "local_step": candidate_config.h / candidate_config.transition_refine_factor,
        },
        "cases": case_results,
        "timing_scope": "one converged solver call on a fresh model; excludes model construction, process startup, and JIT after warmup",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for case, result in case_results.items():
        row = result["summary"]
        print(
            f"{case}: baseline={row['baseline_median_ms']:.3f} ms, "
            f"candidate={row['candidate_median_ms']:.3f} ms, "
            f"paired ratio={row['paired_median_ratio_candidate_over_baseline']:.4f}, "
            f"95% CI={row['paired_ratio_bootstrap_95_ci']}"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
