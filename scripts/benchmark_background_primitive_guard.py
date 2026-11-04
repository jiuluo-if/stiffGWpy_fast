"""Paired normal-case overhead check for the Fast primitive finite-value guard.

The disabled arm replaces only the private guard with a no-op; equations,
backgrounds, grids, and all other solver code are identical.
中文：配对测量背景原函数有限值检查的正常点开销；关闭组只将私有检查器替换为 no-op。
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
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    "default": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "lowT": dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    "highT": dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "stiff": dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    "low_r": dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    "high_kappa": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def run(params, guard_enabled, config):
    original_guard = FS._first_nonfinite_primitive
    if not guard_enabled:
        FS._first_nonfinite_primitive = lambda _arrays: None
    try:
        model = LCDM_SG(**params)
        started = time.perf_counter_ns()
        result = FS.SGWB_iter_fast(
            model, config=config, tol=1e-6, freq_res=1.0,
            physical_guard=False,
        )
        elapsed_ms = (time.perf_counter_ns() - started) / 1e6
    finally:
        FS._first_nonfinite_primitive = original_guard
    if result is not model:
        raise RuntimeError(f"Fast failed on normal benchmark point: {model.fast_failure_reason}")
    spectrum = np.asarray(model.log10OmegaGW, dtype=float).copy()
    digest = hashlib.sha256(np.ascontiguousarray(spectrum).tobytes()).hexdigest()
    return elapsed_ms, spectrum, float(model.DN_gw[-1]), digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=501)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--affinity", default=None,
                        help="optional comma-separated CPU IDs")
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "docs" / "background_primitive_guard_benchmark_20261104.json",
    )
    args = parser.parse_args()
    if args.pairs < 1 or args.threads < 1:
        parser.error("pairs and threads must be positive")
    process = psutil.Process()
    available = list(process.cpu_affinity())
    selected = ([int(value) for value in args.affinity.split(",")]
                if args.affinity else available[:args.threads])
    if not selected or any(cpu not in available for cpu in selected):
        parser.error(f"affinity must be a nonempty subset of {available}")
    process.cpu_affinity(selected)
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)
    rows = []
    for name, params in CASES.items():
        config = FS.resolve_config("fast", z_tail=5.0, threads=args.threads)
        for enabled in (False, True):
            for _ in range(2):
                run(params, enabled, config)
        times = {False: [], True: []}
        digests = {False: [], True: []}
        outputs = {}
        for pair in range(args.pairs):
            order = (False, True) if pair % 2 == 0 else (True, False)
            for enabled in order:
                elapsed, spectrum, dn_gw, digest = run(params, enabled, config)
                times[enabled].append(elapsed)
                digests[enabled].append(digest)
                outputs[enabled] = (spectrum, dn_gw)
        baseline50, guard50 = (statistics.median(times[key]) for key in (False, True))
        baseline95, guard95 = (float(np.percentile(times[key], 95)) for key in (False, True))
        base_spectrum, base_dn = outputs[False]
        guard_spectrum, guard_dn = outputs[True]
        rows.append({
            "case": name,
            "baseline_p50_p95_ms": [baseline50, baseline95],
            "guard_p50_p95_ms": [guard50, guard95],
            "guard_over_baseline_p50_p95": [guard50 / baseline50, guard95 / baseline95],
            "dn_gw_abs_delta": abs(guard_dn - base_dn),
            "spectrum_max_abs_delta_dex": float(np.max(np.abs(guard_spectrum - base_spectrum))),
            "baseline_deterministic": len(set(digests[False])) == 1,
            "guard_deterministic": len(set(digests[True])) == 1,
            "baseline_elapsed_ms": times[False],
            "guard_elapsed_ms": times[True],
        })
    record = {
        "source_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "candidate_diff": "pre-kernel finite-value check for Phi/S2 primitives",
        "pairs": args.pairs,
        "threads": args.threads,
        "affinity": selected,
        "resource_after_mib": process.memory_info().rss / (1024.0 ** 2),
        "rows": rows,
    }
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps({
        "pairs": args.pairs,
        "rows": [{key: row[key] for key in (
            "case", "baseline_p50_p95_ms", "guard_p50_p95_ms",
            "guard_over_baseline_p50_p95", "dn_gw_abs_delta",
            "spectrum_max_abs_delta_dex")}
            for row in rows],
        "artifact": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()
