"""Paired sample11 full-root check for continuous-sigma primitive dispatch.

The default variant follows the transition-refined API default. Passing
``--keep-node-primitive`` gives the previous explicit fast-node path.
中文：配对运行 sample11 完整自洽根，比较转折细化的连续 sigma 原函数与旧节点快速原函数。
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

from scripts.benchmark_sample11_reference_map import PARAMS  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _run(sigma_exact, config, reference_freqs, reference_log):
    model = LCDM_SG(**PARAMS, DN_eff=0.0)
    started = time.perf_counter_ns()
    result = FS.SGWB_iter_fast(
        model,
        config=config,
        freq_res=1.0,
        tol=1e-6,
        transition_refine=True,
        kink_split=True,
        sigma_exact=sigma_exact,
    )
    elapsed_ms = (time.perf_counter_ns() - started) / 1e6
    if result is not model:
        raise RuntimeError("sample11 Fast failed: %s" % model.fast_failure_reason)

    frequencies = np.asarray(model.f, dtype=float)
    order = np.argsort(frequencies)
    indexes = np.searchsorted(frequencies[order], reference_freqs)
    if (np.any(indexes >= frequencies.size)
            or not np.array_equal(frequencies[order][indexes], reference_freqs)):
        raise RuntimeError("Fast and Reference native frequency grids differ")
    error = np.abs(np.asarray(model.log10OmegaGW)[order][indexes] - reference_log)
    return {
        "elapsed_ms": elapsed_ms,
        "dn_eff": float(model.cosmo_param["DN_eff"]),
        "dn_gw_relative_residual": float(model.dn_converged_delta),
        "background_nodes": int(model.Nv.size),
        "spectrum_error_median_p95_max_dex": np.percentile(
            error, [50, 95, 100]
        ).tolist(),
        "modes_over_5e-4_dex": int(np.count_nonzero(error > 5e-4)),
        "spectrum_sha256": hashlib.sha256(
            np.ascontiguousarray(np.asarray(model.log10OmegaGW)[order][indexes])
            .tobytes()
        ).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--pairs", type=int, default=15)
    parser.add_argument("--affinity", default=None,
                        help="optional comma-separated CPU IDs")
    parser.add_argument(
        "--reference",
        type=Path,
        default=ROOT / "docs" / "sample11_h02_reference_map_20261104.npz",
    )
    parser.add_argument(
        "--reference-map",
        type=Path,
        default=ROOT / "docs" / "sample11_h02_reference_map_20261104.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs" / "sample11_sigma_exact_benchmark_20261104.json",
    )
    args = parser.parse_args()
    if args.threads < 1 or args.pairs < 1:
        parser.error("threads and pairs must be positive")

    process = psutil.Process()
    available = list(process.cpu_affinity())
    affinity = ([int(value) for value in args.affinity.split(",")]
                if args.affinity else available[:args.threads])
    if not affinity or any(cpu not in available for cpu in affinity):
        parser.error("affinity must be a nonempty subset of %s" % available)
    process.cpu_affinity(affinity)
    threads = len(affinity)

    reference = np.load(args.reference)
    reference_map = json.loads(args.reference_map.read_text(encoding="utf-8"))
    reference_freqs = np.asarray(reference["freqs"], dtype=float)
    reference_log = np.asarray(reference["reference_log"], dtype=float)
    order = np.argsort(reference_freqs)
    reference_freqs = reference_freqs[order]
    reference_log = reference_log[order]

    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    config = FS.resolve_config(
        "fast", freq_grid="grid_independent", h=0.02, z_tail=8.0,
        threads=threads, transition_refine_radius=0.0,
        transition_refine_factor=32,
    )

    cold = {
        "node_primitive": _run(False, config, reference_freqs, reference_log),
        "continuous_default": _run(None, config, reference_freqs, reference_log),
    }
    samples = {"node_primitive": [], "continuous_default": []}
    for pair in range(args.pairs):
        order = (["node_primitive", "continuous_default"] if pair % 2 == 0
                 else ["continuous_default", "node_primitive"])
        for name in order:
            sigma_exact = False if name == "node_primitive" else None
            samples[name].append(_run(
                sigma_exact, config, reference_freqs, reference_log
            ))

    variants = {}
    for name, rows in samples.items():
        variants[name] = {
            "warm_p50_ms": statistics.median(row["elapsed_ms"] for row in rows),
            "warm_p95_ms": float(np.percentile(
                [row["elapsed_ms"] for row in rows], 95
            )),
            **{key: value for key, value in rows[-1].items()
               if key != "elapsed_ms"},
            "raw_elapsed_ms": [row["elapsed_ms"] for row in rows],
            "deterministic": len({row["spectrum_sha256"] for row in rows}) == 1,
        }
    record = {
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "reference_artifact": str(args.reference.resolve()),
        "reference_map_record": str(args.reference_map.resolve()),
        "reference_input_dn_eff": reference_map.get("dn_eff"),
        "reference_dn_map": reference_map.get("reference_dn_map"),
        "threads": threads,
        "affinity": affinity,
        "pairs": args.pairs,
        "continuous_default_sigma_exact_argument": None,
        "cold_calls": cold,
        "variants": variants,
    }
    record["continuous_over_node_p50_ratio"] = (
        variants["continuous_default"]["warm_p50_ms"]
        / variants["node_primitive"]["warm_p50_ms"]
    )
    record["continuous_over_node_p95_ratio"] = (
        variants["continuous_default"]["warm_p95_ms"]
        / variants["node_primitive"]["warm_p95_ms"]
    )
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
