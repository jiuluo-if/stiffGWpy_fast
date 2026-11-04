"""Measure DOP853 dense-output storage overhead on fixed reference modes.

The baseline stores an interpolant after every accepted step. The candidate
keeps only the event state returned by ``y_events`` and disables the global
dense-output trajectory. The physics, tolerances, background, and event are
identical.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)
DN_EFF = 0.133732083004284
CASES = {"no_tail": -14.6013, "tail": -12.6723}


class _PeakRss:
    def __init__(self):
        try:
            import psutil
        except ImportError:
            self.process = None
            self.peak = None
            return
        self.process = psutil.Process()
        self.peak = self.process.memory_info().rss
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._sample, daemon=True)

    def _sample(self):
        while not self.stop.wait(0.02):
            self.peak = max(self.peak, self.process.memory_info().rss)

    def __enter__(self):
        if self.process is not None:
            self.thread.start()
        return self

    def __exit__(self, *_):
        if self.process is not None:
            self.stop.set()
            self.thread.join()
            self.peak = max(self.peak, self.process.memory_info().rss)

    def delta(self):
        if self.process is None:
            return None
        return int(self.peak - self.process.memory_info().rss)


def _run_variant(model, freq, z_tail, rtol, dense_output):
    original_solve_ivp = REF.solve_ivp

    def configured_solve_ivp(*args, **kwargs):
        kwargs["dense_output"] = dense_output
        return original_solve_ivp(*args, **kwargs)

    REF.solve_ivp = configured_solve_ivp
    try:
        start = time.perf_counter()
        result = REF.solve_reference_mode(
            model, freq, DN_EFF, z_tail=z_tail, rtol=rtol,
            backend="cartesian")
        runtime = time.perf_counter() - start
    finally:
        REF.solve_ivp = original_solve_ivp
    return result, runtime


def _percentiles(values):
    return {"p50": float(np.percentile(values, 50)),
            "p95": float(np.percentile(values, 95))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=3)
    parser.add_argument("--rtol", type=float, default=1e-11)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" /
                        "reference_dense_output_benchmark_20261104.json")
    args = parser.parse_args(argv)
    if args.pairs < 1:
        parser.error("--pairs must be positive")

    reference_path = ROOT / "stiffgwpy_fast" / "reference.py"
    script_digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    reference_digest = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    records = []
    summaries = {}
    for case, freq in CASES.items():
        model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
        rows = []
        for pair in range(args.pairs):
            order = ((True, False) if pair % 2 == 0 else (False, True))
            pair_outputs = {}
            for dense_output in order:
                variant = "dense" if dense_output else "nodense"
                with _PeakRss() as memory:
                    result, runtime = _run_variant(
                        model, freq, args.z_tail, args.rtol, dense_output)
                row = {
                    "pair": pair,
                    "variant": variant,
                    "runtime_s": runtime,
                    "nfev": int(result["n_steps"]),
                    "used_tail": bool(result["used_tail"]),
                    "output": [float(result[key]) for key in
                               ("Ogw_today", "Oj_today", "Opgw_today")],
                    "rss_peak_delta_bytes": memory.delta(),
                }
                rows.append(row)
                pair_outputs[variant] = row["output"]
            if not np.array_equal(pair_outputs["dense"],
                                  pair_outputs["nodense"]):
                raise RuntimeError(
                    "%s pair %d changed Cartesian outputs with dense_output"
                    % (case, pair))

        dense_rows = [row for row in rows if row["variant"] == "dense"]
        nodense_rows = [row for row in rows if row["variant"] == "nodense"]
        expected_tail = case == "tail"
        if any(row["used_tail"] is not expected_tail for row in rows):
            raise RuntimeError("%s no longer exercises its intended tail path" % case)
        dense_times = [row["runtime_s"] for row in dense_rows]
        nodense_times = [row["runtime_s"] for row in nodense_rows]
        summaries[case] = {
            "freq_log10_hz": freq,
            "used_tail": expected_tail,
            "dense_runtime_s": _percentiles(dense_times),
            "nodense_runtime_s": _percentiles(nodense_times),
            "median_speedup": statistics.median(dense_times)
            / statistics.median(nodense_times),
            "dense_nfev": dense_rows[-1]["nfev"],
            "nodense_nfev": nodense_rows[-1]["nfev"],
            "nfev_reduction": dense_rows[-1]["nfev"]
            - nodense_rows[-1]["nfev"],
            "dense_rss_peak_delta_bytes": max(
                (row["rss_peak_delta_bytes"] for row in dense_rows
                 if row["rss_peak_delta_bytes"] is not None), default=None),
            "nodense_rss_peak_delta_bytes": max(
                (row["rss_peak_delta_bytes"] for row in nodense_rows
                 if row["rss_peak_delta_bytes"] is not None), default=None),
            "bitwise_output_match": True,
        }
        records.extend({"case": case, **row} for row in rows)

    payload = {
        "schema_version": 1,
        "experiment": "cartesian_dop853_dense_output_storage_overhead",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "source_script_sha256": script_digest,
        "reference_module_sha256": reference_digest,
        "parameters": {**PARAMS, "DN_eff": DN_EFF},
        "z_tail": args.z_tail,
        "rtol": args.rtol,
        "atol": [1e-12, 1e-22, 1e-22],
        "pairs_per_case": args.pairs,
        "resources": telemetry(workers=1, threads=1),
        "summaries": summaries,
        "samples": records,
        "method": "same process alternating DOP853 dense_output True versus False; y_events supplies event state",
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
