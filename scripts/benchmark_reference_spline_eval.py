"""Paired Reference benchmark for compiled scalar cubic-spline evaluation.

The SciPy variant calls the exact same ``CubicSpline`` objects as the previous
implementation. The candidate evaluates their stored piecewise-cubic
coefficients with a Numba-compiled scalar Horner evaluator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)
DN_EFF = 0.133732083004284
CASES = {"no_tail": -14.6013, "tail": -12.6723}


def _percentiles(values):
    return {"p50": float(np.percentile(values, 50)),
            "p95": float(np.percentile(values, 95))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=int, default=3)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--reference-archive", type=Path,
                        default=ROOT / "docs" /
                        "res4_sigma_exact_reference_map_20261104.npz")
    parser.add_argument("--archive-indices", type=int, nargs="+",
                        default=[0, 65, 130, 259])
    parser.add_argument("--archive-rtol", type=float, default=1e-11)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" /
                        "reference_spline_eval_benchmark_20261104.json")
    args = parser.parse_args(argv)
    if args.pairs < 1:
        parser.error("--pairs must be positive")

    splines = [gp.spl_rho, gp.spl_rhop, REF._FD_RHO_CS, REF._FD_P_CS]
    scipy_splines = {(spline.x.ctypes.data, spline.c.ctypes.data): spline
                     for spline in splines}

    def scipy_eval(x, c, point):
        return float(scipy_splines[(x.ctypes.data, c.ctypes.data)](point))

    compiled_eval = REF._cubic_spline_eval
    compile_start = time.perf_counter()
    compiled_eval(gp.spl_rho.x, gp.spl_rho.c, float(gp.spl_rho.x[0]))
    compile_seconds = time.perf_counter() - compile_start
    reference_path = ROOT / "stiffgwpy_fast" / "reference.py"
    script_digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    reference_digest = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    rows = []
    summaries = {}
    for case, freq in CASES.items():
        case_rows = []
        for pair in range(args.pairs):
            order = (("scipy", "numba") if pair % 2 == 0
                     else ("numba", "scipy"))
            pair_outputs = {}
            for method in order:
                REF._cubic_spline_eval = (
                    scipy_eval if method == "scipy" else compiled_eval)
                start = time.perf_counter()
                result = REF.solve_reference_mode(
                    model, freq, DN_EFF, z_tail=args.z_tail,
                    rtol=args.rtol, backend="cartesian")
                runtime = time.perf_counter() - start
                output = [float(result[key]) for key in
                          ("Ogw_today", "Oj_today", "Opgw_today")]
                pair_outputs[method] = output
                row = {
                    "case": case,
                    "pair": pair,
                    "method": method,
                    "runtime_s": runtime,
                    "nfev": int(result["n_steps"]),
                    "used_tail": bool(result["used_tail"]),
                    "output": output,
                }
                case_rows.append(row)
            scipy_output = pair_outputs["scipy"]
            numba_output = pair_outputs["numba"]
            abs_diff = [abs(a - b) for a, b in zip(scipy_output, numba_output)]
            rel_diff = [abs(a - b) / max(abs(a), abs(b), 1e-300)
                        for a, b in zip(scipy_output, numba_output)]
            case_rows[-1]["paired_max_abs_output_diff"] = max(abs_diff)
            case_rows[-1]["paired_max_rel_output_diff"] = max(rel_diff)
            if any(not math.isfinite(value) for value in rel_diff):
                raise FloatingPointError("non-finite paired spline output")
        dense_rows = [row for row in case_rows if row["method"] == "scipy"]
        compiled_rows = [row for row in case_rows if row["method"] == "numba"]
        dense_times = [row["runtime_s"] for row in dense_rows]
        compiled_times = [row["runtime_s"] for row in compiled_rows]
        ratios = [old / new for old, new in zip(dense_times, compiled_times)]
        expected_tail = case == "tail"
        if any(row["used_tail"] is not expected_tail for row in case_rows):
            raise RuntimeError("%s no longer exercises its expected path" % case)
        summaries[case] = {
            "freq_log10_hz": freq,
            "used_tail": expected_tail,
            "scipy_runtime_s": _percentiles(dense_times),
            "numba_runtime_s": _percentiles(compiled_times),
            "median_speedup": statistics.median(dense_times)
            / statistics.median(compiled_times),
            "paired_speedup_p50": float(np.percentile(ratios, 50)),
            "paired_speedup_p95": float(np.percentile(ratios, 95)),
            "scipy_nfev": [row["nfev"] for row in dense_rows],
            "numba_nfev": [row["nfev"] for row in compiled_rows],
            "max_pair_relative_output_diff": max(
                row["paired_max_rel_output_diff"] for row in case_rows
                if "paired_max_rel_output_diff" in row),
        }
        rows.extend(case_rows)
    REF._cubic_spline_eval = compiled_eval
    with np.load(args.reference_archive) as archive:
        reference_freqs = np.asarray(archive["freqs"], dtype=float)
        reference_log = np.asarray(archive["reference_log"], dtype=float)
    if reference_freqs.shape != reference_log.shape:
        raise ValueError("Reference archive frequency and spectrum shapes differ")
    if any(index < 0 or index >= len(reference_freqs)
           for index in args.archive_indices):
        parser.error("archive indices must be inside the Reference grid")
    archive_rows = []
    for index in args.archive_indices:
        freq = float(reference_freqs[index])
        start = time.perf_counter()
        result = REF.solve_reference_mode(
            model, freq, DN_EFF, z_tail=args.z_tail,
            rtol=args.archive_rtol, backend="cartesian")
        runtime = time.perf_counter() - start
        omega = float(result["Ogw_today"] - result["Oj_today"])
        if not math.isfinite(omega) or omega <= 0.0:
            raise FloatingPointError(
                "nonpositive/nonfinite candidate Omega at archive index %d"
                % index)
        log_omega = math.log10(omega)
        archive_rows.append({
            "index": int(index),
            "freq_log10_hz": freq,
            "runtime_s": runtime,
            "nfev": int(result["n_steps"]),
            "used_tail": bool(result["used_tail"]),
            "candidate_log10_omega": log_omega,
            "reference_log10_omega": float(reference_log[index]),
            "abs_dex_error": abs(log_omega - float(reference_log[index])),
        })
    payload = {
        "schema_version": 1,
        "experiment": "reference_numba_scalar_cubic_spline_evaluation",
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
        "reference_archive": str(args.reference_archive.resolve()),
        "archive_rtol": args.archive_rtol,
        "archive_indices": args.archive_indices,
        "jit_warmup_compile_seconds": compile_seconds,
        "resources": telemetry(workers=1, threads=1),
        "summaries": summaries,
        "samples": rows,
        "archive_check": {
            "max_abs_dex_error": max(row["abs_dex_error"]
                                      for row in archive_rows),
            "rows": archive_rows,
        },
        "method": "alternating SciPy CubicSpline scalar call versus Numba Horner evaluation of the same stored coefficients",
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
