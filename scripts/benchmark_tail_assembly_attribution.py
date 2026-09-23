"""Amdahl screen for common tail-assembly factor reuse.

The production tail path repeatedly forms ``ev_minus[kk]`` and
``fp_minus[kk]`` for each mode and output slot.  This diagnostic measures the
real mode/slot workload and a candidate that loads their precomputed product.
It is not a production twin: it must first show enough kernel headroom to
justify reproducing the complete solve kernel and then pass its exactness
gates.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def _tail_factors(ev_minus, fp_minus):
    return np.asarray(ev_minus, dtype=np.float64) * np.asarray(fp_minus, dtype=np.float64)


@njit(parallel=True, cache=False)
def _tail_baseline(Nv, S2, fp_freq, P_t, ev_minus, fp_minus,
                   coeff, eNz, col_step, Ogw, Oj, Opgw):
    for mode in prange(len(fp_freq)):
        for slot in range(Ogw.shape[1]):
            kk = col_step * slot
            if slot == Ogw.shape[1] - 1 or kk >= len(Nv):
                kk = len(Nv) - 1
            th = coeff[mode] * eNz[mode] * ev_minus[kk]
            xf = th * fp_freq[mode] * fp_minus[kk]
            oj = -th * th / 3.0 * P_t[mode]
            op = xf * xf / 36.0 * P_t[mode]
            Oj[mode, slot] = oj
            Opgw[mode, slot] = op
            Ogw[mode, slot] = 3.0 * op + oj


@njit(parallel=True, cache=False)
def _tail_factor_candidate(Nv, S2, fp_freq, P_t, tail_factor,
                           coeff, eNz, col_step, Ogw, Oj, Opgw):
    for mode in prange(len(fp_freq)):
        for slot in range(Ogw.shape[1]):
            kk = col_step * slot
            if slot == Ogw.shape[1] - 1 or kk >= len(Nv):
                kk = len(Nv) - 1
            th = coeff[mode] * eNz[mode] * tail_factor[kk]
            xf = th * fp_freq[mode]
            oj = -th * th / 3.0 * P_t[mode]
            op = xf * xf / 36.0 * P_t[mode]
            Oj[mode, slot] = oj
            Opgw[mode, slot] = op
            Ogw[mode, slot] = 3.0 * op + oj


def _inputs(common):
    Nv, _, _, S2, _, _, _, P_t, ev_minus, fp_minus, fp_freq, _, n_coarse, col_step, *_ = common
    modes = len(fp_freq)
    coeff = np.linspace(0.5, 1.5, modes, dtype=np.float64)
    eNz = np.linspace(0.8, 1.2, modes, dtype=np.float64)
    shape = (modes, n_coarse)
    return (
        Nv, S2, fp_freq, P_t, ev_minus, fp_minus, coeff, eNz, col_step,
        np.zeros(shape), np.zeros(shape), np.zeros(shape),
    )


def _clear(outputs):
    for output in outputs:
        output.fill(0.0)


def run_case(name, repeats=30, threads=2):
    model, common = _prepared(name, threads)
    Nv, S2, fp_freq, P_t, ev_minus, fp_minus, coeff, eNz, col_step, *_ = _inputs(common)
    baseline_outputs = tuple(np.zeros_like(_inputs(common)[-3]) for _ in range(3))
    candidate_outputs = tuple(np.zeros_like(_inputs(common)[-3]) for _ in range(3))
    baseline_args = (Nv, S2, fp_freq, P_t, ev_minus, fp_minus, coeff, eNz,
                     col_step, *baseline_outputs)
    candidate_args = (Nv, S2, fp_freq, P_t, _tail_factors(ev_minus, fp_minus),
                      coeff, eNz, col_step, *candidate_outputs)
    _tail_baseline(*baseline_args)
    _tail_factor_candidate(*candidate_args)
    full_args = __import__("scripts.benchmark_phase_recurrence", fromlist=["_make_args"])._make_args(common)
    full_times = []
    base_times = []
    cand_times = []
    for _ in range(int(repeats)):
        _clear(baseline_outputs)
        _clear(candidate_outputs)
        start = time.perf_counter()
        FS.solve_kernel(*full_args)
        full_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _tail_baseline(*baseline_args)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _tail_factor_candidate(*candidate_args)
        cand_times.append(time.perf_counter() - start)
    bitwise_equal = all(np.array_equal(a, b) for a, b in zip(baseline_outputs, candidate_outputs))
    full_median = statistics.median(full_times)
    base_median = statistics.median(base_times)
    cand_median = statistics.median(cand_times)
    return {
        "case": name,
        "repeats": int(repeats),
        "full_kernel_median_ms": full_median * 1e3,
        "tail_baseline_median_ms": base_median * 1e3,
        "tail_candidate_median_ms": cand_median * 1e3,
        "tail_candidate_over_baseline": cand_median / base_median,
        "tail_share_of_full_kernel": base_median / full_median,
        "free_tail_upper_bound_fraction": base_median / full_median,
        "bitwise_equal": bitwise_equal,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="default,lowT,highT,stiff,high_kappa")
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--out", default="docs/tail_assembly_attribution_round44_20260923.json")
    args = parser.parse_args()
    rows = [run_case(name, args.repeats, args.threads) for name in args.cases.split(",")]
    payload = {
        "experiment": "round44_tail_assembly_factor_amdahl_screen",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "resources": {"numba_threads": args.threads, "blas_threads": 1, "workers": 1},
        "settings": vars(args),
        "records": rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
