"""Standalone grouped-SoA propagation spike.

The production kernel parallelizes one frequency mode per ``prange`` item.
This diagnostic groups nearby horizon starts into narrow buckets and advances
the Cartesian transfer states in lockstep with structure-of-arrays storage.
Each lane still calls the exact production ``_phase_segment`` and assembly
helpers; only execution layout changes.  The wrapper reports bucket work
overhead because inactive lanes are the central SIMD trade-off.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
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

from scripts.benchmark_phase_recurrence import CASES, _make_args, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

DEFAULT_BUCKET_WIDTH = 32


def _bucket_layout(j0s, bucket_width):
    j0s = np.asarray(j0s, dtype=np.int64)
    if bucket_width < 1:
        raise ValueError("bucket_width must be positive")
    order = np.argsort(j0s, kind="stable").astype(np.int64)
    keys = j0s[order] // int(bucket_width)
    starts = [0]
    for index in range(1, len(order)):
        if keys[index] != keys[index - 1]:
            starts.append(index)
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.empty(starts.size, dtype=np.int64)
    ends[:-1] = starts[1:]
    ends[-1] = len(order)
    return order, starts, ends


def _work_overhead(j0s, z0s, phi_grid, bucket_width, z_tail):
    order, starts, ends = _bucket_layout(j0s, bucket_width)
    j0s = np.asarray(j0s, dtype=np.int64)
    z0s = np.asarray(z0s, dtype=np.float64)
    phi_grid = np.asarray(phi_grid, dtype=np.float64)
    tail_end = np.empty(len(j0s), dtype=np.int64)
    for mode, j0 in enumerate(j0s):
        hit = np.flatnonzero(z0s[mode] + phi_grid - phi_grid[j0] >= z_tail)
        tail_end[mode] = int(hit[0]) if hit.size else len(phi_grid) - 1
    actual = np.maximum(tail_end - j0s, 0).astype(np.int64)
    actual_work = int(np.sum(actual))
    lockstep_work = 0
    for begin, end in zip(starts, ends):
        members = order[begin:end]
        lockstep_work += len(members) * int(np.max(tail_end[members]) - np.min(j0s[members]))
    return {
        "bucket_width": int(bucket_width),
        "bucket_count": int(starts.size),
        "max_bucket_lanes": int(max(ends - starts)),
        "actual_work": actual_work,
        "lockstep_work": int(lockstep_work),
        "work_overhead": float(lockstep_work / max(actual_work, 1)),
        "j0_span": int(np.max(j0s) - np.min(j0s)),
        "tail_end_span": int(np.max(tail_end) - np.min(tail_end)),
    }


@njit(parallel=True, cache=True)
def _solve_grouped_soa(
    Nv,
    Phi_grid,
    Phi_mid,
    S2,
    S2inv,
    j0s,
    z0s,
    P_t,
    ev_minus,
    fp_minus,
    fp_freq,
    assemble,
    n_coarse,
    col_step,
    h,
    z_tail,
    Ogw,
    Oj,
    Opgw,
    h_arr=None,
    Sv=None,
    phase_max=0.0,
    handoff_eps=None,
    kink_index=-1,
    kink_fraction=0.0,
    phi_re=0.0,
    mode_order=None,
    bucket_starts=None,
    bucket_ends=None,
):
    """Run exact Cartesian lanes in grouped lockstep; diagnostic only."""
    nv = len(Nv)
    for bucket in prange(len(bucket_starts)):
        begin = bucket_starts[bucket]
        end = bucket_ends[bucket]
        first_mode = mode_order[begin]
        min_k = j0s[first_mode]
        for lane in range(begin, end):
            mode = mode_order[lane]
            if j0s[mode] < min_k:
                min_k = j0s[mode]
        count = end - begin
        x_state = np.empty(count, dtype=np.float64)
        y_state = np.empty(count, dtype=np.float64)
        z_state = np.empty(count, dtype=np.float64)
        last_z_state = np.empty(count, dtype=np.float64)
        lx_state = np.empty(count, dtype=np.float64)
        ly_state = np.empty(count, dtype=np.float64)
        kend_state = np.empty(count, dtype=np.int64)
        active = np.zeros(count, dtype=np.bool_)
        done = np.zeros(count, dtype=np.bool_)
        needs_tail = np.zeros(count, dtype=np.bool_)

        for local in range(count):
            mode = mode_order[begin + local]
            j0 = j0s[mode]
            z0 = z0s[mode]
            xh = 0.0
            yh = math.exp(z0) * S2inv[j0]
            zz = z0
            x_state[local] = xh
            y_state[local] = yh
            z_state[local] = zz
            last_z_state[local] = zz
            lx_state[local] = 0.0
            ly_state[local] = yh
            kend_state[local] = j0
            if j0 % col_step == 0:
                if assemble:
                    slot = j0 // col_step
                    if slot >= n_coarse - 1:
                        slot = n_coarse - 1
                    FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[j0], xh, yh, zz, P_t[mode])
            elif j0 == nv - 1:
                FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1, S2[j0], xh, yh, zz, P_t[mode])
            if j0 >= nv - 1:
                done[local] = True
            elif zz >= z_tail:
                done[local] = True
                needs_tail[local] = True
                kend_state[local] = j0
            else:
                active[local] = True

        for step_k in range(min_k, nv - 1):
            any_active = False
            for local in range(count):
                mode = mode_order[begin + local]
                if done[local] or j0s[mode] > step_k:
                    continue
                any_active = True
                xh = x_state[local]
                yh = y_state[local]
                zz = z_state[local]
                h_step = h_arr[step_k] if h_arr is not None else h
                z_node = zz
                z_mid_step = z0s[mode] + Phi_mid[step_k] - Phi_grid[j0s[mode]]
                if step_k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0s[mode] + phi_re - Phi_grid[j0s[mode]]
                    z_end = z0s[mode] + Phi_grid[step_k + 1] - Phi_grid[j0s[mode]]
                    h_left = h_step * kink_fraction
                    h_right = h_step - h_left
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(xh, yh, z_break, z_end, h_right, phase_max)
                else:
                    z_end = 2.0 * z_mid_step - z_node
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_end, h_step, phase_max)
                step_next = step_k + 1
                zz = z0s[mode] + Phi_grid[step_next] - Phi_grid[j0s[mode]]
                if step_next % col_step == 0:
                    if assemble:
                        slot = step_next // col_step
                        if slot >= n_coarse - 1:
                            slot = n_coarse - 1
                        FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[step_next], xh, yh, zz, P_t[mode])
                elif step_next == nv - 1:
                    FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1, S2[step_next], xh, yh, zz, P_t[mode])
                if zz < z_tail:
                    lx_state[local] = xh
                    ly_state[local] = yh
                    last_z_state[local] = zz
                else:
                    done[local] = True
                    needs_tail[local] = True
                    kend_state[local] = step_k
                x_state[local] = xh
                y_state[local] = yh
                z_state[local] = zz
            if not any_active:
                break

        for local in range(count):
            if not needs_tail[local]:
                continue
            mode = mode_order[begin + local]
            kend = kend_state[local]
            if kend >= nv - 1:
                continue
            s2k = S2[kend]
            e_z = math.exp(-last_z_state[local])
            if Sv is not None:
                gamma = FS._tail_match_gamma(Sv[kend])
                amp2 = (lx_state[local] * lx_state[local]
                        + ly_state[local] * ly_state[local] * (1.0 + gamma * gamma * e_z * e_z)
                        + 2.0 * gamma * lx_state[local] * ly_state[local] * e_z)
            else:
                amp2 = lx_state[local] * lx_state[local] + ly_state[local] * ly_state[local]
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z_state[local])
            if handoff_eps is not None:
                if Sv is not None:
                    handoff_eps[mode] = abs(1.5 * Sv[kend] - 1.0) * e_z
                else:
                    handoff_eps[mode] = 0.0
            slot_start = kend // col_step
            if slot_start >= n_coarse - 1:
                slot_start = n_coarse - 1
            while slot_start < n_coarse:
                kk2 = col_step * slot_start
                if slot_start == n_coarse - 1:
                    kk2 = nv - 1
                if kk2 > kend:
                    break
                slot_start += 1
            for slot in range(slot_start, n_coarse):
                if not (assemble or slot == n_coarse - 1):
                    continue
                kk2 = col_step * slot
                if slot == n_coarse - 1:
                    kk2 = nv - 1
                FS.assemble_tail(
                    Ogw, Oj, Opgw, mode, slot, kk2, coeff, eNz, fp_freq[mode],
                    P_t[mode], ev_minus, fp_minus)


def solve_kernel_grouped_soa(*args):
    """Production-signature adapter that supplies a cached bucket layout."""
    j0s = np.asarray(args[5], dtype=np.int64)
    bucket_width = getattr(solve_kernel_grouped_soa, "bucket_width", DEFAULT_BUCKET_WIDTH)
    layout = _bucket_layout(j0s, bucket_width)
    return _solve_grouped_soa(*args, mode_order=layout[0], bucket_starts=layout[1], bucket_ends=layout[2])


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _first_divergence(base, candidate):
    for name, index in (("Ogw", 16), ("Oj", 17), ("Opgw", 18), ("handoff_eps", 22)):
        lhs = np.ascontiguousarray(base[index], dtype=np.float64)
        rhs = np.ascontiguousarray(candidate[index], dtype=np.float64)
        if not np.array_equal(lhs, rhs):
            where = np.argwhere(lhs != rhs)
            pos = tuple(int(value) for value in where[0])
            return {"variable": name, "index": pos, "baseline": float(lhs[pos]), "candidate": float(rhs[pos])}
    return None


def _run_kernel(name, repeats=30, threads=2, bucket_width=DEFAULT_BUCKET_WIDTH):
    solve_kernel_grouped_soa.bucket_width = int(bucket_width)
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    overhead = _work_overhead(common[5], common[6], common[1], bucket_width, common[15])
    FS.solve_kernel(*baseline)
    solve_kernel_grouped_soa(*candidate)
    first_divergence = _first_divergence(baseline, candidate)
    baseline_times = []
    candidate_times = []
    for _ in range(repeats):
        base_args = _make_args(common)
        candidate_args = _make_args(common)
        started = time.perf_counter()
        FS.solve_kernel(*base_args)
        baseline_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        solve_kernel_grouped_soa(*candidate_args)
        candidate_times.append(time.perf_counter() - started)
    return {
        "case": name,
        "threads": threads,
        "repeats": repeats,
        "bucket_width": int(bucket_width),
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "bitwise_equal": first_divergence is None,
        "first_divergence": first_divergence,
        "digest_Ogw_equal": _digest(baseline[16]) == _digest(candidate[16]),
        "digest_Oj_equal": _digest(baseline[17]) == _digest(candidate[17]),
        "digest_Opgw_equal": _digest(baseline[18]) == _digest(candidate[18]),
        "digest_handoff_eps_equal": _digest(baseline[22]) == _digest(candidate[22]),
        **overhead,
    }


def _run_outer(name, repeats=25, threads=2, bucket_width=DEFAULT_BUCKET_WIDTH):
    solve_kernel_grouped_soa.bucket_width = int(bucket_width)
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    original = FS.solve_kernel
    baseline_times = []
    candidate_times = []
    baseline_model = candidate_model = None
    try:
        for _ in range(repeats):
            baseline_model = LCDM_SG(**CASES[name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(baseline_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
            baseline_times.append(time.perf_counter() - started)
            FS.solve_kernel = solve_kernel_grouped_soa
            candidate_model = LCDM_SG(**CASES[name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(candidate_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
            candidate_times.append(time.perf_counter() - started)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original
    return {
        "case": name,
        "threads": threads,
        "repeats": repeats,
        "bucket_width": int(bucket_width),
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "digest_spectrum_equal": _digest(baseline_model.log10OmegaGW) == _digest(candidate_model.log10OmegaGW),
        "digest_DN_gw_equal": _digest(baseline_model.DN_gw) == _digest(candidate_model.DN_gw),
        "digest_g2_equal": _digest(baseline_model.g2) == _digest(candidate_model.g2),
        "digest_w2_equal": _digest(baseline_model.w2) == _digest(candidate_model.w2),
        "failure_equal": getattr(baseline_model, "fast_failure_reason", None) == getattr(candidate_model, "fast_failure_reason", None),
        "converged_equal": getattr(baseline_model, "SGWB_converge", False) == getattr(candidate_model, "SGWB_converge", False),
        "spectrum_dex_max": float(np.max(np.abs(np.asarray(candidate_model.log10OmegaGW) - np.asarray(baseline_model.log10OmegaGW)))),
        "DN_gw_relative": abs(float(candidate_model.DN_gw[-1]) - float(baseline_model.DN_gw[-1])) / max(abs(float(baseline_model.DN_gw[-1])), 1e-300),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--bucket-width", type=int, default=DEFAULT_BUCKET_WIDTH)
    parser.add_argument("--outer", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    os.environ.setdefault("FAST_THREADS", str(args.threads))
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)
    cases = args.case or list(CASES)
    runner = _run_outer if args.outer else _run_kernel
    payload = {
        "candidate": "grouped_soa_cartesian_transfer",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "bucket_width": int(args.bucket_width),
        "scope": "full_outer" if args.outer else "kernel",
        "cases": [runner(case, args.repeats, args.threads, args.bucket_width) for case in cases],
    }
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
