"""Standalone P1 spike replacing per-step modulo/division assembly checks.

Only the assembly scheduling state changes.  Tensor propagation arithmetic,
phase map, kink handling, tail crossing, tail matching and assembly formulas
are the production versions.  The candidate records assembly slots so the
production schedule can be compared independently of floating-point outputs.
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

from scripts.benchmark_phase_recurrence import (  # noqa: E402
    CASES,
    _make_args,
    _prepared,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(parallel=True, cache=True)
def reference_assembly_trace(Nv, Phi_grid, j0s, z0s, assemble, n_coarse, col_step, z_tail, traces):
    """Record the slots selected by the production modulo/division schedule."""
    nv = len(Nv)
    for m in prange(len(j0s)):
        count = 0
        j0 = j0s[m]
        z0 = z0s[m]
        Phi0 = Phi_grid[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        if k % col_step == 0:
            if assemble:
                slot = k // col_step
                if slot >= n_coarse - 1:
                    slot = n_coarse - 1
                traces[m, count] = slot
                count += 1
        elif k == nv - 1:
            traces[m, count] = n_coarse - 1
            count += 1
        if zz < z_tail:
            while k < nv - 1 and zz < z_tail:
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0:
                    if assemble:
                        slot = k // col_step
                        if slot >= n_coarse - 1:
                            slot = n_coarse - 1
                        traces[m, count] = slot
                        count += 1
                elif k == nv - 1:
                    traces[m, count] = n_coarse - 1
                    count += 1
            if zz < z_tail:
                kend = nv - 1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        if kend < nv - 1:
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
                traces[m, count] = slot
                count += 1
        traces[m, count] = -2


@njit(parallel=True, cache=True)
def solve_kernel_counted_assembly(
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
    traces=None,
):
    """Production solve kernel with counted main-loop output state."""
    nv = len(Nv)
    for m in prange(len(j0s)):
        trace_count = 0
        j0 = j0s[m]
        z0 = z0s[m]
        Pt = P_t[m]
        fp_i = fp_freq[m]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0) * S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh = 0.0
        lyh = yh
        last_z = zz
        if k % col_step == 0:
            output_slot = k // col_step
            if output_slot >= n_coarse - 1:
                output_slot = n_coarse - 1
            next_output_k = k + col_step
            if assemble:
                if traces is not None:
                    traces[m, trace_count] = output_slot
                    trace_count += 1
                FS.assemble_main(Ogw, Oj, Opgw, m, output_slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            output_slot = n_coarse - 1
            next_output_k = nv
            if traces is not None:
                traces[m, trace_count] = output_slot
                trace_count += 1
            FS.assemble_main(Ogw, Oj, Opgw, m, output_slot, S2[k], xh, yh, zz, Pt)
        else:
            output_slot = k // col_step
            next_output_k = (output_slot + 1) * col_step
        if zz < z_tail:
            while k < nv - 1 and zz < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = zz
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    z_end = z0 + Phi_grid[k + 1] - Phi0
                    h_left = h_step * kink_fraction
                    h_right = h_step - h_left
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(xh, yh, z_break, z_end, h_right, phase_max)
                else:
                    z_end = 2.0 * z_mid_step - z_node
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_end, h_step, phase_max)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k == next_output_k:
                    output_slot += 1
                    if output_slot >= n_coarse - 1:
                        output_slot = n_coarse - 1
                    next_output_k += col_step
                    if assemble:
                        if traces is not None:
                            traces[m, trace_count] = output_slot
                            trace_count += 1
                        FS.assemble_main(Ogw, Oj, Opgw, m, output_slot, S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    if traces is not None:
                        traces[m, trace_count] = n_coarse - 1
                        trace_count += 1
                    FS.assemble_main(Ogw, Oj, Opgw, m, n_coarse - 1, S2[k], xh, yh, zz, Pt)
                if zz < z_tail:
                    lxh = xh
                    lyh = yh
                    last_z = zz
            if zz < z_tail:
                kend = nv - 1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        if kend < nv - 1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            if Sv is not None:
                gamma = FS._tail_match_gamma(Sv[kend])
                amp2 = lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z) + 2.0 * gamma * lxh * lyh * e_z
            else:
                amp2 = lxh * lxh + lyh * lyh
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                if Sv is not None:
                    handoff_eps[m] = abs(1.5 * Sv[kend] - 1.0) * e_z
                else:
                    handoff_eps[m] = 0.0
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
                if traces is not None:
                    traces[m, trace_count] = slot
                    trace_count += 1
                FS.assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus)
        if traces is not None:
            traces[m, trace_count] = -2


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _first_divergence(base, candidate):
    for name, index in (("Ogw", 16), ("Oj", 17), ("Opgw", 18), ("handoff_eps", 22)):
        lhs = np.ascontiguousarray(base[index], dtype=np.float64)
        rhs = np.ascontiguousarray(candidate[index], dtype=np.float64)
        if not np.array_equal(lhs, rhs):
            where = np.argwhere(lhs != rhs)
            pos = tuple(int(x) for x in where[0])
            return {"variable": name, "index": pos, "baseline": float(lhs[pos]), "candidate": float(rhs[pos])}
    return None


def _run_kernel(name, repeats=30, threads=2):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    trace_shape = (len(common[5]), common[12] + 2)
    expected_trace = np.full(trace_shape, -1, dtype=np.int64)
    candidate_trace = np.full(trace_shape, -1, dtype=np.int64)
    reference_assembly_trace(
        common[0], common[1], common[5], common[6], common[11], common[12], common[13], common[15], expected_trace
    )
    candidate = candidate + (candidate_trace,)
    FS.solve_kernel(*baseline)
    solve_kernel_counted_assembly(*candidate)
    first_divergence = _first_divergence(baseline, candidate)
    trace_equal = np.array_equal(expected_trace, candidate_trace)
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    base_dn = FS.integrate_frequency_pchip(model.f, base_obs)
    cand_dn = FS.integrate_frequency_pchip(model.f, cand_obs)
    baseline_times = []
    candidate_times = []
    for _ in range(repeats):
        args = _make_args(common)
        started = time.perf_counter()
        FS.solve_kernel(*args)
        baseline_times.append(time.perf_counter() - started)
        args = _make_args(common) + (np.full(trace_shape, -1, dtype=np.int64),)
        started = time.perf_counter()
        solve_kernel_counted_assembly(*args)
        candidate_times.append(time.perf_counter() - started)
    return {
        "case": name,
        "threads": threads,
        "repeats": repeats,
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "bitwise_equal": first_divergence is None,
        "assembly_nodes_equal": trace_equal,
        "first_divergence": first_divergence,
        "DN_gw_relative": abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
    }


def _run_outer(name, repeats=25, threads=2):
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
            FS.solve_kernel = solve_kernel_counted_assembly
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
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "digest_spectrum_equal": _digest(baseline_model.log10OmegaGW) == _digest(candidate_model.log10OmegaGW),
        "digest_DN_gw_equal": _digest(baseline_model.DN_gw) == _digest(candidate_model.DN_gw),
        "failure_equal": getattr(baseline_model, "fast_failure_reason", None)
        == getattr(candidate_model, "fast_failure_reason", None),
        "converged_equal": getattr(baseline_model, "SGWB_converge", False)
        == getattr(candidate_model, "SGWB_converge", False),
        "spectrum_dex_max": float(
            np.max(np.abs(np.asarray(candidate_model.log10OmegaGW) - np.asarray(baseline_model.log10OmegaGW)))
        ),
        "DN_gw_relative": abs(float(candidate_model.DN_gw[-1]) - float(baseline_model.DN_gw[-1]))
        / max(abs(float(baseline_model.DN_gw[-1])), 1e-300),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--outer", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    os.environ.setdefault("FAST_THREADS", str(args.threads))
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)
    cases = args.case or list(CASES)
    runner = _run_outer if args.outer else _run_kernel
    payload = {
        "candidate": "counted_assembly_state",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "scope": "full_outer" if args.outer else "kernel",
        "cases": [runner(case, args.repeats, args.threads) for case in cases],
    }
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
