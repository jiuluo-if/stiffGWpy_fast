"""Standalone fourth-order commutator-free Magnus block spike.

This prototype applies the two-exponential Gaussian-node CF4 composition over
four native background intervals.  A cheap commutator defect gate decides
whether a block is eligible; otherwise the production Cartesian transfer is
used exactly.  It is diagnostic only and never changes the production path.
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

SQRT3 = math.sqrt(3.0)
A_PLUS = 0.25 + SQRT3 / 6.0
A_MINUS = 0.25 - SQRT3 / 6.0
C_LEFT = 0.5 - SQRT3 / 6.0
C_RIGHT = 0.5 + SQRT3 / 6.0


@njit(inline="always")
def _combined_step(xh, yh, w1, w2, c1, c2, h):
    """Apply exp(h*(c1*A(w1)+c2*A(w2))) to one Cartesian state."""
    d = c1 + c2
    w = c1 * w1 + c2 * w2
    q = w * w - d * d
    if q >= 0.0:
        omega = math.sqrt(q)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega if omega != 0.0 else h
    else:
        x = math.sqrt(-q) * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return (c - d * si) * xh - w * si * yh, w * si * xh + (c + d * si) * yh


@njit(inline="always")
def _cf4_block(xh, yh, z_start, z_end, block_h):
    delta = z_end - z_start
    z1 = z_start + C_LEFT * delta
    z2 = z_start + C_RIGHT * delta
    w1 = math.exp(z1)
    w2 = math.exp(z2)
    xh, yh = _combined_step(xh, yh, w1, w2, A_PLUS, A_MINUS, block_h)
    return _combined_step(xh, yh, w1, w2, A_MINUS, A_PLUS, block_h)


@njit(inline="always")
def _cf4_defect(z_start, z_end, block_h, phase_cap):
    z_mid = 0.5 * (z_start + z_end)
    w_mid = math.exp(z_mid)
    w_delta = abs(math.exp(z_end) - math.exp(z_start))
    commutator = w_delta * (1.0 + w_mid)
    phase = block_h * max(w_mid, 1.0)
    return block_h * block_h * commutator, phase


@njit(parallel=True, cache=True)
def solve_kernel_cf4(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re, defect_threshold, block_phase_cap, accepted_counts,
):
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Pt = P_t[mode]
        fp_i = fp_freq[mode]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0) * S2inv[j0]
        k = j0
        zz = z0
        lxh, lyh, last_z = 0.0, yh, zz
        if k % col_step == 0 and assemble:
            slot = k // col_step
            if slot >= n_coarse - 1:
                slot = n_coarse - 1
            FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1,
                             S2[k], xh, yh, zz, Pt)
        if zz < z_tail:
            while k < nv - 1 and zz < z_tail:
                can_block = (k + 4 < nv and k != kink_index
                             and k + 1 != kink_index and k + 2 != kink_index
                             and k + 3 != kink_index
                             and k % col_step <= col_step - 4
                             and not (0.0 < kink_fraction < 1.0
                                      and k <= kink_index < k + 4))
                if can_block:
                    hs = h_arr[k] if h_arr is not None else h
                    block_h = hs * 4.0
                    z_end = z0 + Phi_grid[k + 4] - Phi0
                    defect, phase = _cf4_defect(zz, z_end, block_h, block_phase_cap)
                    if (defect <= defect_threshold and phase <= block_phase_cap
                            and z_end < z_tail):
                        xh, yh = _cf4_block(xh, yh, zz, z_end, block_h)
                        k += 4
                        zz = z0 + Phi_grid[k] - Phi0
                        accepted_counts[mode] += 1
                    else:
                        z_mid_step = z0 + Phi_mid[k] - Phi0
                        z_end_step = 2.0 * z_mid_step - zz
                        xh, yh = FS._phase_segment(
                            xh, yh, zz, z_end_step, hs, phase_max)
                        k += 1
                        zz = z0 + Phi_grid[k] - Phi0
                else:
                    hs = h_arr[k] if h_arr is not None else h
                    z_mid_step = z0 + Phi_mid[k] - Phi0
                    z_end_step = 2.0 * z_mid_step - zz
                    if k == kink_index and 0.0 < kink_fraction < 1.0:
                        z_break = z0 + phi_re - Phi0
                        z_end_step = z0 + Phi_grid[k + 1] - Phi0
                        h_left = hs * kink_fraction
                        xh, yh = FS._phase_segment(
                            xh, yh, zz, z_break, h_left, phase_max)
                        xh, yh = FS._phase_segment(
                            xh, yh, z_break, z_end_step,
                            hs - h_left, phase_max)
                    else:
                        xh, yh = FS._phase_segment(
                            xh, yh, zz, z_end_step, hs, phase_max)
                    k += 1
                    zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0 and assemble:
                    slot = k // col_step
                    if slot >= n_coarse - 1:
                        slot = n_coarse - 1
                    FS.assemble_main(
                        Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(
                        Ogw, Oj, Opgw, mode, n_coarse - 1,
                        S2[k], xh, yh, zz, Pt)
                if zz < z_tail:
                    lxh, lyh, last_z = xh, yh, zz
            kend = nv - 1 if zz < z_tail else k - 1
            if kend < j0:
                kend = j0
        else:
            kend = j0

        if kend < nv - 1:
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend]) if Sv is not None else 0.0
            amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                    + 2.0 * gamma * lxh * lyh * e_z)
            coeff = math.sqrt(0.5 * S2[kend] * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                handoff_eps[mode] = (abs(1.5 * Sv[kend] - 1.0) * e_z
                                      if Sv is not None else 0.0)
            slot_start = kend // col_step
            if slot_start >= n_coarse - 1:
                slot_start = n_coarse - 1
            while slot_start < n_coarse:
                kk = col_step * slot_start
                if slot_start == n_coarse - 1:
                    kk = nv - 1
                if kk > kend:
                    break
                slot_start += 1
            for slot in range(slot_start, n_coarse):
                if not (assemble or slot == n_coarse - 1):
                    continue
                kk = col_step * slot
                if slot == n_coarse - 1:
                    kk = nv - 1
                FS.assemble_tail(
                    Ogw, Oj, Opgw, mode, slot, kk, coeff, eNz,
                    fp_i, Pt, ev_minus, fp_minus)
        else:
            FS.assemble_main(
                Ogw, Oj, Opgw, mode, n_coarse - 1, S2[kend],
                xh, yh, zz, Pt)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _candidate_args(common, threshold, phase_cap):
    args = _make_args(common)
    accepted = np.zeros(len(args[5]), dtype=np.int64)
    return args + (float(threshold), float(phase_cap), accepted), accepted


def _run_kernel(case_name, threshold, repeats=30, threads=2, phase_cap=1.0):
    model, common = _prepared(case_name, threads)
    baseline = _make_args(common)
    candidate, accepted = _candidate_args(common, threshold, phase_cap)
    FS.solve_kernel(*baseline)
    solve_kernel_cf4(*candidate)
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    rel = np.abs(cand_obs - base_obs) / np.maximum(np.abs(base_obs), 1e-300)
    baseline_times, candidate_times = [], []
    for _ in range(repeats):
        for arr in baseline[16:19]:
            arr.fill(0.0)
        baseline[22].fill(-1.0)
        for arr in candidate[16:19]:
            arr.fill(0.0)
        candidate[22].fill(-1.0)
        accepted.fill(0)
        started = time.perf_counter()
        FS.solve_kernel(*baseline)
        baseline_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        solve_kernel_cf4(*candidate)
        candidate_times.append(time.perf_counter() - started)
    bmed = statistics.median(baseline_times)
    cmed = statistics.median(candidate_times)
    return {
        "case": case_name,
        "threshold": float(threshold),
        "phase_cap": float(phase_cap),
        "threads": threads,
        "repeats": repeats,
        "baseline_median_ms": bmed * 1e3,
        "candidate_median_ms": cmed * 1e3,
        "candidate_over_baseline": cmed / bmed,
        "baseline_p95_ms": float(np.percentile(baseline_times, 95) * 1e3),
        "candidate_p95_ms": float(np.percentile(candidate_times, 95) * 1e3),
        "accepted_blocks": int(np.sum(accepted)),
        "bitwise_equal": all(
            np.array_equal(baseline[i][:, -1], candidate[i][:, -1])
            for i in (16, 17, 18)
        ),
        "full_output_bitwise_equal": all(
            np.array_equal(baseline[i], candidate[i]) for i in (16, 17, 18)
        ),
        "spectrum_rel_p95": float(np.percentile(rel, 95)),
        "spectrum_rel_max": float(np.max(rel)),
        "digest_Ogw_equal": _digest(baseline[16][:, -1]) == _digest(candidate[16][:, -1]),
        "digest_Oj_equal": _digest(baseline[17][:, -1]) == _digest(candidate[17][:, -1]),
        "digest_Opgw_equal": _digest(baseline[18][:, -1]) == _digest(candidate[18][:, -1]),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--threshold", type=float, nargs="+", default=[0.0, 1e-3, 3e-3])
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--phase-cap", type=float, default=1.0)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = args.case or list(CASES)
    rows = [_run_kernel(case, threshold, args.repeats, args.threads, args.phase_cap)
            for threshold in args.threshold for case in cases]
    payload = {
        "candidate": "cf4_gaussian_commutator_free_magnus",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "rows": rows,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
