"""Standalone probe for an adiabaticity-triggered WKB handoff.

The candidate keeps the production phase propagation and tail formula, but stops
at the first native node where |omega'/omega^2| is below a selected threshold.
It is diagnostic only and never changes the formal fast path.
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
from numba import njit, prange

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import (  # noqa: E402
    CASES,
    _make_args,
    _prepared,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(inline='always')
def _adiabaticity(sigma_value, z_value):
    """Return |omega'/omega^2| for omega=exp(z)."""
    return abs(1.5 * sigma_value - 1.0) * math.exp(-z_value)


@njit(parallel=True, cache=True)
def solve_kernel_adiabatic(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re, epsilon_trigger,
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
        zz = z0 + Phi_grid[k] - Phi0
        lxh, lyh, last_z = 0.0, yh, zz
        if k % col_step == 0 and assemble:
            slot = k // col_step
            if slot >= n_coarse - 1:
                slot = n_coarse - 1
            FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1,
                             S2[k], xh, yh, zz, Pt)

        trigger_now = (_adiabaticity(Sv[k], zz) <= epsilon_trigger
                       if Sv is not None else False)
        if zz < z_tail and not trigger_now:
            while k < nv - 1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                z_node = z0 + Phi_grid[k] - Phi0
                if _adiabaticity(Sv[k], z_node) <= epsilon_trigger:
                    break
                h_step = h_arr[k] if h_arr is not None else h
                z_end = z0 + Phi_grid[k + 1] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(
                        xh, yh, z_break, z_end, h_step - h_left, phase_max)
                else:
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_end, h_step, phase_max)
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
                if zz < z_tail and _adiabaticity(Sv[k], zz) > epsilon_trigger:
                    lxh, lyh, last_z = xh, yh, zz
            kend = nv - 1 if zz < z_tail and k == nv - 1 else k - 1
            if kend < j0:
                kend = j0
        else:
            kend = j0

        if kend < nv - 1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend]) if Sv is not None else 0.0
            amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                    + 2.0 * gamma * lxh * lyh * e_z)
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                handoff_eps[mode] = _adiabaticity(Sv[kend], last_z)
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
                    Ogw, Oj, Opgw, mode, slot, kk2, coeff, eNz, fp_i, Pt,
                    ev_minus, fp_minus)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _run_kernel(case_name, threshold, repeats, threads):
    model, common = _prepared(case_name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common) + (threshold,)
    FS.solve_kernel(*baseline)
    solve_kernel_adiabatic(*candidate)
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    rel = np.abs(cand_obs - base_obs) / np.maximum(np.abs(base_obs), 1e-300)
    base_dn = FS.integrate_frequency_pchip(model.f, base_obs)
    cand_dn = FS.integrate_frequency_pchip(model.f, cand_obs)
    baseline_times, candidate_times = [], []
    for _ in range(repeats):
        for array in baseline[16:19]:
            array.fill(0.0)
        baseline[22].fill(-1.0)
        for array in candidate[16:19]:
            array.fill(0.0)
        candidate[22].fill(-1.0)
        start = time.perf_counter()
        FS.solve_kernel(*baseline)
        baseline_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        solve_kernel_adiabatic(*candidate)
        candidate_times.append(time.perf_counter() - start)
    bmed = statistics.median(baseline_times)
    cmed = statistics.median(candidate_times)
    return {
        'case': case_name,
        'threshold': threshold,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_ms': bmed * 1e3,
        'candidate_median_ms': cmed * 1e3,
        'candidate_over_baseline': cmed / bmed,
        'spectrum_rel_p50': float(np.percentile(rel, 50)),
        'spectrum_rel_p95': float(np.percentile(rel, 95)),
        'spectrum_rel_max': float(np.max(rel)),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
        'digest_equal': all(_digest(baseline[i]) == _digest(candidate[i])
                            for i in (16, 17, 18)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', type=float, nargs='+',
                        default=[2e-2, 1e-2, 5e-3, 1e-3])
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    rows = [_run_kernel(name, threshold, args.repeats, args.threads)
            for threshold in args.threshold
            for name in CASES]
    payload = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                          cwd=ROOT, text=True).strip(),
        'candidate': 'adiabaticity_trigger_handoff',
        'production_path_changed': False,
        'settings': vars(args),
        'rows': rows,
    }
    (ROOT / args.out).write_text(json.dumps(payload, indent=2) + '\n',
                                 encoding='utf-8')
    print(json.dumps(payload, indent=2))


if __name__ == '__main__':
    main()
