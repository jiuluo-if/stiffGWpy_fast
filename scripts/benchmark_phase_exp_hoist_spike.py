"""Standalone spike for hoisting one repeated phase ``exp(z_mid)``.

In the production phase segment, ``_phase_substeps`` evaluates ``exp(z_mid)``
to choose ``n_sub`` and the following ``scaled_step`` evaluates the same value
again when ``n_sub == 1``.  This twin reuses that scalar only in the one-step
case.  Transfer-map formulas, phase subdivision, kink splitting, tail matching
and assembly are otherwise unchanged; production code is not modified.
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
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(inline='always', fastmath=True, cache=True)
def _scaled_step_with_w(xh, yh, w, h):
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


@njit(inline='always', fastmath=True, cache=True)
def _phase_segment_hoist(xh, yh, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    w_mid = math.exp(z_mid)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * w_mid / phase_max))
        if n_sub < 1:
            n_sub = 1
    if n_sub == 1:
        return _scaled_step_with_w(xh, yh, w_mid, h_step)
    h_sub = h_step / n_sub
    dz_half = z_mid - z_start
    for sub in range(n_sub):
        zs = z_start + dz_half * (2.0 * sub + 1.0) / n_sub
        xh, yh = FS.scaled_step(xh, yh, zs, h_sub)
    return xh, yh


@njit(parallel=True, fastmath=True, cache=True)
def solve_kernel_exp_hoist(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re,
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
        if k % col_step == 0:
            if assemble:
                slot = k // col_step
                if slot >= n_coarse - 1:
                    slot = n_coarse - 1
                FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1,
                             S2[k], xh, yh, zz, Pt)
        if zz < z_tail:
            while k < nv - 1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = z0 + Phi_grid[k] - Phi0
                z_end = z0 + Phi_grid[k + 1] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    xh, yh = _phase_segment_hoist(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = _phase_segment_hoist(
                        xh, yh, z_break, z_end, h_step - h_left, phase_max)
                else:
                    xh, yh = _phase_segment_hoist(
                        xh, yh, z_node, z_end, h_step, phase_max)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0:
                    if assemble:
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
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend]) if Sv is not None else 0.0
            amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                    + 2.0 * gamma * lxh * lyh * e_z)
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                handoff_eps[mode] = (
                    abs(1.5 * Sv[kend] - 1.0) * e_z if Sv is not None else 0.0)
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


def _run_kernel(name, repeats, threads):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_exp_hoist(*candidate)
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    rel = np.abs(cand_obs - base_obs) / np.maximum(np.abs(base_obs), 1e-300)
    base_dn = FS.integrate_frequency_pchip(model.f, base_obs)
    cand_dn = FS.integrate_frequency_pchip(model.f, cand_obs)
    base_times, cand_times = [], []
    for _ in range(repeats):
        args = _make_args(common)
        t0 = time.perf_counter()
        FS.solve_kernel(*args)
        base_times.append(time.perf_counter() - t0)
        args = _make_args(common)
        t0 = time.perf_counter()
        solve_kernel_exp_hoist(*args)
        cand_times.append(time.perf_counter() - t0)
    return {
        'case': name,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_s': statistics.median(base_times),
        'candidate_median_s': statistics.median(cand_times),
        'median_ratio': statistics.median(cand_times) / statistics.median(base_times),
        'baseline_p95_s': float(np.percentile(base_times, 95)),
        'candidate_p95_s': float(np.percentile(cand_times, 95)),
        'digest_equal': all(_digest(baseline[i]) == _digest(candidate[i]) for i in (16, 17, 18)),
        'spectrum_rel_p50': float(np.percentile(rel, 50)),
        'spectrum_rel_p95': float(np.percentile(rel, 95)),
        'spectrum_rel_max': float(np.max(rel)),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
    }


def _run_outer(name, repeats, threads):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    original = FS.solve_kernel
    base_times, cand_times = [], []
    base_model = cand_model = None
    try:
        for _ in range(repeats):
            base_model = LCDM_SG(**CASES[name])
            t0 = time.perf_counter()
            FS.SGWB_iter_fast(base_model, kink_split=True, freq_grid='goal',
                              frequency_quadrature='pchip')
            base_times.append(time.perf_counter() - t0)
            FS.solve_kernel = solve_kernel_exp_hoist
            cand_model = LCDM_SG(**CASES[name])
            t0 = time.perf_counter()
            FS.SGWB_iter_fast(cand_model, kink_split=True, freq_grid='goal',
                              frequency_quadrature='pchip')
            cand_times.append(time.perf_counter() - t0)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original
    base = np.asarray(base_model.log10OmegaGW, dtype=np.float64)
    cand = np.asarray(cand_model.log10OmegaGW, dtype=np.float64)
    delta = np.abs(cand - base)
    rel = np.abs(10.0 ** cand - 10.0 ** base) / np.maximum(10.0 ** base, 1e-300)
    bdn = float(np.asarray(base_model.DN_gw)[-1])
    cdn = float(np.asarray(cand_model.DN_gw)[-1])
    return {
        'case': name,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_s': statistics.median(base_times),
        'candidate_median_s': statistics.median(cand_times),
        'median_ratio': statistics.median(cand_times) / statistics.median(base_times),
        'baseline_p95_s': float(np.percentile(base_times, 95)),
        'candidate_p95_s': float(np.percentile(cand_times, 95)),
        'spectrum_dex_p50': float(np.percentile(delta, 50)),
        'spectrum_dex_p95': float(np.percentile(delta, 95)),
        'spectrum_dex_max': float(np.max(delta)),
        'spectrum_rel_p50': float(np.percentile(rel, 50)),
        'spectrum_rel_p95': float(np.percentile(rel, 95)),
        'spectrum_rel_max': float(np.max(rel)),
        'DN_gw_relative': abs(cdn - bdn) / max(abs(bdn), 1e-300),
        'digest_spectrum_equal': _digest(base) == _digest(cand),
        'digest_DN_gw_equal': _digest(np.asarray(base_model.DN_gw)) == _digest(
            np.asarray(cand_model.DN_gw)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', action='append', choices=sorted(CASES))
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--outer', action='store_true')
    parser.add_argument('--output')
    args = parser.parse_args()
    cases = args.case or list(CASES)
    runner = _run_outer if args.outer else _run_kernel
    result = {
        'candidate': 'phase_exp_hoist_fastmath',
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'production_unchanged': True,
        'scope': 'full_outer' if args.outer else 'kernel',
        'cases': [runner(case, args.repeats, args.threads) for case in cases],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        Path(args.output).write_text(text + '\n', encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
