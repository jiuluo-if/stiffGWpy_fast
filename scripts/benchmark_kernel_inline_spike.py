"""Standalone benchmark for inlining the scalar tensor transfer map.

This diagnostic keeps the production phase subdivision, kink split, tail
matching, and assembly unchanged.  The only change is that the small
``scaled_step`` expression is inlined into the phase loop, avoiding a Numba
call boundary.  It must not be used as a production solver without a fresh
accuracy audit.
"""
from __future__ import annotations

import argparse
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
    _digest,
    _make_args,
    _metrics,
    _prepared,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(inline='always', cache=True)
def _scaled_step_inline(xh, yh, z_mid, h):
    w = math.exp(z_mid)
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


@njit(inline='always', cache=True)
def _phase_segment_inline(xh, yh, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * math.exp(z_mid) / phase_max))
        if n_sub < 1:
            n_sub = 1
    if n_sub == 1:
        return _scaled_step_inline(xh, yh, z_mid, h_step)
    h_sub = h_step / n_sub
    dz_half = z_mid - z_start
    for sub in range(n_sub):
        zs = z_start + dz_half * (2.0 * sub + 1.0) / n_sub
        xh, yh = _scaled_step_inline(xh, yh, zs, h_sub)
    return xh, yh


@njit(parallel=True, cache=True)
def solve_kernel_inline(
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
        lxh = 0.0
        lyh = yh
        last_z = zz
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
                    xh, yh = _phase_segment_inline(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = _phase_segment_inline(
                        xh, yh, z_break, z_end, h_step - h_left, phase_max)
                else:
                    xh, yh = _phase_segment_inline(
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


def _run_case(name, repeats, threads):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_inline(*candidate)
    metrics = _metrics(model, baseline, candidate)
    base_times = []
    cand_times = []
    for _ in range(repeats):
        args = _make_args(common)
        t0 = time.perf_counter()
        FS.solve_kernel(*args)
        base_times.append(time.perf_counter() - t0)
        args = _make_args(common)
        t0 = time.perf_counter()
        solve_kernel_inline(*args)
        cand_times.append(time.perf_counter() - t0)
    metrics.update({
        'case': name,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_s': statistics.median(base_times),
        'candidate_median_s': statistics.median(cand_times),
        'median_ratio': statistics.median(cand_times) / statistics.median(base_times),
        'baseline_p95_s': float(np.percentile(base_times, 95)),
        'candidate_p95_s': float(np.percentile(cand_times, 95)),
        'baseline_digest': [_digest(a) for a in baseline[16:19]],
        'candidate_digest': [_digest(a) for a in candidate[16:19]],
    })
    return metrics


def _run_outer_case(name, repeats, threads):
    """Run the same full outer solve with the kernel twin monkeypatched."""
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    original = FS.solve_kernel
    baseline_times = []
    candidate_times = []
    baseline_model = None
    candidate_model = None
    try:
        for _ in range(repeats):
            baseline_model = LCDM_SG(**CASES[name])
            t0 = time.perf_counter()
            FS.SGWB_iter_fast(
                baseline_model, kink_split=True, freq_grid='goal',
                frequency_quadrature='pchip')
            baseline_times.append(time.perf_counter() - t0)
            FS.solve_kernel = solve_kernel_inline
            candidate_model = LCDM_SG(**CASES[name])
            t0 = time.perf_counter()
            FS.SGWB_iter_fast(
                candidate_model, kink_split=True, freq_grid='goal',
                frequency_quadrature='pchip')
            candidate_times.append(time.perf_counter() - t0)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original
    base_spectrum = np.asarray(baseline_model.log10OmegaGW, dtype=np.float64)
    cand_spectrum = np.asarray(candidate_model.log10OmegaGW, dtype=np.float64)
    base_dn = float(np.asarray(baseline_model.DN_gw)[-1])
    cand_dn = float(np.asarray(candidate_model.DN_gw)[-1])
    delta = np.abs(cand_spectrum - base_spectrum)
    rel = np.abs(10.0 ** cand_spectrum - 10.0 ** base_spectrum) / np.maximum(
        np.abs(10.0 ** base_spectrum), 1e-300)
    return {
        'case': name,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_s': statistics.median(baseline_times),
        'candidate_median_s': statistics.median(candidate_times),
        'median_ratio': statistics.median(candidate_times) / statistics.median(baseline_times),
        'baseline_p95_s': float(np.percentile(baseline_times, 95)),
        'candidate_p95_s': float(np.percentile(candidate_times, 95)),
        'spectrum_dex_p50': float(np.percentile(delta, 50)),
        'spectrum_dex_p95': float(np.percentile(delta, 95)),
        'spectrum_dex_max': float(np.max(delta)),
        'spectrum_rel_p50': float(np.percentile(rel, 50)),
        'spectrum_rel_p95': float(np.percentile(rel, 95)),
        'spectrum_rel_max': float(np.max(rel)),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
        'digest_spectrum_equal': _digest(base_spectrum) == _digest(cand_spectrum),
        'digest_DN_gw_equal': _digest(np.asarray(baseline_model.DN_gw)) == _digest(
            np.asarray(candidate_model.DN_gw)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', action='append', choices=sorted(CASES))
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--output', type=str)
    parser.add_argument('--outer', action='store_true')
    args = parser.parse_args()
    cases = args.case or list(CASES)
    runner = _run_outer_case if args.outer else _run_case
    result = {
        'candidate': 'inline_scaled_step',
        'commit': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'production_unchanged': True,
        'scope': 'full_outer' if args.outer else 'kernel',
        'cases': [runner(name, args.repeats, args.threads) for name in cases],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as handle:
            handle.write(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
