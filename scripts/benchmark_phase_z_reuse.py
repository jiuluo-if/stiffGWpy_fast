"""Standalone exact loop-state reuse A/B for the tensor propagation kernel.

The candidate reuses the already-computed ``zz`` state for the next loop
condition and interval start.  It removes repeated arithmetic only; production
code and the transfer map are unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

from scripts.benchmark_phase_recurrence import (  # noqa: E402
    CASES,
    _make_args,
    _metrics,
    _prepared,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def phase_state_reuse(z0, phi0, phi_grid):
    """Reference-level helper for the exact repeated expression contract."""
    return z0 + np.asarray(phi_grid, dtype=np.float64) - phi0


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


@njit(parallel=True, cache=True)
def solve_kernel_z_reuse(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re,
):
    """Twin of ``solve_kernel`` using ``zz`` as the loop-state cache."""
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
                FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh,
                                 zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1, S2[k], xh, yh,
                             zz, Pt)
        if zz < z_tail:
            while k < nv - 1 and zz < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = zz
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(
                        xh, yh, z_break,
                        z0 + Phi_grid[k + 1] - Phi0,
                        h_step - h_left, phase_max)
                else:
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node,
                        2.0 * z_mid_step - z_node,
                        h_step, phase_max)
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
                        Ogw, Oj, Opgw, mode, n_coarse - 1, S2[k], xh, yh,
                        zz, Pt)
                if zz < z_tail:
                    lxh, lyh, last_z = xh, yh, zz
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
                amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                        + 2.0 * gamma * lxh * lyh * e_z)
            else:
                amp2 = lxh * lxh + lyh * lyh
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
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
                    Ogw, Oj, Opgw, mode, slot, kk2, coeff, eNz, fp_i, Pt,
                    ev_minus, fp_minus)


def _run_case(name, repeats, threads):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_z_reuse(*candidate)
    accuracy = _metrics(model, baseline, candidate)
    base_times = []
    cand_times = []
    for _ in range(repeats):
        for array in baseline[16:19]:
            array.fill(0.0)
        baseline[22].fill(-1.0)
        for array in candidate[16:19]:
            array.fill(0.0)
        candidate[22].fill(-1.0)
        start = time.perf_counter()
        FS.solve_kernel(*baseline)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        solve_kernel_z_reuse(*candidate)
        cand_times.append(time.perf_counter() - start)
    base = statistics.median(base_times)
    cand = statistics.median(cand_times)
    return {
        'case': name,
        'n_freq': len(model.f),
        'repeats': repeats,
        'baseline_median_ms': base * 1e3,
        'candidate_median_ms': cand * 1e3,
        'baseline_p95_ms': float(np.percentile(base_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(cand_times, 95) * 1e3),
        'candidate_over_baseline': cand / base,
        'accuracy': accuracy,
    }


def _run_outer_case(name, repeats, threads):
    """Run the full outer solve with the standalone kernel twin."""
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    original = FS.solve_kernel
    base_times = []
    cand_times = []
    baseline_model = None
    candidate_model = None
    try:
        for _ in range(repeats):
            baseline_model = LCDM_SG(**CASES[name])
            start = time.perf_counter()
            FS.SGWB_iter_fast(
                baseline_model, kink_split=True, freq_grid='goal',
                frequency_quadrature='pchip')
            base_times.append(time.perf_counter() - start)
            FS.solve_kernel = solve_kernel_z_reuse
            candidate_model = LCDM_SG(**CASES[name])
            start = time.perf_counter()
            FS.SGWB_iter_fast(
                candidate_model, kink_split=True, freq_grid='goal',
                frequency_quadrature='pchip')
            cand_times.append(time.perf_counter() - start)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original
    base_spectrum = np.asarray(baseline_model.log10OmegaGW, dtype=np.float64)
    cand_spectrum = np.asarray(candidate_model.log10OmegaGW, dtype=np.float64)
    delta = np.abs(cand_spectrum - base_spectrum)
    base_dn = float(np.asarray(baseline_model.DN_gw)[-1])
    cand_dn = float(np.asarray(candidate_model.DN_gw)[-1])
    return {
        'case': name,
        'threads': threads,
        'repeats': repeats,
        'baseline_median_s': statistics.median(base_times),
        'candidate_median_s': statistics.median(cand_times),
        'median_ratio': statistics.median(cand_times) /
        statistics.median(base_times),
        'baseline_p95_s': float(np.percentile(base_times, 95)),
        'candidate_p95_s': float(np.percentile(cand_times, 95)),
        'spectrum_dex_p50': float(np.percentile(delta, 50)),
        'spectrum_dex_p95': float(np.percentile(delta, 95)),
        'spectrum_dex_max': float(np.max(delta)),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
        'digest_spectrum_equal': _digest(base_spectrum) == _digest(cand_spectrum),
        'digest_DN_gw_equal': _digest(np.asarray(baseline_model.DN_gw)) ==
        _digest(np.asarray(candidate_model.DN_gw)),
        'baseline_failure': getattr(baseline_model, 'fast_failure_reason', None),
        'candidate_failure': getattr(candidate_model, 'fast_failure_reason', None),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), action='append')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--outer', action='store_true')
    parser.add_argument('--out', default='docs/phase_z_reuse_round20_20260917.json')
    args = parser.parse_args(argv)
    cases = args.case or list(CASES)
    runner = _run_outer_case if args.outer else _run_case
    records = [runner(case, args.repeats, args.threads) for case in cases]
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_phase_loop_state_reuse',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'repeats': args.repeats, 'cases': cases,
                     'threads': args.threads, 'scope':
                     'full_outer' if args.outer else 'kernel'},
        'records': records,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
