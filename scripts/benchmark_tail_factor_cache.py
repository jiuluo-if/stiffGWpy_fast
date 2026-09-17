"""Standalone tail-assembly factor cache A/B.

The candidate precomputes ``exp(-Nv) * exp(-f_hor * ln10)`` once per grid
node and reuses it in the tail-column assembly.  Production code is unchanged;
the benchmark is valid only if the output digest, accuracy and failure state
remain unchanged across the target regimes.
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


def tail_factor(ev_minus, fp_minus):
    """Build the grid-only tail factor used by the standalone candidate."""
    return np.asarray(ev_minus, dtype=np.float64) * np.asarray(
        fp_minus, dtype=np.float64)


@njit(cache=True)
def _assemble_tail_cached(Ogw, Oj, Opgw, mode, slot, coeff, eNz, fp_i, Pt,
                          ev_value, tail_static):
    Th = coeff * eNz * ev_value
    oj = -Th * Th / 3.0 * Pt
    xf = coeff * eNz * fp_i * tail_static
    op = xf * xf / 36.0 * Pt
    Oj[mode, slot] = oj
    Opgw[mode, slot] = op
    Ogw[mode, slot] = 3.0 * op + oj


@njit(parallel=True, cache=True)
def solve_kernel_tail_factor(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
    Ogw, Oj, Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index,
    kink_fraction, phi_re, tail_static,
):
    """Production-kernel twin with only the tail factor expression changed."""
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
            while k < nv - 1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = z0 + Phi_grid[k] - Phi0
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
                if k % col_step == 0:
                    if assemble:
                        slot = k // col_step
                        if slot >= n_coarse - 1:
                            slot = n_coarse - 1
                        FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh,
                                         yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1, S2[k],
                                     xh, yh, zz, Pt)
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
                    abs(1.5 * Sv[kend] - 1.0) * e_z
                    if Sv is not None else 0.0)
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
                _assemble_tail_cached(
                    Ogw, Oj, Opgw, mode, slot, coeff, eNz, fp_i, Pt,
                    ev_minus[kk2], tail_static[kk2])


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _run_case(name, repeats, threads):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    tail_static = tail_factor(baseline[8], baseline[9])
    candidate_args = candidate + (tail_static,)
    FS.solve_kernel(*baseline)
    solve_kernel_tail_factor(*candidate_args)
    accuracy = _metrics(model, baseline, candidate_args)
    base_times = []
    cand_times = []
    for _ in range(repeats):
        for array in baseline[16:19]:
            array.fill(0.0)
        baseline[22].fill(-1.0)
        for array in candidate_args[16:19]:
            array.fill(0.0)
        candidate_args[22].fill(-1.0)
        start = time.perf_counter()
        FS.solve_kernel(*baseline)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        solve_kernel_tail_factor(*candidate_args)
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
        'tail_factor_digest': _digest(tail_static),
    }


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), action='append')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', default='docs/tail_factor_cache_round19_20260917.json')
    args = parser.parse_args(argv)
    cases = args.case or list(CASES)
    records = [_run_case(case, args.repeats, args.threads) for case in cases]
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_tail_factor_cache',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'repeats': args.repeats, 'cases': cases,
                     'threads': args.threads},
        'records': records,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
