"""Standalone full-kernel twin for phase-increment sin/cos recurrence.

Only the high-frequency phase factors in the production transfer map are
changed.  Subdivision, z endpoints, kink handling, assembly and tail matching
remain the production expressions.  This module is diagnostic until all
fixed-DN and full outer gates pass.
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


@njit(cache=True, inline='always')
def _phase_segment_increment(
    xh, yh, z_start, z_end, h_step, phase_max, reanchor,
    phase_age, phase_c, phase_s, angle_previous,
):
    """Advance one production segment with periodic exact phase anchors."""
    z_mid = 0.5 * (z_start + z_end)
    n_sub = FS._phase_substeps(h_step, z_mid, phase_max)
    if n_sub == 1:
        xh, yh = FS.scaled_step(xh, yh, z_mid, h_step)
        return xh, yh, reanchor, 1.0, 0.0, 0.0

    if reanchor < 1:
        reanchor = 1
    h_sub = h_step / n_sub
    dz_half = z_mid - z_start
    for sub in range(n_sub):
        zs = z_start + dz_half * (2.0 * sub + 1.0) / n_sub
        w = math.exp(zs)
        w2 = w * w
        if w2 < 1.0:
            omega = math.sqrt(1.0 - w2)
            x = omega * h_sub
            c = 1.0 + 0.5 * x * x
            si = h_sub * (1.0 + x * x / 6.0)
            xh, yh = (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh
            phase_age = reanchor
            phase_c = 1.0
            phase_s = 0.0
            angle_previous = 0.0
            continue

        omega = math.sqrt(w2 - 1.0)
        angle = omega * h_sub
        if phase_age >= reanchor:
            phase_c = math.cos(angle)
            phase_s = math.sin(angle)
            phase_age = 1
        else:
            delta = angle - angle_previous
            delta2 = delta * delta
            cd = 1.0 - 0.5 * delta2
            sd = delta - delta2 * delta / 6.0
            phase_c, phase_s = (
                phase_c * cd - phase_s * sd,
                phase_s * cd + phase_c * sd,
            )
            phase_age += 1
        angle_previous = angle
        si = phase_s / omega
        xh, yh = (
            (phase_c - si) * xh - w * si * yh,
            w * si * xh + (phase_c + si) * yh,
        )
    return xh, yh, phase_age, phase_c, phase_s, angle_previous


@njit(parallel=True, cache=True)
def solve_kernel_phase_increment(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re, reanchor,
):
    """Full production-kernel twin with only phase sin/cos recurrence changed."""
    nv = len(Nv)
    for m in prange(len(j0s)):
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
        phase_age = reanchor
        phase_c = 1.0
        phase_s = 0.0
        angle_previous = 0.0
        if k % col_step == 0:
            if assemble:
                slot = k // col_step
                if slot >= n_coarse - 1:
                    slot = n_coarse - 1
                FS.assemble_main(Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, m, n_coarse - 1,
                             S2[k], xh, yh, zz, Pt)

        if zz < z_tail:
            while k < nv - 1 and zz < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = zz
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    h_right = h_step - h_left
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(
                        xh, yh, z_break, z0 + Phi_grid[k + 1] - Phi0,
                        h_right, phase_max)
                    phase_age = reanchor
                    phase_c = 1.0
                    phase_s = 0.0
                    angle_previous = 0.0
                else:
                    z_end = 2.0 * z_mid_step - z_node
                    (xh, yh, phase_age, phase_c, phase_s,
                     angle_previous) = _phase_segment_increment(
                        xh, yh, z_node, z_end, h_step, phase_max,
                        reanchor, phase_age, phase_c, phase_s,
                        angle_previous)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0:
                    if assemble:
                        slot = k // col_step
                        if slot >= n_coarse - 1:
                            slot = n_coarse - 1
                        FS.assemble_main(
                            Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(
                        Ogw, Oj, Opgw, m, n_coarse - 1,
                        S2[k], xh, yh, zz, Pt)
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
                amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                        + 2.0 * gamma * lxh * lyh * e_z)
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
                FS.assemble_tail(
                    Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt,
                    ev_minus, fp_minus)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _candidate_args(common):
    args = list(_make_args(common))
    return tuple(args)


def run_case(case_name, repeats=25, threads=2, reanchor=32):
    model, common = _prepared(case_name, threads)
    baseline = _candidate_args(common)
    candidate = _candidate_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_phase_increment(*candidate, reanchor)
    accuracy = _metrics(model, baseline, candidate)
    accuracy['handoff_eps_max'] = float(np.max(np.abs(candidate[22])))
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
        solve_kernel_phase_increment(*candidate, reanchor)
        cand_times.append(time.perf_counter() - start)
    bm = statistics.median(base_times)
    cm = statistics.median(cand_times)
    return {
        'case': case_name,
        'reanchor': reanchor,
        'repeats': repeats,
        'n_freq': len(model.f),
        'baseline_median_ms': bm * 1e3,
        'candidate_median_ms': cm * 1e3,
        'baseline_p95_ms': float(np.percentile(base_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(cand_times, 95) * 1e3),
        'candidate_over_baseline': cm / bm,
        'accuracy': accuracy,
        'baseline_digests': {
            'Ogw': _digest(baseline[16]), 'Oj': _digest(baseline[17]),
            'Opgw': _digest(baseline[18]),
        },
        'candidate_digests': {
            'Ogw': _digest(candidate[16]), 'Oj': _digest(candidate[17]),
            'Opgw': _digest(candidate[18]),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), action='append')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--reanchor', type=int, default=32)
    parser.add_argument('--out', default='docs/phase_increment_kernel_twin_round47_20260923.json')
    args = parser.parse_args()
    cases = args.case or list(CASES)
    records = [run_case(case, args.repeats, args.threads, args.reanchor)
               for case in cases]
    payload = {
        'schema_version': 1,
        'experiment': 'round47_phase_increment_full_kernel_twin',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(threads=args.threads),
        'settings': {'repeats': args.repeats, 'cases': cases,
                     'threads': args.threads, 'reanchor': args.reanchor},
        'records': records,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
