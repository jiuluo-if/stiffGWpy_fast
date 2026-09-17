"""Benchmark a phase-substep exponential recurrence as a standalone prototype.

The production solver evaluates ``exp(z)`` independently at every phase
substep.  This prototype keeps the production transfer map, kink split,
assembly and tail matching unchanged, but advances the substep frequency with
one initial exponential and a fixed recurrence.  It is diagnostic only.
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

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.exact_background import fast_phi_s2_split  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}


@njit(inline='always', cache=True)
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


@njit(inline='always', cache=True)
def _phase_segment_recurrence(xh, yh, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * math.exp(z_mid) / phase_max))
        if n_sub < 1:
            n_sub = 1
    if n_sub == 1:
        return FS.scaled_step(xh, yh, z_mid, h_step)

    h_sub = h_step / n_sub
    dz = (z_end - z_start) / n_sub
    w = math.exp(z_start + 0.5 * dz)
    w_ratio = math.exp(dz)
    for _ in range(n_sub):
        xh, yh = _scaled_step_with_w(xh, yh, w, h_sub)
        w *= w_ratio
    return xh, yh


@njit(parallel=True, cache=True)
def solve_kernel_recurrence(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re,
):
    """Full-assembly twin of ``solve_kernel`` with only phase exp recurrence changed."""
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
                    xh, yh = _phase_segment_recurrence(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = _phase_segment_recurrence(
                        xh, yh, z_break, z_end, h_step - h_left, phase_max)
                else:
                    xh, yh = _phase_segment_recurrence(
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
                    lxh = xh
                    lyh = yh
                    last_z = zz
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
                handoff_eps[mode] = (abs(1.5 * Sv[kend] - 1.0) * e_z
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
                FS.assemble_tail(
                    Ogw, Oj, Opgw, mode, slot, kk2, coeff, eNz, fp_i, Pt,
                    ev_minus, fp_minus)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _prepared(case_name, threads):
    FS.apply_accuracy_mode('fast')
    FS.set_threads(threads)
    model = LCDM_SG(**CASES[case_name])
    FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal',
                      frequency_quadrature='pchip')
    Nv = np.asarray(model.Nv, dtype=np.float64)
    freqs = np.asarray(model.f, dtype=np.float64)
    dn_eff = float(model.cosmo_param['DN_eff'])
    Phi_grid, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = (
        fast_phi_s2_split(model, Nv, dn_eff, sigma_nodes=model.sigma))
    _, _, j0s, z0s, fp_minus = FS.prep_frequency_only(model, Nv, freqs)
    ev_minus = np.exp(-Nv)
    fp_freq = np.power(10.0, freqs)
    col_step = 8
    idx_out = np.unique(np.append(np.arange(0, len(Nv), col_step), len(Nv) - 1))
    P_t = (model.derived_param['A_t']
           * np.power((10.0 ** freqs) / FS.gp.f_piv,
                      model.derived_param['nt']))
    common = (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
        fp_minus, fp_freq, 1, len(idx_out), col_step, 0.005, 5.0,
        None, model.sigma, 0.25, kink_index, kink_fraction, phi_re,
    )
    return model, common


def _make_args(common):
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus, fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, _, Sv, phase_max, kink_index, kink_fraction, phi_re = common
    shape = (len(fp_freq), n_coarse)
    return (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
        fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
        np.zeros(shape), np.zeros(shape), np.zeros(shape), None, Sv,
        phase_max, np.full(len(fp_freq), -1.0), kink_index, kink_fraction,
        phi_re,
    )


def _metrics(model, baseline, candidate):
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    scale = max(float(np.max(np.abs(base_obs))), 1e-300)
    mask = np.abs(base_obs) > 1e-12 * scale
    rel = np.abs(cand_obs[mask] - base_obs[mask]) / np.maximum(
        np.abs(base_obs[mask]), 1e-300)
    base_dn = FS.integrate_frequency_pchip(model.f, base_obs)
    cand_dn = FS.integrate_frequency_pchip(model.f, cand_obs)
    return {
        'digest_equal': all(_digest(baseline[i]) == _digest(candidate[i])
                            for i in (16, 17, 18)),
        'max_abs_output': max(float(np.max(np.abs(candidate[i] - baseline[i])))
                            for i in (16, 17, 18)),
        'spectrum_rel_p50': float(np.percentile(rel, 50)),
        'spectrum_rel_p95': float(np.percentile(rel, 95)),
        'spectrum_rel_max': float(np.max(rel)),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
    }


def run_case(case_name, repeats, threads):
    model, common = _prepared(case_name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_recurrence(*candidate)
    accuracy = _metrics(model, baseline, candidate)
    baseline_times = []
    candidate_times = []
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
        solve_kernel_recurrence(*candidate)
        candidate_times.append(time.perf_counter() - start)
    base_median = statistics.median(baseline_times)
    cand_median = statistics.median(candidate_times)
    return {
        'case': case_name,
        'n_freq': len(model.f),
        'repeats': repeats,
        'baseline_median_ms': base_median * 1e3,
        'candidate_median_ms': cand_median * 1e3,
        'baseline_p95_ms': float(np.percentile(baseline_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
        'candidate_over_baseline': cand_median / base_median,
        'accuracy': accuracy,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), action='append')
    parser.add_argument('--repeats', type=int, default=30)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', default='docs/phase_recurrence_round_20260917.json')
    args = parser.parse_args()
    cases = args.case or list(CASES)
    records = [run_case(case, args.repeats, args.threads) for case in cases]
    payload = {
        'schema_version': 1,
        'experiment': 'standalone_phase_exp_recurrence',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': telemetry(workers=1, threads=args.threads,
                              concurrent_processes=1),
        'settings': {'repeats': args.repeats, 'cases': cases,
                     'threads': args.threads},
        'records': records,
    }
    output = ROOT / args.out
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
