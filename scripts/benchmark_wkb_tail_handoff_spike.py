"""Standalone WKB-corrected early-tail handoff prototype.

The candidate stops Cartesian propagation at an earlier ``z_tail`` and applies
Oracle C's first-order boundary correction to the frozen-tail transfer.  It is
an evidence-only kernel twin; production ``solve_kernel`` is not modified.
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

import numpy as np
from numba import njit, prange

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from scripts.benchmark_phase_recurrence import _make_args, _prepared  # noqa: E402
from scripts.benchmark_prufer_oracle import spectrum_prufer  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


@njit(parallel=True, cache=True)
def solve_kernel_wkb_tail(args, correction_enabled):
    (Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
     fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
     Ogw, Oj, Opgw, h_arr, Sv, phase_max, handoff_eps,
     kink_index, kink_fraction, phi_re) = args
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Pt = P_t[mode]
        fp_i = fp_freq[mode]
        Phi0 = Phi_grid[j0]
        xh = 0.0
        yh = math.exp(z0) * S2inv[j0]
        k = j0
        zz = z0
        lxh = 0.0
        lyh = yh
        last_z = zz
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
                # This prototype uses the formal uniform kink-split grid;
                # ``h_arr`` is intentionally None on that path.
                h_step = h
                z_node = z0 + Phi_grid[k] - Phi0
                z_end = z0 + Phi_grid[k + 1] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_break,
                                                h_left, phase_max)
                    xh, yh = FS._phase_segment(xh, yh, z_break, z_end,
                                                h_step - h_left, phase_max)
                else:
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_end,
                                                h_step, phase_max)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0 and assemble:
                    slot = k // col_step
                    if slot >= n_coarse - 1:
                        slot = n_coarse - 1
                    FS.assemble_main(Ogw, Oj, Opgw, mode, slot,
                                     S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1,
                                     S2[k], xh, yh, zz, Pt)
                if zz < z_tail:
                    lxh = xh
                    lyh = yh
                    last_z = zz
        if k >= nv - 1:
            kend = nv - 1
        else:
            kend = max(j0, k - 1)
        if kend < nv - 1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend]) if Sv is not None else 0.0
            amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                    + 2.0 * gamma * lxh * lyh * e_z)
            coeff = math.sqrt(0.5 * s2k * amp2)
            if correction_enabled:
                omega = math.exp(last_z)
                theta = math.atan2(lxh, lyh)
                correction = 1.0 + math.sin(2.0 * theta) / omega
                if correction > 0.0:
                    coeff *= math.sqrt(correction)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                handoff_eps[mode] = (abs(1.5 * Sv[kend] - 1.0) * e_z
                                      if Sv is not None else 0.0)
            slot_start = kend // col_step
            if slot_start >= n_coarse - 1:
                slot_start = n_coarse - 1
            for slot in range(slot_start, n_coarse):
                if not (assemble or slot == n_coarse - 1):
                    continue
                kk2 = col_step * slot
                if slot == n_coarse - 1:
                    kk2 = nv - 1
                if kk2 < kend:
                    continue
                FS.assemble_tail(Ogw, Oj, Opgw, mode, slot, kk2, coeff,
                                 eNz, fp_i, Pt, ev_minus, fp_minus)


def digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def clear_outputs(args):
    for index in (16, 17, 18):
        args[index].fill(0.0)
    args[22].fill(-1.0)


def stats(values):
    array = np.asarray(values, dtype=float)
    return {'p50': float(np.percentile(array, 50)),
            'p95': float(np.percentile(array, 95)),
            'max': float(np.max(array)), 'n': int(array.size)}


def run_case(name, z_tail, repeats, threads, oracle_z_tail):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate_common = list(common)
    candidate_common[15] = float(z_tail)
    candidate = _make_args(tuple(candidate_common))
    FS.solve_kernel(*baseline)
    solve_kernel_wkb_tail(candidate, True)
    base_omega = baseline[16][:, -1] - baseline[17][:, -1]
    cand_omega = candidate[16][:, -1] - candidate[17][:, -1]
    rel = np.abs(cand_omega - base_omega) / np.maximum(np.abs(base_omega), 1e-300)
    dn_base = FS.integrate_frequency_pchip(model.f, base_omega)
    dn_cand = FS.integrate_frequency_pchip(model.f, cand_omega)
    oracle = spectrum_prufer(
        model, np.asarray(model.f, dtype=float),
        float(model.cosmo_param['DN_eff']), z_tail=oracle_z_tail, rtol=1e-10)
    oracle_omega = oracle['Ogw'] - oracle['Oj']
    base_oracle_rel = np.abs(base_omega - oracle_omega) \
        / np.maximum(np.abs(oracle_omega), 1e-300)
    cand_oracle_rel = np.abs(cand_omega - oracle_omega) \
        / np.maximum(np.abs(oracle_omega), 1e-300)
    dn_oracle = FS.integrate_frequency_pchip(model.f, oracle_omega)
    base_times = []
    cand_times = []
    for _ in range(repeats):
        clear_outputs(baseline)
        clear_outputs(candidate)
        start = time.perf_counter()
        FS.solve_kernel(*baseline)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        solve_kernel_wkb_tail(candidate, True)
        cand_times.append(time.perf_counter() - start)
    return {
        'case': name, 'z_tail': float(z_tail), 'repeats': repeats,
        'baseline_median_ms': statistics.median(base_times) * 1e3,
        'candidate_median_ms': statistics.median(cand_times) * 1e3,
        'candidate_over_baseline': statistics.median(cand_times) / statistics.median(base_times),
        'spectrum_relative': stats(rel),
        'baseline_vs_oracle': stats(base_oracle_rel),
        'candidate_vs_oracle': stats(cand_oracle_rel),
        'DN_base': float(dn_base), 'DN_candidate': float(dn_cand),
        'DN_relative': abs(dn_cand - dn_base) / max(abs(dn_base), 1e-300),
        'DN_oracle': float(dn_oracle),
        'DN_base_vs_oracle_relative': abs(dn_base - dn_oracle) / max(abs(dn_oracle), 1e-300),
        'DN_candidate_vs_oracle_relative': abs(dn_cand - dn_oracle) / max(abs(dn_oracle), 1e-300),
        'baseline_digest': digest(base_omega),
        'candidate_digest': digest(cand_omega),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--z-tail', type=float, nargs='+', default=[3.0, 3.5, 4.0, 4.5])
    parser.add_argument('--cases', default='default,highT,stiff,high_kappa,lowT')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--oracle-z-tail', type=float, default=10.0)
    parser.add_argument('--out', default='docs/wkb_tail_handoff_round13_20260917.json')
    args = parser.parse_args()
    rows = []
    for z_tail in args.z_tail:
        for name in args.cases.split(','):
            row = run_case(name, z_tail, args.repeats, args.threads,
                           args.oracle_z_tail)
            rows.append(row)
            print(row)
    payload = {
        'experiment': 'standalone_wkb_corrected_early_tail_handoff',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': {'numba_threads': args.threads, 'blas_threads': 1, 'workers': 1},
        'settings': vars(args), 'records': rows,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps({'commit': payload['commit'], 'out': args.out}, ensure_ascii=False))


if __name__ == '__main__':
    main()
