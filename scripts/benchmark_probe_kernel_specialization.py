"""Benchmark a no-assembly specialization for the first outer probe."""
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

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}


@njit(parallel=True, cache=True)
def solve_probe(Nv, Phi_grid, Phi_mid, S2, S2inv,
                j0s, z0s, P_t, ev_minus, fp_minus, fp_freq,
                n_coarse, col_step, h, z_tail, Ogw, Oj, Opgw,
                h_arr=None, Sv=None, phase_max=0.0, handoff_eps=None,
                kink_index=-1, kink_fraction=0.0, phi_re=0.0):
    """Exact transfer path for assemble=0; only the final tail slot is emitted."""
    nv = len(Nv)
    for m in prange(len(j0s)):
        j0 = j0s[m]
        z0 = z0s[m]
        Pt = P_t[m]
        fp_i = fp_freq[m]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0)*S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh = 0.0
        lyh = yh
        last_z = zz
        if zz < z_tail:
            while k < nv-1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_node = z0 + Phi_grid[k] - Phi0
                    z_break = z0 + phi_re - Phi0
                    z_end = z0 + Phi_grid[k + 1] - Phi0
                    h_left = h_step * kink_fraction
                    h_right = h_step - h_left
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = FS._phase_segment(
                        xh, yh, z_break, z_end, h_right, phase_max)
                else:
                    z_node = z0 + Phi_grid[k] - Phi0
                    z_end = 2.0 * z_mid_step - z_node
                    xh, yh = FS._phase_segment(
                        xh, yh, z_node, z_end, h_step, phase_max)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if zz < z_tail:
                    lxh = xh
                    lyh = yh
                    last_z = zz
            if zz < z_tail:
                kend = nv-1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        if kend < nv-1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            if Sv is not None:
                gamma = FS._tail_match_gamma(Sv[kend])
                amp2 = (lxh*lxh + lyh*lyh*(1.0 + gamma*gamma*e_z*e_z)
                        + 2.0*gamma*lxh*lyh*e_z)
            else:
                amp2 = lxh*lxh + lyh*lyh
            coeff = math.sqrt(0.5*s2k*amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                if Sv is not None:
                    handoff_eps[m] = abs(1.5*Sv[kend] - 1.0)*e_z
                else:
                    handoff_eps[m] = 0.0
            FS.assemble_tail(Ogw, Oj, Opgw, m, n_coarse-1, nv-1,
                             coeff, eNz, fp_i, Pt, ev_minus, fp_minus)


def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def run(case, repeats, specialized):
    original = FS.solve_kernel
    if specialized:
        def wrapped(*args, **kwargs):
            if args[11] == 0:
                probe_args = list(args)
                probe_args.pop(11)
                # solve_probe omits the assemble flag but keeps all remaining args.
                return solve_probe(*probe_args, **kwargs)
            return original(*args, **kwargs)
        FS.solve_kernel = wrapped
    rows = []
    try:
        for _ in range(repeats):
            model = LCDM_SG(**CASES[case])
            start = time.perf_counter()
            result = FS.SGWB_iter_fast(model, kink_split=True)
            rows.append({
                'total_s': time.perf_counter() - start,
                'converged': result is model and bool(model.SGWB_converge),
                'digest_f': digest(model.f),
                'digest_spectrum': digest(model.log10OmegaGW),
                'digest_DN_gw': digest(model.DN_gw),
                'DN_gw_last': float(np.asarray(model.DN_gw)[-1]),
            })
    finally:
        FS.solve_kernel = original
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    baseline = run(args.case, args.repeats, specialized=False)
    candidate = run(args.case, args.repeats, specialized=True)
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'case_id': args.case,
        'threads': int(FS._THREADS),
        'repeats': args.repeats,
        'baseline_median_s': statistics.median(row['total_s'] for row in baseline),
        'candidate_median_s': statistics.median(row['total_s'] for row in candidate),
        'baseline_p95_s': float(np.percentile([row['total_s'] for row in baseline], 95)),
        'candidate_p95_s': float(np.percentile([row['total_s'] for row in candidate], 95)),
        'candidate_over_baseline': statistics.median(row['total_s'] for row in candidate) /
        statistics.median(row['total_s'] for row in baseline),
        'all_converged': all(row['converged'] for row in baseline + candidate),
        'digest_equal': all(
            left[key] == right[key]
            for left, right in zip(baseline, candidate)
            for key in ('digest_f', 'digest_spectrum', 'digest_DN_gw')
        ),
        'dn_max_abs': max(abs(left['DN_gw_last'] - right['DN_gw_last'])
                          for left, right in zip(baseline, candidate)),
    }
    with open(args.json, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
