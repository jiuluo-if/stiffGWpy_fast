# -*- coding: utf-8 -*-
"""Standalone A/B for removing the unused Psi preparation buffer."""
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
from numba import njit

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.fast_sgwb import (  # noqa: E402
    ln10,
    prep_fast,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


@njit(cache=True)
def prep_kernel_without_psi(Nv, Sv, f_hor, freqs, h, ln10v,
                            Phi_grid, Phi_mid, S2, S2inv,
                            j0s, z0s, fp_minus):
    nv = len(Nv)
    nseg = nv - 1
    M = np.empty(nv)
    inv_h2 = 1.0 / (h * h)
    M[1] = (Sv[0] - 2.0 * Sv[1] + Sv[2]) * inv_h2
    M[nv - 2] = (Sv[nv - 3] - 2.0 * Sv[nv - 2] + Sv[nv - 1]) * inv_h2
    m = nv - 4
    aa = np.empty(m)
    bb = np.empty(m)
    cc = np.empty(m)
    dd = np.empty(m)
    for i in range(m):
        k = i + 2
        aa[i] = 1.0
        bb[i] = 4.0
        cc[i] = 1.0
        dd[i] = 6.0 * (Sv[k - 1] - 2.0 * Sv[k] + Sv[k + 1]) * inv_h2
    dd[0] -= M[1]
    dd[m - 1] -= M[nv - 2]
    for i in range(1, m):
        w = aa[i] / bb[i - 1]
        bb[i] -= w * cc[i - 1]
        dd[i] -= w * dd[i - 1]
    M[m + 1] = dd[m - 1] / bb[m - 1]
    for i in range(m - 2, -1, -1):
        M[i + 2] = (dd[i] - cc[i] * M[i + 3]) / bb[i]
    M[0] = 2.0 * M[1] - M[2]
    M[nv - 1] = 2.0 * M[nv - 2] - M[nv - 3]
    F = np.empty(nv)
    F[0] = 0.0
    for i in range(nseg):
        Mi = M[i]
        Mi1 = M[i + 1]
        a = (Mi1 - Mi) / (6.0 * h)
        b = Mi / 2.0
        c = (Sv[i + 1] - Sv[i]) / h - h * (2.0 * Mi + Mi1) / 6.0
        d = Sv[i]
        F[i + 1] = F[i] + (a / 4.0 * h**4 + b / 3.0 * h**3
                            + c / 2.0 * h**2 + d * h)
    N0 = Nv[0]
    for i in range(nv):
        Phi_grid[i] = 1.5 * F[i] - Nv[i] + N0
        S2[i] = math.exp(3.0 * F[i] - 4.0 * Nv[i])
        S2inv[i] = math.exp(-0.5 * (3.0 * F[i] - 4.0 * Nv[i]))
    for i in range(nv):
        if i < nseg:
            seg = i
            dx = 0.5 * h
        else:
            seg = nseg - 1
            dx = 1.5 * h
        Mi = M[seg]
        Mi1 = M[seg + 1]
        a = (Mi1 - Mi) / (6.0 * h)
        b = Mi / 2.0
        c = (Sv[seg + 1] - Sv[seg]) / h - h * (2.0 * Mi + Mi1) / 6.0
        d = Sv[seg]
        prim = F[seg] + (a / 4.0 * dx**4 + b / 3.0 * dx**3
                         + c / 2.0 * dx**2 + d * dx)
        xm = Nv[i] + 0.5 * h
        Phi_mid[i] = 1.5 * prim - xm + N0
    for j in range(nv):
        fp_minus[j] = math.exp(-f_hor[j] * ln10v)
    Nf = len(freqs)
    for mm in range(Nf):
        freq3 = freqs[mm] + 3.0
        lo = 0
        hi = nv
        while lo < hi:
            mid = (lo + hi) // 2
            if f_hor[mid] >= freq3:
                lo = mid + 1
            else:
                hi = mid
        j0 = lo - 1
        if j0 < 0:
            j0 = 0
        if j0 > nv - 1:
            j0 = nv - 1
        j0s[mm] = j0
        z0s[mm] = (freqs[mm] - f_hor[j0]) * ln10v


def digest(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def run(case, reps):
    model = LCDM_SG(**case)
    FS.gen_fast(model, 0.005)
    Nv = np.asarray(model.Nv, dtype=np.float64)
    model.construct_f(1.0)
    freqs = np.asarray(model.f, dtype=np.float64)
    h = 0.005
    nv = len(Nv)
    baseline = prep_fast(model, Nv, freqs, h, variable_grid=False)
    candidate = [np.empty(nv) for _ in range(4)]
    candidate += [np.empty(len(freqs), dtype=np.int64), np.empty(len(freqs)), np.empty(nv)]
    prep_kernel_without_psi(Nv, model.sigma, model.f_hor, freqs, h, ln10, *candidate)
    baseline_arrays = [baseline[i] for i in (2, 3, 5, 6, 7, 8, 9)]
    equal = all(np.array_equal(a, b) for a, b in zip(baseline_arrays, candidate))
    base_times = []
    cand_times = []
    for _ in range(reps):
        start = time.perf_counter()
        prep_fast(model, Nv, freqs, h, variable_grid=False)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        prep_kernel_without_psi(Nv, model.sigma, model.f_hor, freqs, h, ln10, *candidate)
        cand_times.append(time.perf_counter() - start)
    return {
        'digest_equal': equal,
        'baseline_median_s': statistics.median(base_times),
        'candidate_median_s': statistics.median(cand_times),
        'ratio': statistics.median(cand_times) / statistics.median(base_times),
        'digests': [digest(x) for x in candidate],
        'n_v': nv,
        'n_freq': len(freqs),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reps', type=int, default=50)
    parser.add_argument('--json', type=str, required=True)
    parser.add_argument('--case', choices=('lowT', 'highT', 'stiff', 'high_kappa'),
                        default='default')
    args = parser.parse_args()
    cases = {
        'default': {'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1e-2},
        'lowT': {'r': 1e-2, 'cr': 1, 'T_re': 1e1, 'kappa10': 1e-2},
        'highT': {'r': 1e-2, 'cr': 1, 'T_re': 1e4, 'kappa10': 1e-2},
        'stiff': {'r': 1e-1, 'cr': 1, 'T_re': 2e3, 'kappa10': 1e-2},
        'high_kappa': {'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1.0},
    }
    case = cases[args.case]
    result = {'commit': os.popen('git rev-parse HEAD').read().strip(),
              'reps': args.reps, 'case_id': args.case, 'case': case,
              'result': run(case, args.reps)}
    with open(args.json, 'w', encoding='utf-8') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
