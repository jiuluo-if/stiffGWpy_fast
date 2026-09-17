"""Benchmark an allocation-reduced exact primitive prototype."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}


def fused_primitive(m, Nv, DN_eff, sigma_nodes=None):
    """Use the same formulas while reusing the midpoint work buffer."""
    Nv = np.asarray(Nv, dtype=float)
    h_arr = np.diff(Nv).astype(np.float64)
    if Nv.size < 2:
        return (np.zeros_like(Nv), np.zeros_like(Nv), np.ones_like(Nv),
                np.ones_like(Nv), -1, 0.0, 0.0)
    nodes_left, nodes_right = EB._sigma_node_limits(Nv, m, DN_eff, sigma_nodes)
    edge_sum = nodes_left[:-1] + nodes_right[1:]
    mid_sigma = 0.5 * edge_sum
    integral = h_arr * (nodes_left[:-1] + 4.0 * mid_sigma + nodes_right[1:]) / 6.0

    d = m.derived_param
    n_re = float(d['N_inf'] - d['N_re'])
    kink_index = int(np.searchsorted(Nv, n_re, side='right') - 1)
    kink_fraction = 0.0
    left_integral = None
    if 0 <= kink_index < Nv.size - 1:
        left = float(Nv[kink_index])
        right = float(Nv[kink_index + 1])
        kink_fraction = (n_re - left) / (right - left)
        if 0.0 < kink_fraction < 1.0:
            probes = np.array([
                left + 0.5 * (n_re - left), n_re,
                n_re + 0.5 * (right - n_re),
            ], dtype=float)
            sig_left, sig_re, sig_right = EB.sigma_vec(probes, m, DN_eff)
            left_h = n_re - left
            right_h = right - n_re
            left_integral = left_h * (nodes_left[kink_index] +
                                      4.0 * sig_left + 1.0) / 6.0
            right_integral = right_h * (sig_re + 4.0 * sig_right +
                                        nodes_right[kink_index + 1]) / 6.0
            integral[kink_index] = left_integral + right_integral
        else:
            kink_index = -1
            kink_fraction = 0.0

    F_nodes = np.concatenate(([0.0], np.cumsum(integral)))
    if kink_index >= 0 and left_integral is not None:
        phi_re = 1.5 * (F_nodes[kink_index] + left_integral) - n_re + Nv[0]
    else:
        phi_re = 0.0
    # Reuse the midpoint buffer for the quarter-point interpolation.
    mid_sigma = 0.75 * nodes_left[:-1] + 0.25 * nodes_right[1:]
    F_mid = F_nodes[:-1] + h_arr * (
        nodes_left[:-1] + 4.0 * mid_sigma + 0.5 * edge_sum) / 12.0
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - 0.5 * (Nv[:-1] + Nv[1:]) + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64),
            kink_index, float(kink_fraction), float(phi_re))


def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=50)
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    model = LCDM_SG(**CASES[args.case])
    FS.apply_accuracy_mode('fast')
    FS.gen_fast(model, kink_split=True)
    Nv = model.Nv.astype(np.float64)
    sigma = model.sigma.copy()
    dn = model.cosmo_param['DN_eff']
    baseline = EB.fast_phi_s2_split
    for _ in range(3):
        baseline(model, Nv, dn, sigma_nodes=sigma)
        fused_primitive(model, Nv, dn, sigma_nodes=sigma)
    timings = {'baseline': [], 'fused': []}
    outputs = {}
    for name, fn in (('baseline', baseline), ('fused', fused_primitive)):
        for _ in range(args.repeats):
            start = time.perf_counter()
            value = fn(model, Nv, dn, sigma_nodes=sigma)
            timings[name].append(time.perf_counter() - start)
            outputs[name] = value
    base = outputs['baseline']
    cand = outputs['fused']
    arrays = [0, 1, 2, 3]
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'case_id': args.case,
        'threads': int(FS._THREADS),
        'repeats': args.repeats,
        'median_s': {name: statistics.median(values) for name, values in timings.items()},
        'p95_s': {name: float(np.percentile(values, 95)) for name, values in timings.items()},
        'ratio_fused_over_baseline': statistics.median(timings['fused']) /
        statistics.median(timings['baseline']),
        'digest_equal': all(digest(base[i]) == digest(cand[i]) for i in arrays),
        'max_abs': max(float(np.max(np.abs(base[i] - cand[i]))) for i in arrays),
        'max_rel': max(float(np.max(np.abs(base[i] - cand[i]) /
                              np.maximum(np.abs(base[i]), 1e-300))) for i in arrays),
    }
    with open(args.json, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
