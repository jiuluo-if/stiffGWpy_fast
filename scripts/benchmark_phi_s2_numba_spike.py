"""Standalone Numba prototype for the kink-split Phi/S2 primitive.

This module is intentionally not imported by production code.  It keeps the
existing Python/NumPy preparation and moves only the array arithmetic and
prefix integration into a Numba helper so that bitwise and end-to-end costs
can be measured independently.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time

import numpy as np
from numba import njit

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
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                    T_re=1e3, kappa10=1e-3),
    'positive_tilt': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                         T_re=2e3, kappa10=1e-2),
    'negative_tilt': dict(r=1e-2, cr=1, n_t=-0.4,
                         T_re=2e3, kappa10=1e-2),
    'sobol_000': dict(r=3.4845505440304265e-06, n_t=-0.12252988945692778,
                     cr=1.0, T_re=34989.571710260374,
                     DN_re=18.787082182243466, kappa10=5.28549585835966e-06),
    'sobol_002': dict(r=0.008513729433124471, n_t=-0.3924837075173855,
                     cr=1.0, T_re=118660.64687585819,
                     DN_re=11.174977160990238, kappa10=0.0017982905763510582),
    'sobol_006': dict(r=0.002986405620277968, n_t=-0.16131426580250263,
                     cr=0.0, T_re=2832.147926998191,
                     DN_re=3.932188209146261, kappa10=0.01617670894806547),
    'sobol_010': dict(r=0.08346659409929896, n_t=-0.024900889955461025,
                     cr=0.0, T_re=4977.754805490329,
                     DN_re=26.626055845990777, kappa10=0.15852849886128825),
    'sobol_015': dict(r=1.5975020128853773e-06, n_t=0.18561450485140085,
                     cr=0.0, T_re=114872.20639509047,
                     DN_re=5.934050939977169, kappa10=2.79610302157714e-05),
}


@njit(cache=True)
def numba_primitive_from_nodes(Nv, nodes_left, nodes_right,
                               kink_index, kink_fraction,
                               sig_left, sig_re, sig_right):
    """Evaluate the existing fast primitive formulas in one compiled loop."""
    n = len(Nv)
    h_arr = np.diff(Nv)
    mid_sigma = 0.5 * (nodes_left[:-1] + nodes_right[1:])
    integral = h_arr * (
        nodes_left[:-1] + 4.0 * mid_sigma + nodes_right[1:]) / 6.0
    left_integral = 0.0
    if 0 <= kink_index < n - 1 and 0.0 < kink_fraction < 1.0:
        left_h = h_arr[kink_index] * kink_fraction
        right_h = h_arr[kink_index] - left_h
        left_integral = left_h * (
            nodes_left[kink_index] + 4.0 * sig_left + 1.0) / 6.0
        right_integral = right_h * (
            sig_re + 4.0 * sig_right + nodes_right[kink_index + 1]) / 6.0
        integral[kink_index] = left_integral + right_integral

    F_nodes = np.empty(n, dtype=np.float64)
    F_nodes[0] = 0.0
    for i in range(n - 1):
        F_nodes[i + 1] = F_nodes[i] + integral[i]
    F_mid = F_nodes[:-1] + h_arr * (
        nodes_left[:-1] + 4.0 * (0.75 * nodes_left[:-1]
                                  + 0.25 * nodes_right[1:])
        + mid_sigma) / 12.0
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - 0.5 * (Nv[:-1] + Nv[1:]) + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    return Phi_grid, Phi_mid, S2, S2inv, F_nodes, left_integral


def numba_fast_phi_s2_split(m, Nv, DN_eff, sigma_nodes=None):
    """Python adapter matching ``exact_background.fast_phi_s2_split``."""
    Nv = np.asarray(Nv, dtype=float)
    if Nv.size < 2:
        zeros = np.zeros_like(Nv)
        ones = np.ones_like(Nv)
        return zeros, zeros, ones, ones, -1, 0.0, 0.0
    nodes_left, nodes_right = EB._sigma_node_limits(
        Nv, m, DN_eff, sigma_nodes)
    d = m.derived_param
    n_re = float(d['N_inf'] - d['N_re'])
    kink_index = int(np.searchsorted(Nv, n_re, side='right') - 1)
    kink_fraction = 0.0
    sig_left = sig_re = sig_right = 0.0
    if 0 <= kink_index < Nv.size - 1:
        left = float(Nv[kink_index])
        right = float(Nv[kink_index + 1])
        kink_fraction = (n_re - left) / (right - left)
        if 0.0 < kink_fraction < 1.0:
            probes = np.array([
                left + 0.5 * (n_re - left), n_re,
                n_re + 0.5 * (right - n_re)], dtype=float)
            sig_left, sig_re, sig_right = EB.sigma_vec(
                probes, m, DN_eff)
        else:
            kink_index = -1
            kink_fraction = 0.0
    (Phi_grid, Phi_mid, S2, S2inv, F_nodes,
     left_integral) = numba_primitive_from_nodes(
         Nv, nodes_left, nodes_right, kink_index, kink_fraction,
         sig_left, sig_re, sig_right)
    if kink_index >= 0 and left_integral != 0.0:
        phi_re = 1.5 * (F_nodes[kink_index] + left_integral) - n_re + Nv[0]
    else:
        phi_re = 0.0
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64),
            kink_index, float(kink_fraction), float(phi_re))


def _digest(value):
    return hashlib.sha256(np.ascontiguousarray(
        np.asarray(value, dtype=np.float64)).tobytes()).hexdigest()


def _primitive_ab(model, repeats):
    FS.gen_fast(model, kink_split=True)
    Nv = model.Nv.astype(np.float64)
    sigma = model.sigma.copy()
    dn = model.cosmo_param['DN_eff']
    baseline = EB.fast_phi_s2_split
    # Compile and warm both paths before timing.
    baseline(model, Nv, dn, sigma_nodes=sigma)
    numba_fast_phi_s2_split(model, Nv, dn, sigma_nodes=sigma)
    timings = {'baseline': [], 'candidate': []}
    outputs = {}
    for name, fn in (('baseline', baseline), ('candidate', numba_fast_phi_s2_split)):
        for _ in range(repeats):
            start = time.perf_counter()
            outputs[name] = fn(model, Nv, dn, sigma_nodes=sigma)
            timings[name].append(time.perf_counter() - start)
    base = outputs['baseline']
    cand = outputs['candidate']
    delta = [np.abs(np.asarray(cand[i]) - np.asarray(base[i])) for i in range(4)]
    return {
        'baseline_median_ms': statistics.median(timings['baseline']) * 1e3,
        'candidate_median_ms': statistics.median(timings['candidate']) * 1e3,
        'candidate_over_baseline': statistics.median(timings['candidate']) /
        statistics.median(timings['baseline']),
        'digest_equal': all(_digest(base[i]) == _digest(cand[i]) for i in range(4)),
        'max_abs': max(float(np.max(value)) for value in delta),
        'max_rel': max(float(np.max(value / np.maximum(np.abs(base[i]), 1e-300)))
                    for i, value in enumerate(delta)),
    }


def _outer_ab(name, repeats):
    original = EB.fast_phi_s2_split
    timings = {'baseline': [], 'candidate': []}
    outputs = {'baseline': [], 'candidate': []}
    try:
        for _ in range(repeats):
            order = ('baseline', 'candidate') if _ % 2 == 0 else ('candidate', 'baseline')
            for mode in order:
                EB.fast_phi_s2_split = (original if mode == 'baseline'
                                        else numba_fast_phi_s2_split)
                model = LCDM_SG(**CASES[name])
                start = time.perf_counter()
                result = FS.SGWB_iter_fast(model, kink_split=True)
                timings[mode].append(time.perf_counter() - start)
                outputs[mode].append({
                    'result_is_model': result is model,
                    'converged': bool(getattr(model, 'SGWB_converge', False)),
                    'failure': getattr(model, 'fast_failure_reason', None),
                    'f': _digest(model.f),
                    'spectrum': _digest(model.log10OmegaGW),
                    'DN': _digest(model.DN_gw),
                    'g2': _digest(model.g2),
                    'w2': _digest(model.w2),
                    'DN_last': float(np.asarray(model.DN_gw)[-1]),
                })
    finally:
        EB.fast_phi_s2_split = original
    base = outputs['baseline'][-1]
    cand = outputs['candidate'][-1]
    return {
        'baseline_median_ms': statistics.median(timings['baseline']) * 1e3,
        'candidate_median_ms': statistics.median(timings['candidate']) * 1e3,
        'candidate_over_baseline': statistics.median(timings['candidate']) /
        statistics.median(timings['baseline']),
        'digest_equal': all(base[key] == cand[key]
                            for key in ('f', 'spectrum', 'DN', 'g2', 'w2')),
        'status_equal': (base['result_is_model'] == cand['result_is_model']
                         and base['converged'] == cand['converged']
                         and base['failure'] == cand['failure']),
        'dn_relative': abs(cand['DN_last'] - base['DN_last']) /
        max(abs(base['DN_last']), 1e-300),
        'baseline_failure': base['failure'],
        'candidate_failure': cand['failure'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(int(os.environ.get('FAST_THREADS', '2')))
    model = LCDM_SG(**CASES[args.case])
    primitive = _primitive_ab(model, args.repeats)
    outer = _outer_ab(args.case, args.repeats)
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'standalone_numba_phi_s2_primitive',
        'production_path_changed': False,
        'case': args.case,
        'repeats': args.repeats,
        'threads': int(FS._THREADS),
        'primitive': primitive,
        'outer': outer,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
