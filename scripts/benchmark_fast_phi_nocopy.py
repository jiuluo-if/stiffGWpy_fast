"""Standalone no-copy primitive spike for the formal kink-split path."""
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


def no_copy_primitive(m, Nv, DN_eff, sigma_nodes=None):
    """Evaluate the production primitive without redundant node copies.

    The formal kink-split grid normally has no node exactly at ``N_re``.  In
    that case both one-sided node arrays in ``_sigma_node_limits`` are equal,
    so the read-only input can serve both views.  The exact-node branch keeps
    the original copies and one-sided convention.
    """
    Nv = np.asarray(Nv, dtype=float)
    h_arr = np.diff(Nv).astype(np.float64)
    if Nv.size < 2:
        return (np.zeros_like(Nv), np.zeros_like(Nv), np.ones_like(Nv),
                np.ones_like(Nv), -1, 0.0, 0.0)
    if sigma_nodes is None:
        nodes = EB.sigma_vec(Nv, m, DN_eff)
    else:
        nodes = np.asarray(sigma_nodes, dtype=float)
        if nodes.shape != Nv.shape:
            raise ValueError('sigma_nodes must match Nv shape')
    n_re = float(m.derived_param['N_inf'] - m.derived_param['N_re'])
    at_re = np.isclose(Nv, n_re, rtol=0.0, atol=1e-12)
    if np.any(at_re):
        nodes_left = nodes.copy()
        nodes_right = nodes.copy()
        nodes_left[at_re] = 1.0
        nodes_right[at_re] = EB.sigma_vec(
            np.array([n_re]), m, DN_eff)[0]
    else:
        nodes_left = nodes
        nodes_right = nodes

    mid_sigma = 0.5 * (nodes_left[:-1] + nodes_right[1:])
    quarter_sigma = 0.75 * nodes_left[:-1] + 0.25 * nodes_right[1:]
    integral = h_arr * (nodes_left[:-1] + 4.0 * mid_sigma + nodes_right[1:]) / 6.0

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
    F_mid = F_nodes[:-1] + h_arr * (
        nodes_left[:-1] + 4.0 * quarter_sigma + mid_sigma) / 12.0
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - 0.5 * (Nv[:-1] + Nv[1:]) + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64),
            kink_index, float(kink_fraction), float(phi_re))


def _digest(value):
    return hashlib.sha256(np.ascontiguousarray(
        np.asarray(value, dtype=np.float64)).tobytes()).hexdigest()


def _run_case(name, candidate):
    original = EB.fast_phi_s2_split
    if candidate:
        EB.fast_phi_s2_split = no_copy_primitive
    try:
        model = LCDM_SG(**CASES[name])
        start = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True)
        elapsed = (time.perf_counter() - start) * 1e3
        return model, elapsed
    finally:
        EB.fast_phi_s2_split = original


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    name = args.case
    _run_case(name, False)
    _run_case(name, True)
    timings = {'baseline': [], 'candidate': []}
    outputs = {'baseline': [], 'candidate': []}
    for repeat in range(args.repeats):
        order = (False, True) if repeat % 2 == 0 else (True, False)
        for candidate in order:
            key = 'candidate' if candidate else 'baseline'
            model, elapsed = _run_case(name, candidate)
            timings[key].append(elapsed)
            outputs[key].append(model)
    baseline = outputs['baseline'][-1]
    candidate = outputs['candidate'][-1]
    fields = ('f', 'log10OmegaGW', 'DN_gw', 'g2', 'w2')
    baseline_digests = {field: _digest(getattr(baseline, field)) for field in fields}
    candidate_digests = {field: _digest(getattr(candidate, field)) for field in fields}
    spec_delta = np.abs(np.asarray(candidate.log10OmegaGW) -
                        np.asarray(baseline.log10OmegaGW))
    dn_base = float(np.asarray(baseline.DN_gw)[-1])
    dn_candidate = float(np.asarray(candidate.DN_gw)[-1])
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'fast_phi_s2_no_copy_standalone_spike',
        'production_path_changed': False,
        'case': name,
        'repeats': args.repeats,
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'baseline_median_ms': statistics.median(timings['baseline']),
        'candidate_median_ms': statistics.median(timings['candidate']),
        'baseline_p95_ms': float(np.percentile(timings['baseline'], 95)),
        'candidate_p95_ms': float(np.percentile(timings['candidate'], 95)),
        'candidate_over_baseline': statistics.median(timings['candidate']) /
        statistics.median(timings['baseline']),
        'digest_equal': baseline_digests == candidate_digests,
        'baseline_digests': baseline_digests,
        'candidate_digests': candidate_digests,
        'spectrum_abs_dex_p50': float(np.percentile(spec_delta, 50)),
        'spectrum_abs_dex_p95': float(np.percentile(spec_delta, 95)),
        'spectrum_abs_dex_max': float(np.max(spec_delta)),
        'dn_relative': abs(dn_candidate - dn_base) / max(abs(dn_base), 1e-300),
        'baseline_failure': getattr(baseline, 'fast_failure_reason', None),
        'candidate_failure': getattr(candidate, 'fast_failure_reason', None),
        'baseline_converged': bool(getattr(baseline, 'SGWB_converge', False)),
        'candidate_converged': bool(getattr(candidate, 'SGWB_converge', False)),
    }
    with open(args.json, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
