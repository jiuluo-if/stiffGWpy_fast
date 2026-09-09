# -*- coding: utf-8 -*-
"""Isolated Phase-A comparison: plain fixed grid versus exact kink split."""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast._metrics import dex_abs, rel_linear_omega, signal_mask  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASE = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
FREQ_SUBSET = np.array([
    -18.4, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, -15.0, -14.5,
    -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0, -7.0, -6.0,
    -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
    6.0,
])


def solve(kink_split):
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='construct', threads=1)
    model = LCDM_SG(**CASE)
    start = time.perf_counter()
    FS.SGWB_iter_fast(model, tol=1e-6, config=cfg, kink_split=kink_split)
    return model, time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reps', type=int, default=5)
    parser.add_argument('--json', default=None)
    parser.add_argument('--subset-reference', action='store_true',
                        help='use the historical sparse reference grid (diagnostic only)')
    args = parser.parse_args()

    # Warm the two paths before collecting timing samples.
    solve(False)
    solve(True)
    rows = {}
    for name, flag in (('plain', False), ('kink_split', True)):
        samples = [solve(flag)[1] for _ in range(args.reps)]
        model, _ = solve(flag)
        rows[name] = {
            'warm_ms': [x * 1e3 for x in samples],
            'median_ms': statistics.median(samples) * 1e3,
            'p95_ms': float(np.percentile(samples, 95)),
            'n_nodes': int(len(model.Nv)),
            'n_freq': int(len(model.f)),
            'DN_gw': float(model.DN_gw[-1]),
            'DN_eff': float(model.cosmo_param['DN_eff']),
            'f_re': float(model.f_re),
        }

    reference_dn = {}
    for name in ('plain', 'kink_split'):
        # Keep the reference solve independent and evaluate it at the candidate's
        # converged DN_eff so only the transfer path is compared here.  The full
        # candidate frequency grid is the default; a sparse subset is retained
        # only as an explicitly labelled diagnostic because it changes DN_gw.
        candidate = solve(name == 'kink_split')[0]
        reference_freqs = (FREQ_SUBSET if args.subset_reference
                           else np.asarray(candidate.f, dtype=float))
        model = LCDM_SG(**CASE)
        ref = REF.run_reference(model, dn_eff=rows[name]['DN_eff'], freq_res=1.0,
                                z_tail=5.0, rtol=1e-11,
                                freq_subset=reference_freqs, self_consistent=False)
        candidate = solve(name == 'kink_split')[0]
        lo_ref = np.asarray(ref['log10OmegaGW'])
        if args.subset_reference:
            lo_test = np.interp(reference_freqs, np.asarray(candidate.f)[::-1],
                                np.asarray(candidate.log10OmegaGW)[::-1])
        else:
            lo_test = np.asarray(candidate.log10OmegaGW)
        mask = signal_mask(lo_ref)
        lin = rel_linear_omega(lo_ref[mask], lo_test[mask]) if mask.any() else np.array([0.0])
        all_rel = np.abs(10.0 ** lo_test - 10.0 ** lo_ref) / np.maximum(10.0 ** lo_ref, 1e-300)
        transition = np.abs(reference_freqs - rows[name]['f_re']) <= 0.5
        reference_dn[name] = {
            'reference_DN_gw': float(ref['DN_gw']),
            'DN_gw_rel': abs(rows[name]['DN_gw'] - float(ref['DN_gw'])) /
                         max(abs(float(ref['DN_gw'])), 1e-300),
            'spectrum_dex_median': float(np.median(dex_abs(lo_ref, lo_test))),
            'spectrum_dex_p95': float(np.percentile(dex_abs(lo_ref, lo_test), 95)),
            'spectrum_dex_max': float(np.max(dex_abs(lo_ref, lo_test))),
            'spectrum_rel_signal_median': float(np.median(lin)),
            'spectrum_rel_signal_max': float(np.max(lin)),
            'spectrum_rel_all_max': float(np.max(all_rel)),
            'transition_rel_max': float(np.max(all_rel[transition])) if transition.any() else None,
        }
    result = {'case': CASE, 'config': {'h': 0.02, 'z_tail': 5.0,
                                        'phase_max': 0.0, 'col_step': 8,
                                        'freq_grid': 'construct', 'threads': 1},
              'rows': rows, 'reference': reference_dn,
              'reference_grid': ('subset' if args.subset_reference else 'candidate_full')}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)


if __name__ == '__main__':
    main()
