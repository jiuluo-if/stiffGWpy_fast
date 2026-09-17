"""A/B benchmark for the opt-in fast derived-parameter cache."""
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


def _digest(value):
    return hashlib.sha256(np.ascontiguousarray(
        np.asarray(value, dtype=np.float64)).tobytes()).hexdigest()


def _uncached_frequency_only(model, Nv, freqs):
    fp_minus = np.empty(len(Nv), dtype=np.float64)
    j0s = np.empty(len(freqs), dtype=np.int64)
    z0s = np.empty(len(freqs), dtype=np.float64)
    FS.prep_frequency_kernel(model.f_hor, freqs, FS.ln10, j0s, z0s,
                             fp_minus)
    return model.sigma, model.f_hor, j0s, z0s, fp_minus


def _uncached_snapshot(Nv, sigma, f_hor):
    return Nv.copy(), sigma.copy(), f_hor.copy()


def _run_case(name, repeats, threads):
    original_frequency_only = FS.prep_frequency_only
    original_snapshot = FS._retain_outer_background_snapshot

    timings = {'uncached': [], 'cached': []}
    outputs = {'uncached': [], 'cached': []}
    try:
        for repeat in range(repeats):
            order = ('uncached', 'cached') if repeat % 2 == 0 else ('cached', 'uncached')
            for mode in order:
                FS.prep_frequency_only = (
                    _uncached_frequency_only if mode == 'uncached'
                    else original_frequency_only)
                FS._retain_outer_background_snapshot = (
                    _uncached_snapshot if mode == 'uncached'
                    else original_snapshot)
                model = LCDM_SG(**CASES[name])
                model._fast_derived_cache_disabled = (mode == 'uncached')
                start = time.perf_counter()
                FS.SGWB_iter_fast(model, kink_split=True)
                timings[mode].append((time.perf_counter() - start) * 1e3)
                outputs[mode].append({
                    'f': _digest(model.f),
                    'spectrum': _digest(model.log10OmegaGW),
                    'DN': _digest(model.DN_gw),
                    'g2': _digest(model.g2),
                    'w2': _digest(model.w2),
                    'spectrum_values': np.asarray(model.log10OmegaGW, dtype=float).copy(),
                    'DN_last': float(np.asarray(model.DN_gw)[-1]),
                    'failure': getattr(model, 'fast_failure_reason', None),
                    'converged': bool(getattr(model, 'SGWB_converge', False)),
                })
    finally:
        FS.prep_frequency_only = original_frequency_only
        FS._retain_outer_background_snapshot = original_snapshot
    base = outputs['uncached'][-1]
    candidate = outputs['cached'][-1]
    delta = np.abs(candidate['spectrum_values'] - base['spectrum_values'])
    return {
        'uncached_median_ms': statistics.median(timings['uncached']),
        'cached_median_ms': statistics.median(timings['cached']),
        'cached_over_uncached': statistics.median(timings['cached']) /
        statistics.median(timings['uncached']),
        'uncached_p95_ms': float(np.percentile(timings['uncached'], 95)),
        'cached_p95_ms': float(np.percentile(timings['cached'], 95)),
        'digest_equal': all(base[key] == candidate[key]
                            for key in ('f', 'spectrum', 'DN', 'g2', 'w2')),
        'status_equal': (base['failure'] == candidate['failure']
                         and base['converged'] == candidate['converged']),
        'spectrum_abs_dex_p50': float(np.percentile(delta, 50)),
        'spectrum_abs_dex_p95': float(np.percentile(delta, 95)),
        'spectrum_abs_dex_max': float(np.max(delta)),
        'dn_relative': abs(candidate['DN_last'] - base['DN_last']) /
        max(abs(base['DN_last']), 1e-300),
        'uncached_failure': base['failure'],
        'cached_failure': candidate['failure'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    result = _run_case(args.case, args.repeats, args.threads)
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'derived_param_cache_and_frequency_workspace',
        'production_path_changed': True,
        'case': args.case,
        'repeats': args.repeats,
        'threads': args.threads,
        'outer': result,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
