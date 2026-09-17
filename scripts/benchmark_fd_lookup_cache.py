"""Standalone cache spike for repeated FD interpolator lookup."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402

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

_FD_CACHE = EB._fd_from_ref()


def cached_fd_lookup():
    """Return the already-loaded FD interpolation callables by identity."""
    return _FD_CACHE


def _digest(value):
    return hashlib.sha256(np.ascontiguousarray(
        np.asarray(value, dtype=np.float64)).tobytes()).hexdigest()


def _run_case(name, candidate):
    original = EB._fd_from_ref
    if candidate:
        EB._fd_from_ref = cached_fd_lookup
    try:
        model = LCDM_SG(**CASES[name])
        start = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True)
        elapsed = (time.perf_counter() - start) * 1e3
        return model, elapsed
    finally:
        EB._fd_from_ref = original


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
    baseline_digest_unique = {
        field: len({_digest(getattr(model, field)) for model in outputs['baseline']})
        for field in fields
    }
    candidate_digest_unique = {
        field: len({_digest(getattr(model, field)) for model in outputs['candidate']})
        for field in fields
    }
    spec_delta = np.abs(np.asarray(candidate.log10OmegaGW) -
                        np.asarray(baseline.log10OmegaGW))
    dn_base = float(np.asarray(baseline.DN_gw)[-1])
    dn_candidate = float(np.asarray(candidate.DN_gw)[-1])
    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'fd_lookup_cache_standalone_spike',
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
        'baseline_digest_unique': baseline_digest_unique,
        'candidate_digest_unique': candidate_digest_unique,
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
