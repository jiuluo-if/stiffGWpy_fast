"""Standalone A/B for caching the continuous post-kink sigma endpoint.

The candidate only changes the diagnostic implementation of
``_sigma_node_limits``.  It caches the exact ``sigma_vec([N_re])`` result for
one solver call and reuses it on the next outer iteration.  Production code
is never modified.
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

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}
DIGEST_FIELDS = ('f', 'log10OmegaGW', 'DN_gw', 'g2', 'w2')


def _digest(model, name):
    array = np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64))
    return hashlib.sha256(array.tobytes()).hexdigest()


def _cached_sigma_factory():
    original = EB.sigma_vec
    cache = {}
    calls = {'probe': 0}

    def cached(N, model, dn_eff):
        array = np.asarray(N, dtype=float)
        # Only cache the fixed three-probe call in fast_phi_s2_split.  All
        # other sigma_vec calls retain the production implementation exactly.
        if array.shape == (3,):
            key = (id(model), float(dn_eff), array.tobytes())
            if key not in cache:
                cache[key] = original(array, model, dn_eff)
                calls['probe'] += 1
            return cache[key].copy()
        return original(array, model, dn_eff)

    return original, cached, calls


def run_once(case, candidate):
    original = EB.sigma_vec
    calls = None
    if candidate:
        original, replacement, calls = _cached_sigma_factory()
        EB.sigma_vec = replacement
    try:
        model = LCDM_SG(**case)
        start = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True)
        elapsed = time.perf_counter() - start
        row = {'seconds': elapsed, 'converged': bool(getattr(model, 'SGWB_converge', False)),
               'reason': getattr(model, 'fast_failure_reason', None)}
        row.update({'digest_' + name: _digest(model, name) for name in DIGEST_FIELDS})
        if calls is not None:
            row['probe_sigma_vec_calls'] = calls['probe']
        return row
    finally:
        EB.sigma_vec = original


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeats', type=int, default=50)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--cases', default='default,highT,stiff,high_kappa,lowT')
    parser.add_argument('--out', default='docs/kink_sigma_cache_round12_20260917.json')
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    FS.set_threads(args.threads)
    records = []
    for name in args.cases.split(','):
        case = CASES[name]
        run_once(case, False)  # compile/warm-up
        run_once(case, True)
        baseline = [run_once(case, False) for _ in range(args.repeats)]
        candidate = [run_once(case, True) for _ in range(args.repeats)]
        base_times = [row['seconds'] for row in baseline]
        cand_times = [row['seconds'] for row in candidate]
        base_last = baseline[-1]
        cand_last = candidate[-1]
        records.append({
            'case': name, 'repeats': args.repeats,
            'baseline_median_ms': statistics.median(base_times) * 1e3,
            'candidate_median_ms': statistics.median(cand_times) * 1e3,
            'baseline_p95_ms': float(np.percentile(base_times, 95)) * 1e3,
            'candidate_p95_ms': float(np.percentile(cand_times, 95)) * 1e3,
            'candidate_over_baseline': statistics.median(cand_times) / statistics.median(base_times),
            'digest_equal': all(base_last['digest_' + field] == cand_last['digest_' + field]
                                for field in DIGEST_FIELDS),
            'status_equal': (base_last['converged'] == cand_last['converged']
                             and base_last['reason'] == cand_last['reason']),
            'probe_sigma_vec_calls_candidate': cand_last.get('probe_sigma_vec_calls'),
            'digests': {field: {'baseline': base_last['digest_' + field],
                                'candidate': cand_last['digest_' + field]}
                        for field in DIGEST_FIELDS},
        })
        print(name, records[-1])
    payload = {
        'experiment': 'standalone_kink_sigma_endpoint_cache',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'resource': {'threads': args.threads, 'workers': 1, 'blas': '1'},
        'settings': vars(args), 'records': records,
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps({'commit': payload['commit'], 'out': args.out}, ensure_ascii=False))


if __name__ == '__main__':
    main()
