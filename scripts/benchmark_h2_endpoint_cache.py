# -*- coding: utf-8 -*-
"""验证 goal-grid 内 H2 端点缓存的逐位输出与 runtime 影响。"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

POINTS = ['default', 'lowT', 'highT', 'stiff', 'high_kappa',
          'sobol_000', 'sobol_002', 'sobol_006', 'sobol_010', 'sobol_015',
          'edge_r_hi', 'edge_dnre_hi', 'edge_kap_hi']
DIGEST_FIELDS = ('f', 'log10OmegaGW', 'DN_gw', 'g2', 'w2')


def _digest(value):
    array = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    return hashlib.sha256(array.tobytes()).hexdigest()


def _old_grid_builder():
    source = subprocess.check_output(
        ['git', 'show', 'HEAD:stiffgwpy_fast/freq_adaptive.py'],
        text=True, encoding='utf-8',
    )
    start = source.index('def grid_independent_freqs(')
    end = source.index('\ndef goal_oriented_freqs(', start)
    method = textwrap.dedent(source[start:end])
    namespace = vars(__import__('stiffgwpy_fast.freq_adaptive', fromlist=['*']))
    exec(method, namespace)
    return namespace['grid_independent_freqs']


def _run(name, builder):
    original = FA.grid_independent_freqs
    FA.grid_independent_freqs = builder
    try:
        model = LCDM_SG(**CASES[name])
        start = time.perf_counter()
        FS.apply_accuracy_mode('fast')
        FS.set_threads(2)
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
        elapsed = time.perf_counter() - start
    finally:
        FA.grid_independent_freqs = original
    return {
        'elapsed_s': elapsed,
        'digests': {field: _digest(getattr(model, field)) for field in DIGEST_FIELDS},
        'DN_gw': float(np.asarray(model.DN_gw)[-1]),
        'n_freq': len(model.f),
        'converged': bool(getattr(model, 'SGWB_converge', False)),
        'failure': getattr(model, 'fast_failure_reason', None),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reps', type=int, default=25)
    parser.add_argument('--json', type=Path, required=True)
    args = parser.parse_args()
    old_builder = _old_grid_builder()
    new_builder = FA.grid_independent_freqs
    rows = []
    for name in POINTS:
        _run(name, old_builder)
        _run(name, new_builder)
        baseline = [_run(name, old_builder) for _ in range(args.reps)]
        candidate = [_run(name, new_builder) for _ in range(args.reps)]
        base_median = statistics.median(x['elapsed_s'] for x in baseline)
        candidate_median = statistics.median(x['elapsed_s'] for x in candidate)
        rows.append({
            'case': name,
            'reps': args.reps,
            'digest_equal': baseline[0]['digests'] == candidate[0]['digests'],
            'baseline_median_ms': base_median * 1e3,
            'candidate_median_ms': candidate_median * 1e3,
            'speed_ratio_candidate_over_baseline': candidate_median / base_median,
            'baseline': baseline[0],
            'candidate': candidate[0],
        })
    payload = {
        'experiment': 'goal_grid_h2_endpoint_cache_ab',
        'generated_commit': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True,
        ).strip(),
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'points': POINTS,
        'rows': rows,
    }
    args.json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({
        'commit': payload['generated_commit'],
        'digest_mismatch': sum(not row['digest_equal'] for row in rows),
        'speed_ratios': {row['case']: row['speed_ratio_candidate_over_baseline'] for row in rows},
        'out': str(args.json),
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
