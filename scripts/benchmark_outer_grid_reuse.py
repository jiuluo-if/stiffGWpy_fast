"""验证 outer iteration 是否可以安全复用首轮 goal frequency grid。"""

from __future__ import annotations

import hashlib
import json
import os
import statistics
import sys
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


def _run_once(name, reuse_grid):
    original = FA.goal_oriented_freqs
    cached = [None]
    calls = [0]

    def goal(*args, **kwargs):
        calls[0] += 1
        if reuse_grid and cached[0] is not None:
            return cached[0].copy()
        value = np.asarray(original(*args, **kwargs), dtype=np.float64).copy()
        if reuse_grid:
            cached[0] = value
        return value

    FA.goal_oriented_freqs = goal
    try:
        model = LCDM_SG(**CASES[name])
        FS.apply_accuracy_mode('fast')
        FS.set_threads(2)
        start = time.perf_counter()
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
        elapsed = time.perf_counter() - start
    finally:
        FA.goal_oriented_freqs = original
    if not getattr(model, 'SGWB_converge', False):
        return {'converged': False, 'failure': getattr(model, 'fast_failure_reason', None),
                'seconds': elapsed, 'goal_calls': calls[0]}
    return {
        'converged': True,
        'failure': getattr(model, 'fast_failure_reason', None),
        'seconds': elapsed,
        'goal_calls': calls[0],
        'DN_gw': float(np.asarray(model.DN_gw)[-1]),
        'digests': {field: _digest(getattr(model, field)) for field in DIGEST_FIELDS},
        '_values': {field: np.asarray(getattr(model, field), dtype=float).copy()
                    for field in DIGEST_FIELDS},
    }


def run_case(name, repeats):
    baseline = _run_once(name, False)
    candidate = _run_once(name, True)
    if not baseline.get('converged') or not candidate.get('converged'):
        return {'case': name, 'baseline': baseline, 'candidate': candidate,
                'digest_equal': False, 'speed_ratio': None}
    base_times = []
    candidate_times = []
    # 预热两条路径，然后交替采样，隔离进程初始化和 JIT 影响。
    for _ in range(repeats):
        start = time.perf_counter()
        _run_once(name, False)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _run_once(name, True)
        candidate_times.append(time.perf_counter() - start)
    base_median = statistics.median(base_times)
    candidate_median = statistics.median(candidate_times)
    base_values = baseline.pop('_values')
    candidate_values = candidate.pop('_values')
    field_deltas = {}
    for field in DIGEST_FIELDS:
        delta = np.abs(candidate_values[field] - base_values[field])
        scale = np.maximum(np.abs(base_values[field]), 1e-300)
        field_deltas[field] = {
            'max_abs': float(np.max(delta)),
            'max_rel': float(np.max(delta / scale)),
        }
    return {
        'case': name,
        'baseline': baseline,
        'candidate': candidate,
        'digest_equal': baseline['digests'] == candidate['digests'],
        'field_deltas': field_deltas,
        'dn_abs_delta': abs(candidate['DN_gw'] - baseline['DN_gw']),
        'baseline_median_ms': base_median * 1e3,
        'candidate_median_ms': candidate_median * 1e3,
        'baseline_p95_ms': float(np.percentile(base_times, 95) * 1e3),
        'candidate_p95_ms': float(np.percentile(candidate_times, 95) * 1e3),
        'speed_ratio': candidate_median / base_median,
    }


def main():
    repeats = int(os.environ.get('OUTER_GRID_REUSE_REPEATS', '25'))
    payload = {
        'experiment': 'outer_goal_grid_reuse_ab',
        'production_path_changed': False,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'repeats': repeats,
        'points': POINTS,
        'rows': [run_case(name, repeats) for name in POINTS],
        'semantics': 'Reuse first outer goal grid only; all other production inputs unchanged.',
    }
    output = ROOT / 'docs' / 'outer_goal_grid_reuse_round_20260917.json'
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n',
                      encoding='utf-8')
    print(json.dumps({
        'commit': payload['generated_commit'],
        'digest_mismatch': sum(not row['digest_equal'] for row in payload['rows']),
        'speed_ratios': {row['case']: row['speed_ratio'] for row in payload['rows']},
        'out': str(output),
    }, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
