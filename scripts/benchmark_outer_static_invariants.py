"""审计 outer iteration 中背景、频率网格和准备数组的跨轮不变量。"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import exact_background as EB  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def _digest(value):
    array = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    return hashlib.sha256(array.tobytes()).hexdigest()


def _delta(first, second):
    if first is None or second is None:
        return None
    a = np.asarray(first, dtype=float)
    b = np.asarray(second, dtype=float)
    if a.shape != b.shape:
        return {'shape_changed': True, 'shape_first': a.shape, 'shape_second': b.shape}
    return {'shape_changed': False,
            'max_abs': float(np.max(np.abs(a - b))) if a.size else 0.0,
            'digest_equal': _digest(a) == _digest(b)}


def run_case(name):
    original = (FS.gen_fast, FA.goal_oriented_freqs, FS.prep_frequency_only,
                EB.fast_phi_s2_split)
    observations = {'gen': [], 'goal': [], 'prep': [], 'primitive': []}

    def gen(*args, **kwargs):
        result = original[0](*args, **kwargs)
        observations['gen'].append({
            'Nv': np.asarray(args[0].Nv, dtype=float).copy(),
            'sigma': np.asarray(args[0].sigma, dtype=float).copy(),
            'f_hor': np.asarray(args[0].f_hor, dtype=float).copy(),
        })
        return result

    def goal(*args, **kwargs):
        result = original[1](*args, **kwargs)
        observations['goal'].append(np.asarray(result, dtype=float).copy())
        return result

    def prep(*args, **kwargs):
        result = original[2](*args, **kwargs)
        observations['prep'].append(tuple(
            np.asarray(value, dtype=float).copy() if isinstance(value, np.ndarray)
            else value for value in result))
        return result

    def primitive(*args, **kwargs):
        result = original[3](*args, **kwargs)
        observations['primitive'].append(tuple(
            np.asarray(value, dtype=float).copy() if isinstance(value, np.ndarray)
            else value for value in result))
        return result

    FS.gen_fast, FA.goal_oriented_freqs = gen, goal
    FS.prep_frequency_only, EB.fast_phi_s2_split = prep, primitive
    try:
        model = LCDM_SG(**CASES[name])
        FS.apply_accuracy_mode('fast')
        FS.set_threads(2)
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid='goal')
    finally:
        (FS.gen_fast, FA.goal_oriented_freqs, FS.prep_frequency_only,
         EB.fast_phi_s2_split) = original

    row = {'case': name, 'converged': bool(getattr(model, 'SGWB_converge', False)),
           'failure': getattr(model, 'fast_failure_reason', None),
           'calls': {key: len(value) for key, value in observations.items()}}
    for key, values in observations.items():
        if not values:
            continue
        if key == 'gen':
            fields = ('Nv', 'sigma', 'f_hor')
            row[key] = {field: _delta(values[0][field], values[-1][field])
                         for field in fields}
        elif key in ('goal',):
            row[key] = _delta(values[0], values[-1])
        else:
            row[key] = [_delta(values[0][index], values[-1][index])
                        for index in range(len(values[0]))]
    return row


def main():
    names = sys.argv[1:] or list(CASES)
    payload = {
        'experiment': 'outer_static_invariant_probe',
        'production_path_changed': False,
        'generated_commit': os.popen('git rev-parse HEAD').read().strip(),
        'resources': {'numba_threads': 2, 'blas_threads': 1, 'workers': 1},
        'rows': [run_case(name) for name in names],
        'semantics': 'Compare first and last captured outer inputs; diagnostic only.',
    }
    output = ROOT / 'docs' / 'outer_static_invariants_round_20260917.json'
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n',
                      encoding='utf-8')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
