# -*- coding: utf-8 -*-
"""A/B test the local derived-parameter cache in gen_expansion."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import LCDM_SN  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def _digest(value):
    array = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    return hashlib.sha256(array.tobytes()).hexdigest()


def _old_gen_expansion():
    source = subprocess.check_output(
        ['git', 'show', 'HEAD:stiffgwpy_fast/LCDM_stiff_Neff.py'], text=True,
    )
    start = source.index('    def gen_expansion(')
    end = source.index('\ndef input_nt(', start)
    method = textwrap.dedent(source[start:end]).replace(
        'def gen_expansion(', 'def old_gen_expansion(', 1)
    namespace = vars(__import__('stiffgwpy_fast.LCDM_stiff_Neff', fromlist=['*']))
    exec(method, namespace)
    return namespace['old_gen_expansion']


def _run(case, old_method, new_method):
    LCDM_SN.gen_expansion = old_method if old_method else new_method
    model = LCDM_SG(**CASES[case])
    start = time.perf_counter()
    result = FS.SGWB_iter_fast(model)
    elapsed = time.perf_counter() - start
    return {
        'elapsed_s': elapsed,
        'digest_Nv': _digest(model.Nv),
        'digest_sigma': _digest(model.sigma),
        'digest_f_hor': _digest(model.f_hor),
        'digest_spectrum': _digest(model.log10OmegaGW),
        'digest_DN_gw': _digest(model.DN_gw),
        'DN_gw': float(np.asarray(model.DN_gw)[-1]),
        'n_freq': len(model.f),
        'result_type': type(result).__name__,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', type=Path, default=None)
    args = parser.parse_args()
    reps = 25
    old_method = _old_gen_expansion()
    new_method = LCDM_SN.gen_expansion
    FS.apply_accuracy_mode('fast')
    FS.set_threads(int(os.environ.get('FAST_THREADS', '2')))
    rows = []
    for case in CASES:
        # 先预热两个实现，再采集计时。
        _run(case, old_method, new_method)
        _run(case, None, new_method)
        baseline = [_run(case, old_method, new_method) for _ in range(reps)]
        candidate = [_run(case, None, new_method) for _ in range(reps)]
        rows.append({
            'case': case,
            'reps': reps,
            'baseline_median_s': statistics.median(x['elapsed_s'] for x in baseline),
            'candidate_median_s': statistics.median(x['elapsed_s'] for x in candidate),
            'speed_ratio_candidate_over_baseline': (
                statistics.median(x['elapsed_s'] for x in candidate)
                / statistics.median(x['elapsed_s'] for x in baseline)
            ),
            'digest_equal': all(
                baseline[0][name] == candidate[0][name]
                for name in (
                    'digest_Nv', 'digest_sigma', 'digest_f_hor',
                    'digest_spectrum', 'digest_DN_gw',
                )
            ),
            'baseline': baseline[0],
            'candidate': candidate[0],
        })
    output = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'threads': FS._THREADS,
        'rows': rows,
    }
    rendered = json.dumps(output, ensure_ascii=False, indent=2)
    print(rendered)
    if args.json is not None:
        args.json.write_text(rendered + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
