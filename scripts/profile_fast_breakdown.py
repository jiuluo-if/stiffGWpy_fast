# -*- coding: utf-8 -*-
"""对 fast/plain-grid 做可重复的阶段级计时，不改变求解器行为。"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time
from collections import defaultdict

import numpy as np

ROOT = os.environ.get('STIFFGWPY_ROOT',
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from stiffgwpy import fast_sgwb as FS  # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG  # noqa: E402

CASE = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)


def p95(values):
    return float(np.percentile(np.asarray(values), 95))


def timed_call(totals, name, fn, *args, **kwargs):
    start = time.perf_counter()
    result = fn(*args, **kwargs)
    totals[name].append(time.perf_counter() - start)
    return result


def run_once():
    totals = defaultdict(list)
    original = {
        'gen_fast': FS.gen_fast,
        'prep_fast': FS.prep_fast,
        'solve_kernel': FS.solve_kernel,
        'int_SGWB_W': FS.int_SGWB_W,
        'construct_f': LCDM_SG.construct_f,
    }
    FS.gen_fast = lambda *a, **k: timed_call(totals, 'expansion_background', original['gen_fast'], *a, **k)
    FS.prep_fast = lambda *a, **k: timed_call(totals, 'kernel_prepare', original['prep_fast'], *a, **k)
    FS.solve_kernel = lambda *a, **k: timed_call(totals, 'tensor_solve_kernel', original['solve_kernel'], *a, **k)
    FS.int_SGWB_W = lambda *a, **k: timed_call(totals, 'column_integration', original['int_SGWB_W'], *a, **k)
    LCDM_SG.construct_f = lambda *a, **k: timed_call(totals, 'frequency_grid', original['construct_f'], *a, **k)
    try:
        model = LCDM_SG(**CASE)
        start = time.perf_counter()
        result = FS.SGWB_iter_fast(model)
        totals['total'].append(time.perf_counter() - start)
        totals['outer_iteration'].append(len(getattr(model, 'DN_gw', [])))
        return result, totals, model
    finally:
        FS.gen_fast = original['gen_fast']
        FS.prep_fast = original['prep_fast']
        FS.solve_kernel = original['solve_kernel']
        FS.int_SGWB_W = original['int_SGWB_W']
        LCDM_SG.construct_f = original['construct_f']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reps', type=int, default=5)
    parser.add_argument('--json', default=None)
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    # 保留基准命令显式指定的线程数，避免 preset 默认值遮蔽 scaling 点。
    if os.environ.get('FAST_THREADS'):
        FS.set_threads(int(os.environ['FAST_THREADS']))
    records = []
    for _ in range(args.reps):
        _, totals, model = run_once()
        row = {'total_s': totals['total'][0],
               'iterations': len(totals['tensor_solve_kernel'])}
        for name, values in totals.items():
            if name in ('total', 'outer_iteration'):
                continue
            row[name + '_s'] = float(sum(values))
            row[name + '_calls'] = len(values)
        records.append(row)
    summary = {'threads': FS._THREADS, 'case': CASE, 'reps': args.reps, 'records': records,
               'median_s': {}, 'p95_s': {}}
    names = sorted(k for k in records[0] if k.endswith('_s'))
    for name in names:
        vals = [r[name] for r in records]
        summary['median_s'][name] = statistics.median(vals)
        summary['p95_s'][name] = p95(vals)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as handle:
            json.dump(summary, handle, ensure_ascii=False, indent=2)


if __name__ == '__main__':
    main()
