# -*- coding: utf-8 -*-
"""对 fast/plain-grid 做可重复的阶段级计时，不改变求解器行为。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time
from collections import defaultdict

import numpy as np
from numba import set_parallel_chunksize

ROOT = os.environ.get('STIFFGWPY_ROOT',
                      os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'A': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'B': dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
}
CASE = CASES['A']


def p95(values):
    return float(np.percentile(np.asarray(values), 95))


def timed_call(totals, name, fn, *args, **kwargs):
    start = time.perf_counter()
    result = fn(*args, **kwargs)
    totals[name].append(time.perf_counter() - start)
    return result


def run_once(case):
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
    def profiled_solve(*args, **kwargs):
        j0s = np.asarray(args[5])
        totals['j0_cv'].append(float(np.std(j0s) / max(np.mean(j0s), 1.0)))
        totals['j0_span'].append(float(np.max(j0s) - np.min(j0s)))
        solve_args = args
        if os.environ.get('PROFILE_ASSEMBLE') == '0':
            # 仅用于估算 assembly 的理论上限；该模式的完整输出不作为
            # 数值结果，正式 benchmark 必须保持默认 assemble=1。
            solve_args = list(args)
            solve_args[11] = 0
            solve_args = tuple(solve_args)
        return timed_call(totals, 'tensor_solve_kernel', original['solve_kernel'], *solve_args, **kwargs)
    FS.solve_kernel = profiled_solve
    FS.int_SGWB_W = lambda *a, **k: timed_call(totals, 'column_integration', original['int_SGWB_W'], *a, **k)
    LCDM_SG.construct_f = lambda *a, **k: timed_call(totals, 'frequency_grid', original['construct_f'], *a, **k)
    try:
        model = LCDM_SG(**case)
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
    parser.add_argument('--case', choices=sorted(CASES), default='A')
    args = parser.parse_args()
    FS.apply_accuracy_mode('fast')
    # 保留基准命令显式指定的线程数，避免 preset 默认值遮蔽 scaling 点。
    if os.environ.get('FAST_THREADS'):
        FS.set_threads(int(os.environ['FAST_THREADS']))
    if os.environ.get('FAST_CHUNKSIZE'):
        set_parallel_chunksize(int(os.environ['FAST_CHUNKSIZE']))
    records = []
    for _ in range(args.reps):
        _, totals, model = run_once(CASES[args.case])
        def digest(name):
            arr = np.ascontiguousarray(np.asarray(getattr(model, name), dtype=np.float64))
            return hashlib.sha256(arr.tobytes()).hexdigest()
        row = {'total_s': totals['total'][0],
               'iterations': len(totals['tensor_solve_kernel']),
               'n_freq': len(model.f),
               'DN_eff': float(model.cosmo_param['DN_eff']),
               'DN_gw_last': float(np.asarray(model.DN_gw)[-1]),
               'digest_f': digest('f'),
               'digest_spectrum': digest('log10OmegaGW'),
               'digest_DN_gw': digest('DN_gw'),
               'digest_g2': digest('g2'),
               'digest_w2': digest('w2')}
        for name, values in totals.items():
            if name in ('total', 'outer_iteration'):
                continue
            if name == 'j0_cv':
                row[name] = float(values[-1])
                continue
            if name == 'j0_span':
                row[name] = float(values[-1])
                continue
            row[name + '_s'] = float(sum(values))
            row[name + '_calls'] = len(values)
        records.append(row)
    summary = {'threads': FS._THREADS, 'case': CASES[args.case], 'case_id': args.case,
               'reps': args.reps, 'records': records,
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
