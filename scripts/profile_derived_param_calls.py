"""Profile repeated derived-parameter property evaluation.

The report separates property-access time from the complete fast solve.

中文：统计 fast 求解期间 `derived_param` 的访问次数和耗时，并与总求解时间并列输出。
本脚本用于定位开销，不改变正式求解器实现。
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.LCDM_stiff_Neff import LCDM_SN  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'stiff': dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    'high_kappa': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0,
                    T_re=1e3, kappa10=1e-3),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=sorted(CASES), default='default')
    parser.add_argument('--repeats', type=int, default=25)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be at least 1')

    original_property = LCDM_SN.derived_param
    original = original_property.fget
    stats = {'calls': 0, 'seconds': 0.0}

    def wrapped(self):
        start = time.perf_counter()
        result = original(self)
        stats['calls'] += 1
        stats['seconds'] += time.perf_counter() - start
        return result

    LCDM_SN.derived_param = property(wrapped)
    try:
        FS.apply_accuracy_mode('fast')
        FS.set_threads(2)
        timings = []
        call_samples = []
        prop_samples = []
        for _ in range(args.repeats):
            before = stats['calls']
            before_seconds = stats['seconds']
            model = LCDM_SG(**CASES[args.case])
            start = time.perf_counter()
            FS.SGWB_iter_fast(model, kink_split=True)
            elapsed = time.perf_counter() - start
            timings.append(elapsed)
            call_samples.append(stats['calls'] - before)
            prop_samples.append(stats['seconds'] - before_seconds)
    finally:
        LCDM_SN.derived_param = original_property

    def median(values):
        return statistics.median(values)

    payload = {
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'experiment': 'derived_param_call_attribution',
        'case': args.case,
        'repeats': args.repeats,
        'resources': {
            'numba_threads': 2,
            'workers': 1,
            'thread_environment': {
                name: os.environ.get(name)
                for name in ('NUMBA_NUM_THREADS', 'NUMBA_THREADING_LAYER',
                             'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                             'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
            },
        },
        'median_total_ms': 1000.0 * median(timings),
        'median_property_ms': 1000.0 * median(prop_samples),
        'median_calls': median(call_samples),
        'all_call_counts': sorted(set(call_samples)),
    }
    with open(args.out, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write('\n')
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
