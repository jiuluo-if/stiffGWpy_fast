# -*- coding: utf-8 -*-
"""在固定 CPU affinity 与 Numba 线程层下运行分层 profiler。"""
from __future__ import annotations

import argparse
import os
import sys

try:
    from scripts._resource_budget import apply_environment, limit_affinity, telemetry
except ImportError:
    from _resource_budget import apply_environment, limit_affinity, telemetry

apply_environment()

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--affinity-count', type=int, default=2)
    parser.add_argument('--reps', type=int, default=7)
    parser.add_argument('--case', choices=('A', 'B'), default='A')
    parser.add_argument('--kink-split', action='store_true')
    parser.add_argument('--disable-outer-reuse', action='store_true')
    parser.add_argument('--json', default=None)
    args = parser.parse_args()
    if args.threads > args.affinity_count:
        raise SystemExit('threads must not exceed affinity-count')

    try:
        import psutil
    except ImportError as exc:
        raise SystemExit('固定 affinity 需要 psutil') from exc

    process = psutil.Process()
    selected = limit_affinity(process, args.affinity_count)

    os.environ['FAST_THREADS'] = str(args.threads)
    from numba import threading_layer
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import profile_fast_breakdown as profile

    profile_args = ['--reps', str(args.reps), '--case', args.case]
    if args.kink_split:
        profile_args.append('--kink-split')
    if args.disable_outer_reuse:
        profile_args.append('--disable-outer-reuse')
    if args.json:
        profile_args.extend(['--json', args.json])
    sys.argv = ['profile_fast_breakdown.py', *profile_args]
    profile.main()
    print('PROFILE_META affinity=%s threading_layer=%s threads=%s' %
          (selected, threading_layer(), args.threads))
    print('PROFILE_RESOURCES %s' % telemetry(
        process, threads=args.threads).get('blas_threads'))


if __name__ == '__main__':
    main()
