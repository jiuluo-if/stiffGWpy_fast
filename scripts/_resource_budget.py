"""统一的低压力 benchmark 资源预算与运行环境记录。"""

from __future__ import annotations

import os
import platform
import sys


def nested_thread_budget(workers: int = 1, threads: int = 2) -> int:
    """Return safe inner threads for a possibly parallel outer sweep."""
    workers = int(workers)
    threads = int(threads)
    if workers < 1 or threads < 1:
        raise ValueError('workers and threads must be positive integers')
    return 1 if workers > 1 else threads


def apply_environment(numba_threads: int = 2, fast_threads: int | None = None) -> None:
    """在导入数值库前设置保守的线程环境。"""
    fast_value = numba_threads if fast_threads is None else fast_threads
    defaults = {
        'NUMBA_NUM_THREADS': str(numba_threads),
        'NUMBA_THREADING_LAYER': 'workqueue',
        'FAST_THREADS': str(fast_value),
        'OMP_NUM_THREADS': '1',
        'OPENBLAS_NUM_THREADS': '1',
        'MKL_NUM_THREADS': '1',
        'VECLIB_MAXIMUM_THREADS': '1',
        'NUMEXPR_NUM_THREADS': '1',
    }
    for name, value in defaults.items():
        os.environ.setdefault(name, value)


def limit_affinity(process, requested: int | None = None) -> list[int]:
    """限制当前进程 affinity；不可用时返回空列表而不改变求解。"""
    try:
        available = list(process.cpu_affinity())
    except (AttributeError, OSError):
        return []
    count = len(available) if requested is None else max(1, int(requested))
    selected = available[:min(count, len(available))]
    try:
        process.cpu_affinity(selected)
    except (AttributeError, OSError):
        return available
    return selected


def telemetry(process=None, workers: int = 1, threads: int | None = None,
              concurrent_processes: int = 1) -> dict:
    """返回可比较的资源、版本和提交信息。"""
    try:
        import numpy as np
    except ImportError:
        np = None
    try:
        import numba
        numba_threads = int(numba.get_num_threads())
        threading_layer = numba.threading_layer()
        numba_version = numba.__version__
    except Exception:
        numba_threads = threads
        threading_layer = os.environ.get('NUMBA_THREADING_LAYER')
        numba_version = None
    affinity = []
    if process is not None:
        try:
            affinity = list(process.cpu_affinity())
        except (AttributeError, OSError):
            pass
    try:
        import scipy
        scipy_version = scipy.__version__
    except ImportError:
        scipy_version = None
    payload = {
        'logical_cpu_count': os.cpu_count(),
        'affinity': affinity,
        'numba_threads': numba_threads,
        'numba_threading_layer': threading_layer,
        'blas_threads': {
            name: os.environ.get(name)
            for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                         'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS',
                         'NUMEXPR_NUM_THREADS')
        },
        'workers': int(workers),
        'concurrent_processes': int(concurrent_processes),
        'python': sys.version.split()[0],
        'platform': platform.platform(),
        'numpy': None if np is None else np.__version__,
        'scipy': scipy_version,
        'numba': numba_version,
        'threads_argument': threads,
    }
    return payload
