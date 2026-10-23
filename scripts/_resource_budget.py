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
    """在导入 Numba 前设置保守线程环境；导入后不再改写进程配置。"""
    if 'numba' in sys.modules:
        return
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


def _affinity_errors() -> tuple[type[BaseException], ...]:
    """Return optional psutil affinity errors alongside built-in errors.

    中文：合并 psutil 可选异常与 Python 内置异常，统一处理 affinity 的查询和设置。
    """
    try:
        import psutil
    except ImportError:
        return AttributeError, OSError
    return AttributeError, OSError, psutil.Error


def limit_affinity(process, requested: int | None = None) -> list[int]:
    """Best-effort affinity restriction; return selected CPU IDs or an empty list.

    中文：尽力限制当前进程的 CPU affinity；设置失败时返回空列表，实际 affinity
    可由 ``telemetry()`` 单独读取。
    """
    affinity_errors = _affinity_errors()
    try:
        available = list(process.cpu_affinity())
    except affinity_errors:
        return []
    count = len(available) if requested is None else max(1, int(requested))
    selected = available[:min(count, len(available))]
    try:
        process.cpu_affinity(selected)
    except affinity_errors:
        return []
    return selected


def telemetry(process=None, workers: int = 1, threads: int | None = None,
              concurrent_processes: int = 1) -> dict:
    """Return comparable resource metadata without selecting a Numba backend.

    中文：收集可比较的资源、版本和提交信息；Numba 尚未初始化时只记录配置，不触发后端选择。
    """
    try:
        import numpy as np
    except ImportError:
        np = None
    try:
        import numba
        numba_version = numba.__version__
    except Exception:
        numba_threads = threads
        numba_threads_configured = None
        threading_layer = None
        threading_layer_requested = os.environ.get('NUMBA_THREADING_LAYER')
        numba_version = None
    else:
        numba_threads_configured = int(numba.config.NUMBA_NUM_THREADS)
        threading_layer_requested = (
            os.environ.get('NUMBA_THREADING_LAYER') or numba.config.THREADING_LAYER)
        try:
            threading_layer = numba.threading_layer()
        except ValueError:
            numba_threads = None
            threading_layer = None
        else:
            try:
                numba_threads = int(numba.get_num_threads())
            except Exception:
                numba_threads = None
    affinity = []
    if process is not None:
        affinity_errors = _affinity_errors()
        try:
            affinity = list(process.cpu_affinity())
        except affinity_errors:
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
        'numba_threads_configured': numba_threads_configured,
        'numba_threading_layer': threading_layer,
        'numba_threading_layer_requested': threading_layer_requested,
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
