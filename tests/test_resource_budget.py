"""验证 benchmark 资源预算的默认契约。"""

import os
import subprocess
import sys

from scripts._resource_budget import telemetry


def test_nested_budget_reduces_inner_threads():
    from scripts._resource_budget import nested_thread_budget

    assert nested_thread_budget(workers=1, threads=4) == 4
    assert nested_thread_budget(workers=2, threads=4) == 1


def test_default_environment_caps_blas_and_numba(monkeypatch):
    code = (
        'import os; '
        'from scripts._resource_budget import apply_environment; '
        'apply_environment(); '
        "names = ('NUMBA_NUM_THREADS', 'NUMBA_THREADING_LAYER', "
        "'FAST_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', "
        "'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'); "
        "print('\\n'.join(os.environ[name] for name in names))"
    )
    clean_env = os.environ.copy()
    for name in (
        'NUMBA_NUM_THREADS', 'NUMBA_THREADING_LAYER', 'FAST_THREADS',
        'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
    ):
        clean_env.pop(name, None)
    result = subprocess.run(
        [sys.executable, '-c', code],
        check=True,
        capture_output=True,
        text=True,
        env=clean_env,
    )
    assert result.stdout.splitlines() == [
        '2', 'workqueue', '2', '1', '1', '1', '1', '1'
    ]


def test_telemetry_records_resource_contract(monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    payload = telemetry(workers=1, threads=2)
    assert payload['workers'] == 1
    assert payload['concurrent_processes'] == 1
    assert payload['threads_argument'] == 2
    assert payload['blas_threads']['OMP_NUM_THREADS'] == '1'
    assert payload['logical_cpu_count'] >= 1
