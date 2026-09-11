"""验证 benchmark 资源预算的默认契约。"""

import os

from scripts._resource_budget import apply_environment, telemetry


def test_default_environment_caps_blas_and_numba(monkeypatch):
    names = (
        'NUMBA_NUM_THREADS', 'NUMBA_THREADING_LAYER', 'FAST_THREADS',
        'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
    )
    for name in names:
        monkeypatch.delenv(name, raising=False)
    apply_environment()
    assert os.environ['NUMBA_NUM_THREADS'] == '2'
    assert os.environ['FAST_THREADS'] == '2'
    assert os.environ['NUMBA_THREADING_LAYER'] == 'workqueue'
    for name in names[3:]:
        assert os.environ[name] == '1'


def test_telemetry_records_resource_contract(monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS', '1')
    payload = telemetry(workers=1, threads=2)
    assert payload['workers'] == 1
    assert payload['concurrent_processes'] == 1
    assert payload['threads_argument'] == 2
    assert payload['blas_threads']['OMP_NUM_THREADS'] == '1'
    assert payload['logical_cpu_count'] >= 1
