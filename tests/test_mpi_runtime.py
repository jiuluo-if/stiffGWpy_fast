import builtins
import sys
from types import ModuleType, SimpleNamespace

import pytest

from stiffgwpy_fast.reference import _sync_size


def _install_mpi(monkeypatch, mpi):
    module = ModuleType("mpi4py")
    module.MPI = mpi
    monkeypatch.setitem(sys.modules, "mpi4py", module)


def test_sync_size_defaults_to_one_when_mpi4py_is_missing(monkeypatch):
    real_import = builtins.__import__

    def import_without_mpi(name, *args, **kwargs):
        if name == "mpi4py":
            raise ModuleNotFoundError("mpi4py is not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_mpi)

    assert _sync_size() == 1


@pytest.mark.parametrize(
    ("initialized", "finalized"),
    [(False, False), (True, True)],
)
def test_sync_size_defaults_to_one_outside_active_mpi(
    monkeypatch, initialized, finalized
):
    comm = SimpleNamespace(Get_size=lambda: pytest.fail("inactive MPI was queried"))
    mpi = SimpleNamespace(
        Is_initialized=lambda: initialized,
        Is_finalized=lambda: finalized,
        COMM_WORLD=comm,
    )
    _install_mpi(monkeypatch, mpi)

    assert _sync_size() == 1


def test_sync_size_returns_active_world_size(monkeypatch):
    mpi = SimpleNamespace(
        Is_initialized=lambda: True,
        Is_finalized=lambda: False,
        COMM_WORLD=SimpleNamespace(Get_size=lambda: 4),
    )
    _install_mpi(monkeypatch, mpi)

    assert _sync_size() == 4


def test_sync_size_propagates_errors_from_active_communicator(monkeypatch):
    def fail_get_size():
        raise RuntimeError("MPI communicator failure")

    mpi = SimpleNamespace(
        Is_initialized=lambda: True,
        Is_finalized=lambda: False,
        COMM_WORLD=SimpleNamespace(Get_size=fail_get_size),
    )
    _install_mpi(monkeypatch, mpi)

    with pytest.raises(RuntimeError, match="MPI communicator failure"):
        _sync_size()
