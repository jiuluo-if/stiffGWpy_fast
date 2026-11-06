"""FastSolverConfig controls TBB chunk size without leaking thread-local state."""

import numpy as np
import pytest
from numba import get_parallel_chunksize, set_parallel_chunksize

from stiffgwpy_fast import LCDM_SG
from stiffgwpy_fast import fast_sgwb as FS


def test_parallel_chunksize_is_applied_and_restored_after_failed_solve(monkeypatch):
    previous = set_parallel_chunksize(3)
    observed = []

    def nonfinite_solve(
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
        fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
        Ogw, Oj, Opgw, h_arr=None, Sv=None, phase_max=0.0,
        handoff_eps=None, kink_index=-1, kink_fraction=0.0,
        phi_re=0.0, h_final=0.0,
    ):
        observed.append(get_parallel_chunksize())
        Ogw[...] = np.nan
        Oj[...] = 0.0
        Opgw[...] = 0.0

    monkeypatch.setattr(FS, "solve_kernel", nonfinite_solve)
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    config = FS.resolve_config("fast", threads=1, parallel_chunksize=1)
    try:
        assert FS.SGWB_iter_fast(model, config=config) is None
        assert observed and all(value == 1 for value in observed)
        assert get_parallel_chunksize() == 3
    finally:
        set_parallel_chunksize(previous)


def test_parallel_chunksize_is_restored_when_solve_raises(monkeypatch):
    previous = set_parallel_chunksize(3)
    observed = []

    def raising_solve(*args, **kwargs):
        observed.append(get_parallel_chunksize())
        raise RuntimeError("injected solve failure")

    monkeypatch.setattr(FS, "solve_kernel", raising_solve)
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    config = FS.resolve_config("fast", threads=1, parallel_chunksize=1)
    try:
        with pytest.raises(RuntimeError, match="injected solve failure"):
            FS.SGWB_iter_fast(model, config=config)
        assert observed == [1]
        assert get_parallel_chunksize() == 3
    finally:
        set_parallel_chunksize(previous)
