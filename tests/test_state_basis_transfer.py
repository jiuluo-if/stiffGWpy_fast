import numpy as np

from scripts.benchmark_state_basis_transfer import (
    _basis_step,
    _cartesian_step,
)


def test_state_basis_transfer_is_finite_and_reconstructs():
    x, y = 0.125, 0.875
    u, v = x + y, y - x
    bx, by = _basis_step(u, v, 1.2, 0.005)
    cx, cy = _cartesian_step(x, y, 1.2, 0.005)
    rx, ry = 0.5 * (bx - by), 0.5 * (bx + by)
    assert np.isfinite((rx, ry)).all()
    assert np.max(np.abs(np.asarray((rx, ry)) - np.asarray((cx, cy)))) < 1e-14


def test_state_basis_long_chain_is_finite():
    x, y = 0.1, 1.0
    u, v = x + y, y - x
    for _ in range(256):
        u, v = _basis_step(u, v, 0.2, 0.01)
    assert np.isfinite((u, v)).all()
