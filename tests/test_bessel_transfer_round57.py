"""Contract tests for the standalone variable-coefficient Bessel map."""

import numpy as np

from scripts.benchmark_bessel_transfer_round57 import (
    _fundamental,
    _ode_transfer,
    bessel_transfer,
)


def test_bessel_fundamental_is_finite_on_representative_q():
    matrix = _fundamental(1.25, 0.5, 0.01)
    assert matrix.shape == (2, 2)
    assert np.isfinite(matrix).all()
    assert abs(np.linalg.det(matrix)) > 0.0


def test_bessel_transfer_matches_independent_dop853():
    candidate = bessel_transfer(0.1, 1.0, 0.4, 1.1, 0.005)
    oracle = _ode_transfer(0.1, 1.0, 0.4, 1.1, 0.005)
    assert np.allclose(candidate, oracle, rtol=2e-9, atol=2e-11)
