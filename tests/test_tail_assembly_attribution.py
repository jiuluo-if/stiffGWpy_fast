"""Contract tests for the tail-assembly Amdahl screen."""

import numpy as np

from scripts.benchmark_tail_assembly_attribution import _tail_factors


def test_tail_factor_screen_preserves_shape_and_finiteness():
    ev_minus = np.exp(-np.arange(5, dtype=np.float64))
    fp_minus = np.linspace(1.0, 2.0, 5)
    factors = _tail_factors(ev_minus, fp_minus)
    assert factors.shape == ev_minus.shape
    assert np.isfinite(factors).all()
    np.testing.assert_array_equal(factors, ev_minus * fp_minus)
