"""TDD contract for the standalone tail-factor cache spike."""

import numpy as np

from scripts.benchmark_tail_factor_cache import tail_factor


def test_tail_factor_preserves_elementwise_product():
    ev_minus = np.array([1.0, 0.5, 0.125], dtype=np.float64)
    fp_minus = np.array([2.0, 3.0, 5.0], dtype=np.float64)
    expected = ev_minus * fp_minus

    actual = tail_factor(ev_minus, fp_minus)

    np.testing.assert_array_equal(actual, expected)
