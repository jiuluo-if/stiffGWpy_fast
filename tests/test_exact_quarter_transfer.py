"""Tests for exact-background quarter-sample tensor transfers."""

import numpy as np

from stiffgwpy_fast.fast_sgwb import _phase_segment, _phase_segment_quarters


def test_quarter_transfer_matches_four_const_z_segments():
    state = (0.25, -0.4)
    z = (-0.3, 0.05, 0.4, 0.85, 1.2)
    width = 0.08
    actual = _phase_segment_quarters(*state, *z, width, 0.0)

    expected = state
    for left, right in zip(z[:-1], z[1:]):
        expected = _phase_segment(*expected, left, right, 0.25 * width, 0.0)

    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2e-15)
    coarse = _phase_segment(*state, z[0], z[-1], width, 0.0)
    assert not np.allclose(actual, coarse, rtol=1e-10, atol=1e-12)
