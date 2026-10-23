"""TDD contract for the Round 53 scaled Riccati prototype."""

import numpy as np


def test_riccati_constant_coefficient_step_is_finite_and_close():
    from scripts.benchmark_scaled_riccati_round53 import (
        _cartesian_from_state,
        _riccati_step,
        _state_from_cartesian,
    )
    from stiffgwpy_fast.fast_sgwb import scaled_step

    x0, y0 = 0.37, -0.81
    state = _state_from_cartesian(x0, y0)
    candidate = _cartesian_from_state(_riccati_step(state, 0.8, 0.005))
    reference = scaled_step(x0, y0, np.log(0.8), 0.005)
    assert np.all(np.isfinite(candidate))
    assert np.max(np.abs(np.asarray(candidate) - np.asarray(reference))) < 2e-4


def test_riccati_chart_switch_round_trip_is_finite():
    from scripts.benchmark_scaled_riccati_round53 import (
        _cartesian_from_state,
        _state_from_cartesian,
    )

    state = _state_from_cartesian(1.0e-12, 1.0)
    x, y = _cartesian_from_state(state)
    assert np.isfinite(x) and np.isfinite(y)
    assert abs(x - 1.0e-12) < 1e-24
    assert abs(y - 1.0) < 1e-14


def test_scaled_step_zero_frequency_boundary_has_finite_limit():
    from stiffgwpy_fast.fast_sgwb import _phase_segment, scaled_step

    x0, y0, h = 0.37, -0.81, 0.005
    expected = ((1.0 - h) * x0 - h * y0, h * x0 + (1.0 + h) * y0)

    direct = scaled_step(x0, y0, z_mid=0.0, h=h)
    segmented = _phase_segment(x0, y0, z_start=-0.1, z_end=0.1, h_step=h, phase_max=0.0)

    assert np.all(np.isfinite(direct))
    np.testing.assert_allclose(direct, expected, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(segmented, expected, rtol=0.0, atol=1e-15)
