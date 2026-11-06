"""Regression tests for the public frequency-quadrature method semantics."""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from stiffgwpy_fast import fast_sgwb as FS


def test_natural_cubic_quadrature_preserves_legacy_not_a_knot_boundary():
    freqs = np.array([0.0, 0.2, 0.5, 1.1, 1.7])
    integrand = np.array([1.0, 1.3, 0.9, 1.5, 1.2])
    spline = CubicSpline(freqs, integrand)
    expected = float(spline.integrate(freqs[0], freqs[-1]))

    actual = FS.integrate_frequency_quadrature(
        freqs, integrand, method="natural_cubic"
    )

    assert actual == pytest.approx(expected, rel=1e-13, abs=1e-13)


def test_natural_cubic_local_estimator_uses_legacy_interval_integrals():
    freqs = np.array([0.0, 0.2, 0.5, 1.1, 1.7])
    integrand = np.array([1.0, 1.3, 0.9, 1.5, 1.2])
    spline = CubicSpline(freqs, integrand)
    expected = np.array([
        spline.integrate(left, right)
        for left, right in zip(freqs[:-1], freqs[1:])
    ])

    _error, candidate, _baseline = FS.estimate_frequency_quadrature_local(
        freqs, integrand, method="natural_cubic"
    )

    np.testing.assert_allclose(candidate, expected, rtol=1e-13, atol=1e-13)
