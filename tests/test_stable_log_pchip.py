import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator

from scripts.stable_log_pchip import LogPchip


def test_log_pchip_matches_scipy_for_representable_positive_data():
    x = np.array([0.0, 0.01, 0.2, 0.21, 1.3, 3.0])
    y = np.array([1e-200, 3e-100, 2e-80, 1e-60, 1e-20, 1.0])
    points = np.linspace(x[0], x[-1], 301)
    reference = PchipInterpolator(x, y)
    candidate = LogPchip(x, np.log(y))

    np.testing.assert_allclose(candidate.evaluate_log(points), np.log(reference(points)), rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(candidate.integrate_log(x[0], x[-1]), np.log(reference.integrate(x[0], x[-1])), rtol=0.0, atol=2e-12)


def test_log_pchip_preserves_sign_changes_and_endpoint_shape_clipping():
    x = np.array([-1.0, -0.7, 0.0, 0.1, 3.0])
    y = np.array([1e-50, 3e-50, 2e-50, 5e-50, 1e-50])
    points = np.linspace(x[0], x[-1], 301)
    reference = PchipInterpolator(x, y)
    candidate = LogPchip(x, np.log(y))

    np.testing.assert_allclose(candidate.evaluate_log(points), np.log(reference(points)), rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(candidate.integrate_log(x[0], x[-1]), np.log(reference.integrate(x[0], x[-1])), rtol=0.0, atol=2e-12)


def test_log_pchip_keeps_finite_nodes_and_integral_beyond_linear_range():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    log_y = np.array([0.0, -500.0, -900.0, -1000.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        candidate = LogPchip(x, log_y)
        values = candidate.evaluate_log(np.array([0.5, 1.5, 2.5]))
        integral = candidate.integrate_log(x[0], x[-1])

    np.testing.assert_array_equal(candidate.evaluate_log(x), log_y)
    assert np.isfinite(values).all()
    assert np.isfinite(integral)
