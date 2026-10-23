import numpy as np

from stiffgwpy_fast.cobaya.likelihoods.PTA._spectrum import primordial_rho


def test_primordial_rho_returns_zero_when_no_signal_band_knots_exist():
    work_freqs = np.array([1e-8, 1e-7], dtype=np.float64)

    rho = primordial_rho([-4.0, -3.0], [-9.0, -8.0], work_freqs, 0.67, 1e8)

    np.testing.assert_array_equal(rho, np.zeros_like(work_freqs))


def test_primordial_rho_does_not_extrapolate_outside_theory_support():
    work_freqs = np.array([1e-8, 1e-7, 1e-6, 1e-5], dtype=np.float64)

    rho = primordial_rho([-7.0, -6.0], [-9.0, -8.0], work_freqs, 0.67, 1e8)

    assert rho[0] == 0.0
    assert rho[-1] == 0.0
    assert np.all(rho[1:-1] > 0.0)


def test_primordial_rho_single_knot_only_supports_an_exact_query():
    work_freqs = np.array([1e-8, 2e-8], dtype=np.float64)

    rho = primordial_rho([-8.0], [-9.0], work_freqs, 0.67, 1e8)

    expected = 10.0**-9 * 0.67**2 / (8.0 * np.pi**4 * work_freqs[0] ** 5 * 1e8)
    np.testing.assert_allclose(rho, [expected, 0.0], rtol=1e-14, atol=0.0)
