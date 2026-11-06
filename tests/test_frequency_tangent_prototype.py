"""Contracts for the standalone derivative of Fast's discrete transfer map."""

import numpy as np

from scripts.prototype_frequency_tangent import (
    phase_segment_tangent,
    propagate_frequency_grid,
    propagate_mode_with_tangent,
    scaled_step_tangent,
    terminal_log_spectrum_slope,
)
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.fast_sgwb import _phase_segment, scaled_step


def _production_terminal_spectrum(Nv, Phi, Phi_mid, S2, S2inv, Sv, z0, pt, fp, z_tail):
    nv = Nv.size
    arrays = [np.zeros((1, 1), dtype=np.float64) for _ in range(3)]
    FS.solve_kernel(
        Nv, Phi, Phi_mid, S2, S2inv,
        np.array([0], dtype=np.int64), np.array([z0]),
        np.array([pt]), np.exp(-Nv), np.ones(nv), np.array([fp]),
        1, 1, nv, float(Nv[1] - Nv[0]), z_tail,
        *arrays, np.diff(Nv), Sv, 0.25, np.full(1, -1.0), -1, 0.0, 0.0, 0.0,
    )
    return float(arrays[0][0, -1] - arrays[1][0, -1])


def test_scaled_step_tangent_matches_fixed_branch_finite_difference():
    h = 0.013
    state = np.array([0.37, -0.81])
    zero_tangent = np.zeros(2)
    for z_mid in (-1.0, -1e-8, 0.0, 1e-8, 0.3, 2.0):
        _, tangent = scaled_step_tangent(*state, *zero_tangent, z_mid, h)
        eps = 1e-6
        plus = np.array(scaled_step(*state, z_mid + eps, h))
        minus = np.array(scaled_step(*state, z_mid - eps, h))
        finite_difference = (plus - minus) / (2.0 * eps)
        np.testing.assert_allclose(tangent, finite_difference, rtol=2e-7, atol=2e-9)


def test_phase_segment_tangent_matches_shifted_segment_finite_difference():
    state = np.array([0.23, 0.91])
    tangent = np.array([-0.07, 0.91])
    z_start, z_end, h_step, phase_max = 0.2, 0.4, 0.2, 0.1
    result, result_tangent = phase_segment_tangent(
        *state, *tangent, z_start, z_end, h_step, phase_max
    )
    eps = 1e-6
    plus = np.array(_phase_segment(
        *(state + eps * tangent), z_start + eps, z_end + eps, h_step, phase_max
    ))
    minus = np.array(_phase_segment(
        *(state - eps * tangent), z_start - eps, z_end - eps, h_step, phase_max
    ))
    finite_difference = (plus - minus) / (2.0 * eps)
    np.testing.assert_allclose(result_tangent, finite_difference, rtol=3e-7, atol=3e-9)
    assert np.isfinite(result).all()


def test_terminal_log_spectrum_slope_includes_transfer_and_primordial_tilt():
    x, y, vx, vy, nt = 0.7, -0.2, -0.1, 0.35, 0.04
    expected = nt + 2.0 * (x * vx + y * vy) / (x * x + y * y)

    actual = terminal_log_spectrum_slope(x, y, vx, vy, nt)

    assert actual == expected


def test_full_mode_tangent_matches_production_exact_endpoint_spectrum():
    Nv = np.linspace(0.0, 0.2, 21)
    Phi = 0.5 * Nv
    Phi_mid = np.zeros_like(Nv)
    Phi_mid[:-1] = 0.5 * (Phi[:-1] + Phi[1:])
    S2 = np.ones_like(Nv)
    S2inv = np.ones_like(Nv)
    sigma = np.ones_like(Nv)
    ev_minus = np.exp(-Nv)
    fp_minus = np.ones_like(Nv)
    h_arr = np.diff(Nv)
    z0, nt, freq = 0.4, 0.07, 2.0
    pt = freq**nt

    _, slope, *_ = propagate_mode_with_tangent(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr,
        0, z0, pt, freq, ev_minus, fp_minus, 0.25, np.inf, -1, 0.0, 0.0, nt,
    )
    eps = 1e-6
    center = _production_terminal_spectrum(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, z0, pt, freq, np.inf
    )
    plus = _production_terminal_spectrum(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, z0 + eps,
        pt * np.exp(nt * eps), freq * np.exp(eps), np.inf,
    )
    minus = _production_terminal_spectrum(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, z0 - eps,
        pt * np.exp(-nt * eps), freq * np.exp(-eps), np.inf,
    )
    finite_difference = (np.log(plus) - np.log(minus)) / (2.0 * eps)

    assert center > 0.0
    np.testing.assert_allclose(slope, finite_difference, rtol=2e-6, atol=2e-7)


def test_full_mode_tangent_differentiates_analytic_tail_assembly():
    Nv = np.linspace(0.0, 1.0, 21)
    Phi = 0.5 * Nv
    Phi_mid = np.zeros_like(Nv)
    Phi_mid[:-1] = 0.5 * (Phi[:-1] + Phi[1:])
    S2 = np.ones_like(Nv)
    S2inv = np.ones_like(Nv)
    sigma = np.full_like(Nv, 2.0 / 3.0)
    ev_minus = np.exp(-Nv)
    fp_minus = np.ones_like(Nv)
    h_arr = np.diff(Nv)
    z0, nt, freq = 4.83, -0.03, 1.7
    pt = freq**nt

    _, slope, *_ = propagate_mode_with_tangent(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr,
        0, z0, pt, freq, ev_minus, fp_minus, 0.25, 5.0, -1, 0.0, 0.0, nt,
    )
    eps = 1e-6
    plus = _production_terminal_spectrum(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, z0 + eps,
        pt * np.exp(nt * eps), freq * np.exp(eps), 5.0,
    )
    minus = _production_terminal_spectrum(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, z0 - eps,
        pt * np.exp(-nt * eps), freq * np.exp(-eps), 5.0,
    )
    finite_difference = (np.log(plus) - np.log(minus)) / (2.0 * eps)

    assert plus > 0.0 and minus > 0.0
    np.testing.assert_allclose(slope, finite_difference, rtol=2e-5, atol=2e-6)


def test_parallel_frequency_grid_tangent_matches_independent_mode_calls():
    Nv = np.linspace(0.0, 0.04, 5)
    Phi = 0.5 * Nv
    Phi_mid = np.zeros_like(Nv)
    Phi_mid[:-1] = 0.5 * (Phi[:-1] + Phi[1:])
    S2 = np.ones_like(Nv)
    S2inv = np.ones_like(Nv)
    sigma = np.ones_like(Nv)
    h_arr = np.diff(Nv)
    z0 = np.array([0.2, 0.4])
    pt = np.array([1.0, 0.9])
    fp = np.array([1.3, 1.7])

    spectra, slopes = propagate_frequency_grid(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr,
        np.array([0, 0], dtype=np.int64), z0, pt, fp,
        np.exp(-Nv), np.ones_like(Nv), 0.25, np.inf, -1, 0.0, 0.0,
        np.array([0.0, 0.1]),
    )
    independent = [
        propagate_mode_with_tangent(
            Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr,
            0, z0[i], pt[i], fp[i], np.exp(-Nv), np.ones_like(Nv),
            0.25, np.inf, -1, 0.0, 0.0, np.array([0.0, 0.1])[i],
        )
        for i in range(2)
    ]

    np.testing.assert_array_equal(spectra, np.array([row[0] for row in independent]))
    np.testing.assert_array_equal(slopes, np.array([row[1] for row in independent]))
