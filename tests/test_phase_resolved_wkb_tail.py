import math

import numpy as np

from stiffgwpy_fast.fast_sgwb import (
    _blend_phase_tail_moments,
    _phase_after_tail_entry,
    _phase_integral_suffix,
    _phase_resolved_tail_weight,
    _phase_resolved_tail_weights,
    _phase_resolved_wkb_tail_state,
    _phase_tail_effective_weight,
    _phase_tail_state_for_weights,
    _requires_phase_resolved_tail,
    _terminal_base_signal,
    solve_kernel,
)


def test_phase_tail_blends_all_terminal_moments_with_one_weight():
    average = (4.0, -2.0, 3.0)
    phase = (12.0, 6.0, -5.0)
    weight = 0.25

    actual = _blend_phase_tail_moments(*average, *phase, weight)

    np.testing.assert_allclose(actual, (6.0, 0.0, 1.0), rtol=0.0, atol=0.0)


def test_phase_resolved_tail_reconstructs_exact_radiation_mode():
    """For radiation, h=e^-N sin(e^N) is an exact oscillatory solution."""
    n_h = 2.0
    omega_h = 150.0
    omega_end = 165.0
    delta_n = math.log(omega_end / omega_h)
    delta_phase = omega_end - omega_h

    h_h = math.exp(-n_h) * math.sin(omega_h)
    x_h = math.exp(-n_h) * (-math.sin(omega_h) + omega_h * math.cos(omega_h))
    y_h = omega_h * h_h

    x_end, y_end = _phase_resolved_wkb_tail_state(
        x_h,
        y_h,
        math.log(omega_h),
        math.log(omega_end),
        delta_n,
        delta_phase,
    )

    n_end = n_h + delta_n
    h_expected = math.exp(-n_end) * math.sin(omega_end)
    x_expected = math.exp(-n_end) * (
        -math.sin(omega_end) + omega_end * math.cos(omega_end)
    )
    y_expected = omega_end * h_expected

    np.testing.assert_allclose(x_end, x_expected, rtol=2e-14, atol=2e-15)
    np.testing.assert_allclose(y_end, y_expected, rtol=2e-14, atol=2e-15)


def test_tail_selector_resolves_only_when_ripple_bound_exceeds_spectrum_budget():
    assert _requires_phase_resolved_tail(math.log(165.0))
    assert _requires_phase_resolved_tail(math.log(434.0))
    assert not _requires_phase_resolved_tail(math.log(435.0))


def test_phase_tail_weight_tapers_continuously_to_phase_average():
    relative_bound = math.expm1(math.log(10.0) * 5e-4)
    inv_q_cut = 2.0 * relative_bound / (1.0 + math.sqrt(1.0 + 4.0 * relative_bound))
    q_cut = 1.0 / inv_q_cut
    z_below = math.log(q_cut * (1.0 - 1e-9))
    z_above = math.log(q_cut * (1.0 + 1e-9))

    w_below = _phase_resolved_tail_weight(z_below)
    w_above = _phase_resolved_tail_weight(z_above)

    assert 0.0 < w_below < 1e-6
    assert w_above == 0.0
    assert 0.0 < _phase_resolved_tail_weight(math.log(165.0)) < 1.0


def test_phase_tail_weight_is_continuous_at_analytic_tail_entry():
    z_tail = 5.0
    leading = _phase_resolved_tail_weights(
        np.array([z_tail + 1e-8]), z_tail, enabled=True
    )[0]
    below = _phase_tail_effective_weight(leading, 0.0)
    above = _phase_tail_effective_weight(leading, 1e-8)

    assert below == 0.0
    assert 0.0 < above < 1e-6


def test_tail_phase_after_boundary_vanishes_at_terminal_entry():
    nv = np.array([0.0, 1.0])
    phi = np.array([0.0, 1.0])
    phase_suffix = _phase_integral_suffix(nv, phi)
    z_tail = 1.0 - 1e-9

    theta = _phase_after_tail_entry(
        nv, phi, 0, 0.0, z_tail, 0, 0, phase_suffix
    )

    assert 0.0 < theta < 3e-9


def test_analytic_tail_recovers_positive_energy_before_cancellation():
    ogw = np.array([10.0, 1.0])
    oj = np.array([2.0, 1.0])
    opgw = np.array([3.0, 0.1])
    handoff_eps = np.array([-1.0, 0.0])

    actual = _terminal_base_signal(ogw, oj, opgw, handoff_eps)

    np.testing.assert_allclose(actual, np.array([8.0, 0.3]), rtol=0.0, atol=1e-16)


def test_first_probe_writes_terminal_main_sample_on_coarse_boundary():
    nv = np.arange(5, dtype=np.float64) * 0.01
    phi = np.arange(5, dtype=np.float64) * 0.001
    phi_mid = np.r_[0.5 * (phi[:-1] + phi[1:]), phi[-1]]
    s2 = np.ones(5, dtype=np.float64)
    s2inv = np.ones(5, dtype=np.float64)
    ogw = np.zeros((1, 3), dtype=np.float64)
    oj = np.zeros_like(ogw)
    opgw = np.zeros_like(ogw)

    solve_kernel(
        nv,
        phi,
        phi_mid,
        s2,
        s2inv,
        np.array([1], dtype=np.int64),
        np.array([-0.5], dtype=np.float64),
        np.ones(1, dtype=np.float64),
        np.exp(-nv),
        np.ones(5, dtype=np.float64),
        np.ones(1, dtype=np.float64),
        0,
        3,
        2,
        0.01,
        100.0,
        ogw,
        oj,
        opgw,
        np.diff(nv),
        None,
        0.0,
        np.full(1, -1.0),
    )

    assert ogw[0, -1] - oj[0, -1] > 0.0


def test_solve_kernel_skips_stale_tail_state_when_phase_weight_is_none():
    nv = np.arange(5, dtype=np.float64) * 0.01
    phi = np.arange(5, dtype=np.float64) * 0.001
    phi_mid = np.r_[0.5 * (phi[:-1] + phi[1:]), phi[-1]]
    s2 = np.ones(5, dtype=np.float64)
    s2inv = np.ones(5, dtype=np.float64)
    ogw = np.zeros((1, 3), dtype=np.float64)
    oj = np.zeros_like(ogw)
    opgw = np.zeros_like(ogw)
    tail_state = np.zeros((1, 4), dtype=np.float64)
    tail_state[:, 3] = -1.0
    tail_state_for_solve = _phase_tail_state_for_weights(None, tail_state)

    solve_kernel(
        nv,
        phi,
        phi_mid,
        s2,
        s2inv,
        np.array([1], dtype=np.int64),
        np.array([-0.001], dtype=np.float64),
        np.ones(1, dtype=np.float64),
        np.exp(-nv),
        np.ones(5, dtype=np.float64),
        np.ones(1, dtype=np.float64),
        0,
        3,
        2,
        0.01,
        0.0,
        ogw,
        oj,
        opgw,
        np.diff(nv),
        None,
        0.0,
        np.full(1, -1.0),
        phase_weight=None,
        tail_state=tail_state_for_solve,
    )

    np.testing.assert_array_equal(tail_state[:, 3], np.array([-1.0]))


def test_phase_tail_state_is_not_passed_without_phase_weights():
    stale_state = np.zeros((2, 4), dtype=np.float64)
    weights = np.zeros(2, dtype=np.float64)

    assert _phase_tail_state_for_weights(None, stale_state) is None
    assert _phase_tail_state_for_weights(weights, stale_state) is stale_state


def test_phase_suffix_integrates_exponential_linear_frequency_exactly():
    n = np.array([0.0, 0.1, 0.2])
    phi = n.copy()
    expected = np.array([1.0 - math.exp(-0.2), 1.0 - math.exp(-0.1), 0.0])

    np.testing.assert_allclose(_phase_integral_suffix(n, phi), expected, rtol=2e-14, atol=2e-15)
