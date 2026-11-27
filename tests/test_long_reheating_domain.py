import numpy as np
import pytest

from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.fast_sgwb import (
    SGWB_iter_fast,
    _defer_nonadiabatic_wkb_tail,
    _defer_nonadiabatic_wkb_tails_kernel,
    resolve_config,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_wkb_tail_defers_past_a_future_nonadiabatic_background_segment():
    z_grid = np.array([4.8, 5.1, 5.2, 5.3, 6.0])
    sigma = np.array([1.0, 1.0, 100.0, 1.0, 1.0])

    deferred = _defer_nonadiabatic_wkb_tail(
        z_grid, sigma, z_tail=5.0, epsilon_limit=0.02
    )

    assert 5.3 < deferred <= 6.0


def test_wkb_tail_checks_the_state_before_the_first_crossing_node():
    z_grid = np.array([4.8, 4.9, 5.1, 5.3, 6.0])
    sigma = np.array([1.0, 100.0, 1.0, 1.0, 1.0])

    deferred = _defer_nonadiabatic_wkb_tail(
        z_grid, sigma, z_tail=5.0, epsilon_limit=0.02
    )

    assert 5.1 < deferred <= 5.3


def test_wkb_tail_keeps_the_existing_handoff_for_safe_and_phase_active_modes():
    z_grid = np.array([4.8, 5.1, 6.0])
    sigma = np.ones_like(z_grid)

    assert _defer_nonadiabatic_wkb_tail(
        z_grid, sigma, z_tail=5.0, epsilon_limit=0.02
    ) == 5.0


def test_fused_wkb_tail_screen_matches_the_scalar_definition():
    phi = np.array([0.0, 0.2, 0.3, 0.4, 1.0])
    sigma = np.array([1.0, 100.0, 1.0, 1.0, 1.0])
    z0s = np.array([4.8, 4.8])
    j0s = np.array([0, 0], dtype=np.int64)
    tails = np.array([5.0, 5.0])
    phases = np.array([0.0, 1.0])

    deferred_count, unresolved_count = _defer_nonadiabatic_wkb_tails_kernel(
        tails, z0s, j0s, phi, sigma, phases, True, 0.02
    )

    expected = _defer_nonadiabatic_wkb_tail(
        z0s[0] + phi - phi[0], sigma, z_tail=5.0, epsilon_limit=0.02
    )
    assert deferred_count == 1
    assert unresolved_count == 0
    assert tails[0] == expected
    assert tails[1] == 5.0


def test_wkb_tail_screen_accepts_an_equivalent_precomputed_background_summary():
    phi = np.array([0.0, 0.2, 0.3, 0.4, 1.0])
    sigma = np.array([1.0, 100.0, 1.0, 1.0, 1.0])
    z0s = np.array([4.8, 4.8])
    j0s = np.array([0, 0], dtype=np.int64)
    phases = np.array([0.0, 1.0])
    scanned_tails = np.array([5.0, 5.0])
    summarized_tails = scanned_tails.copy()

    scanned = _defer_nonadiabatic_wkb_tails_kernel(
        scanned_tails, z0s, j0s, phi, sigma, phases, True, 0.02
    )
    summarized = _defer_nonadiabatic_wkb_tails_kernel(
        summarized_tails, z0s, j0s, phi, sigma, phases, True, 0.02,
        True, True, 0.6, 148.5,
    )

    assert summarized == scanned
    np.testing.assert_array_equal(summarized_tails, scanned_tails)


def test_fused_wkb_screen_marks_a_nonadiabatic_segment_that_reaches_the_terminal_node():
    phi = np.array([0.0, 0.1, 0.2])
    sigma = np.array([1.0, 1.0, 100.0])
    tails = np.array([5.0])

    deferred_count, unresolved_count = _defer_nonadiabatic_wkb_tails_kernel(
        tails, np.array([4.9]), np.array([0], dtype=np.int64),
        phi, sigma, np.empty(0), False, 0.02,
    )

    assert deferred_count == 1
    assert unresolved_count == 1
    assert tails[0] > 5.0


def test_positive_h2_s2_overflow_matches_the_reference_fixed_point_map():
    model = LCDM_SG(
        Omega_bh2=1.0e5,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=0.01,
        n_t=0.0,
        cr=0.0,
        T_re=2000.0,
        DN_re=10.0,
        kappa10=0.01,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        tol=1e-6,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert model.DN_gw[-1] == pytest.approx(0.012394368479831894, rel=2e-4)
    assert np.isfinite(model.log10OmegaGW).all()
    assert model._fast_nonadiabatic_handoff_selected > 0


def test_steep_positive_h2_background_preserves_a_previously_successful_point():
    model = LCDM_SG(
        Omega_bh2=3.0e5,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=0.01,
        n_t=0.0,
        cr=0.0,
        T_re=4000.0,
        DN_re=10.0,
        kappa10=0.01,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        tol=1e-6,
        physical_guard=True,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert model.DN_gw[-1] == pytest.approx(0.041791685648698645, rel=2e-4)
    assert np.isfinite(model.log10OmegaGW).all()
    assert model._fast_h2_identity_used
    assert _defer_nonadiabatic_wkb_tail(
        np.array([4.8, 5.1, 5.2, 5.3, 6.0]),
        np.array([1.0, 1.0, 100.0, 1.0, 1.0]),
        z_tail=5.0,
        epsilon_limit=0.02,
        phase_active=True,
    ) == 5.0


@pytest.mark.parametrize("duration", [680.0, 800.0])
def test_inactive_matter_prefix_preserves_complete_tensor_solution(duration):
    # The blue cutoff and every initial condition lie after reheating here.
    params = dict(cr=0.0, T_re=2000.0, kappa10=0.01, n_t=2.0, r=0.01)
    models = [LCDM_SG(DN_re=d, **params) for d in (16.0, duration)]
    for model in models:
        assert model.SGWB_iter(engine="fast", accuracy_mode="fast", threads=2,
                               fallback=False, physical_guard=False) is not None
        assert np.all(np.isfinite(model.log10OmegaGW))
        assert model.f.max() + 3.0 < model.f_re
    short, long = models
    np.testing.assert_allclose(long.f, short.f, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(long.log10OmegaGW, short.log10OmegaGW, atol=1e-7, rtol=0.0)
    offset = round((duration - 16.0) / 0.005)
    for name in ("DN_gw", "g2", "w2"):
        values = getattr(long, name)
        assert np.all(values[:offset] == 0.0)
        np.testing.assert_allclose(values[offset:], getattr(short, name),
                                   atol=1e-18, rtol=1e-5)


@pytest.mark.parametrize("T_re", [10.0, 2000.0, 1.0e5])
@pytest.mark.parametrize("kappa10", [0.001, 0.01, 0.1])
def test_long_reheating_cutoff_slice_remains_evaluable(T_re, kappa10):
    seed = LCDM_SG(cr=0.0, T_re=T_re, kappa10=kappa10, r=0.01,
                   n_t=1.0, A_s=2.1e-9, DN_re=0.0)
    N_re = float(seed.derived_param["N_re"])
    H_over_2pi = float(seed.derived_param["H_0"]) / (2.0 * np.pi)
    amplitude = (gp.f_piv / (1.01 * H_over_2pi)) / 0.01
    model = LCDM_SG(
        cr=0.0,
        T_re=T_re,
        kappa10=kappa10,
        r=0.01,
        n_t=1.0,
        A_s=amplitude,
        DN_re=740.0 - N_re,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        tol=1e-6,
    )

    assert result is model
    assert model.derived_param["N_inf"] == pytest.approx(740.0)
    assert model.fast_failure_reason is None
    assert np.all(np.isfinite(model.log10OmegaGW))

