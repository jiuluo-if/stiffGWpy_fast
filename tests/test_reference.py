# -*- coding: utf-8 -*-
"""Tests for the independent high-accuracy reference pipeline."""

import math
import warnings

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_background_analytic_limits_md_and_radiation():
    """sigma = 1 in the MD reheating branch, ~4/3 just after reheating (no stiff)."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=0.0, DN_eff=0.0)
    n_re_abs = m.derived_param['N_inf'] - m.derived_param['N_re']
    # One e-fold before reheating ends (MD).
    # 再加热结束前一个 e-fold，此时为物质主导阶段。
    _, s_md = REF.background_at(m, n_re_abs - 1.0, 0.0)
    # Just after reheating (radiation dominated, w=1/3 -> sigma = 1 + w = 4/3).
    # 再加热结束后为辐射主导，w=1/3，因此 sigma=1+w=4/3。
    _, s_rad = REF.background_at(m, n_re_abs + 1.0, 0.0)
    assert s_md == pytest.approx(1.0, rel=1e-10)
    assert s_rad == pytest.approx(4.0 / 3.0, rel=2e-3)


def test_background_analytic_limit_stiff():
    """Deep in the stiff era (w=1) sigma tends to 2."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1.0, DN_eff=0.0)
    n_re_abs = m.derived_param['N_inf'] - m.derived_param['N_re']
    _, s_stiff = REF.background_at(m, n_re_abs + 0.05, 0.0)
    assert s_stiff == pytest.approx(2.0, rel=2e-3)


def test_background_matches_fast_grid_away_from_kink():
    """Continuous background reproduces the fast grid sigma away from the kink."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    DN = m.cosmo_param['DN_eff']
    FS.gen_fast(m, 0.01)
    n_re_abs = m.derived_param['N_inf'] - m.derived_param['N_re']
    for N_abs in [n_re_abs + 2.0, n_re_abs + 8.0, m.Nv[-1] - 10.0]:
        i = int(np.argmin(np.abs(m.Nv - N_abs)))
        _, s_ref = REF.background_at(m, float(m.Nv[i]), DN)
        assert abs(s_ref - m.sigma[i]) < 5e-3


def test_f_hor_abs_matches_public_horizon_frequency():
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    dn_eff = m.cosmo_param['DN_eff']
    background = REF._Background(m, dn_eff)

    for N in (0.0, background.N_re_abs, background.N_inf):
        assert REF.f_hor_abs(background, N) == pytest.approx(
            REF.f_hor_at(m, N, dn_eff), abs=1e-14
        )


@pytest.mark.parametrize('freq', [-18.4, -15.0, -12.0, -7.0, 0.0])
def test_find_start_n_reuses_horizon_anchors_and_stops_at_float_resolution(
        model, monkeypatch, freq):
    """Start search must preserve the old root while avoiding redundant H2 work."""
    dn_eff = model.cosmo_param['DN_eff']
    background = REF._Background(model, dn_eff)

    def legacy_search():
        lo, hi = 0.0, background.N_inf

        def z_at(N):
            return (freq - REF.f_hor_abs(background, N)) * math.log(10.0)

        if z_at(lo) >= -12.0:
            return 0.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if z_at(mid) < -12.0:
                lo = mid
            else:
                hi = mid
        return hi

    expected = legacy_search()
    original = REF._H2_and_sigma
    calls = 0

    def counted_h2(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(REF, '_H2_and_sigma', counted_h2)
    actual = REF._find_start_N(background, freq, -12.0)

    assert actual == expected
    assert calls < 100


def test_numba_h2_sigma_matches_python_at_every_background_branch(model):
    """Compiled continuous background preserves exact values and branch edges."""
    dn_eff = model.cosmo_param['DN_eff']
    background = REF._Background(model, dn_eff)
    R_func = REF._H2_and_sigma_python
    REF._warmup_numba_reference_background(background)

    points = [
        0.0,
        np.nextafter(background.N_re_abs, -np.inf),
        background.N_re_abs,
        np.nextafter(background.N_re_abs, np.inf),
        background.N_inf - background.N_fin,
        np.nextafter(background.N_inf - background.N_fin, -np.inf),
        np.nextafter(background.N_inf - background.N_fin, np.inf),
        background.N_inf - background.N_max,
        np.nextafter(background.N_inf - background.N_max, -np.inf),
        np.nextafter(background.N_inf - background.N_max, np.inf),
        background.N_inf,
    ]
    for nu_edge in (100.0, 0.1):
        n_edge = background.N_inf + math.log(nu_edge / background.nu_today)
        points.extend((
            np.nextafter(n_edge, -np.inf),
            n_edge,
            np.nextafter(n_edge, np.inf),
        ))
    points.extend(np.random.default_rng(20261105).uniform(
        0.0, background.N_inf, 2000
    ))

    for N in points:
        expected = R_func(background, float(N))
        actual = REF._H2_and_sigma(background, float(N))
        assert actual == expected


@pytest.mark.parametrize('params, dn_eff', [
    ({'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1e-2},
     0.133732083004284),
    ({'r': 1e-2, 'cr': 1, 'T_re': 1e4, 'kappa10': 1e-2}, 0.2),
    ({'r': 1e-2, 'cr': 1, 'T_re': 1e1, 'kappa10': 1e-2}, 0.2),
    ({'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1.0}, 0.2),
    ({'r': 0.19649188544593907, 'cr': 0,
      'T_re': 1761.314039059193, 'kappa10': 0.011598483151630803},
     2.291411867022825),
])
def test_numba_h2_sigma_matches_python_across_cosmologies(params, dn_eff):
    """The compiled evaluator agrees across reheating, stiffness and cr branches."""
    model = LCDM_SG(**params, DN_eff=dn_eff)
    background = REF._Background(model, dn_eff)
    REF._warmup_numba_reference_background(background)

    points = [
        0.0,
        np.nextafter(background.N_re_abs, -np.inf),
        background.N_re_abs,
        np.nextafter(background.N_re_abs, np.inf),
        background.N_inf - background.N_fin,
        background.N_inf - background.N_max,
        background.N_inf,
    ]
    for nu_edge in (100.0, 0.1):
        edge = background.N_inf + math.log(nu_edge / background.nu_today)
        points.extend((np.nextafter(edge, -np.inf), edge,
                       np.nextafter(edge, np.inf)))
    rng = np.random.default_rng(20261105)
    points.extend(rng.uniform(0.0, background.N_inf, 256))

    for N in points:
        assert REF._H2_and_sigma(background, float(N)) == \
            REF._H2_and_sigma_python(background, float(N))


@pytest.mark.parametrize('params, dn_eff', [
    ({'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1e-2},
     0.133732083004284),
    ({'r': 1e-2, 'cr': 1, 'T_re': 1e4, 'kappa10': 1e-2}, 0.2),
    ({'r': 1e-2, 'cr': 1, 'T_re': 1e1, 'kappa10': 1e-2}, 0.2),
    ({'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1.0}, 0.2),
    ({'r': 0.19649188544593907, 'cr': 0,
      'T_re': 1761.314039059193, 'kappa10': 0.011598483151630803},
     2.291411867022825),
])
def test_numba_tensor_rhs_matches_python_on_piecewise_boundaries(params, dn_eff):
    """A fused Numba tensor RHS preserves the Python derivatives exactly."""
    model = LCDM_SG(**params, DN_eff=dn_eff)
    background = REF._Background(model, dn_eff)
    REF._warmup_numba_reference_background(background)
    rhs = REF._make_numba_tensor_rhs(background)

    points = [background.N_re_abs, background.N_inf - background.N_fin,
              background.N_inf - background.N_max, background.N_inf]
    for nu_edge in (100.0, 0.1):
        edge = background.N_inf + math.log(nu_edge / background.nu_today)
        points.extend((np.nextafter(edge, -np.inf), edge,
                       np.nextafter(edge, np.inf)))
    rng = np.random.default_rng(20261105)
    points.extend(rng.uniform(background.N_re_abs, background.N_inf, 128))
    states = [np.array([z, x, y]) for z, x, y in (
        (-12.0, 0.0, 1.0), (-1.0, 1e-8, -2e-8),
        (0.0, 1.0, 0.5), (8.0, -0.2, 0.3),
    )]

    for N in points:
        for state in states:
            expected = REF._tensor_orig(float(N), state, background)
            actual = rhs(float(N), state, background)
            assert actual == expected


def test_numba_tensor_rhs_preserves_python_overflow_failures(model):
    """Extreme exponents retain the existing Python OverflowError behavior."""
    dn_eff = model.cosmo_param['DN_eff']
    background = REF._Background(model, dn_eff)
    rhs = REF._make_numba_tensor_rhs(background)
    background.N_re_abs = -1000.0
    with pytest.raises(OverflowError):
        rhs(background.N_inf - 711.0, np.array([0.0, 0.0, 1.0]), background)
    with pytest.raises(OverflowError):
        rhs(background.N_inf, np.array([711.0, 0.0, 1.0]), background)


@pytest.mark.parametrize('freqs, expected_warmups', [
    (np.array([-12.0]), 0),
    (np.linspace(-18.0, -7.0, 32), 1),
])
def test_spectrum_reference_jit_warmup_is_batched(model, monkeypatch, freqs,
                                                    expected_warmups):
    """Do not pay the H2 JIT compile cost for a one-frequency reference call."""
    h2_calls = []
    rhs_calls = []
    monkeypatch.setattr(REF, '_NUMBA_H2_SIGMA_READY', False)
    monkeypatch.setattr(REF, '_NUMBA_TENSOR_RHS_READY', False)
    monkeypatch.setattr(REF, '_warmup_numba_spline', lambda: None)

    def warm_h2(bg):
        h2_calls.append(bg)
        monkeypatch.setattr(REF, '_NUMBA_H2_SIGMA_READY', True)

    monkeypatch.setattr(
        REF, '_warmup_numba_reference_background', warm_h2,
    )
    monkeypatch.setattr(
        REF, '_warmup_numba_reference_tensor_rhs',
        lambda bg: rhs_calls.append(bg),
    )
    monkeypatch.setattr(
        REF, 'solve_reference_mode',
        lambda *args, **kwargs: {
            'Ogw_today': 1.0, 'Oj_today': 0.0, 'Opgw_today': 1.0,
            'used_tail': False,
        },
    )

    REF.spectrum_reference(
        model, freqs, model.cosmo_param['DN_eff'], workers=1
    )

    assert len(h2_calls) == expected_warmups
    assert len(rhs_calls) == expected_warmups


def test_spectrum_reference_prufer_batch_skips_cartesian_rhs_warmup(model,
                                                                    monkeypatch):
    calls = []
    monkeypatch.setattr(REF, '_NUMBA_H2_SIGMA_READY', True)
    monkeypatch.setattr(REF, '_NUMBA_TENSOR_RHS_READY', False)
    monkeypatch.setattr(REF, '_warmup_numba_spline', lambda: None)
    monkeypatch.setattr(
        REF, '_warmup_numba_reference_tensor_rhs',
        lambda bg: calls.append(bg),
    )
    monkeypatch.setattr(
        REF, 'solve_reference_mode',
        lambda *args, **kwargs: {
            'Ogw_today': 1.0, 'Oj_today': 0.0, 'Opgw_today': 1.0,
            'used_tail': False,
        },
    )

    REF.spectrum_reference(
        model, np.linspace(-18.0, -7.0, 32), model.cosmo_param['DN_eff'],
        workers=1, backend='prufer',
    )

    assert calls == []


def test_numba_h2_sigma_preserves_cr_positive_mode_solutions(monkeypatch):
    """The compiled background also preserves a positive-cr Reference solve."""
    dn_eff = 0.133732083004284
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2, DN_eff=dn_eff)
    background = REF._Background(m, dn_eff)
    REF._warmup_numba_reference_background(background)
    compiled = REF._H2_and_sigma
    freqs = (-14.6013, -12.6723)

    for freq in freqs:
        actual = REF.solve_reference_mode(
            m, freq, dn_eff, z_tail=10.0, rtol=1e-11
        )
        monkeypatch.setattr(REF, '_H2_and_sigma', REF._H2_and_sigma_python)
        expected = REF.solve_reference_mode(
            m, freq, dn_eff, z_tail=10.0, rtol=1e-11
        )
        monkeypatch.setattr(REF, '_H2_and_sigma', compiled)

        for key in ('Ogw_today', 'Oj_today', 'Opgw_today', 'used_tail',
                    'event_N', 'n_steps', 'z0'):
            assert actual[key] == expected[key]


def test_numba_h2_sigma_preserves_python_overflow_and_nonfinite_behavior(model):
    """The compiled path keeps Python's explicit overflow failure behavior."""
    dn_eff = model.cosmo_param['DN_eff']
    background = REF._Background(model, dn_eff)
    REF._warmup_numba_reference_background(background)

    with np.errstate(all='ignore'):
        expected = REF._H2_and_sigma_python(
            background, background.N_inf + 1000.0
        )
        actual = REF._H2_and_sigma(
            background, background.N_inf + 1000.0
        )
    assert np.isinf(expected[0]) and expected[1] == 0.0
    assert np.isinf(actual[0]) and actual[1] == 0.0

    background.N_re_abs = -1000.0
    overflow_N = background.N_inf - 1000.0
    with pytest.raises(OverflowError):
        REF._H2_and_sigma_python(background, overflow_N)
    with pytest.raises(OverflowError):
        REF._H2_and_sigma(background, overflow_N)

    singular = REF._Background(model, dn_eff)
    singular.N_inf = 0.0
    singular.N_re_abs = -1.0
    singular.nu_today = 1000.0
    singular.Omh2 = 1.0
    singular.Omega_mnuh2 = 0.0
    singular.Osh2 = 0.0
    singular.Oerh2 = 0.0
    singular.Omega_ph2 = 0.0
    singular.Omega_nh2 = 0.0
    singular.OLh2 = np.float64(-1.0)
    with np.errstate(all='ignore'):
        expected_zero = REF._H2_and_sigma_python(singular, 0.0)
        actual_zero = REF._H2_and_sigma(singular, 0.0)
    assert expected_zero == actual_zero == (0.0, math.inf)
    with np.errstate(all='ignore'):
        expected_rhs = REF._tensor_orig(0.0, np.array([0.0, 1.0, 1.0]), singular)
        actual_rhs = REF._make_numba_tensor_rhs(singular)(
            0.0, np.array([0.0, 1.0, 1.0]), singular)
    assert actual_rhs == expected_rhs


def test_integrate_spectrum_exponential():
    """integrate_spectrum reproduces a simple analytic integral."""
    freqs = np.linspace(-18.0, 6.0, 401)
    Ogw = 10.0 ** (-14.0 - 0.3 * freqs)   # 平滑幂律谱。
    Oj = np.zeros_like(Ogw)
    g2, qerr, ierr = REF.integrate_spectrum(freqs, Ogw, Oj)
    # Analytical: integral of 10^{-14-0.3 x} over x in [a,b], times ln10.
    # 解析值为区间 [a,b] 上 10^(-14-0.3x) 的积分，再乘 ln(10)。
    a, b = freqs[0], freqs[-1]
    ln10 = math.log(10.0)
    exact = ln10 * (-1.0 / (0.3 * ln10)) * (
        (10.0 ** (-14.0 - 0.3 * b)) - (10.0 ** (-14.0 - 0.3 * a)))
    assert g2 == pytest.approx(exact, rel=1e-3)
    assert qerr >= 0.0


def test_integrate_spectrum_rescales_subnormal_dynamic_range():
    """Tiny physical spectra should not overflow PCHIP's reciprocal slopes."""
    freqs = np.linspace(0.0, 10.0, 257)
    unit_scale_spectrum = np.exp(-48.0 * freqs)
    scale = 1e-106
    tiny_spectrum = scale * unit_scale_spectrum

    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        actual, qerr, ierr = REF.integrate_spectrum(
            freqs, tiny_spectrum, np.zeros_like(tiny_spectrum)
        )
    expected, expected_qerr, expected_ierr = REF.integrate_spectrum(
        freqs, unit_scale_spectrum, np.zeros_like(unit_scale_spectrum)
    )

    assert actual == pytest.approx(scale * expected, rel=2e-13, abs=0.0)
    assert qerr == pytest.approx(scale * expected_qerr, rel=2e-13, abs=0.0)
    assert ierr == pytest.approx(scale * expected_ierr, rel=2e-13, abs=0.0)


def test_integrate_spectrum_scales_quad_budget_with_breakpoint_count():
    """A dense reference grid must not exceed quad's fixed interval budget."""
    freqs = np.linspace(-18.0, -7.0, 1600)
    Ogw = 10.0 ** (-14.0 - 0.3 * freqs)

    g2, qerr, _ = REF.integrate_spectrum(freqs, Ogw, np.zeros_like(Ogw))

    assert np.isfinite(g2)
    assert qerr is not None and qerr >= 0.0


def test_integrate_spectrum_interpolation_error_matches_spline_integrals():
    freqs = np.linspace(-20.0, 7.0, 247)
    spectrum = np.exp(-0.23 * (freqs + 4.0)) * (
        1.0
        + 0.17 * np.exp(-0.5 * ((freqs + 5.7) / 0.21) ** 2)
        + 0.08 * np.exp(-0.5 * ((freqs - 0.6) / 0.48) ** 2)
        + 0.012 * np.sin(3.7 * freqs)
    )
    _, _, interpolation_error = REF.integrate_spectrum(
        freqs, spectrum, np.zeros_like(spectrum)
    )

    half = freqs[::2]
    expected = math.log(10.0) * abs(
        PchipInterpolator(freqs, spectrum).integrate(freqs[0], freqs[-1])
        - PchipInterpolator(half, spectrum[::2]).integrate(freqs[0], freqs[-1])
    )
    assert interpolation_error == pytest.approx(expected, rel=1e-10, abs=1e-11)


def test_integrate_spectrum_even_grid_includes_last_coarse_endpoint():
    freqs = np.linspace(-22.0, 8.0, 240)
    spectrum = np.exp(-0.18 * (freqs + 2.0)) * (
        1.0
        + 0.2 * np.exp(-0.5 * ((freqs + 8.0) / 0.3) ** 2)
        + 0.05 * np.sin(2.1 * freqs)
    )
    _, _, interpolation_error = REF.integrate_spectrum(
        freqs, spectrum, np.zeros_like(spectrum)
    )

    coarse_indices = np.arange(0, freqs.size, 2)
    coarse_indices = np.append(coarse_indices, freqs.size - 1)
    coarse_freqs = freqs[coarse_indices]
    coarse_spectrum = spectrum[coarse_indices]
    expected = math.log(10.0) * abs(
        PchipInterpolator(freqs, spectrum).integrate(freqs[0], freqs[-1])
        - PchipInterpolator(coarse_freqs, coarse_spectrum).integrate(
            freqs[0], freqs[-1]
        )
    )
    assert interpolation_error == pytest.approx(expected, rel=1e-10, abs=1e-11)


def test_summarize_tail_convergence_reports_observed_bound():
    rows = [
        {'z_tail': 5.0, 'DN_gw': 1.0300},
        {'z_tail': 6.0, 'DN_gw': 1.0120},
        {'z_tail': 7.0, 'DN_gw': 1.0040},
        {'z_tail': 8.0, 'DN_gw': 1.0015},
        {'z_tail': 10.0, 'DN_gw': 1.0005},
    ]
    summary = REF.summarize_tail_convergence(rows)
    assert summary['central_z_tail'] == 10.0
    assert summary['central_DN_gw'] == pytest.approx(1.0005)
    assert summary['systematic_uncertainty_abs'] == pytest.approx(0.0295)
    assert summary['systematic_uncertainty_rel'] == pytest.approx(0.0295 / 1.0005)
    assert summary['empirical_exponential_decay_rate'] is not None
    assert summary['convergence_monotone_to_central']


def test_reference_mode_reports_tail_phase_and_adiabaticity():
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    sol = REF.solve_reference_mode(m, -8.0, m.cosmo_param['DN_eff'],
                                   z_tail=5.0, rtol=1e-8)
    for name in ('phase_handoff', 'amplitude_handoff', 'omega_handoff',
                 'omega_prime_over_omega2', 'omega_second_over_omega3',
                 'event_N'):
        assert name in sol
        assert math.isfinite(sol[name])
    assert sol['omega_handoff'] == pytest.approx(math.exp(5.0), rel=1e-8)
    assert sol['amplitude_handoff'] > 0.0


@pytest.mark.parametrize(
    ('freq', 'expected_tail'),
    [(-8.0, True), (-20.0, False)],
)
def test_prufer_reference_backend_matches_cartesian_mode(freq, expected_tail):
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    dn_eff = m.cosmo_param['DN_eff']
    cartesian = REF.solve_reference_mode(
        m, freq, dn_eff, z_tail=5.0, rtol=1e-9, backend='cartesian')
    prufer = REF.solve_reference_mode(
        m, freq, dn_eff, z_tail=5.0, rtol=1e-9, backend='prufer')

    assert cartesian['used_tail'] is expected_tail
    assert prufer['used_tail'] is expected_tail
    for name in ('Ogw_today', 'Oj_today', 'Opgw_today'):
        assert prufer[name] == pytest.approx(cartesian[name], rel=1e-4, abs=1e-30)


def test_prufer_reference_backend_fails_closed_on_ode_failure(monkeypatch):
    from types import SimpleNamespace

    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    monkeypatch.setattr(
        REF, 'solve_ivp',
        lambda *args, **kwargs: SimpleNamespace(success=False,
                                                  message='mocked solver failure'))

    with pytest.raises(RuntimeError, match='mocked solver failure'):
        REF.solve_reference_mode(
            m, -8.0, m.cosmo_param['DN_eff'], z_tail=5.0, rtol=1e-8,
            backend='prufer')


def test_cartesian_tail_uses_event_state_without_dense_output(monkeypatch):
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    original_solve_ivp = REF.solve_ivp
    dense_output_values = []

    def record_solve_ivp(*args, **kwargs):
        dense_output_values.append(kwargs.get('dense_output'))
        return original_solve_ivp(*args, **kwargs)

    monkeypatch.setattr(REF, 'solve_ivp', record_solve_ivp)
    result = REF.solve_reference_mode(
        m, -8.0, m.cosmo_param['DN_eff'], z_tail=5.0, rtol=1e-9)

    assert result['used_tail']
    assert dense_output_values == [False]
    assert math.isfinite(result['Ogw_today'])


@pytest.mark.parametrize(
    'name', ('spl_rho', 'spl_rhop', '_FD_RHO_CS', '_FD_P_CS'))
def test_reference_compiled_cubic_spline_matches_scipy(name):
    spline = (getattr(REF.gp, name) if name.startswith('spl_')
              else getattr(REF, name))
    points = np.concatenate((
        np.array([spline.x[0] - 0.1, spline.x[0], spline.x[-1],
                  spline.x[-1] + 0.1, np.nan]),
        np.linspace(spline.x[0], spline.x[-1], 51),
    ))
    expected = np.asarray(spline(points), dtype=float)
    for evaluator in (REF._cubic_spline_eval,
                      REF._cubic_spline_eval_python):
        actual = np.array([
            evaluator(spline.x, spline.c, float(point))
            for point in points
        ])
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=3e-16,
                                   equal_nan=True)


def test_reference_spline_dispatch_avoids_cold_jit_for_cheap_modes():
    assert REF._spline_evaluator_for_frequency(
        -18.5, compiled_ready=False) is REF._cubic_spline_eval_python
    assert REF._spline_evaluator_for_frequency(
        -14.0, compiled_ready=False) is REF._cubic_spline_eval
    assert REF._spline_evaluator_for_frequency(
        -18.5, compiled_ready=True) is REF._cubic_spline_eval


def test_spectrum_reference_warms_compiled_spline_before_mode_dispatch(monkeypatch):
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    events = []

    def fake_mode(*args, **kwargs):
        events.append('mode')
        return {'Ogw_today': 1.0, 'Oj_today': 0.0, 'Opgw_today': 1.0,
                'used_tail': False}

    monkeypatch.setattr(REF, '_sync_size', lambda: 1)
    monkeypatch.setattr(REF, '_warmup_numba_spline',
                        lambda: events.append('warmup'))
    monkeypatch.setattr(REF, 'solve_reference_mode', fake_mode)

    REF.spectrum_reference(
        m, np.array([-18.0, -14.0, np.nan]),
        m.cosmo_param['DN_eff'], workers=2)

    assert events == ['warmup', 'mode', 'mode', 'mode']


def test_exact_zero_spectrum_does_not_warm_numba(monkeypatch):
    m = LCDM_SG(r=0.0, cr=0.0, T_re=2e3, kappa10=1e-2)

    def unexpected_warmup():
        raise AssertionError('exact zero tensor should skip spline JIT warmup')

    monkeypatch.setattr(REF, '_warmup_numba_spline', unexpected_warmup)
    monkeypatch.setattr(REF, '_sync_size', lambda: 1)
    monkeypatch.setattr(
        REF, 'solve_reference_mode',
        lambda *args, **kwargs: {
            'Ogw_today': 0.0, 'Oj_today': 0.0, 'Opgw_today': 0.0,
            'used_tail': False,
        })

    result = REF.spectrum_reference(
        m, np.array([-14.0]), m.cosmo_param['DN_eff'], workers=1)

    assert np.array_equal(result[0], np.array([0.0]))


def test_prufer_rhs_rejects_frequency_exponential_underflow():
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    bg = REF._Background(m, m.cosmo_param['DN_eff'])

    with pytest.raises(FloatingPointError, match='frequency underflow'):
        REF._tensor_prufer(bg.N_inf, (-1000.0, 0.0, 0.0), bg)


def test_run_reference_forwards_optional_mode_backend():
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    freqs = np.array([-8.0, -15.0, -20.0, -21.0])
    dn_eff = m.cosmo_param['DN_eff']

    result = REF.run_reference(
        m, dn_eff=dn_eff, freq_subset=freqs, z_tail=5.0, rtol=1e-9,
        self_consistent=False, backend='prufer')

    assert result['mode_backend'] == 'prufer'
    assert result['DN_eff'] == pytest.approx(dn_eff)
    assert result['n_freq'] == len(freqs)


def test_run_reference_forwards_z_start_to_spectrum_modes(monkeypatch):
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    calls = []

    def fake_mode(*args, **kwargs):
        calls.append(kwargs)
        return {'Ogw_today': 1.0, 'Oj_today': 0.0, 'Opgw_today': 1.0,
                'used_tail': False}

    monkeypatch.setattr(REF, 'solve_reference_mode', fake_mode)
    REF.run_reference(
        m, dn_eff=m.cosmo_param['DN_eff'], freq_subset=np.array([-8.0]),
        self_consistent=False, z_start=-9.0, backend='prufer')

    assert calls[0]['z_start'] == -9.0


def test_tail_diagnostic_cache_key_binds_tail_and_frequency():
    from scripts.benchmark_tail_diagnostics import _cache_key

    base = _cache_key('default', {'r': 1e-2}, [-2.0, 0.0], 0.002, 5.0,
                      1e-10)
    assert base != _cache_key('default', {'r': 1e-2}, [-2.0, 0.0], 0.002,
                              6.0, 1e-10)
    assert base != _cache_key('default', {'r': 1e-2}, [-2.0, 1.0], 0.002,
                              5.0, 1e-10)


def test_phase_window_observable_uses_later_handoff():
    from scripts.benchmark_phase_averaged_oracle import _phase_window_observable

    result = _phase_window_observable(
        z_handoff=5.0, event_N=1.0, x_handoff=2.0, y_handoff=3.0,
        n_inf=4.0, today_f_hor=-2.0, freq=-1.0, tensor_power=1.0)

    assert result['observable'] == pytest.approx(100.0 * 6.5 * math.exp(-16.0) / 36.0)
    assert result['phase_averaged'] is True


def test_prufer_derivatives_preserve_subhorizon_amplitude_equations():
    from scripts.benchmark_prufer_oracle import _prufer_derivatives

    z_prime, log_amplitude_prime, phase_prime = _prufer_derivatives(
        z=5.0, sigma=2.0, phase=0.0)

    assert z_prime == pytest.approx(2.0)
    assert log_amplitude_prime == pytest.approx(2.0)
    assert phase_prime == pytest.approx(-math.exp(5.0))


def test_prufer_phase_averaged_power_uses_half_amplitude_square():
    from scripts.benchmark_prufer_oracle import _phase_averaged_power

    assert _phase_averaged_power(4.0) == pytest.approx(8.0)


def test_prufer_phase_averaged_today_observable_matches_tail_formula():
    from scripts.benchmark_prufer_oracle import _phase_averaged_today

    result = _phase_averaged_today(
        coefficient_squared=6.5, z_handoff=5.0, event_N=1.0,
        n_inf=4.0, today_f_hor=-2.0, freq=-1.0, tensor_power=1.0)
    expected_x_squared = 6.5 * math.exp(-16.0) * 100.0

    assert result['Opgw_today'] == pytest.approx(expected_x_squared / 36.0)
    assert result['Oj_today'] == pytest.approx(-6.5 * math.exp(-16.0) / 3.0)
    assert result['Ogw_today'] == pytest.approx(
        3.0 * result['Opgw_today'] + result['Oj_today'])


def test_prufer_outer_convergence_metric_matches_reference_definition():
    from scripts.benchmark_prufer_oracle import _outer_convergence_metric
    from stiffgwpy_fast import global_param as gp

    assert _outer_convergence_metric(2.0, 0.0024, 0.0020) == pytest.approx(
        (gp.Neff0 + 2.0024) / (gp.Neff0 + 2.0020) - 1.0)


def test_oracle_checkpoint_key_and_schema_are_reproducible(tmp_path):
    from scripts.benchmark_oracle_tail_convergence import (
        _cache_key,
        _load_checkpoint,
        _new_checkpoint,
        _save_checkpoint,
    )

    parameters = {'r': 1e-2, 'cr': 1, 'T_re': 2e3, 'kappa10': 1e-2}
    frequencies = [-2.0, 0.0, 1.0]
    first = _cache_key('default', parameters, frequencies, 0.002, 5.0, 1e-10)
    second = _cache_key('default', parameters, frequencies, 0.002, 6.0, 1e-10)
    assert first != second
    path = str(tmp_path / 'oracle.checkpoint.json')
    checkpoint = _new_checkpoint()
    checkpoint['records']['default'] = {'5.0': {'cache_key': first}}
    _save_checkpoint(path, checkpoint)
    restored = _load_checkpoint(path)
    assert restored['records']['default']['5.0']['cache_key'] == first


@pytest.mark.slow
def test_reference_mode_vs_fast_mid_frequency():
    """Reference ODE agrees with the fast solver in the resolved mid-frequency region."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.apply_accuracy_mode('production')
    FS.SGWB_iter_fast(m, tol=1e-7)
    dn = m.cosmo_param['DN_eff']
    logf = -8.0
    # m.f / m.log10OmegaGW are stored descending; sort with a shared
    # argsort so the (f, log10Omega) pairing survives interpolation.
    # m.f 和 m.log10OmegaGW 均按降序保存；用同一个 argsort 排序，保持插值时频率与谱值配对。
    o = np.argsort(np.asarray(m.f, dtype=float))
    lo_fast = np.interp(logf, np.asarray(m.f)[o],
                        np.asarray(m.log10OmegaGW)[o])
    sol = REF.solve_reference_mode(m, logf, dn, z_tail=5.0, rtol=1e-11)
    lo_ref = math.log10(sol['Ogw_today'] - sol['Oj_today'])
    assert abs(lo_ref - lo_fast) < 5e-2
    assert sol['used_tail']
    assert math.isfinite(sol['Ogw_today'])


def test_apply_reference_to_model(monkeypatch):
    """apply_reference_to_model exposes the reference result on the model."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    freqs = np.array([6.0, 4.0, 2.0, 0.0])
    logO = np.array([-15.6, -15.6, -11.6, -7.9])
    fake = dict(freqs=freqs, log10OmegaGW=logO, DN_eff=2.27e-3, DN_gw=2.27e-3,
                kappa_r=1.94e-3, g2=2.81e-8)
    monkeypatch.setattr(REF, 'run_reference', lambda *a, **k: fake)
    REF.apply_reference_to_model(m)
    assert m.cosmo_param['DN_eff'] == pytest.approx(2.27e-3)
    assert m.DN_gw[-1] == pytest.approx(2.27e-3)
    assert m.kappa_r == pytest.approx(1.94e-3)
    assert np.allclose(m.log10OmegaGW, logO)
    assert m.SGWB_converge


def test_reference_physical_guard_is_optional_for_large_fixed_point(monkeypatch):
    """The physical total-Neff guard must not cap numeric exploration by default."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2, DN_eff=0.0)
    freqs = np.array([-1.0, 0.0, 1.0, 2.0])
    omega_nu = REF.gp.Omega_nh2 / m.derived_param['h']**2
    target_dn = 12.0
    target_g2 = target_dn * omega_nu / REF.gp.Neff0

    def fake_spectrum(model, frequencies, dn_eff, **kwargs):
        count = len(frequencies)
        return (np.ones(count), np.zeros(count), np.zeros(count),
                np.zeros(count, dtype=bool))

    monkeypatch.setattr(REF, 'spectrum_reference', fake_spectrum)
    monkeypatch.setattr(REF, 'integrate_spectrum',
                        lambda *args, **kwargs: (target_g2, 0.0, 0.0))

    result = REF.run_reference(
        m, freq_subset=freqs, z_tail=5.0, rtol=1e-9, tol=1e-8)

    assert result['DN_gw'] == pytest.approx(target_dn)
    assert result['DN_eff'] == pytest.approx(target_dn)
    with pytest.raises(RuntimeError, match='physical guard'):
        REF.run_reference(
            m, freq_subset=freqs, z_tail=5.0, rtol=1e-9, tol=1e-8,
            physical_guard=True)


def test_reference_zero_tensor_limit_is_exact_and_grid_safe():
    """r=0, cr=0 has a zero spectrum independent of the input tilt."""
    from stiffgwpy_fast.freq_adaptive import grid_independent_freqs

    m = LCDM_SG(r=0.0, cr=0.0, n_t=0.5, T_re=2e3,
                kappa10=1e-2, DN_re=10.0, DN_eff=0.0)
    grid, _, _ = grid_independent_freqs(m, 1.0)
    assert np.isfinite(grid).all()

    freqs = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
    for backend in ('cartesian', 'prufer'):
        result = REF.run_reference(
            m, dn_eff=0.0, freq_subset=freqs, z_tail=5.0, rtol=1e-9,
            self_consistent=False, backend=backend)

        assert result['DN_gw'] == 0.0
        assert result['mode_backend'] == backend
        assert np.all(result['Ogw'] == 0.0)
        assert np.all(result['Oj'] == 0.0)
        assert np.all(result['Opgw'] == 0.0)
        assert np.all(np.isneginf(result['log10OmegaGW']))

    single_field = LCDM_SG(r=0.0, cr=1.0, T_re=2e3, kappa10=1e-2)
    with pytest.raises(ValueError, match='positive r'):
        REF.run_reference(
            single_field, dn_eff=0.0, freq_subset=freqs,
            self_consistent=False)

    underflow = LCDM_SG(
        r=1e-300, cr=0.0, n_t=0.0, A_s=1e-60, T_re=2e3,
        kappa10=1e-2, DN_re=10.0, DN_eff=0.0,
    )
    with pytest.raises(FloatingPointError, match='underflow'):
        REF.run_reference(
            underflow, dn_eff=0.0, freq_subset=freqs,
            self_consistent=False)
