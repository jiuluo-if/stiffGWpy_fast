# -*- coding: utf-8 -*-
"""Tests for the independent high-accuracy reference pipeline."""

import math

import numpy as np
import pytest

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_background_analytic_limits_md_and_radiation():
    """sigma = 1 in the MD reheating branch, ~4/3 just after reheating (no stiff)."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=0.0, DN_eff=0.0)
    n_re_abs = m.derived_param['N_inf'] - m.derived_param['N_re']
    # One e-fold before reheating ends (MD).
    _, s_md = REF.background_at(m, n_re_abs - 1.0, 0.0)
    # Just after reheating (radiation dominated, w=1/3 -> sigma = 1 + w = 4/3).
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


def test_integrate_spectrum_exponential():
    """integrate_spectrum reproduces a simple analytic integral."""
    freqs = np.linspace(-18.0, 6.0, 401)
    Ogw = 10.0 ** (-14.0 - 0.3 * freqs)   # smooth power law
    Oj = np.zeros_like(Ogw)
    g2, qerr, ierr = REF.integrate_spectrum(freqs, Ogw, Oj)
    # Analytical: integral of 10^{-14-0.3 x} over x in [a,b], times ln10.
    a, b = freqs[0], freqs[-1]
    ln10 = math.log(10.0)
    exact = ln10 * (-1.0 / (0.3 * ln10)) * (
        (10.0 ** (-14.0 - 0.3 * b)) - (10.0 ** (-14.0 - 0.3 * a)))
    assert g2 == pytest.approx(exact, rel=1e-3)
    assert qerr >= 0.0


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
