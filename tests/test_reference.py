# -*- coding: utf-8 -*-
"""Tests for the independent high-accuracy reference pipeline."""

import math

import numpy as np
import pytest

from stiffgwpy import fast_sgwb as FS
from stiffgwpy import reference as REF
from stiffgwpy.stiff_SGWB import LCDM_SG


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


@pytest.mark.slow
def test_reference_mode_vs_fast_mid_frequency():
    """Reference ODE agrees with the fast solver in the resolved mid-frequency region."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.apply_accuracy_mode('production')
    FS.SGWB_iter_fast(m, tol=1e-7)
    dn = m.cosmo_param['DN_eff']
    logf = -8.0
    lo_fast = np.interp(logf, np.sort(m.f), np.sort(m.log10OmegaGW))
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
