# -*- coding: utf-8 -*-
"""Physics/analytic-limit consistency tests for the SGWB pipeline."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def _z_of_N(m, freq, N, DN):
    return (freq - REF.f_hor_at(m, N, DN)) * math.log(10.0)


def test_deep_subhorizon_wkb_frequency():
    """In the deep sub-horizon the tensor mode oscillates with d(theta)/dN = e^z
    (the WKB / horizon-crossing frequency).  We count mode sign changes and
    compare with the integral of e^z over the window."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    dn = m.cosmo_param['DN_eff']
    freq = -8.0

    def rhs(N, s):
        z, xh, yh = s
        _, sig = REF.background_at(m, N, dn)
        ez = math.exp(z)
        return (1.5 * sig - 1.0, -xh - ez * yh, yh + ez * xh)

    # Start at z = 3 (already deep sub-horizon) and integrate a few e-folds.
    N0 = None
    # z increases with N; bisection for z(N) = 3.
    N_inf = m.derived_param['N_inf']
    lo, hi = 0.0, N_inf
    if _z_of_N(m, freq, lo, dn) >= 3.0:
        N0 = 0.0
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if _z_of_N(m, freq, mid, dn) < 3.0:
                lo = mid
            else:
                hi = mid
        N0 = hi
    Z0 = _z_of_N(m, freq, N0, dn)
    s0 = (Z0, 1.0, 0.0)
    span = 1.0
    res = solve_ivp(rhs, (N0, N0 + span), s0, method='DOP853', rtol=1e-12,
                    atol=[1e-14, 1e-20, 1e-20], dense_output=True)
    Ns = np.linspace(N0, N0 + span, 20000)
    S = res.sol(Ns)
    z, xh, yh = S
    # Predicted accumulated phase (in cycles): integral e^z dN / 2pi.
    dN = Ns[1] - Ns[0]
    predicted_cycles = np.sum(np.exp(z)) * dN / (2.0 * math.pi)
    # Measured sign changes of the carrier oscillation (2 crossings per cycle).
    signs = np.sign(xh)
    crossings = int(np.sum(signs[1:] != signs[:-1]))
    measured_cycles = crossings / 2.0
    assert measured_cycles > 3.0
    rel = abs(measured_cycles - predicted_cycles) / predicted_cycles
    assert rel < 0.1
    assert np.all(np.isfinite(z)) and z.max() > 3.9


@pytest.mark.slow
def test_superhorizon_mode_frozen_until_today():
    """A very-low-frequency mode should stay super-horizon to today (no tail)."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    dn = 2.24e-3
    sol = REF.solve_reference_mode(m, -18.4, dn, z_tail=5.0, rtol=1e-11)
    assert sol['used_tail'] is False
    assert math.isfinite(sol['Ogw_today'])
    # Physical GW density is Ogw - Oj (the shear/j bookkeeping term is
    # negative at super-horizon, so Ogw itself need not be positive).
    assert sol['Ogw_today'] - sol['Oj_today'] > 0.0
    assert sol['Oj_today'] < 0.0  # super-horizon negative Omega_j contribution


@pytest.mark.slow
def test_energy_consistency_dn_gw_definition():
    """DN_gw = Neff0 * g2 / Omega_nu and kappa_r are self-consistent."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    ref = REF.run_reference(m, dn_eff=2.24e-3, freq_subset=[-14.0, -8.0, 0.0],
                            z_tail=5.0, rtol=1e-11, self_consistent=False)
    Omega_nu = gp.Omega_nh2 / m.derived_param['h'] ** 2
    assert ref['DN_gw'] == pytest.approx(gp.Neff0 * ref['g2'] / Omega_nu, rel=1e-12)
    assert ref['kappa_r'] == pytest.approx(
        ref['DN_eff'] * 7 / 8 * (4 / 11) ** (4 / 3) * gp.z_ratio ** 4, rel=1e-12)


def test_float_robustness_extreme_params():
    """Extreme parameters must not produce NaN/inf in the fast spectrum.

    Covers tiny tensor-to-scalar ratio, negligible stiff density, low reheating
    temperature and the deterministic >5 Delta N_eff guard.  The guard rejection
    is a shared fast/LSODA physical rejection, not a numerical failure.
    """
    from stiffgwpy_fast import fast_sgwb as FS

    saved = FS.get_settings()
    try:
        FS.apply_accuracy_mode('production')
        ok = [
            dict(r=1e-8, cr=1, T_re=2e3, kappa10=1e-2),
            dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-8),
            dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
        ]
        for kw in ok:
            m = LCDM_SG(**kw)
            r = FS.SGWB_iter_fast(m, tol=1e-7)
            assert r is not None
            lo = np.asarray(m.log10OmegaGW)
            assert np.all(np.isfinite(lo))
            assert lo.min() > -320.0  # no underflow to -inf
        # Huge stiff drives total N_eff > 5 -> deterministic physical rejection.
        m_bad = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e3)
        r = FS.SGWB_iter_fast(m_bad, tol=1e-7)
        assert r is None
        assert m_bad.fast_failure_reason == 'shared_Neff_guard'
    finally:
        FS.set_threads(saved['threads'])
        FS.set_col_step(saved['col_step'])
        FS.set_h(saved['h'])
        FS.set_z_tail(saved['z_tail'])


def test_sigma_exact_reduces_model_bias():
    """With the continuous present-day anchor, fast Delta N_eff is close to the
    continuous-sigma reference (the dominant historical error was the
    grid-quantised present-day anchor, not the sigma kink)."""
    from stiffgwpy_fast import fast_sgwb as FS

    saved = FS.get_settings()
    try:
        FS.apply_accuracy_mode('production')
        FS.set_z_tail(5.0)
        m_grid = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
        FS.SGWB_iter_fast(m_grid, tol=1e-7)
        m_exact = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
        FS.SGWB_iter_fast(m_exact, tol=1e-7, sigma_exact=True)
        ref = 0.00227081
        d_grid = m_grid.cosmo_param['DN_eff']
        d_exact = m_exact.cosmo_param['DN_eff']
        # Both engines are now within 0.5% of the continuous-sigma reference;
        # the previously-dominant grid-anchor bias is removed.
        assert abs(d_grid - ref) / ref < 5e-3
        assert abs(d_exact - ref) / ref < 5e-3
    finally:
        FS.set_threads(saved['threads'])
        FS.set_col_step(saved['col_step'])
        FS.set_h(saved['h'])
        FS.set_z_tail(saved['z_tail'])


def test_stiff_enhances_high_frequency_spectrum():
    """The stiff era boosts high-frequency Omega_GW (kappa10 -> 0 recovers
    the standard LCDM tensor spectrum, which is much smaller at high f)."""
    for kappa10 in (1.0e-2, 1.0e-6):
        m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=kappa10)
        dn = 2.24e-3
        sol = REF.solve_reference_mode(m, 0.0, dn, z_tail=5.0, rtol=1e-10)
        if kappa10 == 1.0e-2:
            high_s = sol['Ogw_today'] - sol['Oj_today']
        else:
            high_n = sol['Ogw_today'] - sol['Oj_today']
    assert high_s > 2.0 * high_n


def test_no_stiff_low_frequency_limit_similar():
    """Below the stiff-affected band the spectrum is (nearly) kappa10-independent."""
    vals = []
    for kappa10 in (1.0e-2, 1.0e-6):
        m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=kappa10)
        dn = 2.24e-3
        sol = REF.solve_reference_mode(m, -14.0, dn, z_tail=5.0, rtol=1e-10)
        vals.append(math.log10(max(sol['Ogw_today'] - sol['Oj_today'], 1e-300)))
    assert abs(vals[0] - vals[1]) < 0.5


def test_rd_plateau_scaling_slope_equals_nt():
    """In the pure-radiation sub-horizon plateau, Omega_GW ~ f^{n_t}."""
    m = LCDM_SG(r=1e-2, cr=0, T_re=2e3, kappa10=0.0, n_t=-0.4, DN_re=0)
    dn = 2.24e-3
    freqs = np.array([-6.5, -6.25, -6.0, -5.75, -5.5])
    Ogw, Oj, _, _ = REF.spectrum_reference(m, freqs, dn, z_tail=5.0,
                                           rtol=1e-10, workers=1)
    lo = np.log10(Ogw - Oj)
    # Linear least-squares slope of log10(Omega) vs log10(f): Omega ~ f^{n_t}.
    slope = np.polyfit(freqs, lo, 1)[0]
    assert abs(slope - (-0.4)) / 0.4 < 0.1


def test_stiff_plateau_slope_unity():
    """In the stiff-enhanced band the SGWB slope is Omega_GW ~ f^{1}."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1.0)
    dn = 2.24e-3
    freqs = np.array([-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0])
    Ogw, Oj, _, _ = REF.spectrum_reference(m, freqs, dn, z_tail=5.0,
                                           rtol=1e-10, workers=1)
    lo = np.log10(Ogw - Oj)
    slope = np.polyfit(freqs, lo, 1)[0]
    assert abs(slope - 1.0) < 0.05
