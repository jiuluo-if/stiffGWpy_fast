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
    # 从 z=3（已深度进入视界）开始，向后积分若干个 e-fold。
    N0 = None
    # z increases with N; bisection for z(N) = 3.
    # z 随 N 增大；用二分法求 z(N)=3 的起点。
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
    # 理论累计相位（周期数）为积分 ∫e^z dN/(2π)。
    dN = Ns[1] - Ns[0]
    predicted_cycles = np.sum(np.exp(z)) * dN / (2.0 * math.pi)
    # Measured sign changes of the carrier oscillation (2 crossings per cycle).
    # 通过载波振荡的符号变化估计周期数；每个周期有两次过零。
    signs = np.sign(xh)
    crossings = int(np.sum(signs[1:] != signs[:-1]))
    measured_cycles = crossings / 2.0
    assert measured_cycles > 3.0
    rel = abs(measured_cycles - predicted_cycles) / predicted_cycles
    assert rel < 0.1
    assert np.all(np.isfinite(z)) and z.max() > 3.9


def test_negative_inflation_reheating_efolds_report_actual_condition(capsys):
    """A negative inferred inflation interval should not claim a false energy comparison."""
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e40)

    assert model.derived_param['N_inf'] is None
    output = capsys.readouterr().out
    assert 'negative (Delta_N < 0)' in output
    assert 'V_inf smaller than rho' not in output


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
    # 物理引力波密度为 Ogw-Oj；超视界时剪切记账项为负，因此 Ogw 本身不一定为正。
    assert sol['Ogw_today'] - sol['Oj_today'] > 0.0
    assert sol['Oj_today'] < 0.0  # 超视界阶段 Omega_j 的负贡献。


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
    temperature, plus the opt-in total Delta N_eff guard behavior.
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
            assert lo.min() > -320.0  # 确认数值没有下溢为 -inf。
        # The high-level fast path has the physical guard disabled by default.
        m_bad = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e3)
        r = FS.SGWB_iter_fast(m_bad, tol=1e-7)
        assert r is m_bad
        assert m_bad.cosmo_param['DN_eff'] > 5.0

        # Callers can still opt into the previous physical rejection.
        m_guarded = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e3)
        r = FS.SGWB_iter_fast(m_guarded, tol=1e-7, physical_guard=True)
        assert r is None
        assert m_guarded.fast_failure_reason == 'shared_Neff_guard'
    finally:
        FS.set_threads(saved['threads'])
        FS.set_col_step(saved['col_step'])
        FS.set_h(saved['h'])
        FS.set_z_tail(saved['z_tail'])
        FS.set_phase_max(saved['phase_max'])
        FS.set_freq_grid(saved['freq_grid'])
        FS._KINK_SPLIT = saved['kink_split']


def test_gen_fast_avoids_present_day_cancellation_for_large_dn_eff():
    """Large extra radiation must not cancel the finite present-day H2 anchor."""
    from stiffgwpy_fast import fast_sgwb as FS
    from stiffgwpy_fast.exact_background import H2_vec, sigma_vec

    model = LCDM_SG(
        r=0.026636847157969202,
        n_t=2.9357513057156224,
        kappa10=6.370761587894848,
        T_re=8349.551883276852,
        DN_re=5.104113570393838,
        Omega_bh2=0.02493600106013432,
        Omega_ch2=0.3458984920924526,
        H0=52.30419572699695,
        A_s=2.190471632694714e-9,
        cr=0.0,
        DN_eff=2.0648360336399695e53,
    )

    FS.gen_fast(model, h=0.005)

    assert np.all(np.isfinite(model.f_hor))
    assert np.all(np.isfinite(model.sigma))
    today = np.array([model.derived_param['N_inf']])
    h2_today = H2_vec(today, model, model.cosmo_param['DN_eff'])[0]
    sigma_today = sigma_vec(today, model, model.cosmo_param['DN_eff'])[0]
    assert h2_today == pytest.approx(model.derived_param['h'] ** 2, rel=1e-14)
    assert np.isfinite(sigma_today)


def test_fast_guard_off_reports_nonfinite_background_instead_of_index_error():
    """Runaway parameters must fail as a recorded background error, not leak IndexError."""
    model = LCDM_SG(
        r=2.981401117884715e-6,
        n_t=2.433506923077321,
        kappa10=0.005956142774765002,
        T_re=152.38081588419843,
        DN_re=22.053096877437447,
        Omega_bh2=0.014198405074411415,
        Omega_ch2=0.3483299468903473,
        H0=69.1484144351007,
        A_s=3.2897061708711526e-9,
        cr=0.0,
        DN_eff=0.0,
    )

    result = model.SGWB_iter(
        engine='fast', fallback=False, tol=1e-6, z_tail=5.0,
        freq_res=1.0, threads=2, accuracy_mode='fast', physical_guard=False,
    )

    assert result is None
    assert model.fast_failure_reason == 'nonfinite_background'


def test_sigma_exact_reduces_model_bias():
    """Continuous-sigma construction remains closer to the reference oracle."""
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
        assert abs(d_exact - ref) / ref < 5e-3
        assert abs(d_exact - ref) < abs(d_grid - ref)
    finally:
        FS.set_threads(saved['threads'])
        FS.set_col_step(saved['col_step'])
        FS.set_h(saved['h'])
        FS.set_z_tail(saved['z_tail'])
        FS.set_phase_max(saved['phase_max'])
        FS.set_freq_grid(saved['freq_grid'])
        FS._KINK_SPLIT = saved['kink_split']


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
    # 对 log10(Omega) 与 log10(f) 做线性最小二乘，斜率对应 Omega∝f^n_t。
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


def test_fast_outer_reuse_recomputes_when_frequency_grid_changes():
    from stiffgwpy_fast import fast_sgwb as FS

    FS.apply_accuracy_mode('fast')
    model = LCDM_SG(
        r=2.5157032591313454e-5,
        n_t=1.3755212701484574,
        kappa10=0.1653015270125515,
        T_re=0.0485083021944659,
        DN_re=10.0,
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        A_s=2.100549e-9,
        cr=0.0,
        DN_eff=0.0,
    )

    result = model.SGWB_iter(
        engine='fast', fallback=False, tol=1e-6, z_tail=5.0,
        freq_res=1.0, threads=2, accuracy_mode='fast', physical_guard=False,
    )

    assert result is model
    assert np.isfinite(model.log10OmegaGW).all()
    assert np.all(model.Ogw_today - model.Oj_today > 0.0)
    assert model.DN_gw[-1] > 0.0

def test_fast_reports_nonpositive_background_h2_without_leaking_division_error():
    model = LCDM_SG(
        Omega_bh2=-1e100,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-2,
        n_t=0.5,
        cr=0.0,
        T_re=2e3,
        DN_re=10.0,
        kappa10=1e-2,
    )

    result = model.SGWB_iter(
        engine='fast', fallback=False, tol=1e-6, z_tail=5.0,
        freq_res=1.0, threads=2, accuracy_mode='fast', physical_guard=False,
    )

    assert result is None
    assert model.fast_failure_reason == 'nonpositive_h2'
