# -*- coding: utf-8 -*-
"""Tests for the curvature-adaptive frequency refinement."""

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from stiffgwpy_fast import exact_background as EB
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import freq_adaptive as FA
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def _sharp(x):
    """Steep step-like feature around x=0.5 (the analogue of a spectral knee)."""
    return -10.0 * np.tanh((np.asarray(x, dtype=float) - 0.5) / 0.02)


def test_adaptive_refines_sharp_feature():
    """The adaptive grid keeps the PCHIP interpolant below the target error and
    concentrates points around the feature."""
    target = 5e-2
    logf, vals, _ = FA.adapt_refine_grid(np.linspace(0.0, 1.0, 12), _sharp,
                                         target_dex=target, min_dlogf=1e-3,
                                         max_iter=12)
    # Interpolation error on a dense probe grid.
    xp = np.linspace(0.0, 1.0, 5000)
    spl = PchipInterpolator(logf, vals)
    err = np.max(np.abs(spl(xp) - _sharp(xp)))
    assert err < 2.0 * target
    # Points cluster around the feature (0.5), not in the flat outer regions.
    density_mid = np.sum((logf > 0.4) & (logf < 0.6)) / 0.2
    density_edge = np.sum((logf > 0.0) & (logf < 0.2)) / 0.2
    assert density_mid > 2.0 * density_edge


def test_adaptive_is_cheap_on_smooth_function():
    """A smooth quadratic needs almost no refinement beyond the seed grid."""
    f0 = np.linspace(0.0, 1.0, 20)
    logf, vals, n_solves = FA.adapt_refine_grid(f0,
                                                lambda x: (x - 0.5) ** 2,
                                                target_dex=1e-2, min_dlogf=1e-3)
    # A pure quadratic has constant second derivative; a handful of refinements
    # should suffice and it should not blow up to thousands of points.
    assert logf.size <= 4 * f0.size


def test_grid_independent_freqs_invariant_to_sigma_grid():
    """The frequency set must not depend on the sigma-grid resolution."""
    FS.apply_accuracy_mode('production')
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.gen_fast(m, 0.01)
    gf, _, _ = FA.grid_independent_freqs(m, 1.0)
    m2 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    EB.build_transition_grid(m2, 0.01)
    gf2, _, _ = FA.grid_independent_freqs(m2, 1.0)
    assert gf.shape == gf2.shape
    assert np.allclose(gf, gf2, atol=0.0)
    # The independent grid is down to the CMB pivot and covers the UV cutoff.
    assert gf[0] > gf[-1]
    assert gf.size > 100


def test_breakpoint_phi_s2_is_accurate_without_dense_subgrid():
    """The Phase-A primitive uses the breakpoint grid without a global subgrid."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    EB.build_transition_grid(m, 0.02)
    fast = EB.exact_phi_s2_breakpoint(m, m.Nv, m.cosmo_param['DN_eff'])
    dense = EB.exact_phi_s2_grid(m, m.Nv, m.cosmo_param['DN_eff'])
    assert np.max(np.abs(fast[0] - dense[0])) < 1.5e-3
    assert np.max(np.abs(fast[1] - dense[1])) < 5e-3
    assert np.allclose(fast[2], dense[2], rtol=2e-3, atol=1e-8)
    assert np.allclose(fast[3], dense[3], rtol=2e-3, atol=1e-6)


def test_split_primitive_uses_left_limit_at_reheating_kink():
    """The pre-reheating Simpson panel must use sigma=1 at its right limit."""
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    h = 0.02
    n_re = float(m.derived_param['N_inf'] - m.derived_param['N_re'])
    Nv = np.arange(0.0, float(m.derived_param['N_inf']) + 1e-12, h)
    split = EB.exact_phi_s2_split(m, Nv, m.cosmo_param['DN_eff'])
    idx = int(np.searchsorted(Nv, n_re, side='right') - 1)
    left = float(Nv[idx])
    left_h = n_re - left
    sig_left = float(EB.sigma_vec(np.array([left]), m, m.cosmo_param['DN_eff'])[0])
    sig_mid = float(EB.sigma_vec(np.array([(left + n_re) / 2.0]), m,
                                 m.cosmo_param['DN_eff'])[0])
    integrals = []
    sig_nodes = EB.sigma_vec(Nv, m, m.cosmo_param['DN_eff'])
    sig_mid_all = EB.sigma_vec((Nv[:-1] + Nv[1:]) / 2.0, m,
                               m.cosmo_param['DN_eff'])
    integrals.extend(h * (sig_nodes[:-1] + 4.0 * sig_mid_all + sig_nodes[1:]) / 6.0)
    expected_f_re = (np.sum(integrals[:idx])
                     + left_h * (sig_left + 4.0 * sig_mid + 1.0) / 6.0)
    expected_phi_re = 1.5 * expected_f_re - n_re + Nv[0]
    assert split[5] == pytest.approx((n_re - left) / h)
    assert split[6] == pytest.approx(expected_phi_re, rel=2e-12, abs=2e-12)
