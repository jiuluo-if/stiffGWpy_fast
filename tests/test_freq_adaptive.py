# -*- coding: utf-8 -*-
"""Tests for the curvature-adaptive frequency refinement."""

import numpy as np
from scipy.interpolate import PchipInterpolator

from stiffgwpy import exact_background as EB
from stiffgwpy import fast_sgwb as FS
from stiffgwpy import freq_adaptive as FA
from stiffgwpy.stiff_SGWB import LCDM_SG


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
