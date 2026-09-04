# -*- coding: utf-8 -*-
"""freq_adaptive.py -- curvature-adaptive frequency sampling for Omega_GW.

The default empirical frequency grid under-samples the steep low-frequency tail
and the spectral knee (a ~0.16 dex error at f ~ 1e-18 Hz on the default point).
Rather than a hand-tuned uniform grid, this module refines the ``log10 f`` grid
where the spectrum is locally curved, so that the pointwise interpolated
``Omega_GW(f)`` stays under a target dex error while avoiding redundant solves
in smooth regions.

The refinement is engine-agnostic: ``evaluate`` (a callable mapping a batch of
``log10 f`` values to ``log10 Omega_GW``) is supplied by the caller (the
reference pipeline below, or the fast solver later).  New points are solved
incrementally; existing ones are never re-solved.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ['adapt_refine_grid', 'adaptive_spectrum_reference',
           'grid_independent_freqs']


def _local_curvature(x, y):
    """Second derivative of the PCHIP interpolant on the grid points."""
    spl = PchipInterpolator(x, y)
    return np.asarray(spl(x, 2))


def _merge_cached(grid, values, new_grid, new_values):
    """Merge (grid, values) with (new_grid, new_values), re-sorting."""
    g = np.concatenate((grid, new_grid))
    v = np.concatenate((values, new_values))
    order = np.argsort(g)
    return np.asarray(g)[order], np.asarray(v)[order]


def adapt_refine_grid(logf_init, evaluate, target_dex=1e-3, min_dlogf=0.02,
                      max_iter=10, max_points=3000, return_detail=False):
    """Refine a ``log10 f`` grid to keep the PCHIP interpolant under ``target_dex``.

    ``evaluate`` maps an array of ``log10 f`` to ``log10 Omega_GW``.  New
    midpoints are solved incrementally; existing points are cached.

    Returns ``(logf, logOmega, n_solves)`` (or a detail dict with
    ``local_error_max`` when ``return_detail`` is set).
    """
    logf = np.unique(np.asarray(logf_init, dtype=float))
    vals = None
    n_solves = 0
    local_err = 0.0
    for _ in range(max_iter):
        if vals is None:
            vals = np.asarray(evaluate(logf), dtype=float)
            n_solves += int(logf.size)
        c2 = _local_curvature(logf, vals)
        to_add = []
        local_err = 0.0
        for i in range(logf.size - 1):
            h = float(logf[i + 1] - logf[i])
            if h <= min_dlogf:
                continue
            c = abs(float(c2[i]))
            err = c * h * h / 8.0
            local_err = max(local_err, err)
            if err > target_dex:
                to_add.append(0.5 * (logf[i] + logf[i + 1]))
        if not to_add or logf.size >= max_points:
            break
        new = np.unique(np.asarray(to_add, dtype=float))
        nv = np.asarray(evaluate(new), dtype=float)
        n_solves += int(new.size)
        logf, vals = _merge_cached(logf, vals, new, nv)
    if return_detail:
        return logf, vals, n_solves, {'local_error_max': float(local_err),
                                      'n_grid': int(logf.size)}
    return logf, vals, n_solves


def adaptive_spectrum_reference(m, dn_eff, fmin, fmax, target_dex=1e-3,
                                z_tail=5.0, rtol=1e-11, workers=8,
                                init_n=48, **kw):
    """Curvature-adaptive Omega_GW spectrum using the high-accuracy reference."""
    from . import reference as REF

    logf_init = np.linspace(fmin, fmax, init_n)

    def evaluate(grid):
        grid = np.asarray(grid, dtype=float)
        Ogw, Oj, Opgw, used = REF.spectrum_reference(m, grid, dn_eff,
                                                     z_tail=z_tail, rtol=rtol,
                                                     workers=workers)
        return np.log10(Ogw - Oj)

    logf, lo, n_solves = adapt_refine_grid(logf_init, evaluate,
                                           target_dex=target_dex, **kw)
    return logf, lo, n_solves


def grid_independent_freqs(m, freq_res=1.0):
    """A log-frequency grid built only from continuous background quantities.

    ``construct_f`` reads ``m.f_hor`` (the *grid* array) so its sampling shifts
    whenever the sigma-grid resolution changes, which pollutes the bolometric
    integral (see docs/audit_reference.md §7.2).  This builder instead derives
    ``fmax`` (horizon at inflation start), ``fmin`` (horizon today) and the
    reheating feature directly from the continuous background, so the frequency
    set is invariant to the sigma-grid resolution.
    """
    import math as _m

    from . import global_param as gp
    from .exact_background import H2_vec

    d = m.derived_param
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    dn = m.cosmo_param['DN_eff']

    def f_hor_cont(N):
        # Continuous log10(aH/(2pi)/Hz) at absolute N (grid-independent).
        N = float(N)
        H2 = float(H2_vec(np.array([N]), m, dn)[0])
        raw = -0.5 * N + 0.5 * _m.log(H2)
        H2_last = float(H2_vec(np.array([N_inf]), m, dn)[0])
        raw_last = -0.5 * N_inf + 0.5 * _m.log(H2_last)
        H2_re = float(H2_vec(np.array([N_re_abs]), m, dn)[0])
        raw_re = -0.5 * N_re_abs + 0.5 * _m.log(H2_re)
        if N < N_re_abs:
            raw = raw_re - 0.5 * (N - N_re_abs)
        Delta_f = _m.log(2.0 * _m.pi / d['H_0'])
        return (raw - raw_last - Delta_f) / _m.log(10.0)

    fmax = f_hor_cont(0.0)
    fmin = f_hor_cont(N_inf)
    fcmb = _m.log10(gp.f_piv)
    if d['nt'] > 0:
        fmax = min(fmax, (-_m.log10(d['A_t'])) / d['nt'] + _m.log10(gp.f_piv))
    # Uniform-in-log10 frequency grid scaled by freq_res (density tuned to be
    # comparable to the model's empirical grid at freq_res=1), capped at the
    # CMB pivot.  Chosen so freq_res=1 gives ~246 points (same order as
    # construct_f) while staying grid-independent.
    n = max(2, int((fmax - fcmb) * freq_res * 10)) + int((fcmb - fmin) * freq_res * 5)
    logf = np.linspace(fmax, fmin, n + 1)
    return logf, fmin, fmax
