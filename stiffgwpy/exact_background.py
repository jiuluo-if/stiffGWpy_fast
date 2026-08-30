# -*- coding: utf-8 -*-
"""exact_background.py -- continuous-sigma expansion integrals for the fast path.

The fixed-step fast solver builds ``F = integral(sigma dN)`` (and hence the
rescaled ``Phi``/``S2`` amplitudes) from a cubic spline of ``sigma`` on a
uniform grid.  That spline smooths the instantaneous-reheating ``sigma`` kink
(``sigma=1`` -> ``4/3`` or ``2``), producing the ~1% ``model_bias`` documented in
``docs/audit_reference.md``.

This module recomputes ``F`` from the *continuous piecewise-exact* ``sigma``
evaluator, treating the reheating boundary as an exact breakpoint, so the kink
bias is removed.  The callers (``stiffgwpy.fast_sgwb``) can then hand the exact
``Phi``/``Phi_mid``/``S2``/``S2inv`` arrays to the existing numba stepping kernel
without changing the ODE integration itself.
"""

import math

import numpy as np

from . import global_param as gp

__all__ = ['sigma_vec', 'H2_vec', 'exact_phi_s2', 'build_transition_grid',
           'build_kink_refined_grid', 'exact_phi_s2_grid']

ln10 = math.log(10.0)


def _fd_from_ref():
    from .reference import _FD_P_CS, _FD_RHO_CS
    return _FD_RHO_CS, _FD_P_CS


def sigma_vec(N, m, DN_eff):
    """Vectorized exact ``sigma`` at an array of absolute ``N`` values."""
    N = np.asarray(N, dtype=float)
    d = m.derived_param
    Omh2 = d['Omega_mh2']
    Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
    Otrh2 = gp.Omega_orh2 + Oerh2
    Otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + Oerh2
    OLh2 = d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 \
        - gp.Omega_ph2 - Oerh2 - Osh2
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    fd_rho, fd_p = _fd_from_ref()

    sigma = np.ones_like(N)
    running = N >= N_re_abs
    Nr = N[running]
    eN = np.exp(N_inf - Nr)
    e3N = eN * eN * eN
    nu = gp.nu_today / eN
    Omega_ph2 = gp.Omega_ph2
    Omega_nh2 = gp.Omega_nh2
    Omega_mnuh2 = gp.Omega_mnuh2

    m1 = nu > 100.0
    out = np.empty(running.sum())
    if m1.any():
        H2 = (Omh2 + Omega_mnuh2 + (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN[m1]
              + Osh2 * e3N[m1] + OLh2 / e3N[m1])
        out[m1] = (Omh2 + Omega_mnuh2
                   + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN[m1]
                   + 2.0 * Osh2 * e3N[m1]) / H2
    m2 = (~m1) & (nu >= 0.1)
    if m2.any():
        lognu = np.log10(nu[m2])
        rho_nu = fd_rho(lognu)
        p_nu = fd_p(lognu)
        H2 = (Omh2 + (Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * Omega_nh2 + Oerh2) * eN[m2]
              + Osh2 * e3N[m2] + OLh2 / e3N[m2])
        out[m2] = (Omh2
                   + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN[m2]
                   + (rho_nu + p_nu) * Omega_nh2 / 3.0 * eN[m2]
                   + 2.0 * Osh2 * e3N[m2]) / H2
    m3 = (~m1) & (~m2) & (Nr > N_inf - gp.N_fin)
    if m3.any():
        H2 = Omh2 + Otrh2 * eN[m3] + Osh2 * e3N[m3] + OLh2 / e3N[m3]
        out[m3] = (Omh2 + 4.0 / 3.0 * Otrh2 * eN[m3] + 2.0 * Osh2 * e3N[m3]) / H2
    m4 = (~m1) & (~m2) & (~m3) & (Nr >= N_inf - gp.N_max)
    if m4.any():
        Nl = N_inf - Nr[m4]
        rho_i = gp.spl_rho(Nl)
        rhop_i = gp.spl_rhop(Nl)
        H2 = Omh2 + (Omega_ph2 * rho_i + Oerh2) * eN[m4] + Osh2 * e3N[m4] + OLh2 / e3N[m4]
        out[m4] = (Omh2 + (Omega_ph2 * rhop_i + 4.0 / 3.0 * Oerh2) * eN[m4]
                   + 2.0 * Osh2 * e3N[m4]) / H2
    m5 = (~m1) & (~m2) & (~m3) & (~m4)
    if m5.any():
        H2 = Omh2 + Otreh2 * eN[m5] + Osh2 * e3N[m5] + OLh2 / e3N[m5]
        out[m5] = (Omh2 + 4.0 / 3.0 * Otreh2 * eN[m5] + 2.0 * Osh2 * e3N[m5]) / H2
    sigma[running] = out
    return sigma


def H2_vec(N, m, DN_eff):
    """Vectorized exact ``H^2`` at an array of absolute ``N`` (mirrors sigma_vec)."""
    N = np.asarray(N, dtype=float)
    d = m.derived_param
    Omh2 = d['Omega_mh2']
    Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
    Otrh2 = gp.Omega_orh2 + Oerh2
    Otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + Oerh2
    OLh2 = d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 \
        - gp.Omega_ph2 - Oerh2 - Osh2
    N_inf = d['N_inf']
    fd_rho, _ = _fd_from_ref()
    out = np.empty_like(N)
    for i in range(N.size):
        Ni = float(N[i])
        eN = math.exp(N_inf - Ni)
        e3N = eN * eN * eN
        nu = gp.nu_today / eN
        if nu > 100.0:
            out[i] = (Omh2 + gp.Omega_mnuh2
                      + (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2 + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        elif nu >= 0.1:
            lognu = math.log10(nu)
            rho_nu = float(fd_rho(lognu))
            out[i] = (Omh2
                      + (gp.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * gp.Omega_nh2 + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        elif Ni > N_inf - gp.N_fin:
            out[i] = Omh2 + Otrh2 * eN + Osh2 * e3N + OLh2 / e3N
        elif Ni >= N_inf - gp.N_max:
            Nl = N_inf - Ni
            out[i] = (Omh2 + (gp.Omega_ph2 * float(gp.spl_rho(Nl)) + Oerh2) * eN
                      + Osh2 * e3N + OLh2 / e3N)
        else:
            out[i] = Omh2 + Otreh2 * eN + Osh2 * e3N + OLh2 / e3N
    return out


def exact_phi_s2(m, Nv, DN_eff, h):
    """Exact ``Phi``/``Phi_mid``/``S2``/``S2inv`` from the continuous sigma.

    ``F = integral sigma dN`` is evaluated on a fine sub-grid with the reheating
    boundary as an exact breakpoint, so the kink is not smoothed.  Returns the
    arrays expected by ``fast_sgwb.solve_kernel`` (which uses ``N0=0``).
    """
    Nv = np.asarray(Nv, dtype=float)
    # Fine sub-grid (a few points per base step) plus its midpoints for Phi_mid.
    half = 0.5 * h
    frac = 0.125 * h      # fine sub-grid spacing for the F/Phi integral
    sub = np.concatenate([
        np.arange(Nv[0], Nv[-1] + 1e-12, frac),
        np.arange(Nv[0] + half, Nv[-1] + 1e-12, h),
    ])
    sub = np.unique(np.sort(np.concatenate((sub, Nv))))
    d = m.derived_param
    N_re_abs = d['N_inf'] - d['N_re']
    if (N_re_abs > sub[0]) and (N_re_abs < sub[-1]):
        sub = np.unique(np.concatenate((sub, [N_re_abs])))
    sig = sigma_vec(sub, m, DN_eff)
    # Cumulative integral F on the sub-grid (shape-preserving trapezoid; the kink
    # is a breakpoint, so sigma is smooth within each sub-interval).
    F = np.concatenate(([0.0], np.cumsum(0.5 * (sig[1:] + sig[:-1]) * np.diff(sub))))
    # Interpolate F onto grid nodes and midpoints.
    F_grid = np.interp(Nv, sub, F)
    F_mid = np.interp(Nv + half, sub, F)
    N0 = Nv[0]
    Phi_grid = 1.5 * F_grid - Nv + N0
    Phi_mid = 1.5 * F_mid - (Nv + half) + N0
    Psi = 3.0 * F_grid - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    return Phi_grid.astype(np.float64), Phi_mid.astype(np.float64), \
        S2.astype(np.float64), S2inv.astype(np.float64)


def build_transition_grid(m, h):
    """Variable N-grid with a node exactly at the reheating kink.

    The kink becomes a grid node so that no Magnus step straddles the
    ``sigma`` jump (removing the O(h) frozen-midpoint phase error).  Returns
    ``(Nv, sigma, f_hor, index_re)`` and sets the model grid attributes.
    """
    d = m.derived_param
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    len_inf = math.floor(N_inf / h) + 1
    Nv = np.arange(0, len_inf) * h
    present = Nv[-1]
    Nv = np.unique(np.sort(np.concatenate((Nv, [N_re_abs, present]))))
    dn = m.cosmo_param['DN_eff']
    sigma = sigma_vec(Nv, m, dn)
    # Exact f_hor (log10 aH/(2pi)/Hz) with the MD extension for N < N_re_abs,
    # normalised to the grid's present node (same convention as gen_fast).
    H2 = H2_vec(Nv, m, dn)
    raw = np.empty_like(Nv)
    raw_re = -0.5 * N_re_abs + 0.5 * math.log(float(H2_vec(np.array([N_re_abs]), m, dn)[0]))
    for i in range(Nv.size):
        if Nv[i] < N_re_abs:
            raw[i] = raw_re - 0.5 * (Nv[i] - N_re_abs)
        else:
            raw[i] = -0.5 * Nv[i] + 0.5 * math.log(float(H2[i]))
    Delta_f = math.log(2.0 * math.pi / d['H_0'])
    f_hor = (raw - raw[-1] - Delta_f) / math.log(10.0)
    index_re = int(np.argmin(np.abs(Nv - N_re_abs)))
    m.Nv = Nv
    m.N = Nv - Nv[-1]
    m.sigma = sigma
    m.f_hor = f_hor
    m.f_re = f_hor[index_re]
    return Nv, sigma, f_hor, index_re


def build_kink_refined_grid(m, h, refine_radius=2.0, refine_factor=8):
    """Kink-aware non-uniform grid: uniform ``h`` plus a refined zone around the
    reheating boundary ``N_re_abs`` (step ``h/refine_factor`` within
    ``+/- refine_radius`` e-folds), with the kink as an exact node.

    Returns ``(Nv, sigma, f_hor, index_re)`` and sets the model grid attributes.
    """
    d = m.derived_param
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    len_inf = math.floor(N_inf / h) + 1
    Nv = np.arange(0, len_inf) * h
    present = N_inf
    lo = N_re_abs - refine_radius
    hi = N_re_abs + refine_radius
    if lo < 0.0:
        lo = 0.0
    if hi > present:
        hi = present
    refined = np.arange(lo, hi + 1e-12, h / refine_factor)
    Nv = np.unique(np.sort(np.concatenate((Nv, refined, [N_re_abs, present]))))
    dn = m.cosmo_param['DN_eff']
    sigma = sigma_vec(Nv, m, dn)
    H2 = H2_vec(Nv, m, dn)
    raw = np.empty_like(Nv)
    raw_re = -0.5 * N_re_abs + 0.5 * math.log(float(H2_vec(np.array([N_re_abs]), m, dn)[0]))
    for i in range(Nv.size):
        if Nv[i] < N_re_abs:
            raw[i] = raw_re - 0.5 * (Nv[i] - N_re_abs)
        else:
            raw[i] = -0.5 * Nv[i] + 0.5 * math.log(float(H2[i]))
    Delta_f = math.log(2.0 * math.pi / d['H_0'])
    f_hor = (raw - raw[-1] - Delta_f) / math.log(10.0)
    index_re = int(np.argmin(np.abs(Nv - N_re_abs)))
    m.Nv = Nv
    m.N = Nv - Nv[-1]
    m.sigma = sigma
    m.f_hor = f_hor
    m.f_re = f_hor[index_re]
    return Nv, sigma, f_hor, index_re


def exact_phi_s2_grid(m, Nv, DN_eff):
    """Exact ``Phi``/``midpoint-Phi``/``S2``/``S2inv``/``h_arr`` on a variable grid."""
    Nv = np.asarray(Nv, dtype=float)
    d = m.derived_param
    N_re_abs = d['N_inf'] - d['N_re']
    # Fine sub-grid (roughly half the local spacing) for the F integral, with the
    # kink as a breakpoint.  Use the *median* spacing so that a single inserted
    # kink node cannot collapse the step and blow up the sub-grid.
    h = float(np.median(np.abs(np.diff(Nv))))
    step = 0.5 * h
    sub = np.concatenate((np.arange(Nv[0], Nv[-1] + 1e-12, step), Nv))
    if (N_re_abs > sub[0]) and (N_re_abs < sub[-1]):
        sub = np.concatenate((sub, [N_re_abs]))
    sub = np.unique(np.sort(sub))
    sig = sigma_vec(sub, m, DN_eff)
    F = np.concatenate(([0.0], np.cumsum(0.5 * (sig[1:] + sig[:-1]) * np.diff(sub))))
    mid = 0.5 * (Nv[:-1] + Nv[1:])
    F_nodes = np.interp(Nv, sub, F)
    F_mid = np.interp(mid, sub, F)
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - mid + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    h_arr = np.diff(Nv).astype(np.float64)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64), h_arr)
