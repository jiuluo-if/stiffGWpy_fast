# -*- coding: utf-8 -*-
"""
fast_sgwb.py -- experimental approximate fast solver for LCDM_SG.SGWB_iter().

The original SGWB_iter() solves the tensor-mode Boltzmann equations with
scipy.integrate.solve_ivp (LSODA) per frequency channel and integrates the
resulting spectrum with scipy.integrate.simpson, repeating both inside a
bisection loop on Delta N_eff.  A typical call chain takes ~7-24 s.

This module implements a *different, approximate* numerical scheme for the
same physical equations and the same outer bisection target:

  * numba JIT kernels for the expansion history and the ODE stepping,
  * a fixed-step analytic-rotation (Magnus-type) solver (h = 0.01) instead of
    the adaptive LSODA solver,
  * an analytic deep-subhorizon tail beyond z = 5 (the original code uses
    such a tail as well),
  * a precomputed Simpson weight matrix instead of per-column scipy calls,
  * PCHIP refinement of the bolometric integrals onto the fine N grid,
  * OpenMP parallelism over frequency channels.

Because the ODE solver and the time-column integration scheme differ from the
original, results are close but NOT bit-identical to SGWB_iter(): on the
12-case spot validation the final Delta N_eff agrees to ~5e-5-8e-5 relative,
the spectrum agrees to ~4e-4 dex (linear-Omega relative difference
~8e-4-1e-3), while the full DN_gw(N) evolution can differ by up to ~1%-37%
in the early near-zero region.  Treat it as an experimental fast solver;
keep the LSODA path for cross-checks and fallback.

Usage
-----
    from stiffgwpy import LCDM_SG
    from stiffgwpy import fast_sgwb

    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    fast_sgwb.SGWB_iter_fast(m)     # fills the fast-path output attributes
    # or, to keep a single API with automatic LSODA fallback:
    m.SGWB_iter(engine='fast', fallback=True)

Optional tuning (read before importing this module):
    os.environ['FAST_THREADS']  = '8'    # OpenMP threads; default = numba default
    os.environ['FAST_COL_STEP'] = '4'    # output-column stride (1-8, default 4)

The module is deterministic; its numba kernels are cache-compiled on first use.
"""
import math
import os as _os

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads
from scipy import interpolate

from . import global_param as gp
from .functions import int_FD

__all__ = ['SGWB_iter_fast', 'gen_fast', 'set_threads', 'set_col_step', 'set_h',
           'set_z_tail', 'get_settings', 'apply_accuracy_mode', 'ACCURACY_MODES']

# Default OpenMP threads: numba's own default (no more than the detected core
# count).  We do NOT force a fixed number at import time -- that previously
# raised ValueError on machines with fewer than 32 cores -- and only call
# set_num_threads() when FAST_THREADS is explicitly set.
_MAX_THREADS = get_num_threads()
_THREADS = _MAX_THREADS
_fast_threads_env = _os.environ.get('FAST_THREADS')
if _fast_threads_env is not None:
    _THREADS = int(_fast_threads_env)
    if not 1 <= _THREADS <= _MAX_THREADS:
        raise ValueError('FAST_THREADS must be an integer in [1, %d], got %r'
                         % (_MAX_THREADS, _fast_threads_env))
    set_num_threads(_THREADS)

_COL_STEP = 4
_fast_col_step_env = _os.environ.get('FAST_COL_STEP')
if _fast_col_step_env is not None:
    _COL_STEP = int(_fast_col_step_env)
    if not 1 <= _COL_STEP <= 8:
        raise ValueError('FAST_COL_STEP must be an integer in [1, 8], got %r'
                         % _fast_col_step_env)

# Fixed step size / Nv grid spacing.  Both the expansion grid and the Magnus
# step use this h; 0.01 is the default.  Tunable for convergence studies.
_FAST_H = 0.01
_fast_h_env = _os.environ.get('FAST_H')
if _fast_h_env is not None:
    _FAST_H = float(_fast_h_env)
    if not 1e-4 <= _FAST_H <= 0.1:
        raise ValueError('FAST_H must be in [1e-4, 0.1], got %r' % _fast_h_env)

# Deep-subhorizon analytic-tail threshold (z = ln(k/aH) at which the numerical
# stepping hands over to the analytic tail).  The original LSODA path uses the
# same z_tail value when configured (see functions.solve_SGWB).
_Z_TAIL = 5.0
_fast_ztail_env = _os.environ.get('FAST_Z_TAIL')
if _fast_ztail_env is not None:
    _Z_TAIL = float(_fast_ztail_env)
    if not 2.0 <= _Z_TAIL <= 15.0:
        raise ValueError('FAST_Z_TAIL must be in [2.0, 15.0], got %r'
                         % _fast_ztail_env)

MAX_ITER = 60            # cap on the outer bisection loop
ln10 = math.log(10.0)


def set_threads(n):
    """Set the number of OpenMP threads used by the frequency-parallel kernel.

    `n` must be an integer in [1, numba's detected thread count].
    """
    global _THREADS
    n = int(n)
    if not 1 <= n <= _MAX_THREADS:
        raise ValueError('thread count must be an integer in [1, %d], got %d'
                         % (_MAX_THREADS, n))
    _THREADS = n
    set_num_threads(_THREADS)


def set_col_step(n):
    """Set the output-column stride (1..8); 4 is a good speed/accuracy trade-off."""
    global _COL_STEP
    n = int(n)
    if not 1 <= n <= 8:
        raise ValueError('column stride must be an integer in [1, 8], got %d' % n)
    _COL_STEP = n


def set_h(h):
    """Set the fixed step size / expansion-grid spacing (1e-4 .. 0.1); 0.01 default."""
    global _FAST_H
    h = float(h)
    if not 1e-4 <= h <= 0.1:
        raise ValueError('step size must be in [1e-4, 0.1], got %r' % h)
    _FAST_H = h


def set_z_tail(z):
    """Set the analytic-tail threshold z_tail (2.0 .. 15.0); 5.0 default."""
    global _Z_TAIL
    z = float(z)
    if not 2.0 <= z <= 15.0:
        raise ValueError('z_tail must be in [2.0, 15.0], got %r' % z)
    _Z_TAIL = z


def get_settings():
    """Snapshot the current fast-solver module settings (threads/col_step/h/z_tail)."""
    return dict(threads=_THREADS, col_step=_COL_STEP, h=_FAST_H, z_tail=_Z_TAIL)


# Named accuracy presets (audit phase "three recommended modes").
# Values come from the phase-2 convergence study: engine-vs-LSODA difference
# at h=0.01 is ~1e-5 on Delta N_eff while the shared sigma-grid bias vs a
# deep reference is ~0.73% (h=0.01) / ~0.33% (h=0.005); z_tail=7 reduces the
# analytic-tail error to ~2e-5 and z_tail=10 to ~2.4e-7; col_step has <1e-9
# effect on the final Delta N_eff (it only shapes early small-value curves);
# freq_res=2.0 halves the low-frequency-tail undersampling error.
ACCURACY_MODES = {
    'reference': dict(h=0.00125, col_step=1, z_tail=10.0, freq_res=2.0,
                      tol=1e-8, threads=8),
    'production': dict(h=0.01, col_step=4, z_tail=7.0, freq_res=1.0,
                       tol=1e-7, threads=8),
    'ultra-fast': dict(h=0.01, col_step=8, z_tail=5.0, freq_res=1.0,
                       tol=1e-6, threads=16),
}


def apply_accuracy_mode(name):
    """Apply a named accuracy preset to the module settings and return its table.

    ``name`` must be a key of :data:`ACCURACY_MODES`.  The preset's
    threads/col_step/h/z_tail are applied through the usual setters (threads is
    clamped to the numba-detected maximum); the preset's freq_res/tol are
    returned in the table so the caller can forward them to
    :func:`SGWB_iter_fast`.  Module settings are process-global, as documented
    for the setters.
    """
    if name not in ACCURACY_MODES:
        raise ValueError('unknown accuracy mode %r; choose from %s'
                         % (name, sorted(ACCURACY_MODES)))
    cfg = dict(ACCURACY_MODES[name])
    set_col_step(cfg['col_step'])
    set_h(cfg['h'])
    set_z_tail(cfg['z_tail'])
    set_threads(min(cfg['threads'], _MAX_THREADS))
    return cfg


# ================= module-level tables (once) =================
def _build_fd_table():
    """Fermi-Dirac rho/p table, loaded from a cached copy when available.

    A precomputed table is shipped inside the package (``fd_table.npz``); when
    that is missing (e.g. an editable install without data files) we fall back
    to the user cache and finally to recomputation, saving the result to the
    user cache so the next import is cheap.
    """
    table_file = None
    pkg_file = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'fd_table.npz')
    if _os.path.exists(pkg_file):
        table_file = pkg_file
    else:
        cache_file = _os.path.join(_os.path.expanduser('~'), '.cache', 'stiffgwpy', 'fd_table.npz')
        if _os.path.exists(cache_file):
            table_file = cache_file
    if table_file is not None:
        try:
            data = np.load(table_file)
            nu, vals = data['nu'], data['vals']
            if nu.shape == (3001,) and vals.shape == (3001, 2):
                return nu, vals
        except (IOError, ValueError, KeyError):
            pass
    nu = np.logspace(-1.0, 2.0, 3001)
    vals = np.array([int_FD(u) for u in nu])
    try:
        cache_dir = _os.path.join(_os.path.expanduser('~'), '.cache', 'stiffgwpy')
        _os.makedirs(cache_dir, exist_ok=True)
        np.savez(_os.path.join(cache_dir, 'fd_table.npz'), nu=nu, vals=vals)
    except OSError:
        pass
    return nu, vals


_FD_NU, _FD_VALS = _build_fd_table()
_FD_RHO = interpolate.CubicSpline(np.log10(_FD_NU), _FD_VALS[:, 0])
_FD_P = interpolate.CubicSpline(np.log10(_FD_NU), _FD_VALS[:, 1])
_FD_X0 = float(np.log10(_FD_NU[0]))
_FD_DX = float(np.log10(_FD_NU[1]) - _FD_X0)
_FD_C_RHO = np.ascontiguousarray(_FD_RHO.c)
_FD_C_P = np.ascontiguousarray(_FD_P.c)
_TH_X = np.ascontiguousarray(gp.spl_rho.x)
_TH_C_RHO = np.ascontiguousarray(gp.spl_rho.c)
_TH_C_RHOP = np.ascontiguousarray(gp.spl_rhop.c)

@njit(cache=True)
def cubic_uniform(x, x0, dx, c):
    nseg = c.shape[1]
    i = int((x - x0)/dx)
    if i < 0: i = 0
    if i > nseg-1: i = nseg-1
    t = x - (x0 + i*dx)
    return ((c[0,i]*t + c[1,i])*t + c[2,i])*t + c[3,i]

@njit(cache=True)
def cubic_any(x, breaks, c):
    n = len(breaks)
    lo = 0; hi = n-1
    while hi - lo > 1:
        mid = (lo + hi)//2
        if breaks[mid] <= x: lo = mid
        else: hi = mid
    nseg = c.shape[1]
    if lo < 0: lo = 0
    if lo > nseg-1: lo = nseg-1
    t = x - breaks[lo]
    return ((c[0,lo]*t + c[1,lo])*t + c[2,lo])*t + c[3,lo]

# ================= gen_expansion in numba =================
@njit(cache=True)
def gen_kernel(Nv, Sv, f_hor, index_re, Omh2, Osh2, Oerh2, Otrh2, Otreh2, OLh2,
               Omega_mnuh2, Omega_ph2, Omega_nh2, nu_today, N_fin, N_max,
               FD_X0, FD_DX, FD_C_RHO, FD_C_P, TH_X, TH_C_RHO, TH_C_RHOP,
               Delta_f, ln10v):
    n = len(Nv)
    Nlast = Nv[n-1]
    for i in range(index_re, n):
        eN = math.exp(Nlast - Nv[i])
        e3N = eN*eN*eN
        nu = nu_today / eN
        if nu > 100.0:
            H2 = Omh2 + Omega_mnuh2 + (Omega_ph2 + 2.0/3.0*Omega_nh2 + Oerh2)*eN + Osh2*e3N + OLh2/e3N
            Sv[i] = (Omh2 + Omega_mnuh2 + 4.0/3.0*(Omega_ph2 + 2.0/3.0*Omega_nh2 + Oerh2)*eN + 2.0*Osh2*e3N)/H2
        elif nu >= 0.1:
            lnnu = math.log10(nu)
            rho_nu = cubic_uniform(lnnu, FD_X0, FD_DX, FD_C_RHO)
            p_nu = cubic_uniform(lnnu, FD_X0, FD_DX, FD_C_P)
            H2 = Omh2 + (Omega_ph2 + (2.0/3.0 + rho_nu/3.0)*Omega_nh2 + Oerh2)*eN + Osh2*e3N + OLh2/e3N
            Sv[i] = (Omh2 + 4.0/3.0*(Omega_ph2 + 2.0/3.0*Omega_nh2 + Oerh2)*eN + (rho_nu + p_nu)*Omega_nh2/3.0*eN + 2.0*Osh2*e3N)/H2
        elif Nv[i] > Nlast - N_fin:
            H2 = Omh2 + Otrh2*eN + Osh2*e3N + OLh2/e3N
            Sv[i] = (Omh2 + 4.0/3.0*Otrh2*eN + 2.0*Osh2*e3N)/H2
        elif Nv[i] >= Nlast - N_max:
            Nl = Nlast - Nv[i]
            rho_i = cubic_any(Nl, TH_X, TH_C_RHO)
            rhop_i = cubic_any(Nl, TH_X, TH_C_RHOP)
            H2 = Omh2 + (Omega_ph2*rho_i + Oerh2)*eN + Osh2*e3N + OLh2/e3N
            Sv[i] = (Omh2 + (Omega_ph2*rhop_i + 4.0/3.0*Oerh2)*eN + 2.0*Osh2*e3N)/H2
        else:
            H2 = Omh2 + Otreh2*eN + Osh2*e3N + OLh2/e3N
            Sv[i] = (Omh2 + 4.0/3.0*Otreh2*eN + 2.0*Osh2*e3N)/H2
        f_hor[i] = -0.5*Nv[i] + 0.5*math.log(H2)
    fv_re = f_hor[index_re]
    Nv_re = Nv[index_re]
    for i in range(index_re):
        Sv[i] = 1.0
        f_hor[i] = fv_re - 0.5*(Nv[i] - Nv_re)
    f0 = f_hor[n-1] + Delta_f
    for i in range(n):
        f_hor[i] = (f_hor[i] - f0)/ln10v

def gen_fast(m, h=0.01):
    d = m.derived_param
    p = m.cosmo_param
    Omh2 = d['Omega_mh2']; Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2*7/8*(4/11)**(4/3)*p['DN_eff']
    Otrh2 = gp.Omega_orh2 + Oerh2
    Otreh2 = gp.Omega_ph2*gp.rho_th[-1] + Oerh2
    OLh2 = d['h']**2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2*2/3 - gp.Omega_ph2 - Oerh2 - Osh2
    len_inf = math.floor(d['N_inf']/h)+1
    Nv = np.arange(0, len_inf)*h
    index_re = len_inf-1 - math.floor(d['N_re']/h)
    Sv = np.empty(len_inf); f_hor = np.empty(len_inf)
    Delta_f = math.log(2*math.pi/d['H_0'])
    gen_kernel(Nv, Sv, f_hor, index_re, Omh2, Osh2, Oerh2, Otrh2, Otreh2, OLh2,
               gp.Omega_mnuh2, gp.Omega_ph2, gp.Omega_nh2, gp.nu_today, gp.N_fin, gp.N_max,
               _FD_X0, _FD_DX, _FD_C_RHO, _FD_C_P, _TH_X, _TH_C_RHO, _TH_C_RHOP,
               Delta_f, ln10)
    m.Nv = Nv; m.N = Nv - Nv[-1]; m.sigma = Sv; m.f_hor = f_hor
    m.f_re = f_hor[index_re]
    return len_inf, index_re

# ================= prep in numba (spline + primitive + phi/psi/s2 + j0s/z0s) =================
@njit(cache=True)
def prep_kernel(Nv, Sv, f_hor, freqs, h, ln10v,
                Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus):
    nv = len(Nv); nseg = nv - 1
    M = np.empty(nv)
    inv_h2 = 1.0/(h*h)
    M[1] = (Sv[0] - 2.0*Sv[1] + Sv[2])*inv_h2
    M[nv-2] = (Sv[nv-3] - 2.0*Sv[nv-2] + Sv[nv-1])*inv_h2
    m = nv - 4
    aa = np.empty(m); bb = np.empty(m); cc = np.empty(m); dd = np.empty(m)
    for i in range(m):
        k = i + 2
        aa[i] = 1.0; bb[i] = 4.0; cc[i] = 1.0
        dd[i] = 6.0*(Sv[k-1] - 2.0*Sv[k] + Sv[k+1])*inv_h2
    dd[0] -= M[1]
    dd[m-1] -= M[nv-2]
    for i in range(1, m):
        w = aa[i]/bb[i-1]
        bb[i] -= w*cc[i-1]
        dd[i] -= w*dd[i-1]
    M[m+1] = dd[m-1]/bb[m-1]
    for i in range(m-2, -1, -1):
        M[i+2] = (dd[i] - cc[i]*M[i+3])/bb[i]
    M[0] = 2.0*M[1] - M[2]
    M[nv-1] = 2.0*M[nv-2] - M[nv-3]
    F = np.empty(nv)
    F[0] = 0.0
    for i in range(nseg):
        Mi = M[i]; Mi1 = M[i+1]
        a = (Mi1 - Mi)/(6.0*h)
        b = Mi/2.0
        c = (Sv[i+1] - Sv[i])/h - h*(2.0*Mi + Mi1)/6.0
        d = Sv[i]
        F[i+1] = F[i] + (a/4.0*h**4 + b/3.0*h**3 + c/2.0*h**2 + d*h)
    N0 = Nv[0]
    for i in range(nv):
        Phi_grid[i] = 1.5*F[i] - Nv[i] + N0
        Psi[i] = 3.0*F[i] - 4.0*Nv[i]
        S2[i] = math.exp(Psi[i])
        S2inv[i] = math.exp(-0.5*Psi[i])
    for i in range(nv):
        if i < nseg:
            seg = i; dx = 0.5*h
        else:
            seg = nseg-1; dx = 1.5*h
        Mi = M[seg]; Mi1 = M[seg+1]
        a = (Mi1 - Mi)/(6.0*h)
        b = Mi/2.0
        c = (Sv[seg+1] - Sv[seg])/h - h*(2.0*Mi + Mi1)/6.0
        d = Sv[seg]
        prim = F[seg] + (a/4.0*dx**4 + b/3.0*dx**3 + c/2.0*dx**2 + d*dx)
        xm = Nv[i] + 0.5*h
        Phi_mid[i] = 1.5*prim - xm + N0
    for j in range(nv):
        fp_minus[j] = math.exp(-f_hor[j]*ln10v)
    Nf = len(freqs)
    for mm in range(Nf):
        freq3 = freqs[mm] + 3.0
        lo = 0; hi = nv
        while lo < hi:
            mid = (lo + hi)//2
            if f_hor[mid] >= freq3:
                lo = mid + 1
            else:
                hi = mid
        j0 = lo - 1
        if j0 < 0: j0 = 0
        if j0 > nv-1: j0 = nv-1
        j0s[mm] = j0
        z0s[mm] = (freqs[mm] - f_hor[j0])*ln10v

def prep_fast(m, Nv, freqs, h):
    Sv = m.sigma; f_hor = m.f_hor
    nv = len(Nv)
    Phi_grid = np.empty(nv); Phi_mid = np.empty(nv); Psi = np.empty(nv)
    S2 = np.empty(nv); S2inv = np.empty(nv); fp_minus = np.empty(nv)
    j0s = np.empty(len(freqs), dtype=np.int64); z0s = np.empty(len(freqs))
    prep_kernel(Nv, Sv, f_hor, freqs, h, ln10, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus)
    return Sv, f_hor, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus

# ================= stepping / assembly / kernel (from proto11) =================
@njit(cache=True)
def scaled_step(xh, yh, z_mid, h):
    w = math.exp(z_mid)
    w2 = w*w
    if w2 >= 1.0:
        Om = math.sqrt(w2 - 1.0)
        c = math.cos(Om*h); si = math.sin(Om*h)/Om
    else:
        Om = math.sqrt(1.0 - w2)
        x = Om*h
        c = 1.0 + 0.5*x*x
        si = h*(1.0 + x*x/6.0)
    return (c-si)*xh - w*si*yh, w*si*xh + (c+si)*yh

@njit(cache=True)
def assemble_main(Ogw, Oj, Opgw, m, slot, s2, xh, yh, zz, Pt):
    ss = math.sqrt(s2)
    P = s2*(xh*xh + yh*yh)
    Th = ss*yh*math.exp(-zz)
    oj = ss*xh*Th/3.0*Pt
    Oj[m,slot] = oj
    Ogw[m,slot] = P/24.0*Pt + oj
    Opgw[m,slot] = s2*(-5.0*xh*xh + 7.0*yh*yh)/72.0*Pt

@njit(cache=True)
def assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus):
    Th = coeff*eNz*ev_minus[kk2]
    xf = Th*fp_i*fp_minus[kk2]
    oj = -Th*Th/3.0*Pt
    Op = xf*xf/36.0*Pt
    Oj[m,slot] = oj
    Opgw[m,slot] = Op
    Ogw[m,slot] = 3.0*Op + oj

@njit(parallel=True, cache=True)
def solve_kernel(Nv, Phi_grid, Phi_mid, S2, S2inv,
                 j0s, z0s, P_t, ev_minus, fp_minus, fp_freq,
                 assemble, n_coarse, col_step, h, z_tail, Ogw, Oj, Opgw):
    nv = len(Nv)
    for m in prange(len(j0s)):
        j0 = j0s[m]; z0 = z0s[m]
        Pt = P_t[m]; fp_i = fp_freq[m]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0)*S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh = 0.0; lyh = yh; last_z = zz
        if k % col_step == 0:
            if assemble:
                slot = k//col_step
                if slot >= n_coarse-1: slot = n_coarse-1
                assemble_main(Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv-1:
            assemble_main(Ogw, Oj, Opgw, m, n_coarse-1, S2[k], xh, yh, zz, Pt)
        if zz < z_tail:
            while k < nv-1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                xh, yh = scaled_step(xh, yh, z0 + Phi_mid[k] - Phi0, h)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0:
                    if assemble:
                        slot = k//col_step
                        if slot >= n_coarse-1: slot = n_coarse-1
                        assemble_main(Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
                elif k == nv-1:
                    assemble_main(Ogw, Oj, Opgw, m, n_coarse-1, S2[k], xh, yh, zz, Pt)
                if zz < z_tail:
                    lxh = xh; lyh = yh; last_z = zz
            if zz < z_tail:
                kend = nv-1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        if kend < nv-1:
            s2k = S2[kend]
            coeff = math.sqrt(0.5*s2k*(lxh*lxh + lyh*lyh))
            eNz = math.exp(Nv[kend] - last_z)
            slot_start = kend//col_step
            if slot_start >= n_coarse-1: slot_start = n_coarse-1
            while slot_start < n_coarse:
                kk2 = col_step*slot_start
                if slot_start == n_coarse-1: kk2 = nv-1
                if kk2 > kend: break
                slot_start += 1
            for slot in range(slot_start, n_coarse):
                if not (assemble or slot == n_coarse-1):
                    continue
                kk2 = col_step*slot
                if slot == n_coarse-1: kk2 = nv-1
                assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus)

# ================= simpson weights (precomputed matrix) =================
@njit(cache=True)
def simpson_row(Xf, h, a, b, W):
    m = b - a + 1
    if m < 2:
        return
    if m == 2:
        W[a] += 0.5*h[a]
        W[a+1] += 0.5*h[a]
        return
    s = a
    while s + 2 <= b:
        h0 = h[s]; h1 = h[s+1]
        hsum = h0+h1; r = h0/h1; t = hsum/6.0
        W[s] += t*(2.0 - 1.0/r)
        W[s+1] += t*hsum*hsum/(h0*h1)
        W[s+2] += t*(2.0 - r)
        s += 2
    if (m % 2 == 0) and (m >= 4):
        i3 = b-2; i2 = b-1; i1 = b
        h0 = h[i3]; h1 = h[i2]
        alpha = (2*h1*h1 + 3*h0*h1)/(6.0*(h0+h1))
        beta = (h1*h1 + 3*h0*h1)/(6.0*h0)
        eta = h1*h1*h1/(6.0*h0*(h0+h1))
        W[i1] += alpha
        W[i2] += beta
        W[i3] -= eta

@njit(cache=True)
def build_Wmat(Nf, Xf, h, Wmat):
    for jh in range(Nf):
        a = Nf-1-jh
        W = np.zeros(Nf)
        simpson_row(Xf, h, a, Nf-1, W)
        for p in range(Nf):
            Wmat[jh, p] = W[p]

@njit(parallel=True, cache=True)
def int_SGWB_W(Nf, n_coarse, j_hi, Wmat, Ogw, Oj, Opgw, g2c, w2c, ln10v):
    for c in prange(n_coarse):
        jh = j_hi[c]
        if jh < 0:
            g2c[c] = 0.0; w2c[c] = 0.0
            continue
        acc_g = 0.0; acc_w = 0.0
        for p in range(Nf):
            jj = Nf-1-p
            acc_g += Wmat[jh, p]*(Ogw[jj, c] - Oj[jj, c])
            acc_w += Wmat[jh, p]*Opgw[jj, c]
        g2c[c] = acc_g*ln10v
        w2c[c] = acc_w*ln10v

# ================= PCHIP in numba =================
@njit(cache=True)
def _edge_case(h0, h1, m0, m1):
    d = ((2.0*h0 + h1)*m0 - h0*m1)/(h0 + h1)
    def sgn(v):
        if v > 0.0: return 1.0
        if v < 0.0: return -1.0
        return 0.0
    if sgn(d) != sgn(m0):
        return 0.0
    if (sgn(m0) != sgn(m1)) and (abs(d) > 3.0*abs(m0)):
        return 3.0*m0
    return d

@njit(cache=True)
def pchip_fine(idx_out, y, nv, out):
    nc = len(idx_out)
    hk = np.empty(nc-1); mk = np.empty(nc-1)
    for i in range(nc-1):
        hk[i] = idx_out[i+1] - idx_out[i]
        mk[i] = (y[i+1] - y[i])/hk[i]
    dk = np.empty(nc)
    for i in range(1, nc-1):
        s0 = 1.0 if mk[i-1] > 0.0 else (-1.0 if mk[i-1] < 0.0 else 0.0)
        s1 = 1.0 if mk[i] > 0.0 else (-1.0 if mk[i] < 0.0 else 0.0)
        if (s1 != s0) or (mk[i] == 0.0) or (mk[i-1] == 0.0):
            dk[i] = 0.0
        else:
            w1 = 2.0*hk[i] + hk[i-1]
            w2 = hk[i] + 2.0*hk[i-1]
            whmean = (w1/mk[i-1] + w2/mk[i])/(w1 + w2)
            dk[i] = 1.0/whmean
    dk[0] = _edge_case(hk[0], hk[1], mk[0], mk[1])
    dk[nc-1] = _edge_case(hk[nc-2], hk[nc-3], mk[nc-2], mk[nc-3])
    c = 0
    for p in range(nv):
        while (c+1 < nc-1) and (p >= idx_out[c+1]):
            c += 1
        x0 = idx_out[c]; x1 = idx_out[c+1]
        hh = x1 - x0
        dx = p - x0
        y0 = y[c]; y1 = y[c+1]; d0 = dk[c]; d1 = dk[c+1]
        c0 = 2.0*(y0 - y1)/(hh*hh*hh) + (d0 + d1)/(hh*hh)
        c1 = 3.0*(y1 - y0)/(hh*hh) - (2.0*d0 + d1)/hh
        out[p] = ((c0*dx + c1)*dx + d0)*dx + y0

# ================= full fast SGWB_iter =================
def SGWB_iter_fast(m, tol=1e-4, freq_res=1.0):
    """Accelerated (approximate) self-consistent SGWB iteration.

    Solves the same physics as ``LCDM_SG.SGWB_iter()`` with a fixed-step
    approximate solver.  On success fills the fast-path output attributes and
    returns the model; on failure restores ``m.cosmo_param['DN_eff']`` and
    ``m.DN_eff_orig`` and returns ``None``.  Stale per-channel full-evolution
    attributes from a previous original ``SGWB_iter()`` run are removed so they
    cannot be confused with fast-path results.  ``tol`` is the outer Delta
    N_eff self-consistency stopping criterion (default 1e-4) and ``freq_res``
    scales the frequency-grid density (audit-only; 1.0 = default grid).
    """
    if m.cosmo_param['r'] <= 0:
        print('Must set a positive r to calculate the inflationary GWs!')
        return None
    if m.derived_param['N_inf'] is None:
        print('High-end cutoff frequency has not been set properly.')
        return None
    if getattr(m, 'SGWB_converge', False):
        return m

    # The fast path does not produce the per-channel full-evolution outputs of
    # the original run_SGWB(); drop any leftovers from an earlier LSODA run so
    # stale data cannot be mixed with a fast-path result.
    for stale in ('N_hc', 'Th', 'Oj', 'Ogw', 'Opgw'):
        if hasattr(m, stale):
            try:
                delattr(m, stale)
            except AttributeError:
                pass

    col_step = _COL_STEP        # snapshot; a mid-run set_col_step() must not
    Omega_nu = gp.Omega_nh2/m.derived_param['h']**2    # change array layouts
    DN_eff_orig = m.cosmo_param['DN_eff']
    DN_gw_list = [0.0]; DN_gw_new = 0.0; DN_gw_min = 0.0; DN_gw_max = 10.0
    converged = False
    h = _FAST_H
    z_tail = _Z_TAIL
    if not 0.0 < tol < 1.0:
        raise ValueError('outer tol must be in (0, 1), got %r' % tol)
    freqs = None; Nf = 0; nv = 0; Nv = None
    idx_out = None; n_coarse = 0
    ev_minus = P_t = fp_freq = Wmat = W_last = None
    Ogw = Oj = Opgw = None
    first = True
    try:
        for _iter in range(MAX_ITER):
            gen_fast(m, h)
            m.construct_f(freq_res)
            freqs_new = m.f.astype(np.float64)
            nv_new = len(m.Nv)
            # Reuse grid-dependent quantities across bisection iterations when
            # the frequency/e-fold grids are unchanged (they drift only at the
            # 1e-13 level); this is bit-safe at the 1e-9 tolerance and saves ~10%.
            grid_same = (not first and nv_new == nv and len(freqs_new) == Nf
                         and np.max(np.abs(freqs_new - freqs)) < 1e-9
                         and np.max(np.abs(m.Nv - Nv)) < 1e-9)
            if not grid_same:
                Nf = len(freqs_new); nv = nv_new
                freqs = freqs_new
                Nv = m.Nv.astype(np.float64)
                ev_minus = np.exp(-Nv)
                idx_out = np.unique(np.append(np.arange(0, nv, col_step), nv-1))
                n_coarse = len(idx_out)
                P_t = m.derived_param['A_t']*np.power((10**freqs)/gp.f_piv, m.derived_param['nt'])
                fp_freq = np.power(10.0, freqs)
                Xf = np.flip(freqs); hf = np.diff(Xf)
                Wmat = np.zeros((Nf, Nf))
                build_Wmat(Nf, Xf, hf, Wmat)
                Wmat = np.ascontiguousarray(Wmat)
                W_last = Wmat[Nf-1].copy()
                Ogw = Oj = Opgw = None
            first = False
            Sv, f_hor, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus = prep_fast(m, Nv, freqs, h)
            if Ogw is None or Ogw.shape[0] != Nf or Ogw.shape[1] != n_coarse:
                # Zero-fill, not np.empty: solve_kernel starts each channel at
                # j0 = horizon-crossing + 3 decades, so the early columns
                # (N < N_j0) are never written.  Those cells are read by
                # int_SGWB_W and physically equal 0 (modes still far outside
                # the horizon); np.empty made the result depend on stale
                # heap contents and could occasionally leak NaN into g2/w2.
                Ogw = np.zeros((Nf, n_coarse)); Oj = np.zeros((Nf, n_coarse)); Opgw = np.zeros((Nf, n_coarse))
            solve_kernel(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus, fp_minus, fp_freq,
                         1, n_coarse, col_step, h, z_tail, Ogw, Oj, Opgw)
            g2_last = np.dot(W_last, (Ogw[:, -1] - Oj[:, -1])[::-1]) * ln10
            DN_gw_new = gp.Neff0 * g2_last / Omega_nu
            if not math.isfinite(DN_gw_new):
                print('SGWB_iter_fast: non-finite DN_gw_new, aborting.')
                break
            if DN_eff_orig + DN_gw_new > 5:
                print('SGWB_iter_fast: total N_eff too large, aborting.')
                break
            if abs((gp.Neff0+DN_eff_orig+DN_gw_new)/(gp.Neff0+DN_eff_orig+DN_gw_list[-1]) - 1) < tol:
                converged = True
                break
            if DN_gw_new > DN_gw_list[-1] > DN_gw_min and DN_gw_max >= DN_gw_list[-1]:
                DN_gw_min = DN_gw_list[-1]
            elif DN_gw_new < DN_gw_list[-1] < DN_gw_max and DN_gw_min <= DN_gw_list[-1]:
                DN_gw_max = DN_gw_list[-1]
            if 0 < DN_gw_min <= DN_gw_max < 10:
                DN_gw_new = (DN_gw_min + DN_gw_max)/2
            m.cosmo_param['DN_eff'] = DN_eff_orig + DN_gw_new
            DN_gw_list.append(DN_gw_new)
        else:
            print('SGWB_iter_fast: did not converge within %d iterations.' % MAX_ITER)
    finally:
        if not converged:
            m.cosmo_param['DN_eff'] = DN_eff_orig
            m.DN_eff_orig = None
            m.SGWB_converge = False
    if not converged:
        return None
    m.cosmo_param['DN_eff'] = DN_eff_orig + DN_gw_new
    m.DN_eff_orig = DN_eff_orig
    m.SGWB_converge = True
    m.hubble = math.log10(2*math.pi) + f_hor + (Nv[-1]-Nv)/gp.ln10
    j_hi = np.searchsorted(j0s, idx_out, side='right') - 1
    g2c = np.zeros(n_coarse); w2c = np.zeros(n_coarse)
    int_SGWB_W(Nf, n_coarse, j_hi.astype(np.int64), Wmat, Ogw, Oj, Opgw, g2c, w2c, ln10)
    g2_fine = np.empty(nv); w2_fine = np.empty(nv)
    pchip_fine(idx_out.astype(np.float64), g2c, nv, g2_fine)
    pchip_fine(idx_out.astype(np.float64), w2c, nv, w2_fine)
    m.g2 = g2_fine; m.w2 = w2_fine
    m.DN_gw = gp.Neff0 * np.multiply(g2_fine, np.exp(2*(f_hor-f_hor[-1])*gp.ln10 + 2*(Nv-Nv[-1]))) / Omega_nu
    m.Ogw_today = Ogw[:, -1].copy(); m.Opgw_today = Opgw[:, -1].copy(); m.Oj_today = Oj[:, -1].copy()
    m.log10OmegaGW = np.log10(m.Ogw_today - m.Oj_today)
    m.kappa_r = m.cosmo_param['DN_eff'] * 7/8*(4/11)**(4/3) * gp.z_ratio**4
    return m
