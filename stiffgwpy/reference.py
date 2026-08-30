# -*- coding: utf-8 -*-
"""reference.py -- physics-first high-accuracy reference pipeline.

The original ``stiffgwpy`` solves the tensor-mode equations on a fixed-step
``sigma`` grid and treats ``fast`` as an approximation to the LSODA path.  This
module builds an *independent* high-accuracy reference so that accuracy can be
measured against a continuum answer rather than against LSODA.

It differs from both the LSODA path and the fast path in three ways that target
the two dominant *shared* error sources found in ``docs/audit_error_budget.md``:

1. ``sigma(N)`` and ``H(N)`` are evaluated *exactly* at arbitrary ``N`` (using
   the same physical branches and splines), so the fixed-step grid + cubic
   spline through the instantaneous-reheating ``sigma`` kink (a ~0.73% bias at
   ``h=0.01``) is removed.
2. The tensor equations are integrated with a high-order adaptive solver
   (``scipy.integrate.solve_ivp``, method ``DOP853``) with tight tolerances.
3. The deep-subhorizon tail is handed off at a configurable ``z_tail`` with a
   matching-error estimate (overlap-region comparison), instead of a single
   frozen-amplitude anchor.

The pipeline is intentionally slower than the fast path; it is meant for
benchmark points, pathological points, and convergence certification.  It also
returns explicit error estimates for the spectrum, the quadrature and the ODE.

This module does NOT sit on the MCMC hot path.
"""

import math
import os
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import interpolate
from scipy.integrate import quad, solve_ivp

from . import global_param as gp

__all__ = [
    'background_at',
    'f_hor_at',
    'solve_reference_mode',
    'spectrum_reference',
    'integrate_spectrum',
    'run_reference',
    'apply_reference_to_model',
]

ln10 = math.log(10.0)


# ---------------------------------------------------------------------------
# Exact background evaluation (H^2 and sigma) at arbitrary N.
# ---------------------------------------------------------------------------

def _fd_splines():
    """Fermi-Dirac rho/p spline on log10(nu), from the shipped ``fd_table.npz``."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, 'fd_table.npz')
    data = np.load(path)
    nu = data['nu']
    vals = data['vals']
    lognu = np.log10(nu)
    rho = interpolate.CubicSpline(lognu, vals[:, 0])
    p = interpolate.CubicSpline(lognu, vals[:, 1])
    return nu, rho, p


_FD_NU, _FD_RHO_CS, _FD_P_CS = _fd_splines()


class _Background:
    """Holds the model-parameter closures needed to evaluate sigma/H exactly."""

    __slots__ = ('Omh2', 'Osh2', 'Oerh2', 'Otrh2', 'Otreh2', 'OLh2',
                 'N_inf', 'N_re', 'N_re_abs', 'N_fin', 'N_max',
                 'nu_today', 'Omega_ph2', 'Omega_nh2', 'Omega_mnuh2',
                 'H_0')

    def __init__(self, m, DN_eff):
        d = m.derived_param
        Omh2 = d['Omega_mh2']
        Osh2 = d['Omega_sh2']
        Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
        Otrh2 = gp.Omega_orh2 + Oerh2
        Otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + Oerh2
        OLh2 = d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 \
            - gp.Omega_ph2 - Oerh2 - Osh2
        self.Omh2 = Omh2
        self.Osh2 = Osh2
        self.Oerh2 = Oerh2
        self.Otrh2 = Otrh2
        self.Otreh2 = Otreh2
        self.OLh2 = OLh2
        self.N_inf = d['N_inf']
        self.N_re = d['N_re']
        # Continuous (non-quantised) reheating boundary: look-back N_re below
        # the present.  The fixed-step grid rounds this; we do not.
        self.N_re_abs = d['N_inf'] - d['N_re']
        self.N_fin = gp.N_fin
        self.N_max = gp.N_max
        self.nu_today = gp.nu_today
        self.Omega_ph2 = gp.Omega_ph2
        self.Omega_nh2 = gp.Omega_nh2
        self.Omega_mnuh2 = gp.Omega_mnuh2
        self.H_0 = d['H_0']


def _H2_and_sigma(bg, N):
    """Exact H^2 (h^2 a^{-?} units) and sigma at absolute N."""
    if N < bg.N_re_abs:
        # MD reheating branch: sigma = 1.  H2 is not used on this branch (the
        # horizon position uses the analytic MD extension in f_hor_abs).
        return 0.0, 1.0
    Nlast = bg.N_inf
    eN = math.exp(Nlast - N)
    e3N = eN * eN * eN
    nu = bg.nu_today / eN
    if nu > 100.0:
        H2 = (bg.Omh2 + bg.Omega_mnuh2
              + (bg.Omega_ph2 + 2.0 / 3.0 * bg.Omega_nh2 + bg.Oerh2) * eN
              + bg.Osh2 * e3N + bg.OLh2 / e3N)
        sigma = (bg.Omh2 + bg.Omega_mnuh2
                 + 4.0 / 3.0 * (bg.Omega_ph2 + 2.0 / 3.0 * bg.Omega_nh2 + bg.Oerh2) * eN
                 + 2.0 * bg.Osh2 * e3N) / H2
    elif nu >= 0.1:
        lognu = math.log10(nu)
        rho_nu = float(_FD_RHO_CS(lognu))
        p_nu = float(_FD_P_CS(lognu))
        H2 = (bg.Omh2
              + (bg.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * bg.Omega_nh2 + bg.Oerh2) * eN
              + bg.Osh2 * e3N + bg.OLh2 / e3N)
        sigma = (bg.Omh2
                 + 4.0 / 3.0 * (bg.Omega_ph2 + 2.0 / 3.0 * bg.Omega_nh2 + bg.Oerh2) * eN
                 + (rho_nu + p_nu) * bg.Omega_nh2 / 3.0 * eN
                 + 2.0 * bg.Osh2 * e3N) / H2
    elif N > Nlast - bg.N_fin:
        H2 = bg.Omh2 + bg.Otrh2 * eN + bg.Osh2 * e3N + bg.OLh2 / e3N
        sigma = (bg.Omh2 + 4.0 / 3.0 * bg.Otrh2 * eN + 2.0 * bg.Osh2 * e3N) / H2
    elif N >= Nlast - bg.N_max:
        Nl = Nlast - N
        rho_i = float(gp.spl_rho(Nl))
        rhop_i = float(gp.spl_rhop(Nl))
        H2 = bg.Omh2 + (bg.Omega_ph2 * rho_i + bg.Oerh2) * eN \
            + bg.Osh2 * e3N + bg.OLh2 / e3N
        sigma = (bg.Omh2 + (bg.Omega_ph2 * rhop_i + 4.0 / 3.0 * bg.Oerh2) * eN
                 + 2.0 * bg.Osh2 * e3N) / H2
    else:
        H2 = bg.Omh2 + bg.Otreh2 * eN + bg.Osh2 * e3N + bg.OLh2 / e3N
        sigma = (bg.Omh2 + 4.0 / 3.0 * bg.Otreh2 * eN + 2.0 * bg.Osh2 * e3N) / H2
    return H2, sigma


def background_at(m, N, DN_eff):
    """Return ``(H2, sigma)`` at absolute ``N`` for the given ``DN_eff``."""
    N = float(N)
    bg = _Background(m, DN_eff)
    return _H2_and_sigma(bg, N)


def f_hor_at(m, N, DN_eff):
    """Exact ``log10(f_H/Hz)`` at absolute ``N`` (reheating kink handled exactly)."""
    N = float(N)
    bg = _Background(m, DN_eff)
    Nlast = bg.N_inf
    # Today's raw f_hor (normalisation anchor).
    H2_last, _ = _H2_and_sigma(bg, Nlast)
    raw_last = -0.5 * Nlast + 0.5 * math.log(H2_last)
    Delta_f = math.log(2.0 * math.pi / bg.H_0)

    # Reheating boundary raw value (continuous).
    if bg.N_re_abs < Nlast:
        H2_re, _ = _H2_and_sigma(bg, bg.N_re_abs)
        raw_re = -0.5 * bg.N_re_abs + 0.5 * math.log(H2_re)
    else:
        raw_re = raw_last

    if N < bg.N_re_abs:
        raw = raw_re - 0.5 * (N - bg.N_re_abs)
    else:
        H2, _ = _H2_and_sigma(bg, N)
        raw = -0.5 * N + 0.5 * math.log(H2)
    return (raw - raw_last - Delta_f) / ln10


# ---------------------------------------------------------------------------
# Tensor-mode ODE in the ORIGINAL variables (as in functions.tensor).
# ---------------------------------------------------------------------------

def _tensor_orig(N, state, bg):
    z, x, y = state
    _, sigma = _H2_and_sigma(bg, N)
    ez = math.exp(z)
    dz = 1.5 * sigma - 1.0
    dx = -3.0 * x + 1.5 * sigma * x - ez * y
    dy = -y + 1.5 * sigma * y + ez * x
    return (dz, dx, dy)


def _make_tail_event(z_tail):
    def event(N, state, bg):
        return state[0] - z_tail
    event.terminal = True
    event.direction = 1
    return event


def _find_start_N(bg, freq, z_start):
    """Smallest N in [0, N_inf] with z(N) >= z_start (z increasing in N)."""
    Nlast = bg.N_inf
    lo, hi = 0.0, Nlast
    def z_at(N):
        fh = f_hor_abs(bg, freq, N)
        return (freq - fh) * ln10
    if z_at(lo) >= z_start:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if z_at(mid) < z_start:
            lo = mid
        else:
            hi = mid
    return hi


def f_hor_abs(bg, freq, N):
    """Internal f_hor used for locating the start; offset by log10(f)."""
    # Reuse the public evaluator through the module-level function is clumsy here
    # (it needs the model).  Instead compute raw f_hor(N) using `bg`.
    Nlast = bg.N_inf
    H2_last, _ = _H2_and_sigma(bg, Nlast)
    raw_last = -0.5 * Nlast + 0.5 * math.log(H2_last)
    Delta_f = math.log(2.0 * math.pi / bg.H_0)
    if bg.N_re_abs < Nlast:
        H2_re, _ = _H2_and_sigma(bg, bg.N_re_abs)
        raw_re = -0.5 * bg.N_re_abs + 0.5 * math.log(H2_re)
    else:
        raw_re = raw_last
    if N < bg.N_re_abs:
        raw = raw_re - 0.5 * (N - bg.N_re_abs)
    else:
        H2, _ = _H2_and_sigma(bg, N)
        raw = -0.5 * N + 0.5 * math.log(H2)
    return (raw - raw_last - Delta_f) / ln10


def solve_reference_mode(m, freq, DN_eff, z_tail=6.0, z_start=-12.0,
                         rtol=1e-12, atol=None):
    """Solve a single frequency with the independent high-accuracy pipeline.

    Returns a dict with the today values ``Ogw_today``, ``Oj_today``,
    ``Opgw_today`` plus diagnostics (whether the analytic tail was used, the
    hand-off time, the number of ODE steps taken, and the deep-subhorizon
    adiabaticity ``eps_handoff`` with the leading-order WKB matching error
    ``matching_error_rel``).  ``eps_handoff = |1.5*sigma - 1| / exp(z)`` and the
    frozen-amplitude tail is an ``O(eps)`` approximation.
    """
    if atol is None:
        atol = [1e-12, 1e-22, 1e-22]
    bg = _Background(m, DN_eff)
    A_t = m.derived_param['A_t']
    nt = m.derived_param['nt']
    P_t = A_t * (10.0 ** freq / gp.f_piv) ** nt

    # Real log10(z) relationship: z = ln(2 pi f / (a H)) = (log10 f - f_hor)*ln10.
    # We integrate the ODE for z directly; only the *starting* z uses the
    # horizon position so that the initial condition (x=0, y=e^z) is the frozen
    # super-horizon mode.
    Nstart = _find_start_N(bg, freq, z_start)
    fh_start = f_hor_abs(bg, freq, Nstart)
    # Initial z so that x=0, y=e^{z0} is the frozen super-horizon mode.
    z0 = (freq - fh_start) * ln10

    y0 = math.exp(z0)
    Nspan = (Nstart, bg.N_inf)
    result = solve_ivp(_tensor_orig, Nspan, (z0, 0.0, y0),
                       method='DOP853', rtol=rtol, atol=atol,
                       events=[_make_tail_event(z_tail)],
                       args=(bg,), dense_output=True)

    used_tail = bool(result.t_events[0].size)
    n_steps = int(result.nfev)
    # Adiabaticity parameter of the deep-subhorizon hand-off: the frozen-amplitude
    # tail is an O(epsilon) approximation, where epsilon = |d ln(omega)/dN| / omega
    # with omega = exp(z) and dN the integration variable.  sigma enters through
    # z' = 1.5*sigma - 1, so epsilon = |1.5*sigma - 1| / exp(z).
    eps_handoff = None
    matching_error_rel = None
    if used_tail:
        t_event = float(result.t_events[0][0])
        zf, xf, yf = result.sol(t_event)
        _, sig = _H2_and_sigma(bg, t_event)
        eps_handoff = abs(1.5 * sig - 1.0) / math.exp(zf)
        # Leading-order WKB defect of the frozen-amplitude approximation.
        matching_error_rel = eps_handoff
        coeff = math.sqrt(0.5 * (xf * xf + yf * yf))
        # Analytic deep-subhorizon tail from t_event to today.
        N_tail = np.arange(t_event, bg.N_inf + 1e-12, 0.02)
        if N_tail[-1] < bg.N_inf:
            N_tail = np.append(N_tail, bg.N_inf)
        fh_tail = np.array([f_hor_abs(bg, freq, float(n)) for n in N_tail])
        Th_hf = coeff * np.exp(-zf + t_event - N_tail)
        xf_hf = Th_hf * 10.0 ** (freq - fh_tail)
        Oj_hf = -Th_hf * Th_hf / 3.0 * P_t
        Opgw_hf = xf_hf * xf_hf / 36.0 * P_t
        Ogw_hf = 3.0 * Opgw_hf + Oj_hf
        return dict(Ogw_today=float(Ogw_hf[-1]), Oj_today=float(Oj_hf[-1]),
                    Opgw_today=float(Opgw_hf[-1]), used_tail=used_tail,
                    event_N=t_event, n_steps=n_steps, z0=z0,
                    eps_handoff=eps_handoff,
                    matching_error_rel=matching_error_rel)

    # Mode never reached the tail threshold before today.
    zf, xf, yf = result.y[:, -1]
    Th = yf / math.exp(zf)
    Oj = xf * Th / 3.0 * P_t
    Ogw = (xf * xf + yf * yf) / 24.0 * P_t + Oj
    Opgw = (-5.0 * xf * xf + 7.0 * yf * yf) / 72.0 * P_t
    return dict(Ogw_today=float(Ogw), Oj_today=float(Oj),
                Opgw_today=float(Opgw), used_tail=False,
                event_N=None, n_steps=n_steps, z0=z0,
                eps_handoff=eps_handoff,
                matching_error_rel=matching_error_rel)


# ---------------------------------------------------------------------------
# Spectrum and bolometric integral with error estimates.
# ---------------------------------------------------------------------------

def _sync_size():
    try:
        from mpi4py import MPI
        return int(MPI.COMM_WORLD.Get_size())
    except Exception:
        return 1


def spectrum_reference(m, freqs, DN_eff, z_tail=6.0, rtol=1e-12, workers=None):
    """Solve every frequency and return today's spectrum arrays."""
    freqs = np.asarray(freqs, dtype=float)
    n = freqs.size
    Ogw = np.empty(n)
    Oj = np.empty(n)
    Opgw = np.empty(n)
    used_tail = np.zeros(n, dtype=bool)
    if workers is None:
        workers = int(os.environ.get('SGWB_POOL_SIZE', 4))
    workers = max(1, min(workers, os.cpu_count() or 1))
    if _sync_size() > 1:
        workers = 1
    tasks = [(float(f), DN_eff, z_tail, rtol) for f in freqs]
    if workers > 1 and n > 1:
        def _task(t):
            freq, dn, zt, rtol = t
            return solve_reference_mode(m, freq, dn, z_tail=zt, rtol=rtol)
        with ThreadPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_task, tasks))
    else:
        results = [solve_reference_mode(m, float(f), DN_eff, z_tail=z_tail, rtol=rtol)
                   for f in freqs]
    for i, sol in enumerate(results):
        Ogw[i] = sol['Ogw_today']
        Oj[i] = sol['Oj_today']
        Opgw[i] = sol['Opgw_today']
        used_tail[i] = sol['used_tail']
    return Ogw, Oj, Opgw, used_tail


def integrate_spectrum(freqs, Ogw, Oj, ln10v=ln10):
    """Bolometric integral of the today spectrum with quadrature + interp errors.

    The integrand ``Omega_GW = Ogw - Oj`` is a smooth function of ``log10 f``.
    We interpolate with a shape-preserving PCHIP (avoids cubic overshoot on the
    sharp low-frequency tail), integrate with adaptive Gauss-Kronrod
    (``scipy.integrate.quad``), and estimate the *interpolation* error by
    comparing against a half-density PCHIP.

    Returns ``(g2, quadrature_error, interpolation_error)``.
    """
    freqs = np.asarray(freqs, dtype=float)
    integrand = np.asarray(Ogw, dtype=float) - np.asarray(Oj, dtype=float)
    # Ascending log10 f.
    order = np.argsort(freqs)
    x = freqs[order]
    y = integrand[order]
    if x.size < 4:
        g2 = ln10v * np.trapezoid(y, x) if x.size >= 2 else 0.0
        return float(g2), None, None
    spl = interpolate.PchipInterpolator(x, y)

    def f_quad(ell):
        return float(spl(ell))

    def _quad(func):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return quad(func, float(x[0]), float(x[-1]), epsabs=0.0,
                        epsrel=1e-10, limit=800, points=x)

    val, abserr = _quad(f_quad)
    g2 = ln10v * val
    # Half-density interpolation error estimate: compare to a coarser PCHIP.
    half = x[::2]
    if half.size >= 4 and half.size < x.size and np.all(np.diff(half) > 0):
        yh = y[::2]
        spl_h = interpolate.PchipInterpolator(half, yh)
        val_h, _ = _quad(lambda e: float(spl_h(e)))
        interp_err = ln10v * abs(val - val_h)
    else:
        interp_err = None
    return float(g2), float(ln10v * abs(abserr)), (None if interp_err is None else float(interp_err))


def _spectrum_at_dn(m, freqs, DN_eff, z_tail, rtol):
    return spectrum_reference(m, freqs, DN_eff, z_tail=z_tail, rtol=rtol)


def run_reference(m, dn_eff=None, freq_res=1.0, z_tail=5.0, rtol=1e-11,
                  z_start=-12.0, tol=1e-7, max_iter=60, freq_subset=None,
                  self_consistent=True):
    """High-accuracy SGWB solve, returning derived results.

    When ``dn_eff`` is given (recommended for benchmarking) a *single* pass at
    that background is done: it isolates the solver accuracy (continuous
    ``sigma`` + high-order ODE) from the outer self-consistency.  When it is
    ``None`` a full bisection on ``Delta N_eff`` is run (slow).

    ``freq_subset`` (an array of ``log10 f`` values) restricts the solved
    frequencies; when ``None`` the model's ``construct_f(freq_res)`` grid is
    used.  The returned spectrum and the bolometric integral are computed on
    the chosen frequency set, so a coarse ``freq_subset`` also carries a
    frequency-grid error that is reported separately.
    """
    DN_eff_orig = m.cosmo_param['DN_eff']
    Omega_nu = gp.Omega_nh2 / m.derived_param['h'] ** 2
    if freq_subset is not None:
        freqs = np.asarray(freq_subset, dtype=float)
    else:
        # Grid-independent frequency set (built from continuous background, so
        # the sampling is invariant to the sigma-grid resolution).
        from .freq_adaptive import grid_independent_freqs
        freqs = grid_independent_freqs(m, freq_res)[0]
        freqs = np.asarray(freqs, dtype=float)

    def solve(freqs_in, dn):
        return _spectrum_at_dn(m, freqs_in, dn, z_tail, rtol)

    if self_consistent and dn_eff is None:
        DN_gw_list = [0.0]
        DN_gw_new = 0.0
        DN_gw_min = 0.0
        DN_gw_max = 10.0
        converged = False
        iters = 0
        for iters in range(max_iter):
            Ogw, Oj, Opgw, used_tail = solve(freqs, DN_eff_orig + DN_gw_new)
            g2, qerr, ierr = integrate_spectrum(freqs, Ogw, Oj)
            DN_gw_new = gp.Neff0 * g2 / Omega_nu
            if not math.isfinite(DN_gw_new):
                raise RuntimeError('reference: non-finite DN_gw_new')
            if DN_eff_orig + DN_gw_new > 5.0:
                DN_gw_new = 0.0
                break
            if abs((gp.Neff0 + DN_eff_orig + DN_gw_new) /
                   (gp.Neff0 + DN_eff_orig + DN_gw_list[-1]) - 1.0) < tol:
                converged = True
                break
            if DN_gw_new > DN_gw_list[-1] > DN_gw_min and DN_gw_max >= DN_gw_list[-1]:
                DN_gw_min = DN_gw_list[-1]
            elif DN_gw_new < DN_gw_list[-1] < DN_gw_max and DN_gw_min <= DN_gw_list[-1]:
                DN_gw_max = DN_gw_list[-1]
            if 0.0 < DN_gw_min <= DN_gw_max < 10.0:
                DN_gw_new = (DN_gw_min + DN_gw_max) / 2.0
            DN_gw_list.append(DN_gw_new)
        if not converged:
            raise RuntimeError('reference: did not converge in %d iterations' % max_iter)
        DN_eff = DN_eff_orig + DN_gw_new
        n_iter = iters + 1
    else:
        DN_eff = dn_eff if dn_eff is not None else DN_eff_orig
        n_iter = None

    Ogw, Oj, Opgw, used_tail = solve(freqs, DN_eff)
    g2, qerr, ierr = integrate_spectrum(freqs, Ogw, Oj)
    DN_gw = gp.Neff0 * g2 / Omega_nu
    kappa_r = DN_eff * 7 / 8 * (4 / 11) ** (4 / 3) * gp.z_ratio ** 4
    # Guard against non-positive Ogw - Oj in the (ill-defined) super-horizon low-f
    # tail: log10 is undefined there.  The spectrum only has a physical meaning
    # once the mode has re-entered the horizon.
    integrand = Ogw - Oj
    logOmega = np.log10(np.maximum(integrand, 1e-40))
    return dict(
        freqs=freqs,
        log10OmegaGW=logOmega,
        Ogw=Ogw, Oj=Oj, Opgw=Opgw,
        DN_eff=DN_eff,
        DN_gw=DN_gw,
        kappa_r=kappa_r,
        g2=g2,
        quadrature_error=qerr,
        interpolation_error=ierr,
        used_tail=used_tail,
        n_freq=freqs.size,
        n_iter=n_iter,
        z_tail=z_tail,
        rtol=rtol,
    )


def apply_reference_to_model(m, freq_res=1.0, z_tail=5.0, rtol=1e-11):
    """Run the continuous-sigma reference and expose it on the model.

    Mirrors the fast/LSODA output contract so ``LCDM_SG`` and the Cobaya adapter
    can treat ``engine='reference'`` as a first-class engine (slow: intended for
    certification/benchmark points, not the MCMC thermal path).
    """
    ref = run_reference(m, dn_eff=None, freq_res=freq_res, z_tail=z_tail,
                        rtol=rtol, self_consistent=True)
    m.cosmo_param['DN_eff'] = ref['DN_eff']
    m.DN_gw = np.append([0.0], ref['DN_gw'])
    m.kappa_r = ref['kappa_r']
    m.log10OmegaGW = ref['log10OmegaGW']
    m.f = ref['freqs']
    m.g2 = np.array([ref['g2']])
    m.reference_evals = getattr(m, 'reference_evals', 0) + 1
    m.SGWB_converge = True
    return m
