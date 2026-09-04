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
    from stiffgwpy_fast import LCDM_SG
    from stiffgwpy_fast import fast_sgwb

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
from contextlib import nullcontext
from threading import RLock

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads
from scipy import interpolate

from . import global_param as gp
from ._resources import package_path
from .config import FastSolverConfig
from .functions import int_FD

__all__ = ['SGWB_iter_fast', 'gen_fast', 'set_threads', 'set_col_step', 'set_h',
           'set_z_tail', 'get_settings', 'apply_accuracy_mode', 'ACCURACY_MODES',
           'USER_FAST_PROFILES', 'FAST_PROFILES', 'normalize_accuracy_mode',
           'is_validation_mode', 'MODE_ROLE', 'FastSolverConfig', 'get_config',
           'resolve_config', 'max_threads']

# Default OpenMP threads: numba's own default (no more than the detected core
# count).  We do NOT force a fixed number at import time -- that previously
# raised ValueError on machines with fewer than 32 cores -- and only call
# set_num_threads() when FAST_THREADS is explicitly set.
_MAX_THREADS = get_num_threads()
_THREADS = _MAX_THREADS
_NUMBA_CONFIG_LOCK = RLock()
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

# Maximum phase increment per (sub-)step, dTheta = e^z * dh <= _PHASE_MAX,
# used by the horizon-crossing adaptive step control in solve_kernel.  0.0
# disables sub-stepping (pure fixed-step Magnus on the grid).
_PHASE_MAX = 0.0
_fast_phase_env = _os.environ.get('FAST_PHASE_MAX')
if _fast_phase_env is not None:
    _PHASE_MAX = float(_fast_phase_env)
    if not 0.0 <= _PHASE_MAX <= 10.0:
        raise ValueError('FAST_PHASE_MAX must be in [0, 10], got %r'
                         % _fast_phase_env)

# Default frequency-grid builder for SGWB_iter_fast (overridable per call via
# the freq_grid argument; apply_accuracy_mode sets it from the preset).
_FREQ_GRID = 'construct'

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


def max_threads():
    """返回当前 Numba 进程可用的最大线程数。"""
    return _MAX_THREADS


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


def set_phase_max(pm):
    """Set the max phase increment per (sub-)step (0 disables sub-stepping)."""
    global _PHASE_MAX
    pm = float(pm)
    if not 0.0 <= pm <= 10.0:
        raise ValueError('phase_max must be in [0, 10], got %r' % pm)
    _PHASE_MAX = pm


def set_freq_grid(name):
    """Set the default frequency-grid builder ('construct'/'grid_independent'/'adaptive')."""
    global _FREQ_GRID
    if name not in ('construct', 'grid_independent', 'adaptive'):
        raise ValueError('freq_grid must be construct/grid_independent/adaptive, got %r' % name)
    _FREQ_GRID = name


def get_settings():
    """Snapshot legacy module settings, including the selected grid builder."""
    return dict(threads=_THREADS, col_step=_COL_STEP, h=_FAST_H, z_tail=_Z_TAIL,
                phase_max=_PHASE_MAX, freq_grid=_FREQ_GRID)


# Named accuracy presets (audit phase "three recommended modes").
# Values come from the phase-2 convergence study: engine-vs-LSODA difference
# at h=0.01 is ~1e-5 on Delta N_eff while the shared sigma-grid bias vs a
# deep reference is ~0.73% (h=0.01) / ~0.33% (h=0.005); z_tail=7 reduces the
# analytic-tail error to ~2e-5 and z_tail=10 to ~2.4e-7; col_step has <1e-9
# effect on the final Delta N_eff (it only shapes early small-value curves);
# freq_res=2.0 halves the low-frequency-tail undersampling error.
ACCURACY_MODES = {
    'debug': dict(h=0.005, col_step=1, z_tail=10.0, freq_res=2.0,
                  tol=1e-8, threads=8, transition_refine=True, phase_max=0.25),
    'fast': dict(h=0.02, col_step=8, z_tail=5.0, freq_res=1.0,
                 tol=1e-6, threads=16, transition_refine=False, phase_max=0.0),
    # Backward-compatible alias of 'fast' (identical settings).
    'ultra-fast': dict(h=0.02, col_step=8, z_tail=5.0, freq_res=1.0,
                       tol=1e-6, threads=16, transition_refine=False, phase_max=0.0),
    'reference': dict(h=0.00125, col_step=1, z_tail=10.0, freq_res=2.0,
                      tol=1e-8, threads=8, transition_refine=True, phase_max=0.1),
    'production': dict(h=0.01, col_step=4, z_tail=8.0, freq_res=1.0,
                       tol=1e-7, threads=8, transition_refine=True, phase_max=0.5,
                       freq_grid='adaptive'),
    # Deep-tail variant of production: the frozen-tail (WKB handoff) error at
    # z_tail=8 is ~2.5e-4 relative per mode vs the deep limit, so production
    # already keeps the tail term subdominant in the error budget; this mode
    # pushes it to ~3e-5 for the validation chain.
    'deep': dict(h=0.01, col_step=4, z_tail=10.0, freq_res=1.0,
                 tol=1e-7, threads=8, transition_refine=True, phase_max=0.25,
                 freq_grid='adaptive'),
}

# 档位说明：对外只提供两个快速档位（速度优先的 plain-grid 和精度优先的
# transition-refine），并保留向后兼容的别名。其余档位（debug/deep/reference）
# 是同一内部求解器的验证或基准变体，可用于认证，但不作为第三、第四个生产档位，
# 也不用于 MCMC 热路径；真正的精度锚点是连续 sigma 的
# ``stiffgwpy_fast.reference`` 流程（engine='reference'）。
USER_FAST_PROFILES = ('fast', 'production')

# Alias -> canonical mode name (accepts the human-facing names used in docs).
FAST_PROFILE_ALIASES = {
    'plain_grid': 'fast',
    'plain-grid': 'fast',
    'plain': 'fast',
    'transition_refine': 'production',
    'transition-refine': 'production',
    'tr': 'production',
}

# Role per canonical mode: 'fast' = user-facing fast profile, 'validation' =
# certification/benchmark variant (not a production tier).
MODE_ROLE = {
    'fast': 'fast',
    'ultra-fast': 'fast',       # alias of 'fast'
    'production': 'fast',
    'debug': 'validation',
    'deep': 'validation',
    'reference': 'validation',
}

# The two user-facing fast profiles, described as delivered (plain-grid vs
# transition-refine).  ``profile`` is the human-facing name in docs.
FAST_PROFILES = {
    'fast': dict(ACCURACY_MODES['fast'], profile='plain-grid',
                 role='user_fast'),
    'production': dict(ACCURACY_MODES['production'], profile='transition-refine',
                       role='user_fast'),
}


def normalize_accuracy_mode(name):
    """Resolve an accuracy-mode name/alias to a canonical :data:`ACCURACY_MODES` key.

    Accepts the user-facing names ``plain_grid`` / ``plain-grid`` (alias of
    ``fast``) and ``transition_refine`` / ``transition-refine`` (alias of
    ``production``), plus the historical keys.  ``None`` is returned unchanged.
    """
    if name is None:
        return None
    name = str(name)
    name = FAST_PROFILE_ALIASES.get(name, name)
    if name not in ACCURACY_MODES:
        raise ValueError(
            'unknown accuracy mode %r; choose from %s '
            '(user-facing fast profiles: %s)' % (name, sorted(ACCURACY_MODES),
                                                 list(USER_FAST_PROFILES)))
    return name


def is_validation_mode(name):
    """True if ``name`` is a validation/benchmark variant, not a user-facing fast profile."""
    return MODE_ROLE.get(normalize_accuracy_mode(name), 'fast') == 'validation'


# Calibrated error budgets per accuracy mode, from the physics-first benchmark
# (see docs/audit_reference.md).  ``model_bias`` is the dominant *shared*
# continuous-sigma-vs-fixed-grid bias (~1% at h=0.01) that a fast-vs-fast
# convergence check cannot detect; ``ode``/``quadrature`` are engine + grid
# convergence terms; ``tail`` is the analytic-tail error at the mode's z_tail;
# ``spectrum_dex`` is the pointwise log10(Omega_GW) error (dominated by the
# frequency-grid resolution of the spectral features).
ERROR_BUDGET = {
    # model_bias is the continuous-sigma (reference.py) vs fixed-grid sigma bias
    # at the mode's h, taken from the measured h-convergence curve (default point):
    # h=0.02 -> +3.9%, h=0.01 -> +1.3%, h=0.005 -> +0.55%,
    # h=0.0025 -> +0.46%, h=0.00125 -> +0.22%, h~0 -> +0.10%.
    'fast': dict(model_bias=3.9e-2, ode=4.0e-5, quadrature=1.0e-3,
                 tail=3.8e-3, spectrum_dex=0.20),
    'ultra-fast': dict(model_bias=3.9e-2, ode=4.0e-5, quadrature=1.0e-3,
                       tail=3.8e-3, spectrum_dex=0.20),
    'production': dict(model_bias=1.3e-2, ode=1.0e-5, quadrature=1.0e-4,
                       tail=2.0e-5, spectrum_dex=0.07),
    'deep': dict(model_bias=1.3e-2, ode=2.5e-6, quadrature=1.0e-4,
                 tail=2.4e-7, spectrum_dex=0.02),
    'debug': dict(model_bias=5.5e-3, ode=5.0e-6, quadrature=5.0e-5,
                  tail=2.4e-7, spectrum_dex=0.02),
    'reference': dict(model_bias=2.2e-3, ode=2.6e-7, quadrature=1.0e-5,
                      tail=2.4e-7, spectrum_dex=0.002),
}

# Calibrated coefficient for the ODE/horizon-crossing phase-truncation error:
# the constant-z Magnus sub-stepping leaves an O(phase_max^2) relative
# amplitude/phase error per sub-step (validated by a phase_max Richardson
# study, see docs/audit_reference.md); the coefficient is the relative
# Delta N_eff error per phase_max^2 at the default point.
_ODE_PHASE_COEF = 0.0   # filled by calibrate_ode_phase_coef() from measurement






def apply_accuracy_mode(name):
    """Apply a named accuracy preset to the module settings and return its table.

    ``name`` must be a key of :data:`ACCURACY_MODES`.  The preset's
    threads/col_step/h/z_tail are applied through the usual setters (threads is
    clamped to the numba-detected maximum); the preset's freq_res/tol are
    returned in the table so the caller can forward them to
    :func:`SGWB_iter_fast`.  Module settings are process-global, as documented
    for the setters.
    """
    name = normalize_accuracy_mode(name)
    cfg = dict(ACCURACY_MODES[name])
    set_col_step(cfg['col_step'])
    set_h(cfg['h'])
    set_z_tail(cfg['z_tail'])
    set_phase_max(cfg.get('phase_max', 0.0))
    set_freq_grid(cfg.get('freq_grid', 'construct'))
    set_threads(min(cfg['threads'], _MAX_THREADS))
    return cfg


def get_config():
    """Return an immutable snapshot of the legacy module settings."""
    return FastSolverConfig(h=_FAST_H, col_step=_COL_STEP, z_tail=_Z_TAIL,
                            phase_max=_PHASE_MAX, freq_grid=_FREQ_GRID,
                            threads=_THREADS)


def resolve_config(name=None, base=None, **overrides):
    """Resolve a named preset and explicit overrides without mutating globals.

    ``name=None`` snapshots the compatibility settings. A named preset is
    resolved from :data:`ACCURACY_MODES`; only keys with non-``None`` values in
    ``overrides`` replace it. Preset thread counts are clamped to the current
    Numba process budget, matching the historical preset behavior.
    """
    valid_keys = {'h', 'col_step', 'z_tail', 'phase_max', 'freq_grid', 'threads'}
    unknown = sorted(set(overrides) - valid_keys)
    if unknown:
        raise TypeError('unknown FastSolverConfig override(s): %s' % ', '.join(unknown))
    if base is not None:
        if not isinstance(base, FastSolverConfig):
            raise TypeError("base must be a FastSolverConfig")
        values = dict(h=base.h, col_step=base.col_step, z_tail=base.z_tail,
                      phase_max=base.phase_max, freq_grid=base.freq_grid,
                      threads=base.threads)
    elif name is None:
        values = dict(h=_FAST_H, col_step=_COL_STEP, z_tail=_Z_TAIL,
                      phase_max=_PHASE_MAX, freq_grid=_FREQ_GRID,
                      threads=_THREADS)
    else:
        canonical = normalize_accuracy_mode(name)
        preset = ACCURACY_MODES[canonical]
        values = dict(h=preset['h'], col_step=preset['col_step'],
                      z_tail=preset['z_tail'], phase_max=preset.get('phase_max', 0.0),
                      freq_grid=preset.get('freq_grid', 'construct'),
                      threads=min(int(preset['threads']), _MAX_THREADS))
    for key in values:
        if key in overrides and overrides[key] is not None:
            values[key] = overrides[key]
    return FastSolverConfig(**values)

def estimate_error(name='production'):
    """Return the calibrated error budget for the named accuracy mode.

    The returned dict carries the per-stage relative errors used by the
    production error gate (``DN_gw_error`` / ``spectrum_error`` /
    ``quadrature_error`` / ``integration_error``).  The engine terms
    (``ode``, ``quadrature``) are converged to ~1e-5; the physically dominant
    ``model_bias`` is the continuous-sigma vs fixed-grid bias that a fast-only
    convergence study cannot detect, so it is derived from the reference
    benchmark rather than from a fast-vs-fast comparison.
    """
    name = normalize_accuracy_mode(name)
    if name not in ERROR_BUDGET:
        raise ValueError('unknown accuracy mode %r; choose from %s'
                         % (name, sorted(ERROR_BUDGET)))
    b = ERROR_BUDGET[name]
    integration = max(b['model_bias'], b['ode'], b['quadrature'], b['tail'])
    return dict(
        accuracy_mode=name,
        DN_gw_error=integration,
        spectrum_error=b['spectrum_dex'],
        quadrature_error=b['quadrature'],
        integration_error=integration,
        ODE_error=b['ode'],
        tail_error=b['tail'],
        model_bias_error=b['model_bias'],
    )


# Measured anchors for :func:`estimate_local_error` (default point, this
# session; see docs/audit_reference.md and docs/reference/deep_oracle_default.json):
#   * phase_max Richardson: |dDeltaN(pm 0.5 -> 0.125)| = 7.5e-6 relative,
#     so the frozen-z Magnus phase-truncation error at pm=0.5 is ~8e-6 (pm^2 scaling).
#   * kink-breakpoint (transition_refine) residual envelope: ~1e-4.
#   * plain-grid sigma-kink residual: ~2.4e-3; sigma_exact residual at z8: ~4.2e-3.
#   * frozen-tail aggregate: 2e-5 at z_tail=8, 2.4e-7 at z_tail=10 (z-decay e^-(z-8)).
#   * uniform grid_independent quadrature bias: 5.6e-4 at freq_res=1 (242 pt),
#     ~0 at freq_res>=2; construct-grid bias ~1.6e-3 (246 pt).
#   * col_step PCHIP interpolation: <1e-9 on Delta N_eff.
_LOCAL_ODE_PM05 = 8.0e-6          # relative DN error at phase_max=0.5
_LOCAL_ODE_NO_SUBSTEP = 1.0e-5    # h=0.01 fixed-step (calibrated engine term)
_LOCAL_TRANSITION_REFINED = 1.0e-4
_LOCAL_TRANSITION_PLAIN = 2.4e-3
_LOCAL_TRANSITION_SIGMA_EXACT = 4.2e-3
_LOCAL_TAIL_ANCHOR = 2.0e-5       # relative DN error at z_tail=8
_LOCAL_TAIL_FLOOR = 2.4e-7        # z_tail=10 limit
_LOCAL_GRID_UNIFORM_1 = 5.6e-4    # grid_independent freq_res=1 (242 pt)
_LOCAL_GRID_UNIFORM_2 = 2.0e-5    # grid_independent freq_res>=2 (converged)
_LOCAL_GRID_CONSTRUCT = 1.6e-3    # construct grid (246 pt)
_LOCAL_INTERP = 1.0e-9            # col_step PCHIP interpolation on DN_eff
_LOCAL_SELF_CONS_FLOOR = 1.0e-5   # bisection floor when bracket is not recorded


def estimate_local_error(m):
    """A-posteriori, point-local error budget for the last fast solve on ``m``.

    Returns a dict with one entry per physics error category plus the combined
    ``DN_gw_error`` (relative error on the integrated Delta N_eff) and
    ``Delta_Neff_abs_error`` (absolute).  Each category carries ``value`` and
    ``kind``:

    * ``'local'`` -- computed a-posteriori from this solve's telemetry
      (``handoff_eps``, ``freq_grid_error``, ``phase_max_used``, ``z_tail_used``,
      ``quadrature_error_local``, ``dn_bracket``, ``cancellation_ratio``);
    * ``'calibrated'`` -- measured anchors (default-point fast-vs-reference,
      see docs/audit_reference.md) scaled to the solve's settings.

    Categories (physics meaning): ``background_model`` (continuous-sigma vs
    fixed-grid / present-day-anchor bias), ``sigma_transition`` (reheating-kink
    grid-phase residual), ``ode_integration`` (frozen-z Magnus phase truncation),
    ``horizon_crossing`` (the z~0 crossing-band share of that truncation),
    ``wkb_handoff`` (adiabaticity defect at the WKB handoff node, per-mode max
    and weighted aggregate), ``interpolation`` (PCHIP fine-grid, col_step),
    ``frequency_grid`` (sampling of the spectral features), ``quadrature``
    (Simpson-rule integral error on the frequency grid), ``tail_approximation``
    (frozen analytic tail beyond z_tail), ``floating_point`` (Ogw-Oj cancellation
    plus accumulation rounding), ``self_consistency`` (outer bisection bracket).
    The combined error is systematic (max of model/transition) plus the RSS of
    the remaining independent terms.
    """
    dn_gw = np.asarray(getattr(m, 'DN_gw', [0.0]))
    dn = abs(float(dn_gw[-1])) if dn_gw.size else 0.0
    eps64 = np.finfo(float).eps
    he = getattr(m, 'handoff_eps', None)
    if he is None:
        handoff = np.array([0.0])
    else:
        handoff = np.asarray(he, dtype=float)
        handoff = handoff[handoff >= 0.0]
    eps_max = float(np.max(handoff)) if handoff.size else 0.0
    eps_mean = float(np.mean(handoff)) if handoff.size else 0.0
    z_tail = float(getattr(m, 'z_tail_used', 8.0))
    pm = float(getattr(m, 'phase_max_used', 0.0))
    fg = getattr(m, 'freq_grid_used', 'construct')
    tr = bool(getattr(m, 'transition_refine_used', False))
    se = bool(getattr(m, 'sigma_exact_used', False))
    nf = int(getattr(m, 'freq_grid_n', 0))

    # 1. background/model error (present-day anchor + continuous-sigma bias)
    if se:
        model = _LOCAL_TRANSITION_SIGMA_EXACT
    elif tr:
        model = _LOCAL_TRANSITION_REFINED
    else:
        model = _LOCAL_TRANSITION_PLAIN
    # 2. sigma-transition error (reheating-kink grid-phase residual)
    transition = _LOCAL_TRANSITION_REFINED if tr else _LOCAL_TRANSITION_PLAIN
    # 3. ODE integration error: frozen-z Magnus phase truncation (local)
    ode = _LOCAL_ODE_PM05*(pm/0.5)**2 if pm > 0.0 else _LOCAL_ODE_NO_SUBSTEP
    # 4. horizon-crossing error: the crossing band carries the leading
    #    observable (phase-averaged amplitude) effect of the same truncation.
    horizon = ode
    # 5. WKB handoff error: per-mode adiabaticity defect at the handoff node
    #    (local, verifiable per mode).  The aggregate Delta N_eff effect is
    #    reported under tail_approximation (the frozen-tail weight share).
    wkb_max = eps_max
    wkb_agg = eps_mean
    # 6. interpolation error: col_step PCHIP fine-grid on Delta N_eff
    interp = _LOCAL_INTERP
    # 7. frequency-grid error: adaptive -> local criterion; else calibrated
    if fg == 'adaptive':
        freq_grid = max(float(getattr(m, 'freq_grid_error', 0.0)), 0.0)
    elif fg == 'grid_independent':
        fr = float(getattr(m, 'freq_res_used', 1.0) or 1.0)
        freq_grid = _LOCAL_GRID_UNIFORM_2 if fr >= 2.0 else _LOCAL_GRID_UNIFORM_1
    else:  # 'construct'
        freq_grid = _LOCAL_GRID_CONSTRUCT
    # 8. quadrature error: local Richardson estimate on the final grid
    quadrature = max(float(getattr(m, 'quadrature_error_local', 0.0)), 0.0)
    # 9. tail approximation error: calibrated anchor scaled by z_tail
    tail = max(_LOCAL_TAIL_FLOOR, _LOCAL_TAIL_ANCHOR*math.exp(-(z_tail - 8.0)))
    # 10. floating-point/cancellation error (local)
    fp_err = float(getattr(m, 'floating_point_error', 0.0))
    if fp_err <= 0.0:
        fp_err = eps64*math.sqrt(max(nf, 1))
    # 11. self-consistency iteration error: converged successive difference
    #     (local; the outer criterion bounds |DN_gw_new - DN_gw_prev|).
    dn_delta = float(getattr(m, 'dn_converged_delta', 0.0))
    if dn > 0.0 and dn_delta > 0.0:
        self_consistency = max(dn_delta/dn, _LOCAL_SELF_CONS_FLOOR)
    else:
        self_consistency = _LOCAL_SELF_CONS_FLOOR

    # Honest certification tagging: each term is either *measured from this
    # solve's telemetry* ('local-measured') or *anchored to a fiducial-point
    # measurement* ('calibrated-at-fiducial').  A calibrated term is NOT a
    # universal / per-point error estimate; it is an empirical bound measured
    # at the reference point and scaled by the solve's settings.
    def _cert(kind, local):
        return 'local-measured' if local else 'calibrated-at-fiducial'

    cats = dict(
        background_model=dict(value=model, kind='calibrated',
                              certification='calibrated-at-fiducial'),
        sigma_transition=dict(value=transition, kind='calibrated',
                              certification='calibrated-at-fiducial'),
        ode_integration=dict(value=ode, kind='local',
                             certification='calibrated-at-fiducial'),
        horizon_crossing=dict(value=horizon, kind='local',
                              certification='calibrated-at-fiducial'),
        wkb_handoff=dict(value=wkb_max, aggregate=wkb_agg, kind='local',
                         certification='local-measured'),
        interpolation=dict(value=interp, kind='calibrated',
                           certification='calibrated-at-fiducial'),
        frequency_grid=dict(value=freq_grid,
                            kind='local' if fg == 'adaptive' else 'calibrated',
                            certification=_cert('local' if fg == 'adaptive' else 'calibrated',
                                                fg == 'adaptive')),
        quadrature=dict(value=quadrature, kind='local',
                        certification='local-measured'),
        tail_approximation=dict(value=tail, kind='calibrated',
                                certification='calibrated-at-fiducial'),
        floating_point=dict(value=fp_err,
                            cancellation_ratio=float(getattr(m, 'cancellation_ratio', 0.0)),
                            kind='local', certification='local-measured'),
        self_consistency=dict(value=self_consistency, kind='local',
                              certification='local-measured'),
    )
    systematic = max(model, transition)
    rss = math.sqrt(sum(float(cats[k]['value'])**2 for k in cats
                        if k not in ('background_model', 'sigma_transition',
                                     'wkb_handoff')))
    dn_gw_error = systematic + rss
    local_terms = [k for k, c in cats.items() if c['kind'] == 'local']
    calibrated_terms = [k for k, c in cats.items()
                        if c['certification'] == 'calibrated-at-fiducial']
    # Has this solve produced the telemetry needed to certify a per-point
    # budget at all?  Without it the numbers are defaults, not a measurement.
    has_telemetry = (hasattr(m, 'z_tail_used') and hasattr(m, 'handoff_eps')
                     and hasattr(m, 'freq_grid_used'))
    if not has_telemetry:
        certification_status = 'uncertified'
    else:
        # The dominant systematic terms (background_model / sigma_transition)
        # are always fiducial-calibrated, so the budget is point-local only
        # in the random terms (RSS), never fully certified per point.
        certification_status = 'certified-fiducial-calibrated'
    return dict(
        categories=cats,
        DN_gw_error=dn_gw_error,
        Delta_Neff_abs_error=dn_gw_error*dn,
        systematic_error=systematic,
        random_rss=rss,
        handoff_eps_max=wkb_max,
        handoff_eps_mean=wkb_agg,
        z_tail_used=z_tail,
        phase_max_used=pm,
        freq_grid_used=fg,
        # Honest certification metadata (see module docstring):
        certification_status=certification_status,
        local_terms=local_terms,
        calibrated_terms=calibrated_terms,
        explicit_uncertainty='uncertainty is a mix of point-local and '
                             'fiducial-calibrated terms; calibrated terms are '
                             'NOT universal per-point error estimates',
    )


# ================= module-level tables (once) =================
def _build_fd_table():
    """Fermi-Dirac rho/p table, loaded from a cached copy when available.

    A precomputed table is shipped inside the package (``fd_table.npz``); when
    that is missing (e.g. an editable install without data files) we fall back
    to the user cache and finally to recomputation, saving the result to the
    user cache so the next import is cheap.
    """
    try:
        with package_path('stiffgwpy_fast', 'fd_table.npz') as pkg_file:
            with np.load(pkg_file) as data:
                nu, vals = data['nu'], data['vals']
                if nu.shape == (3001,) and vals.shape == (3001, 2):
                    return nu, vals
    except (IOError, OSError, ValueError, KeyError, FileNotFoundError):
        pass
    cache_file = _os.path.join(_os.path.expanduser('~'), '.cache', 'stiffgwpy_fast', 'fd_table.npz')
    if _os.path.exists(cache_file):
        try:
            with np.load(cache_file) as data:
                nu, vals = data['nu'], data['vals']
                if nu.shape == (3001,) and vals.shape == (3001, 2):
                    return nu, vals
        except (IOError, OSError, ValueError, KeyError, FileNotFoundError):
            pass
    nu = np.logspace(-1.0, 2.0, 3001)
    vals = np.array([int_FD(u) for u in nu])
    try:
        cache_dir = _os.path.join(_os.path.expanduser('~'), '.cache', 'stiffgwpy_fast')
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
    # Use the *continuous* present-day anchor N_inf (not the grid-quantised
    # floor(N_inf/h)*h) so the fast path is anchored identically to the
    # continuous-sigma reference.  The floor(N_inf/h)*h anchor attenuates the
    # analytic-tail amplitude by exp(N_inf - floor(N_inf/h)*h) ~ 0.996 -> a ~0.4%
    # systematic bias (grid-anchor error, not a numerical-order issue).
    len_inf = math.floor(d['N_inf']/h)+1
    Nv = np.arange(0, len_inf)*h
    if Nv[-1] != d['N_inf']:
        Nv = np.append(Nv, d['N_inf'])
    index_re = int(np.argmin(np.abs(Nv - (d['N_inf'] - d['N_re']))))
    Sv = np.empty(len(Nv)); f_hor = np.empty(len(Nv))
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
                 assemble, n_coarse, col_step, h, z_tail, Ogw, Oj, Opgw,
                 h_arr=None, Sv=None, phase_max=0.0, handoff_eps=None):
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
                h_step = h_arr[k] if h_arr is not None else h
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if phase_max > 0.0 and z_mid_step > 0.0:
                    n_sub = int(math.ceil(h_step*math.exp(z_mid_step)/phase_max))
                    if n_sub > 1:
                        z_node = z0 + Phi_grid[k] - Phi0
                        dPhi = z_mid_step - z_node
                        h_sub = h_step/n_sub
                        for _s in range(n_sub):
                            zs = z_node + dPhi*(2.0*_s + 1.0)/n_sub
                            xh, yh = scaled_step(xh, yh, zs, h_sub)
                    else:
                        xh, yh = scaled_step(xh, yh, z_mid_step, h_step)
                else:
                    xh, yh = scaled_step(xh, yh, z_mid_step, h_step)
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
            e_z = math.exp(-last_z)
            if Sv is not None:
                gamma = (3.0 - 1.5*Sv[kend])*0.5
                amp2 = (lxh*lxh + lyh*lyh*(1.0 + gamma*gamma*e_z*e_z)
                        + 2.0*gamma*lxh*lyh*e_z)
            else:
                amp2 = lxh*lxh + lyh*lyh
            coeff = math.sqrt(0.5*s2k*amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                if Sv is not None:
                    handoff_eps[m] = abs(1.5*Sv[kend] - 1.0)*e_z
                else:
                    handoff_eps[m] = 0.0
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

def _SGWB_iter_fast_impl(m, tol=1e-4, freq_res=1.0, sigma_exact=False,
                   transition_refine=False, freq_grid=None, config=None,
                   freq_grid_target=3e-4, freq_grid_max_points=1500, eval_freqs=None):
    """Accelerated (approximate) self-consistent SGWB iteration.

    Solves the same physics as ``LCDM_SG.SGWB_iter()`` with a fixed-step
    approximate solver.  On success fills the fast-path output attributes and
    returns the model; on failure restores ``m.cosmo_param['DN_eff']`` and
    ``m.DN_eff_orig`` and returns ``None``.  Stale per-channel full-evolution
    attributes from a previous original ``SGWB_iter()`` run are removed so they
    cannot be confused with fast-path results.  ``tol`` is the outer Delta
    N_eff self-consistency stopping criterion (default 1e-4) and ``freq_res``
    scales the frequency-grid density (audit-only; 1.0 = default grid).  When
    ``sigma_exact`` is set, ``F``/``Phi``/``S2`` are recomputed from the
    continuous piecewise-exact ``sigma`` (reheating kink as an exact breakpoint)
    instead of the fixed-grid cubic spline, removing the ~1% continuous-sigma-vs-
    grid model bias (see ``stiffgwpy_fast.exact_background``).

    ``transition_refine`` (the ``production`` / ``transition-refine`` profile)
    treats the reheating transition as an ODE integration breakpoint: the
    kink-aware grid from ``stiffgwpy_fast.exact_background.build_kink_refined_grid``
    keeps the instantaneous-reheating kink inside a refined sub-step so it is
    never crossed by a spline/grid, and (with ``phase_max > 0``) horizon
    crossing uses phase-aware sub-stepping.  This is the default scientific
    path.  When ``False`` (the ``fast`` / ``plain-grid`` profile) the plain
    fixed-step grid is used — faster but with a larger sigma-kink bias that is
    only certified for the coarse exploratory envelope.

    ``freq_grid`` selects the frequency sampling: ``'construct'`` (the model's
    empirical grid), ``'grid_independent'`` (built from continuous background
    quantities, invariant to the sigma-grid resolution) or ``'adaptive'``
    (seeded from the grid-independent grid and refined where the PCHIP
    curvature ``|y''| h^2 / 8`` of ``log10 Omega_GW`` exceeds
    ``freq_grid_target`` dex, up to ``freq_grid_max_points``; see
    ``freq_adaptive``).  ``None`` uses the module default set by
    :func:`apply_accuracy_mode` / :func:`set_freq_grid`.

    ``eval_freqs`` optionally force-adds log10(f/Hz) values to the solve grid
    as native nodes (so point evaluations at steep spectral features, e.g.
    likelihood bins, do not inherit interpolation error).  Values outside the
    solved frequency range are ignored; ``'construct'`` grids are left
    untouched.
    """
    # Machine-readable reason used by the engine wrapper to distinguish a
    # deterministic physical rejection from a recoverable numerical failure.
    m.fast_failure_reason = None
    if m.cosmo_param['r'] <= 0:
        m.fast_failure_reason = 'invalid_r'
        print('Must set a positive r to calculate the inflationary GWs!')
        return None
    if m.derived_param['N_inf'] is None:
        m.fast_failure_reason = 'invalid_cutoff'
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

    if config is None:
        config = get_config()
    elif not isinstance(config, FastSolverConfig):
        raise TypeError('config must be a FastSolverConfig')
    col_step = config.col_step  # a mid-run compatibility setter must not
    Omega_nu = gp.Omega_nh2/m.derived_param['h']**2    # change array layouts
    # Resolve the frequency-grid builder: per-call argument wins, otherwise the
    # module default set by apply_accuracy_mode()/set_freq_grid().
    if freq_grid is None:
        freq_grid = config.freq_grid
    if freq_grid not in ('construct', 'grid_independent', 'adaptive'):
        raise ValueError('freq_grid must be construct/grid_independent/adaptive, got %r' % freq_grid)
    DN_eff_orig = m.cosmo_param['DN_eff']
    DN_gw_list = [0.0]; DN_gw_new = 0.0; DN_gw_min = 0.0; DN_gw_max = 10.0
    converged = False
    h = config.h
    z_tail = config.z_tail
    phase_max = config.phase_max
    if not 0.0 < tol < 1.0:
        raise ValueError('outer tol must be in (0, 1), got %r' % tol)
    freqs = None; Nf = 0; nv = 0; Nv = None
    idx_out = None; n_coarse = 0
    ev_minus = P_t = fp_freq = Wmat = W_last = None
    Ogw = Oj = Opgw = None
    adaptive_grid = None
    adaptive_done = freq_grid != 'adaptive'
    freq_grid_error = 0.0
    first = True
    thread_before = get_num_threads()
    if config.threads is not None:
        if config.threads > _MAX_THREADS:
            raise ValueError('thread count must be in [1, %d], got %d'
                             % (_MAX_THREADS, config.threads))
        if config.threads != thread_before:
            set_num_threads(config.threads)
    try:
        for _iter in range(MAX_ITER):
            if transition_refine:
                from .exact_background import build_kink_refined_grid
                build_kink_refined_grid(m, h)
            else:
                gen_fast(m, h)
            if freq_grid == 'grid_independent':
                from .freq_adaptive import grid_independent_freqs
                m.f = grid_independent_freqs(m, freq_res)[0]
            elif freq_grid == 'adaptive':
                from .freq_adaptive import grid_independent_freqs
                if adaptive_grid is None:
                    adaptive_grid = grid_independent_freqs(m, freq_res)[0]
                m.f = np.sort(np.asarray(adaptive_grid, dtype=float))[::-1]
            else:
                m.construct_f(freq_res)
            if eval_freqs is not None and freq_grid != 'construct':
                ef = np.asarray(eval_freqs, dtype=float)
                fmin_g, fmax_g = float(np.min(m.f)), float(np.max(m.f))
                ef = ef[(ef >= fmin_g) & (ef <= fmax_g)]
                if ef.size:
                    m.f = np.sort(np.concatenate(
                        (np.asarray(m.f, dtype=float), ef)))[::-1]
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
            if transition_refine:
                from .exact_background import exact_phi_s2_grid
                Phi_grid, Phi_mid, S2, S2inv, h_arr = exact_phi_s2_grid(
                    m, Nv, m.cosmo_param['DN_eff'])
            elif sigma_exact:
                from .exact_background import exact_phi_s2
                Phi_grid, Phi_mid, S2, S2inv = exact_phi_s2(
                    m, Nv, m.cosmo_param['DN_eff'], h)
                h_arr = np.diff(Nv)
            else:
                h_arr = np.diff(Nv)
            if Ogw is None or Ogw.shape[0] != Nf or Ogw.shape[1] != n_coarse:
                # Zero-fill, not np.empty: solve_kernel starts each channel at
                # j0 = horizon-crossing + 3 decades, so the early columns
                # (N < N_j0) are never written.  Those cells are read by
                # int_SGWB_W and physically equal 0 (modes still far outside
                # the horizon); np.empty made the result depend on stale
                # heap contents and could occasionally leak NaN into g2/w2.
                Ogw = np.zeros((Nf, n_coarse)); Oj = np.zeros((Nf, n_coarse)); Opgw = np.zeros((Nf, n_coarse))
                handoff_eps = np.full(Nf, -1.0)
            # phase_max caps the per-(sub-)step phase increment dTheta = e^z dh
            # (horizon-crossing adaptive step control); Sv supplies sigma at the
            # handoff node for the damping-corrected WKB amplitude, handoff_eps
            # receives the per-mode adiabaticity error |1.5*sigma-1|*e^{-z}.
            solve_kernel(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus, fp_minus, fp_freq,
                         1, n_coarse, col_step, h, z_tail, Ogw, Oj, Opgw, h_arr,
                         m.sigma, phase_max, handoff_eps)
            g2_last = np.dot(W_last, (Ogw[:, -1] - Oj[:, -1])[::-1]) * ln10
            DN_gw_new = gp.Neff0 * g2_last / Omega_nu
            if not math.isfinite(DN_gw_new):
                m.fast_failure_reason = 'nonfinite'
                print('SGWB_iter_fast: non-finite DN_gw_new, aborting.')
                break
            if DN_eff_orig + DN_gw_new > 5:
                m.fast_failure_reason = 'shared_Neff_guard'
                print('SGWB_iter_fast: total N_eff too large, aborting.')
                break
            if freq_grid == 'adaptive' and not adaptive_done:
                # Quadrature-weighted frequency refinement: PCHIP second
                # derivative of log10(Omega_GW) estimates the local interval
                # interpolation error |y''| h^2 / 8; weighting by the local
                # integrand share 0.5(Om_i+Om_{i+1}) h_i / sum(w) makes the
                # criterion sensitive to the Delta N_eff integral, not to
                # weightless low-amplitude tail structure.
                Om_cur = np.maximum(Ogw[:, -1] - Oj[:, -1], 1e-300)
                lo_cur = np.log10(Om_cur)
                srt = np.argsort(freqs)
                x = freqs[srt]; y = lo_cur[srt]; Om = Om_cur[srt]
                from scipy.interpolate import PchipInterpolator
                c2 = np.abs(np.asarray(PchipInterpolator(x, y).derivative(2)(x)))
                hh = np.diff(x)
                err_i = c2[:-1]*hh*hh/8.0
                w_i = 0.5*(Om[:-1] + Om[1:])*hh
                wsum = float(np.sum(w_i))
                frac_i = err_i*w_i/wsum if wsum > 0 else np.zeros_like(err_i)
                freq_grid_error = float(np.max(frac_i)) if frac_i.size else 0.0
                if freq_grid_error > freq_grid_target and adaptive_grid.size < freq_grid_max_points:
                    mids = 0.5*(x[:-1] + x[1:])
                    add = mids[frac_i > freq_grid_target]
                    if add.size:
                        room = freq_grid_max_points - int(adaptive_grid.size)
                        if room <= 0:
                            adaptive_done = True
                        else:
                            add = np.sort(add)
                            if add.size > room:
                                add = add[np.linspace(0, add.size-1, room).astype(int)]
                            adaptive_grid = np.unique(np.concatenate([adaptive_grid, add]))
                            continue
                adaptive_done = True
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
            m.fast_failure_reason = 'max_iter'
            print('SGWB_iter_fast: did not converge within %d iterations.' % MAX_ITER)
    finally:
        if config.threads is not None and config.threads != thread_before:
            set_num_threads(thread_before)
        if not converged:
            m.cosmo_param['DN_eff'] = DN_eff_orig
            m.DN_eff_orig = None
            m.SGWB_converge = False
    if not converged:
        if m.fast_failure_reason is None:
            m.fast_failure_reason = 'failed'
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
    # Local a-posteriori quadrature error: Richardson estimate from the
    # Simpson-vs-trapezoid difference on the final frequency grid,
    # |I_simp - I_trap|/15 (standard result for a smooth integrand).
    _Om = m.Ogw_today - m.Oj_today
    _rev = _Om[::-1]
    _hf = np.diff(np.flip(freqs))
    _wt = np.empty(Nf)
    if Nf == 2:
        _wt[0] = 0.5*_hf[0]; _wt[1] = 0.5*_hf[0]
    else:
        _wt[0] = 0.5*_hf[0]; _wt[-1] = 0.5*_hf[-1]
        _wt[1:-1] = 0.5*(_hf[:-1] + _hf[1:])
    _I_simp = float(np.dot(W_last, _rev))*ln10
    _I_trap = float(np.dot(_wt, _rev))*ln10
    m.quadrature_error_local = abs(_I_simp - _I_trap)/15.0/max(abs(_I_simp), 1e-300)
    # Local floating-point/cancellation: Ogw = Oj + remainder near the peak;
    # eps64 * max|Oj|/|Ogw-Oj| bounds the relative subtraction error.
    _num = float(np.max(np.abs(m.Oj_today)/np.maximum(np.abs(_Om), 1e-300)))
    m.cancellation_ratio = _num
    m.floating_point_error = 2.2e-16*_num + 2.2e-16*math.sqrt(Nf)
    # Physics telemetry for the Cobaya adapter / error budget:
    # handoff_eps = per-mode adiabaticity error at the WKB handoff node,
    # z_tail_used/phase_max_used = the frozen-tail threshold and the max
    # phase increment that were actually used, freq_grid_* = the adaptive
    # frequency-grid status (used builder, final quadrature-weighted error
    # estimate, point count).
    m.handoff_eps = np.asarray(handoff_eps)
    m.z_tail_used = z_tail
    m.phase_max_used = phase_max
    m.freq_grid_used = freq_grid
    m.freq_grid_error = float(freq_grid_error)
    m.freq_grid_n = int(Nf)
    m.freq_res_used = float(freq_res)
    m.sigma_exact_used = bool(sigma_exact)
    m.transition_refine_used = bool(transition_refine)
    m.dn_converged_delta = abs(DN_gw_new - DN_gw_list[-1])
    return m


def SGWB_iter_fast(m, tol=1e-4, freq_res=1.0, sigma_exact=False,
                   transition_refine=False, freq_grid=None, config=None,
                   freq_grid_target=3e-4, freq_grid_max_points=1500, eval_freqs=None):
    """Run the fast solver with an isolated configuration snapshot.

    Numba's thread-count setter is process-wide. Calls that use the legacy
    snapshot or a configured thread count therefore hold a short-lived lock
    for the complete solve, preventing concurrent calls from restoring one
    another's thread setting. A ``FastSolverConfig(threads=None)`` call does
    not touch that process-wide state and can proceed without this lock.
    """
    if config is not None and not isinstance(config, FastSolverConfig):
        raise TypeError('config must be a FastSolverConfig')
    lock = (_NUMBA_CONFIG_LOCK if config is None or config.threads is not None
            else nullcontext())
    with lock:
        return _SGWB_iter_fast_impl(
            m, tol=tol, freq_res=freq_res, sigma_exact=sigma_exact,
            transition_refine=transition_refine, freq_grid=freq_grid,
            config=config, freq_grid_target=freq_grid_target,
            freq_grid_max_points=freq_grid_max_points, eval_freqs=eval_freqs)
