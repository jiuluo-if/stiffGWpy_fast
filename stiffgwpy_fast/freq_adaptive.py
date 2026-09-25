# -*- coding: utf-8 -*-
"""Frequency-grid builders for curvature-adaptive and goal-oriented sampling.

The adaptive builder refines ``log10(f)`` where an evaluated spectrum is
curved; its ``evaluate`` callback is supplied by the calling engine, and
already-solved nodes are reused. The formal fast profile uses the separate
goal-oriented builder, which places native nodes around relevant features and
can include requested ``eval_freqs``. See ``docs/numerical_method.md`` and
``docs/accuracy.md`` for the current grid contract and scoped evidence.

中文说明：本模块提供曲率自适应网格和面向目标的 fast 网格构造器。自适应流程由调用方传入
频谱求值函数，并复用已求解节点；正式 fast 档位使用独立的 goal 网格，在关注的频谱特征附近
布点，同时保留显式请求的 ``eval_freqs`` 原生节点。当前网格契约和有限范围的精度证据见
``docs/numerical_method.md`` 与 ``docs/accuracy.md``。
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ['adapt_refine_grid', 'adaptive_spectrum_reference',
           'grid_independent_freqs', 'goal_oriented_freqs']


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
                                z_tail=5.0, rtol=1e-11, workers=1,
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

    ``construct_f`` reads ``m.f_hor`` (the *grid* array), so its sampling can
    shift when the sigma-grid resolution changes. This builder instead derives
    ``fmax`` (horizon at inflation start), ``fmin`` (horizon today), and the
    reheating feature directly from the continuous background. Its frequency
    set is therefore invariant to sigma-grid resolution; see
    ``docs/numerical_method.md`` and the scoped evidence in ``docs/accuracy.md``.

    中文：本函数从连续背景计算频率端点和再加热特征，使同一物理参数下的频率集合不随
    sigma 网格分辨率变化。端点缓存只在本次网格构造调用中有效，不跨外层 `DN_eff` 迭代复用。
    """
    import math as _m

    from . import global_param as gp
    from .exact_background import H2_vec

    d = m.derived_param
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    dn = m.cosmo_param['DN_eff']
    # 与 N 无关的连续背景量只求值一次；原实现每次调用 f_hor_cont 都重算它们。
    h2_cache = {}

    def h2_at(N):
        N = float(N)
        if N not in h2_cache:
            h2_cache[N] = float(H2_vec(np.array([N]), m, dn)[0])
        return h2_cache[N]

    H2_last = h2_at(N_inf)
    raw_last = -0.5 * N_inf + 0.5 * _m.log(H2_last)
    H2_re = h2_at(N_re_abs)
    raw_re = -0.5 * N_re_abs + 0.5 * _m.log(H2_re)
    Delta_f = _m.log(2.0 * _m.pi / d['H_0'])
    ln10v = _m.log(10.0)

    def f_hor_cont(N):
        # Continuous log10(aH/(2pi)/Hz) at absolute N (grid-independent).
        N = float(N)
        H2 = h2_at(N)
        raw = -0.5 * N + 0.5 * _m.log(H2)
        if N < N_re_abs:
            raw = raw_re - 0.5 * (N - N_re_abs)
        return (raw - raw_last - Delta_f) / ln10v

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


def goal_oriented_freqs(m, freq_res=1.0, seed_n=64, max_points=120,
                        eval_freqs=None, feature_half_width=1.0,
                        feature_n=17):
    """Construct a sparse grid around the reheating spectral feature.

    The broad seed is deliberately small for the MCMC path.  A fixed local
    reserve around ``m.f_re`` protects the transition knee, and caller-supplied
    likelihood nodes are mandatory native points.  This builder only chooses
    coordinates; it never evaluates or interpolates the spectrum.
    """
    seed_n = max(2, int(seed_n))
    max_points = max(seed_n, int(max_points))
    feature_n = max(3, int(feature_n))
    base, fmin, fmax = grid_independent_freqs(m, freq_res)
    if base.size <= seed_n:
        seed = base
    else:
        idx = np.rint(np.linspace(0, base.size - 1, seed_n)).astype(int)
        seed = base[np.unique(idx)]

    f_re = float(getattr(m, 'f_re', 0.5 * (fmin + fmax)))
    lo = max(float(fmin), f_re - float(feature_half_width))
    hi = min(float(fmax), f_re + float(feature_half_width))
    feature = np.linspace(hi, lo, feature_n)

    if eval_freqs is None:
        eval = np.empty(0, dtype=float)
    else:
        eval = np.asarray(eval_freqs, dtype=float).reshape(-1)
        eval = eval[np.isfinite(eval)]
        eval = eval[(eval >= fmin) & (eval <= fmax)]

    mandatory = np.unique(np.concatenate((feature, eval)))
    if mandatory.size >= max_points:
        # Preserve every caller node；保留调用方节点，仅在 unusually large likelihood node set
        # consumes the budget 时减少可选 feature reserve。
        selected = np.unique(np.concatenate((eval, feature)))
        if selected.size > max_points and eval.size < max_points:
            selected = np.unique(np.concatenate((eval, feature[:max_points - eval.size])))
        return np.sort(selected)[::-1]

    optional = seed[(seed < lo) | (seed > hi)]
    room = max_points - mandatory.size
    if optional.size > room:
        idx = np.rint(np.linspace(0, optional.size - 1, room)).astype(int)
        optional = optional[np.unique(idx)]
    selected = np.unique(np.concatenate((mandatory, optional)))
    return np.sort(selected)[::-1]
