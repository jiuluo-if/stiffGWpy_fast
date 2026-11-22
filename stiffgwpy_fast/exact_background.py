# -*- coding: utf-8 -*-
"""Continuous-background primitives used while preparing the fast solver.

This module evaluates ``sigma(N)`` and ``H^2(N)`` from the piecewise physical
background, and constructs the corresponding ``Phi``/``S2`` arrays. Reheating
boundaries are treated as explicit breakpoints so integration does not smooth
across the instantaneous-reheating kink. The resulting arrays feed the existing
Numba tensor-stepping kernel; this module does not replace that ODE stepper.
See ``docs/accuracy.md`` for the measured, scope-limited accuracy evidence.

中文说明：本模块在 fast 求解准备阶段计算连续背景的 ``sigma(N)``、``H^2(N)`` 及
``Phi``/``S2`` 原语；积分会把再加热边界作为显式断点，不跨越瞬时再加热 kink 平滑插值。
生成的数组仍交给现有 Numba 张量传播内核使用；这里不替换 ODE 步进算法。精度结论见
``docs/accuracy.md``，须按其中的实验范围解读。
"""

import math

import numpy as np
from numba import njit

from . import global_param as gp

__all__ = ['sigma_vec', 'H2_vec', 'exact_phi_s2', 'build_transition_grid',
           'build_kink_refined_grid', 'exact_phi_s2_grid',
           'exact_phi_s2_breakpoint', 'exact_phi_s2_split',
           'fast_phi_s2_split']

ln10 = math.log(10.0)
_LOG_FLOAT_MAX = math.log(float(np.finfo(np.float64).max))
_LOG_FLOAT_MIN_SUBNORMAL = math.log(float(np.nextafter(0.0, 1.0)))
_PSI_MIN_FLOAT_RANGE = max(_LOG_FLOAT_MIN_SUBNORMAL, -2.0 * _LOG_FLOAT_MAX)
_PSI_MAX_FLOAT_RANGE = min(_LOG_FLOAT_MAX, -2.0 * _LOG_FLOAT_MIN_SUBNORMAL)
_H2_IDENTITY_PANEL_VARIATION_LIMIT = 5.0e-2
_H2_CLOSURE_REL_ROUNDOFF_LIMIT = math.sqrt(np.finfo(np.float64).eps)
_LOG_FLOAT64_ROUNDOFF_MARGIN = -math.log(np.finfo(np.float64).eps)


def _h2_closure_needs_scaling(Omh2, Osh2, h2_today_anchor):
    """Estimate whether matter/stiff closure subtraction loses sqrt(eps).

    中文：按物质和刚性组分闭合相消的条件数估算，只有舍入误差界超过 sqrt(eps)
    时才选择稳定恒等式，避免在可信热路径增加指数计算。
    """
    cancellation_scale = abs(Omh2) + abs(Osh2)
    return (np.finfo(np.float64).eps * cancellation_scale
            > _H2_CLOSURE_REL_ROUNDOFF_LIMIT * abs(h2_today_anchor))


def _stable_closed_h2_base(Omh2, Osh2, h2_today_anchor, delta_N):
    """Evaluate matter + stiff + Lambda without subtracting large closures.

    For ``x = N_inf - N`` and ``q = exp(-3x)``, the exact identity is
    ``Om*(1-q) + Osh*(exp(3x)-q) + anchor*q``. ``expm1`` and ``sinh`` retain
    the small-x differences needed at the present-day closure.
    """
    exponent = 3.0 * delta_N
    decay = np.exp(-exponent)
    with np.errstate(over="ignore", invalid="ignore"):
        stiff_difference = 2.0 * Osh2 * np.sinh(exponent)
    overflow = ~np.isfinite(stiff_difference)
    if np.any(overflow):
        stiff_term = _scaled_exp_product(Osh2, exponent[overflow])
        stiff_difference[overflow] = stiff_term - Osh2 * decay[overflow]
    return (Omh2 * (-np.expm1(-3.0 * delta_N))
            + stiff_difference
            + h2_today_anchor * decay)


def _scaled_exp_product(coefficient, exponent):
    """Evaluate ``coefficient * exp(exponent)`` without an overflowing exp.

    The logarithmic path is only used when the exponential itself is outside
    float64 range; a finite product is retained when the coefficient rescales it.
    """
    exponent = np.asarray(exponent, dtype=np.float64)
    if coefficient == 0.0:
        return np.zeros_like(exponent)
    log_magnitude = math.log(abs(float(coefficient))) + exponent
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        magnitude = np.exp(log_magnitude)
    return np.copysign(magnitude, coefficient)


def _stiff_h2_terms(Osh2, e3N, delta_N):
    """Return the stiff H2 and sigma terms, scaling overflowing exp(3x)."""
    with np.errstate(over="ignore", invalid="ignore"):
        stiff_h2 = Osh2 * e3N
        stiff_sigma = 2.0 * Osh2 * e3N
    overflow = ~np.isfinite(e3N)
    if np.any(overflow):
        stiff_h2[overflow] = _scaled_exp_product(Osh2, 3.0 * delta_N[overflow])
        with np.errstate(over="ignore", invalid="ignore"):
            stiff_sigma[overflow] = 2.0 * stiff_h2[overflow]
    return stiff_h2, stiff_sigma


def _sigma_over_h2(other_numerator, stiff_h2, stiff_sigma, h2,
                   delta_N, Osh2, Omh2, h2_today_anchor, OLh2_base):
    """Form ``(other + 2*stiff)/H2`` without overflowing a finite ratio."""
    other_numerator = np.asarray(other_numerator, dtype=np.float64)
    stiff_h2 = np.asarray(stiff_h2, dtype=np.float64)
    stiff_sigma = np.asarray(stiff_sigma, dtype=np.float64)
    h2 = np.asarray(h2, dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        numerator = other_numerator + stiff_sigma
        ratio = numerator / h2
    recover = (
        ~np.isfinite(numerator)
        & np.isfinite(other_numerator)
        & np.isfinite(stiff_h2)
        & np.isfinite(h2)
        & (h2 > 0.0)
    )
    if np.any(recover):
        ratio[recover] = (
            other_numerator[recover] / h2[recover]
            + 2.0 * (stiff_h2[recover] / h2[recover])
        )
    unrepresented = ~np.isfinite(ratio) & ~np.isfinite(h2)
    if Osh2 > 0.0:
        for index in np.flatnonzero(unrepresented):
            delta = float(delta_N[index])
            if delta < 0.0:
                continue
            log_stiff = math.log(float(Osh2)) + 3.0 * delta
            terms = [
                float(other_numerator[index]), float(Omh2),
                float(h2_today_anchor) * math.exp(-3.0 * delta),
                float(OLh2_base) * math.exp(-3.0 * delta),
            ]
            log_other = max(
                (math.log(abs(term)) for term in terms if term != 0.0),
                default=-math.inf,
            )
            if log_stiff - log_other > _LOG_FLOAT64_ROUNDOFF_MARGIN:
                ratio[index] = 2.0
    return ratio


def _fd_from_ref():
    from .reference import _FD_P_CS, _FD_RHO_CS
    return _FD_RHO_CS, _FD_P_CS


def sigma_vec(N, m, DN_eff):
    """Vectorized exact ``sigma`` at an array of absolute ``N`` values.

    中文：对一组绝对 e-fold 坐标 ``N`` 向量化计算连续背景的精确 ``sigma``。
    """
    N = np.asarray(N, dtype=float)
    d = m.derived_param
    Omh2 = d['Omega_mh2']
    Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
    Otrh2_base = gp.Omega_orh2
    Otreh2_base = gp.Omega_ph2 * gp.rho_th[-1]
    OLh2_base = d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 \
        - gp.Omega_ph2 - Osh2
    h2_today_anchor = d['h'] ** 2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 - gp.Omega_ph2
    stable_closure = _h2_closure_needs_scaling(Omh2, Osh2, h2_today_anchor)
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    fd_rho, fd_p = _fd_from_ref()

    sigma = np.ones_like(N)
    running = N >= N_re_abs
    Nr = N[running]
    delta_N = N_inf - Nr
    eN = np.exp(delta_N)
    with np.errstate(over="ignore", invalid="ignore"):
        e3N = eN * eN * eN
    stiff_h2, stiff_sigma = _stiff_h2_terms(Osh2, e3N, delta_N)
    overflowed_e3 = ~np.isfinite(e3N)
    decay = np.exp(-3.0 * delta_N)
    with np.errstate(divide="ignore", invalid="ignore"):
        inverse_e3 = 1.0 / e3N
        lambda_term = OLh2_base / e3N
    if np.any(overflowed_e3):
        inverse_e3[overflowed_e3] = decay[overflowed_e3]
        lambda_term[overflowed_e3] = OLh2_base * decay[overflowed_e3]
    if stable_closure:
        closed_h2_base = _stable_closed_h2_base(
            Omh2, Osh2, h2_today_anchor, delta_N
        )
    extra_radiation_h2 = Oerh2 * (eN - inverse_e3)
    nu = gp.nu_today / eN
    Omega_ph2 = gp.Omega_ph2
    Omega_nh2 = gp.Omega_nh2
    Omega_mnuh2 = gp.Omega_mnuh2

    m1 = nu > 100.0
    out = np.empty(running.sum())
    if m1.any():
        if stable_closure:
            H2 = (closed_h2_base[m1] + Omega_mnuh2
                  + (Omega_ph2 + 2.0 / 3.0 * Omega_nh2) * eN[m1]
                  + extra_radiation_h2[m1])
        else:
            H2 = (Omh2 + Omega_mnuh2 + (Omega_ph2 + 2.0 / 3.0 * Omega_nh2) * eN[m1]
                  + stiff_h2[m1] + extra_radiation_h2[m1] + lambda_term[m1])
        sigma_other = (Omh2 + Omega_mnuh2
                       + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN[m1])
        out[m1] = _sigma_over_h2(
            sigma_other, stiff_h2[m1], stiff_sigma[m1], H2,
            delta_N[m1], Osh2, Omh2, h2_today_anchor, OLh2_base,
        )
    m2 = (~m1) & (nu >= 0.1)
    if m2.any():
        lognu = np.log10(nu[m2])
        rho_nu = fd_rho(lognu)
        p_nu = fd_p(lognu)
        if stable_closure:
            H2 = (closed_h2_base[m2]
                  + (Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * Omega_nh2) * eN[m2]
                  + extra_radiation_h2[m2])
        else:
            H2 = (Omh2 + (Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * Omega_nh2) * eN[m2]
                  + stiff_h2[m2] + extra_radiation_h2[m2] + lambda_term[m2])
        sigma_other = (Omh2
                       + 4.0 / 3.0 * (Omega_ph2 + 2.0 / 3.0 * Omega_nh2 + Oerh2) * eN[m2]
                       + (rho_nu + p_nu) * Omega_nh2 / 3.0 * eN[m2])
        out[m2] = _sigma_over_h2(
            sigma_other, stiff_h2[m2], stiff_sigma[m2], H2,
            delta_N[m2], Osh2, Omh2, h2_today_anchor, OLh2_base,
        )
    m3 = (~m1) & (~m2) & (Nr > N_inf - gp.N_fin)
    if m3.any():
        if stable_closure:
            H2 = closed_h2_base[m3] + Otrh2_base * eN[m3] + extra_radiation_h2[m3]
        else:
            H2 = (Omh2 + Otrh2_base * eN[m3] + stiff_h2[m3]
                  + extra_radiation_h2[m3] + lambda_term[m3])
        sigma_other = Omh2 + 4.0 / 3.0 * (Otrh2_base + Oerh2) * eN[m3]
        out[m3] = _sigma_over_h2(
            sigma_other, stiff_h2[m3], stiff_sigma[m3], H2,
            delta_N[m3], Osh2, Omh2, h2_today_anchor, OLh2_base,
        )
    m4 = (~m1) & (~m2) & (~m3) & (Nr >= N_inf - gp.N_max)
    if m4.any():
        Nl = N_inf - Nr[m4]
        rho_i = gp.spl_rho(Nl)
        rhop_i = gp.spl_rhop(Nl)
        if stable_closure:
            H2 = (closed_h2_base[m4] + Omega_ph2 * rho_i * eN[m4]
                  + extra_radiation_h2[m4])
        else:
            H2 = (Omh2 + Omega_ph2 * rho_i * eN[m4] + stiff_h2[m4]
                  + extra_radiation_h2[m4] + lambda_term[m4])
        sigma_other = Omh2 + (Omega_ph2 * rhop_i + 4.0 / 3.0 * Oerh2) * eN[m4]
        out[m4] = _sigma_over_h2(
            sigma_other, stiff_h2[m4], stiff_sigma[m4], H2,
            delta_N[m4], Osh2, Omh2, h2_today_anchor, OLh2_base,
        )
    m5 = (~m1) & (~m2) & (~m3) & (~m4)
    if m5.any():
        if stable_closure:
            H2 = closed_h2_base[m5] + Otreh2_base * eN[m5] + extra_radiation_h2[m5]
        else:
            H2 = (Omh2 + Otreh2_base * eN[m5] + stiff_h2[m5]
                  + extra_radiation_h2[m5] + lambda_term[m5])
        sigma_other = Omh2 + 4.0 / 3.0 * (Otreh2_base + Oerh2) * eN[m5]
        out[m5] = _sigma_over_h2(
            sigma_other, stiff_h2[m5], stiff_sigma[m5], H2,
            delta_N[m5], Osh2, Omh2, h2_today_anchor, OLh2_base,
        )
    sigma[running] = out
    return sigma


def H2_vec(N, m, DN_eff):
    """Vectorized exact ``H^2`` at an array of absolute ``N``.

    中文：对一组绝对坐标 ``N`` 向量化计算精确 ``H^2``，输入与 ``sigma_vec`` 相同。
    """
    N = np.asarray(N, dtype=float)
    d = m.derived_param
    Omh2 = d['Omega_mh2']
    Osh2 = d['Omega_sh2']
    Oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * DN_eff
    Otrh2_base = gp.Omega_orh2
    Otreh2_base = gp.Omega_ph2 * gp.rho_th[-1]
    OLh2_base = d['h'] ** 2 - Omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 \
        - gp.Omega_ph2 - Osh2
    h2_today_anchor = d['h'] ** 2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3 - gp.Omega_ph2
    stable_closure = _h2_closure_needs_scaling(Omh2, Osh2, h2_today_anchor)
    N_inf = d['N_inf']
    fd_rho, _ = _fd_from_ref()
    # Preserve the legacy scalar ``math.exp`` result across platforms while
    # batching the much more expensive spline evaluation below.
    # 保留旧版逐点 ``math.exp`` 的跨平台结果，仅批量处理耗时更高的样条插值。
    eN = np.fromiter(
        (math.exp(N_inf - float(Ni)) for Ni in N),
        dtype=float,
        count=N.size,
    )
    delta_N = N_inf - N
    with np.errstate(over="ignore", invalid="ignore"):
        e3N = eN * eN * eN
    stiff_h2, _ = _stiff_h2_terms(Osh2, e3N, delta_N)
    overflowed_e3 = ~np.isfinite(e3N)
    decay = np.exp(-3.0 * delta_N)
    with np.errstate(divide="ignore", invalid="ignore"):
        inverse_e3 = 1.0 / e3N
        lambda_term = OLh2_base / e3N
    if np.any(overflowed_e3):
        inverse_e3[overflowed_e3] = decay[overflowed_e3]
        lambda_term[overflowed_e3] = OLh2_base * decay[overflowed_e3]
    extra_radiation_h2 = Oerh2 * (eN - inverse_e3)
    if stable_closure:
        closed_h2_base = _stable_closed_h2_base(
            Omh2, Osh2, h2_today_anchor, delta_N
        )
    nu = gp.nu_today / eN
    out = np.empty_like(N)
    # Keep the scalar branch precedence explicit while evaluating each
    # background regime as one array operation. The old per-node PCHIP call
    # dominated variable-grid setup for dense transition-refined grids.
    # 按原有分支优先级构造互斥掩码，并批量计算各背景区间；旧实现逐节点调用 PCHIP，
    # 在密集的转折细化网格中成为主要耗时。
    high_nu = nu > 100.0
    if np.any(high_nu):
        if stable_closure:
            out[high_nu] = (closed_h2_base[high_nu] + gp.Omega_mnuh2
                            + (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2) * eN[high_nu]
                            + extra_radiation_h2[high_nu])
        else:
            out[high_nu] = (Omh2 + gp.Omega_mnuh2
                            + (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2) * eN[high_nu]
                            + stiff_h2[high_nu] + extra_radiation_h2[high_nu]
                            + lambda_term[high_nu])
    thermal_nu = (~high_nu) & (nu >= 0.1)
    if np.any(thermal_nu):
        rho_nu = fd_rho(np.log10(nu[thermal_nu]))
        if stable_closure:
            out[thermal_nu] = (closed_h2_base[thermal_nu]
                               + (gp.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * gp.Omega_nh2)
                               * eN[thermal_nu] + extra_radiation_h2[thermal_nu])
        else:
            out[thermal_nu] = (Omh2
                               + (gp.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * gp.Omega_nh2)
                               * eN[thermal_nu]
                               + stiff_h2[thermal_nu] + extra_radiation_h2[thermal_nu]
                               + lambda_term[thermal_nu])
    relativistic = (~high_nu) & (~thermal_nu) & (N > N_inf - gp.N_fin)
    if np.any(relativistic):
        if stable_closure:
            out[relativistic] = (closed_h2_base[relativistic]
                                 + Otrh2_base * eN[relativistic]
                                 + extra_radiation_h2[relativistic])
        else:
            out[relativistic] = (Omh2 + Otrh2_base * eN[relativistic]
                                 + stiff_h2[relativistic] + extra_radiation_h2[relativistic]
                                 + lambda_term[relativistic])
    thermal_grid = ((~high_nu) & (~thermal_nu) & (~relativistic)
                    & (N >= N_inf - gp.N_max))
    if np.any(thermal_grid):
        Nl = N_inf - N[thermal_grid]
        rho_nu = gp.spl_rho(Nl)
        if stable_closure:
            out[thermal_grid] = (closed_h2_base[thermal_grid]
                                 + gp.Omega_ph2 * rho_nu * eN[thermal_grid]
                                 + extra_radiation_h2[thermal_grid])
        else:
            out[thermal_grid] = (Omh2 + gp.Omega_ph2 * rho_nu * eN[thermal_grid]
                                 + stiff_h2[thermal_grid] + extra_radiation_h2[thermal_grid]
                                 + lambda_term[thermal_grid])
    early = ~(high_nu | thermal_nu | relativistic | thermal_grid)
    if np.any(early):
        if stable_closure:
            out[early] = (closed_h2_base[early] + Otreh2_base * eN[early]
                          + extra_radiation_h2[early])
        else:
            out[early] = (Omh2 + Otreh2_base * eN[early]
                          + stiff_h2[early] + extra_radiation_h2[early]
                          + lambda_term[early])
    return out


def _h2_nonstiff_abs_bound(m, DN_eff, delta):
    """Bound all non-stiff H2 terms by summing their absolute magnitudes."""
    d = m.derived_param
    omh2 = float(d["Omega_mh2"])
    osh2 = float(d["Omega_sh2"])
    anchor = float(d["h"] ** 2 - gp.Omega_mnuh2
                   - gp.Omega_nh2 * 2.0 / 3.0 - gp.Omega_ph2)
    olh2 = float(d["h"] ** 2 - omh2 - gp.Omega_mnuh2
                 - gp.Omega_nh2 * 2.0 / 3.0 - gp.Omega_ph2 - osh2)
    oerh2 = (gp.Omega_ph2 * 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)
             * float(DN_eff))
    decay = math.exp(-3.0 * delta)
    expansion = math.exp(delta)
    bound = (abs(omh2) + abs(osh2) * decay + abs(anchor) * decay
             + abs(olh2) * decay + abs(oerh2) * decay)
    nu = gp.nu_today / expansion
    if nu > 100.0:
        fluid = (abs(gp.Omega_mnuh2)
                 + abs(gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2) * expansion)
    elif nu >= 0.1:
        fd_rho, _ = _fd_from_ref()
        rho_nu = float(fd_rho(math.log10(nu)))
        fluid = abs(gp.Omega_ph2 + (2.0 / 3.0 + rho_nu / 3.0) * gp.Omega_nh2) * expansion
    elif delta < gp.N_fin:
        fluid = abs(gp.Omega_orh2 + oerh2) * expansion
    elif delta <= gp.N_max:
        fluid = abs(gp.Omega_ph2 * float(gp.spl_rho(delta)) + oerh2) * expansion
    else:
        fluid = abs(gp.Omega_ph2 * float(gp.rho_th[-1]) + oerh2) * expansion
    return bound + fluid


def _log_H2_vec(N, m, DN_eff):
    """Return log(H2), retaining positive stiff-dominated values above float64.

    H2 itself is still reported as nonfinite by :func:`H2_vec` when it cannot
    be represented.  Fast background consumers need its logarithm for sigma
    ratios and horizon frequencies, so this helper recovers only the regime
    where a conservative bound proves that the positive stiff term dominates
    every other term by more than float64 roundoff.
    """
    values = np.asarray(N, dtype=np.float64)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        h2 = H2_vec(values, m, DN_eff)
    result = np.full(values.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(h2) & (h2 > 0.0)
    result[valid] = np.log(h2[valid])
    needs_log_scale = ~np.isfinite(h2)
    if not np.any(needs_log_scale):
        return result

    d = m.derived_param
    osh2 = float(d["Omega_sh2"])
    if osh2 <= 0.0:
        return result
    log_roundoff_margin = -math.log(np.finfo(np.float64).eps)
    n_inf = float(d["N_inf"])
    for index in np.flatnonzero(needs_log_scale):
        delta = n_inf - float(values[index])
        if delta < 0.0:
            continue
        log_stiff = math.log(osh2) + 3.0 * delta
        other_bound = _h2_nonstiff_abs_bound(m, DN_eff, delta)
        if not math.isfinite(other_bound) or other_bound < 0.0:
            continue
        log_other_bound = math.log(other_bound) if other_bound > 0.0 else -math.inf
        if log_stiff - log_other_bound > log_roundoff_margin:
            result[index] = log_stiff
    return result


def _sigma_node_limits(Nv, m, DN_eff, sigma_nodes=None):
    """Return continuous sigma at grid nodes with explicit kink limits.

    ``sigma_vec`` evaluates the reheating point using the post-transition
    branch.  Simpson panels ending at that point must instead use the
    pre-transition limit, while panels starting there must use the post-
    transition limit.  Keeping both arrays makes the convention explicit.

    中文：节点值按转变后的分支计算；若 Simpson 面板终止于再加热点，该面板需用转变前极限，
    从该点开始的面板则用转变后极限。分别返回两组单侧极限，明确积分约定。
    """
    Nv = np.asarray(Nv, dtype=float)
    if sigma_nodes is None:
        nodes = sigma_vec(Nv, m, DN_eff)
    else:
        nodes = np.asarray(sigma_nodes, dtype=float)
        if nodes.shape != Nv.shape:
            raise ValueError('sigma_nodes must match Nv shape')
    left = nodes.copy()
    right = nodes.copy()
    d = m.derived_param
    n_re = float(d['N_inf'] - d['N_re'])
    at_re = np.isclose(Nv, n_re, rtol=0.0, atol=1e-12)
    if np.any(at_re):
        left[at_re] = 1.0
        right[at_re] = sigma_vec(np.array([n_re]), m, DN_eff)[0]
    return left, right


def _primitive_scales(psi, relative_scale):
    psi = np.asarray(psi, dtype=np.float64)
    psi_min = float(np.min(psi))
    psi_max = float(np.max(psi))
    auto_scale = bool(
        psi_min < _LOG_FLOAT_MIN_SUBNORMAL or psi_max > _LOG_FLOAT_MAX
    )
    use_relative = bool(relative_scale or auto_scale)
    if use_relative:
        units = np.ones_like(psi)
        return units, units, use_relative, auto_scale and not relative_scale
    return np.exp(psi), np.exp(-0.5 * psi), False, False


def _h2_identity_primitives(m, Nv, DN_eff, midpoints, *, relative_scale=False):
    """Build exact background primitives from the smooth-branch H2 identity.

    On every smooth post-reheating branch, ``sigma = 1 - (d ln H2/dN)/3``.
    Integrate that identity from the physical reheating breakpoint; retain the
    exact pre-reheating matter branch ``sigma=1``. The transition itself is not
    smoothed or shifted.
    """
    nodes = np.asarray(Nv, dtype=np.float64)
    mids = np.asarray(midpoints, dtype=np.float64)
    n0 = float(nodes[0])
    n_re = float(m.derived_param["N_inf"] - m.derived_param["N_re"])
    anchor_n = max(n0, n_re)
    post_nodes = nodes >= n_re
    post_mids = mids >= n_re
    try:
        log_h2_anchor = float(_log_H2_vec(np.asarray([anchor_n]), m, DN_eff)[0])
        log_h2_nodes = _log_H2_vec(nodes[post_nodes], m, DN_eff)
        log_h2_mids = _log_H2_vec(mids[post_mids], m, DN_eff)
    except (OverflowError, ValueError, FloatingPointError):
        log_h2_anchor = float("nan")
        log_h2_nodes = np.full(int(np.count_nonzero(post_nodes)), np.nan)
        log_h2_mids = np.full(int(np.count_nonzero(post_mids)), np.nan)

    valid_anchor = math.isfinite(log_h2_anchor)
    valid_nodes = np.isfinite(log_h2_nodes).all()
    valid_mids = np.isfinite(log_h2_mids).all()
    if not (valid_anchor and valid_nodes and valid_mids):
        failure = "nonfinite_background"

        def is_nonpositive_h2(N):
            try:
                with np.errstate(over="ignore", invalid="ignore"):
                    value = float(H2_vec(np.asarray([N]), m, DN_eff)[0])
            except (OverflowError, ValueError, FloatingPointError):
                return False
            return math.isfinite(value) and value <= 0.0

        if not valid_anchor:
            failure_index = int(np.searchsorted(nodes, anchor_n, side="left"))
            if is_nonpositive_h2(anchor_n):
                failure = "nonpositive_h2"
        elif not valid_nodes:
            bad_index = int(np.flatnonzero(~np.isfinite(log_h2_nodes))[0])
            failure_index = int(np.flatnonzero(post_nodes)[bad_index])
            if is_nonpositive_h2(nodes[failure_index]):
                failure = "nonpositive_h2"
        else:
            bad_index = int(np.flatnonzero(~np.isfinite(log_h2_mids))[0])
            failure_index = max(0, int(np.searchsorted(nodes, mids[bad_index])))
            if is_nonpositive_h2(mids[bad_index]):
                failure = "nonpositive_h2"
        return None, failure, failure_index

    log_anchor = log_h2_anchor
    F_nodes = nodes - n0
    F_mids = mids - n0
    F_nodes[post_nodes] -= (log_h2_nodes - log_anchor) / 3.0
    F_mids[post_mids] -= (log_h2_mids - log_anchor) / 3.0
    phi_nodes = 1.5 * F_nodes - nodes + n0
    phi_mids = 1.5 * F_mids - mids + n0
    psi = 3.0 * F_nodes - 4.0 * nodes
    s2, s2inv, use_relative, _ = _primitive_scales(psi, relative_scale)
    phi_re = 0.5 * (n_re - n0) if n0 < n_re <= nodes[-1] else 0.0
    return (phi_nodes, phi_mids, s2, s2inv, use_relative, float(phi_re)), None, None


def _replace_overflowed_primitive(
    m, Nv, DN_eff, midpoints, *, auto_overflow, relative_scale,
):
    if not auto_overflow:
        return None
    result, failure, failure_index = _h2_identity_primitives(
        m, Nv, DN_eff, midpoints, relative_scale=relative_scale
    )
    if result is None:
        m._fast_background_failure = failure
        m._fast_background_failure_detail = "H2 in sigma-integral identity"
        m._fast_background_failure_index = int(failure_index)
        return None
    phi, phi_mid, s2, s2inv, relative_used, phi_re = result
    m._fast_relative_background_required = bool(relative_used)
    m._fast_h2_identity_used = True
    return phi, phi_mid, s2, s2inv, float(phi_re), bool(relative_used)


def exact_phi_s2(m, Nv, DN_eff, h, *, relative_scale=False):
    """Exact ``Phi``/``Phi_mid``/``S2``/``S2inv`` from the continuous sigma.

    ``F = integral sigma dN`` is evaluated on a fine sub-grid with the reheating
    boundary as an exact breakpoint, so the kink is not smoothed.  Returns the
    arrays expected by ``fast_sgwb.solve_kernel`` (which uses ``N0=0``).

    中文：根据连续 ``sigma`` 计算 ``Phi``、中点 ``Phi``、``S2`` 和 ``S2inv``。在细网格上积分
    ``F = ∫ sigma dN``，并将再加热边界作为精确分段点，避免平滑掉转折；返回值可直接传给
    ``fast_sgwb.solve_kernel``（其坐标起点为 ``N0=0``）。
    """
    Nv = np.asarray(Nv, dtype=float)
    # Fine sub-grid (a few points per base step) plus its midpoints for Phi_mid.
    # 在每个基础步内加密取点，并加入步中点以计算 Phi_mid。
    half = 0.5 * h
    frac = 0.125 * h
    sub = np.concatenate([
        np.arange(Nv[0], Nv[-1] + 1e-12, frac),
        np.arange(Nv[0] + half, Nv[-1] + 1e-12, h),
    ])
    sub = np.unique(np.sort(np.concatenate((sub, Nv))))
    d = m.derived_param
    N_re_abs = d['N_inf'] - d['N_re']
    if (N_re_abs > sub[0]) and (N_re_abs < sub[-1]):
        sub = np.unique(np.concatenate((sub, [N_re_abs])))
    sig_left, sig_right = _sigma_node_limits(sub, m, DN_eff)
    # Cumulative integral F on the sub-grid (shape-preserving trapezoid; the kink
    # is a breakpoint, so sigma is smooth within each sub-interval).
    # 在子网格上用保形梯形法累计积分 F；kink 作为断点，因此每个子区间内的 sigma 连续光滑。
    F = np.concatenate(([0.0], np.cumsum(0.5 * (sig_right[1:] + sig_left[:-1])
                                             * np.diff(sub))))
    # Interpolate F onto grid nodes and midpoints.
    # 将 F 插值到基础网格节点和步中点。
    F_grid = np.interp(Nv, sub, F)
    F_mid = np.interp(Nv + half, sub, F)
    N0 = Nv[0]
    Phi_grid = 1.5 * F_grid - Nv + N0
    Phi_mid = 1.5 * F_mid - (Nv + half) + N0
    Psi = 3.0 * F_grid - 4.0 * Nv
    S2, S2inv, relative_used, auto_overflow = _primitive_scales(Psi, relative_scale)
    if auto_overflow:
        replacement = _replace_overflowed_primitive(
            m, Nv, DN_eff, Nv + half, auto_overflow=True,
            relative_scale=relative_scale,
        )
        if replacement is not None:
            Phi_grid, Phi_mid, S2, S2inv, _, relative_used = replacement
    m._fast_relative_background_required = bool(relative_used)
    return Phi_grid.astype(np.float64), Phi_mid.astype(np.float64), \
        S2.astype(np.float64), S2inv.astype(np.float64)


def build_transition_grid(m, h):
    """Variable N-grid with a node exactly at the reheating kink.

    The kink becomes a grid node so that no Magnus step straddles the
    ``sigma`` jump (removing the O(h) frozen-midpoint phase error).  Returns
    ``(Nv, sigma, f_hor, index_re)`` and sets the model grid attributes.

    中文：构造可变步长 ``N`` 网格，并确保再加热转折 ``N_re`` 精确落在网格节点上，避免 Magnus 步跨过
    ``sigma`` 跳变及其 ``O(h)`` 冻结中点相位误差。返回 ``(Nv, sigma, f_hor, index_re)``，并写入模型网格属性。
    """
    d = m.derived_param
    N_inf = d['N_inf']
    N_re_abs = N_inf - d['N_re']
    len_inf = math.floor(N_inf / h) + 1
    Nv = np.arange(0, len_inf) * h
    # Keep the continuous present-day anchor used by ``gen_fast``.  Rounding
    # the final node down to the nearest h-grid point would add a separate
    # grid-anchor bias to this Phase-A experiment.
    # 保留 gen_fast 使用的连续当前时刻锚点；将末节点向下取整会在本实验中额外引入网格锚定偏差。
    present = N_inf
    Nv = np.unique(np.sort(np.concatenate((Nv, [N_re_abs, present]))))
    dn = m.cosmo_param['DN_eff']
    sigma = sigma_vec(Nv, m, dn)
    # Exact f_hor (log10 aH/(2pi)/Hz) with the MD extension for N < N_re_abs,
    # normalised to the grid's present node (same convention as gen_fast).
    # 计算精确 f_hor；N < N_re_abs 时延拓物质主导分支，并按网格当前时刻节点归一化，与 gen_fast 一致。
    post = Nv >= N_re_abs
    raw = np.empty_like(Nv)
    raw[post] = -0.5 * Nv[post] + 0.5 * np.log(H2_vec(Nv[post], m, dn))
    raw_re = raw[np.searchsorted(Nv, N_re_abs)]
    raw[~post] = raw_re - 0.5 * (Nv[~post] - N_re_abs)
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

    中文：在均匀步长 ``h`` 网格附近，对 ``N_re_abs`` 左右 ``refine_radius`` 范围用
    ``h/refine_factor`` 加密，并将转折点设为精确节点。返回网格、背景量和转折索引，同时写入模型属性。
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
    # ``arange`` can land a few ulps from an explicit boundary or from a node
    # on the other grid.  Retaining both creates near-zero intervals in the
    # variable-step propagation kernel and ambiguous one-sided kink samples.
    # Drop generated copies around the exact physical anchors, merge remaining
    # generated nodes, then restore both anchors even if they are close to each
    # other; a tiny but representable reheating interval must remain explicit.
    node_tol = 1e-12
    nodes = np.concatenate((Nv, refined))
    near_anchor = ((np.abs(nodes - N_re_abs) <= node_tol)
                   | (np.abs(nodes - present) <= node_tol))
    nodes = np.sort(nodes[~near_anchor])
    if nodes.size:
        keep = np.concatenate(([True], np.diff(nodes) > node_tol))
        nodes = nodes[keep]
    anchors = np.unique(np.array([N_re_abs, present], dtype=float))
    Nv = np.sort(np.concatenate((nodes, anchors)))
    dn = m.cosmo_param['DN_eff']
    sigma = sigma_vec(Nv, m, dn)
    post = Nv >= N_re_abs
    raw = np.empty_like(Nv)
    raw[post] = -0.5 * Nv[post] + 0.5 * np.log(H2_vec(Nv[post], m, dn))
    raw_re = raw[np.searchsorted(Nv, N_re_abs)]
    raw[~post] = raw_re - 0.5 * (Nv[~post] - N_re_abs)
    Delta_f = math.log(2.0 * math.pi / d['H_0'])
    f_hor = (raw - raw[-1] - Delta_f) / math.log(10.0)
    index_re = int(np.argmin(np.abs(Nv - N_re_abs)))
    m.Nv = Nv
    m.N = Nv - Nv[-1]
    m.sigma = sigma
    m.f_hor = f_hor
    m.f_re = f_hor[index_re]
    return Nv, sigma, f_hor, index_re


def exact_phi_s2_grid(m, Nv, DN_eff, *, relative_scale=False):
    """Exact ``Phi``/``midpoint-Phi``/``S2``/``S2inv``/``h_arr`` on a variable grid.

    中文：在可变步长网格上精确计算节点和区间中点处的 ``Phi``，以及 ``S2``、``S2inv`` 和步长数组。
    """
    Nv = np.asarray(Nv, dtype=float)
    d = m.derived_param
    N_re_abs = d['N_inf'] - d['N_re']
    # Fine sub-grid (roughly half the local spacing) for the F integral, with the
    # kink as a breakpoint.  Use the *median* spacing so that a single inserted
    # kink node cannot collapse the step and blow up the sub-grid.
    # 用局部步长约一半的子网格积分 F，并把 kink 设为断点；用中位步长避免单个插入节点压缩步长、撑大网格。
    h = float(np.median(np.abs(np.diff(Nv))))
    step = 0.5 * h
    sub = np.concatenate((np.arange(Nv[0], Nv[-1] + 1e-12, step), Nv))
    if (N_re_abs > sub[0]) and (N_re_abs < sub[-1]):
        sub = np.concatenate((sub, [N_re_abs]))
    sub = np.unique(np.sort(sub))
    sig_left, sig_right = _sigma_node_limits(sub, m, DN_eff)
    F = np.concatenate(([0.0], np.cumsum(0.5 * (sig_right[1:] + sig_left[:-1])
                                             * np.diff(sub))))
    mid = 0.5 * (Nv[:-1] + Nv[1:])
    F_nodes = np.interp(Nv, sub, F)
    F_mid = np.interp(mid, sub, F)
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - mid + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2, S2inv, relative_used, auto_overflow = _primitive_scales(Psi, relative_scale)
    if auto_overflow:
        replacement = _replace_overflowed_primitive(
            m, Nv, DN_eff, mid, auto_overflow=True,
            relative_scale=relative_scale,
        )
        if replacement is not None:
            Phi_grid, Phi_mid, S2, S2inv, _, relative_used = replacement
    m._fast_relative_background_required = bool(relative_used)
    h_arr = np.diff(Nv).astype(np.float64)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64), h_arr)


def exact_phi_s2_breakpoint(m, Nv, DN_eff):
    """Build primitives with composite Simpson quadrature on a split grid.

    Unlike :func:`exact_phi_s2_grid`, this evaluates only the interval
    midpoints in addition to the already-computed node values.  Because the
    caller has inserted ``N_re`` as a node, every Simpson interval stays on
    one side of the reheating kink.  The midpoint primitive uses half of the
    symmetric full-interval Simpson integral and avoids changing the transfer
    kernel's stepping logic.

    中文：在已含精确 ``N_re`` 节点的分段网格上使用复合 Simpson 积分。除已有节点值外只计算区间中点；
    每个 Simpson 区间都位于转折点同一侧。中点原函数取对称全区间 Simpson 积分的一半，无需改动传递核的步进逻辑。
    """
    Nv = np.asarray(Nv, dtype=float)
    h_arr = np.diff(Nv).astype(np.float64)
    nodes_left, nodes_right = _sigma_node_limits(Nv, m, DN_eff)
    mid = 0.5 * (Nv[:-1] + Nv[1:])
    mid_sigma = sigma_vec(mid, m, DN_eff)
    integral = h_arr * (nodes_left[:-1] + 4.0 * mid_sigma + nodes_right[1:]) / 6.0
    quarter_sigma = sigma_vec(Nv[:-1] + 0.25 * h_arr, m, DN_eff)
    F_nodes = np.concatenate(([0.0], np.cumsum(integral)))
    half_integral = h_arr * (nodes_left[:-1] + 4.0 * quarter_sigma + mid_sigma) / 12.0
    F_mid = F_nodes[:-1] + half_integral
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - mid + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2 = np.exp(Psi)
    S2inv = np.exp(-0.5 * Psi)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64), h_arr)


def exact_phi_s2_split(m, Nv, DN_eff, sigma_nodes=None, *, relative_scale=False):
    """Build uniform-grid primitives while splitting the one ``N_re`` interval.

    The returned grid remains uniform, so the hot kernel keeps its original
    array layout.  Composite Simpson quadrature is applied independently to
    both sides when the reheating boundary lies inside an interval; the
    resulting ``kink_index``, ``kink_fraction`` and ``Phi_re`` tell the ODE
    kernel to perform two transfer steps for that interval.

    中文：保持原均匀网格布局，只将包含 ``N_re`` 的区间拆成左右两段分别积分。返回 ``kink_index``、
    ``kink_fraction`` 和 ``Phi_re``，供 ODE 核在该网格区间内执行两次传递步。
    """
    Nv = np.asarray(Nv, dtype=float)
    h_arr = np.diff(Nv).astype(np.float64)
    if Nv.size < 2:
        return (np.zeros_like(Nv), np.zeros_like(Nv), np.ones_like(Nv),
                np.ones_like(Nv), -1, 0.0, 0.0)
    nodes_left, nodes_right = _sigma_node_limits(Nv, m, DN_eff, sigma_nodes)
    mid = 0.5 * (Nv[:-1] + Nv[1:])
    mid_sigma = sigma_vec(mid, m, DN_eff)
    integral = h_arr * (nodes_left[:-1] + 4.0 * mid_sigma + nodes_right[1:]) / 6.0
    quarter_sigma = sigma_vec(Nv[:-1] + 0.25 * h_arr, m, DN_eff)

    n_re = float(m.derived_param['N_inf'] - m.derived_param['N_re'])
    kink_index = int(np.searchsorted(Nv, n_re, side='right') - 1)
    kink_fraction = 0.0
    left_integral = None
    if 0 <= kink_index < Nv.size - 1:
        left = float(Nv[kink_index])
        right = float(Nv[kink_index + 1])
        kink_fraction = (n_re - left) / (right - left)
        if 0.0 < kink_fraction < 1.0:
            probes = np.array([
                left + 0.5 * (n_re - left),
                n_re,
                n_re + 0.5 * (right - n_re),
            ], dtype=float)
            sig_left, sig_re, sig_right = sigma_vec(probes, m, DN_eff)
            left_h = n_re - left
            right_h = right - n_re
            left_integral = left_h * (nodes_left[kink_index] +
                                      4.0 * sig_left + 1.0) / 6.0
            right_integral = right_h * (sig_re +
                                         4.0 * sig_right + nodes_right[kink_index + 1]) / 6.0
            integral[kink_index] = left_integral + right_integral
        else:
            kink_index = -1
            kink_fraction = 0.0
    F_nodes = np.concatenate(([0.0], np.cumsum(integral)))
    if kink_index >= 0 and left_integral is not None:
        phi_re = 1.5 * (F_nodes[kink_index] + left_integral) - n_re + Nv[0]
    else:
        phi_re = 0.0
    half_integral = h_arr * (nodes_left[:-1] + 4.0 * quarter_sigma + mid_sigma) / 12.0
    F_mid = F_nodes[:-1] + half_integral
    N0 = Nv[0]
    Phi_grid = 1.5 * F_nodes - Nv + N0
    Phi_mid = 1.5 * F_mid - mid + N0
    Psi = 3.0 * F_nodes - 4.0 * Nv
    S2, S2inv, relative_used, auto_overflow = _primitive_scales(Psi, relative_scale)
    if auto_overflow:
        replacement = _replace_overflowed_primitive(
            m, Nv, DN_eff, mid, auto_overflow=True,
            relative_scale=relative_scale,
        )
        if replacement is not None:
            Phi_grid, Phi_mid, S2, S2inv, phi_re, relative_used = replacement
    m._fast_relative_background_required = bool(relative_used)
    return (Phi_grid.astype(np.float64), Phi_mid.astype(np.float64),
            S2.astype(np.float64), S2inv.astype(np.float64),
            kink_index, float(kink_fraction), float(phi_re))


@njit(cache=True)
def _fast_phi_s2_fill(Nv, nodes_left, nodes_right, kink_index,
                      kink_fraction, sig_left, sig_re, sig_right,
                      h_arr, mid_sigma, integral, F_nodes, F_mid,
                      Phi_grid, Phi_mid, S2, S2inv, relative_scale=False):
    n = len(Nv)
    steep_panel = False
    max_derivative = 0.0
    for i in range(n - 1):
        h_arr[i] = Nv[i + 1] - Nv[i]
        derivative = abs(1.5 * nodes_left[i] - 1.0)
        if derivative > max_derivative:
            max_derivative = derivative
        derivative = abs(1.5 * nodes_right[i + 1] - 1.0)
        if derivative > max_derivative:
            max_derivative = derivative
        if (abs(nodes_right[i + 1] - nodes_left[i]) * h_arr[i]
                > _H2_IDENTITY_PANEL_VARIATION_LIMIT):
            steep_panel = True
        mid_sigma[i] = 0.5 * (nodes_left[i] + nodes_right[i + 1])
        integral[i] = h_arr[i] * (
            nodes_left[i] + 4.0 * mid_sigma[i] + nodes_right[i + 1]) / 6.0
    left_integral = 0.0
    if 0 <= kink_index < n - 1 and 0.0 < kink_fraction < 1.0:
        left_h = h_arr[kink_index] * kink_fraction
        right_h = h_arr[kink_index] - left_h
        left_integral = left_h * (
            nodes_left[kink_index] + 4.0 * sig_left + 1.0) / 6.0
        right_integral = right_h * (
            sig_re + 4.0 * sig_right + nodes_right[kink_index + 1]) / 6.0
        integral[kink_index] = left_integral + right_integral
    F_nodes[0] = 0.0
    compensation = 0.0
    outside_float_range = False
    N0 = Nv[0]
    derivative = abs(1.5 * nodes_left[0] - 1.0)
    if derivative > max_derivative:
        max_derivative = derivative
    derivative = abs(1.5 * nodes_right[0] - 1.0)
    if derivative > max_derivative:
        max_derivative = derivative
    derivative = abs(1.5 * nodes_left[n - 1] - 1.0)
    if derivative > max_derivative:
        max_derivative = derivative
    psi0 = 3.0 * F_nodes[0] - 4.0 * Nv[0]
    if psi0 < _PSI_MIN_FLOAT_RANGE or psi0 > _PSI_MAX_FLOAT_RANGE:
        outside_float_range = True
    Phi_grid[0] = 1.5 * F_nodes[0] - Nv[0] + N0
    monotone_phi = True
    max_dphi = 0.0
    S2[0] = 1.0 if relative_scale else np.exp(psi0)
    S2inv[0] = 1.0 if relative_scale else np.exp(-0.5 * psi0)
    for i in range(n - 1):
        if relative_scale:
            increment = integral[i] - compensation
            F_nodes[i + 1] = F_nodes[i] + increment
            compensation = (F_nodes[i + 1] - F_nodes[i]) - increment
        else:
            F_nodes[i + 1] = F_nodes[i] + integral[i]
        psi = 3.0 * F_nodes[i + 1] - 4.0 * Nv[i + 1]
        if psi < _PSI_MIN_FLOAT_RANGE or psi > _PSI_MAX_FLOAT_RANGE:
            outside_float_range = True
        Phi_grid[i + 1] = 1.5 * F_nodes[i + 1] - Nv[i + 1] + N0
        delta_phi = Phi_grid[i + 1] - Phi_grid[i]
        if delta_phi < 0.0:
            monotone_phi = False
        if delta_phi > max_dphi:
            max_dphi = delta_phi
        S2[i + 1] = 1.0 if relative_scale else np.exp(psi)
        S2inv[i + 1] = 1.0 if relative_scale else np.exp(-0.5 * psi)
        F_mid[i] = F_nodes[i] + h_arr[i] * (
            nodes_left[i] + 4.0 * (0.75 * nodes_left[i]
                                    + 0.25 * nodes_right[i + 1])
            + mid_sigma[i]) / 12.0
        Phi_mid[i] = 1.5 * F_mid[i] - 0.5 * (Nv[i] + Nv[i + 1]) + N0
    relative_used = relative_scale or outside_float_range
    auto_overflow = outside_float_range and not relative_scale
    identity_required = auto_overflow or steep_panel
    if auto_overflow:
        for i in range(n):
            S2[i] = 1.0
            S2inv[i] = 1.0
    has_panel_defect = False
    for i in range(n - 1):
        if identity_required and not (
            i == kink_index and 0.0 < kink_fraction < 1.0
        ):
            defect = Phi_grid[i + 1] - 2.0 * Phi_mid[i] + Phi_grid[i]
            if abs(defect) > 1.0e-5:
                has_panel_defect = True
    return (
        left_integral, relative_used, identity_required, has_panel_defect,
        monotone_phi, max_dphi, max_derivative,
    )


def _fast_phi_s2_workspace(m, n):
    workspace = getattr(m, '_fast_phi_s2_workspace', None)
    if workspace is None or workspace['n'] != n:
        workspace = {'n': n}
        workspace['h_arr'] = np.empty(n - 1, dtype=np.float64)
        workspace['mid_sigma'] = np.empty(n - 1, dtype=np.float64)
        workspace['integral'] = np.empty(n - 1, dtype=np.float64)
        workspace['F_nodes'] = np.empty(n, dtype=np.float64)
        workspace['F_mid'] = np.empty(n - 1, dtype=np.float64)
        workspace['Phi_grid'] = np.empty(n, dtype=np.float64)
        workspace['Phi_mid'] = np.empty(n - 1, dtype=np.float64)
        workspace['S2'] = np.empty(n, dtype=np.float64)
        workspace['S2inv'] = np.empty(n, dtype=np.float64)
        m._fast_phi_s2_workspace = workspace
    return workspace


def fast_phi_s2_split(m, Nv, DN_eff, sigma_nodes=None, *, relative_scale=False):
    """Build a kink-exact primitive with a cheap smooth-background estimate.

    The formal fast path already has continuous background values at every
    expansion node.  Away from the reheating breakpoint, the midpoint and
    quarter-point values are obtained by linear interpolation of those node
    values instead of evaluating the full continuous-sigma spline twice.  The
    interval containing ``N_re`` is still integrated on both sides with exact
    one-sided probes, so no transition step is allowed to cross the kink.
    Returned arrays are model-local work buffers and are overwritten by the
    next call for the same model.  The fast outer iteration consumes each
    primitive before requesting the next one, so this lifetime is intentional.

    中文：构造精确跨越再加热转折的原函数。转折以外的中点和四分点通过相邻连续背景节点线性插值得到，
    避免重复查询完整 ``sigma`` 样条；跨转折区间仍用单侧精确采样分段积分。返回数组是模型局部工作缓冲区，
    同一模型下次调用会覆盖；fast 外层迭代在下一次申请前已消费当前结果。
    """
    Nv = np.asarray(Nv, dtype=np.float64)
    if Nv.size < 2:
        m._fast_wkb_global_summary = None
        return (np.zeros_like(Nv), np.zeros_like(Nv), np.ones_like(Nv),
                np.ones_like(Nv), -1, 0.0, 0.0)
    nodes_left, nodes_right = _sigma_node_limits(Nv, m, DN_eff, sigma_nodes)

    d = m.derived_param
    n_re = float(d['N_inf'] - d['N_re'])
    kink_start = int(np.searchsorted(Nv, n_re - 1.0e-12, side='left'))
    kink_stop = int(np.searchsorted(Nv, n_re + 1.0e-12, side='right'))
    if kink_stop > kink_start:
        # The interval ending at N_re needs the pre-kink left limit; the
        # interval starting there needs the post-kink right limit. Treat a
        # roundoff-close node cluster as the same breakpoint.
        # N_re前一段用左极限、后一段用右极限；将舍入误差内的节点簇视为同一断点。
        sigma_pre = nodes_left[kink_start]
        sigma_post = nodes_right[kink_stop - 1]
        nodes_right[kink_start:kink_stop] = sigma_pre
        nodes_left[kink_start:kink_stop] = sigma_post
    kink_index = int(np.searchsorted(Nv, n_re, side='right') - 1)
    kink_fraction = 0.0
    sig_left = sig_re = sig_right = 0.0
    if 0 <= kink_index < Nv.size - 1:
        left = float(Nv[kink_index])
        right = float(Nv[kink_index + 1])
        kink_fraction = (n_re - left) / (right - left)
        if 0.0 < kink_fraction < 1.0:
            probes = np.array([
                left + 0.5 * (n_re - left),
                n_re,
                n_re + 0.5 * (right - n_re),
            ], dtype=float)
            sig_left, sig_re, sig_right = sigma_vec(probes, m, DN_eff)
        else:
            kink_index = -1
            kink_fraction = 0.0

    workspace = _fast_phi_s2_workspace(m, Nv.size)
    (
        left_integral, relative_used, identity_required, has_panel_defect,
        monotone_phi, max_dphi, max_derivative,
    ) = _fast_phi_s2_fill(
        Nv, nodes_left, nodes_right, kink_index, kink_fraction,
        sig_left, sig_re, sig_right, workspace['h_arr'],
        workspace['mid_sigma'], workspace['integral'], workspace['F_nodes'],
        workspace['F_mid'], workspace['Phi_grid'], workspace['Phi_mid'],
        workspace['S2'], workspace['S2inv'], relative_scale)
    m._fast_wkb_global_summary = None
    if identity_required:
        replacement = _replace_overflowed_primitive(
            m, Nv, DN_eff, 0.5 * (Nv[:-1] + Nv[1:]),
            auto_overflow=True, relative_scale=relative_scale,
        )
        if replacement is not None:
            phi, phi_mid, s2, s2inv, phi_re, relative_used = replacement
            workspace['Phi_grid'][:] = phi
            workspace['Phi_mid'][:phi_mid.size] = phi_mid
            workspace['S2'][:] = s2
            workspace['S2inv'][:] = s2inv
        else:
            phi_re = 0.0
    elif kink_index >= 0 and 0.0 < kink_fraction < 1.0:
        phi_re = (1.5 * (workspace['F_nodes'][kink_index] + left_integral)
                  - n_re + Nv[0])
    else:
        phi_re = 0.0
    m._fast_relative_background_required = bool(relative_used)
    m._fast_panel_curvature_defect = bool(has_panel_defect)
    if not identity_required:
        m._fast_wkb_global_summary = (
            bool(monotone_phi), float(max_dphi), float(max_derivative)
        )
    return (workspace['Phi_grid'], workspace['Phi_mid'],
            workspace['S2'], workspace['S2inv'],
            kink_index, float(kink_fraction), float(phi_re))


def fast_phi_s2_quarters(
    m, Nv, DN_eff, *, phi_grid, kink_index=-1, kink_fraction=0.0,
    sigma_nodes=None,
):
    """Evaluate local ``Phi`` at quarter boundaries on smooth panels.

    The interval integral is assembled from three composite-Simpson quarter
    panels. A panel containing the physical reheating kink is left unused by
    the tensor kernel and is deliberately not sampled across that transition.

    中文：在光滑背景区间上用三个复合 Simpson 子区间计算原函数四分点值。包含真实再加热转折的区间
    由张量核单独精确拆分，因此此函数不跨越该不连续点采样。
    """
    nodes = np.asarray(Nv, dtype=np.float64)
    phi = np.asarray(phi_grid, dtype=np.float64)
    if nodes.ndim != 1 or phi.shape != nodes.shape:
        raise ValueError("Nv and phi_grid must be equal one-dimensional arrays")
    if nodes.size < 2:
        return np.empty((2, 0), dtype=np.float64)
    h = np.diff(nodes)
    if np.any(h <= 0.0) or not np.isfinite(h).all():
        raise ValueError("Nv must be finite and strictly increasing")

    smooth = np.ones(h.size, dtype=bool)
    if 0 <= int(kink_index) < h.size and 0.0 < float(kink_fraction) < 1.0:
        smooth[int(kink_index)] = False
    result = np.zeros((2, h.size), dtype=np.float64)
    active = np.flatnonzero(smooth)
    if active.size == 0:
        return result

    fractions = np.array([0.125, 0.25, 0.375, 0.5, 0.625, 0.75])
    probes = nodes[active, None] + h[active, None] * fractions[None, :]
    sampled_sigma = sigma_vec(probes.reshape(-1), m, DN_eff).reshape(
        active.size, fractions.size
    )
    left, right = _sigma_node_limits(nodes, m, DN_eff, sigma_nodes)
    f_nodes = (phi + nodes - nodes[0]) / 1.5
    width = h[active]
    integral_0_25 = width / 24.0 * (
        left[active] + 4.0 * sampled_sigma[:, 0] + sampled_sigma[:, 1]
    )
    integral_25_50 = width / 24.0 * (
        sampled_sigma[:, 1] + 4.0 * sampled_sigma[:, 2] + sampled_sigma[:, 3]
    )
    integral_50_75 = width / 24.0 * (
        sampled_sigma[:, 3] + 4.0 * sampled_sigma[:, 4] + sampled_sigma[:, 5]
    )
    f_quarter = f_nodes[active] + integral_0_25
    f_three_quarter = f_quarter + integral_25_50 + integral_50_75
    n0 = nodes[0]
    result[0, active] = (
        1.5 * f_quarter - (nodes[active] + 0.25 * width) + n0
    )
    result[1, active] = (
        1.5 * f_three_quarter - (nodes[active] + 0.75 * width) + n0
    )
    return result
