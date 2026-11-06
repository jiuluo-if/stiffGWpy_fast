"""Standalone frequency tangent for Fast's existing discrete phase transfer.

This prototype differentiates the implemented transfer map on a fixed
frequency branch. It does not differentiate integer start-node, phase-substep,
or WKB-handoff changes and is not a production error bound.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from stiffgwpy_fast import fast_sgwb as FS


@njit(cache=True)
def propagate_mode_with_tangent(
    Nv, Phi_grid, Phi_mid, S2, S2inv, Sv, h_arr,
    j0, z0, P_t, fp_freq, ev_minus, fp_minus,
    phase_max, z_tail, kink_index, kink_fraction, phi_re, n_t,
):
    """Propagate one fixed-start mode and its log-frequency tangent.

    The path mirrors the production midpoint/kink transfer and terminal
    analytic tail for one mode. It intentionally omits history assembly and
    the optional pre-event Bessel shortcut. Integer branch changes in `j0`,
    phase subdivision, and tail classification are outside its derivative.
    """
    terminal = len(Nv) - 1
    phi0 = Phi_grid[j0]
    x = 0.0
    y = math.exp(z0) * S2inv[j0]
    vx = 0.0
    vy = y
    k = j0
    zz = z0
    last_x, last_y, last_vx, last_vy, last_z = x, y, vx, vy, zz

    z_tail_mode = z_tail
    phase_max_mode = phase_max
    deferred = False
    if Sv is not None:
        endpoint_z = z0 + Phi_grid[terminal] - phi0
        endpoint_defect = FS._endpoint_wkb_defect(Sv[terminal], endpoint_z)
        if endpoint_defect > FS._WKB_TAIL_ENDPOINT_EPS_LIMIT:
            z_tail_mode = math.inf
            deferred = True
            if 0.25 <= phase_max <= FS._WKB_TAIL_DEFERRED_PHASE_MAX:
                phase_max_mode = FS._WKB_TAIL_DEFERRED_PHASE_MAX

    while k < terminal and zz < z_tail_mode:
        h_step = h_arr[k]
        z_node = zz
        if k == kink_index and 0.0 < kink_fraction < 1.0:
            z_break = z0 + phi_re - phi0
            z_end = z0 + Phi_grid[k + 1] - phi0
            h_left = h_step * kink_fraction
            h_right = h_step - h_left
            (x, y), (vx, vy) = phase_segment_tangent(
                x, y, vx, vy, z_node, z_break, h_left, phase_max_mode
            )
            (x, y), (vx, vy) = phase_segment_tangent(
                x, y, vx, vy, z_break, z_end, h_right, phase_max_mode
            )
        else:
            z_mid = z0 + Phi_mid[k] - phi0
            z_end = 2.0 * z_mid - z_node
            (x, y), (vx, vy) = phase_segment_tangent(
                x, y, vx, vy, z_node, z_end, h_step, phase_max_mode
            )
        k += 1
        zz = z0 + Phi_grid[k] - phi0
        if zz < z_tail_mode:
            last_x, last_y, last_vx, last_vy, last_z = x, y, vx, vy, zz

    if zz < z_tail_mode:
        norm2 = x * x + y * y
        spectrum = S2[terminal] * norm2 * P_t / 24.0
        slope = terminal_log_spectrum_slope(x, y, vx, vy, n_t)
        return spectrum, slope, x, y, vx, vy, deferred

    kend = k - 1
    if kend < j0:
        kend = j0
    x, y, vx, vy, last_z = last_x, last_y, last_vx, last_vy, last_z
    e_z = math.exp(-last_z)
    if Sv is not None:
        gamma = FS._tail_match_gamma(Sv[kend])
        y_factor = 1.0 + gamma * gamma * e_z * e_z
        cross = 2.0 * gamma * e_z
        amp2 = x * x + y * y * y_factor + cross * x * y
        amp2_u = (
            2.0 * x * vx + 2.0 * y * vy * y_factor
            - 2.0 * gamma * gamma * e_z * e_z * y * y
            - cross * x * y + cross * (vx * y + x * vy)
        )
    else:
        amp2 = x * x + y * y
        amp2_u = 2.0 * (x * vx + y * vy)
    if amp2 <= 0.0 or not math.isfinite(amp2):
        return math.nan, math.nan, x, y, vx, vy, deferred
    coeff = math.sqrt(0.5 * S2[kend] * amp2)
    coeff_log_slope = 0.5 * amp2_u / amp2
    e_nz = math.exp(Nv[kend] - last_z)
    th = coeff * e_nz * ev_minus[terminal]
    xf = th * fp_freq * fp_minus[terminal]
    op = xf * xf / 36.0 * P_t
    spectrum = 3.0 * op
    # d ln(eNz)/d ln(f) = -1 cancels d ln(fp_freq)/d ln(f) = +1.
    slope = n_t + 2.0 * coeff_log_slope
    return spectrum, slope, x, y, vx, vy, deferred


@njit(parallel=True, cache=True)
def propagate_frequency_grid(
    Nv, Phi_grid, Phi_mid, S2, S2inv, Sv, h_arr,
    j0s, z0s, P_t, fp_freq, ev_minus, fp_minus,
    phase_max, z_tail, kink_index, kink_fraction, phi_re, n_t,
):
    """Run independent fixed-branch mode tangents in a parallel Numba loop."""
    spectra = np.empty(j0s.size, dtype=np.float64)
    slopes = np.empty(j0s.size, dtype=np.float64)
    for mode in prange(j0s.size):
        spectrum, slope, _, _, _, _, _ = propagate_mode_with_tangent(
            Nv, Phi_grid, Phi_mid, S2, S2inv, Sv, h_arr,
            j0s[mode], z0s[mode], P_t[mode], fp_freq[mode],
            ev_minus, fp_minus, phase_max, z_tail, kink_index,
            kink_fraction, phi_re, n_t[mode],
        )
        spectra[mode] = spectrum
        slopes[mode] = slope
    return spectra, slopes


@njit(cache=True, inline="always")
def scaled_step_tangent(x, y, vx, vy, z_mid, h):
    """Apply `scaled_step` and its exact z-shift derivative to state/tangent."""
    w = math.exp(z_mid)
    w2 = w * w
    q = w2 - 1.0
    if q > 0.0:
        omega = math.sqrt(q)
        theta = omega * h
        c = math.cos(theta)
        s = math.sin(theta) / omega
        if abs(theta) < 0.25:
            t2 = theta * theta
            sinc = 1.0 - t2 / 6.0 + t2 * t2 / 120.0 - t2 * t2 * t2 / 5040.0 + t2 * t2 * t2 * t2 / 362880.0
            c_z = -w2 * h * h * sinc
            s_z = w2 * h * h * h * (
                -1.0 / 3.0 + t2 / 30.0 - t2 * t2 / 840.0
                + t2 * t2 * t2 / 45360.0 - t2 * t2 * t2 * t2 / 3991680.0
            )
        else:
            c_z = -w2 * h * h * math.sin(theta) / theta
            s_z = w2 * h * h * h * (
                theta * math.cos(theta) - math.sin(theta)
            ) / (theta * theta * theta)
    else:
        c = 1.0 + 0.5 * (1.0 - w2) * h * h
        s = h * (1.0 + (1.0 - w2) * h * h / 6.0)
        c_z = -w2 * h * h
        s_z = -w2 * h * h * h / 3.0

    a = c - s
    d = c + s
    b = w * s
    a_z = c_z - s_z
    d_z = c_z + s_z
    b_z = w * (s + s_z)

    x_new = a * x - b * y
    y_new = b * x + d * y
    vx_new = a * vx - b * vy + a_z * x - b_z * y
    vy_new = b * vx + d * vy + b_z * x + d_z * y
    return (x_new, y_new), (vx_new, vy_new)


@njit(cache=True)
def phase_segment_tangent(x, y, vx, vy, z_start, z_end, h_step, phase_max):
    """Differentiate a fixed-subdivision `_phase_segment` path under z shift."""
    z_mid = 0.5 * (z_start + z_end)
    n_sub = FS._phase_substeps(h_step, z_mid, phase_max)
    if n_sub <= 1:
        return scaled_step_tangent(x, y, vx, vy, z_mid, h_step)

    h_sub = h_step / n_sub
    for i in range(n_sub):
        z_sub = z_start + (i + 0.5) * (z_end - z_start) / n_sub
        (x, y), (vx, vy) = scaled_step_tangent(x, y, vx, vy, z_sub, h_sub)
    return (x, y), (vx, vy)


@njit(cache=True, inline="always")
def terminal_log_spectrum_slope(x, y, vx, vy, n_t):
    """Return d ln[(x²+y²) f**n_t]/d ln(f) for nonzero terminal state."""
    norm2 = x * x + y * y
    if norm2 <= 0.0 or not math.isfinite(norm2):
        return math.nan
    return n_t + 2.0 * (x * vx + y * vy) / norm2


@njit(cache=True, inline="always")
def hermite_panel_midpoint(y0, y1, log_slope0, log_slope1, width):
    """Evaluate a cubic Hermite midpoint from endpoint log-spectrum slopes."""
    derivative0 = y0 * log_slope0
    derivative1 = y1 * log_slope1
    return 0.5 * (y0 + y1) + width * (derivative0 - derivative1) / 8.0


@njit(cache=True, inline="always")
def hermite_panel_integral(y0, y1, log_slope0, log_slope1, width):
    """Integrate a cubic Hermite panel from endpoint log-spectrum slopes."""
    derivative0 = y0 * log_slope0
    derivative1 = y1 * log_slope1
    return (
        0.5 * width * (y0 + y1)
        + width * width * (derivative0 - derivative1) / 12.0
    )


def initial_frequency_tangent(x: float, y: float) -> tuple[float, float]:
    """Derivative of the production start state for fixed `j0` and background."""
    return 0.0, float(y)
