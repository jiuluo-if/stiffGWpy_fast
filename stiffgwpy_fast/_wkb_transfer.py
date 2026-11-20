"""局部正频率线性段的嵌入式三阶 WKB 张量传递映射。

中文：仅表达 ``z(N)`` 线性且变换后 ``Omega²`` 全段为正的光滑区间；调用方
根据二阶/三阶映射差估计局部误差，未通过时必须回退既有相位细分。
"""

from __future__ import annotations

import math

from numba import njit

_FLOAT64_MIN_NORMAL = 2.2250738585072014e-308


@njit(cache=True, inline="always")
def _phase_primitive(w: float, alpha: float) -> float:
    omega = math.sqrt(w * w - alpha)
    return omega - math.sqrt(alpha) * math.atan(omega / math.sqrt(alpha))


@njit(cache=True, inline="always")
def _phase2_primitive(omega: float, alpha: float) -> float:
    return -1.0 / omega - 5.0 * alpha / (3.0 * omega * omega * omega)


@njit(cache=True, inline="always")
def _s3_value(w: float, omega: float, q: float, alpha: float) -> float:
    omega2 = omega * omega
    return -q * q * w * w * (w * w + 4.0 * alpha) / (
        16.0 * omega2 * omega2 * omega2
    )


@njit(cache=True, inline="always")
def _s3_derivative(w: float, omega: float, q: float, alpha: float) -> float:
    value = _s3_value(w, omega, q, alpha)
    w2 = w * w
    return value * (2.0 * q * w2) * (
        1.0 / w2 + 1.0 / (w2 + 4.0 * alpha) - 3.0 / (w2 - alpha)
    )


@njit(cache=True, inline="always")
def _basis(
    z: float,
    n_value: float,
    q: float,
    alpha: float,
    phase: float,
    order: int,
) -> tuple[float, float, float, float]:
    w = math.exp(z)
    w2 = w * w
    omega2 = w2 - alpha
    omega = math.sqrt(omega2)
    omega_prime = q * w2 / omega
    theta_prime = omega
    amp_log_prime = 0.5 * q - 0.5 * omega_prime / omega
    amp_log = 0.5 * q * n_value - 0.5 * math.log(omega)
    if order >= 2:
        theta_prime += q * q * w2 * (w2 + 4.0 * alpha) / (
            8.0 * omega * omega2 * omega2
        )
    if order >= 3:
        amp_log += _s3_value(w, omega, q, alpha)
        amp_log_prime += _s3_derivative(w, omega, q, alpha)
    amplitude = math.exp(amp_log)
    c = math.cos(phase)
    s = math.sin(phase)
    u1 = amplitude * c
    u2 = amplitude * s
    up1 = amplitude * (amp_log_prime * c - theta_prime * s)
    up2 = amplitude * (amp_log_prime * s + theta_prime * c)
    y1 = -(up1 + u1) / w
    y2 = -(up2 + u2) / w
    return u1, u2, y1, y2


@njit(cache=True, inline="always")
def _transfer_order(
    x: float,
    y: float,
    z_start: float,
    z_end: float,
    h_step: float,
    order: int,
    alpha: float,
) -> tuple[float, float]:
    q = (z_end - z_start) / h_step
    w_start = math.exp(z_start)
    w_end = math.exp(z_end)
    omega_start = math.sqrt(w_start * w_start - alpha)
    omega_end = math.sqrt(w_end * w_end - alpha)
    phase = (_phase_primitive(w_end, alpha) - _phase_primitive(w_start, alpha)) / q
    if order >= 2:
        phase += q * (
            _phase2_primitive(omega_end, alpha)
            - _phase2_primitive(omega_start, alpha)
        ) / 8.0
    a11, a12, a21, a22 = _basis(z_start, 0.0, q, alpha, 0.0, order)
    b11, b12, b21, b22 = _basis(z_end, h_step, q, alpha, phase, order)
    diagonal_product = a11 * a22
    off_diagonal_product = a12 * a21
    determinant = diagonal_product - off_diagonal_product
    determinant_error = 4.0 * 2.220446049250313e-16 * (
        abs(diagonal_product) + abs(off_diagonal_product)
    )
    if (
        not math.isfinite(determinant)
        or not math.isfinite(determinant_error)
        or abs(determinant) < _FLOAT64_MIN_NORMAL
        or abs(determinant) <= determinant_error
    ):
        # The inverse is not resolved to float64 precision.  Signal an invalid
        # WKB candidate so the caller uses its existing phase-substep transfer.
        return math.nan, math.nan
    c1 = (a22 * x - a12 * y) / determinant
    c2 = (-a21 * x + a11 * y) / determinant
    return b11 * c1 + b12 * c2, b21 * c1 + b22 * c2


@njit(cache=True)
def wkb_transfer_embedded(
    x: float,
    y: float,
    z_start: float,
    z_end: float,
    h_step: float,
    error_tolerance: float,
) -> tuple[float, float, float, bool]:
    """Try WKB3 and estimate local defect from WKB2.

    The segment assumes linear ``z(N)``. False ``accepted`` requires the caller
    to use its established phase-substep transfer; no finite clamp is applied.
    """
    invalid = (math.nan, math.nan, math.inf, False)
    if (
        not math.isfinite(x)
        or not math.isfinite(y)
        or not math.isfinite(z_start)
        or not math.isfinite(z_end)
        or not math.isfinite(h_step)
        or not math.isfinite(error_tolerance)
        or h_step <= 0.0
        or error_tolerance <= 0.0
        or abs(z_start) > 350.0
        or abs(z_end) > 350.0
    ):
        return invalid
    q = (z_end - z_start) / h_step
    if not math.isfinite(q) or abs(q) < 1.0e-10:
        return invalid
    alpha = q + 1.0 + 0.25 * q * q
    if not math.isfinite(alpha) or alpha <= 0.0:
        return invalid
    w_start = math.exp(z_start)
    w_end = math.exp(z_end)
    if w_start * w_start <= alpha or w_end * w_end <= alpha:
        return invalid
    x2_1, y2_1 = _transfer_order(1.0, 0.0, z_start, z_end, h_step, 2, alpha)
    x2_2, y2_2 = _transfer_order(0.0, 1.0, z_start, z_end, h_step, 2, alpha)
    x3_1, y3_1 = _transfer_order(1.0, 0.0, z_start, z_end, h_step, 3, alpha)
    x3_2, y3_2 = _transfer_order(0.0, 1.0, z_start, z_end, h_step, 3, alpha)
    if not (
        math.isfinite(x2_1)
        and math.isfinite(y2_1)
        and math.isfinite(x2_2)
        and math.isfinite(y2_2)
        and math.isfinite(x3_1)
        and math.isfinite(y3_1)
        and math.isfinite(x3_2)
        and math.isfinite(y3_2)
    ):
        return invalid
    difference_norm = math.hypot(
        math.hypot(x3_1 - x2_1, y3_1 - y2_1),
        math.hypot(x3_2 - x2_2, y3_2 - y2_2),
    )
    transfer_norm = math.hypot(
        math.hypot(x3_1, y3_1), math.hypot(x3_2, y3_2)
    )
    if not math.isfinite(difference_norm) or not math.isfinite(transfer_norm):
        return invalid
    scale = max(transfer_norm, math.sqrt(2.0), 1.0e-300)
    defect = difference_norm / scale
    if not math.isfinite(defect):
        return invalid
    x3 = x3_1 * x + x3_2 * y
    y3 = y3_1 * x + y3_2 * y
    if not math.isfinite(x3) or not math.isfinite(y3):
        return invalid
    return x3, y3, defect, defect <= error_tolerance
