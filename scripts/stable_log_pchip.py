"""Stable positive-data PCHIP arithmetic from logarithmic node values."""

from __future__ import annotations

import math

import numpy as np


def _signed_log_add(sign_a: int, log_a: float, sign_b: int, log_b: float) -> tuple[int, float]:
    if sign_a == 0 or log_a == -math.inf:
        return sign_b, log_b
    if sign_b == 0 or log_b == -math.inf:
        return sign_a, log_a
    if sign_a == sign_b:
        hi, lo = (log_a, log_b) if log_a >= log_b else (log_b, log_a)
        return sign_a, hi + math.log1p(math.exp(lo - hi))
    if log_a == log_b:
        return 0, -math.inf
    if log_a > log_b:
        return sign_a, log_a + math.log(-math.expm1(log_b - log_a))
    return sign_b, log_b + math.log(-math.expm1(log_a - log_b))


def _signed_log_sum(terms: list[tuple[int, float]]) -> tuple[int, float]:
    sign = 0
    logabs = -math.inf
    for term_sign, term_log in terms:
        sign, logabs = _signed_log_add(sign, logabs, term_sign, term_log)
    return sign, logabs


def _log_secants(x: np.ndarray, log_y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = np.diff(x)
    delta = np.diff(log_y)
    signs = np.sign(delta).astype(np.int8)
    logs = np.full(delta.shape, -math.inf, dtype=np.float64)
    for i, diff in enumerate(delta):
        if diff > 0.0:
            logs[i] = log_y[i + 1] + math.log(-math.expm1(-diff)) - math.log(h[i])
        elif diff < 0.0:
            logs[i] = log_y[i] + math.log(-math.expm1(diff)) - math.log(h[i])
    return signs, logs


def _endpoint_derivative(
    h0: float,
    h1: float,
    sign0: int,
    log0: float,
    sign1: int,
    log1: float,
) -> tuple[int, float]:
    den = h0 + h1
    coef0 = (2.0 * h0 + h1) / den
    coef1 = -h0 / den
    sign, logabs = _signed_log_sum([
        (sign0, math.log(coef0) + log0 if sign0 else -math.inf),
        (sign1 if coef1 > 0.0 else -sign1,
         math.log(abs(coef1)) + log1 if sign1 else -math.inf),
    ])
    if sign != sign0:
        return 0, -math.inf
    if sign0 != sign1 and logabs > math.log(3.0) + log0:
        return sign0, math.log(3.0) + log0
    return sign, logabs


def _derivatives(x: np.ndarray, log_y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = x.size
    h = np.diff(x)
    secant_sign, secant_log = _log_secants(x, log_y)
    derivative_sign = np.zeros(n, dtype=np.int8)
    derivative_log = np.full(n, -math.inf, dtype=np.float64)
    if n == 2:
        derivative_sign[:] = secant_sign[0]
        derivative_log[:] = secant_log[0]
        return derivative_sign, derivative_log

    derivative_sign[0], derivative_log[0] = _endpoint_derivative(
        h[0], h[1], secant_sign[0], secant_log[0], secant_sign[1], secant_log[1]
    )
    derivative_sign[-1], derivative_log[-1] = _endpoint_derivative(
        h[-1], h[-2], secant_sign[-1], secant_log[-1], secant_sign[-2], secant_log[-2]
    )
    for k in range(1, n - 1):
        sign0, log0 = int(secant_sign[k - 1]), float(secant_log[k - 1])
        sign1, log1 = int(secant_sign[k]), float(secant_log[k])
        if sign0 == 0 or sign1 == 0 or sign0 != sign1:
            continue
        h0, h1 = float(h[k - 1]), float(h[k])
        w1 = 2.0 * h1 + h0
        w2 = h1 + 2.0 * h0
        log_den = float(np.logaddexp(math.log(w1) - log0, math.log(w2) - log1))
        derivative_sign[k] = sign0
        derivative_log[k] = math.log(w1 + w2) - log_den
    return derivative_sign, derivative_log


def _basis(t: float) -> tuple[float, float, float, float]:
    t2 = t * t
    t3 = t2 * t
    return 2.0 * t3 - 3.0 * t2 + 1.0, t3 - 2.0 * t2 + t, -2.0 * t3 + 3.0 * t2, t3 - t2


class LogPchip:
    """Evaluate/integrate the same positive linear PCHIP from log-y nodes: 用 log-y 稳定计算原线性 PCHIP。"""

    def __init__(self, x: np.ndarray, log_y: np.ndarray):
        self.x = np.asarray(x, dtype=np.float64)
        self.log_y = np.asarray(log_y, dtype=np.float64)
        if self.x.ndim != 1 or self.log_y.ndim != 1 or self.x.size != self.log_y.size:
            raise ValueError("x and log_y must be equally sized one-dimensional arrays")
        if self.x.size < 2 or not np.isfinite(self.x).all() or not np.isfinite(self.log_y).all():
            raise ValueError("PCHIP needs at least two finite x and log_y values")
        if np.any(np.diff(self.x) <= 0.0):
            raise ValueError("x must be strictly increasing")
        self.derivative_sign, self.derivative_log = _derivatives(self.x, self.log_y)

    def evaluate_log(self, points: np.ndarray) -> np.ndarray:
        points_arr = np.asarray(points, dtype=np.float64)
        flat = points_arr.ravel()
        out = np.empty_like(flat)
        indices = np.searchsorted(self.x, flat, side="right") - 1
        indices = np.clip(indices, 0, self.x.size - 2)
        for q, (point, i) in enumerate(zip(flat, indices)):
            h = float(self.x[i + 1] - self.x[i])
            t = (float(point) - float(self.x[i])) / h
            h00, h10, h01, h11 = _basis(t)
            terms = []
            if h00 != 0.0:
                terms.append((1, self.log_y[i] + math.log(abs(h00))) if h00 > 0.0 else (-1, self.log_y[i] + math.log(abs(h00))))
            if h01 != 0.0:
                terms.append((1, self.log_y[i + 1] + math.log(abs(h01))) if h01 > 0.0 else (-1, self.log_y[i + 1] + math.log(abs(h01))))
            if h10 != 0.0 and self.derivative_sign[i] != 0:
                sign = int(self.derivative_sign[i]) * (1 if h10 > 0.0 else -1)
                terms.append((sign, math.log(h) + self.derivative_log[i] + math.log(abs(h10))))
            if h11 != 0.0 and self.derivative_sign[i + 1] != 0:
                sign = int(self.derivative_sign[i + 1]) * (1 if h11 > 0.0 else -1)
                terms.append((sign, math.log(h) + self.derivative_log[i + 1] + math.log(abs(h11))))
            sign, logabs = _signed_log_sum(terms)
            if sign <= 0 or not math.isfinite(logabs):
                raise FloatingPointError("positive PCHIP evaluation lost positivity")
            out[q] = logabs
        result = out.reshape(points_arr.shape)
        return result

    def integrate_log(self, a: float | None = None, b: float | None = None) -> float:
        left = float(self.x[0] if a is None else a)
        right = float(self.x[-1] if b is None else b)
        if left != float(self.x[0]) or right != float(self.x[-1]):
            raise ValueError("integrate_log currently integrates the full knot interval")
        h = np.diff(self.x)
        areas = []
        for i, width in enumerate(h):
            log_width = math.log(float(width))
            log_h2 = 2.0 * log_width
            terms = [
                (1, log_width - math.log(2.0) + float(self.log_y[i])),
                (1, log_width - math.log(2.0) + float(self.log_y[i + 1])),
            ]
            if self.derivative_sign[i] != 0:
                terms.append((int(self.derivative_sign[i]), log_h2 - math.log(12.0) + float(self.derivative_log[i])))
            if self.derivative_sign[i + 1] != 0:
                terms.append((-int(self.derivative_sign[i + 1]), log_h2 - math.log(12.0) + float(self.derivative_log[i + 1])))
            sign, log_area = _signed_log_sum(terms)
            if sign <= 0 or not math.isfinite(log_area):
                raise FloatingPointError("positive PCHIP panel integral lost positivity")
            areas.append(log_area)
        return float(np.logaddexp.reduce(np.asarray(areas, dtype=np.float64)))
