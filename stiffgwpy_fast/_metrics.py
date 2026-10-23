# -*- coding: utf-8 -*-
"""Numeric comparison helpers shared by the validation scripts and tests.

Error metrics follow the independent audit recommendations:

* linear-Omega relative error (``rel_linear_omega``) instead of a relative
  error computed on log10(Omega) itself;
* absolute error in dex (``dex_abs``);
* plain relative error with a safe denominator (``rel_abs``), and a
  "signal region" mask so near-zero channels are reported separately.

中文：`Omega_GW` 的线性相对误差由对数谱差通过 ``expm1`` 稳定计算；dex 绝对差单独统计。
近零通道通过 ``signal_mask`` 分区，避免将近零谱值作为信号区相对误差的分母。
"""

import numpy as np

_LN10 = np.log(10.0)

__all__ = ['dex_abs', 'rel_abs', 'rel_linear_omega', 'signal_mask']


def dex_abs(a_log10, b_log10):
    """Absolute difference between two log10 quantities (dex units).

    中文：计算两个以 10 为底的对数值之差的绝对值，单位为 dex。
    """
    a = np.asarray(a_log10, dtype=float)
    b = np.asarray(b_log10, dtype=float)
    return np.abs(a - b)


def rel_abs(a, b):
    """Elementwise relative difference ``|a - b| / |b|``, guarding ``b == 0``.

    中文：逐元素计算 ``|a-b|/|b|``；当 ``b`` 为零时使用安全的正下限作为分母。
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-300)


def rel_linear_omega(a_log10, b_log10):
    """Relative difference of linear Omega, using ``10**b`` as denominator.

    The equivalent log-domain form avoids subtracting nearly equal linear values.

    中文：以候选值 ``10**b`` 为分母，在对数差上用 ``expm1`` 计算，避免线性谱相减时的消减误差。
    """
    a = np.asarray(a_log10, dtype=float)
    b = np.asarray(b_log10, dtype=float)
    return np.abs(np.expm1((a - b) * _LN10))


def signal_mask(b_log10, floor_log10=-30.0):
    """Mask of channels whose reference linear Omega is above ``floor_log10``.

    中文：标记参考谱 ``log10(Omega)`` 不低于 ``floor_log10`` 的通道，用于区分信号区和近零区。
    """
    b = np.asarray(b_log10, dtype=float)
    return b >= floor_log10
