# -*- coding: utf-8 -*-
"""Numeric comparison helpers shared by the validation scripts and tests.

Error metrics follow the independent audit recommendations:

* linear-Omega relative error (``rel_linear_omega``) instead of a relative
  error computed on log10(Omega) itself;
* absolute error in dex (``dex_abs``);
* plain relative error with a safe denominator (``rel_abs``), and a
  "signal region" mask so near-zero channels are reported separately.
"""

import numpy as np

__all__ = ['dex_abs', 'rel_abs', 'rel_linear_omega', 'signal_mask']


def dex_abs(a_log10, b_log10):
    """Absolute difference between two log10 quantities (dex units)."""
    a = np.asarray(a_log10, dtype=float)
    b = np.asarray(b_log10, dtype=float)
    return np.abs(a - b)


def rel_abs(a, b):
    """Elementwise relative difference ``|a - b| / |b|``, guarding ``b == 0``."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-300)


def rel_linear_omega(a_log10, b_log10):
    """Relative difference of *linear* Omega values: ``|10^a - 10^b| / 10^b``."""
    a = np.power(10.0, np.asarray(a_log10, dtype=float))
    b = np.power(10.0, np.asarray(b_log10, dtype=float))
    return np.abs(a - b) / np.maximum(b, 1e-300)


def signal_mask(b_log10, floor_log10=-30.0):
    """Mask of channels whose reference linear Omega is above ``floor_log10``."""
    b = np.asarray(b_log10, dtype=float)
    return b >= floor_log10
