# -*- coding: utf-8 -*-
import numpy as np

from stiffgwpy_fast._metrics import dex_abs, rel_abs, rel_linear_omega, signal_mask


def test_dex_abs():
    assert dex_abs(-10.0, -9.0) == 1.0
    a = np.array([-10.0, -5.0])
    b = np.array([-9.0, -5.0])
    np.testing.assert_allclose(dex_abs(a, b), [1.0, 0.0])


def test_rel_abs_guards_zero():
    assert rel_abs(1.0, 0.0) > 0.0
    assert rel_abs(0.0, 0.0) == 0.0


def test_rel_linear_omega_matches_dex_conversion():
    # A 4e-4 dex difference corresponds to |10^d - 1| in the linear value.
    dex = 4e-4
    rel = rel_linear_omega(np.array([-10.0]), np.array([-10.0 - dex]))
    assert abs(rel[0] - (10.0 ** dex - 1.0)) < 1e-12


def test_signal_mask():
    mask = signal_mask(np.array([-20.0, -40.0, -29.9]), floor_log10=-30.0)
    np.testing.assert_array_equal(mask, [True, False, True])
