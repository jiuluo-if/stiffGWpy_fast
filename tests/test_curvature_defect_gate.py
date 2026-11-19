"""Tests for the local tensor-panel curvature dispatch gate."""

import numpy as np

from stiffgwpy_fast.fast_sgwb import _has_panel_curvature_defect


def test_panel_gate_skips_exactly_split_kink_interval():
    phi = np.array([0.0, 1.0, 1.5])
    phi_mid = np.array([0.5, 1.3])
    assert not _has_panel_curvature_defect(phi, phi_mid, 1, 0.5, 1e-5)
    assert _has_panel_curvature_defect(phi, phi_mid, -1, 0.0, 1e-5)


def test_panel_gate_keeps_flat_background_on_fast_path():
    phi = np.array([0.0, 0.5, 1.0, 1.5])
    phi_mid = np.array([0.25, 0.75, 1.25])
    assert not _has_panel_curvature_defect(phi, phi_mid, -1, 0.0, 1e-5)
