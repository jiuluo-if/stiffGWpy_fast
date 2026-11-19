"""生产张量传递核的局部自适应 WKB3 契约。"""

import numpy as np

from scripts.benchmark_bessel_transfer_round57 import _fundamental
from stiffgwpy_fast.fast_sgwb import _phase_segment_base, _phase_segment_wkb3


def _exact_bessel_map(z_start: float, z_end: float, h_step: float) -> np.ndarray:
    q = (z_end - z_start) / h_step
    left = _fundamental(z_start, q, 0.0)
    right = _fundamental(z_end, q, h_step)
    return right @ np.linalg.solve(left, np.eye(2))


def test_fast_wkb3_replaces_resolved_high_phase_segment():
    z_start = 14.765730539835435
    z_end = 18.896093663906587
    h_step = 0.00125

    x, y, used_wkb, n_sub, defect = _phase_segment_wkb3(
        0.0,
        1.0,
        z_start,
        z_end,
        h_step,
        phase_max=0.25,
        error_tolerance=1.0e-5,
        min_substeps=64,
    )

    assert used_wkb
    assert n_sub == 101987
    assert defect < 1.0e-5
    expected = _exact_bessel_map(z_start, z_end, h_step) @ np.array([0.0, 1.0])
    np.testing.assert_allclose((x, y), expected, rtol=0.0, atol=5.0e-7)


def test_fast_wkb3_defect_estimate_is_state_independent():
    args = (
        14.765730539835435,
        18.896093663906587,
        0.00125,
        1.0e-5,
    )
    defect_x = _phase_segment_wkb3(
        1.0, 0.0, *args[:3], 0.25, args[3], 64
    )[4]
    defect_y = _phase_segment_wkb3(
        0.0, 1.0, *args[:3], 0.25, args[3], 64
    )[4]

    assert defect_x == defect_y


def test_fast_wkb3_uses_exact_phase_fallback_when_defect_is_too_large():
    state = (0.0, 1.0)
    z_start = 14.765730539835435
    z_end = 18.896093663906587
    h_step = 0.00125

    actual = _phase_segment_wkb3(
        *state,
        z_start,
        z_end,
        h_step,
        phase_max=0.25,
        error_tolerance=1.0e-12,
        min_substeps=64,
    )
    expected = _phase_segment_base(*state, z_start, z_end, h_step, 0.25)

    assert not actual[2]
    np.testing.assert_array_equal(actual[:2], expected)


def test_fast_wkb3_preserves_turning_point_fallback():
    state = (0.0, 1.0)
    z_start, z_end, h_step = 0.4, 0.4025, 0.005
    actual = _phase_segment_wkb3(
        *state,
        z_start,
        z_end,
        h_step,
        phase_max=0.25,
        error_tolerance=1.0e-3,
        min_substeps=0,
    )
    expected = _phase_segment_base(*state, z_start, z_end, h_step, 0.25)

    assert not actual[2]
    np.testing.assert_array_equal(actual[:2], expected)
