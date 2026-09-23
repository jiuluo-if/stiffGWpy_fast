"""TDD contract for the Round 51 uniform-integral block screen."""

import math


def test_log_mean_exp_is_symmetric_and_exact_at_zero_delta():
    from scripts.benchmark_uniform_asymptotic_defect_round51 import _log_mean_exp

    assert _log_mean_exp(0.4, 0.4) == math.exp(0.4)
    assert _log_mean_exp(-0.2, 0.7) == _log_mean_exp(0.7, -0.2)


def test_uniform_block_defect_is_zero_for_constant_z():
    from scripts.benchmark_uniform_asymptotic_defect_round51 import _uniform_block_defect

    defect = _uniform_block_defect(0.6, 0.6, 0.01)
    assert defect == 0.0


def test_uniform_block_map_matches_production_constant_z_step():
    import numpy as np

    from scripts.benchmark_uniform_asymptotic_defect_round51 import _uniform_block_map
    from stiffgwpy_fast.fast_sgwb import scaled_step

    x0, y0 = 0.37, -0.81
    x1, y1 = _uniform_block_map(x0, y0, 0.6, 0.6, 0.02)
    x_ref, y_ref = scaled_step(x0, y0, 0.6, 0.02)
    assert math.isfinite(x1) and math.isfinite(y1)
    assert np.array_equal(np.asarray([x1, y1]), np.asarray([x_ref, y_ref]))
