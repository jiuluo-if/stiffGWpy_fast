"""TDD contract for the standalone phase loop-state reuse spike."""

import numpy as np

from scripts.benchmark_phase_z_reuse import phase_state_reuse


def test_phase_state_reuse_matches_repeated_expression():
    z0 = 0.37
    phi0 = -0.2
    grid = np.array([-0.2, 0.1, 0.9], dtype=np.float64)

    expected = np.array([z0 + value - phi0 for value in grid])

    actual = phase_state_reuse(z0, phi0, grid)

    np.testing.assert_array_equal(actual, expected)
