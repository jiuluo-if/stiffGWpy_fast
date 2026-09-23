"""Contract test for the Round 67 transfer perturbation bound."""

import numpy as np

from scripts.benchmark_outer_certificate_round67 import (
    bound_transfer_perturbation,
)


def test_transfer_bound_dominates_direct_endpoint_difference():
    old_z = np.array((0.1, 0.2, 0.3), dtype=np.float64)
    new_z = old_z + 1.0e-5
    h_values = np.array((0.01, 0.01, 0.01), dtype=np.float64)
    old_state = np.array((0.0, 0.8), dtype=np.float64)
    new_state = np.array((0.0, 0.8 + 1.0e-6), dtype=np.float64)
    row = bound_transfer_perturbation(
        old_z, new_z, h_values, old_state, new_state)
    assert row["finite"], row
    assert row["endpoint_relative_actual"] <= (
        row["endpoint_relative_bound"] * (1.0 + 1.0e-12) + 1.0e-15), row
