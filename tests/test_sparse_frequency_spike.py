"""Contract tests for the standalone sparse-frequency propagation screen."""

import numpy as np

from scripts.benchmark_sparse_frequency_spike import (
    _mode_indices,
    _reconstruct,
)


def test_sparse_frequency_indices_retain_endpoints_and_reconstruction_is_exact_on_nodes():
    full = np.linspace(-4.0, 2.0, 11)
    indices = _mode_indices(len(full), stride=2)
    assert indices[0] == 0
    assert indices[-1] == len(full) - 1
    values = np.sin(full[indices])
    reconstructed = _reconstruct(full, full[indices], values)
    np.testing.assert_array_equal(reconstructed[indices], values)
