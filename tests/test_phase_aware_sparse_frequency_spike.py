"""Contract tests for phase-aware sparse-frequency reconstruction."""

import numpy as np

from scripts.benchmark_phase_aware_sparse_frequency_spike import (
    _complex_reconstruct,
    _mode_indices,
)


def test_complex_reconstruction_reproduces_sparse_cartesian_nodes():
    full = np.linspace(-4.0, 2.0, 9)
    indices = _mode_indices(len(full), stride=2)
    x = np.sin(full[indices])
    y = np.cos(full[indices])
    pref = np.ones(indices.size)
    reconstructed = _complex_reconstruct(
        full, full[indices], x, y, pref, np.ones(indices.size),
        np.ones(len(full)),
    )
    np.testing.assert_array_equal(reconstructed[indices], x * x + y * y)
