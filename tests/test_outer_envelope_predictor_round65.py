"""TDD contract for the Round 65 frozen-transfer envelope screen."""

import numpy as np
import pytest

from scripts.benchmark_outer_envelope_predictor_round65 import envelope_factors


def _snapshot(phi, s2, j0, z0, nodes, z_tail=2.0):
    arrays = [np.zeros(3), np.asarray(phi, dtype=float), np.zeros(3),
              np.asarray(s2, dtype=float), np.ones(3), np.asarray(j0, dtype=np.int64),
              np.asarray(z0, dtype=float), np.zeros(1), np.zeros(3), np.zeros(3),
              np.asarray([1.0, 10.0, 100.0]), 1, len(nodes), 1, 0.01, z_tail,
              np.zeros((len(j0), len(nodes))), np.zeros((len(j0), len(nodes))),
              np.zeros((len(j0), len(nodes))), None, np.ones(3), 0.5,
              np.full(len(j0), -1.0), -1, 0.0, 0.0]
    arrays[12] = len(nodes)
    arrays[13] = 1
    return tuple(arrays)


def test_envelope_factor_matches_known_algebra():
    old = _snapshot([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], [0], [0.0], [0, 1, 2])
    new = _snapshot([0.0, 1.0, 2.0], [2.0, 4.0, 8.0], [0], [0.1], [0, 1, 2])
    factors = envelope_factors(old, new)
    # The endpoint S2 ratio is exactly cancelled here by the inverse-S2
    # factor in the initial state, leaving only the horizon-start shift.
    expected = np.full(3, np.exp(0.2))
    np.testing.assert_allclose(factors[0], expected, rtol=0.0, atol=1e-14)


def test_envelope_rejects_discrete_start_change():
    old = _snapshot([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], [0], [0.0], [0, 1, 2])
    new = _snapshot([0.0, 1.0, 2.0], [1.0, 2.0, 4.0], [1], [0.0], [0, 1, 2])
    with pytest.raises(ValueError, match="frequency start"):
        envelope_factors(old, new)
