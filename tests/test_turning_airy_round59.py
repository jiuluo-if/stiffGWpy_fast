"""Contract tests for the standalone turning-point Airy map."""

import numpy as np

from scripts.benchmark_turning_airy_round59 import _oracle, airy_transfer


def test_airy_connection_is_finite_for_simple_turning_interval():
    candidate = airy_transfer(0.0, 1.0, 0.222, 0.2245, 0.005)
    assert np.isfinite(candidate).all()


def test_airy_connection_is_compared_to_independent_oracle():
    candidate = airy_transfer(0.0, 1.0, 0.222, 0.2245, 0.005)
    oracle = _oracle(0.0, 1.0, 0.222, 0.2245, 0.005)
    cp = candidate[0] ** 2 + candidate[1] ** 2
    op = oracle[0] ** 2 + oracle[1] ** 2
    assert abs(cp - op) / max(abs(op), 1e-300) < 0.1
