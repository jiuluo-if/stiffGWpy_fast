"""Contract tests for the standalone native-interval WKB map."""

import math

import numpy as np

from scripts.benchmark_wkb_phase_integral_round58 import (
    _oracle,
    wkb_transfer,
)


def test_wkb_transfer_is_finite_on_representative_interval():
    result = wkb_transfer(0.1, 1.0, 0.4, 0.4025, 0.005)
    assert np.isfinite(result).all()


def test_wkb_transfer_has_independent_oracle_scale():
    candidate = wkb_transfer(0.1, 1.0, 0.4, 0.4025, 0.005)
    oracle = _oracle(0.1, 1.0, 0.4, 0.4025, 0.005)
    assert math.isfinite(candidate[0]) and math.isfinite(candidate[1])
    candidate_power = candidate[0] ** 2 + candidate[1] ** 2
    oracle_power = oracle[0] ** 2 + oracle[1] ** 2
    # This is a feasibility contract, not the production acceptance gate;
    # the actual precision gate remains the project reference/oracle screen.
    assert abs(candidate_power - oracle_power) / oracle_power < 1.0e-2
