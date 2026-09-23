import numpy as np

from scripts.benchmark_phase_increment_kernel_twin import (
    _phase_segment_increment,
)
from stiffgwpy_fast import fast_sgwb as FS


def test_phase_increment_reanchor_one_matches_exact_transfer():
    base = (0.125, 0.875)
    exact = _phase_segment_increment(
        *base, 0.8, 1.6, 0.5, 0.25, 1, 32, 1.0, 0.0, 0.0)
    production = FS._phase_segment(*base, 0.8, 1.6, 0.5, 0.25)
    assert np.array_equal(np.asarray(exact[:2]), np.asarray(production))
    assert np.isfinite(exact[:2]).all()


def test_phase_increment_state_is_finite():
    state = _phase_segment_increment(
        0.0, 1.0, 0.2, 3.0, 0.005, 0.25, 32, 1.0, 0.0, 0.0, 0.0)
    assert np.isfinite(np.asarray(state)).all()
