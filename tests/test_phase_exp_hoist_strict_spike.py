"""Contract tests for the strict repeated-exp standalone candidate."""

from scripts.benchmark_phase_exp_hoist_strict_spike import _run_kernel


def test_strict_repeated_exp_candidate_is_bitwise_equal_on_default():
    row = _run_kernel("default", repeats=1, threads=2)
    assert row["bitwise_equal"] is True
    assert row["first_divergence"] is None
