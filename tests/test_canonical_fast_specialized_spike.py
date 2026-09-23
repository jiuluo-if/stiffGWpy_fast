"""Contract tests for the guarded canonical-fast specialization spike."""

from scripts.benchmark_canonical_fast_specialized_spike import _run_kernel


def test_canonical_specialized_kernel_is_bitwise_equal_on_default():
    row = _run_kernel("default", repeats=1, threads=2)
    assert row["bitwise_equal"] is True
    assert row["first_divergence"] is None
