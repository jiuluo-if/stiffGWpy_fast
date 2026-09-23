"""Contract tests for the standalone grouped SoA propagation spike."""

from scripts.benchmark_grouped_soa_spike import _run_kernel


def test_grouped_soa_kernel_is_bitwise_equal_on_default():
    row = _run_kernel("default", repeats=1, threads=2, bucket_width=32)
    assert row["bitwise_equal"] is True
    assert row["first_divergence"] is None
    assert row["work_overhead"] < 1.10
