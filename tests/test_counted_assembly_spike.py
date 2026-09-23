"""Contract tests for the independent counted assembly-state candidate."""

from scripts.benchmark_counted_assembly_spike import _run_kernel


def test_counted_assembly_preserves_outputs_and_nodes_on_default():
    row = _run_kernel("default", repeats=1, threads=2)
    assert row["bitwise_equal"] is True
    assert row["assembly_nodes_equal"] is True
    assert row["first_divergence"] is None
