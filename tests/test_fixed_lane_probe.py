"""Contract tests for the explicit fixed-width lane feasibility probe."""

from scripts.benchmark_fixed_lane_probe import _lane4_output, _llvm_has_vector_lane


def test_explicit_lane_probe_is_deterministic_and_reports_llvm_evidence():
    first = _lane4_output()
    second = _lane4_output()
    assert first == second
    assert len(first) == 4
    assert all(value == value for value in first)
    assert isinstance(_llvm_has_vector_lane(), bool)
