from scripts.benchmark_phase_increment_recurrence import _probe


def test_phase_increment_recurrence_is_finite_and_accurate():
    row = _probe(0.2, 3.0, 0.01, 256, 32)
    assert row["finite"] is True
    assert row["power_relative_error"] < 1e-8
