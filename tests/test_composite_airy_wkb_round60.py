"""Contract tests for the composite Airy/WKB feasibility screen."""

from scripts.benchmark_composite_airy_wkb_round60 import _probe


def test_composite_probe_reports_a_real_regime():
    row = _probe("default")
    assert row["intervals"] == 4000
    assert row["counts"]["turning_airy"] == 1
    assert row["counts"]["oscillatory_wkb"] > 0
