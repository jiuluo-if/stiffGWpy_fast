"""TDD contract for the Levin interaction-picture phase screen."""

from scripts.benchmark_levin_envelope_spike import _run_case


def test_levin_envelope_reports_residual_and_transfer_metrics():
    result = _run_case("default", mode_count=1)
    assert result["candidate"] == "levin_interaction_envelope"
    assert result["production_unchanged"] is True
    assert result["rows"]
    assert {
        "status",
        "envelope_fit_relative_max",
        "transfer_relative_error",
        "window_z_length",
    } <= result["rows"][0].keys()
