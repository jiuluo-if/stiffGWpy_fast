"""TDD contract for the boundary-conditioned Kummer phase spike."""

from scripts.benchmark_phase_function_boundary_spike import _run_case


def test_boundary_conditioned_phase_reports_certification_metrics():
    result = _run_case("default", mode_count=1)
    assert result["candidate"] == "boundary_conditioned_kummer_phase"
    assert result["production_unchanged"] is True
    assert result["rows"]
    row = result["rows"][0]
    assert {
        "status",
        "rho_min",
        "rho_fit_min",
        "kummer_residual_max",
        "window_z_length",
    } <= row.keys()
