"""Contract test for the standalone Kummer phase-function feasibility spike."""

from scripts.benchmark_phase_function_chebyshev_spike import _run_case


def test_phase_function_chebyshev_screen_detects_insufficient_compression():
    row = _run_case("default", mode_count=4)
    assert row["usable_modes"] > 0
    assert row["positive_modes"] == row["usable_modes"]
    assert row["chebyshev_residual_max"] > 1e-2
