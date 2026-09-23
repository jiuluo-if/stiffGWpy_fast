"""Contract tests for the pole-free projective-coordinate screen."""

from scripts.benchmark_cayley_riccati_spike import _cayley_probe


def test_cayley_probe_is_finite_and_deterministic():
    first = _cayley_probe(0.2, 3.0, 0.01, 256)
    second = _cayley_probe(0.2, 3.0, 0.01, 256)
    for key in ("z0", "z1", "h", "steps", "finite", "power_relative_error"):
        assert first[key] == second[key]
    assert first["finite"] is True
    assert first["power_relative_error"] < 1e-10
