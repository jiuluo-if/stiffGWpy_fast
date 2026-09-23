"""TDD contract for the Round 56 safeguarded scalar secant screen."""


def test_secant_predictor_returns_bracketed_finite_point():
    from scripts.benchmark_safeguarded_outer_secant_round56 import _secant_predict

    value = _secant_predict(0.0, 0.2, 0.5, -0.1, 0.0, 1.0)
    assert 0.0 <= value <= 1.0


def test_secant_predictor_falls_back_to_midpoint_on_degenerate_slope():
    from scripts.benchmark_safeguarded_outer_secant_round56 import _secant_predict

    assert _secant_predict(0.2, 0.1, 0.2, 0.1, 0.0, 1.0) == 0.5


def test_secant_predictor_is_deterministic():
    from scripts.benchmark_safeguarded_outer_secant_round56 import _secant_predict

    args = (0.1, 0.04, 0.4, -0.02, 0.0, 1.0)
    assert _secant_predict(*args) == _secant_predict(*args)
