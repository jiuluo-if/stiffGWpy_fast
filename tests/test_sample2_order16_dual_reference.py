from __future__ import annotations


def test_single_order_candidate_requires_two_converged_references_and_matching_tail():
    try:
        from scripts.validate_sample2_order16_reference import classify_candidate
    except (ImportError, ModuleNotFoundError):
        import pytest

        pytest.fail("dual-reference acceptance classifier is not implemented")

    candidate = {
        "status": "completed",
        "log10_omega": -15.0,
        "tail_used": False,
        "max_collocation_residual_scaled": 1e-11,
        "max_dense_residual_scaled": 2e-10,
    }
    cartesian = {"log10_omega": -15.00000002, "used_tail": False}
    prufer = {"log10_omega": -14.99999998, "used_tail": False}

    assert classify_candidate(candidate, cartesian, prufer)["accepted"] is True
    assert classify_candidate(candidate, cartesian, {**prufer, "used_tail": True})[
        "accepted"
    ] is False
    assert classify_candidate(candidate, cartesian, {**prufer, "log10_omega": -14.8})[
        "accepted"
    ] is False
