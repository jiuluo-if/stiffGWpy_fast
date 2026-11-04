from __future__ import annotations

from scripts import validate_sample2_residual_orders as validation


def test_single_order_candidate_requires_two_converged_references_and_matching_tail():
    try:
        from scripts.validate_sample2_residual_orders import classify_candidate
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
    assert classify_candidate(
        {**candidate, "max_dense_residual_scaled": 2e-9}, cartesian, prufer
    )["accepted"] is False
    assert classify_candidate(
        {**candidate, "max_collocation_residual_scaled": float("nan")}, cartesian, prufer
    )["accepted"] is False


def test_dual_reference_selector_keeps_every_completed_order_only():
    select = getattr(validation, "select_completed_candidates", None)
    assert select is not None, "completed-order selector is not implemented"
    rows = [
        {"order": 16, "status": "completed"},
        {"order": 24, "status": "failed"},
        {"order": 32, "status": "completed"},
    ]
    assert select(rows) == [rows[0], rows[2]]
