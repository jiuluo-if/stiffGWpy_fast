from __future__ import annotations

import importlib

import pytest


def test_best_phase_start_candidate_requires_reference_error_and_tail_match():
    try:
        probe = importlib.import_module("scripts.probe_sample2_phase_start_offsets")
    except ModuleNotFoundError:
        pytest.fail("phase-start candidate selector is not implemented")

    select = getattr(probe, "select_best_candidate", None)
    assert select is not None, "phase-start candidate selector is not implemented"
    trials = [
        {"phase_start": 0.875, "candidate": {"status": "failed"}},
        {
            "phase_start": 1.0,
            "candidate": {
                "status": "completed",
                "log10_omega": -15.0,
                "abs_dex_error_vs_cartesian": 2e-8,
                "tail_used": False,
                "max_collocation_residual_scaled": 1e-11,
                "max_dense_residual_scaled": 2e-10,
            },
        },
        {
            "phase_start": 1.125,
            "candidate": {
                "status": "completed",
                "log10_omega": -15.0,
                "abs_dex_error_vs_cartesian": 3e-8,
                "tail_used": False,
                "max_collocation_residual_scaled": 5e-11,
                "max_dense_residual_scaled": 1e-10,
            },
        },
        {
            "phase_start": 1.25,
            "candidate": {
                "status": "completed",
                "log10_omega": -15.0,
                "abs_dex_error_vs_cartesian": 2e-3,
                "tail_used": False,
                "max_collocation_residual_scaled": 1e-14,
                "max_dense_residual_scaled": 1e-14,
            },
        },
        {
            "phase_start": 1.375,
            "candidate": {
                "status": "completed",
                "log10_omega": -15.0,
                "abs_dex_error_vs_cartesian": 2e-8,
                "tail_used": True,
                "max_collocation_residual_scaled": 1e-14,
                "max_dense_residual_scaled": 1e-14,
            },
        },
    ]

    best = select(trials, reference_used_tail=False)

    assert best["phase_start"] == 1.125
