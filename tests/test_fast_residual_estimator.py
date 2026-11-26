"""Regression for state propagation error in a curved background panel."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from stiffgwpy_fast import fast_sgwb as fast
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_residual_half_step_subdivision_uses_both_endpoints():
    x0, y0 = 0.3, -0.7
    z_start, z_mid, z_end = 2.0, 3.0, 5.0
    h_step, phase_max = 0.2, 0.1
    n_sub = fast._phase_substeps(h_step, z_mid, phase_max)
    assert n_sub > 1
    expected_x, expected_y = x0, y0
    h_sub = h_step / n_sub
    dz_half = 0.5 * (z_end - z_start)
    for index in range(n_sub):
        z_sub_mid = z_start + dz_half * (2 * index + 1) / n_sub
        expected_x, expected_y = fast.scaled_step(
            expected_x, expected_y, z_sub_mid, h_sub
        )

    actual_x, actual_y = fast._phase_segment_with_midpoint(
        x0, y0, z_start, z_mid, z_end, h_step, phase_max
    )
    assert actual_x == expected_x
    assert actual_y == expected_y


def test_fast_controls_internal_cr1_root_error_from_saved_reference():
    root = Path(__file__).resolve().parents[1]
    artifact = root / (
        "LOCAL_RESEARCH_REFERENCE"
        "validate_cr1_t0p004_0p005_r0000000000010641088010743262_"
        "0000000000010675_unseen_holdout_T00475_rlo_z10.json"
    )
    validation = json.loads(artifact.read_text(encoding="utf-8"))
    model = LCDM_SG(**validation["params"])

    assert fast.SGWB_iter_fast(
        model,
        config=fast.resolve_config("fast", threads=2),
        tol=1.0e-6,
        physical_guard=True,
    ) is model

    reference_root = float(validation["reference"]["root"])
    root_error = abs(float(model.DN_gw[-1]) - reference_root) / abs(reference_root)
    assert root_error <= 2.0e-4

    reference_signal = np.asarray(validation["reference"]["Ogw"]) - np.asarray(
        validation["reference"]["Oj"]
    )
    valid = np.isfinite(reference_signal) & (reference_signal > 0.0)
    spectrum_error = np.abs(
        np.asarray(model.log10OmegaGW)[valid] - np.log10(reference_signal[valid])
    )
    assert spectrum_error.size == int(valid.sum())
    assert float(np.max(spectrum_error)) <= 1.0e-3
