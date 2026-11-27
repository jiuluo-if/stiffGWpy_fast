"""Keep every high-density mode when a WKB basis inverse is unresolvable."""

import json
from pathlib import Path

import numpy as np

from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_fast_completes_high_density_modes_after_wkb_inverse_fallback():
    model = LCDM_SG(
        Omega_bh2=3.0e6,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=0.01,
        n_t=0.0,
        cr=0.0,
        T_re=2000.0,
        DN_re=10.0,
        kappa10=0.01,
    )

    solved = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        tol=1.0e-6,
        physical_guard=True,
    )

    assert solved is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.log10OmegaGW).all()
    assert np.isfinite(model.DN_gw).all()
    reference_path = (
        Path(__file__).resolve().parents[1]
        / "tests/data/validate_cr0_3m_wkb_guard_z8_multi.json"
    )
    reference = json.loads(reference_path.read_text(encoding="utf-8"))["cases"][0]
    closed_root = float(reference["closed_reference_root"])
    root_error = abs(float(model.DN_gw[-1]) - closed_root) / abs(closed_root)
    assert root_error <= 2.0e-4

    raw = reference["raw"]
    reference_freqs = np.asarray(raw["freqs_log10_hz"])
    reference_spectrum = np.asarray(raw["reference_log10_OmegaGW_at_fast_root"])
    candidate_freqs = np.asarray(model.f)
    order = np.argsort(reference_freqs)
    reference_at_candidate_freqs = np.interp(
        candidate_freqs,
        reference_freqs[order],
        reference_spectrum[order],
    )
    spectrum_error = np.abs(
        np.asarray(model.log10OmegaGW) - reference_at_candidate_freqs
    )
    assert float(np.max(spectrum_error)) <= 1.0e-3

    linear_spectrum = model.Ogw_today - model.Oj_today
    assert np.all(linear_spectrum > 0.0)
    representable = linear_spectrum > 0.0
    np.testing.assert_allclose(
        model.log10OmegaGW[representable],
        np.log10(linear_spectrum[representable]),
        rtol=0.0,
        atol=2.0e-12,
    )
