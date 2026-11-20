"""Keep every high-density mode when a WKB basis inverse is unresolvable."""

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
    np.testing.assert_allclose(model.DN_gw[-1], 0.3121571825306297, rtol=2.0e-7, atol=0.0)
    linear_spectrum = model.Ogw_today - model.Oj_today
    assert np.all(linear_spectrum > 0.0)
    representable = linear_spectrum > 0.0
    np.testing.assert_allclose(
        model.log10OmegaGW[representable],
        np.log10(linear_spectrum[representable]),
        rtol=0.0,
        atol=2.0e-12,
    )
