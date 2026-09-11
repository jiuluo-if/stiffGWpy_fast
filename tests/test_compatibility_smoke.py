"""跨 Python 版本运行的轻量兼容性 smoke tests。"""

import numpy as np

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_public_import_and_mode_normalization():
    assert FS.normalize_accuracy_mode('production') == 'production'
    assert FS.normalize_accuracy_mode('fast') == 'fast'
    assert FS.get_settings()['threads'] >= 1


def test_quadrature_helper_is_deterministic():
    freqs = np.linspace(-3.0, 1.0, 5)
    values = np.exp(freqs)
    first = FS.integrate_frequency_quadrature(freqs, values, 'simpson')
    second = FS.integrate_frequency_quadrature(freqs, values, 'simpson')
    assert first == second


def test_tiny_fast_smoke_uses_resource_contract(fast_settings):
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.set_threads(1)
    FS.SGWB_iter_fast(model, tol=1e-5, freq_grid='goal')
    assert np.isfinite(model.DN_gw[-1])
    assert model.SGWB_converge
