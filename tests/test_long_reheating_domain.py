import numpy as np
import pytest

from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


@pytest.mark.parametrize("duration", [680.0, 800.0])
def test_inactive_matter_prefix_preserves_complete_tensor_solution(duration):
    # The blue cutoff and every initial condition lie after reheating here.
    params = dict(cr=0.0, T_re=2000.0, kappa10=0.01, n_t=2.0, r=0.01)
    models = [LCDM_SG(DN_re=d, **params) for d in (16.0, duration)]
    for model in models:
        assert model.SGWB_iter(engine="fast", accuracy_mode="fast", threads=2,
                               fallback=False, physical_guard=False) is not None
        assert np.all(np.isfinite(model.log10OmegaGW))
        assert model.f.max() + 3.0 < model.f_re
    short, long = models
    np.testing.assert_allclose(long.f, short.f, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(long.log10OmegaGW, short.log10OmegaGW, atol=1e-7, rtol=0.0)
    offset = round((duration - 16.0) / 0.005)
    for name in ("DN_gw", "g2", "w2"):
        values = getattr(long, name)
        assert np.all(values[:offset] == 0.0)
        np.testing.assert_allclose(values[offset:], getattr(short, name),
                                   atol=1e-18, rtol=1e-5)


@pytest.mark.parametrize("T_re", [10.0, 2000.0, 1.0e5])
@pytest.mark.parametrize("kappa10", [0.001, 0.01, 0.1])
def test_long_reheating_cutoff_slice_remains_evaluable(T_re, kappa10):
    seed = LCDM_SG(cr=0.0, T_re=T_re, kappa10=kappa10, r=0.01,
                   n_t=1.0, A_s=2.1e-9, DN_re=0.0)
    N_re = float(seed.derived_param["N_re"])
    H_over_2pi = float(seed.derived_param["H_0"]) / (2.0 * np.pi)
    amplitude = (gp.f_piv / (1.01 * H_over_2pi)) / 0.01
    model = LCDM_SG(
        cr=0.0,
        T_re=T_re,
        kappa10=kappa10,
        r=0.01,
        n_t=1.0,
        A_s=amplitude,
        DN_re=740.0 - N_re,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2),
        tol=1e-6,
    )

    assert result is model
    assert model.derived_param["N_inf"] == pytest.approx(740.0)
    assert model.fast_failure_reason is None
    assert np.all(np.isfinite(model.log10OmegaGW))

