import inspect

import numpy as np
from scipy.interpolate import PchipInterpolator

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.fast_sgwb import (
    _phase_resolved_tail_weights,
    _rebuild_pchip_terminal_map,
    _scale_frequency_rows_by_log,
)
from stiffgwpy_fast.freq_adaptive import (
    _log_tensor_amplitude,
    _requires_log_tensor_power,
    construct_fast_frequency_grid,
    grid_independent_freqs,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def _underflow_model(**overrides):
    params = dict(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-316,
        n_t=17.38,
        cr=0.0,
        T_re=2e3,
        DN_re=10.0,
        kappa10=1e-2,
    )
    params.update(overrides)
    return LCDM_SG(**params)


def test_fast_recovers_blue_spectrum_after_nonzero_tensor_product_underflows():
    model = _underflow_model()
    assert model.derived_param["A_t"] == 0.0

    result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=5.0,
        freq_res=1.0,
        threads=1,
        accuracy_mode="fast",
        physical_guard=False,
    )

    assert result is model
    assert model.SGWB_converge
    assert np.isfinite(model.cosmo_param["DN_eff"])
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()
    assert len(model.f) > 120
    assert model.freq_grid_requested == "goal"
    assert model.freq_grid_used == "construct"
    assert np.any(model.Ogw_today > 0.0)


def test_phase_tail_blend_is_gated_by_log_amplitude_mode_and_ripple():
    z_today = np.array([4.9, 5.1, 7.0])
    disabled = _phase_resolved_tail_weights(
        z_today, z_tail=5.0, enabled=False
    )
    enabled = _phase_resolved_tail_weights(
        z_today, z_tail=5.0, enabled=True
    )

    assert np.array_equal(disabled, np.zeros_like(z_today))
    assert enabled[0] == 0.0
    assert 0.0 < enabled[1] < 1.0
    assert enabled[2] == 0.0


def test_underflow_log_spectrum_matches_deep_tail_holdout_bin():
    model = _underflow_model()
    model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=5.0,
        freq_res=1.0,
        threads=1,
        accuracy_mode="fast",
        physical_guard=False,
    )

    # Frozen scaled-Reference z8 row; its 235-node grid SHA is f0dae6e2... .
    # 冻结的缩放 Reference z8 频点；235 节点网格哈希为 f0dae6e2...。
    idx = int(np.argmin(np.abs(model.f - (-16.241944809171002))))
    assert abs(float(model.f[idx]) - (-16.241944809171002)) < 1e-9
    assert abs(float(model.log10OmegaGW[idx]) - (-332.2835223133478)) < 1e-3


def test_pchip_terminal_map_is_rebuilt_from_updated_spectrum():
    frequencies = np.array([-2.0, -1.0, 0.0])
    ogw = np.array([[1.0], [2.0], [4.0]])
    oj = np.zeros_like(ogw)
    indices = np.arange(frequencies.size, dtype=np.int64)

    before, _ = _rebuild_pchip_terminal_map(
        ogw, oj, indices, frequencies
    )
    ogw[1, -1] = 3.0
    after, _ = _rebuild_pchip_terminal_map(
        ogw, oj, indices, frequencies
    )
    expected = PchipInterpolator(frequencies, ogw[:, -1]).integrate(
        float(frequencies[0]), float(frequencies[-1])
    )

    assert after != before
    assert np.isclose(after, expected, rtol=1e-14, atol=0.0)


def test_underflow_mode_rejects_nonpositive_transfer_instead_of_returning_minus_inf(
    monkeypatch,
):
    original = FS.solve_kernel
    signature = inspect.signature(original.py_func)

    def force_zero_terminal_transfer(*args, **kwargs):
        result = original(*args, **kwargs)
        bound = signature.bind(*args, **kwargs)
        ogw, oj, opgw = (
            bound.arguments["Ogw"],
            bound.arguments["Oj"],
            bound.arguments["Opgw"],
        )
        ogw[0, -1] = oj[0, -1]
        opgw[0, -1] = 0.0
        return result

    monkeypatch.setattr(FS, "solve_kernel", force_zero_terminal_transfer)
    model = _underflow_model()

    result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=5.0,
        freq_res=1.0,
        threads=1,
        accuracy_mode="fast",
        physical_guard=False,
    )

    assert result is None
    assert model.fast_failure_reason == "nonpositive_log_amplitude_transfer"


def test_log_power_guard_checks_intermediate_and_final_representability():
    freqs = np.array([-18.459378758443304, 2.323238066410573])
    assert _requires_log_tensor_power(np.log(1e-323), 17.38, freqs, gp.f_piv)
    assert not _requires_log_tensor_power(np.log(2.100549e-11), 0.0, freqs, gp.f_piv)
    assert _requires_log_tensor_power(570.0, -18.5, freqs, gp.f_piv)
    assert _requires_log_tensor_power(100.0, 0.1, np.array([-400.0]), gp.f_piv)


def test_fast_local_frequency_grid_accepts_underflowed_amplitude():
    model = _underflow_model()
    model.gen_expansion()
    log_at = np.log(model.cosmo_param["r"]) + np.log(model.cosmo_param["A_s"])

    freqs = construct_fast_frequency_grid(model, freq_res=1.0, log_amplitude=log_at)

    assert np.isfinite(freqs).all()
    assert np.all(np.diff(freqs) < 0.0)
    assert freqs[-1] >= min(model.f_hor)
    assert freqs[-1] - min(model.f_hor) < 0.2


def test_grid_independent_builder_accepts_overflowed_stored_product():
    model = _underflow_model(A_s=1e200, r=1e200, n_t=200.0)
    assert np.isinf(model.derived_param["A_t"])
    assert np.isfinite(_log_tensor_amplitude(model.cosmo_param, model.derived_param))

    model.gen_expansion()
    freqs, _, _ = grid_independent_freqs(model, allow_log_amplitude=True)

    assert np.isfinite(freqs).all()
    assert np.all(np.diff(freqs) < 0.0)


def test_log_scaling_does_not_hide_nonfinite_propagation():
    values = np.array([[np.nan, 2.0]])
    _scale_frequency_rows_by_log(values, np.array([0.0]))
    assert np.isnan(values[0, 0])
    assert values[0, 1] == 2.0
