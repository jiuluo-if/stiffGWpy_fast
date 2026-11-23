import math
import warnings

import numpy as np

from stiffgwpy_fast import exact_background, fast_sgwb
from stiffgwpy_fast.exact_background import H2_vec
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, gen_fast, resolve_config
from stiffgwpy_fast.global_param import N_10, Neff0, Omega_nh2
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def _extreme_reheating_model():
    return LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-300,
        n_t=0.0,
        cr=0.0,
        T_re=1e101,
        DN_re=10.0,
        kappa10=0.01,
    )


def test_fast_keeps_finite_background_ratios_when_early_h2_exceeds_float64():
    model = _extreme_reheating_model()
    config = resolve_config("fast", threads=2)

    today_h2 = H2_vec(np.array([model.derived_param["N_inf"]]), model, 0.0)
    assert np.isfinite(today_h2).all()
    assert np.all(today_h2 > 0.0)
    early_h2 = H2_vec(np.array([0.0]), model, 0.0)
    assert np.isposinf(early_h2[0])

    gen_fast(model, config.h)
    assert np.isfinite(model.sigma).all()
    assert np.isfinite(model.f_hor).all()

    result = SGWB_iter_fast(
        model, config=config, tol=1e-6, physical_guard=True
    )
    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()


def test_h2_identity_preserves_nonpositive_h2_failure_classification(monkeypatch):
    model = _extreme_reheating_model()
    n_re = model.derived_param["N_inf"] - model.derived_param["N_re"]
    nodes = np.array([n_re, n_re + 1.0])
    midpoints = np.array([n_re + 0.5])

    monkeypatch.setattr(
        exact_background,
        "_log_H2_vec",
        lambda values, *_args: np.full(np.asarray(values).shape, np.nan),
    )
    monkeypatch.setattr(
        exact_background,
        "H2_vec",
        lambda values, *_args: -np.ones(np.asarray(values).shape),
    )

    result, failure, _index = exact_background._h2_identity_primitives(
        model, nodes, 0.0, midpoints
    )

    assert result is None
    assert failure == "nonpositive_h2"


def test_stiff_dominance_refuses_nan_remainder_estimates():
    recovered, _log_h2 = fast_sgwb._stiff_dominated_log_h2(
        1.0, 0.0, 800.0, 100.0, np.nan, 1.0, 1.0, 1.0, np.nan
    )

    assert not recovered


def test_scaled_tail_moment_preserves_nan_failure():
    moment = fast_sgwb._square_times_power_numba(np.nan, 1.0, 36.0)

    assert np.isnan(moment)


def test_log_scaled_state_energy_matches_unscaled_quadratic_moments():
    s2 = 2.5
    xh, yh, zz = 0.3, -0.7, 0.2
    state_scale = 1.7
    power = 0.02
    expected = [np.zeros((1, 1)) for _ in range(3)]
    actual = [np.zeros((1, 1)) for _ in range(3)]
    log_base_signal = np.full(1, -math.inf)

    fast_sgwb.assemble_main(
        *expected, 0, 0, s2, state_scale * xh, state_scale * yh, zz, power
    )
    fast_sgwb.assemble_main_log_scaled(
        *actual, 0, 0, 0, s2, xh, yh, zz,
        math.log(state_scale), math.log(power), log_base_signal,
    )

    for observed, reference in zip(actual, expected):
        np.testing.assert_allclose(observed, reference, rtol=2e-15, atol=0.0)
    assert math.isclose(
        log_base_signal[0], math.log(expected[0][0, 0] - expected[1][0, 0]),
        rel_tol=0.0, abs_tol=2e-15,
    )


def test_signed_log_moment_does_not_hide_nan_behind_zero_sign():
    assert np.isnan(fast_sgwb._signed_exp_from_log(math.nan, 0.0))


def test_log_scaled_history_preserves_signed_products_and_nan_scale():
    values = np.array([1e-319, -1e-319, 0.0, 0.0])
    log_scale = np.array([710.0, 710.0, 710.0, math.nan])

    result = fast_sgwb._scale_array_by_log_exponent(values, log_scale, 11.3)

    expected_magnitude = math.exp(math.log(1e-319) + 710.0 + 11.3)
    assert math.isclose(result[0], expected_magnitude, rel_tol=2e-14)
    assert math.isclose(result[1], -expected_magnitude, rel_tol=2e-14)
    assert result[2] == 0.0
    assert np.isnan(result[3])


def test_log_scaled_tail_moment_keeps_representable_final_product():
    moment = fast_sgwb._log_square_times_power_numba(
        math.log(1e200), 2.1e-309, 36.0
    )

    assert np.isfinite(moment)
    assert math.isclose(
        math.log(moment),
        2.0 * math.log(1e200) + math.log(2.1e-309) - math.log(36.0),
        rel_tol=0.0,
        abs_tol=2e-13,
    )


def test_extreme_reheating_keeps_log_of_unrepresentable_reheat_stiff_density():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-300,
        n_t=0.0,
        cr=0.0,
        T_re=1e200,
        DN_re=10.0,
        kappa10=0.01,
    )

    derived = model.derived_param

    assert np.isfinite(derived["N_inf"])
    assert np.isposinf(model.rhostiff_re)
    assert np.isfinite(model.log_rhostiff_re)


def test_fast_scales_relative_background_tail_moment_at_unrepresentable_reheat():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-308,
        n_t=0.0,
        cr=0.0,
        T_re=1e200,
        DN_re=10.0,
        kappa10=1e-200,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()


def test_fast_preserves_subnormal_tail_spectrum_at_unrepresentable_reheat():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-310,
        n_t=0.0,
        cr=0.0,
        T_re=1e200,
        DN_re=10.0,
        kappa10=1e-200,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.log10OmegaGW).all()


def test_fast_scales_tensor_state_before_unrepresentable_transfer_amplitude():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-308,
        n_t=0.0,
        cr=0.0,
        T_re=10.0**307.25,
        DN_re=10.0,
        kappa10=1e-308,
    )
    assert np.isfinite(model.derived_param["N_inf"])
    assert np.isfinite(exact_background._log_H2_vec(np.array([0.0]), model, 0.0)).all()

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()


def test_fast_classifies_log_scaled_final_spectrum_underflow():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-314,
        n_t=0.0,
        cr=0.0,
        T_re=10.0**307.25,
        DN_re=10.0,
        kappa10=1e-308,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is None
    assert model.fast_failure_reason == "nonrepresentable_output"
    assert model.fast_failure_detail == "final_spectrum_underflow"


def test_background_grid_uses_log_h2_when_eN_intermediate_overflows():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-300,
        n_t=0.0,
        cr=0.0,
        T_re=1e300,
        DN_re=10.0,
        kappa10=1e-308,
    )

    log_h2 = exact_background._log_H2_vec(np.array([0.0]), model, 0.0)
    assert np.isfinite(log_h2).all()

    generated = gen_fast(model, h=0.01)
    assert generated is not None
    assert np.isfinite(model.sigma).all()
    assert np.isfinite(model.f_hor).all()

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )
    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()


def test_reheating_stiff_density_recovers_underflowed_prefactor_product():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-300,
        n_t=0.0,
        cr=0.0,
        T_re=1e145,
        DN_re=10.0,
        kappa10=1e-320,
    )

    _derived = model.derived_param
    assert model.rhostiff_re > 0.0
    assert np.isfinite(model.rhostiff_re)
    assert np.isfinite(model.log_rhostiff_re)
    assert np.isclose(np.log(model.rhostiff_re), model.log_rhostiff_re, atol=2e-13)


def test_fast_classifies_final_spectrum_underflow_after_scaled_tail_assembly():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-310,
        n_t=0.0,
        cr=0.0,
        T_re=1e155,
        DN_re=10.0,
        kappa10=1e-100,
    )

    derived = model.derived_param
    assert np.isfinite(derived["N_inf"])
    assert np.isfinite(model.rhostiff_re)
    assert np.isfinite(model.log_rhostiff_re)
    assert 2.0 * (derived["N_re"] - N_10) > np.log(np.finfo(np.float64).max)

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is None
    assert model.fast_failure_reason == "nonrepresentable_output"
    assert model.fast_failure_detail == "final_spectrum_underflow"


def test_fast_solves_with_subnormal_tensor_power_and_representable_spectrum():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-308,
        n_t=0.0,
        cr=0.0,
        T_re=1e155,
        DN_re=10.0,
        kappa10=1e-100,
    )

    result = SGWB_iter_fast(
        model,
        config=resolve_config("fast", threads=2, freq_grid="adaptive"),
        tol=1e-6,
        freq_grid_target=1e-6,
        physical_guard=True,
    )

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()


def test_auto_relative_primitive_scale_is_used_by_tensor_propagation(monkeypatch):
    params = {
        "Omega_bh2": 0.0223828,
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
        "r": 1e-30,
        "n_t": 0.0,
        "cr": 0.0,
        "T_re": 2000.0,
        "DN_re": 10.0,
        "kappa10": 0.01,
    }
    config = resolve_config("fast", threads=2)

    baseline = LCDM_SG(**params)
    assert SGWB_iter_fast(
        baseline, config=config, tol=1e-6, physical_guard=True
    ) is baseline
    baseline_spectrum = baseline.log10OmegaGW.copy()
    baseline_root = float(baseline.DN_gw[-1])

    original_primitive = exact_background.fast_phi_s2_split
    original_kernel = fast_sgwb.solve_kernel
    relative_arguments = []

    def force_relative_primitive(model, *args, **kwargs):
        primitive = original_primitive(model, *args, **kwargs)
        model._fast_relative_background_required = True
        return (primitive[0], primitive[1], np.ones_like(primitive[2]),
                np.ones_like(primitive[3]), *primitive[4:])

    def capture_relative_kernel(*args, **kwargs):
        relative = kwargs.get("relative_background")
        relative_arguments.append(relative is not None)
        return original_kernel(*args, **kwargs)

    monkeypatch.setattr(exact_background, "fast_phi_s2_split", force_relative_primitive)
    monkeypatch.setattr(fast_sgwb, "solve_kernel", capture_relative_kernel)
    scaled = LCDM_SG(**params)
    assert SGWB_iter_fast(
        scaled, config=config, tol=1e-6, physical_guard=True
    ) is scaled

    assert any(relative_arguments)
    assert abs(float(scaled.DN_gw[-1]) - baseline_root) <= 1e-12 * baseline_root
    assert np.array_equal(scaled.f, baseline.f)
    assert np.max(np.abs(scaled.log10OmegaGW - baseline_spectrum)) <= 1e-10


def test_fast_combines_history_scale_before_exponentiating_at_extreme_reheat():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=10.0**-308.5,
        n_t=0.0,
        cr=0.0,
        T_re=1e200,
        DN_re=10.0,
        kappa10=10.0**-91.5,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = SGWB_iter_fast(
            model,
            config=resolve_config("fast", threads=2, freq_grid="adaptive"),
            tol=1e-6,
            freq_grid_target=1e-6,
            physical_guard=True,
        )

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.log10OmegaGW).all()
    assert not any("overflow encountered in exp" in str(item.message) for item in caught)

    scale = 2.0 * (model.f_hor - model.f_hor[-1]) * math.log(10.0)
    scale += 2.0 * (model.Nv - model.Nv[-1])
    overflow_indices = np.flatnonzero(scale >= math.log(np.finfo(np.float64).max))
    assert overflow_indices.size > 0
    omega_nu = Omega_nh2 / model.derived_param["h"] ** 2
    expected_log = (
        np.log(model.g2[overflow_indices])
        + scale[overflow_indices]
        + math.log(Neff0 / omega_nu)
    )
    expected = np.exp(expected_log)
    np.testing.assert_allclose(
        model.DN_gw[overflow_indices], expected, rtol=2e-14, atol=0.0
    )
