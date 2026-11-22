import numpy as np

from stiffgwpy_fast import exact_background, fast_sgwb
from stiffgwpy_fast.exact_background import H2_vec
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, gen_fast, resolve_config
from stiffgwpy_fast.global_param import N_10
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
        1.0, 800.0, np.nan, 1.0, 1.0, 1.0
    )

    assert not recovered


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


def test_fast_solves_when_reheat_exponential_overflows_but_density_is_finite():
    model = LCDM_SG(
        Omega_bh2=0.0223828,
        Omega_ch2=0.1201075,
        H0=67.32117,
        DN_eff=0.0,
        A_s=2.100549e-9,
        r=1e-300,
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

    assert result is model
    assert model.fast_failure_reason is None
    assert np.isfinite(model.DN_gw).all()
    assert np.isfinite(model.f).all()
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
