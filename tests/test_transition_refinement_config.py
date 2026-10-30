import pytest

from stiffgwpy_fast import exact_background as EB
from stiffgwpy_fast import fast_sgwb as FS


def test_transition_refine_config_reaches_grid_builder(monkeypatch, model):
    calls = []
    original = EB.build_kink_refined_grid

    def capture(m, h, refine_radius=2.0, refine_factor=8):
        calls.append((refine_radius, refine_factor))
        return original(m, h, refine_radius, refine_factor)

    monkeypatch.setattr(EB, 'build_kink_refined_grid', capture)
    config = FS.resolve_config(
        'fast', threads=1, h=0.005, transition_refine_radius=0.25,
        transition_refine_factor=8,
    )

    result = FS.SGWB_iter_fast(
        model,
        config=config,
        kink_split=False,
        transition_refine=True,
        freq_grid='goal',
        tol=1e-6,
    )

    assert result is model
    assert calls == [(0.25, 8), (0.25, 8)]


def test_transition_refine_and_kink_split_use_actual_nonuniform_steps(monkeypatch, model):
    """The refined grid passes real local steps to propagation；细化网格须传递真实局部步长。"""
    import numpy as np

    captured = []
    original = FS.solve_kernel

    def capture(*args):
        nv, h_arr = args[0], args[19]
        captured.append((np.asarray(nv).copy(), None if h_arr is None else np.asarray(h_arr).copy()))
        return original(*args)

    monkeypatch.setattr(FS, 'solve_kernel', capture)
    config = FS.resolve_config(
        'fast', threads=1, h=0.02, z_tail=5.0, freq_grid='construct',
        transition_refine_radius=0.25, transition_refine_factor=8,
    )

    assert FS.SGWB_iter_fast(
        model, tol=1e-6, config=config, transition_refine=True, kink_split=True,
    ) is model

    assert captured
    nv, h_arr = captured[-1]
    assert h_arr is not None
    np.testing.assert_allclose(h_arr, np.diff(nv), rtol=0.0, atol=1e-14)
    assert np.unique(np.round(h_arr, 12)).size > 1


def test_transition_refine_controls_survive_base_config_resolution():
    base = FS.resolve_config(
        'fast', transition_refine_radius=0.25, transition_refine_factor=8)
    derived = FS.resolve_config(
        base=base, threads=1, transition_refine_radius=0.5)

    assert derived.transition_refine_radius == 0.5
    assert derived.transition_refine_factor == 8


@pytest.mark.parametrize(
    ('overrides', 'message'),
    [
        ({'transition_refine_radius': -0.1}, 'transition_refine_radius'),
        ({'transition_refine_factor': 0}, 'transition_refine_factor'),
        ({'transition_refine_factor': 129}, 'transition_refine_factor'),
    ],
)
def test_transition_refine_config_rejects_invalid_controls(overrides, message):
    with pytest.raises(ValueError, match=message):
        FS.resolve_config('fast', **overrides)


def test_transition_refine_factor_128_is_allowed_within_grid_budget():
    config = FS.resolve_config(
        'fast', h=0.04, transition_refine_radius=0.125,
        transition_refine_factor=128,
    )

    assert config.transition_refine_factor == 128


def test_transition_refine_config_bounds_combined_grid_size():
    with pytest.raises(ValueError, match='refinement grid'):
        FS.resolve_config(
            'fast', h=1e-4, transition_refine_radius=1.0,
            transition_refine_factor=64,
        )
