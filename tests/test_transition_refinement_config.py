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
    ],
)
def test_transition_refine_config_rejects_invalid_controls(overrides, message):
    with pytest.raises(ValueError, match=message):
        FS.resolve_config('fast', **overrides)


def test_transition_refine_config_bounds_combined_grid_size():
    with pytest.raises(ValueError, match='refinement grid'):
        FS.resolve_config(
            'fast', h=1e-4, transition_refine_radius=1.0,
            transition_refine_factor=64,
        )
