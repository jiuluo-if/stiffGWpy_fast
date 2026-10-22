import json

import numpy as np
import pytest

from scripts import benchmark_fast_true_error as fast_error


@pytest.mark.parametrize(
    ('anchor_dn_eff', 'fast_dn_eff', 'anchor_frequencies',
     'fast_frequencies', 'message'),
    [
        (0.1, 0.2, [0.0], [0.0], 'DN_eff mismatch'),
        (0.2, 0.2, [0.0, 1.0], [0.0], 'frequency grid mismatch'),
    ],
)
def test_fast_error_benchmark_rejects_mismatched_reference(
        monkeypatch, tmp_path, anchor_dn_eff, fast_dn_eff,
        anchor_frequencies, fast_frequencies, message):
    anchor_dir = tmp_path / 'anchors'
    anchor_dir.mkdir()
    anchor = {
        'DN_eff': anchor_dn_eff,
        'frequencies': anchor_frequencies,
        'configuration': {
            'path': 'kink-split',
            'h': 0.005,
            'transition_refine_radius': 2.0,
            'transition_refine_factor': 8,
        },
        'summary': {},
    }
    (anchor_dir / 'oracle_c_wkb_default.json').write_text(
        json.dumps(anchor), encoding='utf-8')

    def solve_with_anchor_mismatch(model, **kwargs):
        model.SGWB_converge = True
        model.DN_gw = np.array([fast_dn_eff])
        model.cosmo_param['DN_eff'] = fast_dn_eff
        model.f = np.asarray(fast_frequencies)
        return model

    monkeypatch.setattr(
        fast_error.FS, 'apply_accuracy_mode',
        lambda mode: {'tol': 1e-6, 'freq_res': 1.0})
    monkeypatch.setattr(
        fast_error.FS, 'SGWB_iter_fast', solve_with_anchor_mismatch)
    out_path = tmp_path / 'must_not_be_written.json'

    with pytest.raises(RuntimeError, match=message):
        fast_error.main([
            '--points', 'default',
            '--anchor-dir', str(anchor_dir),
            '--out', str(out_path),
        ])

    assert not out_path.exists()


def test_fast_error_benchmark_rejects_refinement_mismatch_before_solving(
        monkeypatch, tmp_path):
    anchor_dir = tmp_path / 'anchors'
    anchor_dir.mkdir()
    anchor = {
        'configuration': {
            'path': 'transition-refine',
            'h': 0.0025,
            'transition_refine_radius': 2.0,
            'transition_refine_factor': 8,
        },
        'summary': {},
    }
    (anchor_dir / 'oracle_c_wkb_default.json').write_text(
        json.dumps(anchor), encoding='utf-8')
    monkeypatch.setattr(
        fast_error.FS, 'apply_accuracy_mode',
        lambda mode: {'tol': 1e-6, 'freq_res': 1.0})
    monkeypatch.setattr(
        fast_error.FS, 'SGWB_iter_fast',
        lambda *args, **kwargs: pytest.fail('solver must not run'))
    out_path = tmp_path / 'must_not_be_written.json'

    with pytest.raises(RuntimeError, match='transition_refine_radius mismatch'):
        fast_error.main([
            '--points', 'default',
            '--anchor-dir', str(anchor_dir),
            '--h', '0.0025',
            '--transition-refine',
            '--transition-refine-radius', '0.25',
            '--out', str(out_path),
        ])

    assert not out_path.exists()
