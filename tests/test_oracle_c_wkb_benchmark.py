import pytest

from scripts import benchmark_oracle_c_wkb as oracle_c


def test_oracle_c_does_not_write_reference_for_rejected_fast_point(
        monkeypatch, tmp_path):
    def reject_fast_point(model, **kwargs):
        model.fast_failure_reason = 'shared_Neff_guard'
        model.f = [-1.0]
        return None

    monkeypatch.setattr(oracle_c.FS, 'apply_accuracy_mode', lambda mode: None)
    monkeypatch.setattr(oracle_c.FS, 'SGWB_iter_fast', reject_fast_point)

    with pytest.raises(SystemExit, match='shared_Neff_guard'):
        oracle_c.main([
            '--points', 'default',
            '--out-dir', str(tmp_path),
        ])

    assert list(tmp_path.iterdir()) == []
