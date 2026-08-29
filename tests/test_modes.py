# -*- coding: utf-8 -*-
import pytest

from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb as FS
from stiffgwpy.stiff_SGWB import _sgwb_pool_size


@pytest.fixture
def fast_settings():
    """Snapshot and restore the process-global fast-solver settings."""
    saved = FS.get_settings()
    yield saved
    FS.set_threads(saved['threads'])
    FS.set_col_step(saved['col_step'])
    FS.set_h(saved['h'])
    FS.set_z_tail(saved['z_tail'])


def test_accuracy_modes_valid(fast_settings):
    assert set(FS.ACCURACY_MODES) == {'reference', 'production', 'ultra-fast'}
    for cfg in FS.ACCURACY_MODES.values():
        assert 1 <= cfg['col_step'] <= 8
        assert 1e-4 <= cfg['h'] <= 0.1
        assert 2.0 <= cfg['z_tail'] <= 15.0
        assert cfg['freq_res'] >= 1.0
        assert 0 < cfg['tol'] < 1
        assert 1 <= cfg['threads']
    assert FS.ACCURACY_MODES['reference']['h'] < FS.ACCURACY_MODES['production']['h']
    assert FS.ACCURACY_MODES['production']['z_tail'] > FS.ACCURACY_MODES['ultra-fast']['z_tail']


def test_apply_accuracy_mode_sets_state(fast_settings):
    for name, cfg in FS.ACCURACY_MODES.items():
        assert FS.apply_accuracy_mode(name) == cfg
        st = FS.get_settings()
        assert st['col_step'] == cfg['col_step']
        assert st['h'] == cfg['h']
        assert st['z_tail'] == cfg['z_tail']
        assert st['threads'] == min(cfg['threads'], FS._MAX_THREADS)


def test_apply_accuracy_mode_unknown_raises(fast_settings):
    with pytest.raises(ValueError):
        FS.apply_accuracy_mode('bogus')


def test_engine_fast_forwards_kwargs(monkeypatch, fast_settings):
    captured = {}

    def fake_fast(m, **kw):
        captured.update(kw)
        m.SGWB_converge = True
        return m

    monkeypatch.setattr(FS, 'SGWB_iter_fast', fake_fast)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    r = m.SGWB_iter(engine='fast', accuracy_mode='production', h=0.005,
                    tol=1e-6)
    assert r is m
    assert captured['tol'] == 1e-6
    assert captured['freq_res'] == FS.ACCURACY_MODES['production']['freq_res']
    st = FS.get_settings()
    assert st['h'] == 0.005
    assert st['col_step'] == FS.ACCURACY_MODES['production']['col_step']
    assert st['z_tail'] == FS.ACCURACY_MODES['production']['z_tail']


def test_engine_fast_applies_preset(monkeypatch, fast_settings):
    captured = {}

    def fake_fast(m, **kw):
        captured.update(kw)
        m.SGWB_converge = True
        return m

    monkeypatch.setattr(FS, 'SGWB_iter_fast', fake_fast)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m.SGWB_iter(engine='fast', accuracy_mode='reference')
    cfg = FS.ACCURACY_MODES['reference']
    assert captured['tol'] == cfg['tol']
    assert captured['freq_res'] == cfg['freq_res']
    st = FS.get_settings()
    assert st['h'] == cfg['h']
    assert st['col_step'] == cfg['col_step']
    assert st['z_tail'] == cfg['z_tail']


def test_engine_fast_without_preset_keeps_module_state(monkeypatch,
                                                       fast_settings):
    captured = {}

    def fake_fast(m, **kw):
        captured.update(kw)
        m.SGWB_converge = True
        return m

    FS.set_threads(4)
    FS.set_col_step(4)
    FS.set_h(0.01)
    FS.set_z_tail(5.0)
    monkeypatch.setattr(FS, 'SGWB_iter_fast', fake_fast)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m.SGWB_iter(engine='fast', tol=1e-7)
    assert captured['tol'] == 1e-7
    assert captured['freq_res'] == 1.0
    assert FS.get_settings() == dict(threads=4, col_step=4, h=0.01, z_tail=5.0)


def test_pool_size_env_override(monkeypatch):
    monkeypatch.setenv('SGWB_POOL_SIZE', '2')
    assert _sgwb_pool_size() == 2
    monkeypatch.setenv('SGWB_POOL_SIZE', '0')
    with pytest.raises(ValueError):
        _sgwb_pool_size()


def test_pool_size_mpi_default(monkeypatch):
    monkeypatch.delenv('SGWB_POOL_SIZE', raising=False)
    monkeypatch.setattr('stiffgwpy.stiff_SGWB._mpi_world_size', lambda: 4)
    assert _sgwb_pool_size() == 1
    monkeypatch.setattr('stiffgwpy.stiff_SGWB._mpi_world_size', lambda: 1)
    assert _sgwb_pool_size() == 4
