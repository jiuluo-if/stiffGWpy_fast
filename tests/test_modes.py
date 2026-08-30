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
    assert set(FS.ACCURACY_MODES) == {'debug', 'fast', 'production',
                                      'reference', 'ultra-fast'}
    for cfg in FS.ACCURACY_MODES.values():
        assert 1 <= cfg['col_step'] <= 8
        assert 1e-4 <= cfg['h'] <= 0.1
        assert 2.0 <= cfg['z_tail'] <= 15.0
        assert cfg['freq_res'] >= 1.0
        assert 0 < cfg['tol'] < 1
        assert 1 <= cfg['threads']
    assert FS.ACCURACY_MODES['reference']['h'] < FS.ACCURACY_MODES['production']['h']
    assert FS.ACCURACY_MODES['production']['z_tail'] > FS.ACCURACY_MODES['ultra-fast']['z_tail']
    # The four canonical tiers requested by the audit; ultra-fast is an alias.
    assert FS.ACCURACY_MODES['fast'] == FS.ACCURACY_MODES['ultra-fast']
    assert FS.ACCURACY_MODES['debug']['h'] < FS.ACCURACY_MODES['fast']['h']


def test_error_budget_available(fast_settings):
    for name, cfg in FS.ACCURACY_MODES.items():
        assert name in FS.ERROR_BUDGET
        est = FS.estimate_error(name)
        assert set(est) == {'accuracy_mode', 'DN_gw_error', 'spectrum_error',
                            'quadrature_error', 'integration_error',
                            'ODE_error', 'tail_error', 'model_bias_error'}
        assert est['accuracy_mode'] == name
        # The integration error is the worst per-stage relative error, so it
        # dominates the engine terms and is reported as such.
        assert est['integration_error'] >= max(est['ODE_error'],
                                               est['quadrature_error'],
                                               est['tail_error'])
    with pytest.raises(ValueError):
        FS.estimate_error('bogus')


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


def test_auto_escalate_to_reference_engine(monkeypatch, fast_settings):
    """error-too-large escalates to the continuous-sigma reference engine."""
    import numpy as np

    from stiffgwpy import fast_sgwb as FS
    from stiffgwpy import reference as REF
    from stiffgwpy.stiff_SGWB import LCDM_SG

    def fake_fast(m, **kw):
        m.SGWB_converge = True
        m.cosmo_param['DN_eff'] = 0.002
        return m

    monkeypatch.setattr(FS, 'SGWB_iter_fast', fake_fast)

    def fake_ref(m, **kw):
        m.cosmo_param['DN_eff'] = 0.00227
        m.DN_gw = np.array([0.0, 0.00227])
        m.kappa_r = 1.9e-3
        m.log10OmegaGW = np.array([-15.0])
        m.f = np.array([6.0])
        m.SGWB_converge = True
        m.reference_evals = getattr(m, 'reference_evals', 0) + 1
        return m

    monkeypatch.setattr(REF, 'apply_reference_to_model', fake_ref)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m.SGWB_iter(engine='fast', accuracy_mode='production',
                auto_escalate=True, error_tol=5e-3, escalate_to_reference=True)
    assert m.escalations == 1
    assert m.escalated_from == 'production'
    assert m.reference_evals >= 1
    assert m.cosmo_param['DN_eff'] == pytest.approx(0.00227)
