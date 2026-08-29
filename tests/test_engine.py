# -*- coding: utf-8 -*-
import numpy as np
import pytest

from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb as FS


def test_engine_fast_converges(model):
    r = model.SGWB_iter(engine='fast')
    assert r is model
    assert model.SGWB_converge
    assert model.DN_gw[-1] > 0


def test_engine_fast_deterministic(model):
    m2 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    model.SGWB_iter(engine='fast')
    m2.SGWB_iter(engine='fast')
    assert model.log10OmegaGW.shape == m2.log10OmegaGW.shape
    np.testing.assert_array_equal(model.log10OmegaGW, m2.log10OmegaGW)


def test_engine_invalid_raises(model):
    with pytest.raises(ValueError):
        model.SGWB_iter(engine='bogus')


def test_engine_fast_fallback_reruns_lsoda_on_none(monkeypatch):
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def reset(self):
            calls.append('reset')

        def SGWB_iter(self, engine='lsoda', fallback=False):
            calls.append(engine)
            return self

    f = Fake()
    monkeypatch.setattr(FS, 'SGWB_iter_fast', lambda m, **kw: None)
    result = wrapper(f, engine='fast', fallback=True)
    assert result is f
    assert calls == ['reset', 'lsoda']


def test_engine_fast_no_fallback_returns_none(monkeypatch):
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def SGWB_iter(self, engine='lsoda', fallback=False):
            calls.append(engine)
            return self

    f = Fake()
    monkeypatch.setattr(FS, 'SGWB_iter_fast', lambda m, **kw: None)
    result = wrapper(f, engine='fast', fallback=False)
    assert result is None
    assert calls == []


def test_engine_fast_fallback_on_exception(monkeypatch):
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def reset(self):
            calls.append('reset')

        def SGWB_iter(self, engine='lsoda', fallback=False):
            calls.append(engine)
            return self

    f = Fake()

    def boom(m, **kw):
        raise RuntimeError('jit failure')

    monkeypatch.setattr(FS, 'SGWB_iter_fast', boom)
    result = wrapper(f, engine='fast', fallback=True)
    assert result is f
    assert calls == ['reset', 'lsoda']

    calls.clear()
    with pytest.raises(RuntimeError):
        wrapper(f, engine='fast', fallback=False)


@pytest.mark.slow
def test_engine_lsoda_returns_model(model):
    r = model.SGWB_iter(engine='lsoda')
    assert r is model
    assert model.SGWB_converge


@pytest.mark.slow
def test_engine_fast_close_to_lsoda_on_case_a():
    kw = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    mo = LCDM_SG(**kw)
    mo.SGWB_iter(engine='lsoda')
    mf = LCDM_SG(**kw)
    mf.SGWB_iter(engine='fast')
    assert mo.SGWB_converge and mf.SGWB_converge
    rel = abs(mo.DN_gw[-1] - mf.DN_gw[-1]) / abs(mo.DN_gw[-1])
    assert rel < 1e-3             # measured ~7e-5
    dex = np.abs(np.asarray(mo.log10OmegaGW) - np.asarray(mf.log10OmegaGW)).max()
    assert dex < 1e-2             # measured ~4e-4
