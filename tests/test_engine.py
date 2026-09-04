# -*- coding: utf-8 -*-
import numpy as np
import pytest

from stiffgwpy_fast import LCDM_SG
from stiffgwpy_fast import fast_sgwb as FS


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


def test_fast_state_isolation_across_parameters_and_threads(fast_settings):
    """A→B→A and thread changes must not reuse parameter-dependent results."""
    m_a1 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m_b = LCDM_SG(r=2e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m_a2 = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.set_threads(1)
    m_a1.SGWB_iter(engine='fast')
    FS.set_threads(min(8, FS._MAX_THREADS))
    m_b.SGWB_iter(engine='fast')
    FS.set_threads(1)
    m_a2.SGWB_iter(engine='fast')
    np.testing.assert_array_equal(m_a1.log10OmegaGW, m_a2.log10OmegaGW)
    assert not np.array_equal(m_a1.log10OmegaGW, m_b.log10OmegaGW)


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


def test_failed_fallback_is_attempted_only_once(monkeypatch):
    """A failed LSODA retry must not recurse and double-count fallbacks."""
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def __init__(self):
            self.fast_failure_reason = None
            self.fast_evals = self.fast_failures = 0
            self.lsoda_fallbacks = self.lsoda_evals = 0

        def reset(self):
            calls.append('reset')

        def SGWB_iter(self, engine='lsoda', fallback=False, **kwargs):
            calls.append(engine)
            if engine == 'lsoda':
                return None
            return self

    f = Fake()

    def boom(m, **kwargs):
        raise RuntimeError('jit failure')

    monkeypatch.setattr(FS, 'SGWB_iter_fast', boom)
    assert wrapper(f, engine='fast', fallback=True) is None
    assert calls == ['reset', 'lsoda']
    assert f.lsoda_fallbacks == 1


def test_engine_telemetry_counts_fast_failures_and_fallback(monkeypatch):
    """Numerical failures are observable even when LSODA recovery succeeds."""
    class Fake:
        def __init__(self):
            self.SGWB_converge = False
            self.fast_evals = self.fast_failures = 0
            self.lsoda_evals = self.lsoda_fallbacks = 0
            self.last_engine = None

        def reset(self):
            self.SGWB_converge = False

        def SGWB_iter(self, engine='lsoda', fallback=False, **kwargs):
            if engine == 'lsoda':
                self.lsoda_evals += 1
                self.last_engine = 'lsoda'
                self.SGWB_converge = True
            return self

    f = Fake()
    monkeypatch.setattr(FS, 'SGWB_iter_fast', lambda m, **kw: None)
    result = LCDM_SG.SGWB_iter(f, engine='fast', fallback=True)
    assert result is f
    assert (f.fast_evals, f.fast_failures, f.lsoda_fallbacks,
            f.lsoda_evals) == (1, 1, 1, 1)
    assert f.last_engine == 'lsoda'


def test_shared_guard_rejection_does_not_retry_lsoda(monkeypatch):
    """Deterministic physical guard failures must not trigger slow fallback."""
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def __init__(self):
            self.fast_failure_reason = None
            self.fast_evals = self.fast_failures = 0
            self.fast_guard_rejections = 0
            self.lsoda_fallbacks = self.lsoda_evals = 0

        def reset(self):
            calls.append('reset')

        def SGWB_iter(self, engine='lsoda', fallback=False, **kwargs):
            calls.append(engine)
            if engine == 'lsoda':
                raise AssertionError('shared guard should not retry LSODA')
            return self

    f = Fake()

    def guard(m, **kwargs):
        m.fast_failure_reason = 'shared_Neff_guard'
        return None

    monkeypatch.setattr(FS, 'SGWB_iter_fast', guard)
    result = wrapper(f, engine='fast', fallback=True)
    assert result is None
    assert calls == []
    assert f.fast_guard_rejections == 1
    assert f.fast_failures == 0


@pytest.mark.parametrize('reason', ['invalid_r', 'invalid_cutoff'])
def test_invalid_physical_rejection_does_not_retry_lsoda(monkeypatch, reason):
    """Deterministic input validation failures are not numerical fallbacks."""
    wrapper = LCDM_SG.SGWB_iter
    calls = []

    class Fake:
        def __init__(self):
            self.fast_failure_reason = None
            self.fast_evals = self.fast_failures = 0
            self.fast_physical_rejections = 0
            self.lsoda_fallbacks = self.lsoda_evals = 0

        def SGWB_iter(self, engine='lsoda', fallback=False, **kwargs):
            calls.append(engine)
            if engine == 'lsoda':
                raise AssertionError('input rejection should not retry LSODA')
            return self

    f = Fake()

    def reject(m, **kwargs):
        m.fast_failure_reason = reason
        return None

    monkeypatch.setattr(FS, 'SGWB_iter_fast', reject)
    assert wrapper(f, engine='fast', fallback=True) is None
    assert calls == []
    assert f.fast_physical_rejections == 1
    assert f.lsoda_fallbacks == 0


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
    mf.SGWB_iter(engine='fast', accuracy_mode='production')
    assert mo.SGWB_converge and mf.SGWB_converge
    # 显式 production 档位仍是精度回归档位；保留其基准积分结果作为数值锚点，
    # 独立 reference/oracle 仍是最终精度依据。
    assert mf.DN_gw[-1] == pytest.approx(0.002261731150563835, rel=2e-5)
    rel = abs(mo.DN_gw[-1] - mf.DN_gw[-1]) / abs(mo.DN_gw[-1])
    # LSODA 仅作为回归引擎，不作为精度真值。
    assert rel < 3e-3             # production 与 LSODA 的实测差异约为 2.1e-3
    # 两个引擎使用不同频率网格；先把 LSODA 谱插值到 fast 节点，再在公共网格上比较，
    # 保持与 fast-reference 认证相同的匹配网格纪律。
    o_f = np.argsort(np.asarray(mf.f))
    f_f = np.asarray(mf.f)[o_f]
    lo_f = np.asarray(mf.log10OmegaGW)[o_f]
    o_l = np.argsort(np.asarray(mo.f))
    f_l = np.asarray(mo.f)[o_l]
    lo_l = np.asarray(mo.log10OmegaGW)[o_l]
    d_pt = np.abs(np.interp(f_f, f_l, lo_l) - lo_f)
    # 单模一致性只在物理上已解析的 logf ∈ [-8, +3] 区间断言；两端分别受不可观测的
    # 超视界冻结底和 UV 截止网格节点敏感性影响，对 Delta N_eff 贡献可以忽略。
    m_band = (f_f >= -8.0) & (f_f <= 3.0)
    dex = float(d_pt[m_band].max())
    assert dex < 1e-2, 'dex=%.3e at logf=%.2f' % (
        dex, f_f[m_band][int(np.argmax(d_pt[m_band]))])
