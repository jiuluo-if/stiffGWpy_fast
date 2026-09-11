"""PCHIP 单次拟合共享：全程积分、逐区间积分与局部 estimator 的复用路径。"""
import numpy as np
import pytest
from scipy import interpolate

from stiffgwpy_fast import fast_sgwb as FS


def _random_case(seed, n):
    """构造可复现的随机频率/被积函数样本。"""
    rng = np.random.default_rng(seed)
    freqs = np.sort(rng.uniform(-20.0, 5.0, n))
    values = rng.normal(size=n)
    return freqs, values


def test_pchip_integral_breakdown_matches_scipy_global_integral():
    """全程积分必须等于 scipy 的 PCHIP 全区间积分，并与单点辅助函数一致。"""
    freqs, values = _random_case(20260911, 40)
    total, _ = FS.pchip_integral_breakdown(freqs, values)
    expected = float(interpolate.PchipInterpolator(freqs, values).integrate(
        freqs[0], freqs[-1]))
    assert total == pytest.approx(expected, rel=1e-12)
    assert total == pytest.approx(
        FS.integrate_frequency_pchip(freqs, values), rel=1e-12)


def test_pchip_integral_breakdown_matches_scipy_interval_integrals():
    """逐区间数组必须等于同一拟合的逐区间 scipy 积分，且求和等于全程积分。"""
    freqs, values = _random_case(20260912, 25)
    total, local = FS.pchip_integral_breakdown(freqs, values)
    spline = interpolate.PchipInterpolator(freqs, values)
    expected = np.array([float(spline.integrate(left, right))
                         for left, right in zip(freqs[:-1], freqs[1:])])
    assert local.shape == (freqs.size - 1,)
    assert np.allclose(local, expected, rtol=1e-12, atol=1e-14)
    assert float(np.sum(local)) == pytest.approx(total, rel=1e-12)


def test_pchip_integral_breakdown_sorts_descending_nodes():
    """求解器传入降序频率；排序后必须与升序输入给出完全一致的结果。"""
    freqs, values = _random_case(20260913, 18)
    asc_total, asc_local = FS.pchip_integral_breakdown(freqs, values)
    desc_total, desc_local = FS.pchip_integral_breakdown(
        freqs[::-1], values[::-1])
    assert desc_total == asc_total
    assert np.array_equal(desc_local, asc_local)


def test_local_estimator_reuses_precomputed_candidate_intervals():
    """提供预计算逐区间积分时，局部估计器三路输出必须与标量 PCHIP 路径一致。"""
    freqs = np.array([2.0, 1.6, 1.0, 0.55, 0.2, 0.0])
    integrand = np.exp(-freqs) * (1.0 + 0.2 * np.sin(3.0 * freqs))
    _, local = FS.pchip_integral_breakdown(freqs, integrand)
    reference = FS.estimate_frequency_quadrature_local(freqs, integrand, 'pchip')
    reused = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'pchip', candidate_intervals=local)
    for got, want in zip(reused, reference):
        assert np.allclose(got, want, rtol=1e-12, atol=1e-16)


def test_pchip_integral_breakdown_rejects_non_unique_nodes():
    with pytest.raises(ValueError, match='unique'):
        FS.pchip_integral_breakdown(np.array([0.0, 1.0, 1.0]), np.ones(3))


def test_pchip_integral_breakdown_rejects_nonfinite_values():
    with pytest.raises(ValueError, match='finite'):
        FS.pchip_integral_breakdown(
            np.array([0.0, 1.0]), np.array([1.0, np.nan]))


def test_local_estimator_rejects_mismatched_candidate_intervals():
    """预计算数组长度与区间数不符时必须显式报错，不能静默复用。"""
    freqs = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError, match='interval count'):
        FS.estimate_frequency_quadrature_local(
            freqs, np.ones(3), 'pchip', candidate_intervals=np.zeros(5))
