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


def test_vectorized_pchip_integrals_match_scipy_breakdown():
    """向量化封闭式积分核必须与 scipy 参考实现一致。"""
    for seed, n in ((20260921, 5), (20260922, 17), (20260923, 76), (20260924, 2)):
        freqs, values = _random_case(seed, n)
        values = values * np.exp(values)
        total, local = FS.pchip_integral_breakdown(freqs, values)
        fast = FS._pchip_integrals_vectorized(freqs, values)
        scale = float(np.max(np.abs(local)))
        assert fast.shape == local.shape
        assert np.allclose(fast, local, rtol=1e-9, atol=1e-12 * scale)
        assert float(np.sum(fast)) == pytest.approx(total, rel=1e-12)


def test_vectorized_pchip_integrals_handle_flat_and_sign_changing_data():
    """零斜率与符号翻转（Fritsch-Carlson 置零分支）也必须与 scipy 一致。"""
    cases = (
        np.array([1.0, 1.0, 1.0, 1.0]),
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([1.0, -1.0, 1.0, -1.0]),
        np.array([2.0, 2.0, 0.0, -2.0]),
        np.array([3.0, 5.0]),
    )
    for values in cases:
        freqs = np.arange(values.size, dtype=float)
        total, local = FS.pchip_integral_breakdown(freqs, values)
        fast = FS._pchip_integrals_vectorized(freqs, values)
        assert np.allclose(fast, local, rtol=1e-12, atol=1e-14)
        assert float(np.sum(fast)) == pytest.approx(total, rel=1e-12)


def test_vectorized_pchip_integrals_reject_non_unique_nodes():
    with pytest.raises(ValueError, match='unique'):
        FS._pchip_integrals_vectorized(np.array([0.0, 1.0, 1.0]), np.ones(3))


def _panel_loop_reference(candidate, baseline, allocation):
    """旧逐面板 Python 循环的独立参考副本。"""
    local_error = np.zeros_like(candidate)
    start = 0
    while start + 1 < candidate.size:
        end = start + 2
        panel_error = abs(float(np.sum(candidate[start:end]))
                          - float(np.sum(baseline[start:end])))
        if allocation in ('full_panel', 'panel_envelope'):
            local_error[start:end] = panel_error
        else:
            share = np.abs(candidate[start:end])
            share_sum = float(np.sum(share))
            if share_sum > 0.0:
                local_error[start:end] = panel_error * share / share_sum
            else:
                local_error[start:end] = 0.5 * panel_error
        start += 2
    if start < candidate.size:
        local_error[start] = abs(candidate[start] - baseline[start])
    if allocation == 'panel_envelope':
        panel_error = local_error.copy()
        for i in range(local_error.size):
            local_error[i] = np.max(
                panel_error[max(0, i - 2):min(local_error.size, i + 3)])
    return local_error


def test_vectorized_allocation_matches_panel_loop():
    """向量化分摊必须与逐面板循环逐位一致（三种 allocation）。"""
    for seed, n in ((20260931, 7), (20260932, 26), (20260933, 76)):
        freqs, values = _random_case(seed, n)
        values = values * np.exp(values) * 1e-3
        for allocation in ('weighted', 'full_panel', 'panel_envelope'):
            errors, candidate, baseline = FS.estimate_frequency_quadrature_local(
                freqs, values, 'pchip', allocation=allocation)
            expected = _panel_loop_reference(candidate, baseline, allocation)
            assert np.array_equal(errors, expected), (seed, n, allocation)
