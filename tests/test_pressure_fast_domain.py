"""Tests for domain-pressure classifications at exact limits."""

from scripts.pressure_fast_domain import (
    BASE,
    boundary_cases,
    coupled_cases,
    evaluate,
    qmc_cases,
)


def test_zero_tensor_limit_is_success_with_separate_classification():
    params = {
        'Omega_bh2': 0.0223828,
        'Omega_ch2': 0.1201075,
        'H0': 67.32117,
        'DN_eff': 0.0,
        'A_s': 2.100549e-9,
        'r': 0.0,
        'n_t': 0.5,
        'cr': 0.0,
        'T_re': 2000.0,
        'DN_re': 10.0,
        'kappa10': 0.01,
    }

    result = evaluate({'label': 'zero_tensor', 'params': params},
                      max_n_inf=120.0, threads=2)

    assert result['status'] == 'success'
    assert result['classification'] == 'exact_zero_tensor_spectrum'
    assert result['DN_gw'] == 0.0
    assert result['spectrum_finite_bins'] == 0
    assert result['fast_accuracy_status'] == 'analytic_exact_limit'
    assert result['fast_accuracy_budget'] is None


def test_underflowed_nonzero_tensor_amplitude_is_not_zero_success():
    params = {
        'Omega_bh2': 0.0223828,
        'Omega_ch2': 0.1201075,
        'H0': 67.32117,
        'DN_eff': 0.0,
        'A_s': 1e-60,
        'r': 1e-300,
        'n_t': 0.0,
        'cr': 0.0,
        'T_re': 2000.0,
        'DN_re': 10.0,
        'kappa10': 0.01,
    }

    result = evaluate({'label': 'underflow_tensor', 'params': params},
                      max_n_inf=120.0, threads=2)

    assert result['status'] == 'failed'
    assert result['classification'] == 'floating_point_underflow'
    assert result['reason'] == 'tensor_amplitude_underflow'


def test_finite_success_includes_local_accuracy_budget_without_claiming_certification():
    params = {
        'Omega_bh2': 0.0223828,
        'Omega_ch2': 0.1201075,
        'H0': 67.32117,
        'DN_eff': 0.0,
        'A_s': 2.100549e-9,
        'r': 1e-2,
        'n_t': 0.5,
        'cr': 0.0,
        'T_re': 2000.0,
        'DN_re': 10.0,
        'kappa10': 1e-2,
    }

    result = evaluate({'label': 'finite_with_budget', 'params': params},
                      max_n_inf=120.0, threads=2)

    assert result['status'] == 'success'
    assert result['classification'] == 'finite_positive_spectrum'
    assert result['fast_accuracy_status'] in {
        'uncertified-path-configuration', 'certified-fiducial-calibrated'}
    assert result['fast_accuracy_budget']['DN_gw_error'] >= 0.0
    assert isinstance(result['fast_accuracy_budget']['uncertified_terms'], list)


def test_coupled_suite_skips_redundant_n_t_values_when_consistency_relation_is_active():
    cases = coupled_cases()

    assert len(cases) == 140
    assert {case['params']['cr'] for case in cases} == {0.0, 1.0}
    assert {case['params']['n_t'] for case in cases if case['params']['cr'] == 0.0} == {
        -2.0, 0.0, 1.0, 2.0, 5.0, 10.0}
    assert {case['params']['n_t'] for case in cases if case['params']['cr'] > 0.0} == {
        BASE['n_t']}
    assert {case['params']['kappa10'] for case in cases} == {1e-8, 1e-2, 1.0, 1e2, 1e4}


def test_n_t_boundary_suite_uses_one_canonical_value_when_cr_is_positive():
    cases = [case for case in boundary_cases()
             if case['label'].startswith('boundary_n_t_')]

    assert len([case for case in cases if case['params']['cr'] == 0.0]) == 13
    cr_positive = [case for case in cases if case['params']['cr'] > 0.0]
    assert len(cr_positive) == 1
    assert cr_positive[0]['params']['n_t'] == BASE['n_t']


def test_qmc_canonicalizes_unused_n_t_in_consistency_relation_branch():
    import numpy as np

    points = np.full((4, 10), 0.5)
    points[:, 1] = (0.0, 0.2, 0.8, 1.0)
    cases = qmc_cases(points, 'qmc_test')

    assert [case['params']['cr'] for case in cases] == [0.0, 1.0, 0.0, 1.0]
    assert [case['params']['n_t'] for case in cases if case['params']['cr'] == 0.0] == [-10.0, 6.0]
    assert {case['params']['n_t'] for case in cases if case['params']['cr'] > 0.0} == {BASE['n_t']}


def test_failure_minimizer_skips_n_t_when_consistency_relation_is_active(monkeypatch):
    from scripts import pressure_fast_domain as pressure

    labels = []

    def failed_trial(item, max_n_inf, threads):
        labels.append(item['label'])
        return {'status': 'failed', 'classification': 'nonfinite_numerical_state'}

    monkeypatch.setattr(pressure, 'evaluate', failed_trial)
    item = pressure._case('cr1_failure', cr=1.0, r=0.02)
    pressure.minimize_failure(
        item, {'classification': 'nonfinite_numerical_state'},
        max_n_inf=120.0, threads=1, steps=1,
    )

    assert not any(label.endswith('_n_t') for label in labels)
