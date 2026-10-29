"""Tests for domain-pressure classifications at exact limits."""

from scripts.pressure_fast_domain import evaluate


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
