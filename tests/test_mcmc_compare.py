# -*- coding: utf-8 -*-
import importlib.util
import os
import sys

import numpy as np
import pytest


def _load_mcmc_compare():
    scripts = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'scripts')
    path = os.path.join(scripts, 'mcmc_compare.py')
    spec = importlib.util.spec_from_file_location('mcmc_compare', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['mcmc_compare'] = mod
    spec.loader.exec_module(mod)
    return mod


MC = _load_mcmc_compare()


def _chain(n=1000, seed=0):
    rng = np.random.default_rng(seed)
    sample = np.column_stack([np.ones(n),                    # weight
                              rng.random(n) + 2.0,           # minuslogpost
                              rng.normal(0.0, 1.0, n),       # r
                              rng.normal(5.0, 2.0, n),       # n_t
                              rng.random(n),                 # minuslogprior
                              rng.random(n) + 1.0])          # minusloglike
    sample_params = ['weight', 'minuslogpost', 'r', 'n_t',
                     'minuslogprior', 'minusloglike']
    sample[0, 1] = np.inf                                # one failed sample
    return sample, sample_params


def test_split_sample_columns():
    params, specials = MC.split_sample_columns(
        ['weight', 'minuslogpost', 'r', 'n_t',
         'minuslogprior', 'minusloglike'])
    assert params == ['r', 'n_t']
    assert specials == {'weight': 0, 'minuslogpost': 1,
                        'minuslogprior': 4, 'minusloglike': 5}


def test_split_sample_columns_chi2_and_param_indices():
    cols = ['weight', 'minuslogpost', 'log10r', 'h', 'minuslogprior',
            'minuslogprior__ext', 'chi2', 'chi2__likA', 'chi2__likB']
    params, specials = MC.split_sample_columns(cols)
    assert params == ['log10r', 'h']
    assert specials == {'weight': 0, 'minuslogpost': 1,
                        'minuslogprior': 4, 'minuslogprior__ext': 5,
                        'chi2': 6, 'chi2__likA': 7, 'chi2__likB': 8}
    assert MC.param_indices(cols) == {'log10r': 2, 'h': 3}


def test_chain_stats_matches_numpy():
    sample, sample_params = _chain()
    st = MC.chain_stats(sample, sample_params)
    assert st['n_total'] == 1000
    assert st['n_finite'] == 999
    assert st['failure_rate'] == pytest.approx(1.0 - 999 / 1000)
    ok = np.isfinite(sample[:, 1])
    assert st['stats']['r']['mean'] == np.mean(sample[ok, 2])
    assert st['stats']['n_t']['std'] == np.std(sample[ok, 3])
    i_map = int(np.argmin(sample[ok, 1]))
    assert st['map']['r'] == sample[ok, 2][i_map]
    assert st['map']['n_t'] == sample[ok, 3][i_map]
    assert st['map']['minuslogpost'] == sample[ok, 1][i_map]
    assert st['min_ess'] > 0
    assert len(st['covariance']) == len(st['covariance_params'])


def test_effective_sample_size_constant_and_correlated():
    assert MC.effective_sample_size(np.ones(20)) == 20.0
    rng = np.random.default_rng(5)
    x = np.cumsum(rng.normal(size=200))
    assert 1.0 <= MC.effective_sample_size(x) < 200.0


def test_chain_stats_honors_cobaya_weights():
    sample = np.array([[1.0, 0.0, 0.0], [3.0, 0.0, 10.0]])
    names = ['weight', 'minuslogpost', 'x']
    st = MC.chain_stats(sample, names)
    assert st['stats']['x']['mean'] == pytest.approx(7.5)
    assert st['stats']['x']['median'] == pytest.approx(10.0)


def test_posterior_shift():
    sample, sample_params = _chain()
    ls = MC.chain_stats(sample, sample_params)
    fa_sample = sample.copy()
    ok = np.isfinite(sample[:, 1])
    std_r = np.std(sample[ok, 2])
    fa_sample[ok, 2] += 0.5 * std_r
    fa = MC.chain_stats(fa_sample, sample_params)
    shift = MC.posterior_shift(ls, fa)
    assert np.isclose(shift['r'], 0.5, atol=1e-9)
    assert np.isclose(shift['n_t'], 0.0, atol=1e-9)


def test_max_abs_posterior_shift_handles_missing_and_nonfinite():
    assert MC.max_abs_posterior_shift({'a': 0.05, 'b': -0.2}) == 0.2
    assert MC.max_abs_posterior_shift({'a': np.nan}) == float('inf')
    assert MC.max_abs_posterior_shift({}) == float('inf')


def test_certification_shift_threshold_defaults_to_point_one():
    args = MC.parse_args([])
    assert args.max_posterior_shift == 0.1


def test_delta_logl_stats():
    d = np.array([0.0, 0.1, -0.2, 1.2, -3.0, 2.5])
    st = MC.delta_logl_stats(d)
    assert st['n'] == 6
    assert st['max'] == 2.5
    assert st['max_abs'] == 3.0
    assert st['median'] == np.median(d)
    assert st['frac_abs_gt_1.0'] == 3 / 6
    assert MC.delta_logl_stats(np.array([np.nan, np.inf])) == {'n': 0}


def test_posterior_distances_finite_and_deterministic():
    sample, sample_params = _chain(200, seed=4)
    d = MC.posterior_distances(sample, sample.copy(), sample_params)
    assert d['r']['ks_stat'] == 0.0
    assert d['r']['wasserstein'] == 0.0
    assert d['r']['kl_hist_32'] == 0.0


def test_constant_parameter_ess_does_not_block_certification():
    sample = np.array([
        [1.0, 0.0, 0.0, 3.0],
        [1.0, 0.0, 0.1, 3.0],
        [1.0, 0.0, 0.2, 3.0],
        [1.0, 0.0, 0.3, 3.0],
    ])
    out = MC.chain_stats(sample, ['weight', 'minuslogpost', 'x', 'constant'])
    assert out['stats']['constant']['ess'] == 4.0
    assert out['stats']['x']['ess'] >= 1.0


def test_thin_idx():
    assert np.array_equal(MC.thin_idx(10, 20), np.arange(10))
    idx = MC.thin_idx(100, 5)
    assert len(idx) == 5
    assert idx[0] == 0 and idx[-1] == 99
    assert np.all(np.diff(idx) >= 1)


def test_build_info_overrides_engine(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW:\n    engine: lsoda\n'
        'params:\n  r: {value: 0.01}\n', encoding='utf-8')
    info = MC.build_info(str(yaml_path), 'fast', 'production', 8, 500, 42)
    theory = info['theory']['stiffgwpy_fast.cobaya.stiffGW.stiffGW']
    assert theory['engine'] == 'fast'
    assert theory['accuracy_mode'] == 'production'
    assert theory['fast_threads'] == 8
    assert info['sampler']['mcmc']['max_samples'] == 500
    assert info['sampler']['mcmc']['seed'] == 42


def test_build_info_normalizes_qualified_theory_name(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW:\n    python_path: .\n'
        'params:\n  r: {value: 0.01}\n', encoding='utf-8')
    info = MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1)
    assert list(info['theory']) == ['stiffgwpy_fast.cobaya.stiffGW.stiffGW']
    assert 'python_path' not in info['theory']['stiffgwpy_fast.cobaya.stiffGW.stiffGW']


def test_build_info_applies_explicit_initial_point_to_both_runs(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW: {}\n'
        'params:\n  x:\n    prior: {min: -1, max: 1}\n', encoding='utf-8')
    point = {'x': 0.25}
    info_fast = MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1,
                              initial_point=point)
    info_lsoda = MC.build_info(str(yaml_path), 'lsoda', '', 0, 5, 1,
                               initial_point=point)
    assert info_fast['params']['x']['ref'] == 0.25
    assert info_lsoda['params']['x']['ref'] == 0.25


def test_build_info_rejects_undefined_initial_point(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW: {}\nparams: {}\n', encoding='utf-8')
    with pytest.raises(ValueError, match='not defined'):
        MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1,
                      initial_point={'A_BBH': -16.0})


def test_build_info_common_proposal_scale_override(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW: {}\nparams:\n  x:\n    prior: {min: -1, max: 1}\n',
        encoding='utf-8')
    info = MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1,
                         proposal_scale=0.7)
    assert info['sampler']['mcmc']['proposal_scale'] == 0.7
    with pytest.raises(ValueError, match='positive finite'):
        MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1,
                      proposal_scale=0)


def test_build_info_rejects_derived_initial_point(tmp_path):
    yaml_path = tmp_path / 'run.yaml'
    yaml_path.write_text(
        'theory:\n  stiffGW: {}\n'
        'params:\n  x:\n    derived: "lambda: 1"\n', encoding='utf-8')
    with pytest.raises(ValueError, match='sampled/free'):
        MC.build_info(str(yaml_path), 'fast', 'production', 8, 5, 1,
                      initial_point={'x': 0.25})


def test_yaml_hash_is_content_addressed(tmp_path):
    path = tmp_path / 'run.yaml'
    path.write_text('params: {}\n', encoding='utf-8')
    first = MC._yaml_sha256(path)
    path.write_text('params: {r: {value: 0.01}}\n', encoding='utf-8')
    second = MC._yaml_sha256(path)
    assert first != second
    assert len(first) == 64 and len(second) == 64
