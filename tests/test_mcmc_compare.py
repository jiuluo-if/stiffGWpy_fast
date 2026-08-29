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
    sample = np.column_stack([rng.normal(0.0, 1.0, n),
                              rng.normal(5.0, 2.0, n),
                              rng.random(n),                 # minuslogprior
                              rng.random(n) + 1.0,           # minusloglike
                              rng.random(n) + 2.0,           # minuslogpost
                              np.ones(n)])                   # weight
    sample_params = ['r', 'n_t', 'minuslogprior', 'minusloglike',
                     'minuslogpost', 'weight']
    sample[0] = np.inf                                    # one failed sample
    return sample, sample_params


def test_split_sample_columns():
    params, specials = MC.split_sample_columns(
        ['r', 'n_t', 'minuslogprior', 'minusloglike', 'minuslogpost', 'weight'])
    assert params == ['r', 'n_t']
    assert specials == {'minuslogprior': 2, 'minusloglike': 3,
                        'minuslogpost': 4, 'weight': 5}


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
    ok = np.isfinite(sample[:, 4])
    assert st['stats']['r']['mean'] == np.mean(sample[ok, 0])
    assert st['stats']['n_t']['std'] == np.std(sample[ok, 1])
    i_map = int(np.argmin(sample[ok, 4]))
    assert st['map']['r'] == sample[ok, 0][i_map]
    assert st['map']['minuslogpost'] == sample[ok, 4][i_map]


def test_posterior_shift():
    sample, sample_params = _chain()
    ls = MC.chain_stats(sample, sample_params)
    fa_sample = sample.copy()
    std_r = np.std(sample[np.isfinite(sample[:, 4]), 0])
    fa_sample[:, 0] += 0.5 * std_r
    fa = MC.chain_stats(fa_sample, sample_params)
    shift = MC.posterior_shift(ls, fa)
    assert np.isclose(shift['r'], 0.5, atol=1e-9)
    assert np.isclose(shift['n_t'], 0.0, atol=1e-9)


def test_delta_logl_stats():
    d = np.array([0.0, 0.1, -0.2, 1.2, -3.0, 2.5])
    st = MC.delta_logl_stats(d)
    assert st['n'] == 6
    assert st['max'] == 2.5
    assert st['max_abs'] == 3.0
    assert st['median'] == np.median(d)
    assert st['frac_abs_gt_1.0'] == 3 / 6
    assert MC.delta_logl_stats(np.array([np.nan, np.inf])) == {'n': 0}


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
    assert info['theory']['stiffGW']['engine'] == 'fast'
    assert info['theory']['stiffGW']['accuracy_mode'] == 'production'
    assert info['theory']['stiffGW']['fast_threads'] == 8
    assert info['sampler']['mcmc']['max_samples'] == 500
    assert info['sampler']['mcmc']['seed'] == 42
