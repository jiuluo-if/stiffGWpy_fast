"""Contract tests for the standalone no-copy primitive spike."""

import numpy as np

from scripts.benchmark_fast_phi_nocopy import no_copy_primitive
from stiffgwpy_fast import exact_background as EB
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

CASE = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)


def test_formal_kink_grid_is_bitwise_identical():
    model = LCDM_SG(**CASE)
    FS.gen_fast(model, kink_split=True)
    baseline = EB.fast_phi_s2_split(
        model, model.Nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)
    candidate = no_copy_primitive(
        model, model.Nv, model.cosmo_param['DN_eff'], sigma_nodes=model.sigma)
    for expected, actual in zip(baseline[:4], candidate[:4]):
        assert np.array_equal(expected, actual)
    assert baseline[4:] == candidate[4:]


def test_exact_kink_node_preserves_one_sided_branch():
    model = LCDM_SG(**CASE)
    FS.gen_fast(model, kink_split=True)
    n_re = model.derived_param['N_inf'] - model.derived_param['N_re']
    Nv = np.unique(np.sort(np.concatenate((model.Nv, [n_re]))))
    sigma = EB.sigma_vec(Nv, model, model.cosmo_param['DN_eff'])
    baseline = EB.fast_phi_s2_split(
        model, Nv, model.cosmo_param['DN_eff'], sigma_nodes=sigma)
    candidate = no_copy_primitive(
        model, Nv, model.cosmo_param['DN_eff'], sigma_nodes=sigma)
    for expected, actual in zip(baseline[:4], candidate[:4]):
        assert np.array_equal(expected, actual)
    assert baseline[4:] == candidate[4:]
