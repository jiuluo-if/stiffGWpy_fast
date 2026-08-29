# -*- coding: utf-8 -*-
import numpy as np
import pytest

from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb as FS


def _make_model(**kw):
    defaults = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    defaults.update(kw)
    return LCDM_SG(**defaults)


def test_set_col_step_validates():
    FS.set_col_step(4)
    for bad in (0, 9, -1, 'x'):
        with pytest.raises(ValueError):
            FS.set_col_step(bad)
    FS.set_col_step(1)
    FS.set_col_step(8)


def test_set_threads_validates():
    FS.set_threads(1)
    for bad in (0, -2, 'x'):
        with pytest.raises(ValueError):
            FS.set_threads(bad)
    with pytest.raises(ValueError):
        FS.set_threads(FS._MAX_THREADS + 1)
    FS.set_threads(FS._MAX_THREADS)


def test_r_le_zero_returns_none(model):
    m = _make_model(r=0.0)
    assert FS.SGWB_iter_fast(m) is None


def test_nonfinite_dn_gw_aborts_and_restores(monkeypatch, model):
    orig_dn = model.cosmo_param['DN_eff']

    def bad_solve(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
                  fp_minus, fp_freq, assemble, n_coarse, col_step, Ogw, Oj, Opgw):
        Ogw[...] = np.nan
        Oj[...] = 0.0
        Opgw[...] = 0.0

    monkeypatch.setattr(FS, 'solve_kernel', bad_solve)
    assert FS.SGWB_iter_fast(model) is None
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.DN_eff_orig is None
    assert model.SGWB_converge is False


def test_max_iterations_abort_and_restore(monkeypatch, model):
    orig_dn = model.cosmo_param['DN_eff']
    rng = np.random.default_rng(42)
    monkeypatch.setattr(FS, 'MAX_ITER', 6)

    def fake_solve(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
                   fp_minus, fp_freq, assemble, n_coarse, col_step, Ogw, Oj, Opgw):
        # Produce a random DN_gw_new in (0.05, 4.9) that never converges.
        Omega_nu = FS.gp.Omega_nh2 / model.derived_param['h'] ** 2
        span = np.log10(fp_freq.max()) - np.log10(fp_freq.min())
        val = rng.uniform(0.05, 4.9)
        const = val * Omega_nu / (FS.gp.Neff0 * FS.ln10 * span)
        Ogw[...] = 0.0
        Oj[...] = 0.0
        Opgw[...] = 0.0
        Ogw[:, -1] = const

    monkeypatch.setattr(FS, 'solve_kernel', fake_solve)
    assert FS.SGWB_iter_fast(model) is None
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.DN_eff_orig is None
    assert model.SGWB_converge is False


def test_stale_full_evolution_attrs_cleared(model):
    for name in ('N_hc', 'Th', 'Oj', 'Ogw', 'Opgw'):
        setattr(model, name, [np.array([1.0])])
    assert FS.SGWB_iter_fast(model) is model
    assert model.SGWB_converge
    for name in ('N_hc', 'Th', 'Oj', 'Ogw', 'Opgw'):
        assert not hasattr(model, name), name


def test_success_leaves_dn_eff_orig_for_reset(model):
    orig_dn = model.cosmo_param['DN_eff']
    FS.SGWB_iter_fast(model)
    assert model.SGWB_converge
    assert model.DN_eff_orig == orig_dn
    final_dn = model.cosmo_param['DN_eff']
    assert final_dn > orig_dn
    model.reset()
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.SGWB_converge is False
