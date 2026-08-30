# -*- coding: utf-8 -*-
"""Regression tests pinning the fast kernels against the original/scipy.

These reproduce the component-level agreement reported by the independent
audit (gen_fast vs gen_expansion ~1e-11, Simpson weights ~1e-12, PCHIP ~1e-15).
"""
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator

from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb as FS


def _random_model(rng):
    kw = dict(r=10 ** rng.uniform(-4, -1),
              cr=int(rng.integers(0, 2)),
              T_re=10 ** rng.uniform(1, 4),
              kappa10=10 ** rng.uniform(-3, 0),
              DN_eff=rng.uniform(0, 0.5))
    m = LCDM_SG(**kw)
    return None if m.derived_param['N_inf'] is None else m


def test_gen_fast_matches_original_expansion():
    rng = np.random.default_rng(11)
    worst_sigma, worst_fhor, checked = 0.0, 0.0, 0
    for _ in range(8):
        m = _random_model(rng)
        if m is None:
            continue
        m.gen_expansion()
        sig0, fh0 = m.sigma.copy(), m.f_hor.copy()
        m2 = LCDM_SG(r=m.cosmo_param['r'], cr=m.cosmo_param['cr'],
                     T_re=m.cosmo_param['T_re'], kappa10=m.cosmo_param['kappa10'],
                     DN_eff=m.cosmo_param['DN_eff'])
        FS.gen_fast(m2)
        worst_sigma = max(worst_sigma, float(np.abs(m2.sigma - sig0).max()))
        worst_fhor = max(worst_fhor, float(np.abs(m2.f_hor - fh0).max()))
        checked += 1
    assert checked > 0
    assert worst_sigma < 1e-10    # audit measured ~1.7e-11
    # Both generators are now anchored at the continuous N_inf (present), so the
    # f_hor difference is float-accumulation noise (~1e-10), not a physics gap.
    assert worst_fhor < 1e-9


def test_simpson_row_matches_scipy():
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(6):
        nf = int(rng.integers(8, 40))
        xf = np.sort(rng.uniform(-2, 2, nf))
        hf = np.diff(xf)
        y = rng.normal(size=nf)
        for jh in range(nf):
            a = nf - 1 - jh
            w = np.zeros(nf)
            FS.simpson_row(xf, hf, a, nf - 1, w)
            ref = simpson(y[a:], x=xf[a:])
            mine = np.dot(w[a:], y[a:])
            worst = max(worst, abs(mine - ref) / max(abs(ref), 1e-300))
    assert worst < 1e-9           # audit measured ~1e-12


def test_pchip_fine_matches_scipy():
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(3):
        nv = 400
        idx = np.unique(np.append(np.arange(0, nv, 7), nv - 1))
        y = np.exp(rng.normal(size=len(idx)))
        out = np.empty(nv)
        FS.pchip_fine(idx.astype(np.float64), y, nv, out)
        ref = PchipInterpolator(idx, y)(np.arange(nv))
        worst = max(worst, float(np.abs(out - ref).max()))
    assert worst < 1e-12          # audit measured ~5e-15
