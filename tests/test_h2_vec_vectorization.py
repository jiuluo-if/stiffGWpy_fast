import math

import numpy as np

from stiffgwpy_fast import exact_background as EB
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def _scalar_h2_reference(nodes, model, dn_eff):
    d = model.derived_param
    omh2 = d['Omega_mh2']
    osh2 = d['Omega_sh2']
    oerh2 = gp.Omega_ph2 * 7 / 8 * (4 / 11) ** (4 / 3) * dn_eff
    otrh2 = gp.Omega_orh2 + oerh2
    otreh2 = gp.Omega_ph2 * gp.rho_th[-1] + oerh2
    olh2 = (d['h'] ** 2 - omh2 - gp.Omega_mnuh2 - gp.Omega_nh2 * 2 / 3
            - gp.Omega_ph2 - oerh2 - osh2)
    fd_rho, _ = EB._fd_from_ref()
    result = []
    for n in nodes:
        en = math.exp(d['N_inf'] - n)
        e3n = en * en * en
        nu = gp.nu_today / en
        if nu > 100:
            h2 = (omh2 + gp.Omega_mnuh2
                  + (gp.Omega_ph2 + 2 / 3 * gp.Omega_nh2 + oerh2) * en
                  + osh2 * e3n + olh2 / e3n)
        elif nu >= 0.1:
            rho_nu = float(fd_rho(math.log10(nu)))
            h2 = (omh2
                  + (gp.Omega_ph2 + (2 / 3 + rho_nu / 3) * gp.Omega_nh2
                     + oerh2) * en
                  + osh2 * e3n + olh2 / e3n)
        elif n > d['N_inf'] - gp.N_fin:
            h2 = omh2 + otrh2 * en + osh2 * e3n + olh2 / e3n
        elif n >= d['N_inf'] - gp.N_max:
            n_local = d['N_inf'] - n
            h2 = (omh2 + (gp.Omega_ph2 * gp.spl_rho(n_local) + oerh2) * en
                  + osh2 * e3n + olh2 / e3n)
        else:
            h2 = omh2 + otreh2 * en + osh2 * e3n + olh2 / e3n
        result.append(h2)
    return np.asarray(result)


def test_h2_vec_batches_spline_evaluation_without_changing_values(monkeypatch):
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    dn_eff = 0.00237
    n = np.linspace(0.0, model.derived_param['N_inf'], 257)
    expected = _scalar_h2_reference(n, model, dn_eff)
    original_spline = gp.spl_rho
    calls = []

    def counted_spline(values):
        calls.append(np.asarray(values).ndim)
        return original_spline(values)

    monkeypatch.setattr(gp, 'spl_rho', counted_spline)
    actual = EB.H2_vec(n, model, dn_eff)

    # Stable closure grouping changes only the final rounding order.
    np.testing.assert_allclose(actual, expected, rtol=5e-15, atol=0.0)
    assert calls == [1]
