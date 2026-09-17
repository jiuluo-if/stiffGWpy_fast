from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_fast_derived_param_cache_reuses_only_unchanged_cosmology():
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    model._fast_derived_cache_enabled = True

    first = model.derived_param
    first_rhorad = model.rhorad_re
    second = model.derived_param
    assert first is second

    model.cosmo_param['DN_eff'] += 1e-6
    changed = model.derived_param
    assert changed is not second
    assert model.rhorad_re != first_rhorad
