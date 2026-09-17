import numpy as np

from stiffgwpy_fast import exact_background as EB
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_fast_phi_s2_split_reuses_dead_output_buffers():
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.apply_accuracy_mode('fast')
    FS.gen_fast(model, kink_split=True)
    nv = model.Nv.astype(np.float64)
    sigma = model.sigma.copy()
    dn_eff = model.cosmo_param['DN_eff']

    first = EB.fast_phi_s2_split(model, nv, dn_eff, sigma_nodes=sigma)
    snapshots = tuple(value.copy() for value in first[:4])
    second = EB.fast_phi_s2_split(model, nv, dn_eff, sigma_nodes=sigma)

    for before, after in zip(snapshots, second[:4]):
        np.testing.assert_array_equal(before, after)
    for before, after in zip(first[:4], second[:4]):
        assert before is after
