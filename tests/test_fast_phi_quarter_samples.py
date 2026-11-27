"""Check quarter-point primitives against the same refined Simpson grid."""

import numpy as np

from stiffgwpy_fast import exact_background
from stiffgwpy_fast.exact_background import (
    fast_phi_s2_quarters,
    fast_phi_s2_split,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_quarter_primitives_match_refined_local_integrals():
    model = LCDM_SG(
        Omega_bh2=0.0223828, Omega_ch2=0.1201075, H0=67.32117,
        DN_eff=0.0, A_s=2.100549e-9, r=0.01, n_t=0.0, cr=1.0,
        T_re=2000.0, DN_re=10.0, kappa10=0.01,
    )
    n_re = model.derived_param["N_inf"] - model.derived_param["N_re"]
    nodes = n_re + 0.1 + np.arange(5, dtype=np.float64) * 0.005
    dn_eff = 0.0
    phi, _, _, _, kink_index, kink_fraction, _ = fast_phi_s2_split(
        model, nodes, dn_eff
    )
    quarter_phi = fast_phi_s2_quarters(
        model, nodes, dn_eff, phi_grid=phi,
        kink_index=kink_index, kink_fraction=kink_fraction,
    )
    fine = np.sort(np.concatenate((
        nodes, nodes[:-1] + 0.25 * np.diff(nodes),
        nodes[:-1] + 0.5 * np.diff(nodes),
        nodes[:-1] + 0.75 * np.diff(nodes),
    )))
    phi_fine = fast_phi_s2_split(model, fine, dn_eff)[0]
    for k in range(len(nodes) - 1):
        if k == kink_index and 0.0 < kink_fraction < 1.0:
            continue
        base_i = int(np.searchsorted(fine, nodes[k]))
        q1_i = int(np.searchsorted(fine, nodes[k] + 0.25 * (nodes[k+1]-nodes[k])))
        q3_i = int(np.searchsorted(fine, nodes[k] + 0.75 * (nodes[k+1]-nodes[k])))
        np.testing.assert_allclose(
            quarter_phi[0, k] - phi[k], phi_fine[q1_i] - phi_fine[base_i],
            rtol=0.0, atol=2e-12,
        )
        np.testing.assert_allclose(
            quarter_phi[1, k] - phi[k], phi_fine[q3_i] - phi_fine[base_i],
            rtol=0.0, atol=2e-12,
        )


def test_selected_quarter_panels_match_full_values_and_skip_unused_sigma(
    monkeypatch,
):
    model = LCDM_SG(
        Omega_bh2=0.0223828, Omega_ch2=0.1201075, H0=67.32117,
        DN_eff=0.0, A_s=2.100549e-9, r=0.01, n_t=0.0, cr=1.0,
        T_re=2000.0, DN_re=10.0, kappa10=0.01,
    )
    n_re = model.derived_param["N_inf"] - model.derived_param["N_re"]
    nodes = n_re + 0.1 + np.arange(7, dtype=np.float64) * 0.005
    phi, _, _, _, kink_index, kink_fraction, _ = fast_phi_s2_split(
        model, nodes, 0.0
    )
    sigma_nodes = exact_background.sigma_vec(nodes, model, 0.0)
    mask = np.array([True, False, True, False, True, False], dtype=bool)
    original_sigma_vec = exact_background.sigma_vec
    probe_counts = []

    def counted_sigma_vec(values, *args, **kwargs):
        probe_counts.append(np.asarray(values).size)
        return original_sigma_vec(values, *args, **kwargs)

    monkeypatch.setattr(exact_background, "sigma_vec", counted_sigma_vec)
    full = fast_phi_s2_quarters(
        model, nodes, 0.0, phi_grid=phi, kink_index=kink_index,
        kink_fraction=kink_fraction, sigma_nodes=sigma_nodes,
    )
    full_probe_count = sum(probe_counts)
    probe_counts.clear()

    workspace = model._fast_phi_s2_workspace
    workspace["Phi_mid"][mask] = 0.5 * (phi[:-1] + phi[1:])[mask] - 1.0e-4
    model._fast_quarter_sampling_masked = True
    selected = fast_phi_s2_quarters(
        model, nodes, 0.0, phi_grid=phi, kink_index=kink_index,
        kink_fraction=kink_fraction, sigma_nodes=sigma_nodes,
    )
    del model._fast_quarter_sampling_masked

    expected_probe_count = 6 * int(mask.sum())
    assert sum(probe_counts) == expected_probe_count
    assert expected_probe_count < full_probe_count
    np.testing.assert_array_equal(selected[:, mask], full[:, mask])
    np.testing.assert_array_equal(selected[:, ~mask], 0.0)
