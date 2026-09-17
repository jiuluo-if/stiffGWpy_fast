import numpy as np

from scripts.benchmark_phi_s2_numba_spike import numba_primitive_from_nodes


def test_numba_primitive_matches_reference_formula_without_kink():
    nv = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.float64)
    nodes = np.array([1.0, 1.1, 1.2, 1.3], dtype=np.float64)
    got = numba_primitive_from_nodes(
        nv, nodes, nodes, -1, 0.0, 0.0, 0.0, 0.0)

    h = np.diff(nv)
    mid = 0.5 * (nodes[:-1] + nodes[1:])
    quarter = 0.75 * nodes[:-1] + 0.25 * nodes[1:]
    integral = h * (nodes[:-1] + 4.0 * mid + nodes[1:]) / 6.0
    F = np.concatenate(([0.0], np.cumsum(integral)))
    F_mid = F[:-1] + h * (nodes[:-1] + 4.0 * quarter + mid) / 12.0
    expected = (
        1.5 * F - nv + nv[0],
        1.5 * F_mid - 0.5 * (nv[:-1] + nv[1:]) + nv[0],
        np.exp(3.0 * F - 4.0 * nv),
        np.exp(-0.5 * (3.0 * F - 4.0 * nv)),
    )

    for actual, reference in zip(got[:4], expected):
        np.testing.assert_allclose(actual, reference, rtol=0.0, atol=1e-14)
