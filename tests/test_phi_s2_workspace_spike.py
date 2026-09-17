import numpy as np

from scripts.benchmark_phi_s2_workspace_spike import (
    workspace_primitive_from_nodes,
)


def test_workspace_reuses_output_buffers_without_changing_values():
    nv = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.float64)
    nodes = np.array([1.0, 1.1, 1.2, 1.3], dtype=np.float64)
    workspace = {}
    first = workspace_primitive_from_nodes(
        nv, nodes, nodes, -1, 0.0, 0.0, 0.0, 0.0, workspace)
    first_values = tuple(value.copy() for value in first[:4])
    second = workspace_primitive_from_nodes(
        nv, nodes, nodes, -1, 0.0, 0.0, 0.0, 0.0, workspace)

    for before, after in zip(first_values, second[:4]):
        np.testing.assert_array_equal(before, after)
    for before, after in zip(first[:4], second[:4]):
        assert before is after
