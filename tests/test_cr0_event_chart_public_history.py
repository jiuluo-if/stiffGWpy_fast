"""Focused contracts for the independent public-history validation harness."""

import numpy as np

from scripts.validate_cr0_event_chart_public_history import (
    _assemble_tail_history,
    chart_history_indices,
    public_history_indices,
)


def test_public_history_indices_match_fast_coarse_layout_and_keep_terminal():
    nv = np.arange(19, dtype=np.float64) * 0.005
    idx_out = public_history_indices(nv, col_step=8)

    np.testing.assert_array_equal(idx_out, np.array([0, 8, 16, 18]))


def test_chart_mapping_keeps_pre_event_slots_unmapped_and_post_nodes_exact():
    public_Nv = np.arange(6, dtype=np.float64) * 0.01
    idx_out = np.array([0, 1, 2, 4, 5], dtype=np.int64)
    chart_nodes = np.array([0.02, 0.03, 0.04, 0.05])
    N_star = 0.02

    nodes, mapping = chart_history_indices(public_Nv, idx_out, chart_nodes, N_star)

    np.testing.assert_array_equal(nodes, np.array([0.0, 0.01, 0.02, 0.04, 0.05]))
    np.testing.assert_array_equal(mapping, np.array([-1, -1, 0, 2, 3]))


def test_tail_history_includes_sample_at_the_anchor():
    outputs = [np.zeros((1, 1), dtype=np.float64) for _ in range(3)]

    _assemble_tail_history(
        *outputs,
        0,
        np.array([0], dtype=np.int64),
        np.array([0.0]),
        0.0,
        0.0,
        1.0,
        1.0,
        1.0,
        1.0,
        10.0,
        1.0,
        np.array([1.0]),
        np.array([1.0]),
        1,
    )

    assert all(np.isfinite(channel[0, 0]) for channel in outputs)
    assert any(channel[0, 0] != 0.0 for channel in outputs)
