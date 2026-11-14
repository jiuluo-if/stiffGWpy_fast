"""Shared public-history helpers for the event-chart research harness."""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from stiffgwpy_fast import fast_sgwb as FS


def public_history_indices(Nv: np.ndarray, col_step: int = 8) -> np.ndarray:
    """Return the production coarse history layout, including terminal node."""
    grid = np.asarray(Nv, dtype=np.float64)
    if grid.ndim != 1 or grid.size < 2 or col_step < 1:
        raise ValueError("Nv must be a 1-D grid with at least two nodes and positive col_step")
    return np.unique(np.append(np.arange(0, grid.size, col_step), grid.size - 1)).astype(np.int64)


def chart_history_indices(
    public_Nv: np.ndarray,
    idx_out: np.ndarray,
    chart_N_abs: np.ndarray,
    N_star: float,
    *,
    atol: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """Map post-event public columns to chart nodes; pre-event slots map to -1."""
    nodes = np.asarray(public_Nv, dtype=np.float64)[np.asarray(idx_out, dtype=np.int64)]
    chart_nodes = np.asarray(chart_N_abs, dtype=np.float64)
    mapping = np.full(nodes.size, -1, dtype=np.int64)
    for slot, value in enumerate(nodes):
        if value < N_star:
            continue
        pos = int(np.searchsorted(chart_nodes, value))
        candidates = (pos - 1, pos)
        best = min(
            (q for q in candidates if 0 <= q < chart_nodes.size),
            key=lambda q: abs(float(chart_nodes[q]) - float(value)),
        )
        if abs(float(chart_nodes[best]) - float(value)) > atol:
            raise ValueError(f"public post-event node {value:.17g} is absent from chart")
        mapping[slot] = best
    return nodes, mapping


@njit(inline="always")
def _assemble_tail_history(
    Ogw, Oj, Opgw, mode, sample_chart_indices, sample_tau, anchor_tau,
    last_z, xh, yh, s2k, sigma, fp_freq, P_t, ev_minus, fp_minus,
    stop_slot,
):
    e_z = math.exp(-last_z)
    gamma = FS._tail_match_gamma(sigma)
    amp2 = (
        xh * xh + yh * yh * (1.0 + gamma * gamma * e_z * e_z)
        + 2.0 * gamma * xh * yh * e_z
    )
    coeff = math.sqrt(0.5 * s2k * amp2)
    e_nz = math.exp(anchor_tau - last_z)
    for slot in range(stop_slot):
        chart_index = sample_chart_indices[slot]
        if chart_index >= 0 and sample_tau[slot] > anchor_tau:
            FS.assemble_tail(
                Ogw, Oj, Opgw, mode, slot, chart_index, coeff, e_nz,
                fp_freq, P_t, ev_minus, fp_minus,
            )
