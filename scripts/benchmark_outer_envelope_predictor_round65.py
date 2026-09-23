"""Round 65 standalone outer predictor: frozen transfer, algebraic envelope.

The candidate uses one full propagation at the initial outer point and then
predicts the spectrum at the updated background without a second propagation.
It freezes the transfer dynamics and rescales every assembled column by the
exact algebraic ``S2`` and horizon-start normalization factors.  This is a
diagnostic only; production outer-loop semantics are untouched.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_prufer_oracle import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _relative(a, b):
    return float(abs(a - b) / max(abs(b), 1.0e-300))


def _tail_end(snapshot, mode):
    phi = np.asarray(snapshot[1], dtype=np.float64)
    j0 = int(snapshot[5][mode])
    z0 = float(snapshot[6][mode])
    z_tail = float(snapshot[15])
    phi0 = float(phi[j0])
    k = j0
    zz = z0
    while k < len(phi) - 1 and zz < z_tail:
        k += 1
        zz = z0 + float(phi[k]) - phi0
    if zz < z_tail:
        return len(phi) - 1
    return max(j0, k - 1)


def _column_nodes(snapshot):
    nv = len(snapshot[0])
    n_coarse = int(snapshot[12])
    col_step = int(snapshot[13])
    return np.asarray(
        [nv - 1 if slot == n_coarse - 1 else col_step * slot
         for slot in range(n_coarse)], dtype=np.int64)


def envelope_factors(first_snapshot, updated_snapshot):
    """Return the frozen-transfer algebraic factors for all output columns.

    The first-state normalization is
    ``exp(2*delta_z0) * S2_old[j0]/S2_new[j0]``.  Assuming the propagated
    transfer state is unchanged, each output column then gains the exact
    endpoint factor ``S2_new[k]/S2_old[k]``.  The gate is deliberately strict:
    a frequency-start index or output layout change makes the candidate
    ineligible instead of guessing across a changed discrete representation.
    """
    old_j0 = np.asarray(first_snapshot[5], dtype=np.int64)
    new_j0 = np.asarray(updated_snapshot[5], dtype=np.int64)
    if not np.array_equal(old_j0, new_j0):
        raise ValueError("frequency start indices changed")
    old_s2 = np.asarray(first_snapshot[3], dtype=np.float64)
    new_s2 = np.asarray(updated_snapshot[3], dtype=np.float64)
    old_z0 = np.asarray(first_snapshot[6], dtype=np.float64)
    new_z0 = np.asarray(updated_snapshot[6], dtype=np.float64)
    initial = np.exp(2.0 * (new_z0 - old_z0))
    initial *= old_s2[old_j0] / new_s2[new_j0]
    nodes_old = _column_nodes(first_snapshot)
    nodes_new = _column_nodes(updated_snapshot)
    if not np.array_equal(nodes_old, nodes_new):
        raise ValueError("output column layout changed")
    factors = np.empty((old_j0.size, nodes_old.size), dtype=np.float64)
    for mode, j0 in enumerate(old_j0):
        end_old = _tail_end(first_snapshot, mode)
        end_new = _tail_end(updated_snapshot, mode)
        if end_old != end_new:
            raise ValueError("tail handoff index changed")
        factors[mode] = initial[mode] * (
            new_s2[nodes_old] / old_s2[nodes_old])
    if not np.isfinite(factors).all() or np.any(factors <= 0.0):
        raise FloatingPointError("non-finite envelope factor")
    return factors


def _copy_snapshot(snapshot):
    return tuple(value.copy() if isinstance(value, np.ndarray) else value
                 for value in snapshot)


def _run_capture(name):
    saved = (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
             FS._OUTER_FULL_REUSE_SIGMA_TOL,
             FS._OUTER_FULL_REUSE_FHOR_TOL)
    snapshots = []
    try:
        def gen(*args, **kwargs):
            return saved[0](*args, **kwargs)

        def prep(*args, **kwargs):
            return saved[1](*args, **kwargs)

        def solve(*args, **kwargs):
            result = saved[2](*args, **kwargs)
            snapshots.append(_copy_snapshot(args))
            return result

        FS.gen_fast = gen
        FS.prep_frequency_only = prep
        FS.solve_kernel = solve
        FS.apply_accuracy_mode("fast")
        FS.set_threads(2)
        FS._OUTER_FULL_REUSE_SIGMA_TOL = 0.0
        FS._OUTER_FULL_REUSE_FHOR_TOL = 0.0
        model = LCDM_SG(**CASES[name])
        FS.SGWB_iter_fast(model, kink_split=True, freq_grid="goal")
    finally:
        (FS.gen_fast, FS.prep_frequency_only, FS.solve_kernel,
         FS._OUTER_FULL_REUSE_SIGMA_TOL,
         FS._OUTER_FULL_REUSE_FHOR_TOL) = saved
    if len(snapshots) < 2:
        raise RuntimeError("outer run did not produce two propagation snapshots")
    first = snapshots[0]
    second = snapshots[1]
    full_first = list(_copy_snapshot(first))
    full_first[11] = 1
    for index in (16, 17, 18):
        full_first[index] = np.zeros_like(first[index])
    full_first[22] = np.full_like(first[22], -1.0)
    FS.solve_kernel(*tuple(full_first))
    factors = envelope_factors(full_first, second)
    predicted = [np.asarray(full_first[index], dtype=np.float64) * factors
                 for index in (16, 17, 18)]
    actual = [np.asarray(second[index], dtype=np.float64) for index in (16, 17, 18)]
    predicted_integrand = predicted[0][:, -1] - predicted[1][:, -1]
    actual_integrand = actual[0][:, -1] - actual[1][:, -1]
    freqs = np.log10(np.asarray(second[10], dtype=np.float64))
    g2_pred = float(np.sum(FS._pchip_integrals_vectorized(
        freqs, predicted_integrand)) * FS.ln10)
    g2_actual = float(np.sum(FS._pchip_integrals_vectorized(
        freqs, actual_integrand)) * FS.ln10)
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    dn_pred = float(gp.Neff0 * g2_pred / omega_nu)
    dn_actual = float(gp.Neff0 * g2_actual / omega_nu)
    log_pred = np.log10(np.maximum(predicted_integrand, 1.0e-300))
    log_actual = np.log10(np.maximum(actual_integrand, 1.0e-300))
    return {
        "point": name,
        "status": "ok",
        "frequency_start_indices_stable": True,
        "tail_handoff_stable": True,
        "candidate_dn_gw": dn_pred,
        "oracle_updated_dn_gw": dn_actual,
        "candidate_dn_relative": _relative(dn_pred, dn_actual),
        "spectrum_max_abs_dex": float(np.max(np.abs(log_pred - log_actual))),
        "history_max_abs_relative": float(max(
            _relative(predicted[i], actual[i]) for i in range(3))),
        "candidate_correction_cost_us": _measure_correction_cost(
            full_first, second),
        "production_outer_kernel_calls": len(snapshots),
        "candidate_propagation_calls": 1,
    }


def _measure_correction_cost(first, updated, repeats=100):
    factors = envelope_factors(first, updated)
    start = time.perf_counter()
    for _ in range(repeats):
        _ = [np.asarray(first[index]) * factors for index in (16, 17, 18)]
    elapsed = (time.perf_counter() - start) / repeats
    return elapsed * 1.0e6


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    rows = []
    for name in args.case or ["default", "lowT", "highT", "stiff", "high_kappa"]:
        try:
            rows.append(_run_capture(name))
        except (FloatingPointError, RuntimeError, ValueError) as exc:
            rows.append({"point": name, "status": "ineligible", "reason": str(exc)})
    payload = {
        "schema_version": 1,
        "candidate": "outer_frozen_transfer_envelope_predictor_round65",
        "generated_commit": os.popen("git rev-parse HEAD").read().strip(),
        "resources": telemetry(workers=1, threads=2),
        "production_unchanged": True,
        "rows": rows,
        "semantics": (
            "One full initial propagation plus algebraic S2/horizon-start "
            "rescaling; updated full production propagation is the oracle."),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n",
                                 encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
