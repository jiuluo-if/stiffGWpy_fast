"""Standalone range-reduced exponential transfer twin.

The production phase-subdivision decision remains exact.  Only the ``exp(z)``
value consumed by each transfer step is replaced by a range-reduced Taylor
polynomial; trigonometric calls and all propagation/assembly ordering remain
the production formulas.  This is diagnostic and non-bitwise.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES, _make_args, _prepared  # noqa: E402
from scripts.benchmark_prufer_oracle import CASES as ORACLE_CASES  # noqa: E402
from scripts.benchmark_same_grid_reference import CASES as COVERAGE_CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

EXPERIMENT_CASES = dict(COVERAGE_CASES)
EXPERIMENT_CASES.update({
    name: ORACLE_CASES[name]
    for name in ("edge_r_hi", "edge_tre_hi", "edge_kap_hi")
})

_LN2 = math.log(2.0)
_INV_LN2 = 1.0 / _LN2


@njit(inline="always", cache=True)
def _exp_range(z):
    k = int(math.floor(z * _INV_LN2 + 0.5))
    r = z - k * _LN2
    p = 1.0 / 87178291200.0
    p = 1.0 / 6227020800.0 + r * p
    p = 1.0 / 479001600.0 + r * p
    p = 1.0 / 39916800.0 + r * p
    p = 1.0 / 3628800.0 + r * p
    p = 1.0 / 362880.0 + r * p
    p = 1.0 / 40320.0 + r * p
    p = 1.0 / 5040.0 + r * p
    p = 1.0 / 720.0 + r * p
    p = 1.0 / 120.0 + r * p
    p = 1.0 / 24.0 + r * p
    p = 1.0 / 6.0 + r * p
    p = 0.5 + r * p
    p = 1.0 + r * p
    p = 1.0 + r * p
    return math.ldexp(p, k)


@njit(inline="always", cache=True)
def _scaled_step_range(xh, yh, z_mid, h):
    w = _exp_range(z_mid)
    w2 = w * w
    if w2 >= 1.0:
        omega = math.sqrt(w2 - 1.0)
        c = math.cos(omega * h)
        si = math.sin(omega * h) / omega
    else:
        omega = math.sqrt(1.0 - w2)
        x = omega * h
        c = 1.0 + 0.5 * x * x
        si = h * (1.0 + x * x / 6.0)
    return (c - si) * xh - w * si * yh, w * si * xh + (c + si) * yh


def _local_errors():
    exp_error = 0.0
    transfer_error = 0.0
    for z in np.linspace(-1.0, 5.1, 2001):
        actual = _exp_range(float(z))
        expected = math.exp(float(z))
        exp_error = max(exp_error, abs(actual - expected) / expected)
        base = FS.scaled_step(0.2, -0.1, float(z), 0.005)
        candidate = _scaled_step_range(0.2, -0.1, float(z), 0.005)
        transfer_error = max(
            transfer_error, abs(candidate[0] - base[0]), abs(candidate[1] - base[1]))
    return {
        "exp_relative_max": float(exp_error),
        "transfer_absolute_max": float(transfer_error),
    }


@njit(inline="always", cache=True)
def _phase_segment_range(xh, yh, z_start, z_end, h_step, phase_max):
    z_mid = 0.5 * (z_start + z_end)
    if phase_max <= 0.0 or z_mid <= 0.0:
        n_sub = 1
    else:
        n_sub = int(math.ceil(h_step * math.exp(z_mid) / phase_max))
        if n_sub < 1:
            n_sub = 1
    if n_sub == 1:
        return _scaled_step_range(xh, yh, z_mid, h_step)
    h_sub = h_step / n_sub
    dz = (z_end - z_start) / n_sub
    for sub in range(n_sub):
        z_sub_mid = z_start + (sub + 0.5) * dz
        xh, yh = _scaled_step_range(xh, yh, z_sub_mid, h_sub)
    return xh, yh


@njit(parallel=True, cache=True)
def solve_kernel_range(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re,
):
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Pt = P_t[mode]
        fp_i = fp_freq[mode]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, _exp_range(z0) * S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh = 0.0
        lyh = yh
        last_z = zz
        if k % col_step == 0:
            if assemble:
                slot = k // col_step
                if slot >= n_coarse - 1:
                    slot = n_coarse - 1
                FS.assemble_main(Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, mode, n_coarse - 1,
                             S2[k], xh, yh, zz, Pt)
        if zz < z_tail:
            while k < nv - 1 and (z0 + Phi_grid[k] - Phi0) < z_tail:
                h_step = h_arr[k] if h_arr is not None else h
                z_node = z0 + Phi_grid[k] - Phi0
                z_end = z0 + Phi_grid[k + 1] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    h_left = h_step * kink_fraction
                    xh, yh = _phase_segment_range(
                        xh, yh, z_node, z_break, h_left, phase_max)
                    xh, yh = _phase_segment_range(
                        xh, yh, z_break, z_end, h_step - h_left, phase_max)
                else:
                    xh, yh = _phase_segment_range(
                        xh, yh, z_node, z_end, h_step, phase_max)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % col_step == 0:
                    if assemble:
                        slot = k // col_step
                        if slot >= n_coarse - 1:
                            slot = n_coarse - 1
                        FS.assemble_main(
                            Ogw, Oj, Opgw, mode, slot, S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(
                        Ogw, Oj, Opgw, mode, n_coarse - 1,
                        S2[k], xh, yh, zz, Pt)
                if zz < z_tail:
                    lxh = xh
                    lyh = yh
                    last_z = zz
            kend = nv - 1 if zz < z_tail else k - 1
            if kend < j0:
                kend = j0
        else:
            kend = j0
        if kend < nv - 1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend]) if Sv is not None else 0.0
            amp2 = (lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z)
                    + 2.0 * gamma * lxh * lyh * e_z)
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            if handoff_eps is not None:
                handoff_eps[mode] = (abs(1.5 * Sv[kend] - 1.0) * e_z
                                      if Sv is not None else 0.0)
            slot_start = kend // col_step
            if slot_start >= n_coarse - 1:
                slot_start = n_coarse - 1
            while slot_start < n_coarse:
                kk2 = col_step * slot_start
                if slot_start == n_coarse - 1:
                    kk2 = nv - 1
                if kk2 > kend:
                    break
                slot_start += 1
            for slot in range(slot_start, n_coarse):
                if not (assemble or slot == n_coarse - 1):
                    continue
                kk2 = col_step * slot
                if slot == n_coarse - 1:
                    kk2 = nv - 1
                FS.assemble_tail(
                    Ogw, Oj, Opgw, mode, slot, kk2, coeff, eNz, fp_i, Pt,
                    ev_minus, fp_minus)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _run_kernel(case_name, threads=2, repeats=15):
    model, common = _prepared(case_name, threads, cases=EXPERIMENT_CASES)
    args = _make_args(common)
    FS.solve_kernel(*args)
    solve_kernel_range(*args)
    baseline_times = []
    candidate_times = []
    baseline = None
    candidate = None
    for rep in range(repeats):
        base_args = list(_make_args(common))
        cand_args = list(_make_args(common))
        start = time.perf_counter()
        FS.solve_kernel(*base_args)
        baseline_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        solve_kernel_range(*cand_args)
        candidate_times.append(time.perf_counter() - start)
        if rep == 0:
            baseline = base_args
            candidate = cand_args
    return {
        "case": case_name,
        "threads": int(threads),
        "repeats": int(repeats),
        "baseline_median_ms": statistics.median(baseline_times) * 1000.0,
        "candidate_median_ms": statistics.median(candidate_times) * 1000.0,
        "baseline_p95_ms": float(np.percentile(baseline_times, 95)) * 1000.0,
        "candidate_p95_ms": float(np.percentile(candidate_times, 95)) * 1000.0,
        "candidate_over_baseline": statistics.median(candidate_times) / statistics.median(baseline_times),
        "digest_equal": all(_digest(baseline[i]) == _digest(candidate[i])
                             for i in (16, 17, 18)),
        "max_abs_output_delta": max(float(np.max(np.abs(baseline[i] - candidate[i])))
                                     for i in (16, 17, 18)),
    }


def _oracle_metrics(case_name, threads=2):
    model, common = _prepared(case_name, threads, cases=EXPERIMENT_CASES)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_range(*candidate)
    freqs = np.asarray(model.f, dtype=np.float64)
    dn_eff = float(model.cosmo_param["DN_eff"])
    oracle = REF.spectrum_reference(model, freqs, dn_eff, z_tail=5.0,
                                    rtol=1e-8, workers=1)
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    oracle_g2 = REF.integrate_spectrum(freqs, oracle[0], oracle[1])[0]
    oracle_dn = float(gp.Neff0 * oracle_g2 / omega_nu)
    rows = {}
    for label, args in (("baseline", baseline), ("candidate", candidate)):
        ogw = args[16][:, -1]
        oj = args[17][:, -1]
        g2 = REF.integrate_spectrum(freqs, ogw, oj)[0]
        dn_gw = float(gp.Neff0 * g2 / omega_nu)
        dex = float(np.max(np.abs(
            np.log10(np.maximum(ogw - oj, 1e-300))
            - np.log10(np.maximum(oracle[0] - oracle[1], 1e-300)))))
        rows[label] = {
            "dn_gw": dn_gw,
            "dn_relative_to_oracle": abs(dn_gw - oracle_dn) / max(abs(oracle_dn), 1e-300),
            "spectrum_max_dex_to_oracle": dex,
        }
    return {
        "case": case_name,
        "threads": int(threads),
        "oracle_dn_gw": oracle_dn,
        "baseline": rows["baseline"],
        "candidate": rows["candidate"],
        "candidate_minus_baseline_dn_relative": abs(
            rows["candidate"]["dn_gw"] - rows["baseline"]["dn_gw"]
        ) / max(abs(rows["baseline"]["dn_gw"]), 1e-300),
    }


def _outer_adapter(*args, **kwargs):
    del kwargs
    return solve_kernel_range(*args)


def _run_outer(case_name, repeats=25, threads=2):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    original = FS.solve_kernel
    baseline_times = []
    candidate_times = []
    baseline_model = candidate_model = None
    try:
        for _ in range(repeats):
            baseline_model = LCDM_SG(**CASES[case_name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(
                baseline_model, kink_split=True, freq_grid="goal",
                frequency_quadrature="pchip")
            baseline_times.append(time.perf_counter() - started)
            FS.solve_kernel = _outer_adapter
            candidate_model = LCDM_SG(**CASES[case_name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(
                candidate_model, kink_split=True, freq_grid="goal",
                frequency_quadrature="pchip")
            candidate_times.append(time.perf_counter() - started)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original

    def digest(model, name):
        return _digest(getattr(model, name))

    return {
        "case": case_name,
        "threads": int(threads),
        "repeats": int(repeats),
        "baseline_median_ms": statistics.median(baseline_times) * 1e3,
        "candidate_median_ms": statistics.median(candidate_times) * 1e3,
        "candidate_over_baseline": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_ms": float(np.percentile(baseline_times, 95) * 1e3),
        "candidate_p95_ms": float(np.percentile(candidate_times, 95) * 1e3),
        "f_equal": digest(baseline_model, "f") == digest(candidate_model, "f"),
        "spectrum_equal": digest(baseline_model, "log10OmegaGW") == digest(candidate_model, "log10OmegaGW"),
        "DN_gw_equal": digest(baseline_model, "DN_gw") == digest(candidate_model, "DN_gw"),
        "g2_equal": digest(baseline_model, "g2") == digest(candidate_model, "g2"),
        "w2_equal": digest(baseline_model, "w2") == digest(candidate_model, "w2"),
        "failure_equal": getattr(baseline_model, "fast_failure_reason", None) == getattr(candidate_model, "fast_failure_reason", None),
        "converged_equal": getattr(baseline_model, "SGWB_converge", False) == getattr(candidate_model, "SGWB_converge", False),
        "spectrum_max_dex": float(np.max(np.abs(
            np.asarray(candidate_model.log10OmegaGW)
            - np.asarray(baseline_model.log10OmegaGW)))),
        "DN_gw_relative": abs(float(candidate_model.DN_gw[-1]) - float(baseline_model.DN_gw[-1])) / max(abs(float(baseline_model.DN_gw[-1])), 1e-300),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(EXPERIMENT_CASES))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=15)
    parser.add_argument("--oracle", action="store_true")
    parser.add_argument("--outer", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = args.case or list(CASES)
    payload = {
        "candidate": "range_reduced_exponential_transfer",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "local_errors": _local_errors(),
        "settings": vars(args),
        "cases": ([ _oracle_metrics(name, args.threads) for name in cases ]
                   if args.oracle else
                   [ _run_outer(name, args.repeats, args.threads) for name in cases ]
                   if args.outer else
                   [ _run_kernel(name, args.threads, args.repeats) for name in cases ]),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
