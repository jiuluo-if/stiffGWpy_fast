"""Standalone residual-controlled two-step transfer spike.

The production kernel freezes ``z`` at the midpoint of each native grid
interval.  This diagnostic tries one midpoint exponential over two adjacent
intervals when a cheap local variation/phase defect is below a threshold; all
other intervals use the exact production ``_phase_segment`` fallback.  It is
non-strict and diagnostic only: production code is not changed.
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
for _name in ("edge_r_hi", "edge_tre_hi", "edge_kap_hi"):
    EXPERIMENT_CASES[_name] = ORACLE_CASES[_name]


@njit(inline="always")
def _block_defect(z_start, z_end, width, phase_cap):
    """Cheap midpoint-Magnus variation proxy for one two-step block."""
    z_mid = 0.5 * (z_start + z_end)
    w = math.exp(z_mid)
    phase = width * w
    variation = abs(z_end - z_start)
    return variation * max(1.0, phase) * phase, phase


@njit(parallel=True, cache=True)
def solve_kernel_residual_block(
    Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
    fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, Ogw, Oj,
    Opgw, h_arr, Sv, phase_max, handoff_eps, kink_index, kink_fraction,
    phi_re, defect_threshold, block_phase_cap, accepted_counts, fallback_counts,
):
    """Production-signature twin with guarded two-native-step blocks."""
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        Pt = P_t[mode]
        fp_i = fp_freq[mode]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0) * S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh, lyh, last_z = 0.0, yh, zz
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
            while k < nv - 1 and zz < z_tail:
                can_block = (k + 2 < nv and k != kink_index
                             and k + 1 != kink_index
                             and not (0.0 < kink_fraction < 1.0
                                      and k <= kink_index < k + 2))
                if can_block:
                    h0 = h_arr[k] if h_arr is not None else h
                    h1 = h_arr[k + 1] if h_arr is not None else h
                    z_next = z0 + Phi_grid[k + 2] - Phi0
                    block_h = h0 + h1
                    defect, phase = _block_defect(zz, z_next, block_h, block_phase_cap)
                    if (defect <= defect_threshold and phase <= block_phase_cap
                            and z_next < z_tail):
                        xh, yh = FS.scaled_step(
                            xh, yh, 0.5 * (zz + z_next), block_h)
                        k += 2
                        zz = z_next
                        accepted_counts[mode] += 1
                    else:
                        h_step = h0
                        z_mid_step = z0 + Phi_mid[k] - Phi0
                        z_end = 2.0 * z_mid_step - zz
                        if k == kink_index and 0.0 < kink_fraction < 1.0:
                            z_break = z0 + phi_re - Phi0
                            z_end = z0 + Phi_grid[k + 1] - Phi0
                            h_left = h_step * kink_fraction
                            xh, yh = FS._phase_segment(
                                xh, yh, zz, z_break, h_left, phase_max)
                            xh, yh = FS._phase_segment(
                                xh, yh, z_break, z_end,
                                h_step - h_left, phase_max)
                        else:
                            xh, yh = FS._phase_segment(
                                xh, yh, zz, z_end, h_step, phase_max)
                        k += 1
                        zz = z0 + Phi_grid[k] - Phi0
                        fallback_counts[mode] += 1
                else:
                    h_step = h_arr[k] if h_arr is not None else h
                    z_mid_step = z0 + Phi_mid[k] - Phi0
                    z_end = 2.0 * z_mid_step - zz
                    if k == kink_index and 0.0 < kink_fraction < 1.0:
                        z_break = z0 + phi_re - Phi0
                        z_end = z0 + Phi_grid[k + 1] - Phi0
                        h_left = h_step * kink_fraction
                        xh, yh = FS._phase_segment(
                            xh, yh, zz, z_break, h_left, phase_max)
                        xh, yh = FS._phase_segment(
                            xh, yh, z_break, z_end,
                            h_step - h_left, phase_max)
                    else:
                        xh, yh = FS._phase_segment(
                            xh, yh, zz, z_end, h_step, phase_max)
                    k += 1
                    zz = z0 + Phi_grid[k] - Phi0
                    fallback_counts[mode] += 1

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
                    lxh, lyh, last_z = xh, yh, zz

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
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _candidate_args(common, threshold, phase_cap):
    args = _make_args(common)
    counts = np.zeros(len(args[5]), dtype=np.int64)
    fallbacks = np.zeros(len(args[5]), dtype=np.int64)
    return args + (float(threshold), float(phase_cap), counts, fallbacks), counts, fallbacks


def _run_kernel(case_name, threshold, repeats=30, threads=2, phase_cap=0.5):
    model, common = _prepared(case_name, threads, cases=EXPERIMENT_CASES)
    baseline = _make_args(common)
    candidate, accepted, fallback = _candidate_args(common, threshold, phase_cap)
    FS.solve_kernel(*baseline)
    solve_kernel_residual_block(*candidate)
    first_divergence = None
    for name, index in (("Ogw", 16), ("Oj", 17), ("Opgw", 18), ("handoff_eps", 22)):
        if not np.array_equal(baseline[index], candidate[index]):
            where = np.argwhere(baseline[index] != candidate[index])[0]
            first_divergence = {
                "variable": name,
                "index": tuple(int(v) for v in where),
            }
            break
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    rel = np.abs(cand_obs - base_obs) / np.maximum(np.abs(base_obs), 1e-300)
    base_g2 = float(np.sum(FS._pchip_integrals_vectorized(model.f, base_obs)) * FS.ln10)
    cand_g2 = float(np.sum(FS._pchip_integrals_vectorized(model.f, cand_obs)) * FS.ln10)
    baseline_times, candidate_times = [], []
    for _ in range(repeats):
        for arr in baseline[16:19]:
            arr.fill(0.0)
        baseline[22].fill(-1.0)
        for arr in candidate[16:19]:
            arr.fill(0.0)
        candidate[22].fill(-1.0)
        accepted.fill(0)
        fallback.fill(0)
        started = time.perf_counter()
        FS.solve_kernel(*baseline)
        baseline_times.append(time.perf_counter() - started)
        started = time.perf_counter()
        solve_kernel_residual_block(*candidate)
        candidate_times.append(time.perf_counter() - started)
    bmed = statistics.median(baseline_times)
    cmed = statistics.median(candidate_times)
    return {
        "case": case_name,
        "threshold": float(threshold),
        "phase_cap": float(phase_cap),
        "threads": threads,
        "repeats": repeats,
        "baseline_median_ms": bmed * 1e3,
        "candidate_median_ms": cmed * 1e3,
        "candidate_over_baseline": cmed / bmed,
        "baseline_p95_ms": float(np.percentile(baseline_times, 95) * 1e3),
        "candidate_p95_ms": float(np.percentile(candidate_times, 95) * 1e3),
        "accepted_blocks": int(np.sum(accepted)),
        "fallback_steps": int(np.sum(fallback)),
        "total_candidate_steps": int(2 * np.sum(accepted) + np.sum(fallback)),
        "bitwise_equal": first_divergence is None,
        "first_divergence": first_divergence,
        "spectrum_rel_p95": float(np.percentile(rel, 95)),
        "spectrum_rel_max": float(np.max(rel)),
        "DN_gw_relative": abs(cand_g2 - base_g2) / max(abs(base_g2), 1e-300),
        "digest_Ogw_equal": _digest(baseline[16]) == _digest(candidate[16]),
        "digest_Oj_equal": _digest(baseline[17]) == _digest(candidate[17]),
        "digest_Opgw_equal": _digest(baseline[18]) == _digest(candidate[18]),
    }


def _run_oracle(case_name, threshold, threads=2, phase_cap=0.5):
    """Compare one candidate spectrum with the independent continuous oracle."""
    model, common = _prepared(case_name, threads, cases=EXPERIMENT_CASES)
    baseline = _make_args(common)
    candidate, accepted, fallback = _candidate_args(common, threshold, phase_cap)
    FS.solve_kernel(*baseline)
    solve_kernel_residual_block(*candidate)
    freqs = np.asarray(model.f, dtype=np.float64)
    dn_eff = float(model.cosmo_param["DN_eff"])
    oracle = REF.spectrum_reference(
        model, freqs, dn_eff, z_tail=5.0, rtol=1e-8, workers=1)
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
        "threshold": float(threshold),
        "phase_cap": float(phase_cap),
        "threads": threads,
        "oracle_dn_gw": oracle_dn,
        "accepted_blocks": int(np.sum(accepted)),
        "fallback_steps": int(np.sum(fallback)),
        "baseline": rows["baseline"],
        "candidate": rows["candidate"],
        "candidate_minus_baseline_dn_relative": abs(
            rows["candidate"]["dn_gw"] - rows["baseline"]["dn_gw"]
        ) / max(abs(rows["baseline"]["dn_gw"]), 1e-300),
    }


def _outer_adapter(*args):
    threshold = float(getattr(_outer_adapter, "threshold", 0.0))
    phase_cap = float(getattr(_outer_adapter, "phase_cap", 0.5))
    counts = np.zeros(len(args[5]), dtype=np.int64)
    fallback = np.zeros(len(args[5]), dtype=np.int64)
    return solve_kernel_residual_block(
        *args, threshold, phase_cap, counts, fallback)


def _run_outer(case_name, threshold, repeats=25, threads=2, phase_cap=0.5):
    """Alternating full production/candidate end-to-end comparison."""
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    original = FS.solve_kernel
    baseline_times, candidate_times = [], []
    baseline_model = candidate_model = None
    try:
        _outer_adapter.threshold = float(threshold)
        _outer_adapter.phase_cap = float(phase_cap)
        warm_baseline = LCDM_SG(**CASES[case_name])
        FS.SGWB_iter_fast(
            warm_baseline, kink_split=True, freq_grid="goal",
            frequency_quadrature="pchip")
        FS.solve_kernel = _outer_adapter
        warm_candidate = LCDM_SG(**CASES[case_name])
        FS.SGWB_iter_fast(
            warm_candidate, kink_split=True, freq_grid="goal",
            frequency_quadrature="pchip")
        FS.solve_kernel = original
        for _ in range(repeats):
            baseline_model = LCDM_SG(**EXPERIMENT_CASES[case_name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(
                baseline_model, kink_split=True, freq_grid="goal",
                frequency_quadrature="pchip")
            baseline_times.append(time.perf_counter() - started)
            FS.solve_kernel = _outer_adapter
            candidate_model = LCDM_SG(**EXPERIMENT_CASES[case_name])
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
        "threshold": float(threshold),
        "phase_cap": float(phase_cap),
        "threads": threads,
        "repeats": repeats,
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


def _run_guard_record(case_name, threshold, mode):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    model = LCDM_SG(**EXPERIMENT_CASES[case_name])
    result = FS.SGWB_iter_fast(
        model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
    return {
        "case": case_name,
        "threshold": float(threshold),
        "mode": mode,
        "status": "baseline_physical_guard",
        "failure_reason": getattr(model, "fast_failure_reason", None),
        "candidate_run": False,
        "production_result_is_none": result is None,
    }


def _safe_run(runner, case, threshold, args):
    try:
        if runner == "outer":
            return _run_outer(case, threshold, args.repeats, args.threads, args.phase_cap)
        if runner == "oracle":
            return _run_oracle(case, threshold, args.threads, args.phase_cap)
        return _run_kernel(case, threshold, args.repeats, args.threads, args.phase_cap)
    except RuntimeError as exc:
        if "fast preparation failed" not in str(exc):
            raise
        return _run_guard_record(case, threshold, runner)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(EXPERIMENT_CASES))
    parser.add_argument("--threshold", type=float, nargs="+", default=[0.0, 1e-3, 3e-3])
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--phase-cap", type=float, default=0.5)
    parser.add_argument("--oracle", action="store_true")
    parser.add_argument("--outer", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = args.case or list(CASES)
    if args.outer:
        rows = [_safe_run("outer", case, threshold, args)
                for threshold in args.threshold for case in cases]
    elif args.oracle:
        rows = [_safe_run("oracle", case, threshold, args)
                for threshold in args.threshold for case in cases]
    else:
        rows = [_safe_run("kernel", case, threshold, args)
                for threshold in args.threshold for case in cases]
    payload = {
        "candidate": "residual_controlled_two_step_transfer",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "settings": vars(args),
        "rows": rows,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
