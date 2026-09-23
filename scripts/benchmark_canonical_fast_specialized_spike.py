"""Standalone guarded canonical-fast kernel specialization spike.

This candidate fixes only the already-guarded fast-profile constants and
optional-argument facts: col_step=8, h=0.005, z_tail=5, phase_max=0.25,
assemble=1, h_arr=None and Sv/handoff_eps present.  It intentionally keeps
the production modulo/division assembly schedule, phase segment, kink split,
tail matching and floating-point order unchanged, so it is independent of
the rejected P1 counted-assembly candidate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
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

from scripts.benchmark_phase_recurrence import (  # noqa: E402
    CASES,
    _make_args,
    _prepared,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CANONICAL_COL_STEP = 8
CANONICAL_H = 0.005
CANONICAL_Z_TAIL = 5.0
CANONICAL_PHASE_MAX = 0.25


@njit(parallel=True, cache=True)
def solve_kernel_canonical_specialized(
    Nv,
    Phi_grid,
    Phi_mid,
    S2,
    S2inv,
    j0s,
    z0s,
    P_t,
    ev_minus,
    fp_minus,
    fp_freq,
    assemble,
    n_coarse,
    col_step,
    h,
    z_tail,
    Ogw,
    Oj,
    Opgw,
    h_arr=None,
    Sv=None,
    phase_max=0.0,
    handoff_eps=None,
    kink_index=-1,
    kink_fraction=0.0,
    phi_re=0.0,
):
    """Canonical fast specialization; caller-side guard owns its contract."""
    nv = len(Nv)
    for m in prange(len(j0s)):
        j0 = j0s[m]
        z0 = z0s[m]
        Pt = P_t[m]
        fp_i = fp_freq[m]
        Phi0 = Phi_grid[j0]
        xh, yh = 0.0, math.exp(z0) * S2inv[j0]
        k = j0
        zz = z0 + Phi_grid[k] - Phi0
        lxh = 0.0
        lyh = yh
        last_z = zz
        if k % CANONICAL_COL_STEP == 0:
            slot = k // CANONICAL_COL_STEP
            if slot >= n_coarse - 1:
                slot = n_coarse - 1
            FS.assemble_main(Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
        elif k == nv - 1:
            FS.assemble_main(Ogw, Oj, Opgw, m, n_coarse - 1, S2[k], xh, yh, zz, Pt)
        if zz < CANONICAL_Z_TAIL:
            while k < nv - 1 and zz < CANONICAL_Z_TAIL:
                z_node = zz
                z_mid_step = z0 + Phi_mid[k] - Phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - Phi0
                    z_end = z0 + Phi_grid[k + 1] - Phi0
                    h_left = CANONICAL_H * kink_fraction
                    h_right = CANONICAL_H - h_left
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_break, h_left, CANONICAL_PHASE_MAX)
                    xh, yh = FS._phase_segment(xh, yh, z_break, z_end, h_right, CANONICAL_PHASE_MAX)
                else:
                    z_end = 2.0 * z_mid_step - z_node
                    xh, yh = FS._phase_segment(xh, yh, z_node, z_end, CANONICAL_H, CANONICAL_PHASE_MAX)
                k += 1
                zz = z0 + Phi_grid[k] - Phi0
                if k % CANONICAL_COL_STEP == 0:
                    slot = k // CANONICAL_COL_STEP
                    if slot >= n_coarse - 1:
                        slot = n_coarse - 1
                    FS.assemble_main(Ogw, Oj, Opgw, m, slot, S2[k], xh, yh, zz, Pt)
                elif k == nv - 1:
                    FS.assemble_main(Ogw, Oj, Opgw, m, n_coarse - 1, S2[k], xh, yh, zz, Pt)
                if zz < CANONICAL_Z_TAIL:
                    lxh = xh
                    lyh = yh
                    last_z = zz
            if zz < CANONICAL_Z_TAIL:
                kend = nv - 1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        if kend < nv - 1:
            s2k = S2[kend]
            e_z = math.exp(-last_z)
            gamma = FS._tail_match_gamma(Sv[kend])
            amp2 = lxh * lxh + lyh * lyh * (1.0 + gamma * gamma * e_z * e_z) + 2.0 * gamma * lxh * lyh * e_z
            coeff = math.sqrt(0.5 * s2k * amp2)
            eNz = math.exp(Nv[kend] - last_z)
            handoff_eps[m] = abs(1.5 * Sv[kend] - 1.0) * e_z
            slot_start = kend // CANONICAL_COL_STEP
            if slot_start >= n_coarse - 1:
                slot_start = n_coarse - 1
            while slot_start < n_coarse:
                kk2 = CANONICAL_COL_STEP * slot_start
                if slot_start == n_coarse - 1:
                    kk2 = nv - 1
                if kk2 > kend:
                    break
                slot_start += 1
            for slot in range(slot_start, n_coarse):
                kk2 = CANONICAL_COL_STEP * slot
                if slot == n_coarse - 1:
                    kk2 = nv - 1
                FS.assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _first_divergence(base, candidate):
    for name, index in (("Ogw", 16), ("Oj", 17), ("Opgw", 18), ("handoff_eps", 22)):
        lhs = np.ascontiguousarray(base[index], dtype=np.float64)
        rhs = np.ascontiguousarray(candidate[index], dtype=np.float64)
        if not np.array_equal(lhs, rhs):
            where = np.argwhere(lhs != rhs)
            pos = tuple(int(x) for x in where[0])
            return {"variable": name, "index": pos, "baseline": float(lhs[pos]), "candidate": float(rhs[pos])}
    return None


def _run_kernel(name, repeats=30, threads=2):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    candidate = _make_args(common)
    FS.solve_kernel(*baseline)
    solve_kernel_canonical_specialized(*candidate)
    first_divergence = _first_divergence(baseline, candidate)
    base_obs = baseline[16][:, -1] - baseline[17][:, -1]
    cand_obs = candidate[16][:, -1] - candidate[17][:, -1]
    base_dn = FS.integrate_frequency_pchip(model.f, base_obs)
    cand_dn = FS.integrate_frequency_pchip(model.f, cand_obs)
    baseline_times = []
    candidate_times = []
    for _ in range(repeats):
        args = _make_args(common)
        started = time.perf_counter()
        FS.solve_kernel(*args)
        baseline_times.append(time.perf_counter() - started)
        args = _make_args(common)
        started = time.perf_counter()
        solve_kernel_canonical_specialized(*args)
        candidate_times.append(time.perf_counter() - started)
    return {
        "case": name,
        "threads": threads,
        "repeats": repeats,
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "bitwise_equal": first_divergence is None,
        "first_divergence": first_divergence,
        "DN_gw_relative": abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
    }


def _run_outer(name, repeats=25, threads=2):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    original = FS.solve_kernel
    baseline_times = []
    candidate_times = []
    baseline_model = candidate_model = None
    try:
        for _ in range(repeats):
            baseline_model = LCDM_SG(**CASES[name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(baseline_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
            baseline_times.append(time.perf_counter() - started)
            FS.solve_kernel = solve_kernel_canonical_specialized
            candidate_model = LCDM_SG(**CASES[name])
            started = time.perf_counter()
            FS.SGWB_iter_fast(candidate_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
            candidate_times.append(time.perf_counter() - started)
            FS.solve_kernel = original
    finally:
        FS.solve_kernel = original
    return {
        "case": name,
        "threads": threads,
        "repeats": repeats,
        "baseline_median_s": statistics.median(baseline_times),
        "candidate_median_s": statistics.median(candidate_times),
        "median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
        "baseline_p95_s": float(np.percentile(baseline_times, 95)),
        "candidate_p95_s": float(np.percentile(candidate_times, 95)),
        "digest_spectrum_equal": _digest(baseline_model.log10OmegaGW) == _digest(candidate_model.log10OmegaGW),
        "digest_DN_gw_equal": _digest(baseline_model.DN_gw) == _digest(candidate_model.DN_gw),
        "failure_equal": getattr(baseline_model, "fast_failure_reason", None)
        == getattr(candidate_model, "fast_failure_reason", None),
        "converged_equal": getattr(baseline_model, "SGWB_converge", False)
        == getattr(candidate_model, "SGWB_converge", False),
        "spectrum_dex_max": float(
            np.max(np.abs(np.asarray(candidate_model.log10OmegaGW) - np.asarray(baseline_model.log10OmegaGW)))
        ),
        "DN_gw_relative": abs(float(candidate_model.DN_gw[-1]) - float(baseline_model.DN_gw[-1]))
        / max(abs(float(baseline_model.DN_gw[-1])), 1e-300),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--repeats", type=int, default=30)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--outer", action="store_true")
    parser.add_argument("--output")
    args = parser.parse_args()
    os.environ.setdefault("FAST_THREADS", str(args.threads))
    FS.apply_accuracy_mode("fast")
    FS.set_threads(args.threads)
    cases = args.case or list(CASES)
    runner = _run_outer if args.outer else _run_kernel
    payload = {
        "candidate": "canonical_fast_specialized",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "guard_contract": {
            "h": CANONICAL_H,
            "col_step": CANONICAL_COL_STEP,
            "z_tail": CANONICAL_Z_TAIL,
            "phase_max": CANONICAL_PHASE_MAX,
            "assemble": 1,
            "h_arr": None,
            "Sv": "present",
            "handoff_eps": "present",
        },
        "scope": "full_outer" if args.outer else "kernel",
        "cases": [runner(case, args.repeats, args.threads) for case in cases],
    }
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
