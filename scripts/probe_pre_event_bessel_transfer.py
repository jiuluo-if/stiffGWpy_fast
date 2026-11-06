"""Feasibility screen for exact propagation on the sigma=1 pre-event segment.

This standalone research probe compares the Fast midpoint transfer against
the closed-form order-3/2 Bessel solution at the same coarse output nodes.
It does not modify the production solver or claim an end-to-end speedup.
"""

from __future__ import annotations

import hashlib
import json
import math
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

from numba import njit  # noqa: E402

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.freq_adaptive import goal_oriented_freqs  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

OUTPUT = ROOT / "docs" / "pre_event_bessel_transfer_20261106.json"
H = 0.005
TAIL = 5.0
PHASE_MAX = 0.25
COL_STEP = 8


@njit(cache=True, inline="always")
def _basis(t):
    """Return the basis and N-derivatives for T''+3/2 T'+omega^2 T=0."""
    if t <= 0.0:
        return 0.0, 0.0, 0.0, 0.0
    if t < 0.2:
        t2 = t * t
        t4 = t2 * t2
        t6 = t4 * t2
        t8 = t4 * t4
        t10 = t8 * t2
        r = 1.0 - t2 / 10.0 + t4 / 280.0 - t6 / 15120.0 + t8 / 1330560.0 - t10 / 172972800.0
        rn = -t2 / 10.0 + t4 / 140.0 - t6 / 5040.0 + t8 / 332640.0 - t10 / 34594560.0
        s = -1.0 / (t * t * t) - 1.0 / (2.0 * t) + t / 8.0 - t * t2 / 144.0 + t * t4 / 5760.0 - t * t6 / 403200.0 + t * t8 / 43545600.0
        sn = 1.5 / (t * t * t) + 1.0 / (4.0 * t) + t / 16.0 - t * t2 / 96.0 + t * t4 / 2304.0 - t * t6 / 115200.0 + t * t8 / 9676800.0
        return r, s, rn, sn
    st = math.sin(t)
    ct = math.cos(t)
    f = st - t * ct
    ft = t * st
    g = -(ct + t * st)
    gt = -t * ct
    inv_t2 = 1.0 / (t * t)
    inv_t3 = inv_t2 / t
    r = 3.0 * f * inv_t3
    rt = 3.0 * (ft * inv_t3 - 3.0 * f * inv_t3 / t)
    s = g * inv_t3
    st_derivative = gt * inv_t3 - 3.0 * g * inv_t3 / t
    return r, s, 0.5 * t * rt, 0.5 * t * st_derivative


@njit(cache=True)
def _bessel_coefficients(t0, T0, Tp0):
    r0, s0, rn0, sn0 = _basis(t0)
    determinant = 1.5 / (t0 * t0 * t0)
    c0 = (T0 * sn0 - s0 * Tp0) / determinant
    c1 = (r0 * Tp0 - rn0 * T0) / determinant
    return c0, c1


@njit(cache=True)
def _exact_fast_state_at(N, N0, z0, c0, c1):
    z = z0 + 0.5 * (N - N0)
    t = 2.0 * math.exp(z)
    r, s, rn, sn = _basis(t)
    T = c0 * r + c1 * s
    Tp = c0 * rn + c1 * sn
    scale = math.exp(0.5 * N)
    return scale * Tp, scale * math.exp(z) * T


@njit(cache=True)
def _baseline_prefix(Nv, targets, j0, z0, yh0, phase_max):
    out = np.empty((targets.size, 2), dtype=np.float64)
    xh = 0.0
    yh = yh0
    k = j0
    current_n = Nv[j0]
    for row in range(targets.size):
        target = targets[row]
        while current_n < target:
            end = Nv[k + 1]
            if end > target:
                end = target
            h_step = end - current_n
            z = z0 + 0.5 * (current_n - Nv[j0])
            z_end = z0 + 0.5 * (end - Nv[j0])
            xh, yh = FS._phase_segment(xh, yh, z, z_end, h_step, phase_max)
            current_n = end
            if end >= Nv[k + 1]:
                k += 1
        out[row, 0] = xh
        out[row, 1] = yh
    return out


@njit(cache=True)
def _closed_form_prefix(targets, N0, z0, xh0, yh0):
    T0 = math.exp(-0.5 * N0 - z0) * yh0
    Tp0 = math.exp(-0.5 * N0) * xh0
    t0 = 2.0 * math.exp(z0)
    c0, c1 = _bessel_coefficients(t0, T0, Tp0)
    out = np.empty((targets.size, 2), dtype=np.float64)
    for i in range(targets.size):
        out[i, 0], out[i, 1] = _exact_fast_state_at(
            targets[i], N0, z0, c0, c1
        )
    return out


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _profile(fn, repeats=31):
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return {"p50_s": statistics.median(times), "p95_s": float(np.percentile(times, 95))}


def main() -> int:
    # The small-t expansion switches to the trigonometric form at t=0.2.
    # Check derivative continuity at that boundary before any timing run.
    series_edge = _basis(0.2 - 1e-12)
    trig_edge = _basis(0.2)
    if abs(series_edge[2] - trig_edge[2]) > 1e-12:
        raise AssertionError(
            f"R_N series/trigonometric boundary mismatch: "
            f"{series_edge[2]:.17g} vs {trig_edge[2]:.17g}"
        )
    FS.apply_accuracy_mode("fast")
    FS.set_threads(1)
    model = LCDM_SG(r=0.01, cr=1, T_re=2000.0, kappa10=0.01)
    FS.gen_fast(model, H, kink_split=True)
    freqs = goal_oriented_freqs(model, 1.0, seed_n=64, max_points=120, eval_freqs=None)
    _sigma, _f_hor, j0s, z0s, _fp_minus = FS.prep_frequency_only(
        model, model.Nv, freqs
    )
    n_re = float(model.derived_param["N_inf"] - model.derived_param["N_re"])
    event_idx = int(np.searchsorted(model.Nv, n_re))
    pre_modes = [i for i, j0 in enumerate(j0s) if j0 < event_idx and z0s[i] < TAIL]
    if len(pre_modes) < 8:
        raise RuntimeError(f"screen too small: only {len(pre_modes)} pre-event modes")

    records = []
    errors = []
    speed_ratios = []
    for mode in (pre_modes[0], pre_modes[len(pre_modes) // 2], pre_modes[-1]):
        j0 = int(j0s[mode])
        z0 = float(z0s[mode])
        N0 = float(model.Nv[j0])
        stop_N = min(n_re, N0 + 2.0 * (TAIL - z0))
        target_indices = [
            k for k in range(j0 + 1, event_idx + 1)
            if (k % COL_STEP == 0 and model.Nv[k] <= stop_N)
        ]
        targets = [float(model.Nv[k]) for k in target_indices]
        if not targets or targets[-1] < stop_N:
            targets.append(stop_N)
        targets_arr = np.asarray(targets, dtype=np.float64)
        S2inv0 = math.exp(0.5 * N0)
        xh0 = 0.0
        yh0 = math.exp(z0) * S2inv0
        def baseline_fn():
            return _baseline_prefix(model.Nv, targets_arr, j0, z0, yh0, PHASE_MAX)

        def candidate_fn():
            return _closed_form_prefix(targets_arr, N0, z0, xh0, yh0)
        baseline = baseline_fn()
        candidate = candidate_fn()
        scale = np.maximum(np.linalg.norm(baseline, axis=1), 1e-300)
        state_error = np.linalg.norm(candidate - baseline, axis=1) / scale
        errors.extend(state_error.tolist())
        baseline_timing = _profile(baseline_fn)
        candidate_timing = _profile(candidate_fn)
        ratio = candidate_timing["p50_s"] / baseline_timing["p50_s"]
        speed_ratios.append(ratio)

        # Independent tight ODE checks at the final prefix point.
        def rhs(N, state):
            z = z0 + 0.5 * (N - N0)
            w = math.exp(z)
            x, y = state
            return [-x - w * y, w * x + y]

        reference = solve_ivp(
            rhs, (N0, stop_N), np.array([xh0, yh0]), method="DOP853",
            rtol=1e-11, atol=1e-13,
        )
        independent_error = float(
            np.linalg.norm(candidate[-1] - reference.y[:, -1])
            / max(np.linalg.norm(reference.y[:, -1]), 1e-300)
        )
        if independent_error > 1e-8:
            raise AssertionError(f"Bessel/DOP853 mismatch {independent_error:.3e}")
        records.append({
            "mode": int(mode), "frequency_log10_hz": float(freqs[mode]),
            "j0": j0, "N0": N0, "z0": z0, "N_re_abs": n_re,
            "prefix_stop_N": stop_N, "coarse_targets": len(targets),
            "baseline_vs_closed_form_state_max_rel": float(np.max(state_error)),
            "closed_form_vs_dop853_final_rel": independent_error,
            "baseline_timing": baseline_timing,
            "closed_form_timing": candidate_timing,
            "closed_form_over_baseline_p50": ratio,
            "limitations": [
                "single-mode prefix microbenchmark, not full solve or end-to-end speed",
                "does not yet splice the exact state into post-event kernel or assemble all history columns",
                "baseline is the Fast midpoint transfer, not an independent self-consistent root",
            ],
        })

    record = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "exact order-3/2 Bessel propagation over the sigma=1 prefix",
        "source_commit": "f1568eb1397cf53cd245685dfe567eeb9ff09465",
        "branch": "fast_v0.2",
        "parameters": dict(model.cosmo_param),
        "N_inf": float(model.derived_param["N_inf"]),
        "N_re_abs": n_re,
        "background_nodes": int(model.Nv.size),
        "frequency_nodes": int(freqs.size),
        "pre_event_modes": int(len(pre_modes)),
        "threads": 1,
        "per_mode": records,
        "summary": {
            "max_baseline_state_relative_difference": float(max(errors)),
            "median_candidate_over_baseline_prefix_p50": float(statistics.median(speed_ratios)),
        },
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {
            "fast_sgwb.py": _hash(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "probe_script": _hash(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(record["summary"], indent=2))
    print(f"artifact={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
