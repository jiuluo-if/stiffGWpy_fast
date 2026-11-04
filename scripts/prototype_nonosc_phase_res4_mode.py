"""Standalone Newton–Kantorovich phase-transfer probes for res4 modes.

中文：仅研究参考求解器的单模传播；保留生产 Cartesian 路径，不修改公开后端。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(
    Omega_bh2=0.0223828,
    Omega_ch2=0.1201075,
    H0=67.32117,
    A_s=2.100549e-9,
    r=0.01,
    cr=0.0,
    n_t=2.0,
    T_re=1e4,
    kappa10=0.01,
    DN_re=10.0,
)
DEFAULT_DN = 0.13373681817755848
DEFAULT_FREQ = -12.6723
DEFAULT_PHASE = 5.0
DN = DEFAULT_DN
FREQ = DEFAULT_FREQ
ZSTART = -12.0
ZPHASE = DEFAULT_PHASE
ZTAIL = 10.0
PHASE_GAMMA_MIN = 10.0
model = None
bg = None
NSTART = None
Z0 = None
NCUT = None
NCUT100 = None


def configure_case(dn_eff, freq, phase_start=DEFAULT_PHASE):
    global DN, FREQ, ZPHASE, model, bg, NSTART, Z0, NCUT, NCUT100
    DN = float(dn_eff)
    FREQ = float(freq)
    ZPHASE = float(phase_start)
    model = LCDM_SG(**PARAMS, DN_eff=DN)
    bg = R._Background(model, DN)
    NSTART = R._find_start_N(bg, FREQ, ZSTART)
    Z0 = (FREQ - R.f_hor_abs(bg, NSTART)) * math.log(10.0)
    NCUT = bg.N_inf - math.log(bg.nu_today / 0.1)
    NCUT100 = bg.N_inf - math.log(bg.nu_today / 100.0)


configure_case(DEFAULT_DN, DEFAULT_FREQ, DEFAULT_PHASE)


def p_and_prime(n, regime):
    E = math.exp(bg.N_inf - n)
    E3 = E**3
    inv = 1.0 / E3
    Om, On, Oph, Oer, Osh, OL = (bg.Omh2, bg.Omega_nh2, bg.Omega_ph2, bg.Oerh2, bg.Osh2, bg.OLh2)
    if regime == "radiation" or regime == "high_frequency":
        mass_term = Om
        radiation = bg.Otrh2
        if regime == "high_frequency":
            mass_term += bg.Omega_mnuh2
            radiation = Oph + 2.0 / 3.0 * On + Oer
        D = mass_term + radiation * E + Osh * E3 + OL * inv
        S = mass_term + 4.0 / 3.0 * radiation * E + 2.0 * Osh * E3
        Dp = -radiation * E - 3.0 * Osh * E3 + 3.0 * OL * inv
        Sp = -4.0 / 3.0 * radiation * E - 6.0 * Osh * E3
    else:
        nu = bg.nu_today / E
        ell = math.log10(nu)
        rho = float(R._FD_RHO_CS(ell))
        pressure = float(R._FD_P_CS(ell))
        rhop = float(R._FD_RHO_CS.derivative()(ell)) / math.log(10.0)
        pressurep = float(R._FD_P_CS.derivative()(ell)) / math.log(10.0)
        A = Oph + (2.0 / 3.0 + rho / 3.0) * On + Oer
        Ap = rhop * On / 3.0
        B = 4.0 / 3.0 * (Oph + 2.0 / 3.0 * On + Oer) + (rho + pressure) * On / 3.0
        Bp = (rhop + pressurep) * On / 3.0
        D = Om + A * E + Osh * E3 + OL * inv
        S = Om + B * E + 2.0 * Osh * E3
        Dp = E * (Ap - A) - 3.0 * Osh * E3 + 3.0 * OL * inv
        Sp = E * (Bp - B) - 6.0 * Osh * E3
    sigma = S / D
    sigmap = (Sp * D - S * Dp) / (D * D)
    return 3.0 - 1.5 * sigma, -1.5 * sigmap


def zprime(n, regime):
    return 2.0 - p_and_prime(n, regime)[0]


def dm(k):
    x = np.cos(np.pi * np.arange(k) / (k - 1))
    c = (-1.0) ** np.arange(k)
    c[[0, -1]] *= 2.0
    dx = x[:, None] - x[None, :]
    D = (c[:, None] / c[None, :]) / (dx + np.eye(k))
    D -= np.diag(D.sum(axis=1))
    return x, D


def phase_segment(a, b, regime, k, state, z_left):
    x, D = dm(k)
    half = (b - a) / 2.0
    N = (a + b) / 2.0 + half * x
    Dt = D / half
    zp = np.array([zprime(float(n), regime) for n in N])
    cz = C.chebfit(x, zp, k - 1)
    az = C.chebint(cz)
    z_nodes = z_left + half * (C.chebval(x, az) - C.chebval(-1.0, az))
    z_right = float(z_left + half * (C.chebval(1.0, az) - C.chebval(-1.0, az)))
    p = np.array([p_and_prime(float(n), regime)[0] for n in N])
    pp = np.array([p_and_prime(float(n), regime)[1] for n in N])
    Q = np.exp(2 * z_nodes) - 0.5 * pp - 0.25 * p * p
    if np.any(Q <= 0) or not np.all(np.isfinite(Q)):
        raise FloatingPointError("non-positive/non-finite normal-form Q")
    gamma_min = math.sqrt(float(Q.min())) * (b - a)
    if gamma_min < PHASE_GAMMA_MIN:
        raise ValueError(
            "phase interval is below the high-frequency screen: "
            "gamma=%.6g < %.6g; increase --phase-start or refine the interval" % (gamma_min, PHASE_GAMMA_MIN)
        )
    scale = math.sqrt(math.sqrt(float(Q.min() * Q.max())))
    q = Q / scale**2
    qp = Dt @ q
    r = 1j * scale * np.sqrt(q) - qp / (4.0 * q)
    iterations = 0
    newton_update_converged = False
    newton_residual_converged = False
    last_relative_update = math.inf
    for iterations in range(1, 21):
        resid = Dt @ r + r * r + scale**2 * q
        h = np.linalg.solve(Dt + np.diag(2.0 * r), -resid)
        r += h
        last_relative_update = float(np.max(np.abs(h)) / max(np.max(np.abs(r)), 1.0))
        if last_relative_update < 2e-14:
            newton_update_converged = True
            break
        resid = Dt @ r + r * r + scale**2 * q
        if float(np.max(np.abs(resid)) / scale**2) <= 1e-10:
            newton_residual_converged = True
            break
    resid = Dt @ r + r * r + scale**2 * q
    collocation_residual_scaled = float(np.max(np.abs(resid)) / scale**2)
    if collocation_residual_scaled > 1e-10:
        raise ArithmeticError(
            "collocated Riccati residual %.3g exceeds the prototype gate on %s branch "
            "N=[%.17g, %.17g], last relative update=%.3g after %d steps"
            % (collocation_residual_scaled, regime, a, b, last_relative_update, iterations)
        )
    cr = C.chebfit(x, r, k - 1)
    xd = np.linspace(-1.0, 1.0, 201)
    Nd = (a + b) / 2.0 + half * xd
    rd = C.chebval(xd, cr)
    rpd = C.chebval(xd, C.chebder(cr)) / half
    zd = z_left + half * (C.chebval(xd, az) - C.chebval(-1.0, az))
    pd = np.array([p_and_prime(float(n), regime)[0] for n in Nd])
    ppd = np.array([p_and_prime(float(n), regime)[1] for n in Nd])
    qd = np.exp(2 * zd) - 0.5 * ppd - 0.25 * pd * pd
    dense_res = rpd + rd * rd + qd
    dense_residual_scaled = float(np.max(np.abs(dense_res)) / scale**2)
    if dense_residual_scaled > 1e-9:
        raise ArithmeticError("dense Riccati residual %.3g exceeds the prototype gate" % dense_residual_scaled)
    anti = C.chebint(C.chebfit(x, r.imag, k - 1))
    alpha = half * (C.chebval(1.0, anti) - C.chebval(-1.0, anti))
    ap0, ap1 = float(r.imag[-1]), float(r.imag[0])
    if ap0 <= 0.0 or ap1 <= 0.0:
        raise FloatingPointError("phase derivative is not positive at interval endpoints")
    app0 = -2 * ap0 * float(r.real[-1])
    app1 = -2 * ap1 * float(r.real[0])

    def M(alpha, ap, app):
        c, s = math.cos(alpha), math.sin(alpha)
        rt = math.sqrt(ap)
        return np.array([[c / rt, s / rt], [-rt * s - 0.5 * app * c / (ap * rt), rt * c - 0.5 * app * s / (ap * rt)]])

    transfer = M(alpha, ap1, app1) @ np.linalg.inv(M(0.0, ap0, app0))
    p0 = p_and_prime(a, regime)[0]
    p1 = p_and_prime(b, regime)[0]
    pint = quad(lambda n: p_and_prime(n, regime)[0], a, b, epsabs=1e-13, epsrel=1e-13)[0]
    u0 = np.array([state[0], state[1] + 0.5 * p0 * state[0]])
    u1 = transfer @ u0
    out = math.exp(-pint / 2.0) * np.array([u1[0], u1[1] - 0.5 * p1 * u1[0]])
    if not np.all(np.isfinite(out)) or not math.isfinite(z_right):
        raise FloatingPointError("phase transfer returned a non-finite endpoint state")
    return (
        out,
        z_right,
        {
            "N": [a, b],
            "z": [z_left, z_right],
            "regime": regime,
            "order": k,
            "gamma_min": gamma_min,
            "newton_iterations": iterations,
            "newton_update_converged": newton_update_converged,
            "newton_residual_converged": newton_residual_converged,
            "newton_last_relative_update": last_relative_update,
            "collocation_residual_scaled": collocation_residual_scaled,
            "dense_residual_scaled": dense_residual_scaled,
        },
    )


def _z_between(z0, n0, n1, regime):
    return (
        z0
        + quad(
            lambda value: zprime(value, regime),
            n0,
            n1,
            epsabs=1e-13,
            epsrel=1e-13,
        )[0]
    )


def _invert_z(target, z0, n0, n1, regime):
    return brentq(
        lambda value: _z_between(z0, n0, value, regime) - target,
        n0,
        n1,
        xtol=5e-14,
        rtol=1e-14,
    )


def make_bounds(N5, z5):
    if N5 < NCUT:
        zcut = _z_between(z5, N5, NCUT, "radiation")
    else:
        zcut = z5

    spline_start = max(N5, NCUT)
    spline_zstart = zcut if N5 < NCUT else z5
    if spline_start < NCUT100:
        zcut100 = _z_between(spline_zstart, spline_start, NCUT100, "neutrino_spline")
    else:
        zcut100 = spline_zstart

    high_start = max(N5, NCUT100)
    high_zstart = zcut100 if NCUT100 > N5 else z5
    ztoday = _z_between(high_zstart, high_start, bg.N_inf, "high_frequency")
    started_at_tail = abs(z5 - ZTAIL) <= 1e-10
    if z5 > ZTAIL + 1e-10:
        raise ValueError("phase handoff starts beyond the tail threshold; lower --phase-start")
    tail_used = started_at_tail or ztoday >= ZTAIL

    def invert(target):
        if target < zcut - 1e-12 and N5 < NCUT:
            return _invert_z(target, z5, N5, NCUT, "radiation"), "radiation"
        if target < zcut100 - 1e-12 and spline_start < NCUT100:
            return _invert_z(target, spline_zstart, spline_start, NCUT100, "neutrino_spline"), "neutrino_spline"
        if target <= ztoday:
            return _invert_z(target, high_zstart, high_start, bg.N_inf, "high_frequency"), "high_frequency"
        raise ValueError("requested phase node lies beyond the integrated z range")

    if started_at_tail:
        nend, zend, end_regime = N5, ZTAIL, "phase_start"
    elif tail_used:
        nend, end_regime = invert(ZTAIL)
        zend = ZTAIL
    else:
        nend, zend, end_regime = bg.N_inf, ztoday, "high_frequency"

    bounds = [(N5, z5, "")]
    for target in (6.0, 7.0, 8.0, 9.0):
        if target >= zend - 1e-12 or abs(target - zcut) < 1e-12 or abs(target - zcut100) < 1e-12:
            continue
        n, regime = invert(target)
        bounds.append((n, target, regime))

    if N5 < NCUT < nend:
        bounds.append((NCUT, zcut, "radiation"))
    if N5 < NCUT100 < nend:
        bounds.append((NCUT100, zcut100, "neutrino_spline"))
    bounds.append((nend, zend, end_regime))
    bounds.sort(key=lambda row: row[0])

    deduplicated = []
    for index, bound in enumerate(bounds):
        if not deduplicated or bound[0] - deduplicated[-1][0] > 1e-13:
            deduplicated.append(bound)
        elif index == len(bounds) - 1:
            # 中文：终点即使与分段点几乎重合，也必须覆盖保留为最后边界。
            deduplicated[-1] = bound
    return deduplicated, tail_used, zend


def candidate(k, low_state, N5, z5):
    bounds, tail_used, zend = make_bounds(N5, z5)
    state = np.array(low_state[1:3], dtype=float)
    z = z5
    rows = []
    t = time.perf_counter()
    for (a, _, _), (b, zb, _) in zip(bounds[:-1], bounds[1:]):
        midpoint = 0.5 * (a + b)
        if midpoint < NCUT:
            regime = "radiation"
        elif midpoint < NCUT100:
            regime = "neutrino_spline"
        else:
            regime = "high_frequency"
        state, z, row = phase_segment(a, b, regime, k, state, z)
        row["z_target_right"] = zb
        rows.append(row)
    return state, bounds, rows, time.perf_counter() - t, tail_used, zend


def spectrum_today(state, z):
    T, Tp = state
    x, y = Tp, math.exp(z) * T
    pt = model.derived_param["A_t"] * (10.0**FREQ / gp.f_piv) ** model.derived_param["nt"]
    oj = x * T / 3.0 * pt
    ogw = (x * x + y * y) / 24.0 * pt + oj
    opgw = (-5.0 * x * x + 7.0 * y * y) / 72.0 * pt
    return {"OmegaGW": ogw - oj, "Ogw": ogw, "Oj": oj, "Opgw": opgw}


def spectrum_from_tail(state, Ntail, z):
    T, Tp = state
    omega = math.exp(z)
    x = Tp
    y = omega * T
    coeff = math.sqrt(0.5 * (x * x + y * y))
    transfer = coeff * math.exp(-z + Ntail - bg.N_inf)
    xtoday = transfer * 10.0 ** (FREQ - R.f_hor_abs(bg, bg.N_inf))
    pt = model.derived_param["A_t"] * (10.0**FREQ / gp.f_piv) ** model.derived_param["nt"]
    oj = -transfer * transfer / 3.0 * pt
    op = xtoday * xtoday / 36.0 * pt
    ogw = 3.0 * op + oj
    return {"OmegaGW": ogw - oj, "Ogw": ogw, "Oj": oj, "Opgw": op}


# 生产方程的 z、T、T' 在低频段共同积分，并以 z=5 事件交接。
def low_rhs(n, s):
    z, T, Tp = s
    sig = R._H2_and_sigma(bg, float(n))[1]
    p = 3.0 - 1.5 * sig
    return [1.5 * sig - 1.0, Tp, -p * Tp - math.exp(2 * z) * T]


def low_event(n, s):
    return s[0] - ZPHASE


low_event.terminal = True
low_event.direction = 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--omega-bh2", type=float, default=PARAMS["Omega_bh2"])
    parser.add_argument("--omega-ch2", type=float, default=PARAMS["Omega_ch2"])
    parser.add_argument("--h0", type=float, default=PARAMS["H0"])
    parser.add_argument("--a-s", type=float, default=PARAMS["A_s"])
    parser.add_argument("--r", type=float, default=PARAMS["r"])
    parser.add_argument("--cr", type=float, default=PARAMS["cr"])
    parser.add_argument("--n-t", type=float, default=PARAMS["n_t"])
    parser.add_argument("--t-re", type=float, default=PARAMS["T_re"])
    parser.add_argument("--kappa10", type=float, default=PARAMS["kappa10"])
    parser.add_argument("--dn-re", type=float, default=PARAMS["DN_re"])
    parser.add_argument("--dn-eff", type=float, default=DEFAULT_DN)
    parser.add_argument("--frequency", type=float, default=DEFAULT_FREQ)
    parser.add_argument("--phase-start", type=float, default=DEFAULT_PHASE)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    PARAMS.update(
        Omega_bh2=args.omega_bh2,
        Omega_ch2=args.omega_ch2,
        H0=args.h0,
        A_s=args.a_s,
        r=args.r,
        cr=args.cr,
        n_t=args.n_t,
        T_re=args.t_re,
        kappa10=args.kappa10,
        DN_re=args.dn_re,
    )
    configure_case(args.dn_eff, args.frequency, args.phase_start)
    if args.output is None:
        if FREQ == DEFAULT_FREQ and DN == DEFAULT_DN and ZPHASE == DEFAULT_PHASE:
            args.output = ROOT / "docs" / "nonosc_phase_res4_mode_20261104.json"
        else:
            args.output = ROOT / "docs" / ("nonosc_phase_res4_mode_freq_%+.4f_20261104.json" % FREQ)
    low_start = time.perf_counter()
    lo = solve_ivp(
        low_rhs,
        (NSTART, bg.N_inf),
        (Z0, 1.0, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=[1e-12, 1e-13, 1e-13],
        events=[low_event],
        dense_output=False,
    )
    if not lo.success or not lo.t_events[0].size:
        raise RuntimeError("low-frequency production-z propagation did not reach phase handoff")
    low_seconds = time.perf_counter() - low_start
    N5 = float(lo.t_events[0][0])
    low_state = lo.y_events[0][0].tolist()
    results = []
    for k in (16, 24, 32):
        candidate_start = time.perf_counter()
        state, bounds, rows, phase_s, tail_used, zend = candidate(k, low_state, N5, ZPHASE)
        Nend = bounds[-1][0]
        power = spectrum_from_tail(state, Nend, zend) if tail_used else spectrum_today(state, zend)
        candidate_s = time.perf_counter() - candidate_start + low_seconds
        results.append(
            {
                "order": k,
                "segments": len(rows),
                "state_T_Tprime": state.tolist(),
                "tail_used": tail_used,
                "end_N": Nend,
                "end_z": zend,
                "end_state_spectrum": power,
                "phase_seconds": phase_s,
                "total_candidate_seconds": candidate_s,
                "low_rhs_nfev": lo.nfev,
                "low_seconds": low_seconds,
                "segments_details": rows,
                "gamma_min": min((x["gamma_min"] for x in rows), default=None),
                "dense_residual_max": max(
                    (x["dense_residual_scaled"] for x in rows), default=None
                ),
            }
        )

    # 中文：在相同模态、DN 和尾部约定下运行当前 Cartesian Reference 作为独立对照。
    cartesian_spline_warmup_start = time.perf_counter()
    R._warmup_numba_spline()
    cartesian_spline_warmup_s = time.perf_counter() - cartesian_spline_warmup_start
    t = time.perf_counter()
    ref = R.solve_reference_mode(model, FREQ, DN, z_tail=ZTAIL, rtol=1e-11, backend="cartesian")
    ref_s = time.perf_counter() - t
    refpower = float(ref["Ogw_today"] - ref["Oj_today"])
    ref_log = math.log10(refpower)
    for row in results:
        row["spectrum_abs_dex_error_vs_cartesian"] = abs(math.log10(row["end_state_spectrum"]["OmegaGW"]) - ref_log)
        row["runtime_ratio_cartesian_over_candidate"] = ref_s / row["total_candidate_seconds"]
    if bool(ref["used_tail"]) != bool(results[0]["tail_used"]):
        raise ArithmeticError("phase and Cartesian paths disagree on tail usage")
    source = ROOT / "stiffgwpy_fast" / "reference.py"
    payload = {
        "schema_version": 1,
        "experiment": "newton_kantorovich_phase_transfer_production_z_single_mode",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "prototype_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "params": {**PARAMS, "DN_eff": DN},
        "freq_log10_hz": FREQ,
        "z_start": ZSTART,
        "z_phase_start": ZPHASE,
        "z_tail": ZTAIL,
        "solver_settings": {
            "low_frequency_method": "DOP853",
            "low_frequency_rtol": 1e-12,
            "low_frequency_atol": [1e-12, 1e-13, 1e-13],
            "phase_chebyshev_orders": [16, 24, 32],
            "phase_newton_relative_update_tolerance": 2e-14,
            "phase_newton_max_iterations": 20,
            "phase_high_frequency_gamma_min": PHASE_GAMMA_MIN,
            "phase_dense_residual_points_per_segment": 201,
            "background_segment_splits": ["nu=0.1", "nu=100"],
            "background_branch_derivatives": "analytic, branch-specific p prime",
            "candidate_spline_callback": "Python scalar CubicSpline fallback",
            "cartesian_reference_backend": "cartesian",
            "cartesian_reference_rtol": 1e-11,
            "cartesian_reference_atol": [1e-12, 1e-22, 1e-22],
        },
        "N_start": NSTART,
        "z0": Z0,
        "N_phase_start": N5,
        "low_state_z_T_Tprime": low_state,
        "nu_branch_boundary_N": [NCUT, NCUT100],
        "tail_used": results[0]["tail_used"],
        "N_end_production_z": results[0]["end_N"],
        "z_end_production": results[0]["end_z"],
        "reference_cartesian": {
            "spline_jit_warmup_s": cartesian_spline_warmup_s,
            "runtime_s": ref_s,
            "nfev": ref["n_steps"],
            "log10OmegaGW": ref_log,
            "event_N": ref["event_N"],
            "used_tail": ref["used_tail"],
        },
        "runs": results,
        "resources": telemetry(workers=1, threads=1),
        "limitations": [
            "Standalone single-mode research prototype; not a production backend.",
            "Only one res4 parameter point and three frequencies are checked.",
            "No full-grid integral map, self-consistency root, or parameter-domain validation is included.",
            "The method uses a configurable DOP853 prefix and local Chebyshev phase transfers on smooth positive-Q branches.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
