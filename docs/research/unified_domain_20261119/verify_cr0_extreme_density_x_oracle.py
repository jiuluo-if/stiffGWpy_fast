"""Independent Prüfer check in look-back coordinates for extreme CR0 density.

This research verifier evaluates the unchanged tensor Prüfer equations using
``x = N_inf - N`` and a stable H2 closure identity. The coordinate avoids
subtracting nearly equal absolute N values at the present-day closure. It does
not modify or import a modified Reference implementation.
"""
from __future__ import annotations

import hashlib
import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _background(model: LCDM_SG, dn_eff: float, x: float) -> tuple[float, float]:
    d = model.derived_param
    om = float(d["Omega_mh2"])
    osh = float(d["Omega_sh2"])
    anchor = float(d["h"] ** 2 - gp.Omega_mnuh2
                   - gp.Omega_nh2 * 2.0 / 3.0 - gp.Omega_ph2)
    oer = gp.Omega_ph2 * 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0) * dn_eff
    e = math.exp(x)
    q = math.exp(-3.0 * x)
    e3 = 1.0 / q
    base = om * (-math.expm1(-3.0 * x)) + 2.0 * osh * math.sinh(3.0 * x) + anchor * q
    extra = oer * (e - q)
    nu = gp.nu_today / e
    if nu > 100.0:
        h2 = (base + gp.Omega_mnuh2
              + (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2) * e + extra)
        numerator = (om + gp.Omega_mnuh2
                     + 4.0 / 3.0 * (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2 + oer) * e
                     + 2.0 * osh * e3)
    elif nu >= 0.1:
        rho = float(REF._FD_RHO_CS(math.log10(nu)))
        pressure = float(REF._FD_P_CS(math.log10(nu)))
        h2 = (base + (gp.Omega_ph2 + (2.0 / 3.0 + rho / 3.0) * gp.Omega_nh2) * e + extra)
        numerator = (om + 4.0 / 3.0 * (gp.Omega_ph2 + 2.0 / 3.0 * gp.Omega_nh2 + oer) * e
                     + (rho + pressure) * gp.Omega_nh2 / 3.0 * e + 2.0 * osh * e3)
    elif x < gp.N_fin:
        otr = gp.Omega_orh2 + oer
        h2 = base + otr * e
        numerator = om + 4.0 / 3.0 * otr * e + 2.0 * osh * e3
    elif x <= gp.N_max:
        rho = float(gp.spl_rho(x))
        rhop = float(gp.spl_rhop(x))
        h2 = base + (gp.Omega_ph2 * rho + oer) * e
        numerator = om + (gp.Omega_ph2 * rhop + 4.0 / 3.0 * oer) * e + 2.0 * osh * e3
    else:
        otreh = gp.Omega_ph2 * gp.rho_th[-1] + oer
        h2 = base + otreh * e
        numerator = om + 4.0 / 3.0 * otreh * e + 2.0 * osh * e3
    if h2 <= 0.0 or not math.isfinite(h2):
        raise FloatingPointError(f"nonpositive or nonfinite H2 at x={x}: {h2}")
    return h2, numerator / h2


def _fh(model: LCDM_SG, dn_eff: float, x: float) -> float:
    d = model.derived_param
    nre = float(d["N_re"])
    h2today = float(d["h"] ** 2)
    if x > nre:
        h2re, _ = _background(model, dn_eff, nre)
        raw_delta = 0.5 * x + 0.5 * math.log(h2re / h2today)
    else:
        h2, _ = _background(model, dn_eff, x)
        raw_delta = 0.5 * x + 0.5 * math.log(h2 / h2today)
    return (raw_delta - math.log(2.0 * math.pi / float(d["H_0"]))) / math.log(10.0)


def _h2_coordinate(model: LCDM_SG, dn_eff: float, x: float) -> float:
    """Return the constant H2 used by the analytic pre-reheating horizon extension."""
    nre = float(model.derived_param["N_re"])
    if x <= nre:
        return _background(model, dn_eff, x)[0]
    return _background(model, dn_eff, nre)[0]


def _solve_mode(model: LCDM_SG, freq: float, dn_eff: float, z_tail: float = 8.0) -> dict:
    ninf = float(model.derived_param["N_inf"])
    nre = float(model.derived_param["N_re"])
    z_start = -12.0

    def z_at(x: float) -> float:
        return (freq - _fh(model, dn_eff, x)) * math.log(10.0)

    if z_at(ninf) >= z_start:
        xstart = ninf
    else:
        xstart = brentq(lambda x: z_at(x) - z_start, 0.0, ninf, xtol=2e-13)
    z0 = z_at(xstart)

    def rhs(x: float, state: np.ndarray) -> tuple[float, float]:
        # The background identity is dlog(H2)/dN=3(1-sigma). Thus for
        # b=log(A)+1/2 log(H2)+N/2, d b/dx=-cos(2 theta), exactly.
        phase = float(state[1])
        z = z_at(x)
        return -math.cos(2.0 * phase), math.exp(z) + math.sin(2.0 * phase)

    def event(x: float, _state: np.ndarray) -> float:
        return float(z_at(x) - z_tail)

    event.terminal = True
    event.direction = 1.0
    ninf = float(model.derived_param["N_inf"])
    b0 = z0 + 0.5 * math.log(_h2_coordinate(model, dn_eff, xstart)) \
        + 0.5 * (ninf - xstart)
    state = np.array([b0, 0.0])
    pieces = [(xstart, nre), (nre, 0.0)] if xstart > nre > 0.0 else [(xstart, 0.0)]
    result = None
    x_event = None
    for left, right in pieces:
        if left <= right:
            continue
        result = solve_ivp(rhs, (left, right), state, method="DOP853",
                           rtol=1e-10, atol=1e-12, events=event)
        if not result.success:
            raise RuntimeError(result.message)
        if result.t_events[0].size:
            x_event = float(result.t_events[0][0])
            state = result.y_events[0][0]
            break
        state = result.y[:, -1]
    if result is None:
        raise RuntimeError("mode integration had no interval")
    phase = float(state[1])
    zf = z_at(x_event if x_event is not None else 0.0)
    pt = float(model.derived_param["A_t"]) * (10.0 ** freq / gp.f_piv) ** float(model.derived_param["nt"])
    if x_event is not None:
        logamp = float(state[0]) - 0.5 * math.log(_h2_coordinate(model, dn_eff, x_event)) \
            - 0.5 * (ninf - x_event)
        log_transfer = logamp - 0.5 * math.log(2.0) - zf - x_event
        transfer = math.exp(log_transfer)
        f_today = -math.log(2.0 * math.pi / float(model.derived_param["H_0"])) / math.log(10.0)
        x_today = transfer * 10.0 ** (freq - f_today)
        oj = -transfer * transfer * pt / 3.0
        opgw = x_today * x_today * pt / 36.0
        ogw = 3.0 * opgw + oj
    else:
        logamp = float(state[0]) - 0.5 * math.log(_h2_coordinate(model, dn_eff, 0.0)) \
            - 0.5 * ninf
        amp = math.exp(logamp)
        xx = amp * math.sin(phase)
        yy = amp * math.cos(phase)
        th = yy / math.exp(zf)
        oj = xx * th * pt / 3.0
        ogw = (xx * xx + yy * yy) * pt / 24.0 + oj
        opgw = (-5.0 * xx * xx + 7.0 * yy * yy) * pt / 72.0
    if not all(math.isfinite(v) for v in (ogw, oj, opgw)):
        raise FloatingPointError("nonfinite oracle output")
    return {"Ogw": ogw, "Oj": oj, "Opgw": opgw,
            "nfev": int(result.nfev), "tail": x_event is not None}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--Omega-bh2", type=float, default=1.0e16)
    parser.add_argument("--H0", type=float, default=67.32117)
    parser.add_argument("--T-re", type=float, default=2.0e3)
    parser.add_argument("--n-t", type=float, default=0.0)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    model = LCDM_SG(Omega_bh2=args.Omega_bh2, Omega_ch2=0.1201075, H0=args.H0,
                    DN_eff=0.0, A_s=2.100549e-9, r=1.0e-16, n_t=args.n_t,
                    cr=0.0, T_re=args.T_re, DN_re=10.0, kappa10=1.0e-2)
    config = FS.resolve_config("fast", threads=2, freq_grid="adaptive")
    input_dn_eff = float(model.cosmo_param["DN_eff"])
    started = time.perf_counter()
    FS.SGWB_iter_fast(model, config=config, tol=1e-6,
                      freq_grid_target=1e-6, physical_guard=False)
    fast_seconds = time.perf_counter() - started
    freqs = np.asarray(model.f, dtype=float).copy()
    fast_spectrum = np.asarray(model.log10OmegaGW, dtype=float).copy()
    fast_root = float(np.asarray(model.DN_gw, dtype=float).reshape(-1)[-1])
    dn_eff = input_dn_eff + fast_root
    rows = []
    started = time.perf_counter()
    for i, freq in enumerate(freqs):
        row = _solve_mode(model, float(freq), dn_eff)
        rows.append(row)
        if (i + 1) % 25 == 0 or i + 1 == len(freqs):
            print(f"oracle {i + 1}/{len(freqs)}", flush=True)
    oracle_seconds = time.perf_counter() - started
    ogw = np.asarray([row["Ogw"] for row in rows])
    oj = np.asarray([row["Oj"] for row in rows])
    opgw = np.asarray([row["Opgw"] for row in rows])
    g2, qerr, ierr = REF.integrate_spectrum(freqs, ogw, oj)
    omega_nu = gp.Omega_nh2 / float(model.derived_param["h"] ** 2)
    first_map_root = float(gp.Neff0 * g2 / omega_nu)
    first_map_seconds = oracle_seconds

    # Re-evaluate the self-consistency map at its first independent root
    # estimate. This measures map closure without relying on production Reference.
    dn_eff_closed = input_dn_eff + first_map_root
    rows_closed = []
    started = time.perf_counter()
    for i, freq in enumerate(freqs):
        row = _solve_mode(model, float(freq), dn_eff_closed)
        rows_closed.append(row)
        if (i + 1) % 25 == 0 or i + 1 == len(freqs):
            print(f"closure oracle {i + 1}/{len(freqs)}", flush=True)
    closure_seconds = time.perf_counter() - started
    ogw = np.asarray([row["Ogw"] for row in rows_closed])
    oj = np.asarray([row["Oj"] for row in rows_closed])
    opgw = np.asarray([row["Opgw"] for row in rows_closed])
    g2, qerr, ierr = REF.integrate_spectrum(freqs, ogw, oj)
    oracle_root = float(gp.Neff0 * g2 / omega_nu)
    oracle_spectrum = np.log10(np.maximum(ogw - oj, np.finfo(float).tiny))
    diff = np.abs(fast_spectrum - oracle_spectrum)
    result = {
        "status": "complete",
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {name: _sha(ROOT / "stiffgwpy_fast" / name)
                          for name in ("exact_background.py", "fast_sgwb.py", "reference.py")},
        "inputs": {"Omega_bh2": args.Omega_bh2, "Omega_ch2": 0.1201075, "H0": args.H0,
                   "DN_eff": 0.0, "r": 1.0e-16, "n_t": args.n_t, "cr": 0.0,
                   "T_re": args.T_re, "DN_re": 10.0, "kappa10": 0.01},
        "fast": {"root": fast_root, "seconds": fast_seconds, "frequency_count": len(freqs)},
        "oracle": {"root_map": oracle_root, "seconds": oracle_seconds,
                   "first_map_root": first_map_root,
                   "closure_map_seconds": closure_seconds,
                   "closure_root_delta_rel": abs(oracle_root-first_map_root)/max(abs(oracle_root), 1e-300),
                   "nfev_total": sum(row["nfev"] for row in rows_closed),
                   "tail_fraction": float(np.mean([row["tail"] for row in rows_closed])),
                   "quadrature_error": float(qerr), "interpolation_error": ierr},
        "comparison": {"root_relative_error": abs(oracle_root-fast_root)/max(abs(fast_root), 1e-300),
                       "spectrum_max_error_dex": float(np.max(diff)),
                       "spectrum_p95_error_dex": float(np.quantile(diff, 0.95)),
                       "worst_frequency_index": int(np.argmax(diff)),
                       "worst_frequency": float(freqs[int(np.argmax(diff))])},
        "raw": {"freqs": freqs.tolist(), "fast_log10OmegaGW": fast_spectrum.tolist(),
                "oracle_Ogw": ogw.tolist(), "oracle_Oj": oj.tolist(),
                "oracle_Opgw": opgw.tolist(), "oracle_log10OmegaGW": oracle_spectrum.tolist(),
                "per_mode_nfev": [row["nfev"] for row in rows_closed],
                "per_mode_tail": [row["tail"] for row in rows_closed]},
    }
    out = args.output or (ROOT / "docs/research/unified_domain_20261119/cr0_extreme_density_x_oracle.json")
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("ORACLE_RESULT", json.dumps({k: result[k] for k in ("fast", "oracle", "comparison")}, allow_nan=False))


if __name__ == "__main__":
    main()
