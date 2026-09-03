# -*- coding: utf-8 -*-
"""stiffGWpy vs SageNet+ SGWB emulator: accuracy / error / runtime comparison.

Runs inside the dedicated venv ``F:\\codex\\sagenet_env`` (numpy 1.26 + torch
2.3.1 + astropy 6 + scipy 1.14 + numba 0.60) so the PyTorch DLLs and stiffgwpy
coexist (see the DLL-search preamble below).  For each physical parameter point:

* ``reference`` (stiffgwpy.reference, continuous-sigma DOP853) -- THE truth anchor
* ``production`` (stiffgwpy fast, transition-refine)
* ``plain_grid`` (stiffgwpy fast, plain-grid)
* ``sagenet`` (SageNet+ GWPredictor, one NN model per model_type)

Errors are computed against ``reference`` using (a) the stiffgwpy signal-region
relative/dex metrics and (b) SageNet's own ``metrics`` (area-difference,
SMAPE, per-point dex).  All times are single-point wall-clock seconds.

Usage (run with the venv python):
  F:\\codex\\sagenet_env\\Scripts\\python.exe run_compare.py --phase quick
  F:\\codex\\sagenet_env\\Scripts\\python.exe run_compare.py --phase reference --points 0 1
  F:\\codex\\sagenet_env\\Scripts\\python.exe run_compare.py --phase report
"""

import argparse
import json
import math
import os
import sys
import time

# --- DLL search preamble (torch on Windows needs MKL + torch/lib) ---------------
_TORCH_LIB = r"F:\codex\sagenet_env\Lib\site-packages\torch\lib"
_MKL_BIN = r"C:\miniconda3\Library\bin"
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "True")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ["PATH"] = _TORCH_LIB + os.pathsep + _MKL_BIN + os.pathsep + os.environ.get("PATH", "")
for _p in (_TORCH_LIB, _MKL_BIN):
    try:
        os.add_dll_directory(_p)
    except Exception:
        pass

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAGE_SRC = os.path.join(os.environ.get("TEMP", r"C:\Users\Public"),
                        "sagenet_src2")
if not os.path.isdir(SAGE_SRC):
    SAGE_SRC = r"C:\Users\联想\AppData\Local\Temp\sagenet_src2"
for _p in (REPO, SAGE_SRC):
    sys.path.insert(0, _p)

import numpy as np
from stiffgwpy import fast_sgwb as FS
from stiffgwpy import LCDM_SG
from stiffgwpy import reference as REF
from stiffgwpy.freq_adaptive import grid_independent_freqs

ln10 = math.log(10.0)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError("not JSON serializable: %s" % type(o))


# Representative points inside the SageNet+ valid parameter box (kept clear of
# the stiffgwpy shared-Delta_Neff>5 guard where possible).  ``A_s`` is the
# linear amplitude used by both; SageNet warns outside its box but the box is
# wide enough to cover all of these.
POINTS = [
    dict(label="sage_center", r=3.9585109e-05, n_t=1.0116972, kappa10=110.42477,
         T_re=0.17453859, DN_re=39.366618, Omega_bh2=0.0223828,
         Omega_ch2=0.1201075, H0=67.32117, A_s=2.100549e-9),
    dict(label="neutral", r=1.0e-3, n_t=0.0, kappa10=1.0,
         T_re=1.0e3, DN_re=20.0, Omega_bh2=0.0223828,
         Omega_ch2=0.1201075, H0=67.32117, A_s=2.100549e-9),
    dict(label="red_tilt", r=1.0e-1, n_t=-0.5, kappa10=5.0,
         T_re=5.0e2, DN_re=10.0, Omega_bh2=0.0223828,
         Omega_ch2=0.1201075, H0=67.32117, A_s=2.100549e-9),
    dict(label="high_Tre", r=1.0e-2, n_t=0.2, kappa10=0.05,
         T_re=5.0e4, DN_re=5.0, Omega_bh2=0.0223828,
         Omega_ch2=0.1201075, H0=67.32117, A_s=2.100549e-9),
    dict(label="stiff", r=5.0e-3, n_t=0.0, kappa10=50.0,
         T_re=1.0e2, DN_re=30.0, Omega_bh2=0.0223828,
         Omega_ch2=0.1201075, H0=67.32117, A_s=2.100549e-9),
]


def _env_meta():
    import subprocess
    out = dict(python=sys.version.split()[0], platform=sys.platform,
               numpy=np.__version__)
    for mod in ("torch", "scipy", "numba", "astropy", "sklearn"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            pass
    try:
        out["git_commit"] = subprocess.run(["git", "rev-parse", "HEAD"],
                                           cwd=REPO, capture_output=True,
                                           text=True, timeout=10).stdout.strip()
    except Exception:
        pass
    return out


def _fast_solve(point, profile):
    """One stiffgwpy fast solve; returns dict with spectrum/time or rejection."""
    m = LCDM_SG(r=point["r"], n_t=point["n_t"], kappa10=point["kappa10"],
                T_re=point["T_re"], DN_re=point["DN_re"],
                Omega_bh2=point["Omega_bh2"], Omega_ch2=point["Omega_ch2"],
                H0=point["H0"], A_s=point["A_s"])
    # Use the SAME grid-independent frequency set and z_tail as the oracle so
    # the fast-vs-reference residual isolates engine error (not grid mismatch).
    FS.apply_accuracy_mode("production" if profile == "production" else "fast")
    FS.set_z_tail(8.0)
    t0 = time.perf_counter()
    res = FS.SGWB_iter_fast(
        m, tol=1e-7 if profile == "production" else 1e-6,
        freq_res=1.0,
        transition_refine=(profile == "production"),
        freq_grid="grid_independent")
    dt = time.perf_counter() - t0
    if res is None or not getattr(m, "SGWB_converge", False):
        return {"status": "rejected",
                "reason": getattr(m, "fast_failure_reason", "unknown"),
                "runtime_s": dt}
    return {"status": "ok", "runtime_s": dt, "f": np.asarray(m.f, float),
            "log10OmegaGW": np.asarray(m.log10OmegaGW, float),
            "DN_gw": float(m.DN_gw[-1]), "DN_eff": float(m.cosmo_param["DN_eff"])}


def _reference_solve(point, z_tail=8.0, rtol=1e-6):
    """Independent continuous-sigma oracle on the signal-region frequency subset.

    Full-spectrum reference is prohibitive (~6+ min/point) because the
    sub-horizon-today low-frequency modes integrate the ODE all the way to
    today.  We therefore restrict to ``log10 f >= -8`` (re-entered, signal /
    transition region) so the oracle is the truth anchor for the physics that
    drives ``Omega_GW`` and the likelihood bins.
    """
    m = LCDM_SG(r=point["r"], n_t=point["n_t"], kappa10=point["kappa10"],
                T_re=point["T_re"], DN_re=point["DN_re"],
                Omega_bh2=point["Omega_bh2"], Omega_ch2=point["Omega_ch2"],
                H0=point["H0"], A_s=point["A_s"])
    try:
        G = grid_independent_freqs(m, 1.0)[0]
        G = np.asarray(G, dtype=float)
        # Core signal / transition band: re-entered modes, fast (tail-triggered),
        # and the region that drives the observables.  The sub-horizon-today
        # low-f tail is intentionally excluded (ill-defined static Omega_GW).
        subset = G[(G >= -6.0) & (G <= 1.0)]
        # Common background: use the fast production self-consistent DN_eff so
        # the oracle and fast solve the SAME background (isolates engine error).
        _prod = _fast_solve(point, "production")
        dz = (_prod["DN_eff"] if _prod.get("status") == "ok"
              else float(m.cosmo_param["DN_eff"]))
        t0 = time.perf_counter()
        from stiffgwpy.global_param import Neff0, Omega_nh2
        Omega_nu = Omega_nh2 / m.derived_param["h"] ** 2
        Ogw, Oj, Opgw, used = REF.spectrum_reference(m, subset, dz, z_tail=z_tail,
                                                     rtol=rtol, workers=8)
        g2, qerr, ierr = REF.integrate_spectrum(subset, Ogw, Oj)
        dn_gw = float(Neff0 * g2 / Omega_nu)
        lo = np.log10(np.maximum(Ogw - Oj, 1e-300))
        dt = time.perf_counter() - t0
    except Exception as exc:
        return {"status": "failed", "error": repr(exc)}
    return {"status": "ok", "runtime_s": dt,
            "f": np.asarray(subset, float),
            "log10OmegaGW": np.asarray(lo, float),
            "DN_gw_signal": float(dn_gw),
            "n_freq": int(subset.size),
            "quadrature_error": float(qerr), "interpolation_error": float(ierr),
            "used_tail_fraction": float(np.mean(used))}


def _sage_predict(point, model_type):
    from sagenetgw.classes import GWPredictor
    pred = GWPredictor(model_type=model_type, device="cpu")
    t0 = time.perf_counter()
    out = pred.predict({
        "r": point["r"], "n_t": point["n_t"], "kappa10": point["kappa10"],
        "T_re": point["T_re"], "DN_re": point["DN_re"],
        "Omega_bh2": point["Omega_bh2"], "Omega_ch2": point["Omega_ch2"],
        "H0": point["H0"], "A_s": point["A_s"]})
    dt = time.perf_counter() - t0
    return {"status": "ok", "runtime_s": dt,
            "f": np.asarray(out["f"], float),
            "log10OmegaGW": np.asarray(out["log10OmegaGW"], float)}


def _metrics_vs_ref(ref, test, label):
    """SageNet-author metrics + stiffgwpy signal-region metrics on a common grid."""
    f_ref = np.asarray(ref["f"]); l_ref = np.asarray(ref["log10OmegaGW"])
    f_t = np.asarray(test["f"]); l_t = np.asarray(test["log10OmegaGW"])
    if not (np.isfinite(f_t).all() and np.isfinite(l_t).all() and np.isfinite(f_ref).all()):
        return {"status": "invalid"}
    from scipy.interpolate import PchipInterpolator
    # common log10 f grid covering the intersection
    lo = max(f_ref.min(), f_t.min()); hi = min(f_ref.max(), f_t.max())
    if hi <= lo:
        return {"status": "disjoint"}
    g = np.unique(np.concatenate([f_ref, f_t]))
    g = g[(g >= lo) & (g <= hi)]
    r_ref = PchipInterpolator(np.sort(f_ref), l_ref[np.argsort(f_ref)])
    r_t = PchipInterpolator(np.sort(f_t), l_t[np.argsort(f_t)])
    lo_r = r_ref(g); lo_t = r_t(g)
    # signal mask: stiffgwpy signal-band convention on the reference region
    sig = (g >= -6.0) & (g <= 1.0)
    # local dex error
    dex = np.abs(lo_r - lo_t)
    # Omega-space relative error (linear), guard tiny denominators
    Om_ref = np.power(10.0, lo_r); Om_t = np.power(10.0, lo_t)
    rel = np.abs(Om_t - Om_ref) / np.maximum(Om_ref, 1e-300)
    # SageNet author metrics (need numpy arrays ordered as [f, log10OmegaGW])
    import sagenetgw.metrics as SM
    from numpy import column_stack
    true_coords = column_stack((f_ref, l_ref))
    pred_coords = column_stack((f_t, l_t))
    try:
        area = SM.calculate_area_difference(true_coords, pred_coords)
    except Exception:
        area = float("nan")
    try:
        smape = SM.calculate_smape(true_coords, pred_coords)
    except Exception:
        smape = float("nan")
    return {
        "status": "ok", "n_common": int(g.size),
        "dex_max": float(np.max(dex)),
        "dex_p95": float(np.percentile(dex, 95)),
        "dex_mean": float(np.mean(dex)),
        "rel_max": float(np.max(rel)),
        "rel_p95": float(np.percentile(rel, 95)),
        "rel_mean": float(np.mean(rel)),
        "signal_dex_max": float(np.max(dex[sig])) if sig.any() else float("nan"),
        "signal_rel_max": float(np.max(rel[sig])) if sig.any() else float("nan"),
        "smape": float(smape),
        "area_diff": float(area),
    }


def phase_quick(outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    for pt in POINTS:
        rec = {"label": pt["label"],
               "params": {k: float(v) for k, v in pt.items() if k != "label"}}
        for profile in ("production", "plain_grid"):
            r = _fast_solve(pt, profile)
            rec[profile] = r
        # SageNet Transformer
        try:
            rec["sagenet_transformer"] = _sage_predict(pt, "Transformer")
        except Exception as exc:
            rec["sagenet_transformer"] = {"status": "failed", "error": repr(exc)}
        # internal fast-vs-fast comparison (no reference yet)
        if rec["plain_grid"].get("status") == "ok" and rec["production"].get("status") == "ok":
            rec["plain_vs_prod_dn_gw_rel"] = abs(
                rec["plain_grid"]["DN_gw"] - rec["production"]["DN_gw"]) / abs(
                rec["production"]["DN_gw"])
        rows.append(rec)
        print("done", pt["label"])
    path = os.path.join(outdir, "quick.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"generated": _env_meta(), "points": rows}, fh, indent=2,
                  default=_json_default)
    print("wrote", path)


def phase_reference(outdir, labels):
    """Run the truth anchor (reference) for the chosen points and add to quick.json."""
    path = os.path.join(outdir, "quick.json")
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    by_label = {p["label"]: p for p in data["points"]}
    for label in labels:
        pt = next((p for p in POINTS if p["label"] == label), None)
        if pt is None:
            print("unknown label", label); continue
        ref = _reference_solve(pt)
        if label in by_label:
            by_label[label]["reference"] = ref
            print("ref done", label, ref.get("runtime_s"))
    data["points"] = list(by_label.values())
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, default=_json_default)


def phase_report(outdir):
    path = os.path.join(outdir, "quick.json")
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    rows = []
    for p in data["points"]:
        ref = p.get("reference")
        row = {"label": p["label"], "params": p["params"]}
        if ref and ref.get("status") == "ok":
            for profile in ("production", "plain_grid"):
                if p.get(profile, {}).get("status") == "ok":
                    row[profile] = _metrics_vs_ref(ref, p[profile], profile)
            if p.get("sagenet_transformer", {}).get("status") == "ok":
                row["sagenet_transformer"] = _metrics_vs_ref(
                    ref, p["sagenet_transformer"], "sagenet")
        row["plain_vs_prod_dn_gw_rel"] = p.get("plain_vs_prod_dn_gw_rel")
        rows.append(row)
    out = {"generated": _env_meta(), "rows": rows}
    with open(os.path.join(outdir, "report.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, default=_json_default)
    print("wrote report.json; points:", len(rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="quick",
                    choices=["quick", "reference", "report"])
    ap.add_argument("--points", nargs="*", default=[])
    ap.add_argument("--outdir", default=os.path.join(REPO, "sagenet_compare", "data"))
    args = ap.parse_args()
    if args.phase == "quick":
        phase_quick(args.outdir)
    elif args.phase == "reference":
        phase_reference(args.outdir, args.points)
    elif args.phase == "report":
        phase_report(args.outdir)


if __name__ == "__main__":
    main()
