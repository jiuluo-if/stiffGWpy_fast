"""Standalone coarse fixed-point predictor feasibility experiment.

The production outer loop often evaluates a first map at the original
``DN_eff`` and then a second full propagation at the updated value.  This
prototype replaces that first full map with a sparse Cartesian pre-solve and
uses a fixed contraction extrapolation before one full evaluation.  It is
diagnostic only: no production outer-loop semantics are changed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import CASES  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.exact_background import fast_phi_s2_split  # noqa: E402
from stiffgwpy_fast.freq_adaptive import goal_oriented_freqs  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

# With the production d ln(f) measure restored below, the first fixed-point
# iterate is the only fixed gain justified without an additional contraction
# estimate: x_1 = x_0 + g(x_0).
FIXED_POINT_GAIN = 1.0


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _frequency_subset(freqs, count):
    if count is None or count >= len(freqs):
        return np.asarray(freqs, dtype=np.float64)
    indices = np.unique(np.linspace(0, len(freqs) - 1, count, dtype=int))
    return np.asarray(freqs[indices], dtype=np.float64)


def _evaluate_once(model, dn_eff, coarse_count=None):
    """Evaluate one Cartesian fast solve at a specified total ``DN_eff``."""
    if not np.isfinite(dn_eff) or dn_eff > 5.0:
        return {"status": "physical_guard", "dn_eff": float(dn_eff)}
    model.cosmo_param["DN_eff"] = float(dn_eff)
    settings = FS.get_settings()
    h = float(settings["h"])
    z_tail = float(settings["z_tail"])
    col_step = int(settings["col_step"])
    phase_max = float(settings["phase_max"])
    FS.gen_fast(model, h, kink_split=True)
    kink_index, kink_fraction = FS._correct_kink_background(model)
    freqs = goal_oriented_freqs(model, 1.0, seed_n=64, max_points=120, eval_freqs=None)
    freqs = _frequency_subset(np.asarray(freqs, dtype=np.float64), coarse_count)
    Nv = np.asarray(model.Nv, dtype=np.float64)
    Sv, f_hor, j0s, z0s, fp_minus = FS.prep_frequency_only(model, Nv, freqs)
    Phi_grid, Phi_mid, S2, S2inv, kink_index, kink_fraction, phi_re = fast_phi_s2_split(
        model, Nv, float(dn_eff), sigma_nodes=model.sigma)
    idx_out = np.unique(np.append(np.arange(0, len(Nv), col_step), len(Nv) - 1))
    n_coarse = len(idx_out)
    P_t = model.derived_param["A_t"] * np.power(
        (10.0**freqs) / FS.gp.f_piv, model.derived_param["nt"])
    ev_minus = np.exp(-Nv)
    fp_freq = np.power(10.0, freqs)
    Ogw = np.zeros((len(freqs), n_coarse))
    Oj = np.zeros((len(freqs), n_coarse))
    Opgw = np.zeros((len(freqs), n_coarse))
    handoff_eps = np.full(len(freqs), -1.0)
    solve_args = (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
        fp_minus, fp_freq, 1, n_coarse, col_step, h, z_tail, Ogw, Oj,
        Opgw, None, model.sigma, phase_max, handoff_eps,
        kink_index, kink_fraction, phi_re,
    )
    FS.solve_kernel(*solve_args)
    integrand = Ogw[:, -1] - Oj[:, -1]
    # Match the production PCHIP path exactly: the native coordinate is
    # log10(f), while DN_gw integrates over d ln(f).
    g2 = float(np.sum(FS._pchip_integrals_vectorized(freqs, integrand)) * FS.ln10)
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    dn_gw = float(gp.Neff0 * g2 / omega_nu)
    if not np.isfinite(dn_gw) or not np.isfinite(Ogw).all() or not np.isfinite(Oj).all():
        return {"status": "numerical_failure", "dn_eff": float(dn_eff)}
    return {
        "status": "ok",
        "dn_eff": float(dn_eff),
        "dn_gw": dn_gw,
        "freqs": freqs,
        "Ogw": Ogw[:, -1].copy(),
        "Oj": Oj[:, -1].copy(),
        "Opgw": Opgw[:, -1].copy(),
        "log10OmegaGW": np.log10(Ogw[:, -1] - Oj[:, -1]),
        "handoff_eps": handoff_eps.copy(),
    }


def _run_candidate(name, coarse_count, threads):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    base_model = LCDM_SG(**CASES[name])
    original_dn_eff = float(base_model.cosmo_param["DN_eff"])
    coarse = _evaluate_once(base_model, original_dn_eff, coarse_count=coarse_count)
    if coarse["status"] != "ok":
        return {"status": coarse["status"], "coarse": coarse}
    predicted_dn_eff = original_dn_eff + FIXED_POINT_GAIN * coarse["dn_gw"]
    candidate_model = LCDM_SG(**CASES[name])
    candidate = _evaluate_once(candidate_model, predicted_dn_eff, coarse_count=None)
    return {
        "status": candidate["status"],
        "coarse_dn_gw": coarse["dn_gw"],
        "predicted_dn_eff": predicted_dn_eff,
        "coarse_n_freq": int(len(coarse["freqs"])),
        "candidate": candidate,
    }


def _compare_to_baseline(name, coarse_count, threads):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    baseline_model = LCDM_SG(**CASES[name])
    FS.SGWB_iter_fast(
        baseline_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip"
    )
    start = time.perf_counter()
    candidate_result = _run_candidate(name, coarse_count, threads)
    candidate_elapsed = time.perf_counter() - start
    if candidate_result["status"] != "ok":
        return {
            "case": name,
            "status_equal": False,
            "deterministic_equal": False,
            "candidate_status": candidate_result["status"],
            "candidate_elapsed_s": candidate_elapsed,
        }
    candidate = candidate_result["candidate"]
    base_freqs = np.asarray(baseline_model.f, dtype=np.float64)
    cand_freqs = np.asarray(candidate["freqs"], dtype=np.float64)
    order_base = np.argsort(base_freqs)
    order_cand = np.argsort(cand_freqs)
    baseline_interp = np.interp(
        cand_freqs[order_cand], base_freqs[order_base],
        np.asarray(baseline_model.log10OmegaGW)[order_base]
    )
    candidate_sorted = np.asarray(candidate["log10OmegaGW"])[order_cand]
    spectrum_max_dex = float(np.max(np.abs(candidate_sorted - baseline_interp)))
    baseline_dn = float(np.asarray(baseline_model.DN_gw)[-1])
    candidate_dn = float(candidate["dn_gw"])
    repeat_result = _run_candidate(name, coarse_count, threads)
    deterministic_equal = (
        repeat_result["status"] == "ok"
        and np.array_equal(candidate["freqs"], repeat_result["candidate"]["freqs"])
        and np.array_equal(candidate["log10OmegaGW"], repeat_result["candidate"]["log10OmegaGW"])
        and candidate_dn == float(repeat_result["candidate"]["dn_gw"])
    )
    return {
        "case": name,
        "status_equal": candidate["status"] == "ok" and getattr(baseline_model, "SGWB_converge", False) is True,
        "deterministic_equal": deterministic_equal,
        "baseline_converged": bool(getattr(baseline_model, "SGWB_converge", False)),
        "candidate_status": candidate["status"],
        "candidate_elapsed_s": candidate_elapsed,
        "baseline_dn_gw": baseline_dn,
        "candidate_dn_gw": candidate_dn,
        "candidate_dn_relative": abs(candidate_dn - baseline_dn) / max(abs(baseline_dn), 1e-300),
        "candidate_spectrum_max_dex": spectrum_max_dex,
        "frequency_count_equal": len(base_freqs) == len(cand_freqs),
        "frequency_max_abs": float(np.max(np.abs(base_freqs - cand_freqs))) if len(base_freqs) == len(cand_freqs) else float("inf"),
        "baseline_spectrum_digest": _digest(baseline_model.log10OmegaGW),
        "candidate_spectrum_digest": _digest(candidate["log10OmegaGW"]),
        "predicted_dn_eff": candidate_result["predicted_dn_eff"],
        "coarse_dn_gw": candidate_result["coarse_dn_gw"],
        "coarse_n_freq": candidate_result["coarse_n_freq"],
    }


def _run_case(name, coarse_count=32, repeats=1, threads=2):
    rows = [_compare_to_baseline(name, coarse_count, threads) for _ in range(repeats)]
    first = rows[0]
    return {
        **first,
        "repeats": repeats,
        "coarse_count": coarse_count,
        "threads": threads,
        "median_candidate_s": statistics.median(row["candidate_elapsed_s"] for row in rows),
        "p95_candidate_s": float(np.percentile([row["candidate_elapsed_s"] for row in rows], 95)),
        "all_status_equal": all(row["status_equal"] for row in rows),
        "all_deterministic_equal": all(row["deterministic_equal"] for row in rows),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--coarse-count", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--output")
    args = parser.parse_args()
    cases = args.case or list(CASES)
    payload = {
        "candidate": "coarse_fixed_point_predictor",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "fixed_point_gain": FIXED_POINT_GAIN,
        "cases": [_run_case(case, args.coarse_count, args.repeats, args.threads) for case in cases],
    }
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
