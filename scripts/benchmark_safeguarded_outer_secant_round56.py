"""Round 56 safeguarded scalar outer-secant feasibility screen.

The production outer loop solves ``F(d)=d-d0`` by bracketed updates, where
``F(d)`` is the exact fast propagation result.  This diagnostic uses two exact
full-map evaluations and a bracketed secant prediction, then evaluates the
full spectrum at that predicted point.  No sparse/coarse spectrum is used and
production code is never modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_coarse_fixed_point_predictor import (  # noqa: E402
    CASES,
    _evaluate_once,
)
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _secant_predict(d0, r0, d1, r1, lower, upper):
    denominator = r1 - r0
    if not math.isfinite(denominator) or abs(denominator) <= 1e-14:
        return 0.5 * (lower + upper)
    value = d1 - r1 * (d1 - d0) / denominator
    if not math.isfinite(value):
        return 0.5 * (lower + upper)
    return min(max(value, lower), upper)


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=np.float64).tobytes()).hexdigest()


def _residual(dn_eff, dn_orig, dn_gw):
    return float(dn_gw - (dn_eff - dn_orig))


def _run_case(case_name, threads=2):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(threads)
    baseline_model = LCDM_SG(**CASES[case_name])
    dn_orig = float(baseline_model.cosmo_param["DN_eff"])
    baseline_result = FS.SGWB_iter_fast(
        baseline_model, kink_split=True, freq_grid="goal", frequency_quadrature="pchip")
    if baseline_result is None:
        return {
            "case": case_name,
            "status": "baseline_failure",
            "failure_reason": getattr(baseline_model, "fast_failure_reason", None),
        }

    first_model = LCDM_SG(**CASES[case_name])
    first = _evaluate_once(first_model, dn_orig, coarse_count=None)
    if first["status"] != "ok":
        return {"case": case_name, "status": "first_map_failure", "first": first}
    d0 = dn_orig
    d1 = d0 + float(first["dn_gw"])
    r0 = _residual(d0, dn_orig, first["dn_gw"])
    second_model = LCDM_SG(**CASES[case_name])
    second = _evaluate_once(second_model, d1, coarse_count=None)
    if second["status"] != "ok":
        return {"case": case_name, "status": "second_map_failure", "second": second}
    r1 = _residual(d1, dn_orig, second["dn_gw"])
    predicted = _secant_predict(d0, r0, d1, r1, dn_orig, 5.0)
    third_model = LCDM_SG(**CASES[case_name])
    third = _evaluate_once(third_model, predicted, coarse_count=None)
    if third["status"] != "ok":
        return {"case": case_name, "status": "predicted_map_failure", "predicted_dn_eff": predicted}

    metric = abs((gp.Neff0 + dn_orig + third["dn_gw"])
                 / (gp.Neff0 + predicted) - 1.0)
    base_f = np.asarray(baseline_model.f, dtype=np.float64)
    cand_f = np.asarray(third["freqs"], dtype=np.float64)
    base_log = np.asarray(baseline_model.log10OmegaGW, dtype=np.float64)
    cand_log = np.asarray(third["log10OmegaGW"], dtype=np.float64)
    if np.array_equal(base_f, cand_f):
        spectrum_dex = float(np.max(np.abs(base_log - cand_log)))
    else:
        order = np.argsort(base_f)
        candidate_order = np.argsort(cand_f)
        interp = np.interp(cand_f[candidate_order], base_f[order], base_log[order])
        spectrum_dex = float(np.max(np.abs(interp - cand_log[candidate_order])))
    baseline_dn = float(np.asarray(baseline_model.DN_gw)[-1])
    candidate_final_dn = dn_orig + float(third["dn_gw"])
    return {
        "case": case_name,
        "status": "ok",
        "baseline_converged": bool(getattr(baseline_model, "SGWB_converge", False)),
        "baseline_outer_values": int(len(getattr(baseline_model, "DN_gw", []))),
        "exact_map_evaluations": 3,
        "dn_orig": dn_orig,
        "d0": d0,
        "d1": d1,
        "predicted_dn_eff": predicted,
        "r0": r0,
        "r1": r1,
        "predicted_convergence_metric": float(metric),
        "baseline_dn_gw": baseline_dn,
        "candidate_predicted_dn_gw": float(third["dn_gw"]),
        "candidate_final_dn_eff": candidate_final_dn,
        "candidate_dn_relative_to_baseline": abs(candidate_final_dn - (dn_orig + baseline_dn)) / max(abs(baseline_dn), 1e-300),
        "spectrum_max_dex_to_baseline": spectrum_dex,
        "frequency_equal": bool(np.array_equal(base_f, cand_f)),
        "candidate_spectrum_digest": _digest(cand_log),
        "baseline_spectrum_digest": _digest(base_log),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--case", action="append", choices=sorted(CASES))
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    payload = {
        "experiment": "round56_safeguarded_outer_secant_screen",
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "production_unchanged": True,
        "cases": [_run_case(case, args.threads) for case in (args.case or list(CASES))],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
