"""Paired Fast screen for the H2-identity / nonadiabatic-WKB candidate.

The baseline runs from a temporary source archive of the pinned base commit;
the candidate runs from the current checkout. No Reference or production file
is modified by this screen.
"""

from __future__ import annotations

import hashlib
import importlib
import io
import json
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

exact_background = importlib.import_module("stiffgwpy_fast.exact_background")
fast_sgwb = importlib.import_module("stiffgwpy_fast.fast_sgwb")
LCDM_SG = importlib.import_module("stiffgwpy_fast.stiff_SGWB").LCDM_SG
H2_vec = exact_background.H2_vec
SGWB_iter_fast = fast_sgwb.SGWB_iter_fast
gen_fast = fast_sgwb.gen_fast
resolve_config = fast_sgwb.resolve_config


BASE_SHA = "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45"
BASELINE_TIMEOUT_SECONDS = 15.0
OUTPUT = Path(__file__).with_name("screen_h2_identity_region.json")
BASELINE_PARTIAL = OUTPUT.with_name("screen_h2_identity_region_baseline_partial.json")
BASE_PARAMS = {
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": 0.0,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
RECTANGLE = [
    {"Omega_bh2": ob, "T_re": tr}
    for ob in (1.0e5, 3.0e5, 1.0e6)
    for tr in (1000.0, 2000.0, 4000.0)
]
HOLDOUTS = [
    {"Omega_bh2": 5.0e5, "T_re": 2500.0},
    {"Omega_bh2": 7.5e5, "T_re": 1500.0},
]


BASELINE_RUNNER = r'''
import contextlib, json, math, os, sys, time, warnings
sys.path.insert(0, os.getcwd())
import numpy as np
from stiffgwpy_fast.exact_background import H2_vec
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, gen_fast, resolve_config
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

cases = json.load(sys.stdin)
rows = []
for case in cases:
    params = case["params"]
    preflight = LCDM_SG(**params)
    h = resolve_config("fast", threads=2).h
    gen_fast(preflight, h)
    try:
        h2 = H2_vec(preflight.Nv, preflight, params["DN_eff"])
        h2_ok = bool(np.isfinite(h2).all() and np.all(h2 > 0.0))
        h2_min = float(np.min(h2)) if np.isfinite(h2).all() else None
    except Exception as exc:
        h2_ok, h2_min = False, None
        h2_error = type(exc).__name__ + ": " + str(exc)
    else:
        h2_error = None

    model = LCDM_SG(**params)
    start = time.perf_counter()
    failure = None
    exception = None
    try:
        with contextlib.redirect_stdout(sys.stderr), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = SGWB_iter_fast(
                model, config=resolve_config("fast", threads=2), tol=1e-6,
                physical_guard=True,
            )
        warning_count = len(caught)
        success = result is model
        failure = model.fast_failure_reason
        root = float(model.DN_gw[-1]) if success else None
        finite_spectrum = bool(np.isfinite(model.log10OmegaGW).all()) if success else False
        n_freq = int(len(model.f)) if success else None
    except Exception as exc:
        warning_count = 0
        success, root, finite_spectrum, n_freq = False, None, False, None
        exception = type(exc).__name__ + ": " + str(exc)
    rows.append({
        "case_id": case["case_id"], "params": params,
        "derived_nt": float(preflight.derived_param["nt"]),
        "derived_N_inf": (
            None if preflight.derived_param["N_inf"] is None
            else float(preflight.derived_param["N_inf"])
        ),
        "h2_positive_finite": h2_ok, "h2_min": h2_min,
        "h2_error": h2_error, "success": success,
        "failure_reason": failure, "exception": exception,
        "root": root, "finite_spectrum": finite_spectrum,
        "n_freq": n_freq, "warning_count": warning_count,
        "h2_identity_used": bool(getattr(model, "_fast_h2_identity_used", False)),
        "nonadiabatic_handoff_selected": int(
            getattr(model, "_fast_nonadiabatic_handoff_selected", 0)
        ),
        "runtime_seconds": time.perf_counter() - start,
    })
print(json.dumps(rows))
'''


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cases() -> list[dict[str, object]]:
    rows = []
    for family, specs in (("rectangle", RECTANGLE), ("holdout", HOLDOUTS)):
        for index, override in enumerate(specs):
            params = dict(BASE_PARAMS)
            params.update(override)
            rows.append({
                "case_id": f"{family}_{index:02d}",
                "family": family,
                "params": params,
            })
    return rows


def _run_baseline(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    archive = subprocess.run(
        ["git", "archive", BASE_SHA], cwd=ROOT,
        check=True, stdout=subprocess.PIPE,
    ).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_fast_872a70_baseline_") as tmp:
        baseline_root = Path(tmp)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
            bundle.extractall(baseline_root, filter="data")
        rows = []
        for case in cases:
            start = time.perf_counter()
            try:
                completed = subprocess.run(
                    [sys.executable, "-c", BASELINE_RUNNER],
                    cwd=baseline_root,
                    input=json.dumps([case]),
                    text=True,
                    check=True,
                    capture_output=True,
                    timeout=BASELINE_TIMEOUT_SECONDS,
                )
            except subprocess.TimeoutExpired:
                preflight = LCDM_SG(**case["params"])
                gen_fast(preflight, resolve_config("fast", threads=2).h)
                try:
                    h2 = H2_vec(
                        preflight.Nv, preflight, case["params"]["DN_eff"]
                    )
                    h2_ok = bool(np.isfinite(h2).all() and np.all(h2 > 0.0))
                    h2_min = float(np.min(h2)) if np.isfinite(h2).all() else None
                    h2_error = None
                except Exception as exc:
                    h2_ok, h2_min = False, None
                    h2_error = type(exc).__name__ + ": " + str(exc)
                rows.append({
                    "case_id": case["case_id"], "params": case["params"],
                    "h2_positive_finite": h2_ok, "h2_min": h2_min,
                    "h2_error": h2_error, "success": False,
                    "failure_reason": "resource_timeout",
                    "exception": None, "runtime_seconds": time.perf_counter() - start,
                    "timeout_seconds": BASELINE_TIMEOUT_SECONDS,
                })
                BASELINE_PARTIAL.write_text(
                    json.dumps({"base_sha": BASE_SHA, "rows": rows}, indent=2) + "\n",
                    encoding="utf-8",
                )
                continue
            result_row = json.loads(completed.stdout)[0]
            result_row["process_runtime_seconds"] = time.perf_counter() - start
            rows.append(result_row)
            BASELINE_PARTIAL.write_text(
                json.dumps({"base_sha": BASE_SHA, "rows": rows}, indent=2) + "\n",
                encoding="utf-8",
            )
        return rows


def _preflight_h2(params: dict[str, object]) -> dict[str, object]:
    model = LCDM_SG(**params)
    gen_fast(model, resolve_config("fast", threads=2).h)
    try:
        h2 = H2_vec(model.Nv, model, params["DN_eff"])
        return {
            "h2_positive_finite": bool(np.isfinite(h2).all() and np.all(h2 > 0.0)),
            "h2_min": float(np.min(h2)) if np.isfinite(h2).all() else None,
            "h2_error": None,
        }
    except Exception as exc:
        return {
            "h2_positive_finite": False,
            "h2_min": None,
            "h2_error": type(exc).__name__ + ": " + str(exc),
        }


def _run_candidate(cases: list[dict[str, object]], baseline: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for index, (case, old) in enumerate(zip(cases, baseline)):
        start = time.perf_counter()
        try:
            completed = subprocess.run(
                [sys.executable, "-c", BASELINE_RUNNER],
                cwd=ROOT,
                input=json.dumps([case]),
                text=True,
                check=True,
                capture_output=True,
                timeout=BASELINE_TIMEOUT_SECONDS,
            )
            row = json.loads(completed.stdout)[0]
            row["process_runtime_seconds"] = time.perf_counter() - start
        except subprocess.TimeoutExpired:
            row = {
                "case_id": case["case_id"], "params": case["params"],
                "success": False, "failure_reason": "resource_timeout",
                "exception": None, "runtime_seconds": time.perf_counter() - start,
                "timeout_seconds": BASELINE_TIMEOUT_SECONDS,
            }
        row.update(_preflight_h2(case["params"]))
        rows.append(row)
        partial = {
            "status": "paired_screen_in_progress",
            "base_sha": BASE_SHA,
            "candidate_base_sha": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                capture_output=True, text=True,
            ).stdout.strip(),
            "candidate_source_sha256": {
                "exact_background.py": _sha256(ROOT / "stiffgwpy_fast" / "exact_background.py"),
                "fast_sgwb.py": _sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            },
            "paired_rows": [
                {"input": cases[i], "baseline_fast": baseline[i], "candidate_fast": rows[i]}
                for i in range(index + 1)
            ],
        }
        OUTPUT.with_suffix(".partial.json").write_text(
            json.dumps(partial, indent=2) + "\n", encoding="utf-8"
        )
    return rows


def main() -> None:
    cases = _cases()
    cases.sort(
        key=lambda case: (
            0
            if case["params"]["Omega_bh2"] == 1.0e5
            and case["params"]["T_re"] == 2000.0
            else 1,
            case["case_id"],
        )
    )
    if BASELINE_PARTIAL.exists():
        saved_baseline = json.loads(BASELINE_PARTIAL.read_text(encoding="utf-8"))
    else:
        saved_baseline = {}
    saved_rows = saved_baseline.get("rows", [])
    if (
        saved_baseline.get("base_sha") == BASE_SHA
        and [row.get("case_id") for row in saved_rows]
        == [case["case_id"] for case in cases]
    ):
        baseline = saved_rows
    else:
        baseline = _run_baseline(cases)
    candidate = _run_candidate(cases, baseline)
    result = {
        "status": "fast_only_domain_screen",
        "base_sha": BASE_SHA,
        "candidate_base_sha": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip(),
        "baseline_timeout_seconds": BASELINE_TIMEOUT_SECONDS,
        "params_contract": {
            "cr": float(BASE_PARAMS["cr"]),
            "input_n_t": float(BASE_PARAMS["n_t"]),
            "derived_tilt_rule": "model-derived consistency tilt for cr>0",
            "physical_guard": True,
            "frequency_grid": "goal",
            "threads": 2,
            "tol": 1e-6,
            "h": resolve_config("fast", threads=2).h,
            "z_tail": resolve_config("fast", threads=2).z_tail,
        },
        "candidate_source_sha256": {
            "exact_background.py": _sha256(ROOT / "stiffgwpy_fast" / "exact_background.py"),
            "fast_sgwb.py": _sha256(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        },
        "cases": [
            {"input": case, "baseline_fast": old, "candidate_fast": new}
            for case, old, new in zip(cases, baseline, candidate)
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    baseline_ok = sum(bool(row["success"]) for row in baseline)
    candidate_ok = sum(bool(row["success"]) for row in candidate)
    print(json.dumps({
        "baseline_successes": baseline_ok,
        "candidate_successes": candidate_ok,
        "case_count": len(cases),
        "artifact": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
