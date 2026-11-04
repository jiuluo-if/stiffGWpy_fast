"""Full-grid Prüfer tail-convergence check for QMC sample sobol[146]."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from scripts.pressure_fast_domain import sobol_cases  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20261105
POWER = 8
SAMPLE_INDEX = 146
Z_START = -12.0
TAILS = (5.0, 10.0)
PRUFER_RTOL = 1e-9
CHUNK_NOTICE = 8
QMC_REFERENCE = ROOT / "docs" / "qmc_sobol146_reference_crosscheck_20261105.json"
OUTPUT = ROOT / "docs" / "qmc_sobol146_full_prufer_tail_convergence_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _spectrum_log(ogw: float, oj: float) -> float | None:
    omega = float(ogw - oj)
    return math.log10(omega) if omega > 0.0 and math.isfinite(omega) else None


def _map_summary(freqs, rows, key):
    if not all(row["solutions"].get(key, {}).get("status") == "success" for row in rows):
        return {"complete": False, "DN_gw": None, "max_abs_log10": None}
    ogw = np.array([row["solutions"][key]["Ogw_today"] for row in rows])
    oj = np.array([row["solutions"][key]["Oj_today"] for row in rows])
    omega = ogw - oj
    if not np.all(np.isfinite(omega) & (omega > 0.0)):
        return {"complete": False, "DN_gw": None, "max_abs_log10": None,
                "reason": "nonpositive_or_nonfinite_spectrum"}
    g2, qerr, ierr = R.integrate_spectrum(freqs, ogw, oj)
    model = rows[0]["_model"]
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    return {
        "complete": True,
        "DN_gw": float(gp.Neff0 * g2 / omega_nu),
        "quadrature_error": qerr,
        "interpolation_error": ierr,
        "max_abs_log10": float(np.max(np.abs(np.log10(omega)))),
    }


def main() -> None:
    qmc_row = json.loads(QMC_REFERENCE.read_text(encoding="utf-8"))
    freqs = np.asarray([
        row["frequency_log10_hz"] for row in qmc_row["frequency_records"]
    ], dtype=float)
    item = sobol_cases(POWER, SEED)[SAMPLE_INDEX]
    params = item["params"]
    if item["label"] != qmc_row["label"]:
        raise AssertionError("QMC sample label does not match the archived validation")

    # Reproduce the Fast point and evaluate a candidate with exact background
    # and a deep tail at the same native frequencies.
    FS.apply_accuracy_mode("fast")
    fast_models = {}
    fast_info = {}
    for name, tail, sigma_exact in (
            ("fast_tail5", 5.0, False),
            ("sigma_exact_tail10", 10.0, True)):
        model = LCDM_SG(**params)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            started = time.perf_counter()
            out = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-6, z_tail=tail,
                freq_res=1.0, threads=2, accuracy_mode="fast",
                physical_guard=False, transition_refine=False,
                sigma_exact=sigma_exact,
            )
            elapsed = time.perf_counter() - started
        if out is not model:
            fast_info[name] = {"status": "failed",
                               "reason": model.fast_failure_reason,
                               "runtime_s": elapsed}
            continue
        if not np.array_equal(np.asarray(model.f, dtype=float), freqs):
            raise AssertionError(f"Fast mode {name} changed the native frequency array")
        fast_models[name] = model
        fast_info[name] = {
            "status": "success",
            "z_tail": tail,
            "sigma_exact": sigma_exact,
            "runtime_s": elapsed,
            "DN_eff": float(model.cosmo_param["DN_eff"]),
            "DN_gw": float(model.DN_gw[-1]),
            "local_error_budget": model.local_error_budget,
            "warnings": [str(w.message) for w in caught],
            "log10OmegaGW": np.asarray(model.log10OmegaGW, dtype=float).tolist(),
            "omegaGW_linear": (np.asarray(model.Ogw_today, dtype=float)
                               - np.asarray(model.Oj_today, dtype=float)).tolist(),
            "handoff_eps": np.asarray(model.handoff_eps, dtype=float).tolist(),
        }

    model_ref = LCDM_SG(**params)
    dn_eff = float(fast_models["fast_tail5"].cosmo_param["DN_eff"])
    bg = R._Background(model_ref, dn_eff)
    bg_warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    bg_warm_s = time.perf_counter() - bg_warm_started

    rows = []
    mode_started = time.perf_counter()
    for i, freq in enumerate(freqs):
        row = {"frequency_log10_hz": float(freq), "solutions": {}}
        for tail in TAILS:
            key = f"prufer_tail{int(tail)}"
            started = time.perf_counter()
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    solution = R.solve_reference_mode(
                        model_ref, float(freq), dn_eff, z_tail=tail,
                        z_start=Z_START, rtol=PRUFER_RTOL, backend="prufer",
                    )
                omega = float(solution["Ogw_today"] - solution["Oj_today"])
                row["solutions"][key] = {
                    "status": "success",
                    "elapsed_s": time.perf_counter() - started,
                    "Ogw_today": float(solution["Ogw_today"]),
                    "Oj_today": float(solution["Oj_today"]),
                    "Opgw_today": float(solution["Opgw_today"]),
                    "Omega_GW_today": omega,
                    "log10Omega_GW_today": _spectrum_log(
                        solution["Ogw_today"], solution["Oj_today"]),
                    "used_tail": bool(solution["used_tail"]),
                    "event_N": solution["event_N"],
                    "nfev": int(solution["n_steps"]),
                    "matching_error_rel": solution["matching_error_rel"],
                    "warnings": [str(w.message) for w in caught],
                }
            except Exception as exc:
                row["solutions"][key] = {
                    "status": "exception",
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                    "elapsed_s": time.perf_counter() - started,
                }
        # Runtime metadata needed only for the aggregate integration helper.
        row["_model"] = model_ref
        rows.append(row)
        if (i + 1) % CHUNK_NOTICE == 0 or i + 1 == freqs.size:
            print(f"Pruefer tail maps: {i + 1}/{freqs.size} modes", flush=True)

    total_mode_s = time.perf_counter() - mode_started
    maps = {}
    for tail in TAILS:
        key = f"prufer_tail{int(tail)}"
        maps[key] = _map_summary(freqs, rows, key)
    cartesian_tail5_log = np.array([
        record["reference_cartesian_rtol1e-11_log10_omega"]
        for record in qmc_row["frequency_records"]
    ])
    if maps["prufer_tail5"].get("complete"):
        prufer_tail5_log = np.array([
            row["solutions"]["prufer_tail5"]["log10Omega_GW_today"]
            for row in rows
        ])
        delta_cart_prufer = np.abs(prufer_tail5_log - cartesian_tail5_log)
        maps["prufer_tail5"]["crosscheck_vs_cartesian_tail5"] = {
            "max_abs_spectrum_delta_dex": float(delta_cart_prufer.max()),
            "p95_abs_spectrum_delta_dex": float(
                np.percentile(delta_cart_prufer, 95)),
            "cartesian_DN_gw": qmc_row["reference"][
                "cartesian_rtol_1e-11_DN_gw"],
            "relative_DN_gw_delta": abs(
                maps["prufer_tail5"]["DN_gw"]
                - qmc_row["reference"]["cartesian_rtol_1e-11_DN_gw"]
            ) / max(
                abs(maps["prufer_tail5"]["DN_gw"]),
                abs(qmc_row["reference"]["cartesian_rtol_1e-11_DN_gw"]),
                np.finfo(float).tiny,
            ),
        }

    fast_comparisons = {}
    for fast_name, info in fast_info.items():
        if info.get("status") != "success":
            continue
        fast_log = np.asarray(info["log10OmegaGW"], dtype=float)
        fast_comparisons[fast_name] = {}
        for tail in TAILS:
            key = f"prufer_tail{int(tail)}"
            if not maps[key].get("complete"):
                fast_comparisons[fast_name][key] = {"complete": False}
                continue
            ref_log = np.array([
                row["solutions"][key]["log10Omega_GW_today"] for row in rows
            ])
            delta = np.abs(fast_log - ref_log)
            fast_dn = info["DN_gw"]
            ref_dn = maps[key]["DN_gw"]
            fast_comparisons[fast_name][key] = {
                "complete": True,
                "max_abs_spectrum_error_dex": float(np.max(delta)),
                "p95_abs_spectrum_error_dex": float(np.percentile(delta, 95)),
                "median_abs_spectrum_error_dex": float(np.median(delta)),
                "max_error_frequency_log10_hz": float(freqs[int(np.argmax(delta))]),
                "Fast_DN_gw": fast_dn,
                "Reference_DN_gw": ref_dn,
                "DN_gw_abs_difference": abs(fast_dn - ref_dn),
                "DN_gw_rel_difference": abs(fast_dn - ref_dn)
                / max(abs(fast_dn), abs(ref_dn), np.finfo(float).tiny),
            }

    # Remove non-serializable model references from per-mode rows after integration.
    for row in rows:
        row.pop("_model", None)
    result = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_full_prufer_tail_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_source_file": "docs/fast_domain_pressure_sobol_20261105_256.jsonl",
        "qmc_crosscheck_file": "docs/qmc_sobol146_reference_crosscheck_20261105.json",
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameters": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(
            freqs.astype("<f8").tobytes()
        ).hexdigest(),
        "frequency_range": [float(np.min(freqs)), float(np.max(freqs))],
        "reference_background_DN_eff": dn_eff,
        "reference_backend": "prufer",
        "reference_rtol": PRUFER_RTOL,
        "reference_z_start": Z_START,
        "reference_tail_thresholds": TAILS,
        "reference_background_warm_s": bg_warm_s,
        "reference_total_mode_s": total_mode_s,
        "maps": maps,
        "Fast": fast_info,
        "Fast_vs_independent_reference": fast_comparisons,
        "fast_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "reference_workers": 1,
        },
        "per_mode": rows,
        "limitations": [
            "This uses the QMC point's 77 Fast-native frequencies; frequency-grid convergence remains separate.",
            "Pruefer is an independent propagation formulation, not exact truth; agreement with Cartesian tail5 and per-mode crosschecks are separately archived.",
            "The Fast point remains one uncertified sample, not a domain boundary or general parameter claim.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"maps": maps, "Fast_vs_independent_reference": fast_comparisons,
                      "output": str(OUTPUT)}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
