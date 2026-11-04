"""Targeted cross-backend/reference-tail checks for high-error QMC modes."""

from __future__ import annotations

import hashlib
import json
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
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20261105
POWER = 8
SAMPLE_INDEX = 146
Z_START = -12.0
TAILS = (5.0, 10.0)
FAST_VARIANTS = (
    ("sigma_exact_tail5", 5.0, False, True),
    ("sigma_exact_tail10", 10.0, False, True),
    ("transition_sigma_tail5", 5.0, True, True),
    ("transition_sigma_tail10", 10.0, True, True),
)
MODES = (
    -17.026818276280245,  # largest Fast/reference error; no tail handoff
    -16.50446404464904,   # second largest no-tail error
    -16.086580659344087,  # largest recorded Fast handoff epsilon
    -1.8823780790633404,  # representative tail-mode spectrum error
    -1.0073780790633404,  # opposite-sign tail-mode error
)
OUTPUT = ROOT / "docs" / "qmc_sobol146_mode_tail_probe_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> None:
    item = sobol_cases(POWER, SEED)[SAMPLE_INDEX]
    params = item["params"]
    FS.apply_accuracy_mode("fast")
    model = LCDM_SG(**params)
    with warnings.catch_warnings(record=True) as fast_warnings:
        warnings.simplefilter("always")
        fast_started = time.perf_counter()
        out = model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-6, z_tail=5.0,
            freq_res=1.0, threads=2, accuracy_mode="fast",
            physical_guard=False,
        )
        fast_s = time.perf_counter() - fast_started
    if out is not model:
        raise RuntimeError(f"Fast sample failed: {model.fast_failure_reason}")

    freqs = np.asarray(model.f, dtype=float)
    if not all(np.any(freqs == f) for f in MODES):
        raise AssertionError("target mode missing from regenerated Fast frequency grid")
    index = {float(f): int(np.flatnonzero(freqs == f)[0]) for f in MODES}
    fast_omega = np.asarray(model.Ogw_today) - np.asarray(model.Oj_today)
    fast_log = np.log10(fast_omega)
    fast_handoff = np.asarray(model.handoff_eps, dtype=float)
    dn_eff = float(model.cosmo_param["DN_eff"])

    # Warm only the background spline/H2 implementation; these direct single
    # mode checks deliberately exercise the established Python tensor RHS.
    bg = R._Background(model, dn_eff)
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    reference_warm_s = time.perf_counter() - warm_started

    runs = []
    for freq in MODES:
        mode_runs = {}
        for tail in TAILS:
            for backend, rtol in (("cartesian", 1e-11), ("prufer", 1e-9)):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    started = time.perf_counter()
                    solution = R.solve_reference_mode(
                        model, freq, dn_eff, z_tail=tail, z_start=Z_START,
                        rtol=rtol, backend=backend,
                    )
                    elapsed = time.perf_counter() - started
                omega = float(solution["Ogw_today"] - solution["Oj_today"])
                mode_runs[f"{backend}_tail{int(tail)}"] = {
                    "rtol": rtol,
                    "elapsed_s": elapsed,
                    "Omega_GW_today": omega,
                    "log10Omega_GW_today": math_log10(omega),
                    "used_tail": bool(solution["used_tail"]),
                    "event_N": solution["event_N"],
                    "nfev": solution["n_steps"],
                    "matching_error_rel": solution["matching_error_rel"],
                    "warnings": [str(w.message) for w in caught],
                }

        fi = index[float(freq)]
        fast_log_value = float(fast_log[fi])
        for run in mode_runs.values():
            run["fast_log10Omega_GW"] = fast_log_value
            run["fast_minus_reference_dex"] = (
                fast_log_value - run["log10Omega_GW_today"]
            )
        runs.append({
            "frequency_log10_hz": float(freq),
            "fast_handoff_eps": float(fast_handoff[fi]),
            "runs": mode_runs,
        })

    records = []
    for run in runs:
        r = run["runs"]
        records.append({
            "frequency_log10_hz": run["frequency_log10_hz"],
            "fast_handoff_eps": run["fast_handoff_eps"],
            "fast_reference_error_tail5_dex": r["cartesian_tail5"]["fast_minus_reference_dex"],
            "fast_reference_error_tail10_dex": r["cartesian_tail10"]["fast_minus_reference_dex"],
            "cartesian_prufer_delta_tail5_dex": (
                r["cartesian_tail5"]["log10Omega_GW_today"]
                - r["prufer_tail5"]["log10Omega_GW_today"]),
            "cartesian_tail5_tail10_delta_dex": (
                r["cartesian_tail5"]["log10Omega_GW_today"]
                - r["cartesian_tail10"]["log10Omega_GW_today"]),
            "cartesian_tail5_nfev": r["cartesian_tail5"]["nfev"],
            "cartesian_tail10_nfev": r["cartesian_tail10"]["nfev"],
            "cartesian_tail10_runtime_s": r["cartesian_tail10"]["elapsed_s"],
            "cartesian_tail10_matching_error": r["cartesian_tail10"]["matching_error_rel"],
            "used_tail_tail5": r["cartesian_tail5"]["used_tail"],
            "used_tail_tail10": r["cartesian_tail10"]["used_tail"],
        })

    variant_results = []
    for name, tail, transition_refine, sigma_exact in FAST_VARIANTS:
        variant_model = LCDM_SG(**params)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            started = time.perf_counter()
            variant_out = variant_model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-6, z_tail=tail,
                freq_res=1.0, threads=2, accuracy_mode="fast",
                physical_guard=False, transition_refine=transition_refine,
                sigma_exact=sigma_exact,
            )
            elapsed = time.perf_counter() - started
        if variant_out is not variant_model:
            raise RuntimeError(
                f"Fast variant {name} failed: {variant_model.fast_failure_reason}"
            )
        variant_freqs = np.asarray(variant_model.f, dtype=float)
        if not np.array_equal(variant_freqs, freqs):
            raise AssertionError(f"Fast variant {name} changed the native frequency grid")
        variant_omega = (np.asarray(variant_model.Ogw_today)
                         - np.asarray(variant_model.Oj_today))
        variant_log = np.log10(variant_omega)
        variant_rows = []
        for freq in MODES:
            fi = int(np.flatnonzero(variant_freqs == freq)[0])
            ref_modes = next(
                record["runs"] for record in runs
                if record["frequency_log10_hz"] == freq
            )
            variant_rows.append({
                "frequency_log10_hz": freq,
                "fast_log10Omega_GW": float(variant_log[fi]),
                "fast_minus_cartesian_tail5_dex": (
                    float(variant_log[fi])
                    - ref_modes["cartesian_tail5"]["log10Omega_GW_today"]),
                "fast_minus_cartesian_tail10_dex": (
                    float(variant_log[fi])
                    - ref_modes["cartesian_tail10"]["log10Omega_GW_today"]),
            })
        errors5 = np.abs([row["fast_minus_cartesian_tail5_dex"]
                          for row in variant_rows])
        errors10 = np.abs([row["fast_minus_cartesian_tail10_dex"]
                           for row in variant_rows])
        variant_results.append({
            "name": name,
            "z_tail": tail,
            "transition_refine": transition_refine,
            "sigma_exact": sigma_exact,
            "runtime_s": elapsed,
            "DN_gw": float(variant_model.DN_gw[-1]),
            "DN_eff_final": float(variant_model.cosmo_param["DN_eff"]),
            "warnings": [str(w.message) for w in caught],
            "max_abs_error_vs_cartesian_tail5_dex": float(errors5.max()),
            "max_abs_error_vs_cartesian_tail10_dex": float(errors10.max()),
            "modes": variant_rows,
        })

    result = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_mode_tail_probe",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "params": params,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(
            freqs.astype("<f8").tobytes()
        ).hexdigest(),
        "selected_modes": MODES,
        "z_start": Z_START,
        "tail_thresholds": TAILS,
        "fast": {
            "runtime_s": fast_s,
            "DN_eff_final": dn_eff,
            "DN_gw": float(model.DN_gw[-1]),
            "frequency_grid_error": float(model.freq_grid_error),
            "quadrature_error_local": float(model.quadrature_error_local),
            "local_error_budget": model.local_error_budget,
            "warnings": [str(w.message) for w in fast_warnings],
        },
        "reference_background_warm_s": reference_warm_s,
        "runs": runs,
        "fast_variants": variant_results,
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "reference_workers": 1,
        },
        "summary": records,
        "limitations": [
            "Only five selected modes are compared; this probe cannot certify the full map or root.",
            "The full-grid tail10 Cartesian/Pruefer comparison exceeded the 600 CPU-second screen budget; this targeted follow-up records per-mode tail costs.",
            "The Fast point remains uncertified; finite output is not a boundary claim.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


def math_log10(value: float) -> float | None:
    return float(np.log10(value)) if value > 0.0 and np.isfinite(value) else None


if __name__ == "__main__":
    main()
