"""Cross-check one finite but uncertified QMC Fast point against Reference.

The fixed Fast-native frequency grid isolates per-mode propagation errors;
this does not certify the continuous-frequency integral or parameter domain.
"""

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

QMC_SEED = 20261105
QMC_POWER = 8
SAMPLE_INDEX = 146
REFERENCE_RTOLS = (1e-10, 1e-11)
PRUFER_RTOL = 1e-9
Z_TAIL = 5.0
Z_START = -12.0
OUTPUT = ROOT / "docs" / "qmc_sobol146_reference_crosscheck_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _spectrum_log(ogw, oj):
    omega = np.asarray(ogw, dtype=float) - np.asarray(oj, dtype=float)
    log_omega = np.full(omega.shape, np.nan)
    positive = np.isfinite(omega) & (omega > 0.0)
    log_omega[positive] = np.log10(omega[positive])
    return omega, log_omega, positive


def _reference_map(model, freqs, dn_eff, backend, rtol):
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = R.run_reference(
            model, dn_eff=dn_eff, freq_subset=freqs, z_tail=Z_TAIL,
            z_start=Z_START, rtol=rtol, self_consistent=False, workers=1,
            physical_guard=False, backend=backend,
        )
    return result, time.perf_counter() - started, [str(w.message) for w in caught]


def main() -> None:
    suite = sobol_cases(QMC_POWER, QMC_SEED)
    item = suite[SAMPLE_INDEX]
    params = item["params"]
    model = LCDM_SG(**params)
    FS.apply_accuracy_mode("fast")
    with warnings.catch_warnings(record=True) as fast_warnings:
        warnings.simplefilter("always")
        fast_started = time.perf_counter()
        out = model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-6, z_tail=Z_TAIL,
            freq_res=1.0, threads=2, accuracy_mode="fast",
            physical_guard=False,
        )
        fast_elapsed = time.perf_counter() - fast_started
    if out is not model:
        raise RuntimeError(f"Fast point no longer returns finite output: {model.fast_failure_reason}")

    freqs = np.asarray(model.f, dtype=float).copy()
    fast_omega, fast_log, fast_positive = _spectrum_log(
        model.Ogw_today, model.Oj_today)
    dn_eff = float(model.cosmo_param["DN_eff"])
    fast_dn = float(model.DN_gw[-1])
    background = R._Background(model, dn_eff)
    probe_N = np.linspace(background.N_re_abs, background.N_inf, 513)
    h2 = np.array([
        R._H2_and_sigma_python(background, float(n))[0] for n in probe_N
    ])

    cartesian = {}
    cartesian_times = {}
    cartesian_warnings = {}
    for rtol in REFERENCE_RTOLS:
        result, elapsed, caught = _reference_map(
            model, freqs, dn_eff, "cartesian", rtol)
        omega, log_omega, positive = _spectrum_log(result["Ogw"], result["Oj"])
        cartesian[rtol] = (result, omega, log_omega, positive)
        cartesian_times[rtol] = elapsed
        cartesian_warnings[rtol] = caught

    prufer_result, prufer_elapsed, prufer_warnings = _reference_map(
        model, freqs, dn_eff, "prufer", PRUFER_RTOL)
    prufer_omega, prufer_log, prufer_positive = _spectrum_log(
        prufer_result["Ogw"], prufer_result["Oj"])

    strict_result, strict_omega, strict_log, strict_positive = cartesian[1e-11]
    c10_omega, c10_log, c10_positive = cartesian[1e-10][1:]
    valid_fast_ref = fast_positive & strict_positive
    valid_c10_c11 = c10_positive & strict_positive
    valid_prufer_ref = prufer_positive & strict_positive
    fast_log_delta = np.abs(fast_log[valid_fast_ref] - strict_log[valid_fast_ref])
    ctol_log_delta = np.abs(c10_log[valid_c10_c11] - strict_log[valid_c10_c11])
    prufer_log_delta = np.abs(prufer_log[valid_prufer_ref] - strict_log[valid_prufer_ref])
    ref_dn = float(strict_result["DN_gw"])
    dn_abs = abs(fast_dn - ref_dn)
    reference_dn_rtol_delta = abs(
        float(cartesian[1e-10][0]["DN_gw"]) - ref_dn
    ) / max(abs(ref_dn), np.finfo(float).tiny)
    max_fast_i = int(np.argmax(np.abs(fast_log - strict_log)))
    frecords = []
    handoff_eps = np.asarray(getattr(model, "handoff_eps", np.full(freqs.shape, np.nan)))
    for i, f in enumerate(freqs):
        frecords.append({
            "frequency_log10_hz": float(f),
            "fast_log10_omega": float(fast_log[i]),
            "reference_cartesian_rtol1e-11_log10_omega": float(strict_log[i]),
            "reference_prufer_rtol1e-9_log10_omega": float(prufer_log[i]),
            "fast_minus_reference_dex": float(fast_log[i] - strict_log[i]),
            "cartesian_prufer_delta_dex": float(strict_log[i] - prufer_log[i]),
            "fast_handoff_eps": float(handoff_eps[i]) if i < handoff_eps.size else None,
            "reference_cartesian_used_tail": bool(strict_result["used_tail"][i]),
            "reference_prufer_used_tail": bool(prufer_result["used_tail"][i]),
        })

    result = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_fast_reference_crosscheck",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_source_file": "docs/fast_domain_pressure_sobol_20261105_256.jsonl",
        "qmc_seed": QMC_SEED,
        "qmc_power": QMC_POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameter_sha256": hashlib.sha256(
            json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "parameters": params,
        "fast_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "pressure_script_sha256": _sha(ROOT / "scripts" / "pressure_fast_domain.py"),
        "validation_script_sha256": _sha(Path(__file__)),
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(
            freqs.astype("<f8").tobytes()
        ).hexdigest(),
        "frequency_range": [float(np.min(freqs)), float(np.max(freqs))],
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "fast": {
            "tol": 1e-6,
            "frequency_res": 1.0,
            "threads": 2,
            "runtime_s": fast_elapsed,
            "dn_eff_input": float(params["DN_eff"]),
            "dn_eff_final": dn_eff,
            "DN_gw": fast_dn,
            "spectrum_log10_range": [float(np.min(fast_log)), float(np.max(fast_log))],
            "accuracy_status": getattr(model, "local_error_budget", {}).get(
                "certification_status"),
            "accuracy_budget": getattr(model, "local_error_budget", None),
            "warnings": [str(w.message) for w in fast_warnings],
        },
        "reference": {
            "fixed_background_DN_eff": dn_eff,
            "cartesian_rtol_1e-10_runtime_s": cartesian_times[1e-10],
            "cartesian_rtol_1e-11_runtime_s": cartesian_times[1e-11],
            "prufer_rtol_1e-9_runtime_s": prufer_elapsed,
            "cartesian_rtol_1e-10_DN_gw": float(cartesian[1e-10][0]["DN_gw"]),
            "cartesian_rtol_1e-11_DN_gw": ref_dn,
            "prufer_rtol_1e-9_DN_gw": float(prufer_result["DN_gw"]),
            "cartesian_map_relative_DN_gw_difference_1e-10_vs_1e-11": reference_dn_rtol_delta,
            "prufer_vs_cartesian_DN_gw_relative": abs(
                float(prufer_result["DN_gw"]) - ref_dn
            ) / max(abs(ref_dn), np.finfo(float).tiny),
            "reference_quadrature_error": strict_result["quadrature_error"],
            "reference_interpolation_error": strict_result["interpolation_error"],
            "cartesian_warnings": cartesian_warnings,
            "prufer_warnings": prufer_warnings,
        },
        "comparison": {
            "positive_spectrum_bins": int(valid_fast_ref.sum()),
            "max_abs_spectrum_delta_dex": float(fast_log_delta.max()),
            "p95_abs_spectrum_delta_dex": float(np.percentile(fast_log_delta, 95)),
            "max_fast_error_index": max_fast_i,
            "max_fast_error_frequency_log10_hz": float(freqs[max_fast_i]),
            "max_fast_error_handoff_eps": float(handoff_eps[max_fast_i]),
            "max_cartesian_1e-10_vs_1e-11_dex": float(ctol_log_delta.max()),
            "max_prufer_vs_cartesian_dex": float(prufer_log_delta.max()),
            "DN_gw_absolute_difference_fast_vs_reference": dn_abs,
            "DN_gw_relative_difference_fast_vs_reference": dn_abs / max(
                abs(fast_dn), abs(ref_dn), np.finfo(float).tiny),
            "DN_eff_float_spacing": float(np.spacing(dn_eff)),
            "reference_total_DN_eff_rounds_to_input": float(dn_eff + ref_dn) == dn_eff,
        },
        "sampled_background_check": {
            "points": int(probe_N.size),
            "H2_all_finite_positive": bool(np.all(np.isfinite(h2)) and np.all(h2 > 0.0)),
            "H2_min": float(np.min(h2)),
            "H2_max": float(np.max(h2)),
        },
        "frequency_records": frecords,
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers_fast": 2,
            "workers_reference": 1,
        },
        "limitations": [
            "The Reference maps use the same 77 Fast-native frequencies and fixed final DN_eff; they certify this discrete map, not frequency-grid convergence or a continuous-spectrum root.",
            "The H2 positivity check samples 513 points; it is diagnostic, not a proof over the continuous interval.",
            "The Fast output is intentionally classified uncertified; no input bound or domain claim is inferred from this one point.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
