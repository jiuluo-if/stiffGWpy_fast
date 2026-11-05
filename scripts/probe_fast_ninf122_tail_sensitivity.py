"""Probe Fast tail/sigma/transition choices above the old N_inf=120 pilot cap."""

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
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-2,
    "n_t": 0.5,
    "cr": 0.0,
    "T_re": 1e-3,
    "DN_re": 100.0,
    "kappa10": 1e-2,
}
TARGET = 1e-4
MAX_POINTS = 1500
REFERENCE_RTOL = 1e-9
REFERENCE_WORKERS = 1
REFERENCE_POINTS = 8
OUTPUT = ROOT / "docs" / "fast_ninf122_tail_sensitivity_20261105.json"

CASES = (
    {"name": "default_goal_tabulated_tail5", "freq_grid": None,
     "sigma_exact": False, "z_tail": 5.0, "transition_refine": False},
    {"name": "adaptive_tabulated_tail5", "freq_grid": "adaptive",
     "sigma_exact": False, "z_tail": 5.0, "transition_refine": False},
    {"name": "adaptive_exact_sigma_tail5", "freq_grid": "adaptive",
     "sigma_exact": True, "z_tail": 5.0, "transition_refine": False},
    {"name": "adaptive_exact_sigma_transition_refine_tail5", "freq_grid": "adaptive",
     "sigma_exact": True, "z_tail": 5.0, "transition_refine": True},
    {"name": "adaptive_exact_sigma_tail8", "freq_grid": "adaptive",
     "sigma_exact": True, "z_tail": 8.0, "transition_refine": False},
    {"name": "adaptive_exact_sigma_tail10", "freq_grid": "adaptive",
     "sigma_exact": True, "z_tail": 10.0, "transition_refine": False},
)


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference(model: LCDM_SG, dn_eff: float, freqs: np.ndarray,
               z_tail: float, backend: str) -> tuple[dict, float, list[str]]:
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=freqs,
            z_tail=z_tail,
            rtol=REFERENCE_RTOL,
            self_consistent=False,
            workers=REFERENCE_WORKERS,
            physical_guard=False,
            backend=backend,
        )
    return result, time.perf_counter() - started, [str(w.message) for w in caught]


def main() -> None:
    FS.apply_accuracy_mode("fast")
    results = {}
    for case in CASES:
        model = LCDM_SG(**PARAMS)
        call_args = {}
        if case["freq_grid"] is not None:
            call_args = {
                "freq_grid": case["freq_grid"],
                "freq_grid_target": TARGET,
                "freq_grid_max_points": MAX_POINTS,
            }
        started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = model.SGWB_iter(
                engine="fast",
                fallback=False,
                tol=1e-6,
                z_tail=case["z_tail"],
                freq_res=1.0,
                threads=2,
                accuracy_mode="fast",
                physical_guard=False,
                sigma_exact=case["sigma_exact"],
                transition_refine=case["transition_refine"],
                **call_args,
            )
        fast_elapsed = time.perf_counter() - started
        if out is not model:
            results[case["name"]] = {
                **case,
                "status": "failed",
                "failure_reason": model.fast_failure_reason,
                "runtime_s": fast_elapsed,
                "warnings": [str(w.message) for w in caught],
            }
            continue

        freqs = np.asarray(model.f, dtype=float)
        fast_raw = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
        if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(fast_raw)):
            raise FloatingPointError(f"{case['name']} returned a non-finite grid or spectrum")
        indices = np.unique(np.linspace(0, freqs.size - 1, min(REFERENCE_POINTS, freqs.size), dtype=int))
        ref_freqs = freqs[indices]
        ref_model = LCDM_SG(**PARAMS)
        refs = {}
        budget = getattr(model, "local_error_budget", None)
        for backend in ("cartesian", "prufer"):
            reference, elapsed, ref_warnings = _reference(
                ref_model,
                float(model.cosmo_param["DN_eff"]),
                ref_freqs,
                float(case["z_tail"]),
                backend,
            )
            ref_raw = np.asarray(reference["Ogw"], dtype=float) - np.asarray(reference["Oj"], dtype=float)
            valid = (fast_raw[indices] > 0.0) & (ref_raw > 0.0)
            if not np.all(valid):
                raise FloatingPointError(f"{case['name']} has non-positive compared bins for {backend}")
            errors = np.log10(fast_raw[indices] / ref_raw)
            refs[backend] = {
                "runtime_s": elapsed,
                "DN_gw_on_subset_not_used_as_integral": float(reference["DN_gw"]),
                "omega_gw": ref_raw.tolist(),
                "warnings": ref_warnings,
                "spectrum_error_dex": errors.tolist(),
                "max_abs_spectrum_error_dex": float(np.max(np.abs(errors))),
                "p95_abs_spectrum_error_dex": float(np.quantile(np.abs(errors), 0.95)),
            }
        cart_raw = np.asarray(refs["cartesian"]["omega_gw"], dtype=float)
        prufer_raw = np.asarray(refs["prufer"]["omega_gw"], dtype=float)
        cross_error = np.log10(cart_raw / prufer_raw)
        reference_crosscheck = {
            "cartesian_vs_prufer_max_abs_dex": float(np.max(np.abs(cross_error))),
            "cartesian_vs_prufer_p95_abs_dex": float(np.quantile(np.abs(cross_error), 0.95)),
            "cartesian_vs_prufer_error_dex": cross_error.tolist(),
        }

        results[case["name"]] = {
            **case,
            "status": "complete",
            "fast_runtime_s": fast_elapsed,
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_raw.tolist(),
            "DN_eff_final": float(model.cosmo_param["DN_eff"]),
            "DN_gw": float(model.DN_gw[-1]),
            "frequency_grid_error": float(model.freq_grid_error),
            "frequency_grid_target": TARGET if case["freq_grid"] is not None else None,
            "frequency_grid_max_points": MAX_POINTS if case["freq_grid"] is not None else None,
            "local_error_status": (
                budget.get("certification_status") if isinstance(budget, dict) else None
            ),
            "all_spectrum_bins_positive": bool(np.all(fast_raw > 0.0)),
            "warnings": [str(w.message) for w in caught],
            "reference_indices": indices.tolist(),
            "reference_frequency_sha256": hashlib.sha256(ref_freqs.astype("<f8").tobytes()).hexdigest(),
            "reference": refs,
            "reference_crosscheck": reference_crosscheck,
        }

    record = {
        "schema_version": 1,
        "experiment": "fast_ninf122_adaptive_tail_sigma_transition_sensitivity",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMS,
        "derived_N_inf": float(LCDM_SG(**PARAMS).derived_param["N_inf"]),
        "fast_profile": "fast",
        "fast_threads": 2,
        "reference_rtol": REFERENCE_RTOL,
        "reference_workers": REFERENCE_WORKERS,
        "reference_frequency_samples_per_grid": REFERENCE_POINTS,
        "adaptive_target": TARGET,
        "adaptive_max_points": MAX_POINTS,
        "results": results,
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "hashes": {
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "driver_py": _sha(Path(__file__)),
        },
        "limitations": [
            "Reference checks cover eight stratified frequencies per configuration, not a full spectrum or bolometric Reference map.",
            "Frequency-grid convergence and a self-consistent Reference root have not been certified for this N_inf>120 point.",
            "No production defaults or solver equations are changed by this probe.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "N_inf": record["derived_N_inf"],
        "results": {name: {
            "status": value["status"],
            "frequency_count": value.get("frequency_count"),
            "DN_gw": value.get("DN_gw"),
            "grid_error": value.get("frequency_grid_error"),
            "max_cartesian_abs_dex": value.get("reference", {}).get("cartesian", {}).get("max_abs_spectrum_error_dex"),
            "max_prufer_abs_dex": value.get("reference", {}).get("prufer", {}).get("max_abs_spectrum_error_dex"),
            "cartesian_vs_prufer_max_abs_dex": value.get("reference_crosscheck", {}).get("cartesian_vs_prufer_max_abs_dex"),
        } for name, value in results.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
