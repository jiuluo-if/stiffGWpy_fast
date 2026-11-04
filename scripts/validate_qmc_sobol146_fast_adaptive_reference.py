"""Independently compare the adaptive Fast candidate on its exact frequency map."""

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
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

FAST_RESULT = ROOT / "docs" / "qmc_sobol146_fast_adaptive_grid_20261105.json"
CASE = "adaptive_sigma_exact_tail10_1e-06"
OUTPUT = ROOT / "docs" / "qmc_sobol146_fast_adaptive_reference_20261105.json"
REFERENCE_TAIL = 10.0
REFERENCE_RTOL = 1e-9


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> None:
    source = json.loads(FAST_RESULT.read_text(encoding="utf-8"))
    fast = source["results"][CASE]
    freqs = np.asarray(fast["frequencies_log10_hz"], dtype=float)
    fast_log = np.asarray(fast["log10OmegaGW"], dtype=float)
    fast_dn = float(fast["DN_gw"])
    dn_eff = float(fast["DN_eff_final"])
    if freqs.size != fast["frequency_count"] or fast_log.size != freqs.size:
        raise ValueError("Fast candidate frequency and spectrum sizes do not match")
    frequency_sha = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    if frequency_sha != fast["frequency_sha256"]:
        raise ValueError("Fast candidate frequency array failed its archived hash")

    model = LCDM_SG(**source["parameters"])
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(R._Background(model, dn_eff))
    warm_s = time.perf_counter() - warm_started
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ref = R.run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=freqs,
            z_tail=REFERENCE_TAIL,
            rtol=REFERENCE_RTOL,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    elapsed = time.perf_counter() - started
    ref_freqs = np.asarray(ref["freqs"], dtype=float)
    ref_omega = np.asarray(ref["Ogw"], dtype=float) - np.asarray(ref["Oj"], dtype=float)
    if not np.array_equal(ref_freqs, freqs):
        raise AssertionError("Reference did not preserve the Fast frequency array")
    if not np.all(np.isfinite(ref_omega) & (ref_omega > 0.0)):
        raise FloatingPointError("Reference spectrum is not finite and positive")
    ref_log = np.log10(ref_omega)
    delta_dex = np.abs(fast_log - ref_log)
    ref_dn = float(ref["DN_gw"])
    warnings_rows = [
        {
            "message": str(w.message),
            "filename": w.filename,
            "lineno": w.lineno,
            "category": w.category.__name__,
        }
        for w in caught
    ]
    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_fast_adaptive_same_grid_reference",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "fast_artifact": FAST_RESULT.name,
        "fast_artifact_sha256": _sha(FAST_RESULT),
        "fast_case": CASE,
        "qmc_seed": source["qmc_seed"],
        "qmc_power": source["qmc_power"],
        "sample_index": source["sample_index"],
        "label": source["label"],
        "parameters": source["parameters"],
        "DN_eff_fixed": dn_eff,
        "fast_DN_gw": fast_dn,
        "reference_DN_gw": ref_dn,
        "relative_DN_gw_delta": abs(fast_dn - ref_dn) / max(
            abs(fast_dn), abs(ref_dn), np.finfo(float).tiny
        ),
        "frequency_count": int(freqs.size),
        "frequency_sha256": frequency_sha,
        "frequency_range": [float(freqs.min()), float(freqs.max())],
        "max_abs_spectrum_error_dex": float(np.max(delta_dex)),
        "p95_abs_spectrum_error_dex": float(np.percentile(delta_dex, 95)),
        "median_abs_spectrum_error_dex": float(np.median(delta_dex)),
        "max_error_frequency_log10_hz": float(freqs[int(np.argmax(delta_dex))]),
        "reference_backend": "prufer",
        "reference_rtol": REFERENCE_RTOL,
        "reference_z_tail": REFERENCE_TAIL,
        "reference_workers": 1,
        "reference_warmup_s": warm_s,
        "reference_runtime_s": elapsed,
        "quadrature_error": ref["quadrature_error"],
        "interpolation_error": ref["interpolation_error"],
        "warnings": warnings_rows,
        "fast_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "Fixed-background same-frequency comparison at one extreme QMC point; it does not solve or certify the self-consistent Reference root.",
            "The Fast adaptive frequency grid itself is used for integration; dense independent grid convergence is archived separately.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "fast_DN_gw": fast_dn,
        "reference_DN_gw": ref_dn,
        "relative_DN_gw_delta": record["relative_DN_gw_delta"],
        "max_abs_spectrum_error_dex": record["max_abs_spectrum_error_dex"],
        "p95_abs_spectrum_error_dex": record["p95_abs_spectrum_error_dex"],
        "reference_runtime_s": elapsed,
        "warnings": warnings_rows,
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
