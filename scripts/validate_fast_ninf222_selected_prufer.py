"""Cross-check 12 stratified Fast modes at a second high-N_inf point."""

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

import numba
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf122_tail8_target1e5_validation_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf222_selected_prufer_20261105.json"
SELECTED_MODE_COUNT = 12
Z_TAIL = 8.0
RTOL = 1e-9


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    params = {**parent["parameters"], "DN_re": 200.0}
    FS.apply_accuracy_mode("fast")
    model = LCDM_SG(**params)

    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=Z_TAIL,
            freq_res=1.0,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid="adaptive",
            freq_grid_target=1e-5,
            freq_grid_max_points=1500,
            sigma_exact=True,
        )
    fast_elapsed = time.perf_counter() - started
    if out is not model or not model.SGWB_converge:
        raise RuntimeError(f"Fast failed: {model.fast_failure_reason}")
    if caught:
        raise RuntimeError(f"Fast emitted warnings: {[str(item.message) for item in caught]}")

    freqs = np.asarray(model.f, dtype=float)
    fast_spectrum = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(fast_spectrum)) or not np.all(fast_spectrum > 0.0):
        raise FloatingPointError("Fast returned invalid frequencies or spectrum values")
    parent_freqs = np.asarray(parent["fast"]["frequencies_log10_hz"], dtype=float)
    if not np.array_equal(freqs, parent_freqs):
        raise AssertionError("candidate Fast frequency array differs from parent grid")

    indices = np.linspace(0, freqs.size - 1, SELECTED_MODE_COUNT).astype(int)
    selected_freqs = freqs[indices]
    reference_started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=float(model.cosmo_param["DN_eff"]),
            freq_res=1.0,
            z_tail=Z_TAIL,
            rtol=RTOL,
            freq_subset=selected_freqs,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    reference_elapsed = time.perf_counter() - reference_started
    if caught:
        raise RuntimeError(f"Prüfer Reference emitted warnings: {[str(item.message) for item in caught]}")

    reference_freqs = np.asarray(result["freqs"], dtype=float)
    reference_spectrum = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.array_equal(reference_freqs, selected_freqs):
        raise AssertionError("Reference did not preserve the selected frequency array")
    if not np.all(np.isfinite(reference_spectrum)) or not np.all(reference_spectrum > 0.0):
        raise FloatingPointError("Prüfer Reference returned invalid selected spectrum values")

    fast_error_dex = np.log10(fast_spectrum[indices] / reference_spectrum)
    parent_reference = np.asarray(parent["reference_maps"]["adaptive_fast_grid"]["omega_gw"], dtype=float)
    parent_delta_dex = np.log10(reference_spectrum / parent_reference[indices])
    maximum_fast_error = float(np.max(np.abs(fast_error_dex)))
    maximum_parent_delta = float(np.max(np.abs(parent_delta_dex)))
    if maximum_fast_error > 1e-3:
        raise AssertionError(f"selected-mode Fast spectrum error {maximum_fast_error:g} dex exceeds 1e-3")

    frequency_sha = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    selected_frequency_sha = hashlib.sha256(selected_freqs.astype("<f8").tobytes()).hexdigest()
    record = {
        "schema_version": 1,
        "experiment": "fast_ninf222_selected_prufer_fixed_background_screen",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "N_inf": float(model.derived_param["N_inf"]),
        "parent_artifact": PARENT.name,
        "parent_artifact_sha256": parent_sha,
        "fixed_DN_eff": float(model.cosmo_param["DN_eff"]),
        "fast": {
            "grid": "adaptive",
            "grid_target_dex": 1e-5,
            "z_tail": Z_TAIL,
            "sigma_exact": True,
            "frequency_count": int(freqs.size),
            "frequency_sha256": frequency_sha,
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": fast_spectrum.tolist(),
            "DN_eff_final": float(model.cosmo_param["DN_eff"]),
            "DN_gw": float(model.DN_gw[-1]),
            "frequency_grid_error": float(model.freq_grid_error),
            "runtime_s": fast_elapsed,
        },
        "reference": {
            "backend": "prufer",
            "rtol": RTOL,
            "z_tail": Z_TAIL,
            "workers": 1,
            "selected_indices": indices.tolist(),
            "frequency_sha256": selected_frequency_sha,
            "frequencies_log10_hz": selected_freqs.tolist(),
            "omega_gw": reference_spectrum.tolist(),
            "runtime_s": reference_elapsed,
        },
        "comparison": {
            "selected_mode_count": int(indices.size),
            "fast_vs_prufer_max_abs_dex": maximum_fast_error,
            "fast_vs_prufer_p95_abs_dex": float(np.quantile(np.abs(fast_error_dex), 0.95)),
            "candidate_prufer_vs_parent_prufer_max_abs_dex": maximum_parent_delta,
            "spectrum_gate_1e-3_dex": maximum_fast_error <= 1e-3,
            "same_frequency_grid_as_parent": True,
            "finite_positive_fast_and_reference": True,
            "warnings": [],
        },
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": numba.__version__,
            "fast_threads": 2,
            "reference_workers": 1,
        },
        "hashes": {
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script_py": _sha(Path(__file__)),
        },
        "limitations": [
            "This is a 12-mode fixed-background spectrum screen, not a full-spectrum comparison.",
            "DN_gw integration and the self-consistent root are not independently verified at this point.",
            "The near-identical candidate and parent selected Reference spectra are consistent with an e-fold coordinate shift, but do not prove global invariance.",
        ],
        "status": "selected_mode_screen_passed_fixed_background_only",
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "N_inf": record["N_inf"],
        "frequency_count": record["fast"]["frequency_count"],
        "selected_mode_count": record["comparison"]["selected_mode_count"],
        "fast_vs_prufer_max_abs_dex": maximum_fast_error,
        "candidate_prufer_vs_parent_prufer_max_abs_dex": maximum_parent_delta,
        "DN_gw": record["fast"]["DN_gw"],
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
