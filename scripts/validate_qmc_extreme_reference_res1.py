"""Frequency-grid convergence follow-up for QMC Sobol sample 146."""

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
FREQ_RES = 1.0
Z_TAIL = 5.0
RTOL = 1e-9
OUTPUT = ROOT / "docs" / "qmc_sobol146_reference_res1_prufer_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    item = sobol_cases(POWER, SEED)[SAMPLE_INDEX]
    params = item["params"]
    FS.apply_accuracy_mode("fast")
    fast_model = LCDM_SG(**params)
    fast_out = fast_model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-6, z_tail=Z_TAIL,
        freq_res=1.0, threads=2, accuracy_mode="fast",
        physical_guard=False,
    )
    if fast_out is not fast_model:
        raise RuntimeError(f"Fast sample failed: {fast_model.fast_failure_reason}")
    fast_dn = float(fast_model.DN_gw[-1])
    dn_eff = float(fast_model.cosmo_param["DN_eff"])

    model = LCDM_SG(**params)
    bg = R._Background(model, dn_eff)
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    warm_s = time.perf_counter() - warm_started
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = R.run_reference(
            model, dn_eff=dn_eff, freq_res=FREQ_RES, z_tail=Z_TAIL,
            rtol=RTOL, self_consistent=False, workers=1,
            physical_guard=False, backend="prufer",
        )
    elapsed = time.perf_counter() - started
    ref_dn = float(result["DN_gw"])
    dn_abs = abs(fast_dn - ref_dn)
    cross = json.loads((ROOT / "docs" /
                        "qmc_sobol146_full_prufer_tail_convergence_20261105.json")
                       .read_text(encoding="utf-8"))
    fixed_tail10_dn = float(cross["maps"]["prufer_tail10"]["DN_gw"])
    native77_tail10_dn = float(cross["Fast_vs_independent_reference"]
                               ["sigma_exact_tail10"]["prufer_tail10"]
                               ["Reference_DN_gw"])
    frequencies = np.asarray(result["freqs"], dtype=float)
    ref_omega = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(ref_omega) & (ref_omega > 0.0)):
        raise FloatingPointError("Reference frequency grid contains nonpositive spectrum")

    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_reference_grid_res1",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameters": params,
        "DN_eff_fixed": dn_eff,
        "reference_backend": "prufer",
        "reference_rtol": RTOL,
        "reference_z_tail": Z_TAIL,
        "reference_frequency_res": FREQ_RES,
        "reference_frequency_count": int(frequencies.size),
        "reference_frequency_sha256": hashlib.sha256(
            frequencies.astype("<f8").tobytes()
        ).hexdigest(),
        "reference_frequency_range": [float(np.min(frequencies)),
                                      float(np.max(frequencies))],
        "reference_log10OmegaGW": np.log10(ref_omega).tolist(),
        "fast_native_frequency_count": int(np.asarray(fast_model.f).size),
        "Fast_DN_gw": fast_dn,
        "Reference_res1_DN_gw": ref_dn,
        "Fast_vs_res1_abs_DN_gw_difference": dn_abs,
        "Fast_vs_res1_rel_DN_gw_difference": dn_abs / max(
            abs(fast_dn), abs(ref_dn), np.finfo(float).tiny),
        "Reference_res1_vs_native77_tail10_rel_DN_gw_difference": abs(
            ref_dn - native77_tail10_dn
        ) / max(abs(ref_dn), abs(native77_tail10_dn), np.finfo(float).tiny),
        "Reference_res1_vs_native77_tail10_abs_DN_gw_difference": abs(
            ref_dn - native77_tail10_dn
        ),
        "Reference_res1_vs_fullgrid_tail10_map_relative_difference": abs(
            ref_dn - fixed_tail10_dn
        ) / max(abs(ref_dn), abs(fixed_tail10_dn), np.finfo(float).tiny),
        "reference_quadrature_error": result["quadrature_error"],
        "reference_interpolation_error": result["interpolation_error"],
        "reference_n_iter": result["n_iter"],
        "reference_background_warm_s": warm_s,
        "reference_runtime_s": elapsed,
        "warnings": [str(w.message) for w in caught],
        "fast_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "pressure_script_sha256": _sha(ROOT / "scripts" / "pressure_fast_domain.py"),
        "script_sha256": _sha(Path(__file__)),
        "output": str(OUTPUT),
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "reference_workers": 1,
        },
        "limitations": [
            "This is a fixed-background Reference map on its independent res=1 grid, not a self-consistent root solve.",
            "The frequency grid is denser than the Fast-native grid but a separate res=2 convergence run is not included.",
            "The earlier native-grid Fast-vs-Reference check is the pointwise spectrum test; this record tests integrated frequency-grid sensitivity.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({k: record[k] for k in (
        "reference_frequency_count", "Reference_res1_DN_gw", "Fast_DN_gw",
        "Fast_vs_res1_rel_DN_gw_difference",
        "Reference_res1_vs_native77_tail10_rel_DN_gw_difference",
        "Reference_res1_vs_fullgrid_tail10_map_relative_difference",
        "reference_runtime_s", "warnings", "output")}, indent=2,
        ensure_ascii=False))


if __name__ == "__main__":
    main()
