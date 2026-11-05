"""Resume the native res2 Reference grid check for the N_inf>120 Fast point."""

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
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "fast_ninf122_tail8_fullgrid_20261105.json"
OUTPUT = ROOT / "docs" / "fast_ninf122_tail8_native_res2_20261105.json"
FREQ_RES = 2.0
RTOL = 1e-9
Z_TAIL = 8.0


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(record: dict) -> None:
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": record["status"],
        "frequency_count": record.get("frequency_count"),
        "output": str(OUTPUT),
    }), flush=True)


def main() -> None:
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    if parent["status"] != "adaptive_fast_grid_prufer_complete":
        raise RuntimeError(f"unexpected parent checkpoint status {parent['status']!r}")
    fast = parent["fast"]
    params = parent["parameters"]
    dn_eff = float(fast["fixed_DN_eff"])
    adaptive_dn = float(parent["reference_maps"]["adaptive_fast_grid"]["DN_gw"])
    model = LCDM_SG(**params)
    record = {
        "schema_version": 1,
        "experiment": "fast_ninf122_tail8_native_freq_res2_fixed_background",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "fast_parent_artifact": PARENT.name,
        "fast_parent_artifact_sha256": parent_sha,
        "fast_parent_case": "adaptive_exact_sigma_tail8",
        "fast_frequency_count": int(fast["frequency_count"]),
        "fast_frequency_sha256": fast["frequency_sha256"],
        "fast_DN_eff_final": dn_eff,
        "fast_DN_gw_self_consistent_output": float(fast["DN_gw_self_consistent_output"]),
        "adaptive_grid_prufer_DN_gw": adaptive_dn,
        "reference": {
            "backend": "prufer",
            "freq_res": FREQ_RES,
            "fixed_DN_eff": dn_eff,
            "z_tail": Z_TAIL,
            "rtol": RTOL,
            "workers": 1,
        },
        "status": "native_res2_running",
    }
    _save(record)

    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_res=FREQ_RES,
            z_tail=Z_TAIL,
            rtol=RTOL,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    elapsed = time.perf_counter() - start
    if caught:
        raise RuntimeError(f"Prüfer Reference emitted warnings: {[str(w.message) for w in caught]}")

    freqs = np.asarray(result["freqs"], dtype=float)
    raw = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(raw)):
        raise FloatingPointError("native res2 Reference returned non-finite data")
    ref_dn = float(result["DN_gw"])
    adaptive_delta = abs(ref_dn - adaptive_dn) / max(
        abs(ref_dn), abs(adaptive_dn), np.finfo(float).tiny
    )
    fast_dn = float(fast["DN_gw_self_consistent_output"])
    fast_delta = abs(fast_dn - ref_dn) / max(
        abs(fast_dn), abs(ref_dn), np.finfo(float).tiny
    )
    record.update({
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": raw.tolist(),
        "DN_gw": ref_dn,
        "quadrature_error": float(result["quadrature_error"]),
        "interpolation_error": float(result["interpolation_error"]),
        "used_tail_count": int(np.count_nonzero(result["used_tail"])),
        "runtime_s": elapsed,
        "relative_DN_gw_delta_vs_adaptive_fast_grid_prufer": adaptive_delta,
        "relative_DN_gw_error_vs_fast_output": fast_delta,
        "acceptance": {
            "native_grid_is_denser": int(freqs.size) > int(fast["frequency_count"]),
            "two_grid_DN_gw_delta_le_2e-4": adaptive_delta <= 2e-4,
            "Fast_DN_gw_vs_native_reference_le_2e-4": fast_delta <= 2e-4,
        },
        "hashes": {
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "driver_py": _sha(Path(__file__)),
        },
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": numba.__version__,
            "workers": 1,
        },
        "limitations": [
            "Fixed-background frequency-grid check only; no Reference self-consistent root is solved or certified.",
            "A two-grid relative change below threshold is a convergence screen, not a formal quadrature error bound.",
        ],
        "status": "complete",
    })
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": int(freqs.size),
        "DN_gw": ref_dn,
        "adaptive_grid_delta": adaptive_delta,
        "Fast_vs_native_reference_delta": fast_delta,
        "runtime_s": elapsed,
        "acceptance": record["acceptance"],
    }, indent=2))


if __name__ == "__main__":
    main()
