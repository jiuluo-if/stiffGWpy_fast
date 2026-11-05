"""Confirm the low-T positive-cr Reference grid at native freq_res=8."""

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
BASELINE = ROOT / "docs" / "fast_cr1_lowT_fullgrid_prufer_20261105.json"
PARENT = ROOT / "docs" / "fast_cr1_lowT_res4_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_cr1_lowT_res8_prufer_20261105.json"
FREQ_RES = 8.0
Z_TAIL = 8.0
RTOL = 1e-9
DN_REL_TOL = 2e-4
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    baseline_sha, parent_sha = _sha(BASELINE), _sha(PARENT)
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    model = LCDM_SG(**baseline["parameters_input"])
    dn_eff = float(baseline["fast"]["DN_eff_final"])
    record = {
        "schema_version": 1,
        "experiment": "fast_cr1_lowT_prufer_native_res8_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters_input": baseline["parameters_input"],
        "derived": baseline["derived"],
        "fixed_DN_eff": dn_eff,
        "baseline_artifact": BASELINE.name,
        "baseline_artifact_sha256": baseline_sha,
        "res4_artifact": PARENT.name,
        "res4_artifact_sha256": parent_sha,
        "reference": {"backend": "prufer", "freq_res": FREQ_RES, "z_tail": Z_TAIL, "rtol": RTOL, "workers": 1},
        "status": "native_res8_running",
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": record["status"], "output": str(OUTPUT)}), flush=True)

    started = time.perf_counter()
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
    elapsed = time.perf_counter() - started
    if caught:
        raise RuntimeError(f"res8 Prüfer emitted warnings: {[str(item.message) for item in caught]}")

    freqs = np.asarray(result["freqs"], dtype=float)
    spectrum = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(spectrum)) or not np.all(spectrum > 0.0):
        raise FloatingPointError("res8 Prüfer returned invalid values")
    dn_res8 = float(result["DN_gw"])
    dn_res4 = float(parent["native_res4"]["DN_gw"])
    dn_res2 = float(baseline["reference_maps"]["native_res2"]["DN_gw"])
    dn_adaptive = float(baseline["reference_maps"]["adaptive_exact_fast_grid"]["DN_gw"])
    dn_fast = float(baseline["fast"]["DN_gw"])
    def rel(value: float) -> float:
        return abs(value - dn_res8) / max(abs(value), abs(dn_res8), np.finfo(float).tiny)
    record.update({
        "native_res8": {
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": spectrum.tolist(),
            "DN_gw": dn_res8,
            "quadrature_error": float(result["quadrature_error"]),
            "interpolation_error": float(result["interpolation_error"]),
            "used_tail_count": int(np.count_nonzero(result["used_tail"])),
            "runtime_s": elapsed,
        },
        "comparisons": {
            "res2_to_res8_DN_gw_relative_delta": rel(dn_res2),
            "res4_to_res8_DN_gw_relative_delta": rel(dn_res4),
            "adaptive_to_res8_DN_gw_relative_delta": rel(dn_adaptive),
            "fast_to_res8_DN_gw_relative_error": rel(dn_fast),
        },
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "numba": numba.__version__,
            "reference_workers": 1,
        },
        "hashes": {
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script_py": _sha(Path(__file__)),
        },
        "acceptance": {
            "res4_to_res8_rel_le_2e-4": rel(dn_res4) <= DN_REL_TOL,
            "fast_to_res8_rel_le_2e-4": rel(dn_fast) <= DN_REL_TOL,
        },
        "limitations": [
            "This is fixed-background frequency-grid convergence; it does not solve the independent Reference self-consistent root.",
            "Agreement of res4 and res8 is a convergence screen, not a formal error enclosure.",
        ],
        "status": "native_res8_complete",
    })
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "frequency_count": record["native_res8"]["frequency_count"],
        "DN_gw": dn_res8,
        **record["comparisons"],
        "acceptance": record["acceptance"],
        "runtime_s": elapsed,
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
