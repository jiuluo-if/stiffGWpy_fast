"""Convergence check for the borderline positive-cr low-T frequency map."""

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
PARENT = ROOT / "docs" / "fast_cr1_lowT_fullgrid_prufer_20261105.json"
OUTPUT = ROOT / "docs" / "fast_cr1_lowT_res4_prufer_20261105.json"
FREQ_RES = 4.0
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
    parent_sha = _sha(PARENT)
    parent = json.loads(PARENT.read_text(encoding="utf-8"))
    model = LCDM_SG(**parent["parameters_input"])
    dn_eff = float(parent["fast"]["DN_eff_final"])
    record = {
        "schema_version": 1,
        "experiment": "fast_cr1_lowT_prufer_native_res4_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters_input": parent["parameters_input"],
        "derived": parent["derived"],
        "fixed_DN_eff": dn_eff,
        "parent_artifact": PARENT.name,
        "parent_artifact_sha256": parent_sha,
        "reference": {"backend": "prufer", "freq_res": FREQ_RES, "z_tail": Z_TAIL, "rtol": RTOL, "workers": 1},
        "status": "native_res4_running",
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
        raise RuntimeError(f"res4 Prüfer emitted warnings: {[str(item.message) for item in caught]}")

    freqs = np.asarray(result["freqs"], dtype=float)
    spectrum = np.asarray(result["Ogw"], dtype=float) - np.asarray(result["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(spectrum)) or not np.all(spectrum > 0.0):
        raise FloatingPointError("res4 Prüfer returned invalid values")
    dn_res4 = float(result["DN_gw"])
    dn_adaptive = float(parent["reference_maps"]["adaptive_exact_fast_grid"]["DN_gw"])
    dn_res2 = float(parent["reference_maps"]["native_res2"]["DN_gw"])
    dn_fast = float(parent["fast"]["DN_gw"])
    denominator = max(abs(dn_res4), abs(dn_res2), np.finfo(float).tiny)
    record.update({
        "native_res4": {
            "frequency_count": int(freqs.size),
            "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
            "frequencies_log10_hz": freqs.tolist(),
            "omega_gw": spectrum.tolist(),
            "DN_gw": dn_res4,
            "quadrature_error": float(result["quadrature_error"]),
            "interpolation_error": float(result["interpolation_error"]),
            "used_tail_count": int(np.count_nonzero(result["used_tail"])),
            "runtime_s": elapsed,
        },
        "comparisons": {
            "res2_to_res4_DN_gw_relative_delta": abs(dn_res2 - dn_res4) / denominator,
            "adaptive_to_res4_DN_gw_relative_delta": abs(dn_adaptive - dn_res4) / max(
                abs(dn_adaptive), abs(dn_res4), np.finfo(float).tiny
            ),
            "fast_to_res4_DN_gw_relative_error": abs(dn_fast - dn_res4) / max(
                abs(dn_fast), abs(dn_res4), np.finfo(float).tiny
            ),
        },
        "acceptance": {},
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
        "limitations": [
            "This extends frequency-grid convergence at a fixed background only; no self-consistent Reference root is solved.",
            "A res2-to-res4 agreement screen is not a formal quadrature error bound.",
        ],
    })
    record["acceptance"] = {
        "res2_to_res4_rel_le_2e-4": record["comparisons"]["res2_to_res4_DN_gw_relative_delta"] <= DN_REL_TOL,
        "fast_to_res4_rel_le_2e-4": record["comparisons"]["fast_to_res4_DN_gw_relative_error"] <= DN_REL_TOL,
    }
    record["status"] = "native_res4_complete"
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "frequency_count": record["native_res4"]["frequency_count"],
        "DN_gw": dn_res4,
        **record["comparisons"],
        "acceptance": record["acceptance"],
        "runtime_s": elapsed,
        "output": str(OUTPUT),
    }), flush=True)


if __name__ == "__main__":
    main()
