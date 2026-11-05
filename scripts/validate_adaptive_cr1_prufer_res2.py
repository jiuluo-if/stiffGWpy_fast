"""Check adaptive positive-cr Fast against the denser native Prüfer grid."""

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
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARENT = ROOT / "docs" / "adaptive_cr1_prufer_fullgrid_20261105.json"
OUTPUT = ROOT / "docs" / "adaptive_cr1_prufer_res2_20261105.json"
FREQ_RES = 2.0
Z_TAIL = 10.0
RTOL = 1e-9


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
    if parent["status"] != "partial_resource_limited":
        raise RuntimeError("parent artifact has unexpected status")
    params = parent["parameters"]
    dn_eff = float(parent["fast"]["DN_eff_final"])
    adaptive_dn = float(parent["reference_maps"]["adaptive_fast_grid"]["DN_gw"])
    model = LCDM_SG(**params)
    record = {
        "schema_version": 1,
        "experiment": "adaptive_cr1_prufer_native_freq_res2_fixed_background_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": params,
        "reference": {
            "backend": "prufer",
            "freq_res": FREQ_RES,
            "z_tail": Z_TAIL,
            "rtol": RTOL,
            "workers": 1,
            "fixed_DN_eff": dn_eff,
        },
        "parent_adaptive_artifact": PARENT.name,
        "parent_adaptive_artifact_sha256": parent_sha,
        "parent_adaptive_frequency_count": int(parent["fast"]["frequency_count"]),
        "parent_adaptive_frequency_sha256": parent["fast"]["frequency_sha256"],
        "parent_adaptive_prufer_DN_gw": adaptive_dn,
        "status": "reference_freq_res_2_running",
    }
    _save(record)

    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        output = run_reference(
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
        raise RuntimeError(f"Prüfer Reference emitted warnings: {[str(w.message) for w in caught]}")

    freqs = np.asarray(output["freqs"], dtype=float)
    raw = np.asarray(output["Ogw"], dtype=float) - np.asarray(output["Oj"], dtype=float)
    if not np.all(np.isfinite(freqs)) or not np.all(np.isfinite(raw)):
        raise FloatingPointError("Reference produced a non-finite frequency or spectrum value")
    freq_hash = hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest()
    ref_dn = float(output["DN_gw"])
    grid_delta = abs(ref_dn - adaptive_dn) / max(
        abs(ref_dn), abs(adaptive_dn), np.finfo(float).tiny
    )
    record.update({
        "frequency_count": int(freqs.size),
        "frequency_sha256": freq_hash,
        "frequencies_log10_hz": freqs.tolist(),
        "omega_gw": raw.tolist(),
        "DN_gw": ref_dn,
        "quadrature_error": float(output["quadrature_error"]),
        "interpolation_error": float(output["interpolation_error"]),
        "used_tail_count": int(np.count_nonzero(output["used_tail"])),
        "runtime_s": elapsed,
        "relative_DN_gw_delta_vs_adaptive_256_node_prufer": grid_delta,
        "acceptance": {
            "reference_grid_is_denser_than_adaptive": freqs.size > int(parent["fast"]["frequency_count"]),
            "relative_DN_gw_grid_delta_le_2e-4": grid_delta <= 2e-4,
        },
        "hashes": {
            "fast_sgwb_py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB_py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "reference_py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script_py": _sha(Path(__file__)),
        },
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "workers": 1,
        },
        "limitations": [
            "This compares two frequency grids on the same fixed background; it is not a self-consistent Reference root or a global-domain certification.",
            "The two-grid delta is a convergence screen, not a formal upper bound on quadrature error.",
        ],
        "status": "complete",
    })
    _save(record)
    print(json.dumps({
        "output": str(OUTPUT),
        "frequency_count": int(freqs.size),
        "DN_gw": ref_dn,
        "relative_DN_gw_delta_vs_adaptive": grid_delta,
        "runtime_s": elapsed,
        "acceptance": record["acceptance"],
    }, indent=2))


if __name__ == "__main__":
    main()
