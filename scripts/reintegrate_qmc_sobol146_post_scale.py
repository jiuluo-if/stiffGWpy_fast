"""Reintegrate archived native-grid Prüfer modes with the current integrator."""

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
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SOURCE = ROOT / "docs" / "qmc_sobol146_full_prufer_tail_convergence_20261105.json"
OUTPUT = ROOT / "docs" / "qmc_sobol146_post_scale_reintegration_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    rows = source["per_mode"]
    freqs = np.asarray([row["frequency_log10_hz"] for row in rows], dtype=float)
    if freqs.size != source["frequency_count"]:
        raise ValueError("archived per-mode count does not match frequency_count")
    if hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest() != source["frequency_sha256"]:
        raise ValueError("archived per-mode frequencies do not match their hash")
    model = LCDM_SG(**source["parameters"])
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    results = {}
    warnings_seen = []
    for tail in source["reference_tail_thresholds"]:
        key = f"prufer_tail{int(tail)}"
        ogw = np.asarray([row["solutions"][key]["Ogw_today"] for row in rows], dtype=float)
        oj = np.asarray([row["solutions"][key]["Oj_today"] for row in rows], dtype=float)
        if not np.all(np.isfinite(ogw - oj) & (ogw - oj > 0.0)):
            raise FloatingPointError(f"archived {key} spectrum is not finite and positive")
        started = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            g2, qerr, ierr = R.integrate_spectrum(freqs, ogw, oj)
        elapsed = time.perf_counter() - started
        warning_rows = [
            {
                "message": str(w.message),
                "filename": w.filename,
                "lineno": w.lineno,
                "category": w.category.__name__,
            }
            for w in caught
        ]
        warnings_seen.extend(warning_rows)
        new_dn = float(gp.Neff0 * g2 / omega_nu)
        old_dn = float(source["maps"][key]["DN_gw"])
        results[key] = {
            "frequency_count": int(freqs.size),
            "frequency_sha256": source["frequency_sha256"],
            "old_pre_scale_DN_gw": old_dn,
            "current_DN_gw": new_dn,
            "relative_delta": abs(new_dn - old_dn) / max(
                abs(new_dn), abs(old_dn), np.finfo(float).tiny
            ),
            "quadrature_error": qerr,
            "interpolation_error": ierr,
            "runtime_s": elapsed,
            "warnings": warning_rows,
        }
    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_post_scale_reintegration",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "source_artifact": SOURCE.name,
        "source_artifact_sha256": _sha(SOURCE),
        "source_script_sha256": source["script_sha256"],
        "parameters": source["parameters"],
        "DN_eff_fixed": source["reference_background_DN_eff"],
        "reference_backend": "prufer archived modes, reintegrated only",
        "reference_rtol": source["reference_rtol"],
        "reference_z_start": source["reference_z_start"],
        "reference_tail_thresholds": source["reference_tail_thresholds"],
        "frequency_count": int(freqs.size),
        "frequency_sha256": source["frequency_sha256"],
        "frequency_range": [float(freqs.min()), float(freqs.max())],
        "results": results,
        "warnings": warnings_seen,
        "current_reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "Reuses archived per-mode ODE solutions; it independently reruns only the frequency integration.",
            "This verifies scale-safe integration on the native 77-node grid, not grid convergence, a self-consistent root, or Fast-domain certification.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"results": results, "warnings": warnings_seen, "output": str(OUTPUT)}, indent=2))


if __name__ == "__main__":
    main()
