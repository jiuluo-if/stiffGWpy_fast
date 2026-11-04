"""Cross-check compiled Reference backgrounds on a positive-cr cosmology.

This is an evaluator-equivalence and Cartesian/Pruefer cross-check, not a
full-grid or global-domain certification.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "r": 1e-2,
    "cr": 1,
    "T_re": 2e3,
    "kappa10": 1e-2,
    "DN_eff": 0.133732083004284,
}
FREQUENCIES = [-14.6013, -12.6723]
Z_TAIL = 10.0
Z_START = -12.0
RTOL = 1e-11
OUTPUT = ROOT / "docs" / "reference_h2_numba_cr_positive_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _observable(sol: dict) -> float:
    return float(sol["Ogw_today"] - sol["Oj_today"])


def main() -> None:
    model = LCDM_SG(**PARAMS)
    dn_eff = float(PARAMS["DN_eff"])
    bg = R._Background(model, dn_eff)
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    warm_s = time.perf_counter() - warm_started

    python_eval = R._H2_and_sigma_python
    compiled_eval = R._H2_and_sigma
    records = []
    for index, freq in enumerate(FREQUENCIES):
        methods = [
            ("python", "cartesian"),
            ("compiled", "cartesian"),
            ("python", "prufer"),
            ("compiled", "prufer"),
        ]
        if index % 2:
            methods.reverse()
        runs = {}
        for evaluator, backend in methods:
            R._H2_and_sigma = python_eval if evaluator == "python" else compiled_eval
            started = time.perf_counter()
            solution = R.solve_reference_mode(
                model, freq, dn_eff, z_tail=Z_TAIL, z_start=Z_START,
                rtol=RTOL, backend=backend,
            )
            elapsed = time.perf_counter() - started
            runs[f"{backend}_{evaluator}"] = {
                "elapsed_s": elapsed,
                "Ogw_today": solution["Ogw_today"],
                "Oj_today": solution["Oj_today"],
                "Opgw_today": solution["Opgw_today"],
                "Omega_GW_today": _observable(solution),
                "used_tail": solution["used_tail"],
                "event_N": solution["event_N"],
                "n_steps_nfev": solution["n_steps"],
                "z0": solution["z0"],
            }

        for backend in ("cartesian", "prufer"):
            py = runs[f"{backend}_python"]
            nb = runs[f"{backend}_compiled"]
            for key in ("Ogw_today", "Oj_today", "Opgw_today",
                        "Omega_GW_today", "used_tail", "event_N",
                        "n_steps_nfev", "z0"):
                if py[key] != nb[key]:
                    raise AssertionError(
                        f"{backend} evaluator mismatch at freq={freq}: {key}"
                    )

        cart = runs["cartesian_compiled"]
        pruf = runs["prufer_compiled"]
        a, b = cart["Omega_GW_today"], pruf["Omega_GW_today"]
        relative_delta = abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)
        dex_delta = (abs(math.log10(a) - math.log10(b))
                     if a > 0.0 and b > 0.0 else None)
        records.append({
            "frequency_log10_hz": freq,
            "runs": runs,
            "cartesian_vs_prufer_omega_relative_delta": relative_delta,
            "cartesian_vs_prufer_omega_delta_dex": dex_delta,
        })

    R._H2_and_sigma = python_eval
    result = {
        "schema_version": 1,
        "experiment": "reference_h2_numba_cr_positive_crosscheck",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "parameters": PARAMS,
        "frequency_log10_hz": FREQUENCIES,
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "atol": [1e-12, 1e-22, 1e-22],
        "compiled_warmup_s": warm_s,
        "results": records,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
            "concurrent_processes": 1,
        },
        "limitations": [
            "Two frequencies on one positive-cr cosmology only; not a full-grid certification.",
            "Python and compiled H2/sigma paths are expected to be bitwise equivalent; Prüfer is a separate propagation formulation used as a cross-check, not exact truth.",
            "Per-mode timings are diagnostic observations and are not a performance claim.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
