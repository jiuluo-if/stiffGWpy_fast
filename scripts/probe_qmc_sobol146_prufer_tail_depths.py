"""Check Prüfer tail-depth convergence on the remaining QMC error modes."""

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
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20261105
POWER = 8
SAMPLE_INDEX = 146
FREQUENCIES = (-14.101634579145523,)
TAILS = (10.0, 12.0, 15.0)
BACKENDS = (("cartesian", 1e-11), ("prufer", 1e-9))
OUTPUT = ROOT / "docs" / "qmc_sobol146_prufer_tail_depths_20261105.json"


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
    model = LCDM_SG(**params)
    dn_eff = float(params["DN_eff"])
    bg = R._Background(model, dn_eff)
    warm_started = time.perf_counter()
    R._warmup_numba_reference_background(bg)
    warm_s = time.perf_counter() - warm_started

    rows = []
    total_s = 0.0
    for freq in FREQUENCIES:
        solutions = {}
        for backend, rtol in BACKENDS:
            backend_solutions = {}
            for tail in TAILS:
                started = time.perf_counter()
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    solution = R.solve_reference_mode(
                        model, freq, dn_eff, z_tail=tail, z_start=-12.0,
                        rtol=rtol, backend=backend,
                    )
                elapsed = time.perf_counter() - started
                total_s += elapsed
                omega = float(solution["Ogw_today"] - solution["Oj_today"])
                backend_solutions[f"tail{int(tail)}"] = {
                    "elapsed_s": elapsed,
                    "log10Omega_GW": float(np.log10(omega)) if omega > 0 else None,
                    "used_tail": bool(solution["used_tail"]),
                    "event_N": solution["event_N"],
                    "nfev": int(solution["n_steps"]),
                    "matching_error_rel": solution["matching_error_rel"],
                    "warnings": [str(w.message) for w in caught],
                }
            solutions[backend] = backend_solutions
            d10_12 = abs(backend_solutions["tail10"]["log10Omega_GW"]
                         - backend_solutions["tail12"]["log10Omega_GW"])
            d12_15 = abs(backend_solutions["tail12"]["log10Omega_GW"]
                         - backend_solutions["tail15"]["log10Omega_GW"])
            d10_15 = abs(backend_solutions["tail10"]["log10Omega_GW"]
                         - backend_solutions["tail15"]["log10Omega_GW"])
            solutions[backend]["tail_depth_deltas"] = {
                "tail10_vs_tail12_abs_dex": d10_12,
                "tail12_vs_tail15_abs_dex": d12_15,
                "tail10_vs_tail15_abs_dex": d10_15,
            }
        d15_backend = abs(
            solutions["cartesian"]["tail15"]["log10Omega_GW"]
            - solutions["prufer"]["tail15"]["log10Omega_GW"]
        )
        rows.append({
            "frequency_log10_hz": freq,
            "solutions": solutions,
            "cartesian_vs_prufer_tail15_abs_dex": d15_backend,
        })

    result = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_prufer_tail_depth_convergence",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameters": params,
        "DN_eff_fixed": dn_eff,
        "reference_backends_and_rtol": BACKENDS,
        "tail_thresholds": TAILS,
        "jit_background_warm_s": warm_s,
        "total_mode_runtime_s": total_s,
        "fast_reference_crosscheck_sha256": _sha(
            ROOT / "docs" / "qmc_sobol146_full_prufer_tail_convergence_20261105.json"),
        "module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "workers": 1,
        },
        "modes": rows,
        "limitations": [
            "One selected high-handoff-error mode only; this tests tail convergence, not a full root.",
            "The broader map uses the same 77-node Fast grid; independent frequency-grid convergence is separate.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
