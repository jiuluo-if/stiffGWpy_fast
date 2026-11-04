"""Compare the Fast goal grid with its existing adaptive grid on Sobol[146]."""

from __future__ import annotations

import argparse
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
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20261105
POWER = 8
SAMPLE_INDEX = 146
FREQ_RES = 1.0
Z_TAIL = 5.0
OUTPUT = ROOT / "docs" / "qmc_sobol146_fast_adaptive_grid_20261105.json"
REFERENCE = ROOT / "docs" / "qmc_sobol146_reference_grid_convergence_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(
    params: dict[str, float], grid: str, target: float, max_points: int,
    z_tail: float, sigma_exact: bool,
) -> dict[str, object]:
    model = LCDM_SG(**params)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        output = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=z_tail,
            freq_res=FREQ_RES,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid=grid,
            freq_grid_target=target,
            freq_grid_max_points=max_points,
            sigma_exact=sigma_exact,
        )
        captured = [
            {
                "message": str(w.message),
                "filename": w.filename,
                "lineno": w.lineno,
                "category": w.category.__name__,
            }
            for w in caught
        ]
    elapsed = time.perf_counter() - started
    if output is not model:
        raise RuntimeError(f"Fast {grid} run failed: {model.fast_failure_reason}")
    freqs = np.asarray(model.f, dtype=float)
    omega = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
    if not np.all(np.isfinite(omega) & (omega > 0.0)):
        raise FloatingPointError(f"Fast {grid} spectrum is not finite and positive")
    return {
        "freq_grid": grid,
        "z_tail": z_tail,
        "sigma_exact": sigma_exact,
        "runtime_s": elapsed,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "frequencies_log10_hz": freqs.tolist(),
        "DN_eff_final": float(model.cosmo_param["DN_eff"]),
        "DN_gw": float(model.DN_gw[-1]),
        "frequency_grid_error": float(model.freq_grid_error),
        "frequency_grid_target": target,
        "frequency_grid_max_points": max_points,
        "quadrature_error_local": float(model.quadrature_error_local),
        "quadrature_error_estimator_method": model.quadrature_error_estimator_method,
        "local_error_budget": model.local_error_budget,
        "log10OmegaGW": np.log10(omega).tolist(),
        "warnings": captured,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--adaptive-target", type=float, action="append",
        help="Adaptive grid error targets to compare (defaults: 3e-4, 1e-4, 3e-5, 1e-5).",
    )
    parser.add_argument("--max-points", type=int, default=1500)
    args = parser.parse_args()
    targets = tuple(args.adaptive_target or (3e-4, 1e-4, 3e-5, 1e-5))
    item = sobol_cases(POWER, SEED)[SAMPLE_INDEX]
    params = item["params"]
    FS.apply_accuracy_mode("fast")
    results = {
        "goal_sigma_tabulated_tail5": _run(
            params, "goal", 3e-4, args.max_points, 5.0, False
        )
    }
    for target in targets:
        key = f"adaptive_tabulated_tail5_{target:.0e}"
        results[key] = _run(params, "adaptive", target, args.max_points, 5.0, False)
        exact_key = f"adaptive_sigma_exact_tail10_{target:.0e}"
        results[exact_key] = _run(params, "adaptive", target, args.max_points, 10.0, True)
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    ref_dn = float(reference["maps"]["res8"]["DN_gw"])
    for result in results.values():
        result["relative_DN_gw_delta_vs_reference_res8"] = abs(
            float(result["DN_gw"]) - ref_dn
        ) / max(abs(float(result["DN_gw"])), abs(ref_dn), np.finfo(float).tiny)
        result["relative_DN_gw_error_vs_reference_res8"] = abs(
            float(result["DN_gw"]) - ref_dn
        ) / max(abs(ref_dn), np.finfo(float).tiny)
    record = {
        "schema_version": 1,
        "experiment": "sobol_20261105_0146_fast_goal_vs_adaptive_grid",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "qmc_seed": SEED,
        "qmc_power": POWER,
        "sample_index": SAMPLE_INDEX,
        "label": item["label"],
        "parameters": params,
        "freq_res": FREQ_RES,
        "z_tail_goal": 5.0,
        "threads": 2,
        "outer_tol": 1e-6,
        "frequency_quadrature": "pchip",
        "adaptive_targets": targets,
        "adaptive_max_points": args.max_points,
        "reference_res8_DN_gw": ref_dn,
        "reference_artifact_sha256": _sha(REFERENCE),
        "fast_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
        "model_api_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
        "reference_module_sha256": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        "script_sha256": _sha(Path(__file__)),
        "results": results,
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "This tests existing Fast grid modes on one extreme background; it does not certify the domain or change production defaults.",
            "The res=8 Reference value is a fixed-background converged map, while these Fast values are self-consistent outputs; their DN_gw comparison is diagnostic rather than a same-map error estimate.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "reference_res8_DN_gw": ref_dn,
        "results": {k: {name: v[name] for name in (
            "frequency_count", "DN_gw", "relative_DN_gw_delta_vs_reference_res8",
            "relative_DN_gw_error_vs_reference_res8", "frequency_grid_error",
            "runtime_s", "warnings")}
        for k, v in results.items()},
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
