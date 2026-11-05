"""Screen the existing Fast adaptive frequency grid on the minimized sample2 case."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
import traceback
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-09,
    "r": 0.19649188544593907,
    "cr": 0.0,
    "n_t": 1.0622501865308032,
    "T_re": 1761.314039059193,
    "kappa10": 0.011598483151630803,
    "DN_re": 10.0,
    "DN_eff": 0.0,
}
TARGETS = (None, 1e-4, 1e-5, 1e-6, 1e-7)
MAX_POINTS = 1500
OUTPUT = ROOT / "docs" / "sample2_fast_adaptive_sweep_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _finite_list(values) -> list[float | None]:
    arr = np.asarray(values, dtype=np.float64)
    return [float(value) if np.isfinite(value) else None for value in arr]


def _write(record: dict) -> None:
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    record = {
        "schema_version": 1,
        "experiment": "sample2 Fast goal versus built-in adaptive frequency-grid target sweep",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMS,
        "fast_configuration": {
            "accuracy_mode": "fast",
            "tol": 1e-6,
            "z_tail": 5.0,
            "freq_res": 1.0,
            "threads": 1,
            "physical_guard": False,
            "freq_grid_max_points": MAX_POINTS,
        },
        "target_dex_values": list(TARGETS),
        "runs": [],
        "runner_sha256": _sha(Path(__file__).resolve()),
    }
    _write(record)
    for target in TARGETS:
        model = LCDM_SG(**PARAMS)
        started = time.perf_counter()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                kwargs = {
                    "engine": "fast",
                    "fallback": False,
                    "tol": 1e-6,
                    "z_tail": 5.0,
                    "freq_res": 1.0,
                    "threads": 1,
                    "accuracy_mode": "fast",
                    "physical_guard": False,
                }
                if target is not None:
                    kwargs.update(
                        freq_grid="adaptive",
                        freq_grid_target=float(target),
                        freq_grid_max_points=MAX_POINTS,
                    )
                result = model.SGWB_iter(**kwargs)
                elapsed = time.perf_counter() - started
                row = {
                    "target_dex": target,
                    "frequency_grid": "goal" if target is None else "adaptive",
                    "success": result is model,
                    "runtime_s": elapsed,
                    "failure_reason": getattr(model, "fast_failure_reason", None),
                    "warnings": [str(item.message) for item in caught],
                }
            if result is model:
                freqs = np.asarray(model.f, dtype=np.float64)
                logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
                row.update({
                    "DN_eff_final": float(model.cosmo_param["DN_eff"]),
                    "DN_gw_final": float(model.DN_gw[-1]),
                    "frequency_count": int(freqs.size),
                    "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
                    "frequency_log10_hz": _finite_list(freqs),
                    "log10_OmegaGW": _finite_list(logs),
                    "freq_grid_error": float(model.freq_grid_error),
                    "freq_grid_used": model.freq_grid_used,
                    "phase_max_used": float(model.phase_max_used),
                })
            record["runs"].append(row)
            record["status"] = "running"
            _write(record)
        except Exception as exc:
            record["status"] = "failed"
            record["error"] = f"{type(exc).__name__}: {exc}"
            record["traceback"] = traceback.format_exc()
            record["failed_target_dex"] = target
            record["runs"].append({"target_dex": target, "runtime_s": time.perf_counter() - started, "status": "exception"})
            _write(record)
            raise

    record.update({
        "status": "complete",
        "resources": telemetry(workers=1, threads=1),
        "source_hashes": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB.py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "timing_limit": "Targets run sequentially in one process; these are feasibility timings, not a paired P50/P95 performance comparison.",
        "accuracy_limit": "Fast-only sweep. It does not compare spectra or self-consistent roots to an independent Reference; a best target must be validated separately before any claim.",
    })
    _write(record)
    print(json.dumps({
        "status": record["status"],
        "runs": [{k: row.get(k) for k in ("target_dex", "success", "frequency_count", "DN_gw_final", "runtime_s", "freq_grid_error")} for row in record["runs"]],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
