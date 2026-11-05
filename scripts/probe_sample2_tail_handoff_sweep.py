"""Screen uniform Fast WKB handoff depths on the sample2 adaptive grid."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
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
Z_TAILS = (5.0, 6.0, 7.0, 8.0)
TARGET = 1e-7
MAX_POINTS = 1500
REFERENCE_ROOT = 2.291411867022825
OUTPUT = ROOT / "docs" / "sample2_tail_handoff_sweep_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def _finite_list(values):
    arr = np.asarray(values, dtype=np.float64)
    return [float(value) if np.isfinite(value) else None for value in arr]


def main() -> int:
    record = {
        "schema_version": 1,
        "experiment": "sample2 Fast adaptive grid uniform WKB handoff depth sweep",
        "status": "running",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMS,
        "reference_root_at_z_tail5": REFERENCE_ROOT,
        "fast_configuration": {
            "accuracy_mode": "fast", "freq_res": 1.0, "freq_grid": "adaptive",
            "freq_grid_target": TARGET, "freq_grid_max_points": MAX_POINTS,
            "tol": 1e-6, "threads": 1, "physical_guard": False,
        },
        "z_tail_values": list(Z_TAILS),
        "runs": [],
        "runner_sha256": _sha(Path(__file__).resolve()),
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    for z_tail in Z_TAILS:
        model = LCDM_SG(**PARAMS)
        start = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-6, z_tail=z_tail,
                freq_res=1.0, threads=1, accuracy_mode="fast", physical_guard=False,
                freq_grid="adaptive", freq_grid_target=TARGET,
                freq_grid_max_points=MAX_POINTS,
            )
        runtime = time.perf_counter() - start
        row = {
            "z_tail": z_tail,
            "success": result is model,
            "runtime_s": runtime,
            "warnings": [str(item.message) for item in caught],
            "failure_reason": getattr(model, "fast_failure_reason", None),
        }
        if result is model:
            freqs = np.asarray(model.f, dtype=np.float64)
            logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
            handoff_eps = np.asarray(model.handoff_eps, dtype=np.float64)
            tail_mask = handoff_eps >= 0.0
            row.update({
                "DN_eff_final": float(model.cosmo_param["DN_eff"]),
                "DN_gw_final": float(model.DN_gw[-1]),
                "root_relative_delta_vs_z_tail5_reference": abs(float(model.cosmo_param["DN_eff"]) - REFERENCE_ROOT) / REFERENCE_ROOT,
                "frequency_count": int(freqs.size),
                "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
                "freq_grid_error": float(model.freq_grid_error),
                "tail_mode_count": int(np.count_nonzero(tail_mask)),
                "handoff_eps_max": float(np.max(handoff_eps[tail_mask])) if np.any(tail_mask) else None,
                "handoff_eps_p95": float(np.percentile(handoff_eps[tail_mask], 95)) if np.any(tail_mask) else None,
                "frequency_log10_hz": _finite_list(freqs),
                "log10_OmegaGW": _finite_list(logs),
                "handoff_eps": _finite_list(handoff_eps),
            })
        record["runs"].append(row)
        _write_record = {**record, "status": "running"}
        OUTPUT.write_text(json.dumps(_write_record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    record.update({
        "status": "complete",
        "source_hashes": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB.py": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "runner.py": _sha(Path(__file__).resolve()),
        },
        "resources": telemetry(workers=1, threads=1),
        "timing_limit": "Sequential exploratory timings; the first call can include process/JIT warmup and no paired P50/P95 claim is made.",
        "accuracy_limit": "Compares Fast roots to the z_tail=5 800-point Reference root; no exact-node Reference map at each handoff depth is included in this screen.",
    })
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": record["status"],
        "runs": [{k: row.get(k) for k in ("z_tail", "success", "DN_eff_final", "root_relative_delta_vs_z_tail5_reference", "frequency_count", "handoff_eps_max", "runtime_s")} for row in record["runs"]],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
