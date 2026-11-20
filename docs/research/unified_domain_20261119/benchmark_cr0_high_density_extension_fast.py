"""Warm persistent-process Fast timings for the CR0 high-density screen."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import psutil

ROOT = Path(__file__).resolve().parents[3]
process = psutil.Process()
affinity = process.cpu_affinity()[:2]
process.cpu_affinity(affinity)
for key, value in {
    "NUMBA_NUM_THREADS": "2",
    "NUMBA_THREADING_LAYER": "workqueue",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ[key] = value

sys.path.insert(0, str(ROOT))
from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

REPEATS = 64
BASE_PARAMS = {
    "Omega_bh2": 1.0e6,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = {
    "old_region_lower_edge": {"Omega_bh2": 1.0e6, "n_t": -0.5},
    "new_region_upper_density_blue": {"Omega_bh2": 2.0e6, "n_t": -0.5},
    "new_region_upper_density_zero": {"Omega_bh2": 2.0e6, "n_t": 0.0},
    "new_region_mid_density": {"Omega_bh2": 1.5e6, "n_t": -0.25},
    "new_region_sobol_holdout": {
        "Omega_bh2": 1.5929892305284739e6,
        "n_t": -0.1713485741056502,
    },
    "old_region_holdout": {"Omega_bh2": 7.5e5, "n_t": -0.125},
}
CONFIG = resolve_config("fast", threads=2)
OUTPUT = ROOT / "docs/research/unified_domain_20261119/benchmark_cr0_high_density_extension_fast.json"


def solve(extra: dict[str, float]) -> dict[str, object]:
    params = BASE_PARAMS | extra
    model = LCDM_SG(**params)
    started = time.perf_counter()
    result = SGWB_iter_fast(model, config=CONFIG, tol=1.0e-6, physical_guard=True)
    elapsed = time.perf_counter() - started
    row: dict[str, object] = {
        "success": result is model,
        "seconds": elapsed,
        "failure_reason": model.fast_failure_reason,
        "params": params,
    }
    if result is model:
        freqs = np.asarray(model.f, dtype=np.float64)
        spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
        row.update(
            {
                "root": float(model.DN_gw[-1]),
                "n_freq": int(freqs.size),
                "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
                "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(),
            }
        )
    return row


def main() -> None:
    cold: dict[str, dict[str, object]] = {}
    rows = []
    for name, extra in CASES.items():
        cold[name] = solve(extra)
        warmups = [solve(extra) for _ in range(2)]
        samples = [solve(extra) for _ in range(REPEATS)]
        if not all(sample["success"] for sample in samples):
            raise RuntimeError(f"candidate failed during persistent-process timing: {name}")
        roots = {sample["root"] for sample in samples}
        frequency_hashes = {sample["frequency_sha256"] for sample in samples}
        spectrum_hashes = {sample["spectrum_sha256"] for sample in samples}
        times = [float(sample["seconds"]) for sample in samples]
        row = {
            "case": name,
            "params": BASE_PARAMS | extra,
            "pairs": REPEATS,
            "p50_seconds": float(np.percentile(times, 50)),
            "p95_seconds": float(np.percentile(times, 95)),
            "root_exact_across_repeats": len(roots) == 1,
            "frequency_grid_exact_across_repeats": len(frequency_hashes) == 1,
            "spectrum_exact_across_repeats": len(spectrum_hashes) == 1,
            "cold_first_solve": cold[name],
            "warmup_solvers": warmups,
            "seconds_samples": times,
        }
        rows.append(row)
        print(name, row["p50_seconds"], row["p95_seconds"], flush=True)

    output = {
        "status": "complete_persistent_process_fast_timing",
        "baseline_screen_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "candidate_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "candidate_source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("exact_background.py", "fast_sgwb.py")
        },
        "resources": {
            "affinity": affinity,
            "threads": 2,
            "threading_layer": "workqueue",
            "blas_omp_numexpr_threads": 1,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "warm_repeats": REPEATS,
            "frequency_grid": "goal",
            "frequency_quadrature": "pchip",
            "physical_guard": True,
            "config_h": float(CONFIG.h),
            "config_z_tail": float(CONFIG.z_tail),
        },
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
