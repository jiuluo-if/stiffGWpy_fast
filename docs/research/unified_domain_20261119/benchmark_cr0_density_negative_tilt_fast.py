"""Measure current Fast cost inside the CR0 density x tilt region."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

for _key, _value in {
    "NUMBA_NUM_THREADS": "2",
    "NUMBA_THREADING_LAYER": "workqueue",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ.setdefault(_key, _value)

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).with_name("benchmark_cr0_density_negative_tilt_fast.json")
REPEATS = 64
BASE_PARAMS = {
    "Omega_bh2": 1.0e5,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "n_t": -0.5,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = {
    "low_density_lower_tilt": {"Omega_bh2": 1.0e5, "n_t": -0.5},
    "high_density_lower_tilt": {"Omega_bh2": 1.0e6, "n_t": -0.5},
    "interior_midpoint": {"Omega_bh2": 5.0e5, "n_t": -0.25},
    "independent_holdout": {"Omega_bh2": 7.5e5, "n_t": -0.125},
}


def _worker() -> None:
    import numpy as np

    from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    print(json.dumps({"ready": True}), flush=True)
    for line in sys.stdin:
        request = json.loads(line)
        if request.get("stop"):
            return
        name = request["case"]
        params = dict(BASE_PARAMS, **request["overrides"])
        model = LCDM_SG(**params)
        started = time.perf_counter()
        result = SGWB_iter_fast(
            model,
            config=resolve_config("fast", threads=2),
            tol=1.0e-6,
            physical_guard=True,
        )
        seconds = time.perf_counter() - started
        response: dict[str, object] = {
            "case": name,
            "success": result is model,
            "seconds": seconds,
            "failure_reason": model.fast_failure_reason,
        }
        if result is model:
            spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
            freqs = np.asarray(model.f, dtype=np.float64)
            response.update({
                "root": float(model.DN_gw[-1]),
                "nfreq": int(freqs.size),
                "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(),
                "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
            })
        print(json.dumps(response, allow_nan=False), flush=True)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _solve(process: subprocess.Popen[str], case: str, overrides: dict[str, float]):
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps({"case": case, "overrides": overrides}) + "\n")
    process.stdin.flush()
    return json.loads(process.stdout.readline())


def main() -> None:
    if len(sys.argv) == 2 and sys.argv[1] == "--worker":
        _worker()
        return

    env = os.environ.copy()
    started = time.perf_counter()
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--worker"],
        cwd=ROOT,
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    if not json.loads(process.stdout.readline()).get("ready"):
        raise RuntimeError("Fast benchmark worker failed to initialize")
    startup_seconds = time.perf_counter() - started
    cold_rows = []
    warm_rows = []
    try:
        for index, (name, overrides) in enumerate(CASES.items()):
            first = _solve(process, name, overrides)
            if index == 0:
                cold_rows.append(first)
            for _ in range(2):
                _solve(process, name, overrides)
            samples = [_solve(process, name, overrides) for _ in range(REPEATS)]
            times = [float(item["seconds"]) for item in samples]
            roots = {item.get("root") for item in samples}
            spectra = {item.get("spectrum_sha256") for item in samples}
            grids = {item.get("frequency_sha256") for item in samples}
            if any(not item.get("success") for item in samples):
                raise RuntimeError(f"Fast failed during benchmark: {name}")
            row = {
                "case": name,
                "params": dict(BASE_PARAMS, **overrides),
                "pairs": REPEATS,
                "p50_seconds": float(np_percentile(times, 50)),
                "p95_seconds": float(np_percentile(times, 95)),
                "root_exact_across_repeats": len(roots) == 1,
                "spectrum_exact_across_repeats": len(spectra) == 1,
                "frequency_grid_exact_across_repeats": len(grids) == 1,
                "seconds_samples": times,
            }
            warm_rows.append(row)
            payload = _payload(startup_seconds, cold_rows, warm_rows)
            payload["status"] = "in_progress"
            OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({k: row[k] for k in row if k != "seconds_samples"}), flush=True)
        payload = _payload(startup_seconds, cold_rows, warm_rows)
        payload["status"] = "complete"
        OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"artifact: {OUTPUT}")
    finally:
        if process.stdin is not None:
            process.stdin.write(json.dumps({"stop": True}) + "\n")
            process.stdin.flush()
        process.wait(timeout=10)


def np_percentile(values: list[float], percentile: int) -> float:
    import numpy as np

    return float(np.percentile(np.asarray(values, dtype=np.float64), percentile))


def _payload(startup: float, cold: list[dict[str, object]], warm: list[dict[str, object]]):
    return {
        "status": "in_progress",
        "base_commit": "7c548be4757343137dc3c9e576ee778de8453376",
        "candidate_source_sha256": {
            name: _sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
        "resources": {
            "threads": 2,
            "threading_layer": "workqueue",
            "warm_repeats_per_case": REPEATS,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "worker_startup_seconds": startup,
            "all_full_fast_outputs_retained": True,
        },
        "cold_first_solve": cold,
        "warm_cases": warm,
    }


if __name__ == "__main__":
    main()
