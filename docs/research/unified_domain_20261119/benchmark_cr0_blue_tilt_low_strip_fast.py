"""Measure Fast latency within the CR0 positive-tilt subregion."""

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
    "NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue",
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ.setdefault(_key, _value)

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "benchmark_cr0_blue_tilt_low_strip_fast.json"
REPEATS = 64
FIXED = {
    "Omega_bh2": 1.0e6, "Omega_ch2": 0.1201075, "H0": 67.32117,
    "DN_eff": 0.0, "A_s": 2.100549e-9, "cr": 0.0,
    "T_re": 2000.0, "DN_re": 10.0, "kappa10": 0.01,
}
CASES = {
    "low_r_lower_tilt": {"r": 1.0e-12, "n_t": 0.05},
    "high_r_upper_tilt": {"r": 1.0e-8, "n_t": 0.1},
    "interior": {"r": 1.0e-10, "n_t": 0.075},
    "sobol_holdout": {"r": 3.96526332334898e-9, "n_t": 0.07029469762928785},
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
        params = request["params"]
        model = LCDM_SG(**params)
        started = time.perf_counter()
        solved = SGWB_iter_fast(
            model, config=resolve_config("fast", threads=2),
            tol=1.0e-6, physical_guard=True,
        )
        response: dict[str, object] = {
            "success": solved is model,
            "seconds": time.perf_counter() - started,
            "failure_reason": model.fast_failure_reason,
        }
        if solved is model:
            spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
            freqs = np.asarray(model.f, dtype=np.float64)
            response.update({
                "root": float(model.DN_gw[-1]), "nfreq": int(freqs.size),
                "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(),
                "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
            })
        print(json.dumps(response, allow_nan=False), flush=True)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _solve(process: subprocess.Popen[str], params: dict[str, float]):
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps({"params": params}) + "\n")
    process.stdin.flush()
    return json.loads(process.stdout.readline())


def _percentile(samples: list[float], percentile: int) -> float:
    import numpy as np

    return float(np.percentile(np.asarray(samples, dtype=np.float64), percentile))


def main() -> None:
    if len(sys.argv) == 2 and sys.argv[1] == "--worker":
        _worker()
        return

    started = time.perf_counter()
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--worker"],
        cwd=ROOT, env=os.environ.copy(), stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1,
    )
    assert process.stdout is not None
    if not json.loads(process.stdout.readline()).get("ready"):
        raise RuntimeError("Fast benchmark worker failed to initialize")
    startup = time.perf_counter() - started
    cold: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    try:
        for index, (name, overrides) in enumerate(CASES.items()):
            params = dict(FIXED, **overrides)
            first = _solve(process, params)
            if index == 0:
                cold.append({"case": name, **first})
            for _ in range(2):
                _solve(process, params)
            samples = [_solve(process, params) for _ in range(REPEATS)]
            if any(not row.get("success") for row in samples):
                raise RuntimeError(f"Fast failed during benchmark: {name}")
            row = {
                "case": name, "params": params, "repeats": REPEATS,
                "p50_seconds": _percentile([float(x["seconds"]) for x in samples], 50),
                "p95_seconds": _percentile([float(x["seconds"]) for x in samples], 95),
                "roots_exact": len({x["root"] for x in samples}) == 1,
                "spectra_exact": len({x["spectrum_sha256"] for x in samples}) == 1,
                "grids_exact": len({x["frequency_sha256"] for x in samples}) == 1,
                "seconds_samples": [float(x["seconds"]) for x in samples],
            }
            rows.append(row)
            _write(startup, cold, rows, "in_progress")
            print(json.dumps({k: v for k, v in row.items() if k != "seconds_samples"}), flush=True)
        _write(startup, cold, rows, "complete")
        print(f"artifact: {OUTPUT}")
    finally:
        if process.stdin is not None:
            process.stdin.write(json.dumps({"stop": True}) + "\n")
            process.stdin.flush()
        process.wait(timeout=10)


def _write(startup: float, cold: list[dict[str, object]], rows: list[dict[str, object]], status: str) -> None:
    payload = {
        "status": status,
        "candidate_head": "1a28dc08b237be12f8776c70d2983693cc1a441e",
        "candidate_source_sha256": {
            name: _sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("fast_sgwb.py", "exact_background.py", "reference.py")
        },
        "resources": {
            "threads": 2, "threading_layer": "workqueue",
            "warm_repeats_per_case": REPEATS, "python": sys.version.split()[0],
            "platform": platform.platform(), "worker_startup_seconds": startup,
            "full_spectrum_and_frequency_grid_returned": True,
        },
        "cold_first_solve": cold,
        "warm_cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
