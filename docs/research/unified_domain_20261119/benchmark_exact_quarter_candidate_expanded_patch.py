"""Warm paired performance benchmark at three expanded CR1 density sums."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import tarfile
import tempfile
import time
import warnings
from pathlib import Path

for key, value in {
    "NUMBA_NUM_THREADS": "2", "NUMBA_THREADING_LAYER": "workqueue",
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ.setdefault(key, value)

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).with_name("benchmark_exact_quarter_candidate_expanded_patch.json")
BASE_SHA = "48bbb96a44d4ced3f7e3d80dd2a2deaae440a0a3"
REPEATS = 128
CASES = {
    "sum_1300k": {"Omega_bh2": 650_000.0, "Omega_ch2": 650_000.0},
    "sum_1350k": {"Omega_bh2": 675_000.0, "Omega_ch2": 675_000.0},
    "sum_1400k": {"Omega_bh2": 700_000.0, "Omega_ch2": 700_000.0},
}


def _worker(root: Path) -> None:
    sys.path.insert(0, str(root))
    import numpy as np

    from stiffgwpy_fast.fast_sgwb import SGWB_iter_fast, resolve_config
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    print(json.dumps({"ready": True}), flush=True)
    for line in sys.stdin:
        request = json.loads(line)
        if request.get("stop"):
            return
        name, params = request["case"], request["params"]
        model = LCDM_SG(**params)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            start = time.perf_counter()
            result = SGWB_iter_fast(model, config=resolve_config("fast", threads=2),
                                    tol=1e-6, physical_guard=True)
            seconds = time.perf_counter() - start
        if result is not model:
            response = {"case": name, "success": False,
                        "failure_reason": model.fast_failure_reason, "seconds": seconds}
        else:
            spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
            freqs = np.asarray(model.f, dtype=np.float64)
            response = {"case": name, "success": True, "seconds": seconds,
                "root": float(model.DN_gw[-1]), "nfreq": int(freqs.size),
                "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(),
                "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
                "warning_count": len(caught)}
        print(json.dumps(response, allow_nan=False), flush=True)


def _spawn(root: Path, env: dict[str, str]):
    process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
        "--worker", str(root)], cwd=root, env=env, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1)
    assert process.stdout is not None
    if not json.loads(process.stdout.readline()).get("ready"):
        raise RuntimeError("benchmark worker did not become ready")
    return process


def _solve(process, name: str, params: dict[str, object]) -> dict[str, object]:
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps({"case": name, "params": params}) + "\n")
    process.stdin.flush()
    return json.loads(process.stdout.readline())


def _percentile(values: list[float], quantile: int) -> float:
    import numpy as np
    return float(np.percentile(np.asarray(values, dtype=np.float64), quantile))


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--worker":
        _worker(Path(sys.argv[2]))
        return
    env = os.environ.copy()
    archive = subprocess.run(["git", "archive", BASE_SHA], cwd=ROOT,
        check=True, stdout=subprocess.PIPE).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_expanded_patch_ab_") as tmp:
        base_root = Path(tmp)
        with tarfile.open(fileobj=__import__("io").BytesIO(archive), mode="r:") as tar:
            tar.extractall(base_root, filter="data")
        start = time.perf_counter()
        baseline = _spawn(base_root, env)
        base_start = time.perf_counter() - start
        start = time.perf_counter()
        candidate = _spawn(ROOT, env)
        candidate_start = time.perf_counter() - start
        rows = []
        cold = []
        try:
            for index, (name, density) in enumerate(CASES.items()):
                params = {"H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9,
                    "r": 0.01, "n_t": 0.0, "cr": 1.0, "T_re": 2000.0,
                    "DN_re": 10.0, "kappa10": 0.01, **density}
                order = (baseline, candidate) if index % 2 == 0 else (candidate, baseline)
                first = _solve(order[0], name, params)
                second = _solve(order[1], name, params)
                old, new = (first, second) if order[0] is baseline else (second, first)
                cold.append({"case": name, "baseline": old, "candidate": new})
                for _ in range(2):
                    _solve(baseline, name, params)
                    _solve(candidate, name, params)
                old_times, new_times = [], []
                output = None
                for repeat in range(REPEATS):
                    order = (baseline, candidate) if repeat % 2 == 0 else (candidate, baseline)
                    first = _solve(order[0], name, params)
                    second = _solve(order[1], name, params)
                    old, new = (first, second) if order[0] is baseline else (second, first)
                    old_times.append(float(old["seconds"]))
                    new_times.append(float(new["seconds"]))
                    output = {"baseline_success": old.get("success"),
                        "candidate_success": new.get("success"),
                        "baseline_root": old.get("root"), "candidate_root": new.get("root"),
                        "same_spectrum": old.get("spectrum_sha256") == new.get("spectrum_sha256"),
                        "same_grid": old.get("frequency_sha256") == new.get("frequency_sha256")}
                rows.append({"case": name, "params": params, "pairs": REPEATS,
                    "baseline_p50_s": _percentile(old_times, 50),
                    "baseline_p95_s": _percentile(old_times, 95),
                    "candidate_p50_s": _percentile(new_times, 50),
                    "candidate_p95_s": _percentile(new_times, 95),
                    "p50_ratio": _percentile(new_times, 50)/_percentile(old_times, 50),
                    "p95_ratio": _percentile(new_times, 95)/_percentile(old_times, 95),
                    "output_comparison": output,
                    "baseline_seconds_samples": old_times,
                    "candidate_seconds_samples": new_times})
                OUTPUT.write_text(json.dumps({"status":"in_progress", "base_sha":BASE_SHA,
                    "candidate_source_sha256": hashlib.sha256((ROOT/"stiffgwpy_fast/fast_sgwb.py").read_bytes()).hexdigest(),
                    "exact_background_sha256": hashlib.sha256((ROOT/"stiffgwpy_fast/exact_background.py").read_bytes()).hexdigest(),
                    "paired_cases":rows, "cold_first_solve":cold,
                    "worker_startup_seconds":{"baseline":base_start,"candidate":candidate_start},
                    "resources":{"threads":2,"threading_layer":"workqueue","repeat_pairs":REPEATS,
                        "python":sys.version.split()[0],"platform":platform.platform()}},
                    indent=2)+"\n",encoding="utf-8")
                print(json.dumps(rows[-1]), flush=True)
            old_all = [t for row in rows for t in row["baseline_seconds_samples"]]
            new_all = [t for row in rows for t in row["candidate_seconds_samples"]]
            payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
            payload["status"] = "complete"
            payload["aggregate_p50_ratio"] = _percentile(new_all,50)/_percentile(old_all,50)
            payload["aggregate_p95_ratio"] = _percentile(new_all,95)/_percentile(old_all,95)
            OUTPUT.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
            print(json.dumps({"status":"complete", "aggregate_p50_ratio":payload["aggregate_p50_ratio"],
                "aggregate_p95_ratio":payload["aggregate_p95_ratio"],"artifact":str(OUTPUT)}))
        finally:
            for process in (baseline,candidate):
                try:
                    if process.stdin is not None:
                        process.stdin.write(json.dumps({"stop":True})+"\n")
                        process.stdin.flush()
                    process.wait(timeout=10)
                except Exception:
                    process.kill()


if __name__ == "__main__":
    main()
