"""Paired trusted-region benchmark for the exact-quarter transfer candidate."""

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

import psutil

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
OUTPUT = Path(__file__).resolve().with_name("benchmark_exact_quarter_candidate_trusted_ab_affinity.json")
BASE_SHA = "48bbb96a44d4ced3f7e3d80dd2a2deaae440a0a3"
REPEATS = 128
CASES = {
    "default_cr1": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "lowT_cr1": dict(r=1e-2, cr=1, T_re=10.0, kappa10=1e-2),
    "highT_cr1": dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "highkappa_cr1": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    "stiff_cr0": dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    "highT_cr0": dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
}


def _worker(root: Path) -> None:
    sys.path.insert(0, str(root))
    import numpy as np

    from stiffgwpy_fast.fast_sgwb import (  # noqa: E402
        SGWB_iter_fast,
        resolve_config,
    )
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

    print(json.dumps({"ready": True, "affinity": psutil.Process().cpu_affinity()}), flush=True)
    for line in sys.stdin:
        request = json.loads(line)
        if request.get("stop"):
            return
        name = request["case"]
        params = request["params"]
        model = LCDM_SG(**params)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            start = time.perf_counter()
            result = SGWB_iter_fast(
                model,
                config=resolve_config("fast", threads=2),
                tol=1e-6,
                kink_split=True,
                physical_guard=False,
            )
            seconds = time.perf_counter() - start
        if result is not model:
            response = {
                "case": name,
                "success": False,
                "failure_reason": model.fast_failure_reason,
                "seconds": seconds,
            }
        else:
            spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
            freqs = np.asarray(model.f, dtype=np.float64)
            response = {
                "case": name,
                "success": True,
                "seconds": seconds,
                "root": float(model.DN_gw[-1]),
                "nfreq": int(freqs.size),
                "spectrum_sha256": hashlib.sha256(spectrum.tobytes()).hexdigest(),
                "frequency_sha256": hashlib.sha256(freqs.tobytes()).hexdigest(),
                "warnings": [str(item.message) for item in caught],
                "spectrum": spectrum.tolist(),
            }
        print(json.dumps(response, allow_nan=False), flush=True)


def _spawn(root: Path, env: dict[str, str]):
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--worker", str(root)],
        cwd=root,
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    ready = json.loads(process.stdout.readline())
    if not ready.get("ready") or ready.get("affinity") != [0, 1]:
        raise RuntimeError("benchmark worker failed to initialize")
    return process


def _solve(process, name: str, params: dict[str, float]) -> dict[str, object]:
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps({"case": name, "params": params}) + "\n")
    process.stdin.flush()
    return json.loads(process.stdout.readline())


def _percentile(values: list[float], q: float) -> float:
    import numpy as np

    return float(np.percentile(np.asarray(values, dtype=np.float64), q))


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--worker":
        _worker(Path(sys.argv[2]))
        return

    process = psutil.Process()
    available_cpus = process.cpu_affinity()
    if len(available_cpus) < 2:
        raise RuntimeError("canonical benchmark requires at least two available CPUs")
    affinity = available_cpus[:2]
    process.cpu_affinity(affinity)
    if process.cpu_affinity() != [0, 1]:
        raise RuntimeError(f"canonical affinity [0, 1] unavailable: {process.cpu_affinity()}")

    env = os.environ.copy()
    archive = subprocess.run(
        ["git", "archive", BASE_SHA],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_exact_quarter_base_") as tmp:
        base_root = Path(tmp)
        import io

        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
            tar.extractall(base_root, filter="data")

        started = time.perf_counter()
        base = _spawn(base_root, env)
        base_start = time.perf_counter() - started
        started = time.perf_counter()
        candidate = _spawn(ROOT, env)
        candidate_start = time.perf_counter() - started
        try:
            cold_rows = []
            for name, extra in CASES.items():
                params = {"r": 0.01, "cr": 0, "T_re": 2e3, "kappa10": 0.01, **extra}
                order = (base, candidate) if len(cold_rows) % 2 == 0 else (candidate, base)
                first = _solve(order[0], name, params)
                second = _solve(order[1], name, params)
                pair = (
                    {"baseline": first, "candidate": second}
                    if order[0] is base
                    else {"baseline": second, "candidate": first}
                )
                for row in pair.values():
                    row.pop("spectrum", None)
                cold_rows.append({"case": name, "paired_first_solve": pair})

            rows = []
            for name, extra in CASES.items():
                params = {"r": 0.01, "cr": 0, "T_re": 2e3, "kappa10": 0.01, **extra}
                for _ in range(2):
                    _solve(base, name, params)
                    _solve(candidate, name, params)
                base_times = []
                candidate_times = []
                outputs_exact = True
                candidate_meta = None
                for repeat in range(REPEATS):
                    order = (base, candidate) if repeat % 2 == 0 else (candidate, base)
                    first = _solve(order[0], name, params)
                    second = _solve(order[1], name, params)
                    old, new = (first, second) if order[0] is base else (second, first)
                    base_times.append(float(old["seconds"]))
                    candidate_times.append(float(new["seconds"]))
                    outputs_exact = outputs_exact and (
                        old.get("success") == new.get("success")
                        and old.get("root") == new.get("root")
                        and old.get("nfreq") == new.get("nfreq")
                        and old.get("frequency_sha256") == new.get("frequency_sha256")
                        and old.get("spectrum_sha256") == new.get("spectrum_sha256")
                        and old.get("warnings") == new.get("warnings")
                    )
                    candidate_meta = {
                        "root": new.get("root"),
                        "nfreq": new.get("nfreq"),
                        "frequency_sha256": new.get("frequency_sha256"),
                        "spectrum_sha256": new.get("spectrum_sha256"),
                    }
                rows.append(
                    {
                        "case": name,
                        "params": params,
                        "pairs": REPEATS,
                        "baseline_p50_s": _percentile(base_times, 50),
                        "baseline_p95_s": _percentile(base_times, 95),
                        "candidate_p50_s": _percentile(candidate_times, 50),
                        "candidate_p95_s": _percentile(candidate_times, 95),
                        "p50_ratio": _percentile(candidate_times, 50) / _percentile(base_times, 50),
                        "p95_ratio": _percentile(candidate_times, 95) / _percentile(base_times, 95),
                        "outputs_exact": outputs_exact,
                        "candidate_output": candidate_meta,
                        "baseline_seconds_samples": base_times,
                        "candidate_seconds_samples": candidate_times,
                    }
                )
                payload = {
                    "status": "in_progress",
                    "base_sha": BASE_SHA,
                    "candidate_source_sha256": hashlib.sha256(
                        (ROOT / "stiffgwpy_fast/fast_sgwb.py").read_bytes()
                    ).hexdigest(),
                    "paired_cases": rows,
                    "cold_first_solve": cold_rows,
                    "worker_startup_seconds": {"baseline": base_start, "candidate": candidate_start},
                    "resources": {
                        "affinity": affinity,
                        "threads": 2,
                        "threading_layer": "workqueue",
                        "blas_omp_numexpr_threads": 1,
                        "python": sys.version.split()[0],
                        "platform": platform.platform(),
                        "repeats_per_case": REPEATS,
                        "paired_order": "alternating",
                        "kink_split": True,
                    },
                }
                payload["resources"]["affinity"] = affinity
                payload["resources"]["paired_order"] = "alternating"
                payload["candidate_commit"] = subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
                ).strip()
                OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
                print(json.dumps(rows[-1]), flush=True)

            base_all = [value for row in rows for value in row["baseline_seconds_samples"]]
            candidate_all = [value for row in rows for value in row["candidate_seconds_samples"]]
            payload["status"] = "complete"
            payload["aggregate_p50_ratio"] = _percentile(candidate_all, 50) / _percentile(base_all, 50)
            payload["aggregate_p95_ratio"] = _percentile(candidate_all, 95) / _percentile(base_all, 95)
            OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            print(
                json.dumps(
                    {
                        "status": payload["status"],
                        "p50_ratio": payload["aggregate_p50_ratio"],
                        "p95_ratio": payload["aggregate_p95_ratio"],
                        "artifact": str(OUTPUT),
                    }
                )
            )
        finally:
            for process in (base, candidate):
                try:
                    if process.stdin is not None:
                        process.stdin.write(json.dumps({"stop": True}) + "\n")
                        process.stdin.flush()
                    process.wait(timeout=10)
                except Exception:
                    process.kill()


if __name__ == "__main__":
    main()
