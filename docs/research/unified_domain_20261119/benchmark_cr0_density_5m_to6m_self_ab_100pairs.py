"""Persistent-process full-solver A/B with isolated old/current source trees."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE_SHA = "59a23aec2300f08b50124d7f35dec4b228388d37"
REPEATS = 100
WARMUPS = 3
THREADING_LAYER = os.environ.get("STIFFGW_BENCH_LAYER", "tbb")
OUTPUT = HERE / f"benchmark_cr0_density_5m_to6m_self_ab_{THREADING_LAYER}_100pairs.json"


def _digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def _worker(root: Path) -> None:
    sys.path.insert(0, str(root))
    import numpy as np
    import psutil

    from stiffgwpy_fast import fast_sgwb
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    proc = psutil.Process()
    proc.cpu_affinity(proc.cpu_affinity()[:2])
    config = fast_sgwb.resolve_config("fast", threads=2)
    for line in sys.stdin:
        request = json.loads(line)
        if request.get("shutdown"):
            return
        model = LCDM_SG(**request["params"])
        started = time.perf_counter()
        result = fast_sgwb.SGWB_iter_fast(
            model, config=config, tol=1e-6, physical_guard=True
        )
        elapsed = time.perf_counter() - started
        if result is not model:
            row = {
                "success": False,
                "reason": model.fast_failure_reason,
                "seconds": elapsed,
            }
        else:
            frequencies = np.asarray(model.f, dtype=np.float64)
            spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
            row = {
                "success": True,
                "seconds": elapsed,
                "root": float(model.DN_gw[-1]),
                "frequency_sha256": _digest(frequencies),
                "spectrum_sha256": _digest(spectrum),
                "spectrum": spectrum.tolist(),
            }
        print(json.dumps(row), flush=True)


def _start_worker(root: Path, cache: Path):
    cache.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({
        "NUMBA_NUM_THREADS": "2",
        "NUMBA_THREADING_LAYER": THREADING_LAYER,
        "FAST_THREADS": "2",
        "NUMBA_CACHE_DIR": str(cache),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--worker", str(root)],
        cwd=root,
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    return process


def _request(process, params):
    if process.stdin is None or process.stdout is None:
        raise RuntimeError("solver worker pipes are unavailable")
    process.stdin.write(json.dumps({"params": params}) + "\n")
    process.stdin.flush()
    response = process.stdout.readline()
    if not response:
        stderr = "" if process.stderr is None else process.stderr.read()
        raise RuntimeError(f"solver worker exited {process.poll()}: {stderr}")
    return json.loads(response)


def _stop(process) -> None:
    if process.stdin is not None and process.poll() is None:
        process.stdin.write('{"shutdown":true}\n')
        process.stdin.flush()
    process.wait(timeout=10)


def _cases() -> dict[str, dict[str, float]]:
    lattice = json.loads(
        (HERE / "screen_cr0_density_5m_to6m.json").read_text(encoding="utf-8")
    )
    sobol = json.loads(
        (HERE / "screen_cr0_density_5m_to6m_sobol.json").read_text(encoding="utf-8")
    )
    upper_corner = next(
        row for row in lattice["cases"] if row["input"]["case_id"] == "rectangle_06"
    )["input"]["params"]
    independent_holdout = next(
        row for row in sobol["cases"] if row["input"]["case_id"] == "holdout_03"
    )["input"]["params"]
    return {
        "new_region_upper_density_low_tilt": upper_corner,
        "new_region_independent_sobol": independent_holdout,
        "trusted_physical_cr0_control": {
            "Omega_bh2": 0.022,
            "Omega_ch2": 0.12,
            "H0": 67.32117,
            "DN_eff": 0.0,
            "A_s": 2.100549e-9,
            "r": 0.01,
            "n_t": 0.0,
            "cr": 0.0,
            "T_re": 2000.0,
            "DN_re": 10.0,
            "kappa10": 0.01,
        },
    }


def main() -> None:
    affinity = []
    try:
        import psutil

        affinity = psutil.Process().cpu_affinity()[:2]
        psutil.Process().cpu_affinity(affinity)
    except (ImportError, OSError, AttributeError):
        pass
    payload_cases = []
    with tempfile.TemporaryDirectory(prefix="stiffgwpy_persistent_full_ab_") as tmp:
        temp_root = Path(tmp)
        baseline_root = temp_root / "baseline"
        baseline_root.mkdir()
        archive = subprocess.run(
            ["git", "archive", BASE_SHA], cwd=ROOT, check=True,
            stdout=subprocess.PIPE,
        ).stdout
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
            bundle.extractall(baseline_root, filter="data")
        candidate_root = temp_root / "candidate"
        candidate_root.mkdir()
        shutil.copytree(ROOT / "stiffgwpy_fast", candidate_root / "stiffgwpy_fast")
        base_worker = _start_worker(baseline_root, temp_root / "numba-base")
        current_worker = _start_worker(candidate_root, temp_root / "numba-current")
        try:
            for case_id, params in _cases().items():
                first_base = _request(base_worker, params)
                first_current = _request(current_worker, params)
                if not first_base["success"] or not first_current["success"]:
                    raise RuntimeError(
                        f"cold solve failed for {case_id}: "
                        f"base={first_base.get('reason')}, current={first_current.get('reason')}"
                    )
                for _ in range(WARMUPS):
                    _request(base_worker, params)
                    _request(current_worker, params)
                base_times = []
                current_times = []
                root_deltas = []
                frequency_mismatches = 0
                spectrum_max = 0.0
                pairs = []
                for index in range(REPEATS):
                    order = ("base", "current") if index % 2 == 0 else ("current", "base")
                    first = _request(
                        base_worker if order[0] == "base" else current_worker, params
                    )
                    second = _request(
                        base_worker if order[1] == "base" else current_worker, params
                    )
                    base, current = (first, second) if order[0] == "base" else (second, first)
                    if not base["success"] or not current["success"]:
                        raise RuntimeError(f"warm paired solve failed for {case_id}")
                    base_times.append(base["seconds"])
                    current_times.append(current["seconds"])
                    frequency_equal = base["frequency_sha256"] == current["frequency_sha256"]
                    frequency_mismatches += int(not frequency_equal)
                    if len(base["spectrum"]) != len(current["spectrum"]):
                        spectrum_delta = float("inf")
                    else:
                        spectrum_delta = float(max(
                            abs(x - y) for x, y in zip(base["spectrum"], current["spectrum"])
                        ))
                    spectrum_max = max(spectrum_max, spectrum_delta)
                    root_delta = abs(current["root"] - base["root"]) / max(abs(base["root"]), 1e-300)
                    root_deltas.append(root_delta)
                    pairs.append({
                        "index": index,
                        "first_side": order[0],
                        "base_seconds": base["seconds"],
                        "current_seconds": current["seconds"],
                        "base_root": base["root"],
                        "current_root": current["root"],
                        "root_relative_delta": root_delta,
                        "frequency_digest_equal": frequency_equal,
                        "spectrum_max_abs_delta_dex": spectrum_delta,
                    })
                base_times = np.asarray(base_times)
                current_times = np.asarray(current_times)
                paired = current_times / base_times
                row = {
                    "case_id": case_id,
                    "params": params,
                    "pairs": REPEATS,
                    "first_call": {
                        "fresh_process_compile": not payload_cases,
                        "base_seconds": first_base["seconds"],
                        "current_seconds": first_current["seconds"],
                        "base_root": first_base["root"],
                        "current_root": first_current["root"],
                    },
                    "base_p50_seconds": float(np.median(base_times)),
                    "current_p50_seconds": float(np.median(current_times)),
                    "p50_ratio": float(np.median(current_times) / np.median(base_times)),
                    "base_p95_seconds": float(np.quantile(base_times, 0.95)),
                    "current_p95_seconds": float(np.quantile(current_times, 0.95)),
                    "p95_ratio": float(np.quantile(current_times, 0.95) / np.quantile(base_times, 0.95)),
                    "paired_ratio_p95": float(np.quantile(paired, 0.95)),
                    "root_relative_delta_max": float(max(root_deltas)),
                    "frequency_digest_mismatch_pairs": frequency_mismatches,
                    "spectrum_absdiff_max_dex": spectrum_max,
                    "timing_pairs": pairs,
                }
                payload_cases.append(row)
                print(case_id, row["p50_ratio"], row["p95_ratio"], spectrum_max, flush=True)
        finally:
            _stop(base_worker)
            _stop(current_worker)
    payload = {
        "status": "complete_persistent_process_full_fast_pinned_baseline_ab",
        "base_commit": BASE_SHA,
        "candidate_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "resources": {
            "affinity": affinity,
            "numba_threads": 2,
            "threading_layer": THREADING_LAYER,
            "blas_threads": 1,
            "persistent_processes": 2,
            "alternating_pairs_per_case": REPEATS,
            "warmups_per_side": WARMUPS,
            "separate_numba_caches": True,
            "physical_guard": True,
            "frequency_grid": "goal",
        },
        "cases": payload_cases,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("worker_root", nargs="?")
    args = parser.parse_args()
    if args.worker:
        if args.worker_root is None:
            raise SystemExit("--worker requires a source root")
        _worker(Path(args.worker_root).resolve())
    else:
        main()
