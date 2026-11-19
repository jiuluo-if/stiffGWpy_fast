"""同硬件比较基线提交与候选背景/交接实现的完整 Fast 成本。"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

for key, value in {
    "NUMBA_NUM_THREADS": "2",
    "NUMBA_THREADING_LAYER": "workqueue",
    "FAST_THREADS": "2",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ.setdefault(key, value)

HERE = Path(__file__).resolve().parent
CANDIDATE_ROOT = HERE.parents[2]
OUTPUT = HERE / "benchmark_h2_candidate_ab.json"
REPEATS = 64
CASES = {
    "trusted_default_cr1": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "trusted_stiff_cr0": dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    "trusted_highT_cr0": dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
    "trusted_highkappa_cr1": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
}


def _worker(root: Path, case: str) -> None:
    sys.path.insert(0, str(root))
    import numpy as np

    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    params = CASES[case]

    def solve() -> tuple[dict[str, object], float]:
        model = LCDM_SG(**params)
        with warnings.catch_warnings(record=True) as caught:
            started = time.perf_counter()
            result = model.SGWB_iter(
                engine="fast",
                fallback=False,
                tol=1e-6,
                threads=2,
                accuracy_mode="fast",
                physical_guard=False,
            )
            elapsed = time.perf_counter() - started
        if result is None:
            raise RuntimeError(f"{case} failed: {model.fast_failure_reason}")
        return (
            {
                "root": float(model.DN_gw[-1]),
                "nfreq": int(model.freq_grid_n),
                "spectrum": np.asarray(model.log10OmegaGW, dtype=np.float64),
                "warnings": [str(item.message) for item in caught],
            },
            elapsed,
        )

    first, first_s = solve()
    # 打印并立即 flush，让父进程在初始化后知道 worker 已就绪。
    print(json.dumps({"ready": True, "first_s": first_s}), flush=True)
    for _ in range(2):
        solve()
    for line in sys.stdin:
        if line.strip() != "run":
            continue
        output, elapsed = solve()
        output["spectrum"] = output["spectrum"].tolist()
        print(json.dumps({"seconds": elapsed, **output}, allow_nan=False), flush=True)


def _percentile(values: list[float], q: float) -> float:
    import numpy as np

    return float(np.percentile(np.asarray(values, dtype=np.float64), q))


def main() -> None:
    if len(sys.argv) == 4 and sys.argv[1] == "--worker":
        _worker(Path(sys.argv[2]), sys.argv[3])
        return
    if len(sys.argv) != 2:
        raise SystemExit("usage: benchmark_h2_candidate_ab.py BASELINE_ARCHIVE")

    baseline_root = Path(sys.argv[1]).resolve()
    rows = []
    cold = {}
    for case in CASES:
        workers = {}
        ready = {}
        for label, root in (("baseline", baseline_root), ("candidate", CANDIDATE_ROOT)):
            started = time.perf_counter()
            process = subprocess.Popen(
                [sys.executable, str(Path(__file__).resolve()), "--worker", str(root), case],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                env=os.environ.copy(),
            )
            assert process.stdout is not None
            ready[label] = json.loads(process.stdout.readline())
            ready[label]["startup_to_first_solve_wall_s"] = time.perf_counter() - started
            workers[label] = process
        cold[case] = ready
        times = {"baseline": [], "candidate": []}
        exact = True
        for repeat in range(REPEATS):
            order = ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
            paired = {}
            for label in order:
                process = workers[label]
                assert process.stdin is not None and process.stdout is not None
                process.stdin.write("run\n")
                process.stdin.flush()
                paired[label] = json.loads(process.stdout.readline())
                times[label].append(float(paired[label]["seconds"]))
            left, right = paired["baseline"], paired["candidate"]
            exact = exact and (
                left["root"] == right["root"]
                and left["nfreq"] == right["nfreq"]
                and left["spectrum"] == right["spectrum"]
                and left["warnings"] == right["warnings"]
            )
        for process in workers.values():
            process.terminate()
            process.wait(timeout=10)
        row = {
            "case": case,
            "params": CASES[case],
            "pairs": REPEATS,
            "baseline_p50_s": _percentile(times["baseline"], 50),
            "baseline_p95_s": _percentile(times["baseline"], 95),
            "candidate_p50_s": _percentile(times["candidate"], 50),
            "candidate_p95_s": _percentile(times["candidate"], 95),
            "p50_ratio": _percentile(times["candidate"], 50) / _percentile(times["baseline"], 50),
            "p95_ratio": _percentile(times["candidate"], 95) / _percentile(times["baseline"], 95),
            "outputs_exact": exact,
            "baseline_samples_s": times["baseline"],
            "candidate_samples_s": times["candidate"],
        }
        rows.append(row)
        print(json.dumps(row, allow_nan=False), flush=True)

    baseline_all = [x for row in rows for x in row["baseline_samples_s"]]
    candidate_all = [x for row in rows for x in row["candidate_samples_s"]]
    payload = {
        "status": "complete",
        "baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "resources": {
            "threads": 2,
            "threading_layer": os.environ.get("NUMBA_THREADING_LAYER"),
            "blas_omp_numexpr_threads": 1,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "warm_pairs_per_case": REPEATS,
            "paired_order": "alternating persistent workers",
        },
        "aggregate_p50_ratio": _percentile(candidate_all, 50) / _percentile(baseline_all, 50),
        "aggregate_p95_ratio": _percentile(candidate_all, 95) / _percentile(baseline_all, 95),
        "cold_first_solve": cold,
        "cases": rows,
        "source_sha256": {
            rel: hashlib.sha256((CANDIDATE_ROOT / rel).read_bytes()).hexdigest()
            for rel in ("stiffgwpy_fast/exact_background.py", "stiffgwpy_fast/fast_sgwb.py")
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
