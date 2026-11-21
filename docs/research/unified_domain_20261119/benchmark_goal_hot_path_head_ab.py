"""Compare unchanged default-goal runs on HEAD and candidate source."""
from __future__ import annotations

import io
import hashlib
import json
import statistics
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).resolve().parent / "benchmark_goal_hot_path_head_ab.json"
CASES = {
    "trusted_default": {
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
        "r": 0.01,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 2000.0,
        "DN_re": 10.0,
        "kappa10": 0.01,
    },
    "trusted_highT": {
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
        "r": 0.01,
        "n_t": 0.0,
        "cr": 1.0,
        "T_re": 10000.0,
        "DN_re": 10.0,
        "kappa10": 0.01,
    },
    "cr0_density_edge": {
        "Omega_bh2": 13_890_136.71875,
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
        "r": 0.01,
        "n_t": -0.022705078125,
        "cr": 0.0,
        "T_re": 2000.0,
        "DN_re": 10.0,
        "kappa10": 0.01,
    },
}
CHILD = r'''
import json,os,time
from pathlib import Path
from stiffgwpy_fast import LCDM_SG
from stiffgwpy_fast import fast_sgwb as FS
assert Path(FS.__file__).resolve().is_relative_to(Path.cwd().resolve())
cases=json.loads(os.environ['STIFFGW_BENCH_CASES'])
out={}
for name,params in cases.items():
 times=[];roots=[]
 for _ in range(4):
  m=LCDM_SG(**params);cfg=FS.resolve_config('fast',threads=2,freq_grid='goal')
  started=time.perf_counter();r=FS.SGWB_iter_fast(m,config=cfg,tol=1e-6,physical_guard=True);times.append(time.perf_counter()-started)
  if r is not m:raise RuntimeError((name,getattr(m,'fast_failure_reason',None)))
  roots.append(float(m.DN_gw[-1]))
 out[name]={'seconds':times,'roots':roots}
print(json.dumps(out,sort_keys=True))
'''


def run(cwd: str) -> dict:
    import os

    env = dict(os.environ)
    env["STIFFGW_BENCH_CASES"] = json.dumps(CASES)
    result = subprocess.run(
        [sys.executable, "-c", CHILD],
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    return json.loads(result.stdout.strip().splitlines()[-1])


def main() -> None:
    archive = subprocess.check_output(
        ["git", "archive", "--format=tar", "HEAD"], cwd=ROOT
    )
    with tempfile.TemporaryDirectory(prefix="stiffgwpy-head-goal-ab-") as baseline_dir:
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as source:
            source.extractall(baseline_dir, filter="data")
        samples = {
            variant: {name: {"seconds": [], "roots": []} for name in CASES}
            for variant in ("head", "candidate")
        }
        rounds = 8
        for round_index in range(rounds):
            order = ("head", "candidate") if round_index % 2 == 0 else ("candidate", "head")
            for variant in order:
                cwd = baseline_dir if variant == "head" else str(ROOT)
                result = run(cwd)
                for name in CASES:
                    samples[variant][name]["seconds"].extend(result[name]["seconds"][1:])
                    samples[variant][name]["roots"].extend(result[name]["roots"][1:])

    rows = {}
    for name in CASES:
        base_times = samples["head"][name]["seconds"]
        candidate_times = samples["candidate"][name]["seconds"]
        base_roots = samples["head"][name]["roots"]
        candidate_roots = samples["candidate"][name]["roots"]
        rows[name] = {
            "baseline_warm_seconds": base_times,
            "candidate_warm_seconds": candidate_times,
            "baseline_roots": base_roots,
            "candidate_roots": candidate_roots,
            "baseline_warm_p50": statistics.median(base_times),
            "candidate_warm_p50": statistics.median(candidate_times),
            "baseline_warm_p95_sample_max": max(base_times),
            "candidate_warm_p95_sample_max": max(candidate_times),
            "warm_p50_ratio": statistics.median(candidate_times)
            / statistics.median(base_times),
            "warm_p95_ratio": max(candidate_times) / max(base_times),
            "roots_identical": base_roots == candidate_roots,
        }
    payload = {
        "status": "complete",
        "hardware_threads": 2,
        "paired_process_rounds": 8,
        "warm_sample_count_per_case": 24,
        "cold_note": "Each process discards its first per-case solve; source variants are alternated across process rounds and the extracted HEAD cache is retained between rounds.",
        "baseline_fast_sgwb_sha256": hashlib.sha256(
            subprocess.check_output(
                ["git", "show", "HEAD:stiffgwpy_fast/fast_sgwb.py"], cwd=ROOT
            )
        ).hexdigest(),
        "candidate_fast_sgwb_sha256": hashlib.sha256(
            (ROOT / "stiffgwpy_fast" / "fast_sgwb.py").read_bytes()
        ).hexdigest(),
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
