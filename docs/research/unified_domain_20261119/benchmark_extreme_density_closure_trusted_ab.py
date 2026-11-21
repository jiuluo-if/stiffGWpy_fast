"""Run the established alternating trusted-region HEAD/candidate benchmark."""
from __future__ import annotations

import importlib.util
import io
import json
import os
import subprocess
import statistics
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(__file__).with_name("benchmark_goal_hot_path_head_ab.py")
spec = importlib.util.spec_from_file_location("extreme_density_trusted_ab", SOURCE)
if spec is None or spec.loader is None:
    raise RuntimeError("cannot load the paired Fast benchmark")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)
OUTPUT = Path(__file__).with_name("benchmark_extreme_density_closure_trusted_ab.json")
benchmark.CHILD = r'''
import hashlib,json,os,time
import numpy as np
from pathlib import Path
from stiffgwpy_fast import LCDM_SG
from stiffgwpy_fast import fast_sgwb as FS
assert Path(FS.__file__).resolve().is_relative_to(Path.cwd().resolve())
cases=json.loads(os.environ['STIFFGW_BENCH_CASES'])
out={}
for name,params in cases.items():
 times=[];roots=[];spectra=[];frequencies=[]
 for _ in range(4):
  m=LCDM_SG(**params);cfg=FS.resolve_config('fast',threads=2,freq_grid='goal')
  started=time.perf_counter();r=FS.SGWB_iter_fast(m,config=cfg,tol=1e-6,physical_guard=True);times.append(time.perf_counter()-started)
  if r is not m:raise RuntimeError((name,getattr(m,'fast_failure_reason',None)))
  roots.append(float(m.DN_gw[-1]))
  spectra.append(hashlib.sha256(np.asarray(m.log10OmegaGW,dtype=np.float64).tobytes()).hexdigest())
  frequencies.append(hashlib.sha256(np.asarray(m.f,dtype=np.float64).tobytes()).hexdigest())
 out[name]={'seconds':times,'roots':roots,'spectra':spectra,'frequencies':frequencies}
print(json.dumps(out,sort_keys=True))
'''


def _run_with_cache(cwd: str, cache_dir: str) -> dict:
    env = dict(os.environ)
    env["STIFFGW_BENCH_CASES"] = json.dumps(benchmark.CASES)
    env["NUMBA_CACHE_DIR"] = cache_dir
    result = subprocess.run(
        [sys.executable, "-c", benchmark.CHILD], cwd=cwd, env=env,
        text=True, capture_output=True, check=True, timeout=1200,
    )
    return json.loads(result.stdout.strip().splitlines()[-1])


def main() -> None:
    archive = subprocess.check_output(["git", "archive", "--format=tar", "HEAD"], cwd=ROOT)
    with tempfile.TemporaryDirectory(prefix="stiffgwpy-closure-head-") as baseline_dir, \
            tempfile.TemporaryDirectory(prefix="stiffgwpy-closure-cache-head-") as head_cache, \
            tempfile.TemporaryDirectory(prefix="stiffgwpy-closure-cache-candidate-") as candidate_cache:
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as source:
            source.extractall(baseline_dir, filter="data")
        cold_trials = {"head": [], "candidate": []}
        for _ in range(3):
            with tempfile.TemporaryDirectory(prefix="stiffgwpy-cold-head-") as cold_head_cache, \
                    tempfile.TemporaryDirectory(prefix="stiffgwpy-cold-candidate-") as cold_candidate_cache:
                old = _run_with_cache(baseline_dir, cold_head_cache)
                new = _run_with_cache(str(ROOT), cold_candidate_cache)
                cold_case = next(iter(benchmark.CASES))
                cold_trials["head"].append(old[cold_case]["seconds"][0])
                cold_trials["candidate"].append(new[cold_case]["seconds"][0])
        samples = {
            variant: {name: {"seconds": [], "roots": [], "spectra": [], "frequencies": []}
                      for name in benchmark.CASES}
            for variant in ("head", "candidate")
        }
        cold = {}
        rounds = 16
        for round_index in range(rounds):
            order = ("head", "candidate") if round_index % 2 == 0 else ("candidate", "head")
            for variant in order:
                cwd, cache = ((baseline_dir, head_cache) if variant == "head"
                              else (str(ROOT), candidate_cache))
                result = _run_with_cache(cwd, cache)
                for name in benchmark.CASES:
                    if round_index == 0:
                        cold.setdefault(variant, {})[name] = result[name]["seconds"][0]
                    samples[variant][name]["seconds"].extend(result[name]["seconds"][1:])
                    samples[variant][name]["roots"].extend(result[name]["roots"])
                    samples[variant][name]["spectra"].extend(result[name]["spectra"])
                    samples[variant][name]["frequencies"].extend(result[name]["frequencies"])

    rows = {}
    for name in benchmark.CASES:
        base = samples["head"][name]["seconds"]
        cand = samples["candidate"][name]["seconds"]
        base_roots = samples["head"][name]["roots"]
        candidate_roots = samples["candidate"][name]["roots"]
        base_spectra = samples["head"][name]["spectra"]
        candidate_spectra = samples["candidate"][name]["spectra"]
        base_frequencies = samples["head"][name]["frequencies"]
        candidate_frequencies = samples["candidate"][name]["frequencies"]
        rows[name] = {
            "head_first_solve_seconds": cold["head"][name],
            "candidate_first_solve_seconds": cold["candidate"][name],
            "head_warm_seconds": base,
            "candidate_warm_seconds": cand,
            "head_warm_p50": statistics.median(base),
            "candidate_warm_p50": statistics.median(cand),
            "head_warm_p95_sample_max": max(base),
            "candidate_warm_p95_sample_max": max(cand),
            "warm_p50_ratio": statistics.median(cand) / statistics.median(base),
            "warm_p95_ratio": float(np.quantile(cand, 0.95) / np.quantile(base, 0.95)),
            "warm_sample_max_ratio": max(cand) / max(base),
            "roots_identical": base_roots == candidate_roots,
            "spectra_identical": base_spectra == candidate_spectra,
            "frequencies_identical": base_frequencies == candidate_frequencies,
        }
    payload = {
        "status": "complete",
        "baseline_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "cache_policy": "three cold trials per source use distinct empty NUMBA_CACHE_DIRs; warm A/B uses a separate persistent cache per source",
        "thread_count": 2,
        "cold_process_first_solve": {
            "case": next(iter(benchmark.CASES)),
            "head_trials_seconds": cold_trials["head"],
            "candidate_trials_seconds": cold_trials["candidate"],
            "head_median_seconds": statistics.median(cold_trials["head"]),
            "candidate_median_seconds": statistics.median(cold_trials["candidate"]),
            "candidate_over_head": statistics.median(cold_trials["candidate"])
            / statistics.median(cold_trials["head"]),
        },
        "warm_sample_count_per_case": 48,
        "cases": rows,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
