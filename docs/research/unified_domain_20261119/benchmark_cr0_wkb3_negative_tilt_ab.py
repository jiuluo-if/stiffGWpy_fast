"""Same-source paired A/B for WKB3 overhead on the CR0 negative-tilt region."""

from __future__ import annotations

import hashlib
import importlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import psutil

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
process = psutil.Process()
affinity = process.cpu_affinity()[:2]
process.cpu_affinity(affinity)
for key, value in {
    "NUMBA_NUM_THREADS": "2",
    "NUMBA_THREADING_LAYER": "workqueue",
    "FAST_THREADS": "2",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}.items():
    os.environ[key] = value

import numpy as np  # noqa: E402

BASE_SHA = "63aee9d180f8c570276f07177afe241dc9b8b0cd"
BASE_FAST_SHA256 = hashlib.sha256(
    subprocess.run(
        ["git", "show", f"{BASE_SHA}:stiffgwpy_fast/fast_sgwb.py"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
).hexdigest()
REPEATS = 64
FIXED = {
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
    "old_region_1m_minus_half": {"Omega_bh2": 1.0e6, "n_t": -0.5},
    "new_region_1p5m_minus_half": {"Omega_bh2": 1.5e6, "n_t": -0.5},
    "new_region_2m_minus_half": {"Omega_bh2": 2.0e6, "n_t": -0.5},
    "new_region_2m_zero": {"Omega_bh2": 2.0e6, "n_t": 0.0},
    "new_region_sobol": {"Omega_bh2": 1.5929892305284739e6, "n_t": -0.1713485741056502},
    "old_region_750k_holdout": {"Omega_bh2": 7.5e5, "n_t": -0.125},
}

archive = subprocess.run(
    ["git", "archive", BASE_SHA], cwd=ROOT, check=True, stdout=subprocess.PIPE
).stdout

with tempfile.TemporaryDirectory(prefix="stiffgwpy_cr0_wkb3_ab_") as temp:
    temp_root = Path(temp)
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
        tar.extractall(temp_root, filter="data")
    sys.path.insert(0, str(temp_root))
    from stiffgwpy_fast import fast_sgwb as base_fast
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG as BaseModel

    candidate_root = temp_root / "candidate"
    import shutil

    shutil.copytree(
        ROOT / "stiffgwpy_fast",
        candidate_root / "stiffgwpy_fast",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.nbc", "*.nbi"),
    )
    candidate_name = "stiffgwpy_cr0_wkb3_candidate"
    spec = importlib.util.spec_from_file_location(
        candidate_name,
        candidate_root / "stiffgwpy_fast" / "__init__.py",
        submodule_search_locations=[str(candidate_root / "stiffgwpy_fast")],
    )
    assert spec is not None and spec.loader is not None
    pkg = importlib.util.module_from_spec(spec)
    sys.modules[candidate_name] = pkg
    spec.loader.exec_module(pkg)
    candidate_fast = importlib.import_module(candidate_name + ".fast_sgwb")
    CandidateModel = importlib.import_module(candidate_name + ".stiff_SGWB").LCDM_SG

    base_config = base_fast.resolve_config("fast", threads=2)
    candidate_config = candidate_fast.resolve_config("fast", threads=2)

    def run(side: str, override: dict[str, float]) -> dict[str, object]:
        solver, model_type, config = (
            (base_fast, BaseModel, base_config)
            if side == "base"
            else (candidate_fast, CandidateModel, candidate_config)
        )
        model = model_type(**(FIXED | override))
        started = time.perf_counter()
        result = solver.SGWB_iter_fast(
            model,
            config=config,
            tol=1.0e-6,
            physical_guard=True,
            freq_grid="goal",
            frequency_quadrature="pchip",
        )
        seconds = time.perf_counter() - started
        row: dict[str, object] = {
            "success": result is model,
            "failure_reason": getattr(model, "fast_failure_reason", None),
            "seconds": seconds,
        }
        if result is model:
            freq = np.ascontiguousarray(model.f, dtype=np.float64)
            spec_values = np.ascontiguousarray(model.log10OmegaGW, dtype=np.float64)
            row.update(
                {
                    "root": float(model.DN_gw[-1]),
                    "n_freq": int(freq.size),
                    "frequency_sha256": hashlib.sha256(freq.tobytes()).hexdigest(),
                    "spectrum_sha256": hashlib.sha256(spec_values.tobytes()).hexdigest(),
                    "spectrum": spec_values.tolist(),
                }
            )
        return row

    rows = []
    for name, override in CASES.items():
        params = FIXED | override
        cold = {side: run(side, override) for side in ("base", "candidate")}
        for _ in range(2):
            run("base", override)
            run("candidate", override)
        pairs = []
        for i in range(REPEATS):
            order = ("base", "candidate") if i % 2 == 0 else ("candidate", "base")
            result = {side: run(side, override) for side in order}
            pairs.append({"base": result["base"], "candidate": result["candidate"]})

        base_times = [float(pair["base"]["seconds"]) for pair in pairs]
        candidate_times = [float(pair["candidate"]["seconds"]) for pair in pairs]
        root_diff = 0.0
        spectrum_diff = 0.0
        digest_changes = {key: 0 for key in ("frequency_sha256", "spectrum_sha256")}
        status_diff = 0
        for pair in pairs:
            old, new = pair["base"], pair["candidate"]
            status_diff += old["success"] != new["success"] or old["failure_reason"] != new["failure_reason"]
            if old["success"] and new["success"]:
                root_diff = max(
                    root_diff,
                    abs(old["root"] - new["root"]) / max(abs(old["root"]), 1.0e-300),
                )
                for key in digest_changes:
                    digest_changes[key] += old[key] != new[key]
                spectrum_diff = max(
                    spectrum_diff,
                    float(
                        np.max(
                            np.abs(
                                np.asarray(old["spectrum"], dtype=np.float64)
                                - np.asarray(new["spectrum"], dtype=np.float64)
                            )
                        )
                    ),
                )
            old.pop("spectrum", None)
            new.pop("spectrum", None)

        row = {
            "case": name,
            "params": params,
            "pairs": REPEATS,
            "baseline_status": cold["base"]["failure_reason"] or "success",
            "candidate_status": cold["candidate"]["failure_reason"] or "success",
            "status_mismatch_pairs": int(status_diff),
            "candidate_base_root_relative_difference_max": float(root_diff),
            "candidate_base_spectrum_absdiff_max_dex": float(spectrum_diff),
            "digest_mismatch_pairs": digest_changes,
            "base_p50_seconds": float(np.percentile(base_times, 50)),
            "candidate_p50_seconds": float(np.percentile(candidate_times, 50)),
            "p50_ratio": float(np.percentile(candidate_times, 50) / np.percentile(base_times, 50)),
            "base_p95_seconds": float(np.percentile(base_times, 95)),
            "candidate_p95_seconds": float(np.percentile(candidate_times, 95)),
            "p95_ratio": float(np.percentile(candidate_times, 95) / np.percentile(base_times, 95)),
            "cold_first_calls": cold,
            "pairs_raw": pairs,
        }
        rows.append(row)
        print(name, row["baseline_status"], row["candidate_status"], row["p50_ratio"], row["p95_ratio"], flush=True)

    result = {
        "status": "complete_sameprocess_negative_tilt_ab",
        "baseline_commit": BASE_SHA,
        "baseline_source_sha256": {"fast_sgwb.py": BASE_FAST_SHA256},
        "candidate_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "candidate_source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "_wkb_transfer.py")
        },
        "resources": {
            "affinity": affinity,
            "threads": 2,
            "threading_layer": "workqueue",
            "blas_threads": 1,
            "alternating_pairs_per_case": REPEATS,
            "warmups_per_side": 2,
            "frequency_grid": "goal",
            "frequency_quadrature": "pchip",
            "physical_guard": True,
        },
        "cases": rows,
    }
    output = ROOT / "docs/research/unified_domain_20261119/benchmark_cr0_wkb3_negative_tilt_ab.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(output)
