from __future__ import annotations

import hashlib
import importlib
import importlib.util
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

BASE = subprocess.check_output(
    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
).strip()
archive = subprocess.run(
    ["git", "archive", BASE], cwd=ROOT, check=True, stdout=subprocess.PIPE
).stdout

with tempfile.TemporaryDirectory(prefix="stiffgwpy_wkb3_screen_") as temp_dir:
    base_root = Path(temp_dir)
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
        bundle.extractall(base_root, filter="data")
    sys.path.insert(0, str(base_root))
    from stiffgwpy_fast import fast_sgwb as base_fast
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG as BaseModel

    candidate_name = "stiffgwpy_wkb3_candidate"
    candidate_root = base_root / "candidate_source"
    shutil.copytree(
        ROOT / "stiffgwpy_fast",
        candidate_root / "stiffgwpy_fast",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.nbc", "*.nbi"),
    )
    spec = importlib.util.spec_from_file_location(
        candidate_name,
        candidate_root / "stiffgwpy_fast" / "__init__.py",
        submodule_search_locations=[str(candidate_root / "stiffgwpy_fast")],
    )
    assert spec is not None and spec.loader is not None
    package = importlib.util.module_from_spec(spec)
    sys.modules[candidate_name] = package
    spec.loader.exec_module(package)
    candidate_fast = importlib.import_module(candidate_name + ".fast_sgwb")
    CandidateModel = importlib.import_module(candidate_name + ".stiff_SGWB").LCDM_SG

    config_base = base_fast.resolve_config("fast", threads=2)
    config_candidate = candidate_fast.resolve_config("fast", threads=2)
    fixed = {
        "Omega_bh2": 0.12022,
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": 2.100549e-9,
    }
    cases = {
        "default_cr1": dict(r=0.01, cr=1, T_re=2e3, kappa10=0.01),
        "high_kappa_cr1": dict(r=0.01, cr=1, T_re=2e3, kappa10=1.0),
        "cr0_blue_low": dict(
            Omega_bh2=1e6, r=1e-12, n_t=0.5, cr=0, T_re=2e3, kappa10=0.01
        ),
        "cr0_blue_high": dict(
            Omega_bh2=1e6, r=1e-8, n_t=0.5, cr=0, T_re=2e3, kappa10=0.01
        ),
        "cr0_blue_sobol": dict(
            Omega_bh2=1e6,
            r=1.676304587225033e-9,
            n_t=0.4469903271645308,
            cr=0,
            T_re=2e3,
            kappa10=0.01,
        ),
    }

    def run(side: str, params: dict[str, float]) -> dict[str, object]:
        solver, model_type, config = (
            (base_fast, BaseModel, config_base)
            if side == "base"
            else (candidate_fast, CandidateModel, config_candidate)
        )
        model = model_type(**(fixed | params))
        start = time.perf_counter()
        result = solver.SGWB_iter_fast(
            model,
            config=config,
            tol=1e-6,
            kink_split=True,
            physical_guard=False,
            freq_grid="goal",
            frequency_quadrature="pchip",
        )
        elapsed = time.perf_counter() - start
        names = ("f", "log10OmegaGW", "DN_gw", "g2", "w2")
        digests = {
            name: hashlib.sha256(
                np.ascontiguousarray(
                    np.asarray(getattr(model, name), dtype=np.float64)
                ).tobytes()
            ).hexdigest()
            for name in names
        }
        return {
            "ok": result is model,
            "failure": getattr(model, "fast_failure_reason", None),
            "time_s": elapsed,
            "digests": digests,
            "root": float(model.DN_gw[-1]) if len(model.DN_gw) else None,
            "n_freq": int(model.f.size) if hasattr(model, "f") else 0,
            "spectrum": np.asarray(model.log10OmegaGW, dtype=np.float64).tolist()
            if hasattr(model, "log10OmegaGW")
            else [],
        }

    rows = []
    cold_first_calls = {}
    for name, extra in cases.items():
        params = fixed | extra
        for warmup_index in range(2):
            for side in ("base", "candidate"):
                result = run(side, params)
                if not rows and warmup_index == 0:
                    cold_first_calls[side] = result
        pairs = []
        for index in range(64):
            order = ("base", "candidate") if index % 2 == 0 else ("candidate", "base")
            results = {side: run(side, params) for side in order}
            pairs.append({"base": results["base"], "candidate": results["candidate"]})
        base_times = [pair["base"]["time_s"] for pair in pairs]
        candidate_times = [pair["candidate"]["time_s"] for pair in pairs]
        digest_mismatches = {key: 0 for key in pairs[0]["base"]["digests"]}
        root_mismatches = 0
        grid_size_mismatches = 0
        statuses = 0
        for pair in pairs:
            old = pair["base"]
            new = pair["candidate"]
            statuses += old["ok"] != new["ok"] or old["failure"] != new["failure"]
            for key in digest_mismatches:
                digest_mismatches[key] += old["digests"][key] != new["digests"][key]
            root_mismatches += old["root"] != new["root"]
            grid_size_mismatches += old["n_freq"] != new["n_freq"]
            old_spectrum = np.asarray(old["spectrum"], dtype=np.float64)
            new_spectrum = np.asarray(new["spectrum"], dtype=np.float64)
            pair["max_spectrum_absdiff_dex"] = float(
                np.max(np.abs(old_spectrum - new_spectrum))
            )
            pair["root_relative_difference"] = float(
                abs(old["root"] - new["root"])
                / max(abs(old["root"]), np.finfo(np.float64).tiny)
            )
            pair["base"].pop("spectrum")
            pair["candidate"].pop("spectrum")
        rows.append(
            {
                "case": name,
                "pairs": len(pairs),
                "status_mismatches": int(statuses),
                "digest_mismatch_pairs": digest_mismatches,
                "root_mismatch_pairs": int(root_mismatches),
                "frequency_count_mismatch_pairs": int(grid_size_mismatches),
                "base_p50_s": float(np.percentile(base_times, 50)),
                "candidate_p50_s": float(np.percentile(candidate_times, 50)),
                "p50_ratio": float(np.percentile(candidate_times, 50) / np.percentile(base_times, 50)),
                "base_p95_s": float(np.percentile(base_times, 95)),
                "candidate_p95_s": float(np.percentile(candidate_times, 95)),
                "p95_ratio": float(np.percentile(candidate_times, 95) / np.percentile(base_times, 95)),
                "max_spectrum_absdiff_dex": max(
                    pair["max_spectrum_absdiff_dex"] for pair in pairs
                ),
                "max_root_relative_difference": max(
                    pair["root_relative_difference"] for pair in pairs
                ),
                "pairs_raw": pairs,
            }
        )
        print(
            name,
            "status mismatches",
            statuses,
            "digest mismatches",
            digest_mismatches,
            "p50/p95",
            rows[-1]["p50_ratio"],
            rows[-1]["p95_ratio"],
        )

    output = {
        "status": "complete_sameprocess_fullsolve_candidate_screen",
        "baseline_commit": BASE,
        "candidate_revision": "working tree changes on top of baseline_commit",
        "candidate_source_sha256": {
            name: hashlib.sha256((ROOT / "stiffgwpy_fast" / name).read_bytes()).hexdigest()
            for name in ("fast_sgwb.py", "_wkb_transfer.py")
        },
        "candidate_kernel": "production WKB3 defect <= 1e-5 when phase subdivision >= 32; otherwise exact phase transfer",
        "resources": {
            "affinity": affinity,
            "threads": 2,
            "threading_layer": "workqueue",
            "blas_threads": 1,
            "warmups": 2,
            "alternating_pairs": 64,
            "kink_split": True,
            "frequency_grid": "goal",
            "frequency_quadrature": "pchip",
        },
        "cold_first_calls_including_jit": cold_first_calls,
        "cases": rows,
    }
    output_path = ROOT / "docs/research/unified_domain_20261119/benchmark_wkb3_fullsolve_screen.json"
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(output_path)
