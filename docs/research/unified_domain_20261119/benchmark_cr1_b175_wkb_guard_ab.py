"""Same-process CR1 full-solve A/B for the WKB basis conditioning guard."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).with_name("benchmark_cr1_b175_wkb_guard_ab.json")
BASE_SHA = "4549f5951f781fd275b171a770853cf65407e1c8"
REPEATS = 32
WARMUPS = 2
PARAMS = {
    "Omega_ch2": 7.0e5,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 0.01,
    "cr": 1.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 0.01,
}
CASES = (
    ("cr1_square_high_corner", dict(Omega_bh2=7.0e5, n_t=0.0)),
    ("cr1_square_low_holdout", dict(Omega_bh2=2.0e5, Omega_ch2=2.5e5, n_t=0.0)),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _load_baseline_modules(temp_dir: Path):
    helper_source = subprocess.check_output(
        ["git", "show", f"{BASE_SHA}:stiffgwpy_fast/_wkb_transfer.py"],
        cwd=ROOT,
        text=True,
    )
    helper_path = temp_dir / "_wkb_transfer_baseline.py"
    helper_path.write_text(helper_source, encoding="utf-8")
    helper_name = "stiffgwpy_fast._wkb_transfer_baseline"
    helper_spec = importlib.util.spec_from_file_location(helper_name, helper_path)
    if helper_spec is None or helper_spec.loader is None:
        raise RuntimeError("could not load baseline WKB transfer module")
    helper_module = importlib.util.module_from_spec(helper_spec)
    import sys

    sys.modules[helper_name] = helper_module
    helper_spec.loader.exec_module(helper_module)

    source = (ROOT / "stiffgwpy_fast" / "fast_sgwb.py").read_text(encoding="utf-8")
    original = "from ._wkb_transfer import wkb_transfer_embedded as _wkb_transfer_embedded"
    replacement = (
        "from ._wkb_transfer_baseline import "
        "wkb_transfer_embedded as _wkb_transfer_embedded"
    )
    if source.count(original) != 1:
        raise RuntimeError("expected exactly one WKB transfer import")
    baseline_path = temp_dir / "fast_sgwb_baseline_wkb.py"
    baseline_path.write_text(source.replace(original, replacement), encoding="utf-8")
    baseline_name = "stiffgwpy_fast.fast_sgwb_baseline_wkb"
    baseline_spec = importlib.util.spec_from_file_location(baseline_name, baseline_path)
    if baseline_spec is None or baseline_spec.loader is None:
        raise RuntimeError("could not load baseline Fast module")
    baseline_module = importlib.util.module_from_spec(baseline_spec)
    sys.modules[baseline_name] = baseline_module
    baseline_spec.loader.exec_module(baseline_module)
    return baseline_module


def _solve(module, params, config):
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    model = LCDM_SG(**params)
    started = time.perf_counter()
    result = module.SGWB_iter_fast(
        model,
        config=config,
        tol=1.0e-6,
        physical_guard=True,
    )
    elapsed = time.perf_counter() - started
    if result is not model:
        return {"success": False, "reason": model.fast_failure_reason, "seconds": elapsed}
    return {
        "success": True,
        "seconds": elapsed,
        "root": float(model.DN_gw[-1]),
        "frequencies": np.asarray(model.f, dtype=np.float64),
        "spectrum": np.asarray(model.log10OmegaGW, dtype=np.float64),
    }


def _compact(result):
    return {
        key: value
        for key, value in result.items()
        if key not in ("frequencies", "spectrum")
    } | {
        "frequency_sha256": _digest(result["frequencies"]),
        "spectrum_sha256": _digest(result["spectrum"]),
    }


def main() -> None:
    import sys

    fresh_numba_cache = Path(tempfile.mkdtemp(prefix="stiffgwpy_cr1_numba_cache_"))
    os.environ["NUMBA_CACHE_DIR"] = str(fresh_numba_cache)
    sys.path.insert(0, str(ROOT))
    import stiffgwpy_fast.fast_sgwb as candidate

    try:
        import psutil

        available = psutil.Process().cpu_affinity()
        if 0 in available and 1 in available:
            psutil.Process().cpu_affinity([0, 1])
    except (ImportError, AttributeError, OSError):
        pass

    with tempfile.TemporaryDirectory(prefix="stiffgwpy_wkb_guard_ab_") as temp:
        baseline = _load_baseline_modules(Path(temp))
        config = candidate.resolve_config("fast", threads=2)
        payload_cases = []
        with threadpool_limits(limits=1):
            for case_index, (case_id, override) in enumerate(CASES):
                params = dict(PARAMS, **override)
                first_base = _compact(_solve(baseline, params, config))
                first_candidate = _compact(_solve(candidate, params, config))
                if not first_base["success"] or not first_candidate["success"]:
                    payload_cases.append(
                        {
                            "case_id": case_id,
                            "params": params,
                            "first_call_compilation_cold": case_index == 0,
                            "first_base": first_base,
                            "first_candidate": first_candidate,
                            "status": "baseline_or_candidate_failed",
                        }
                    )
                    continue

                for _ in range(WARMUPS):
                    _solve(baseline, params, config)
                    _solve(candidate, params, config)
                timings = {"baseline_seconds": [], "candidate_seconds": []}
                roots_base = []
                roots_candidate = []
                spectrum_delta_max = 0.0
                frequency_mismatch = 0
                raw_pairs = []
                for index in range(REPEATS):
                    order = (
                        (("baseline", baseline), ("candidate", candidate))
                        if index % 2 == 0
                        else (("candidate", candidate), ("baseline", baseline))
                    )
                    solved = {
                        name: _solve(module, params, config)
                        for name, module in order
                    }
                    b = solved["baseline"]
                    c = solved["candidate"]
                    if not b["success"] or not c["success"]:
                        raise RuntimeError(f"paired solve failed at {case_id} pair {index}")
                    timings["baseline_seconds"].append(b["seconds"])
                    timings["candidate_seconds"].append(c["seconds"])
                    roots_base.append(b["root"])
                    roots_candidate.append(c["root"])
                    frequency_mismatch += int(
                        not np.array_equal(b["frequencies"], c["frequencies"])
                    )
                    if b["frequencies"].shape != c["frequencies"].shape:
                        raise RuntimeError("frequency grid shape changed within paired run")
                    spectrum_delta = float(
                        np.max(np.abs(b["spectrum"] - c["spectrum"]))
                    )
                    spectrum_delta_max = max(spectrum_delta_max, spectrum_delta)
                    raw_pairs.append(
                        {
                            "index": index,
                            "first_side": order[0][0],
                            "baseline_seconds": b["seconds"],
                            "candidate_seconds": c["seconds"],
                            "baseline_root": b["root"],
                            "candidate_root": c["root"],
                            "spectrum_absdiff_max_dex": spectrum_delta,
                        }
                    )

                base_times = np.asarray(timings["baseline_seconds"])
                candidate_times = np.asarray(timings["candidate_seconds"])
                pair_ratios = candidate_times / base_times
                payload_cases.append(
                    {
                        "case_id": case_id,
                        "params": params,
                        "pairs": REPEATS,
                        "first_call_compilation_cold": case_index == 0,
                        "first_base": first_base,
                        "first_candidate": first_candidate,
                        "base_p50_seconds": float(np.median(base_times)),
                        "candidate_p50_seconds": float(np.median(candidate_times)),
                        "p50_ratio": float(np.median(candidate_times) / np.median(base_times)),
                        "base_p95_seconds": float(np.quantile(base_times, 0.95)),
                        "candidate_p95_seconds": float(np.quantile(candidate_times, 0.95)),
                        "p95_ratio": float(
                            np.quantile(candidate_times, 0.95)
                            / np.quantile(base_times, 0.95)
                        ),
                        "paired_ratio_p50": float(np.median(pair_ratios)),
                        "paired_ratio_p95": float(np.quantile(pair_ratios, 0.95)),
                        "root_relative_difference_max": float(
                            np.max(
                                np.abs(np.asarray(roots_candidate) - roots_base)
                                / np.maximum(np.abs(roots_base), 1e-300)
                            )
                        ),
                        "frequency_digest_mismatch_pairs": frequency_mismatch,
                        "spectrum_absdiff_max_dex": spectrum_delta_max,
                        "timing_pairs": raw_pairs,
                    }
                )
                print(
                    case_id,
                    payload_cases[-1]["p50_ratio"],
                    payload_cases[-1]["p95_ratio"],
                    spectrum_delta_max,
                    flush=True,
                )

        payload = {
            "status": "complete_wkb_basis_guard_ab",
            "base_commit": BASE_SHA,
            "candidate_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            "baseline_helper_sha256": hashlib.sha256(
                subprocess.check_output(
                    ["git", "show", f"{BASE_SHA}:stiffgwpy_fast/_wkb_transfer.py"],
                    cwd=ROOT,
                )
            ).hexdigest(),
            "candidate_helper_sha256": _sha256(ROOT / "stiffgwpy_fast" / "_wkb_transfer.py"),
            "resources": {
                "affinity": [0, 1],
                "numba_threads": 2,
                "blas_threads": 1,
                "threading_layer": str(candidate.threading_layer()),
                "warmups_per_side": WARMUPS,
                "alternating_pairs_per_case": REPEATS,
            "numba_cache": "isolated fresh temporary cache for cold-call comparison",
                "frequency_grid": "goal",
                "frequency_quadrature": "pchip",
                "physical_guard": True,
            },
            "cases": payload_cases,
        }
        OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        shutil.rmtree(fresh_numba_cache, ignore_errors=True)
        print(f"artifact: {OUTPUT}")


if __name__ == "__main__":
    main()
