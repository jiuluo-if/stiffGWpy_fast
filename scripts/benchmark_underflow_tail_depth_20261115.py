"""Benchmark z_tail=5/8 at one frozen root and native 937-node grid."""

from __future__ import annotations

import gc
import hashlib
import json
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import underflow_grid_convergence_20261115 as grid_run  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_tail_depth_benchmark_20261115.json"
TAILS = (5.0, 8.0)
THREADS = 1
WARMUPS_PER_TAIL = 2
PAIRED_REPEATS = 12


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rss_bytes(process) -> int | None:
    return None if process is None else int(process.memory_info().rss)


def _run_once(params: dict, freqs: np.ndarray, frozen_dn: float, z_tail: float, process):
    model = LCDM_SG(**params)
    model.cosmo_param["DN_eff"] = frozen_dn
    model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def locked_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = frozen_dn
        return original_gen(*args, **kwargs)

    def fixed_grid(*args, **kwargs):
        return freqs.copy()

    model.gen_expansion = locked_gen
    FA.goal_oriented_freqs = fixed_grid
    FA.construct_fast_frequency_grid = fixed_grid
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=0.99,
            z_tail=z_tail,
            threads=THREADS,
            accuracy_mode="fast",
            freq_grid="goal",
            physical_guard=False,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    if result is None or not np.array_equal(model.f, freqs):
        raise RuntimeError(f"Fast failed at z_tail={z_tail}: {model.fast_failure_reason}")
    return {
        "runtime_s": elapsed,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "map_residual": float(model.DN_gw[-1] - frozen_dn),
        "log10_OmegaGW": np.asarray(model.log10OmegaGW, dtype=np.float64),
        "rss_bytes_after": _rss_bytes(process),
    }


def _summary(rows: list[dict]) -> dict:
    times = np.asarray([row["runtime_s"] for row in rows], dtype=np.float64)
    spectra = np.stack([row["log10_OmegaGW"] for row in rows])
    reference_spectrum = spectra[0]
    return {
        "repeat_count": int(times.size),
        "runtime_p50_s": float(np.quantile(times, 0.50)),
        "runtime_p95_s": float(np.quantile(times, 0.95)),
        "runtime_median_s": float(statistics.median(times.tolist())),
        "peak_rss_bytes": max(
            (row["rss_bytes_after"] for row in rows if row["rss_bytes_after"] is not None),
            default=None,
        ),
        "max_repeat_spectrum_spread_dex": float(
            np.max(np.abs(spectra - reference_spectrum))
        ),
        "map_residual_max_abs": float(max(abs(row["map_residual"]) for row in rows)),
        "log10_OmegaGW": reference_spectrum.tolist(),
    }


def main() -> int:
    data = json.loads(INPUT.read_text(encoding="utf-8"))
    params = data["parameters"]
    frozen_dn = float(data["background"]["common_DN_eff"])
    freqs = np.asarray(data["frequency_grids"]["937"]["log10_hz"], dtype=np.float64)
    expected_grid_sha = data["frequency_grids"]["937"]["sha256_f64le"]
    if grid_run._grid_sha256(freqs) != expected_grid_sha:
        raise AssertionError("937-node frequency grid hash mismatch")
    try:
        import psutil
    except ImportError:
        process = None
        psutil_version = None
    else:
        process = psutil.Process()
        psutil_version = psutil.__version__

    rss_before = _rss_bytes(process)
    cold_calls = {}
    for z_tail in TAILS:
        cold_calls[str(int(z_tail))] = _run_once(
            params, freqs, frozen_dn, z_tail, process
        )
        gc.collect()

    warmup_calls = {str(int(z)): [] for z in TAILS}
    for index in range(WARMUPS_PER_TAIL):
        order = TAILS if index % 2 == 0 else tuple(reversed(TAILS))
        for z_tail in order:
            warmup_calls[str(int(z_tail))].append(
                _run_once(params, freqs, frozen_dn, z_tail, process)
            )
            gc.collect()

    measured = {str(int(z)): [] for z in TAILS}
    for pair in range(PAIRED_REPEATS):
        order = TAILS if pair % 2 == 0 else tuple(reversed(TAILS))
        for z_tail in order:
            measured[str(int(z_tail))].append(
                _run_once(params, freqs, frozen_dn, z_tail, process)
            )
            gc.collect()

    summaries = {key: _summary(rows) for key, rows in measured.items()}
    serializable_summaries = {
        key: {
            **summary,
            "log10_OmegaGW": np.asarray(
                summary["log10_OmegaGW"], dtype=np.float64
            ).tolist(),
        }
        for key, summary in summaries.items()
    }
    spectra = {
        key: np.asarray(summary["log10_OmegaGW"], dtype=np.float64)
        for key, summary in summaries.items()
    }
    delta = spectra["8"] - spectra["5"]
    output = {
        "schema_version": 1,
        "status": "complete_same_root_alternating_tail_depth_benchmark",
        "parameters": params,
        "frozen_DN_eff": frozen_dn,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": expected_grid_sha,
        "threads": THREADS,
        "warmups_per_tail": WARMUPS_PER_TAIL,
        "paired_repeats": PAIRED_REPEATS,
        "cold_first_call": {
            key: {k: v for k, v in row.items() if k != "log10_OmegaGW"}
            for key, row in cold_calls.items()
        },
        "warmup_calls": {
            key: [{k: v for k, v in row.items() if k != "log10_OmegaGW"} for row in rows]
            for key, rows in warmup_calls.items()
        },
        "steady_state": serializable_summaries,
        "paired_runtime_ratio_p50_z8_over_z5": (
            summaries["8"]["runtime_p50_s"] / summaries["5"]["runtime_p50_s"]
        ),
        "z8_minus_z5_spectrum_delta_dex": {
            "max_abs": float(np.max(np.abs(delta))),
            "p95_abs": float(np.quantile(np.abs(delta), 0.95)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "psutil": psutil_version,
            "rss_before_bytes": rss_before,
        },
        "source_sha256": {
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "grid_input_json": _sha(INPUT),
            "grid_runner": _sha(ROOT / "scripts" / "underflow_grid_convergence_20261115.py"),
            "benchmark_script": _sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "runtime_p50_p95_s": {
            key: [row["runtime_p50_s"], row["runtime_p95_s"]]
            for key, row in summaries.items()
        },
        "paired_runtime_ratio_p50_z8_over_z5": output[
            "paired_runtime_ratio_p50_z8_over_z5"
        ],
        "spectrum_delta_max_p95_dex": [
            output["z8_minus_z5_spectrum_delta_dex"]["max_abs"],
            output["z8_minus_z5_spectrum_delta_dex"]["p95_abs"],
        ],
        "repeat_spread_dex": {
            key: row["max_repeat_spectrum_spread_dex"]
            for key, row in summaries.items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
