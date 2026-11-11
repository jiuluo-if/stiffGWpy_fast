"""Compare a 235-mode underflow run with a representable twin: 对比同规模下溢与可表示算例。"""
from __future__ import annotations

import hashlib
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, limit_affinity, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402
import psutil  # noqa: E402

from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402


def _sha(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def _model_params(a_s: float) -> dict:
    return {
        "Omega_bh2": 0.0223828,
        "Omega_ch2": 0.1201075,
        "H0": 67.32117,
        "DN_eff": 0.0,
        "A_s": a_s,
        "r": 1e-316,
        "n_t": 17.38,
        "cr": 0.0,
        "T_re": 2e3,
        "DN_re": 10.0,
        "kappa10": 1e-2,
    }


def _run(params: dict) -> dict:
    model = LCDM_SG(**params)
    start = time.perf_counter()
    FS.SGWB_iter_fast(
        model, kink_split=True, freq_grid="construct", physical_guard=False
    )
    elapsed = time.perf_counter() - start
    return {
        "seconds": elapsed,
        "frequency_count": int(len(model.f)),
        "frequency_sha256": _sha(model.f),
        "root_DN_gw": float(model.DN_gw[-1]),
        "converged": bool(model.SGWB_converge),
        "failure_reason": model.fast_failure_reason,
    }


def main() -> None:
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    process = psutil.Process()
    selected = limit_affinity(process, 2)
    cases = {
        "underflow": _model_params(2.100549e-9),
        # A_s is raised only enough to make r*A_s the smallest positive float64.
        "representable": _model_params(3.0e-8),
    }
    for params in cases.values():
        _run(params)
        _run(params)  # JIT and first-allocation warmup
    pairs = 100
    samples = {name: [] for name in cases}
    for pair in range(pairs):
        order = ("underflow", "representable")
        if pair % 2:
            order = order[::-1]
        for name in order:
            samples[name].append(_run(cases[name]))
    summary = {}
    for name, rows in samples.items():
        times = [row["seconds"] for row in rows]
        summary[name] = {
            "p50_ms": statistics.median(times) * 1e3,
            "p95_ms": float(np.percentile(times, 95)) * 1e3,
            "frequency_counts": sorted({row["frequency_count"] for row in rows}),
            "frequency_sha256": sorted({row["frequency_sha256"] for row in rows}),
            "root_DN_gw": rows[-1]["root_DN_gw"],
            "all_converged": all(row["converged"] for row in rows),
            "failure_reasons": sorted({str(row["failure_reason"]) for row in rows}),
            "samples_ms": [value * 1e3 for value in times],
        }
    payload = {
        "schema_version": 1,
        "experiment": "same-background 235-mode cost comparison",
        "git_head": "a5c7afdce50074ac16a1555c0bc8a204ddae6688",
        "fast_sgwb_sha256": hashlib.sha256(
            (ROOT / "stiffgwpy_fast" / "fast_sgwb.py").read_bytes()
        ).hexdigest(),
        "threads": 2,
        "pairs": pairs,
        "affinity_selected": selected,
        "resources": telemetry(process, threads=2),
        "cases": summary,
        "underflow_over_representable_p50_ratio": (
            summary["underflow"]["p50_ms"] / summary["representable"]["p50_ms"]
        ),
        "underflow_over_representable_p95_ratio": (
            summary["underflow"]["p95_ms"] / summary["representable"]["p95_ms"]
        ),
        "scope_limit": "same parameter background and 235 nodes, but the representable twin has a larger tensor amplitude and therefore a different self-consistent DN_gw root",
    }
    out = ROOT / "docs" / "gated_phase_underflow_same_scale_cost_terminalfix_20261108.json"
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
