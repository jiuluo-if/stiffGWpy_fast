"""Warm-path A/B for the overflow-safe Fast PCHIP slope formula."""

from __future__ import annotations

import hashlib
import json
import platform
import statistics
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

REPEATS = 30
OUTPUT = ROOT / "docs" / "fast_pchip_stable_slopes_ab_20261105.json"
PARAMETERS = {
    "r": 1e-2,
    "cr": 1,
    "T_re": 2e3,
    "kappa10": 1e-2,
    "DN_eff": 0.0,
}


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _old_pchip_slopes(x, y):
    """Original reciprocal-slope implementation for the measured baseline."""
    n = x.size
    hk = x[1:] - x[:-1]
    mk = (y[1:] - y[:-1]) / hk
    if n == 2:
        return np.full(2, mk[0], dtype=np.float64)
    dk = np.zeros(n, dtype=np.float64)
    smk = np.sign(mk)
    condition = (smk[1:] != smk[:-1]) | (mk[1:] == 0.0) | (mk[:-1] == 0.0)
    w1 = 2.0 * hk[1:] + hk[:-1]
    w2 = hk[1:] + 2.0 * hk[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        whmean = (w1 / mk[:-1] + w2 / mk[1:]) / (w1 + w2)
        harmonic = 1.0 / whmean
    dk[1:-1] = np.where(condition, 0.0, harmonic)
    dk[0] = FS._pchip_edge_slope(hk[0], hk[1], mk[0], mk[1])
    dk[-1] = FS._pchip_edge_slope(hk[-1], hk[-2], mk[-1], mk[-2])
    return dk


def _run(slope_function):
    previous = FS._pchip_slopes
    FS._pchip_slopes = slope_function
    model = LCDM_SG(**PARAMETERS)
    caught = []
    started = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as warning_rows:
            warnings.simplefilter("always")
            result = model.SGWB_iter(
                engine="fast",
                fallback=False,
                tol=1e-6,
                z_tail=5.0,
                freq_res=1.0,
                threads=2,
                accuracy_mode="fast",
                physical_guard=False,
            )
        caught = [str(w.message) for w in warning_rows]
    finally:
        FS._pchip_slopes = previous
    if result is not model:
        raise RuntimeError(f"Fast warm-path benchmark failed: {model.fast_failure_reason}")
    elapsed = time.perf_counter() - started
    freqs = np.asarray(model.f, dtype=float)
    omega = np.asarray(model.Ogw_today, dtype=float) - np.asarray(model.Oj_today, dtype=float)
    return {
        "runtime_s": elapsed,
        "DN_gw": float(model.DN_gw[-1]),
        "DN_eff_final": float(model.cosmo_param["DN_eff"]),
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "log10OmegaGW": np.log10(omega).tolist(),
        "warnings": caught,
    }


def main() -> None:
    FS.apply_accuracy_mode("fast")
    warm_started = time.perf_counter()
    warm_old = _run(_old_pchip_slopes)
    warm_new = _run(FS._pchip_slopes)
    warmup_s = time.perf_counter() - warm_started

    times = {"pre_change": [], "stable": []}
    runs = {"pre_change": [], "stable": []}
    for _ in range(REPEATS):
        for key, fn in (("pre_change", _old_pchip_slopes), ("stable", FS._pchip_slopes)):
            result = _run(fn)
            times[key].append(float(result["runtime_s"]))
            runs[key].append(result)

    reference = runs["pre_change"][0]
    max_abs_spectrum_delta = max(
        float(np.max(np.abs(
            np.asarray(old["log10OmegaGW"]) - np.asarray(new["log10OmegaGW"])
        )))
        for old, new in zip(runs["pre_change"], runs["stable"])
    )
    record = {
        "schema_version": 1,
        "experiment": "fast_pchip_stable_slopes_warm_ab",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "parameters": PARAMETERS,
        "fast_settings": {
            "profile": "fast",
            "freq_res": 1.0,
            "z_tail": 5.0,
            "threads": 2,
            "outer_tol": 1e-6,
            "frequency_grid": "goal",
            "frequency_quadrature": "pchip",
        },
        "warmup_s": warmup_s,
        "warmup_outputs": {"pre_change": warm_old, "stable": warm_new},
        "repeat_pairs": REPEATS,
        "frequency_count": reference["frequency_count"],
        "frequency_sha256": reference["frequency_sha256"],
        "timings_s": {
            key: {
                "samples": values,
                "p50": statistics.median(values),
                "p95": float(np.quantile(values, 0.95)),
            }
            for key, values in times.items()
        },
        "max_abs_spectrum_delta_dex": max_abs_spectrum_delta,
        "maximum_relative_DN_gw_delta": max(
            abs(old["DN_gw"] - new["DN_gw"])
            / max(abs(old["DN_gw"]), abs(new["DN_gw"]), np.finfo(float).tiny)
            for old, new in zip(runs["pre_change"], runs["stable"])
        ),
        "warning_counts": {
            key: sum(len(run["warnings"]) for run in values)
            for key, values in runs.items()
        },
        "source_sha256": {
            "fast_module": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "model_api_module": _sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "script": _sha(Path(__file__)),
        },
        "resources": {
            **telemetry(threads=2),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "Warm single-point Fast timing A/B; it does not measure cold process start, batches, or other cosmologies.",
            "The original slope implementation is reproduced inline only as the A/B baseline; the production default is the stable formula.",
        ],
    }
    if any(run["frequency_sha256"] != reference["frequency_sha256"] for group in runs.values() for run in group):
        raise AssertionError("A/B frequency arrays differ")
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "timings_s": {key: {"p50": value["p50"], "p95": value["p95"]}
                      for key, value in record["timings_s"].items()},
        "max_abs_spectrum_delta_dex": max_abs_spectrum_delta,
        "maximum_relative_DN_gw_delta": record["maximum_relative_DN_gw_delta"],
        "warning_counts": record["warning_counts"],
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
