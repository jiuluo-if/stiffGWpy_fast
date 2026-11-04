"""Measure normal-signal cost of the subnormal-safe Reference integrator."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy import integrate, interpolate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402

REPEATS = 24
FREQUENCY_COUNT = 800
OUTPUT = ROOT / "docs" / "reference_integrate_spectrum_scale_ab_20261105.json"


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _old_integrate(freqs, ogw, oj, ln10v=math.log(10.0)):
    """Pre-change implementation used only as an A/B baseline."""
    freqs = np.asarray(freqs, dtype=float)
    integrand = np.asarray(ogw, dtype=float) - np.asarray(oj, dtype=float)
    order = np.argsort(freqs)
    x = freqs[order]
    y = integrand[order]
    if x.size < 4:
        g2 = ln10v * np.trapezoid(y, x) if x.size >= 2 else 0.0
        return float(g2), None, None
    spl = interpolate.PchipInterpolator(x, y)

    def _quad(func, breakpoints):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            limit = max(800, len(breakpoints) + 50)
            return integrate.quad(
                func, float(x[0]), float(x[-1]), epsabs=0.0,
                epsrel=1e-10, limit=limit, points=breakpoints,
            )

    val, abserr = _quad(lambda ell: float(spl(ell)), x)
    half = x[::2]
    yh = y[::2]
    if x.size % 2 == 0:
        half = np.append(half, x[-1])
        yh = np.append(yh, y[-1])
    if half.size >= 4 and half.size < x.size and np.all(np.diff(half) > 0):
        spl_h = interpolate.PchipInterpolator(half, yh)
        val_h, _ = _quad(lambda ell: float(spl_h(ell)), half)
        ierr = ln10v * abs(val - val_h)
    else:
        ierr = None
    return float(ln10v * val), float(ln10v * abs(abserr)), ierr


def main() -> None:
    freqs = np.linspace(-18.5, 12.2, FREQUENCY_COUNT)
    spectrum = 10.0 ** (-112.0 + 0.35 * np.sin(freqs * 0.73) + 0.04 * freqs)
    zero = np.zeros_like(spectrum)
    implementations = {"pre_scale": _old_integrate, "current": R.integrate_spectrum}
    samples: dict[str, list[float]] = {name: [] for name in implementations}
    outputs: dict[str, tuple[float, float, float | None]] = {}

    for _ in range(REPEATS):
        for name, function in implementations.items():
            started = time.perf_counter()
            value = function(freqs, spectrum, zero)
            samples[name].append(time.perf_counter() - started)
            outputs[name] = value

    maximum_relative_output_delta = max(
        abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), np.finfo(float).tiny)
        for a, b in zip(outputs["pre_scale"], outputs["current"])
    )
    record = {
        "schema_version": 1,
        "experiment": "reference_integrate_spectrum_scale_ab",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "repeat_pairs": REPEATS,
        "frequency_count": FREQUENCY_COUNT,
        "frequency_range": [float(freqs[0]), float(freqs[-1])],
        "frequency_sha256": hashlib.sha256(freqs.astype("<f8").tobytes()).hexdigest(),
        "spectrum_definition": "10**(-112 + 0.35*sin(0.73*f) + 0.04*f)",
        "spectrum_min": float(spectrum.min()),
        "spectrum_max": float(spectrum.max()),
        "spectrum_is_subnormal": bool(np.any((spectrum != 0.0) & (spectrum < np.finfo(float).tiny))),
        "zero_j_lambda_sha256": hashlib.sha256(zero.astype("<f8").tobytes()).hexdigest(),
        "outputs": {
            name: [None if value is None else float(value) for value in result]
            for name, result in outputs.items()
        },
        "maximum_relative_output_delta": maximum_relative_output_delta,
        "timings_s": {
            name: {
                "samples": values,
                "p50": statistics.median(values),
                "p95": float(np.quantile(values, 0.95)),
            }
            for name, values in samples.items()
        },
        "source_sha256": {
            "reference_module": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "script": _sha(Path(__file__)),
        },
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "limitations": [
            "Synthetic normal-scale PCHIP integral microbenchmark, not a cosmological full-grid runtime comparison.",
            "Alternating A/B on one process; timings are scoped to this machine and input size.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "timings_s": {key: {"p50": value["p50"], "p95": value["p95"]}
                      for key, value in record["timings_s"].items()},
        "maximum_relative_output_delta": maximum_relative_output_delta,
        "output": str(OUTPUT),
    }, indent=2))


if __name__ == "__main__":
    main()
