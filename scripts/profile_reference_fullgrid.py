"""Collect per-frequency Reference time/nfev and stage attribution.

This profiling run solves the archived sample2 fixed background on all native
res8 nodes. It is not a paired speed benchmark; use the dedicated A/B driver
for performance comparisons.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "r": 0.19649188544593907,
    "cr": 0.0,
    "n_t": 1.0622501865308032,
    "T_re": 1761.314039059193,
    "kappa10": 0.011598483151630803,
    "DN_re": 10.0,
}
DN_EFF = 2.291411867022825
Z_TAIL = 5.0
Z_START = -12.0
RTOL = 1e-9
REFERENCE_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "reference_fullgrid_profile_compiled_20261105.json"


def _git(*args):
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def main():
    archive = np.load(REFERENCE_NPZ)
    freqs = np.asarray(archive["freqs"], dtype=float)
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    R._warmup_numba_reference_background(
        R._Background(model, DN_EFF)
    )
    timings = {"mode_s": [], "start_search_s": 0.0, "tensor_rhs_s": 0.0, "h2_sigma_s": 0.0}
    counts = {"mode": 0, "tensor_rhs": 0, "h2_sigma": 0}
    records = []

    old_start = R._find_start_N
    old_tensor = R._tensor_orig
    old_h2 = R._H2_and_sigma
    old_solve = R.solve_reference_mode

    def timed_start(*args, **kwargs):
        t0 = time.perf_counter()
        try:
            return old_start(*args, **kwargs)
        finally:
            timings["start_search_s"] += time.perf_counter() - t0

    def timed_tensor(*args, **kwargs):
        t0 = time.perf_counter()
        try:
            return old_tensor(*args, **kwargs)
        finally:
            timings["tensor_rhs_s"] += time.perf_counter() - t0
            counts["tensor_rhs"] += 1

    def timed_h2(*args, **kwargs):
        t0 = time.perf_counter()
        try:
            return old_h2(*args, **kwargs)
        finally:
            timings["h2_sigma_s"] += time.perf_counter() - t0
            counts["h2_sigma"] += 1

    def timed_solve(*args, **kwargs):
        t0 = time.perf_counter()
        result = old_solve(*args, **kwargs)
        elapsed = time.perf_counter() - t0
        timings["mode_s"].append(elapsed)
        counts["mode"] += 1
        records.append({
            "index": counts["mode"] - 1,
            "freq": float(args[1] if len(args) > 1 else kwargs["freq"]),
            "runtime_s": elapsed,
            "nfev": int(result["n_steps"]),
            "used_tail": bool(result["used_tail"]),
            "event_N": result["event_N"],
            "eps_handoff": result["eps_handoff"],
        })
        return result

    R._find_start_N = timed_start
    R._tensor_orig = timed_tensor
    R._H2_and_sigma = timed_h2
    R.solve_reference_mode = timed_solve
    run_started = time.perf_counter()
    try:
        Ogw, Oj, Opgw, used_tail = R.spectrum_reference(
            model,
            freqs,
            DN_EFF,
            z_tail=Z_TAIL,
            rtol=RTOL,
            workers=1,
            backend="cartesian",
            z_start=Z_START,
        )
    finally:
        R._find_start_N = old_start
        R._tensor_orig = old_tensor
        R._H2_and_sigma = old_h2
        R.solve_reference_mode = old_solve
    solve_wall_s = time.perf_counter() - run_started

    integral_started = time.perf_counter()
    integral = R.integrate_spectrum(freqs, Ogw, Oj)
    integral_s = time.perf_counter() - integral_started
    mode_times = np.asarray(timings["mode_s"], dtype=float)
    nfev = np.asarray([row["nfev"] for row in records], dtype=int)
    tail = np.asarray(used_tail, dtype=bool)
    if not (len(records) == len(freqs) == 800 and np.all(np.isfinite(Ogw))):
        raise RuntimeError("profile did not produce the expected finite 800-mode grid")

    try:
        import numba
        import scipy
        scipy_version, numba_version = scipy.__version__, numba.__version__
    except ImportError:
        scipy_version = numba_version = None
    result = {
        "schema_version": 1,
        "experiment": "reference_sample2_fullgrid_compiled_stage_profile",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_module_sha256": hashlib.sha256(
            (ROOT / "stiffgwpy_fast" / "reference.py").read_bytes()
        ).hexdigest(),
        "parameters": PARAMS,
        "DN_eff": DN_EFF,
        "frequency_count": int(freqs.size),
        "frequency_sha256": hashlib.sha256(
            np.asarray(freqs, dtype="<f8").tobytes()
        ).hexdigest(),
        "z_tail": Z_TAIL,
        "z_start": Z_START,
        "rtol": RTOL,
        "backend": "cartesian",
        "workers": 1,
        "jit_warmed_before_profile": True,
        "solve_wall_s": solve_wall_s,
        "integral_s": integral_s,
        "integral_g2": float(integral[0]),
        "quadrature_error": integral[1],
        "interpolation_error": integral[2],
        "mode_runtime_s": {
            "p50": float(np.percentile(mode_times, 50)),
            "p95": float(np.percentile(mode_times, 95)),
            "max": float(np.max(mode_times)),
            "sum": float(np.sum(mode_times)),
        },
        "nfev": {
            "sum": int(np.sum(nfev)),
            "p50": float(np.percentile(nfev, 50)),
            "p95": float(np.percentile(nfev, 95)),
            "max": int(np.max(nfev)),
        },
        "stages": {
            "start_search_s": timings["start_search_s"],
            "tensor_rhs_inclusive_s": timings["tensor_rhs_s"],
            "h2_sigma_inclusive_s": timings["h2_sigma_s"],
            "h2_sigma_call_count": counts["h2_sigma"],
            "tensor_rhs_call_count": counts["tensor_rhs"],
            "note": "Nested times are inclusive and must not be summed; solve_wall_s is authoritative total.",
        },
        "tail_count": int(np.count_nonzero(tail)),
        "log10OmegaGW_min": float(np.min(np.log10(Ogw - Oj))),
        "log10OmegaGW_max": float(np.max(np.log10(Ogw - Oj))),
        "per_frequency": records,
        "resources": {
            **telemetry(threads=1),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy_version,
            "numba": numba_version,
            "workers": 1,
        },
        "limitations": [
            "Single fixed cosmology and fixed background; this profile is not a speed comparison.",
            "Stage times include Python wrapper overhead and overlap through nested calls.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "per_frequency"}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
