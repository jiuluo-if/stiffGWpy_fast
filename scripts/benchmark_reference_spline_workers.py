"""Check compiled spline dispatch and Reference output across worker counts."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)
DN_EFF = 0.133732083004284
INDICES = [0, 65, 130, 259]


def _run(model, freqs, workers, z_tail, rtol):
    warmup_seconds = []
    original_warmup = REF._warmup_numba_spline

    def timed_warmup():
        start = time.perf_counter()
        original_warmup()
        warmup_seconds.append(time.perf_counter() - start)

    REF._warmup_numba_spline = timed_warmup
    start = time.perf_counter()
    try:
        result = REF.spectrum_reference(
            model, freqs, DN_EFF, z_tail=z_tail, rtol=rtol,
            workers=workers, backend="cartesian")
    finally:
        REF._warmup_numba_spline = original_warmup
    return result, time.perf_counter() - start, warmup_seconds


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" /
                        "reference_spline_workers_20261104.json")
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be positive")

    archive_path = ROOT / "docs" / "res4_sigma_exact_reference_map_20261104.npz"
    with np.load(archive_path) as archive:
        freqs = np.asarray(archive["freqs"], dtype=float)[INDICES]
        ref_log = np.asarray(archive["reference_log"], dtype=float)[INDICES]
    model = LCDM_SG(**PARAMS, DN_eff=DN_EFF)
    parallel, parallel_s, jit_times = _run(
        model, freqs, args.workers, args.z_tail, args.rtol)
    serial, serial_s, serial_jit_times = _run(
        model, freqs, 1, args.z_tail, args.rtol)
    output_diffs = [float(np.max(np.abs(a - b)))
                    for a, b in zip(parallel[:3], serial[:3])]
    exact_match = all(np.array_equal(a, b)
                      for a, b in zip(parallel[:3], serial[:3]))
    if not exact_match or not np.array_equal(parallel[3], serial[3]):
        raise RuntimeError("Reference outputs differ by worker count")
    omega = parallel[0] - parallel[1]
    if not np.isfinite(omega).all() or np.any(omega <= 0.0):
        raise FloatingPointError("nonpositive/nonfinite spectrum in worker check")
    logs = np.log10(omega)
    errors = np.abs(logs - ref_log)
    reference_path = ROOT / "stiffgwpy_fast" / "reference.py"
    script_path = Path(__file__)
    payload = {
        "schema_version": 1,
        "experiment": "reference_numba_spline_worker_consistency",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "source_script_sha256": hashlib.sha256(
            script_path.read_bytes()).hexdigest(),
        "reference_module_sha256": hashlib.sha256(
            reference_path.read_bytes()).hexdigest(),
        "parameters": {**PARAMS, "DN_eff": DN_EFF},
        "frequency_indices": INDICES,
        "frequencies_log10_hz": freqs.tolist(),
        "z_tail": args.z_tail,
        "rtol": args.rtol,
        "atol": [1e-12, 1e-22, 1e-22],
        "worker_count": args.workers,
        "parallel_cold_or_warm_seconds": parallel_s,
        "parallel_jit_warmup_seconds": jit_times,
        "serial_seconds": serial_s,
        "serial_jit_warmup_seconds": serial_jit_times,
        "bitwise_output_match": exact_match,
        "max_abs_component_diffs": output_diffs,
        "max_abs_spectrum_error_dex_vs_archive": float(errors.max()),
        "resources": telemetry(workers=args.workers, threads=1),
        "semantics": "four selected res4 modes only; worker consistency smoke, not a scaling estimate",
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
