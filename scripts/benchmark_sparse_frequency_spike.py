"""Standalone sparse-frequency propagation and spectrum reconstruction screen.

The production API and frequency grid remain untouched.  This diagnostic
propagates only every second goal-grid mode, retains both endpoints, and
reconstructs the full output spectrum with a monotone cubic interpolant.  It
is intentionally stopped at fixed ``DN_eff``: outer self-consistency is only
worth implementing if this cheaper full-spectrum screen passes first.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from scipy.interpolate import PchipInterpolator  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import _make_args, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def _mode_indices(size, stride=2):
    size = int(size)
    stride = int(stride)
    if size < 2 or stride < 1:
        raise ValueError("size must be >= 2 and stride must be positive")
    indices = np.arange(0, size, stride, dtype=np.int64)
    if indices[-1] != size - 1:
        indices = np.append(indices, size - 1)
    return indices


def _reconstruct(full_x, sparse_x, sparse_values):
    full_x = np.asarray(full_x, dtype=np.float64)
    sparse_x = np.asarray(sparse_x, dtype=np.float64)
    sparse_values = np.asarray(sparse_values, dtype=np.float64)
    order = np.argsort(sparse_x, kind="stable")
    sparse_x = sparse_x[order]
    sparse_values = sparse_values[order]
    interpolator = PchipInterpolator(
        sparse_x,
        sparse_values,
        extrapolate=False,
    )
    reconstructed = np.asarray(interpolator(full_x), dtype=np.float64)
    for sparse_index, sparse_frequency in enumerate(sparse_x):
        exact = np.flatnonzero(full_x == sparse_frequency)
        reconstructed[exact] = sparse_values[sparse_index]
    return reconstructed


def _digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _sparse_args(common, indices):
    (Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
     fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail, _, Sv,
     phase_max, kink_index, kink_fraction, phi_re) = common
    shape = (len(indices), n_coarse)
    return (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s[indices], z0s[indices],
        P_t[indices], ev_minus, fp_minus, fp_freq[indices],
        assemble, n_coarse, col_step, h, z_tail,
        np.zeros(shape), np.zeros(shape), np.zeros(shape), None, Sv,
        phase_max, np.full(len(indices), -1.0), kink_index, kink_fraction,
        phi_re,
    )


def _full_spectrum(model, args, sparse=False, indices=None):
    if not sparse:
        FS.solve_kernel(*args)
        return args[16][:, -1] - args[17][:, -1]
    FS.solve_kernel(*args)
    full_x = np.asarray(model.f, dtype=np.float64)
    sparse_x = full_x[indices]
    return _reconstruct(full_x, sparse_x, args[16][:, -1] - args[17][:, -1])


def _stats(values):
    values = np.asarray(values, dtype=np.float64)
    return {
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def run_case(name, repeats=25, threads=2, stride=2):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    indices = _mode_indices(len(model.f), stride)
    candidate = _sparse_args(common, indices)
    base = _full_spectrum(model, baseline)
    sparse = _full_spectrum(model, candidate, sparse=True, indices=indices)
    positive_base = np.maximum(base, 1e-300)
    positive_sparse = np.maximum(sparse, 1e-300)
    dex = np.abs(np.log10(positive_sparse) - np.log10(positive_base))
    dn_base = FS.integrate_frequency_pchip(model.f, base)
    dn_sparse = FS.integrate_frequency_pchip(model.f, sparse)
    base_times = []
    candidate_times = []
    for _ in range(int(repeats)):
        baseline[16].fill(0.0)
        baseline[17].fill(0.0)
        baseline[18].fill(0.0)
        candidate[16].fill(0.0)
        candidate[17].fill(0.0)
        candidate[18].fill(0.0)
        start = time.perf_counter()
        _full_spectrum(model, baseline)
        base_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _full_spectrum(model, candidate, sparse=True, indices=indices)
        candidate_times.append(time.perf_counter() - start)
    return {
        "case": name,
        "stride": int(stride),
        "full_modes": int(len(model.f)),
        "propagated_modes": int(len(indices)),
        "baseline_median_ms": statistics.median(base_times) * 1e3,
        "baseline_p95_ms": np.percentile(base_times, 95) * 1e3,
        "candidate_median_ms": statistics.median(candidate_times) * 1e3,
        "candidate_p95_ms": np.percentile(candidate_times, 95) * 1e3,
        "candidate_over_baseline": statistics.median(candidate_times)
        / statistics.median(base_times),
        "spectrum_dex": _stats(dex),
        "dn_base": float(dn_base),
        "dn_sparse": float(dn_sparse),
        "dn_relative": abs(dn_sparse - dn_base) / max(abs(dn_base), 1e-300),
        "baseline_digest": _digest(base),
        "candidate_digest": _digest(sparse),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="default,lowT,highT,stiff,high_kappa")
    parser.add_argument("--repeats", type=int, default=25)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--out", default="docs/sparse_frequency_round42_20260923.json")
    args = parser.parse_args()
    rows = []
    for name in args.cases.split(","):
        row = run_case(name, args.repeats, args.threads, args.stride)
        rows.append(row)
        print(json.dumps(row, ensure_ascii=False))
    payload = {
        "experiment": "round42_sparse_frequency_fixed_dn_screen",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "resources": {"numba_threads": args.threads, "blas_threads": 1, "workers": 1},
        "settings": vars(args),
        "records": rows,
    }
    Path(args.out).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
