"""Screen phase-aware sparse-frequency reconstruction at fixed ``DN_eff``.

Round 42 interpolated the already squared spectrum and failed badly.  This
prototype instead interpolates the two real Cartesian handoff components and
forms the physical tail power only after reconstruction.  It remains a
standalone feasibility screen; production output columns and outer semantics
are not changed.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

try:
    from scripts._resource_budget import apply_environment
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment

apply_environment()

import numpy as np  # noqa: E402
from numba import njit, prange  # noqa: E402
from scipy.interpolate import PchipInterpolator  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.benchmark_phase_recurrence import _make_args, _prepared  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402


def _mode_indices(size, stride=2):
    indices = np.arange(0, int(size), int(stride), dtype=np.int64)
    if indices[-1] != int(size) - 1:
        indices = np.append(indices, int(size) - 1)
    return indices


def _interp(full_x, sparse_x, values):
    order = np.argsort(sparse_x, kind="stable")
    x = np.asarray(sparse_x, dtype=np.float64)[order]
    y = np.asarray(values, dtype=np.float64)[order]
    result = np.asarray(
        PchipInterpolator(x, y, extrapolate=False)(
            np.asarray(full_x, dtype=np.float64)
        ),
        dtype=np.float64,
    )
    for index, frequency in enumerate(x):
        result[np.asarray(full_x) == frequency] = y[index]
    return result


def _complex_reconstruct(
    full_x, sparse_x, x_values, y_values, pref_values, last_z_values,
    tail_scale, gamma=0.0,
):
    """Interpolate Cartesian state, then form a tail power observable."""
    x_full = _interp(full_x, sparse_x, x_values)
    y_full = _interp(full_x, sparse_x, y_values)
    pref_full = _interp(full_x, sparse_x, pref_values)
    last_z_full = _interp(full_x, sparse_x, last_z_values)
    e_z = np.exp(-last_z_full)
    amp2 = (
        x_full * x_full
        + y_full * y_full * (1.0 + gamma * gamma * e_z * e_z)
        + 2.0 * gamma * x_full * y_full * e_z
    )
    result = np.asarray(tail_scale, dtype=np.float64) * pref_full * pref_full * amp2
    return result


@njit(parallel=True, cache=False)
def _handoff_state(
    Nv, Phi_grid, Phi_mid, S2inv, j0s, z0s, h, z_tail, Sv, phase_max,
    kink_index, kink_fraction, phi_re, x_out, y_out, z_out, kend_out,
    x_final_out, y_final_out,
):
    """Propagate sparse modes and retain the production tail handoff state."""
    _ = Sv
    nv = len(Nv)
    for mode in prange(len(j0s)):
        j0 = j0s[mode]
        z0 = z0s[mode]
        phi0 = Phi_grid[j0]
        x_value = 0.0
        y_value = math.exp(z0) * S2inv[j0]
        k = j0
        z_value = z0 + Phi_grid[k] - phi0
        last_x = 0.0
        last_y = y_value
        last_z = z_value
        if z_value < z_tail:
            while k < nv - 1 and z_value < z_tail:
                z_mid = z0 + Phi_mid[k] - phi0
                if k == kink_index and 0.0 < kink_fraction < 1.0:
                    z_break = z0 + phi_re - phi0
                    z_end = z0 + Phi_grid[k + 1] - phi0
                    h_left = h * kink_fraction
                    h_right = h - h_left
                    x_value, y_value = FS._phase_segment(
                        x_value, y_value, z_value, z_break, h_left, phase_max
                    )
                    x_value, y_value = FS._phase_segment(
                        x_value, y_value, z_break, z_end, h_right, phase_max
                    )
                else:
                    z_end = 2.0 * z_mid - z_value
                    x_value, y_value = FS._phase_segment(
                        x_value, y_value, z_value, z_end, h, phase_max
                    )
                k += 1
                z_value = z0 + Phi_grid[k] - phi0
                if z_value < z_tail:
                    last_x = x_value
                    last_y = y_value
                    last_z = z_value
            if z_value < z_tail:
                kend = nv - 1
            else:
                kend = k - 1
                if kend < j0:
                    kend = j0
        else:
            kend = j0
        x_out[mode] = last_x
        y_out[mode] = last_y
        z_out[mode] = last_z
        kend_out[mode] = kend
        x_final_out[mode] = x_value
        y_final_out[mode] = y_value


def _stats(values):
    values = np.asarray(values, dtype=np.float64)
    return {
        "median": float(np.median(values)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def run_case(name, threads=2, stride=2):
    model, common = _prepared(name, threads)
    baseline = _make_args(common)
    FS.solve_kernel(*baseline)
    base = baseline[16][:, -1] - baseline[17][:, -1]
    (
        Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus, fp_minus,
        fp_freq, _, _, _, h, z_tail, _, Sv, phase_max, kink_index,
        kink_fraction, phi_re,
    ) = common
    indices = _mode_indices(len(model.f), stride)
    x_values = np.empty(len(indices), dtype=np.float64)
    y_values = np.empty(len(indices), dtype=np.float64)
    x_final = np.empty(len(indices), dtype=np.float64)
    y_final = np.empty(len(indices), dtype=np.float64)
    z_values = np.empty(len(indices), dtype=np.float64)
    kend = np.empty(len(indices), dtype=np.int64)
    _handoff_state(
        Nv, Phi_grid, Phi_mid, S2inv, j0s[indices], z0s[indices], h,
        z_tail, Sv, phase_max, kink_index, kink_fraction, phi_re,
        x_values, y_values, z_values, kend, x_final, y_final,
    )
    pref = np.empty(len(indices), dtype=np.float64)
    for row, k_end in enumerate(kend):
        pref[row] = math.sqrt(0.5 * S2[k_end]) * math.exp(Nv[k_end] - z_values[row])
    tail_scale = (
        np.asarray(fp_freq, dtype=np.float64) ** 2
        * np.asarray(P_t, dtype=np.float64)
        * (ev_minus[-1] * fp_minus[-1]) ** 2 / 12.0
    )
    candidate = _complex_reconstruct(
        np.asarray(model.f), np.asarray(model.f)[indices], x_values, y_values,
        pref, z_values, tail_scale, gamma=1.0,
    )
    for row, k_end in enumerate(kend):
        if k_end == len(Nv) - 1:
            candidate[indices[row]] = (
                S2[-1] * (x_final[row] * x_final[row] + y_final[row] * y_final[row])
                / 24.0 * P_t[indices[row]]
            )
    positive_base = np.maximum(base, 1e-300)
    positive_candidate = np.maximum(candidate, 1e-300)
    dex = np.abs(np.log10(positive_candidate) - np.log10(positive_base))
    dn_base = FS.integrate_frequency_pchip(model.f, base)
    dn_candidate = FS.integrate_frequency_pchip(model.f, candidate)
    return {
        "case": name,
        "stride": int(stride),
        "full_modes": int(len(model.f)),
        "propagated_modes": int(len(indices)),
        "kend_nv_minus_one": int(np.sum(kend == len(Nv) - 1)),
        "spectrum_dex": _stats(dex),
        "dn_base": float(dn_base),
        "dn_candidate": float(dn_candidate),
        "dn_relative": abs(dn_candidate - dn_base) / max(abs(dn_base), 1e-300),
        "finite": bool(np.isfinite(candidate).all()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", default="default,lowT,highT,stiff,high_kappa")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--out", default="docs/phase_aware_sparse_frequency_round43_20260923.json")
    args = parser.parse_args()
    rows = [run_case(name, args.threads, args.stride) for name in args.cases.split(",")]
    payload = {
        "experiment": "round43_phase_aware_sparse_frequency_fixed_dn_screen",
        "production_path_changed": False,
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "resources": {"numba_threads": args.threads, "blas_threads": 1, "workers": 1},
        "settings": vars(args),
        "records": rows,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
