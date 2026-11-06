"""Feasibility screen for exact discrete frequency tangents on low-T Fast."""

from __future__ import annotations

import hashlib
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402

from scripts import prototype_frequency_tangent as tangent  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.exact_background import (  # noqa: E402
    build_transition_grid,
    exact_phi_s2_split,
)
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PANEL = ROOT / "docs" / "cr1_lowT_native_grid_convergence_20261106.json"
REF_ROOT = ROOT / "docs" / "log_natural_cubic_cr1_corner_301_root_20261106.json"
OUTPUT = ROOT / "docs" / "frequency_tangent_lowT_20261106.json"
THREADS = 2
PAIRS = 25


def _sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _production_terminal(freqs, arrays):
    _, model, Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr, kink, fraction, phi_re = arrays
    nf = freqs.size
    j0, z0, fp_minus = FS.prep_frequency_only(model, Nv, freqs)[2:]
    fp_freq = np.power(10.0, freqs)
    P_t = model.derived_param["A_t"] * np.power(
        fp_freq / gp.f_piv, model.derived_param["nt"]
    )
    ev_minus = np.exp(-Nv)
    out = [np.zeros((nf, 2), dtype=np.float64) for _ in range(3)]
    FS.solve_kernel(
        Nv, Phi, Phi_mid, S2, S2inv, j0, z0, P_t, ev_minus,
        fp_minus, fp_freq, 1, 2, Nv.size - 1, float(h_arr[0]), 5.0,
        *out, h_arr, sigma, 0.25, np.full(nf, -1.0, dtype=np.float64),
        kink, fraction, phi_re, float(h_arr[-1]),
    )
    return out[0][:, -1] - out[1][:, -1]


def _candidate_terminal(freqs, arrays):
    _, model, Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr, kink, fraction, phi_re = arrays
    j0, z0, fp_minus = FS.prep_frequency_only(model, Nv, freqs)[2:]
    fp_freq = np.power(10.0, freqs)
    n_t = float(model.derived_param["nt"])
    P_t = model.derived_param["A_t"] * np.power(fp_freq / gp.f_piv, n_t)
    return tangent.propagate_frequency_grid(
        Nv, Phi, Phi_mid, S2, S2inv, sigma, h_arr,
        j0, z0, P_t, fp_freq, np.exp(-Nv), fp_minus,
        0.25, 5.0, kink, fraction, phi_re,
        np.full(freqs.size, n_t, dtype=np.float64),
    )


def _build_arrays(params, dn_eff, base_grid):
    model_params = dict(params)
    model_params["DN_eff"] = dn_eff
    model = LCDM_SG(**model_params)
    h = float(FS.resolve_config("fast", threads=THREADS).h)
    Nv, sigma, _, _ = build_transition_grid(model, h)
    Phi, Phi_mid, S2, S2inv, kink, fraction, phi_re = exact_phi_s2_split(
        model, Nv, dn_eff, sigma_nodes=sigma
    )
    return (
        base_grid, model, Nv, Phi, Phi_mid, S2, S2inv, sigma,
        np.diff(Nv), kink, fraction, phi_re,
    )


def _hermite_screen(base_freq, base_spectrum, base_slope, mid_spectrum):
    order = np.argsort(base_freq)
    u = np.log(np.power(10.0, base_freq[order]))
    y = base_spectrum[order]
    # `base_slope` is d ln(G)/d ln(f); Hermite needs dG/d ln(f).
    m = y * base_slope[order]
    h = np.diff(u)
    predicted_mid = 0.5 * (y[:-1] + y[1:]) + h * (m[:-1] - m[1:]) / 8.0
    panel_integrals = (
        0.5 * h * (y[:-1] + y[1:])
        + h * h * (m[:-1] - m[1:]) / 12.0
    )
    actual_mid = mid_spectrum[1::2][::-1]
    positive = (predicted_mid > 0.0) & (actual_mid > 0.0)
    error = np.full(predicted_mid.size, np.nan)
    error[positive] = np.abs(
        np.log10(predicted_mid[positive]) - np.log10(actual_mid[positive])
    )
    return predicted_mid, panel_integrals, error


def run() -> dict:
    panel = json.loads(PANEL.read_text(encoding="utf-8"))
    ref = json.loads(REF_ROOT.read_text(encoding="utf-8"))
    params = dict(panel["parameters"])
    dn_eff = float(ref["roots"]["reference_301"])
    start_model = LCDM_SG(**params)
    h = float(FS.resolve_config("fast", threads=THREADS).h)
    build_transition_grid(start_model, h)
    base_freq = np.asarray(FA.goal_oriented_freqs(start_model, freq_res=1.0), dtype=np.float64)
    nested_freq = np.sort(np.concatenate((base_freq, 0.5 * (base_freq[:-1] + base_freq[1:]))))[::-1]
    expected_base_hash = panel["grids"]["base_sha256"]
    if hashlib.sha256(base_freq.astype("<f8").tobytes()).hexdigest() != expected_base_hash:
        raise RuntimeError("reconstructed 76-node low-T frequency support differs from the archived screen")
    arrays = _build_arrays(params, dn_eff, base_freq)

    base_actual = _production_terminal(base_freq, arrays)
    base_tangent, base_slope = _candidate_terminal(base_freq, arrays)
    mid_actual = _production_terminal(nested_freq, arrays)
    _, _, midpoint_error = _hermite_screen(base_freq, base_tangent, base_slope, mid_actual)
    predicted_mid, panel_integrals, _ = _hermite_screen(base_freq, base_tangent, base_slope, mid_actual)

    base_relative = np.abs(base_tangent - base_actual) / np.maximum(np.abs(base_actual), 1e-300)
    base_log10_error = np.abs(np.log10(base_tangent) - np.log10(base_actual))
    base_map = FS.integrate_frequency_quadrature(base_freq, base_actual, "pchip") * np.log(10.0)
    nested_map = FS.integrate_frequency_quadrature(nested_freq, mid_actual, "pchip") * np.log(10.0)
    hermite_map = float(np.sum(panel_integrals))
    Omega_nu = gp.Omega_nh2 / arrays[1].derived_param["h"] ** 2
    coeff = gp.Neff0 / Omega_nu
    mapped = {
        "base_76_dn_gw_map": float(coeff * base_map),
        "hermite_76_to_151_dn_gw_map": float(coeff * hermite_map),
        "nested_151_dn_gw_map": float(coeff * nested_map),
        "reference_301_root": dn_eff,
    }
    deferred = []
    j0_base = FS.prep_frequency_only(arrays[1], arrays[2], base_freq)[2]
    j0_nested = FS.prep_frequency_only(arrays[1], arrays[2], nested_freq)[2]
    for i in range(base_freq.size - 1):
        deferred.append((int(j0_base[i]) == int(j0_base[i + 1]), int(j0_nested[2 * i + 1]) in (int(j0_nested[2 * i]), int(j0_nested[2 * i + 2]))))

    def baseline_call():
        _production_terminal(base_freq, arrays)

    def tangent_call():
        _candidate_terminal(base_freq, arrays)

    baseline_call()
    tangent_call()
    times = {"baseline": [], "tangent": []}
    for pair in range(PAIRS):
        order_names = ("baseline", "tangent") if pair % 2 == 0 else ("tangent", "baseline")
        for name in order_names:
            started = time.perf_counter()
            baseline_call() if name == "baseline" else tangent_call()
            times[name].append(time.perf_counter() - started)
    runtime = {
        name: {
            "p50_s": float(statistics.median(samples)),
            "p95_s": float(np.percentile(samples, 95)),
        }
        for name, samples in times.items()
    }
    runtime["tangent_over_baseline_p50"] = runtime["tangent"]["p50_s"] / runtime["baseline"]["p50_s"]
    runtime["tangent_over_baseline_p95"] = runtime["tangent"]["p95_s"] / runtime["baseline"]["p95_s"]
    positive_mid_count = int(np.count_nonzero(predicted_mid > 0.0))
    candidate_pass = (
        positive_mid_count == midpoint_error.size
        and float(np.nanmax(midpoint_error)) <= 1e-3
        and abs(hermite_map - nested_map) / abs(nested_map) <= 2e-4
        and runtime["tangent_over_baseline_p50"] <= 1.05
        and runtime["tangent_over_baseline_p95"] <= 1.05
    )

    record = {
        "schema_version": 1,
        "status": "complete",
        "experiment": "fixed-background low-T discrete frequency tangent versus nested Fast midpoint grid",
        "candidate_decision": "PASS_SCREEN" if candidate_pass else "REJECTED_AT_LOW_T_FEASIBILITY_SCREEN",
        "case": panel["case"],
        "parameters": {**params, "DN_eff": dn_eff},
        "base_frequency_count": int(base_freq.size),
        "nested_frequency_count": int(nested_freq.size),
        "base_frequency_sha256": hashlib.sha256(base_freq.astype("<f8").tobytes()).hexdigest(),
        "nested_frequency_sha256": hashlib.sha256(nested_freq.astype("<f8").tobytes()).hexdigest(),
        "tangent_vs_production_base_spectrum": {
            "max_relative": float(np.max(base_relative)),
            "max_abs_dex": float(np.max(base_log10_error)),
        },
        "hermite_prediction_of_nested_midpoints": {
            "positive_prediction_count": int(np.count_nonzero(predicted_mid > 0.0)),
            "nonpositive_prediction_count": int(np.count_nonzero(predicted_mid <= 0.0)),
            "max_abs_dex_on_positive_predictions": float(np.nanmax(midpoint_error)),
            "p95_abs_dex_on_positive_predictions": float(np.nanpercentile(midpoint_error, 95)),
            "midpoint_count": int(midpoint_error.size),
        },
        "map_comparison": {
            **mapped,
            "base_76_vs_nested_151_relative": abs(base_map - nested_map) / abs(nested_map),
            "hermite_76_vs_nested_151_relative": abs(hermite_map - nested_map) / abs(nested_map),
            "hermite_76_vs_reference_root_relative": abs(mapped["hermite_76_to_151_dn_gw_map"] - dn_eff) / abs(dn_eff),
            "nonpositive_hermite_panel_count": int(np.count_nonzero(panel_integrals <= 0.0)),
        },
        "branch_screen": {
            "base_adjacent_j0_changes": int(np.count_nonzero(np.diff(j0_base))),
            "midpoint_j0_matches_an_endpoint_count": int(sum(ok[1] for ok in deferred)),
            "interval_count": int(len(deferred)),
            "limitations": ["j0 changes are counted but not segmented in the Hermite prediction", "phase-substep ceil changes are not yet counted"],
        },
        "kernel_timing": {"pairs": PAIRS, "threads": THREADS, "alternating": True, **runtime},
        "resources": telemetry(workers=1, threads=THREADS),
        "source_hashes": {
            "fast_sgwb": _sha_file(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "freq_adaptive": _sha_file(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "tangent_prototype": _sha_file(Path(tangent.__file__)),
            "panel": _sha_file(PANEL),
            "reference_root": _sha_file(REF_ROOT),
        },
        "limitations": [
            "One fixed cr=1 low-temperature background; this is only a feasibility screen.",
            "The Hermite map uses exact derivatives of the implemented discrete step on each local fixed path but does not split all j0, substep, or tail branch boundaries.",
            "This is a kernel timing, not full public Fast end-to-end timing or acceptance evidence.",
            "No candidate self-consistent root or independent z_tail=5/deep Reference spectrum was solved here.",
        ],
    }
    return record


def main() -> int:
    FS.apply_accuracy_mode("fast")
    FS.set_threads(THREADS)
    record = run()
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "base_spectrum_parity": record["tangent_vs_production_base_spectrum"],
        "midpoint_prediction": record["hermite_prediction_of_nested_midpoints"],
        "map_comparison": record["map_comparison"],
        "timing": record["kernel_timing"],
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
