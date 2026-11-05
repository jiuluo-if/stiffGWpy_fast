"""Screen small curvature-ranked midpoint budgets on fixed Fast backgrounds.

The driver replays each archived adaptive Fast case, evaluates only selected
actual midpoint frequencies, and compares their PCHIP maps with the retained
exact nested Prüfer midpoint maps. Timing is kernel-only and excludes the
background replay/preparation. It is a feasibility screen, not an error bound.
"""

from __future__ import annotations

import hashlib
import json
import statistics
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import scripts.probe_adaptive_embedded_integral as probe  # noqa: E402
from stiffgwpy_fast import freq_adaptive  # noqa: E402

NESTED_PATH = ROOT / "docs" / "adaptive_nested_reference_midpoints_20261105.json"
OUTPUT = ROOT / "docs" / "adaptive_panel_fraction_budget_20261105.json"
FRACTIONS = (0.05, 0.10, 0.15)
REPEATS = 15
END_TO_END_REPEATS = 10


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)


def _run_end_to_end(params: dict, z_tail: float, target: float,
                    sigma_exact: bool, base_grid: np.ndarray,
                    eval_freqs: np.ndarray | None) -> tuple[float, float, int, float]:
    model = probe.LCDM_SG(**params)
    original_builder = freq_adaptive.grid_independent_freqs
    seed_grid = base_grid.copy()
    if eval_freqs is not None and eval_freqs.size:
        seed_grid = np.sort(np.concatenate((seed_grid, eval_freqs)))[::-1].copy()
    freq_adaptive.grid_independent_freqs = lambda _model, _freq_res=1.0: (
        seed_grid.copy(), float(np.min(seed_grid)), float(np.max(seed_grid))
    )
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=z_tail,
            freq_res=1.0,
            threads=2,
            accuracy_mode="fast",
            physical_guard=False,
            freq_grid="adaptive",
            freq_grid_target=target,
            freq_grid_max_points=int(seed_grid.size),
            eval_freqs=None,
            sigma_exact=sigma_exact,
        )
    finally:
        freq_adaptive.grid_independent_freqs = original_builder
    elapsed = time.perf_counter() - started
    if result is not model:
        raise RuntimeError(f"end-to-end Fast replay failed: {model.fast_failure_reason}")
    dn_gw = float(np.asarray(model.DN_gw, dtype=np.float64).flat[-1])
    return elapsed, float(model.cosmo_param["DN_eff"]), int(np.asarray(model.f).size), dn_gw


def _case(case: dict) -> dict:
    source_path = ROOT / "docs" / case["source_artifact"]
    source = json.loads(source_path.read_text(encoding="utf-8"))
    probe_path = ROOT / "docs" / case["probe_artifact"]
    midpoint_data = json.loads(probe_path.read_text(encoding="utf-8"))[
        "actual_fast_midpoint_samples"
    ]
    params = source.get("parameters", source.get("parameters_input"))
    dn_eff = float(source.get("fixed_DN_eff", source["fast"]["DN_eff_final"]))
    z_tail = float(source["reference"]["z_tail"])
    base_freqs = np.asarray(source["fast"]["frequencies_log10_hz"], dtype=np.float64)
    base_spectrum_archived = np.asarray(source["fast"]["omega_gw"], dtype=np.float64)
    mids = np.asarray(midpoint_data["frequencies_log10_hz"], dtype=np.float64)
    mids_order = np.argsort(mids)
    mids = mids[mids_order]
    nested = case["nested_reference_spectrum"]
    nested_freqs = np.asarray(nested["frequencies_log10_hz"], dtype=np.float64)
    nested_ref = np.asarray(nested["omega_gw"], dtype=np.float64)
    base_x = np.sort(base_freqs)
    base_order = np.argsort(base_freqs)
    base_spectrum_archived = base_spectrum_archived[base_order]

    probe.SOURCE_Z_TAIL = z_tail
    probe.FS.apply_accuracy_mode("fast")
    probe.FS.set_threads(2)
    model = probe.LCDM_SG(**params)
    started = time.perf_counter()
    fast_result = model.SGWB_iter(
        engine="fast",
        fallback=False,
        tol=1e-6,
        z_tail=z_tail,
        freq_res=1.0,
        threads=2,
        accuracy_mode="fast",
        physical_guard=False,
        freq_grid="adaptive",
        freq_grid_target=float(source["fast"].get("target_dex", 1e-4)),
        freq_grid_max_points=int(source["fast"].get("max_points", 1500)),
        sigma_exact=bool(source["fast"].get("sigma_exact", True)),
    )
    replay_seconds = time.perf_counter() - started
    if fast_result is not model:
        raise RuntimeError(f"Fast replay failed for {source_path.name}: {model.fast_failure_reason}")
    if not np.array_equal(np.asarray(model.f), np.sort(base_freqs)[::-1]):
        raise RuntimeError(f"adaptive grid parity failure for {source_path.name}")
    if not np.isclose(float(model.cosmo_param["DN_eff"]), dn_eff, rtol=0.0, atol=2e-12):
        raise RuntimeError(f"fixed-point parity failure for {source_path.name}")

    base_args, _ = probe._prepared(model, base_freqs, dn_eff)
    base_outputs = [base_args[16], base_args[17], base_args[18]]
    base_replay = probe._evaluate(base_args, base_outputs)
    base_replay_sorted = base_replay[np.argsort(base_freqs)]
    base_parity_dex = float(np.max(np.abs(np.log10(base_replay_sorted / base_spectrum_archived))))
    if base_parity_dex > 1e-3:
        raise RuntimeError(
            f"archived spectrum replay mismatch for {source_path.name}: {base_parity_dex:.6g} dex"
        )
    curvature = np.abs(
        np.asarray(PchipInterpolator(base_x, np.log10(base_spectrum_archived)).derivative(2)(base_x))
    )
    widths = np.diff(base_x)
    panel_weights = 0.5 * (base_spectrum_archived[:-1] + base_spectrum_archived[1:]) * widths
    scores = curvature[:-1] * widths * widths / 8.0 * panel_weights
    ranking = np.argsort(scores)[::-1]

    ref_lookup = {float(freq): float(value) for freq, value in zip(nested_freqs, nested_ref)}
    full_ref_integral = float(PchipInterpolator(nested_freqs, nested_ref).integrate(nested_freqs[0], nested_freqs[-1]))
    full_ref_dn = float(case["pchip_DN_gw"]["nested_reference"])
    results = []
    for fraction in FRACTIONS:
        count = max(1, int(np.ceil(fraction * mids.size)))
        selected_indices = np.sort(ranking[:count])
        selected_freqs = np.sort(np.concatenate((base_x, mids[selected_indices])))
        selected_args, _ = probe._prepared(model, selected_freqs[::-1].copy(), dn_eff)
        outputs = [selected_args[16], selected_args[17], selected_args[18]]
        selected_fast = probe._evaluate(selected_args, outputs)[::-1]
        dn_fast = probe._dn_gw(selected_freqs, selected_fast, model)
        times, paired_base_times = [], []
        for repeat in range(REPEATS):
            if repeat % 2 == 0:
                started = time.perf_counter()
                probe._evaluate(base_args, base_outputs)
                paired_base_times.append(time.perf_counter() - started)
                started = time.perf_counter()
                probe._evaluate(selected_args, outputs)
                times.append(time.perf_counter() - started)
            else:
                started = time.perf_counter()
                probe._evaluate(selected_args, outputs)
                times.append(time.perf_counter() - started)
                started = time.perf_counter()
                probe._evaluate(base_args, base_outputs)
                paired_base_times.append(time.perf_counter() - started)
        paired_base_median = statistics.median(paired_base_times)
        selected_median = statistics.median(times)
        ref_values = np.asarray([ref_lookup[float(freq)] for freq in selected_freqs])
        dn_ref_selected = full_ref_dn * float(
            PchipInterpolator(selected_freqs, ref_values).integrate(selected_freqs[0], selected_freqs[-1])
        ) / full_ref_integral
        ref_on_fast_nodes = np.asarray([ref_lookup[float(freq)] for freq in selected_freqs])
        max_spectrum_error = float(np.max(np.abs(np.log10(selected_fast / ref_on_fast_nodes))))
        results.append({
            "selected_fraction": fraction,
            "selected_midpoint_count": int(count),
            "frequency_count": int(selected_freqs.size),
            "kernel_median_s": selected_median,
            "kernel_p95_s": float(np.percentile(times, 95)),
            "kernel_median_ratio_vs_base": selected_median / paired_base_median,
            "fast_DN_gw": dn_fast,
            "nested_reference_selected_DN_gw": dn_ref_selected,
            "nested_reference_full_DN_gw": full_ref_dn,
            "reference_selected_vs_full_relative": _rel(dn_ref_selected, full_ref_dn),
            "fast_selected_vs_full_reference_relative": _rel(dn_fast, full_ref_dn),
            "max_spectrum_abs_dex_vs_exact_nested_reference": max_spectrum_error,
        })
    # Measure the public Fast path with a held-fixed archived coarse grid. The
    # only candidate difference is the same curvature-ranked native midpoint
    # set; this is an end-to-end diagnostic and does not alter production code.
    selected_count = max(1, int(np.ceil(0.10 * mids.size)))
    selected_indices = np.sort(ranking[:selected_count])
    selected_eval_freqs = mids[selected_indices]
    baseline_times, candidate_times = [], []
    baseline_roots, candidate_roots = [], []
    baseline_dn, candidate_dn = [], []
    baseline_counts, candidate_counts = [], []
    for repeat in range(END_TO_END_REPEATS):
        calls = (
            ((None, baseline_times, baseline_roots, baseline_counts, baseline_dn),
             (selected_eval_freqs, candidate_times, candidate_roots, candidate_counts, candidate_dn))
            if repeat % 2 == 0 else
            ((selected_eval_freqs, candidate_times, candidate_roots, candidate_counts, candidate_dn),
             (None, baseline_times, baseline_roots, baseline_counts, baseline_dn))
        )
        for eval_grid, times, roots, counts, dn_values in calls:
            elapsed, root, count, dn_value = _run_end_to_end(
                params,
                z_tail,
                0.999,
                bool(source["fast"].get("sigma_exact", True)),
                base_freqs,
                eval_grid,
            )
            times.append(elapsed)
            roots.append(root)
            counts.append(count)
            dn_values.append(dn_value)
    return {
        "probe_artifact": probe_path.name,
        "source_artifact": source_path.name,
        "source_sha256": _sha(source_path),
        "fixed_DN_eff": dn_eff,
        "base_frequency_count": int(base_x.size),
        "outer_fast_replay_s": replay_seconds,
        "base_spectrum_parity_max_abs_dex": base_parity_dex,
        "candidates": results,
        "end_to_end_10pct_midpoint_candidate": {
            "selected_midpoint_count": int(selected_count),
            "baseline_frequency_counts": baseline_counts,
            "candidate_frequency_counts": candidate_counts,
            "baseline_DN_eff_final": baseline_roots,
            "candidate_DN_eff_final": candidate_roots,
            "baseline_DN_gw": baseline_dn,
            "candidate_DN_gw": candidate_dn,
            "baseline_seconds": baseline_times,
            "candidate_seconds": candidate_times,
            "baseline_median_s": statistics.median(baseline_times),
            "candidate_median_s": statistics.median(candidate_times),
            "candidate_to_baseline_median_ratio": statistics.median(candidate_times) / statistics.median(baseline_times),
            "candidate_median_p95_s": float(np.percentile(candidate_times, 95)),
            "baseline_median_p95_s": float(np.percentile(baseline_times, 95)),
            "root_relative_difference_median": _rel(statistics.median(candidate_roots), statistics.median(baseline_roots)),
            "dn_gw_relative_difference_median": _rel(statistics.median(candidate_dn), statistics.median(baseline_dn)),
            "fixed_grid_diagnostic": True,
        },
    }


def main() -> None:
    nested = json.loads(NESTED_PATH.read_text(encoding="utf-8"))
    cases = [_case(case) for case in nested["cases"]]
    record = {
        "schema_version": 1,
        "experiment": "curvature_ranked_small_panel_fraction_screen",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": probe._git("rev-parse", "HEAD"),
        "nested_reference_artifact": NESTED_PATH.name,
        "nested_reference_sha256": _sha(NESTED_PATH),
        "fractions": list(FRACTIONS),
        "timing_repeats": REPEATS,
        "cases": cases,
        "script_sha256": _sha(Path(__file__)),
        "limitations": [
            "Fixed-background results only; no self-consistent root certification.",
            "Kernel-only timing excludes Fast replay, background preparation, and full adaptive-loop orchestration.",
            "Curvature-ranked fractions are feasibility diagnostics, not a calibrated error estimator or an online policy.",
            "Exact nested Prüfer interpolation comparisons do not constitute a rigorous error bound.",
        ],
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": OUTPUT.name,
        "cases": [
            {
                "probe": case["probe_artifact"],
                "candidates": [
                    {
                        "fraction": item["selected_fraction"],
                        "kernel_ratio": item["kernel_median_ratio_vs_base"],
                        "Fast_to_nested_DN": item["fast_selected_vs_full_reference_relative"],
                        "max_spectrum_dex": item["max_spectrum_abs_dex_vs_exact_nested_reference"],
                    }
                    for item in case["candidates"]
                ],
            }
            for case in cases
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
