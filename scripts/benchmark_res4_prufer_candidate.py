"""Screen the production Prüfer Reference backend against an archived res4 map.

This reuses a same-input Cartesian Reference spectrum; it does not rerun the
expensive Cartesian propagation. Timings cover only candidate propagation.
"""

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import telemetry  # noqa: E402
from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as REF  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PARAMS = dict(r=0.01, cr=0.0, n_t=2.0, T_re=1e4, kappa10=0.01,
              DN_re=10.0)
DEFAULT_INDICES = [0, 65, 100, 130, 160, 195, 220, 259]


def _log10_output_omega(value):
    value = float(value)
    if not math.isfinite(value):
        raise FloatingPointError("non-finite OmegaGW output in Reference benchmark")
    return math.log10(max(value, 1e-40))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "docs" / "res4_sigma_exact_reference_map_20261104.npz")
    parser.add_argument("--dn-eff", type=float, default=0.133732083004284)
    parser.add_argument("--z-tail", type=float, default=10.0)
    parser.add_argument("--rtol", type=float, default=1e-11)
    parser.add_argument("--cartesian-rtol", type=float, default=None,
                        help="optional Cartesian baseline tolerance; defaults to --rtol")
    parser.add_argument("--indices", nargs="+", type=int,
                        default=DEFAULT_INDICES)
    parser.add_argument("--repeats", type=int, default=1,
                        help="alternating A/B repetitions per selected mode")
    parser.add_argument("--candidate-only", action="store_true",
                        help="skip Cartesian reruns and compare against the archive only")
    parser.add_argument("--full-grid", action="store_true",
                        help="evaluate all archived modes and compute the candidate map")
    parser.add_argument("--reference-map", type=Path,
                        default=ROOT / "docs" / "res4_sigma_exact_reference_map_20261104.json")
    parser.add_argument("--skip-map-crosscheck", action="store_true",
                        help="do not compare DN_gw to the default map at a different DN_eff")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs" / "reference_prufer_res4_screen_20261104.json")
    args = parser.parse_args(argv)

    with np.load(args.reference) as archive:
        freqs = np.asarray(archive["freqs"], dtype=float)
        reference_log = np.asarray(archive["reference_log"], dtype=float)
    if freqs.shape != reference_log.shape:
        parser.error("archived Reference frequencies and spectrum must have equal shape")
    indices = list(range(len(freqs))) if args.full_grid else args.indices
    if any(index < 0 or index >= len(freqs) for index in indices):
        parser.error("every mode index must be inside the archived spectrum")
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.full_grid and args.repeats != 1:
        parser.error("--full-grid currently requires exactly one run per mode")

    model = LCDM_SG(**PARAMS, DN_eff=args.dn_eff)
    reference_source = ROOT / "stiffgwpy_fast" / "reference.py"
    script_digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    reference_digest = hashlib.sha256(reference_source.read_bytes()).hexdigest()
    jit_warmup_start = time.perf_counter()
    REF._warmup_numba_spline()
    jit_warmup_seconds = time.perf_counter() - jit_warmup_start
    cartesian_rtol = (args.rtol if args.cartesian_rtol is None
                      else args.cartesian_rtol)
    rows = []
    candidate_ogw = []
    candidate_oj = []
    all_candidate_times = []
    all_cartesian_times = []
    all_speedups = []
    for index in indices:
        freq = float(freqs[index])
        reference_value = float(reference_log[index])
        candidate_rows = []
        candidate_times = []
        cartesian_rows = []
        cartesian_times = []
        speedups = []
        cartesian_archive_errors = []
        for repeat in range(args.repeats):
            cartesian = None
            cartesian_s = None
            if args.candidate_only:
                start = time.perf_counter()
                row = REF.solve_reference_mode(
                    model, freq, args.dn_eff, z_tail=args.z_tail,
                    rtol=args.rtol, backend="prufer")
                candidate_s = time.perf_counter() - start
            elif (index + repeat) % 2:
                start = time.perf_counter()
                cartesian = REF.solve_reference_mode(
                    model, freq, args.dn_eff, z_tail=args.z_tail,
                    rtol=cartesian_rtol)
                cartesian_s = time.perf_counter() - start
                start = time.perf_counter()
                row = REF.solve_reference_mode(
                    model, freq, args.dn_eff, z_tail=args.z_tail,
                    rtol=args.rtol, backend="prufer")
                candidate_s = time.perf_counter() - start
            else:
                start = time.perf_counter()
                row = REF.solve_reference_mode(
                    model, freq, args.dn_eff, z_tail=args.z_tail,
                    rtol=args.rtol, backend="prufer")
                candidate_s = time.perf_counter() - start
                start = time.perf_counter()
                cartesian = REF.solve_reference_mode(
                    model, freq, args.dn_eff, z_tail=args.z_tail,
                    rtol=cartesian_rtol)
                cartesian_s = time.perf_counter() - start

            candidate_rows.append(row)
            candidate_times.append(candidate_s)
            all_candidate_times.append(candidate_s)
            if cartesian is not None:
                cartesian_rows.append(cartesian)
                cartesian_times.append(cartesian_s)
                all_cartesian_times.append(cartesian_s)
                speedup = cartesian_s / candidate_s
                speedups.append(speedup)
                all_speedups.append(speedup)
                cartesian_value = float(
                    cartesian["Ogw_today"] - cartesian["Oj_today"])
                cartesian_log = _log10_output_omega(cartesian_value)
                cartesian_archive_errors.append(
                    abs(cartesian_log - reference_value)
                    if math.isfinite(cartesian_log)
                    and math.isfinite(reference_value) else None)

        row = candidate_rows[-1]
        candidate_ogw.append(float(row["Ogw_today"]))
        candidate_oj.append(float(row["Oj_today"]))
        candidate_log_samples = [
            _log10_output_omega(item["Ogw_today"] - item["Oj_today"])
            for item in candidate_rows]
        candidate_log = candidate_log_samples[-1]
        errors_for_mode = [abs(value - reference_value)
                           for value in candidate_log_samples
                           if math.isfinite(value) and math.isfinite(reference_value)]
        cartesian_log = (
            _log10_output_omega(cartesian_rows[-1]["Ogw_today"]
                                - cartesian_rows[-1]["Oj_today"])
            if cartesian_rows else None)
        rows.append({
            "index": index,
            "freq_log10_hz": freq,
            "candidate_runtime_s": float(np.median(candidate_times)),
            "candidate_runtime_samples_s": candidate_times,
            "candidate_nfev": int(row["n_steps"]),
            "candidate_nfev_samples": [int(item["n_steps"])
                                       for item in candidate_rows],
            "used_tail": bool(row["used_tail"]),
            "cartesian_runtime_s": (float(np.median(cartesian_times))
                                    if cartesian_times else None),
            "cartesian_runtime_samples_s": cartesian_times,
            "cartesian_nfev": (int(cartesian_rows[-1]["n_steps"])
                               if cartesian_rows else None),
            "speedup": float(np.median(speedups)) if speedups else None,
            "speedup_samples": speedups,
            "cartesian_log10_omega": cartesian_log,
            "cartesian_abs_dex_from_archive": (
                max(cartesian_archive_errors) if cartesian_archive_errors else None),
            "reference_log10_omega": reference_value,
            "candidate_log10_omega": candidate_log,
            "candidate_log10_omega_samples": candidate_log_samples,
            "abs_dex_error": max(errors_for_mode) if errors_for_mode else None,
        })

    errors = [row["abs_dex_error"] for row in rows
              if row["abs_dex_error"] is not None]
    payload = {
        "schema_version": 1,
        "experiment": "production_prufer_backend_against_archived_cartesian_res4",
        "source_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_script_sha256": script_digest,
        "reference_module_sha256": reference_digest,
        "worktree_dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "parameters": {**PARAMS, "DN_eff": args.dn_eff},
        "z_tail": args.z_tail,
        "candidate_rtol": args.rtol,
        "cartesian_rtol": cartesian_rtol,
        "repeats_per_mode": args.repeats,
        "atol": [1e-12, 1e-12, 1e-12],
        "spline_jit_warmup_seconds": jit_warmup_seconds,
        "resources": telemetry(workers=1, threads=1),
        "reference_archive": str(args.reference.resolve()),
        "frequency_count_full_grid": int(len(freqs)),
        "selected_indices": indices,
        "max_abs_dex_error": max(errors) if errors else None,
        "median_abs_dex_error": float(np.median(errors)) if errors else None,
        "sum_candidate_seconds": sum(all_candidate_times),
        "sum_cartesian_seconds": (sum(all_cartesian_times)
                                  if all_cartesian_times else None),
        "aggregate_speedup": (
            sum(all_cartesian_times) / sum(all_candidate_times)
            if all_cartesian_times else None),
        "per_mode_speedup_p50": (float(np.percentile(all_speedups, 50))
                                 if all_speedups else None),
        "per_mode_speedup_p95": (float(np.percentile(all_speedups, 95))
                                 if all_speedups else None),
        "rows": rows,
        "semantics": "production Reference backend=prufer; compare to archived Cartesian same-input map",
    }
    if args.full_grid:
        g2, quadrature_error, interpolation_error = REF.integrate_spectrum(
            freqs, np.asarray(candidate_ogw), np.asarray(candidate_oj))
        candidate_dn_gw = (gp.Neff0 * g2
                           / (gp.Omega_nh2 / model.derived_param["h"] ** 2))
        reference_dn_gw = None
        if not args.skip_map_crosscheck:
            archived_map = json.loads(args.reference_map.read_text(encoding="utf-8"))
            reference_dn_gw = float(archived_map["reference_dn_map"])
        errors = [row["abs_dex_error"] for row in rows]
        payload.update({
            "candidate_DN_gw": float(candidate_dn_gw),
            "reference_DN_gw": reference_dn_gw,
            "DN_gw_relative_error_vs_reference": (
                abs(candidate_dn_gw - reference_dn_gw) / abs(reference_dn_gw)
                if reference_dn_gw is not None else None),
            "candidate_map_residual_relative": (
                candidate_dn_gw - args.dn_eff) / args.dn_eff,
            "candidate_quadrature_error": quadrature_error,
            "candidate_interpolation_error": interpolation_error,
            "full_grid_spectrum_error_max_dex": max(
                value for value in errors if value is not None),
            "full_grid_spectrum_error_p95_dex": float(np.percentile(
                [value for value in errors if value is not None], 95)),
        })
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
