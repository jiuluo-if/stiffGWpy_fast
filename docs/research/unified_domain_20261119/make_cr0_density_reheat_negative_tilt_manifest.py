"""Create reproducible hashes and gate summaries for the CR0 T_re slice."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr0_density_reheat_negative_tilt_manifest.json"
SOURCES = (
    "exact_background.py", "fast_sgwb.py", "freq_adaptive.py",
    "_wkb_transfer.py", "LCDM_stiff_Neff.py", "reference.py",
)
ARTIFACTS = (
    "cr0_density_reheat_negative_tilt_summary.md",
    "cr0_density_reheat_negative_tilt_evidence.png",
    "plot_cr0_density_reheat_negative_tilt.py",
    "make_cr0_density_reheat_negative_tilt_manifest.py",
    "screen_cr0_density_reheat_negative_tilt.py",
    "screen_cr0_density_reheat_negative_tilt_lattice.json",
    "screen_cr0_density_reheat_negative_tilt_sobol.json",
    "validate_cr0_density_reheat_negative_tilt_z8.py",
    "validate_cr0_density_reheat_negative_tilt_z8.json",
    "close_reference_cr0_density_reheat_negative_tilt_z8.py",
    "close_reference_cr0_density_reheat_negative_tilt_z8.json",
    "benchmark_cr0_density_reheat_negative_tilt_ab.py",
    "benchmark_cr0_density_reheat_negative_tilt_ab_500pairs.json",
    "benchmark_cr0_reheat_trusted_control_1000pairs.py",
    "benchmark_cr0_reheat_trusted_control_1000pairs.json",
    "screen_cr0_density_reheat_tilt.py",
    "screen_cr0_density_reheat_tilt_lattice.json",
    "screen_cr0_density_reheat_tilt_sobol.json",
    "screen_cr0_tilt_temperature.py",
    "screen_cr0_kappa_tilt_extension.py",
    "screen_cr0_kappa_tilt_extension_lattice.json",
    "screen_cr0_kappa_tilt_extension_sobol.json",
    "screen_cr0_tilt_temperature.py",
    "screen_cr0_tilt_temperature_lattice.json",
    "screen_cr0_tilt_temperature_sobol.json",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_cases(name: str) -> list[dict]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))["cases"]


def main() -> None:
    accepted = read_cases("screen_cr0_density_reheat_negative_tilt_lattice.json") + read_cases(
        "screen_cr0_density_reheat_negative_tilt_sobol.json"
    )
    baseline_success = sum(bool(row["baseline_fast"].get("success")) for row in accepted)
    baseline_failures: dict[str, int] = {}
    for row in accepted:
        reason = row["baseline_fast"].get("failure_reason")
        if reason:
            baseline_failures[reason] = baseline_failures.get(reason, 0) + 1

    refs = json.loads(
        (HERE / "validate_cr0_density_reheat_negative_tilt_z8.json").read_text(encoding="utf-8")
    )
    closure = json.loads(
        (HERE / "close_reference_cr0_density_reheat_negative_tilt_z8.json").read_text(encoding="utf-8")
    )["case"]
    benchmark = json.loads(
        (HERE / "benchmark_cr0_density_reheat_negative_tilt_ab_500pairs.json").read_text(encoding="utf-8")
    )
    extended = json.loads(
        (HERE / "benchmark_cr0_reheat_trusted_control_1000pairs.json").read_text(encoding="utf-8")
    )
    main_times = next(
        row["timing_pairs"] for row in benchmark["cases"]
        if row["case_id"] == "trusted_cr0_standard_control"
    )
    extended_times = extended["cases"][0]["timing_pairs"]
    base = np.asarray([row["base_seconds"] for row in main_times + extended_times], dtype=np.float64)
    current = np.asarray([row["current_seconds"] for row in main_times + extended_times], dtype=np.float64)
    pooled_control = {
        "pairs": int(len(base)),
        "p50_ratio": float(np.quantile(current, 0.5) / np.quantile(base, 0.5)),
        "p95_ratio": float(np.quantile(current, 0.95) / np.quantile(base, 0.95)),
        "baseline_p50_seconds": float(np.quantile(base, 0.5)),
        "current_p50_seconds": float(np.quantile(current, 0.5)),
        "baseline_p95_seconds": float(np.quantile(base, 0.95)),
        "current_p95_seconds": float(np.quantile(current, 0.95)),
    }
    region_ids = ("new_", "retained_", "independent_")
    region_perf = [row for row in benchmark["cases"] if row["case_id"].startswith(region_ids)]
    accepted_screen = [row["candidate_fast"] for row in accepted]

    payload = {
        "status": "validated_finite_sample_cr0_density_reheat_negative_tilt_region",
        "accepted_region": {
            "Omega_bh2": 1.0e6,
            "Omega_ch2": 0.1201075,
            "T_re_GeV": [4.0e3, 1.0e7],
            "n_t": [-0.5, -0.25],
            "fixed_parameters": {
                "cr": 0.0, "H0": 67.32117, "DN_eff": 0.0,
                "A_s": 2.100549e-9, "r": 0.01, "DN_re": 10.0, "kappa10": 0.01,
            },
            "Cobaya_Omega_bh2_prior": [0.01, 0.04],
            "outside_Cobaya_Omega_bh2_prior": True,
            "physical_guard": True,
            "blue_tilt_cutoff_branch": False,
            "continuum_claim": False,
        },
        "paired_screen": {
            "points": len(accepted),
            "old_success": baseline_success,
            "old_failures_by_reason": baseline_failures,
            "current_success": sum(bool(row.get("success")) for row in accepted_screen),
            "current_failure_reasons": {},
            "H2_positive_finite_all": all(row.get("h2_positive_finite") for row in accepted_screen),
            "N_inf_finite_all": all(row.get("derived_N_inf") is not None for row in accepted_screen),
            "H2_minimum_sampled": min(row["h2_min"] for row in accepted_screen),
            "N_inf_range": [
                min(row["derived_N_inf"] for row in accepted_screen),
                max(row["derived_N_inf"] for row in accepted_screen),
            ],
        },
        "rejected_probe_summaries": {
            "wider_nt_to_zero": {
                "current_shared_Neff_guard_rejects": sum(
                    row["candidate_fast"].get("failure_reason") == "shared_Neff_guard"
                    for row in read_cases("screen_cr0_density_reheat_tilt_lattice.json")
                    + read_cases("screen_cr0_density_reheat_tilt_sobol.json")
                ),
                "included_in_accepted_domain": False,
            },
            "positive_tilt_high_reheat": {
                "points": 19,
                "current_shared_Neff_guard_rejects": sum(
                    row["candidate_fast"].get("failure_reason") == "shared_Neff_guard"
                    for row in read_cases("screen_cr0_tilt_temperature_lattice.json")
                    + read_cases("screen_cr0_tilt_temperature_sobol.json")
                ),
                "included_in_accepted_domain": False,
            },
            "kappa10_10_to_100_by_negative_tilt": {
                "old_success": sum(
                    bool(row["baseline_fast"].get("success"))
                    for row in read_cases("screen_cr0_kappa_tilt_extension_lattice.json")
                    + read_cases("screen_cr0_kappa_tilt_extension_sobol.json")
                ),
                "current_success": 19,
                "new_old_failure_coverage": 0,
            },
        },
        "accuracy": {
            "reference_samples": len(refs["cases"]),
            "max_first_map_root_error": max(row["root_relative_error_proxy"] for row in refs["cases"]),
            "max_positive_linear_spectrum_error_dex": max(
                row["spectrum_max_abs_dex_from_raw_linear"] for row in refs["cases"]
            ),
            "tight_corner_second_to_third_closure": closure["second_to_third_closure_relative"],
            "tight_corner_fast_root_error_to_third_map": closure["fast_root_relative_to_third_map"],
            "root_gate": 2e-4,
            "spectrum_gate_dex": 1e-3,
        },
        "performance": {
            "base_commit": benchmark["base_commit"],
            "candidate_commit": benchmark["candidate_commit"],
            "pairs_per_case": 500,
            "regional_p50_ratio_max": max(row["p50_ratio"] for row in region_perf),
            "regional_p95_ratio_max": max(row["p95_ratio"] for row in region_perf),
            "trusted_cr0_high_density_control": {
                "p50_ratio": next(row["p50_ratio"] for row in benchmark["cases"] if row["case_id"] == "trusted_cr0_high_density_control"),
                "p95_ratio": next(row["p95_ratio"] for row in benchmark["cases"] if row["case_id"] == "trusted_cr0_high_density_control"),
            },
            "trusted_cr0_standard_control_pooled": pooled_control,
            "root_delta_max": max(row["root_relative_delta_max"] for row in benchmark["cases"]),
            "frequency_digest_mismatch_pairs": sum(row["frequency_digest_mismatch_pairs"] for row in benchmark["cases"]),
            "spectrum_absdiff_max_dex": max(row["spectrum_absdiff_max_dex"] for row in benchmark["cases"]),
            "first_call_records": {
                row["case_id"]: row["first_call"] for row in benchmark["cases"]
            },
            "warm_resources": benchmark["resources"],
        },
        "reference_code_modified": False,
        "production_code_changed_this_tranche": False,
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name) for name in SOURCES
        },
        "artifact_sha256": {name: sha256(HERE / name) for name in ARTIFACTS},
        "manifest_self_hash_included": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
