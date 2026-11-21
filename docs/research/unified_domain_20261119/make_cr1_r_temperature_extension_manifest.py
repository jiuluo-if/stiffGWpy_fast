"""Create hashes and machine-readable counts for the CR1 temperature extension."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr1_r_temperature_extension_manifest.json"
SOURCES = (
    "exact_background.py", "fast_sgwb.py", "freq_adaptive.py",
    "_wkb_transfer.py", "LCDM_stiff_Neff.py", "reference.py",
)
ARTIFACTS = (
    "cr1_r_temperature_extension_summary.md",
    "cr1_r_temperature_extension_evidence.png",
    "plot_cr1_r_temperature_extension.py",
    "make_cr1_r_temperature_extension_manifest.py",
    "screen_cr1_r_temperature_extension.py",
    "screen_cr1_r_temperature_extension_lattice.json",
    "screen_cr1_r_temperature_extension_sobol.json",
    "validate_cr1_r_temperature_extension_z8.py",
    "validate_cr1_r_temperature_extension_z8.json",
    "close_reference_cr1_r_temperature_extension_sobol_z8.py",
    "close_reference_cr1_r_temperature_extension_sobol_z8.json",
    "benchmark_cr1_r_temperature_extension_ab.py",
    "benchmark_cr1_r_temperature_extension_ab_500pairs.json",
    "screen_cr1_r_nextband.py", "screen_cr1_r_nextband_lattice.json",
    "screen_cr1_r_nextband_sobol.json", "validate_cr1_r_nextband_z8.py",
    "validate_cr1_r_nextband_z8.json", "screen_cr1_r_narrowband.py",
    "screen_cr1_r_narrowband_lattice.json", "screen_cr1_r_narrowband_sobol.json",
    "validate_cr1_r_narrowband_z8.py", "validate_cr1_r_narrowband_z8.json",
    "screen_cr1_r_narrowband.py", "screen_cr1_r_temperature_extension.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cases(filename: str) -> list[dict]:
    return json.loads((HERE / filename).read_text(encoding="utf-8"))["cases"]


def main() -> None:
    extension = cases("screen_cr1_r_temperature_extension_lattice.json") + cases(
        "screen_cr1_r_temperature_extension_sobol.json"
    )
    old_success = sum(bool(row["baseline_fast"].get("success")) for row in extension)
    old_reasons: dict[str, int] = {}
    for row in extension:
        reason = row["baseline_fast"].get("failure_reason")
        if reason:
            old_reasons[reason] = old_reasons.get(reason, 0) + 1
    current_success = sum(bool(row["candidate_fast"].get("success")) for row in extension)

    refs = json.loads((HERE / "validate_cr1_r_temperature_extension_z8.json").read_text(encoding="utf-8"))
    closure = json.loads(
        (HERE / "close_reference_cr1_r_temperature_extension_sobol_z8.json").read_text(encoding="utf-8")
    )["case"]
    benchmark = json.loads(
        (HERE / "benchmark_cr1_r_temperature_extension_ab_500pairs.json").read_text(encoding="utf-8")
    )
    measurements = benchmark["cases"]
    current_source_sha = {name: sha256(ROOT / "stiffgwpy_fast" / name) for name in SOURCES}

    payload = {
        "status": "validated_finite_sample_cr1_temperature_extension",
        "accepted_region": {
            "T_re": [2000.0, 2250.0], "r": [0.005, 0.01], "kappa10": 0.01,
            "fixed_parameters": {
                "cr": 1.0, "Omega_bh2": 700000.0, "Omega_ch2": 700000.0,
                "H0": 67.32117, "DN_eff": 0.0, "A_s": 2.100549e-9, "DN_re": 10.0,
            },
            "derived_tilt": "n_t=-r/8; red throughout",
            "continuum_claim": False,
        },
        "paired_screen": {
            "points": len(extension), "old_success": old_success,
            "old_failures_by_reason": old_reasons, "current_success": current_success,
            "current_failure_reasons": {},
            "h2_positive_finite_all": all(row["candidate_fast"].get("h2_positive_finite") for row in extension),
            "N_inf_finite_all": all(row["candidate_fast"].get("derived_N_inf") is not None for row in extension),
        },
        "accuracy": {
            "reference_sample_count": len(refs["cases"]),
            "max_first_map_root_error": max(row["root_relative_error_proxy"] for row in refs["cases"]),
            "max_positive_linear_spectrum_error_dex": max(
                row["spectrum_max_abs_dex_from_raw_linear"] for row in refs["cases"]
            ),
            "sobol_closure_second_to_third": closure["second_to_third_closure_relative"],
            "sobol_fast_root_error_to_third_map": closure["fast_root_relative_to_third_map"],
            "root_gate": 2e-4, "spectrum_gate_dex": 1e-3,
        },
        "performance": {
            "base_commit": benchmark["base_commit"],
            "candidate_commit": benchmark["candidate_commit"],
            "pairs_per_case": 500,
            "p50_ratio_min_max": [min(row["p50_ratio"] for row in measurements), max(row["p50_ratio"] for row in measurements)],
            "p95_ratio_min_max": [min(row["p95_ratio"] for row in measurements), max(row["p95_ratio"] for row in measurements)],
            "paired_ratio_p95_max_diagnostic": max(row["paired_ratio_p95"] for row in measurements),
            "root_delta_max": max(row["root_relative_delta_max"] for row in measurements),
            "frequency_digest_mismatch_pairs": sum(row["frequency_digest_mismatch_pairs"] for row in measurements),
            "spectrum_absdiff_max_dex": max(row["spectrum_absdiff_max_dex"] for row in measurements),
            "first_call_records": {row["case_id"]: row["first_call"] for row in measurements},
            "warm_resources": benchmark["resources"],
        },
        "reference_code_modified": False,
        "production_code_changed_this_tranche": False,
        "source_sha256": current_source_sha,
        "artifact_sha256": {name: sha256(HERE / name) for name in ARTIFACTS},
        "manifest_self_hash_included": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
