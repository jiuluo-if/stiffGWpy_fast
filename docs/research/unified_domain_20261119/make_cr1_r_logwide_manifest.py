"""Create source and evidence hashes for the CR1 log-wide r extension."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr1_r_logwide_manifest.json"
SOURCE_FILES = [
    "exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "_wkb_transfer.py",
    "LCDM_stiff_Neff.py", "reference.py",
]
ARTIFACTS = [
    "screen_cr1_r_logwide.py", "screen_cr1_r_logwide_lattice.json", "screen_cr1_r_logwide_sobol.json",
    "validate_cr1_r_logwide_z8.py", "validate_cr1_r_logwide_z8.json",
    "validate_cr1_r_cutoff_probe_z8.py", "validate_cr1_r_cutoff_probe_z8.json",
    "validate_cr1_r_log_tail_probe_z8.py", "validate_cr1_r_log_tail_probe_z8.json",
    "audit_cr1_r_old_success_reference.py", "audit_cr1_r_old_success_reference.json",
    "close_reference_cr1_r_lowstrip_z8.py", "close_reference_cr1_r_lowstrip_z8.json",
    "benchmark_cr1_r_logwide_ab.py", "benchmark_cr1_r_logwide_ab_500pairs.json",
    "plot_cr1_r_logwide.py", "cr1_r_logwide_evidence.png", "cr1_r_logwide_summary.md",
    "make_cr1_r_logwide_manifest.py",
    "screen_cr1_r_lowstrip.py", "screen_cr1_r_lowstrip_lattice.json", "screen_cr1_r_lowstrip_sobol.json",
    "validate_cr1_r_lowstrip_upper_probe_z8.py", "validate_cr1_r_lowstrip_upper_probe_z8.json",
    "validate_cr1_r_lowstrip_z8.py", "validate_cr1_r_lowstrip_z8.json",
    "screen_cr1_r_reheat.py", "screen_cr1_r_reheat_lattice.json", "screen_cr1_r_reheat_sobol.json",
    "validate_cr1_r_reheat_z8.py", "validate_cr1_r_reheat_z8.json",
    "validate_cr1_r_reheat_upper_probe_z8.py", "validate_cr1_r_reheat_upper_probe_z8.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    payload = {
        "status": "validated_finite_sample_cr1_logwide_rectangle",
        "region": {
            "cr": 1.0, "Omega_bh2": 700000.0, "Omega_ch2": 700000.0,
            "T_re": [1750.0, 2000.0], "r": [4e-43, 0.005], "kappa10": 0.01,
            "derived_nt": "-r/8", "sobol_measure": "uniform in T_re and log10(r)",
            "cutoff_check_r_min_max_101_T_samples": 3.602014777069887e-43,
        },
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "performance_baseline_commit": "d4ee1b9c8f9fad401c5e22c2316094a608074727",
        "reference_code_modified": False,
        "production_code_changed_this_tranche": False,
        "continuum_claim": False,
        "screen_counts": {
            "paired_points": 19, "old_fast_success": 1, "old_fast_failures": 18,
            "old_failure_reason_counts": {"nonfinite_background_primitive": 18},
            "current_fast_success": 19, "independent_log_sobol": 8,
        },
        "accuracy": {
            "reference_points": 7,
            "max_root_relative_to_third_map": 0.00012885708587717442,
            "max_spectrum_abs_error_dex_from_raw_linear": 0.0008799070981417856,
            "old_success_point_old_root_reference_relative_error": 0.05303642586737549,
            "old_success_point_current_root_reference_relative_error": 0.00001171860112758007,
        },
        "performance_500_pairs": {
            "new_region_p50_ratio_min": 0.9829206931676476,
            "new_region_p50_ratio_max": 1.0230819713613586,
            "new_region_p95_ratio_max": 1.0250009891276477,
            "trusted_cr0_p50_ratio": 0.9993449427168174,
            "trusted_cr0_p95_ratio": 0.980970154214752,
        },
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name) for name in SOURCE_FILES
        },
        "artifact_sha256": {name: sha256(HERE / name) for name in ARTIFACTS},
        "manifest_self_hash_included": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
