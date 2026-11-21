"""Create source and evidence hashes for the nearby CR1 extension."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr1_nearby_reheat_kappa_manifest.json"
SOURCE_FILES = [
    "exact_background.py", "fast_sgwb.py", "freq_adaptive.py", "_wkb_transfer.py",
    "LCDM_stiff_Neff.py", "reference.py",
]
ARTIFACTS = [
    "screen_cr1_nearby_reheat_kappa.py",
    "screen_cr1_nearby_reheat_kappa_lattice.json",
    "screen_cr1_nearby_reheat_kappa_sobol.json",
    "validate_cr1_nearby_reheat_kappa_z8.py",
    "validate_cr1_nearby_reheat_kappa_z8.json",
    "close_reference_cr1_nearby_reheat_kappa_z8.py",
    "close_reference_cr1_nearby_reheat_kappa_z8.json",
    "benchmark_cr1_nearby_reheat_kappa_ab.py",
    "benchmark_cr1_nearby_reheat_kappa_ab_500pairs.json",
    "plot_cr1_nearby_reheat_kappa.py",
    "make_cr1_nearby_reheat_kappa_manifest.py",
    "cr1_nearby_reheat_kappa_evidence.png",
    "cr1_nearby_reheat_kappa_summary.md",
    "screen_cr1_low_reheat_kappa.py",
    "screen_cr1_low_reheat_kappa_lattice.json",
    "screen_cr1_low_reheat_kappa_sobol.json",
    "validate_cr1_low_reheat_kappa_z8.py",
    "validate_cr1_low_reheat_kappa_z8.json",
    "close_reference_cr1_low_reheat_kappa_z8.py",
    "close_reference_cr1_low_reheat_kappa_z8.json",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    payload = {
        "status": "validated_finite_sample_adjacent_cr1_rectangle",
        "accepted_region": {
            "cr": 1.0, "r": 0.01, "derived_nt": -0.00125,
            "T_re": [1687.5, 1750.0], "kappa10": [0.0125, 0.015],
            "fixed_densities": {"Omega_bh2": 700000.0, "Omega_ch2": 700000.0},
        },
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "performance_baseline_commit": "d4ee1b9c8f9fad401c5e22c2316094a608074727",
        "reference_code_modified": False,
        "production_code_changed_this_tranche": False,
        "continuum_claim": False,
        "sample_counts": {"paired_old_success": 0, "paired_current_success": 19,
                          "reference_points": 3, "independent_sobol_holdouts": 8},
        "max_root_relative_to_third_reference_map": 0.0001942824847016016,
        "max_spectrum_abs_error_dex": 0.0009223931316419964,
        "performance_500_pairs": {
            "new_region_p50_ratio_max": 0.967568735844922,
            "new_region_p95_ratio_max": 0.9771221044934549,
            "trusted_cr0_p50_ratio": 0.9963020173479681,
            "trusted_cr0_p95_ratio": 0.9901314618289662,
        },
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name) for name in SOURCE_FILES
        },
        "artifact_sha256": {
            name: sha256(HERE / name) for name in ARTIFACTS
        },
        "manifest_self_hash_included": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
