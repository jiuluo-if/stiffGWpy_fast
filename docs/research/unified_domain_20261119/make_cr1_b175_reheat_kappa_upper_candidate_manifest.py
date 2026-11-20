from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr1_b175_reheat_kappa_upper_candidate_manifest.json"
FAST_SOURCE_FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "screen_cr1_b175_reheat_kappa_upper_lattice.json",
    "screen_cr1_b175_reheat_kappa_upper_sobol.json",
    "screen_cr1_b175_reheat_kappa_upper_candidate_replay.py",
    "screen_h2_identity_region.py",
    "screen_cr1_b175_reheat_kappa_upper_candidate_lattice.json",
    "screen_cr1_b175_reheat_kappa_upper_candidate_sobol.json",
    "benchmark_cr1_upper_candidate_vs_fda_replay.py",
    "benchmark_cr1_upper_candidate_vs_fda_replay_tbb.json",
    "benchmark_cr1_upper_vs_fda_candidate_region_100pairs.py",
    "benchmark_cr1_upper_vs_fda_candidate_region_tbb_100pairs.json",
    "benchmark_cr1_upper_vs_fda_head_process_ab.py",
    "benchmark_cr1_upper_vs_fda_head_process_tbb_500pairs.json",
    "validate_cr1_b175_reheat_kappa_upper_z8.py",
    "validate_cr1_b175_reheat_kappa_upper_z8.json",
    "close_reference_cr1_b175_reheat_kappa_upper_z8.py",
    "close_reference_cr1_b175_reheat_kappa_upper_z8.json",
    "audit_cr1_upper_old_success_point_z8.json",
    "benchmark_cr1_b175_reheat_kappa_upper_wkb_100pairs.json",
    "validate_cr1_b175_reheat_kappa_narrow_z8.py",
    "validate_cr1_b175_reheat_kappa_narrow_z8.json",
    "validate_cr1_b175_reheat_kappa_z8.py",
    "validate_cr1_b175_reheat_kappa_z8.json",
    "plot_cr1_b175_reheat_kappa_upper.py",
    "cr1_b175_reheat_kappa_upper_evidence.png",
    "cr1_b175_reheat_kappa_upper_summary.md",
    "make_cr1_b175_reheat_kappa_upper_candidate_manifest.py",
    "tests/test_long_reheating_domain.py",
    "tests/test_pressure_fast_domain.py",
    "progress.md",
    "findings.md",
    "task_plan.md",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    lattice = json.loads((HERE / "screen_cr1_b175_reheat_kappa_upper_candidate_lattice.json").read_text(encoding="utf-8"))
    sobol = json.loads((HERE / "screen_cr1_b175_reheat_kappa_upper_candidate_sobol.json").read_text(encoding="utf-8"))
    replay = json.loads((HERE / "benchmark_cr1_upper_candidate_vs_fda_replay_tbb.json").read_text(encoding="utf-8"))
    region_ab = json.loads((HERE / "benchmark_cr1_upper_vs_fda_candidate_region_tbb_100pairs.json").read_text(encoding="utf-8"))
    trusted_ab = json.loads((HERE / "benchmark_cr1_upper_vs_fda_head_process_tbb_500pairs.json").read_text(encoding="utf-8"))
    reference = json.loads((HERE / "validate_cr1_b175_reheat_kappa_upper_z8.json").read_text(encoding="utf-8"))
    closure = json.loads((HERE / "close_reference_cr1_b175_reheat_kappa_upper_z8.json").read_text(encoding="utf-8"))

    baseline_source = {}
    for filename in FAST_SOURCE_FILES:
        data = subprocess.run(
            ["git", "show", f"fda343d661e74de911276efc647a16a40d3b5441:{filename}"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
        ).stdout
        baseline_source[Path(filename).name] = sha256(data)

    rows = lattice["cases"] + sobol["cases"]
    baseline_success = sum(row["baseline_fast"]["success"] for row in rows)
    candidate_success = sum(row["candidate_fast"]["success"] for row in rows)
    reference_rows = reference["cases"]
    closure_rows = closure["cases"]
    source_files = FAST_SOURCE_FILES + (
        "tests/test_long_reheating_domain.py",
        "tests/test_pressure_fast_domain.py",
    )
    payload = {
        "status": "validated_finite_sample_cr1_adjacent_region_and_equivalent_fast_optimization",
        "source_base_commit": "fda343d661e74de911276efc647a16a40d3b5441",
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "reference_code_modified": False,
        "source_hashes": {
            "candidate": {Path(name).name: sha256((ROOT / name).read_bytes()) for name in source_files},
            "base_fda343d": baseline_source,
            "reference_reused_source": reference["source_sha256"],
        },
        "screen": {
            "parameter_region": {"Omega_bh2": 700000.0, "Omega_ch2": 700000.0, "cr": 1.0, "T_re": [1750.0, 2000.0], "kappa10": [0.0125, 0.015]},
            "derived_tilt": -0.00125,
            "old_successes": baseline_success,
            "candidate_successes": candidate_success,
            "point_count": len(rows),
            "old_failure_reasons": sorted({row["baseline_fast"]["failure_reason"] for row in rows if not row["baseline_fast"]["success"]}),
            "candidate_full_spectrum_points": sum(row["candidate_fast"]["finite_spectrum"] and row["candidate_fast"]["n_freq"] == 76 for row in rows),
            "source_equivalence_to_fda": {
                "paired_points": len(replay["cases"]),
                "root_relative_delta_max": max(row["root_relative_delta_max"] for row in replay["cases"]),
                "frequency_digest_mismatch_pairs": sum(row["frequency_digest_mismatch_pairs"] for row in replay["cases"]),
                "spectrum_absdiff_max_dex": max(row["spectrum_absdiff_max_dex"] for row in replay["cases"]),
            },
        },
        "reference_reuse": {
            "artifact": "validate_cr1_b175_reheat_kappa_upper_z8.json",
            "reference_code_sha256": reference["source_sha256"]["reference.py"],
            "max_root_error_proxy": max(row["root_relative_error_proxy"] for row in reference_rows),
            "max_spectrum_error_dex": max(row["spectrum_max_abs_dex"] for row in reference_rows),
            "third_map_closure_artifact": "close_reference_cr1_b175_reheat_kappa_upper_z8.json",
            "max_third_map_root_error": max(row["fast_root_relative_to_third_map"] for row in closure_rows),
            "max_second_to_third_map_closure": max(row["second_to_third_closure_relative"] for row in closure_rows),
            "reuse_basis": "All 19 current Fast roots, frequency digests, and spectra were bitwise replayed against fda343d; reused Reference inputs therefore remain unchanged.",
        },
        "performance": {
            "resources": trusted_ab["resources"],
            "trusted_and_old_success_500_pairs": [
                {key: row[key] for key in ("case_id", "pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex", "first_call")}
                for row in trusted_ab["cases"]
            ],
            "new_region_100_pairs": [
                {key: row[key] for key in ("case_id", "pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex")}
                for row in region_ab["cases"]
            ],
            "independent_numba_caches": True,
            "alternating_order": True,
            "cold_first_solve_includes_jit": True,
        },
        "artifact_sha256": {
            name: sha256(((HERE / name) if (HERE / name).exists() else (ROOT / name)).read_bytes())
            for name in ARTIFACT_FILES
        },
        "manifest_self_hash_included": False,
        "continuum_claim": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
