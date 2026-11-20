from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr0_density_13m_to13_75m_manifest.json"
SOURCE_FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "screen_h2_identity_region.py",
    "screen_cr0_density_13m_to13_75m.py",
    "screen_cr0_density_13m_to13_75m.json",
    "screen_cr0_density_13m_to13_75m_sobol.py",
    "screen_cr0_density_13m_to13_75m_sobol.json",
    "screen_cr0_density_13m_to14m.py",
    "screen_cr0_density_13m_to14m.json",
    "screen_cr0_density_13m_to14m_sobol.py",
    "screen_cr0_density_13m_to14m_sobol.json",
    "validate_cr0_density_13m_to13_75m_z8.py",
    "validate_cr0_density_13m_to13_75m_z8.json",
    "validate_cr0_density_13m_to14m_interior_z8.py",
    "validate_cr0_density_13m_to14m_interior_z8.json",
    "validate_cr0_density_13m_to14m_upper_interior_z8.py",
    "validate_cr0_density_13m_to14m_upper_interior_z8.json",
    "validate_cr0_density_13m_to14m_z8.py",
    "validate_cr0_density_13m_to14m_z8_prufer.py",
    "diagnose_cr0_density_14m_reference_step_underflow.py",
    "diagnose_cr0_density_14m_reference_step_underflow.json",
    "benchmark_cr0_density_13m_to13_75m_self_ab_500pairs.py",
    "benchmark_cr0_density_13m_to13_75m_self_ab_tbb_500pairs.json",
    "plot_cr0_density_13m_to13_75m.py",
    "cr0_density_13m_to13_75m_domain.png",
    "cr0_density_13m_to13_75m_summary.md",
    "make_cr0_density_13m_to13_75m_manifest.py",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def artifact_path(name: str) -> Path:
    local = HERE / name
    return local if local.exists() else ROOT / name


def main() -> None:
    lattice = read_json("screen_cr0_density_13m_to13_75m.json")
    sobol = read_json("screen_cr0_density_13m_to13_75m_sobol.json")
    broad = read_json("screen_cr0_density_13m_to14m.json")
    broad_sobol = read_json("screen_cr0_density_13m_to14m_sobol.json")
    reference = read_json("validate_cr0_density_13m_to13_75m_z8.json")
    interior = read_json("validate_cr0_density_13m_to14m_interior_z8.json")
    timing = read_json("benchmark_cr0_density_13m_to13_75m_self_ab_tbb_500pairs.json")
    failure = read_json("diagnose_cr0_density_14m_reference_step_underflow.json")
    paired = lattice["cases"] + sobol["cases"]
    broad_rows = broad["cases"] + broad_sobol["cases"]
    trusted = next(r for r in timing["cases"] if r["case_id"] == "trusted_physical_cr0_control")
    current_source = {name.rsplit("/", 1)[-1]: sha256((ROOT / name).read_bytes()) for name in SOURCE_FILES}
    baseline_source = {}
    for name in SOURCE_FILES:
        blob = subprocess.run(
            ["git", "show", f"59a23aec2300f08b50124d7f35dec4b228388d37:{name}"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
        ).stdout
        baseline_source[name.rsplit("/", 1)[-1]] = sha256(blob)
    refs = reference["cases"] + interior["cases"]
    payload = {
        "status": "validated_finite_sample_cr0_density_extension_13m_to13_75m_with_14m_reference_limit",
        "source_base_commit": "59a23aec2300f08b50124d7f35dec4b228388d37",
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "reference_code_modified": False,
        "production_source_changed": False,
        "source_sha256": {
            "current_worktree": current_source,
            "source_base_59a23ae": baseline_source,
            "reference_code_used": reference["source_sha256"]["reference.py"],
        },
        "accepted_region": {
            "Omega_bh2": [13000000.0, 13750000.0], "Omega_ch2": 0.1201075,
            "cr": 0.0, "n_t": [-0.5, 0.0], "r": 0.01,
            "T_re": 2000.0, "kappa10": 0.01,
            "tilt_contract": "CR0 input n_t is independent",
        },
        "paired_screen": {
            "point_count": len(paired),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in paired),
            "candidate_successes": sum(bool(r["candidate_fast"]["success"]) for r in paired),
            "old_failure_reasons": sorted({r["baseline_fast"]["failure_reason"] for r in paired if not r["baseline_fast"]["success"]}),
            "candidate_positive_finite_h2_all": all(r["candidate_fast"]["h2_positive_finite"] for r in paired),
            "candidate_finite_cutoff_all": all(r["candidate_fast"]["derived_N_inf"] is not None for r in paired),
            "candidate_frequency_count_all": all(r["candidate_fast"]["n_freq"] == 75 for r in paired),
            "candidate_warning_count_max": max(r["candidate_fast"]["warning_count"] for r in paired),
            "sobol_seed": sobol["sampling_design"]["seed"],
            "sobol_point_count": len(sobol["cases"]),
        },
        "broader_exploratory_13m_to14m_screen": {
            "point_count": len(broad_rows),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in broad_rows),
            "current_fast_successes": sum(bool(r["candidate_fast"]["success"]) for r in broad_rows),
            "reference_validated_region_excludes_14m_face": True,
        },
        "reference_checks": {
            "z_tail": 8.0, "rtol": 1.0e-9,
            "accepted_region_cases": [
                {key: row[key] for key in ("case_id", "params", "fast_root", "closed_reference_root", "root_relative_error", "map_closure_residual", "spectrum_max_abs_dex", "spectrum_valid_count")}
                for row in refs
            ],
            "max_root_error": max(r["root_relative_error"] for r in refs),
            "max_spectrum_error_dex": max(r["spectrum_max_abs_dex"] for r in refs),
            "max_map_closure": max(r["map_closure_residual"] for r in refs),
        },
        "reference_boundary_failure": failure,
        "performance": {
            "resources": timing["resources"],
            "500_pair_full_solver": [
                {key: row[key] for key in ("case_id", "pairs", "base_p50_seconds", "current_p50_seconds", "p50_ratio", "base_p95_seconds", "current_p95_seconds", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex", "first_call")}
                for row in timing["cases"]
            ],
            "trusted_control_500_pair": {key: trusted[key] for key in ("pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex")},
        },
        "artifact_sha256": {name: sha256(artifact_path(name).read_bytes()) for name in ARTIFACT_FILES},
        "manifest_self_hash_included": False,
        "continuum_claim": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
