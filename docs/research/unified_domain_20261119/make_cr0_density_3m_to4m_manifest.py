from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr0_density_3m_to4m_manifest.json"
SOURCE_FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "screen_h2_identity_region.py",
    "screen_cr0_density_3m_to4m.py",
    "screen_cr0_density_3m_to4m.json",
    "screen_cr0_density_3m_to4m_sobol.py",
    "screen_cr0_density_3m_to4m_sobol.json",
    "replay_cr0_3m_reference_boundary_current.py",
    "replay_cr0_3m_reference_boundary_current.json",
    "close_cr0_3m_reference_origin_corrected.py",
    "close_cr0_3m_reference_origin_corrected.json",
    "validate_cr0_density_3m_to4m_z8.py",
    "validate_cr0_density_3m_to4m_z8.json",
    "benchmark_cr0_density_3m_to4m_self_ab_100pairs.py",
    "benchmark_cr0_density_3m_to4m_self_ab_tbb_100pairs.json",
    "validate_cr0_3m_wkb_guard_z8.py",
    "validate_cr0_3m_wkb_guard_z8_multi.json",
    "cr0_2m_to3m_wkb_basis_guard_summary.md",
    "benchmark_cr1_upper_vs_fda_head_process_ab.py",
    "benchmark_cr1_upper_vs_fda_head_process_tbb_500pairs.json",
    "plot_cr0_density_3m_to4m.py",
    "cr0_density_3m_to4m_domain.png",
    "cr0_density_3m_to4m_summary.md",
    "make_cr0_density_3m_to4m_manifest.py",
    "tests/test_long_reheating_domain.py",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def artifact_path(name: str) -> Path:
    local = HERE / name
    return local if local.exists() else ROOT / name


def main() -> None:
    lattice = read_json("screen_cr0_density_3m_to4m.json")
    sobol = read_json("screen_cr0_density_3m_to4m_sobol.json")
    reference = read_json("validate_cr0_density_3m_to4m_z8.json")
    closure = read_json("close_cr0_3m_reference_origin_corrected.json")
    replay = read_json("replay_cr0_3m_reference_boundary_current.json")
    perf = read_json("benchmark_cr0_density_3m_to4m_self_ab_tbb_100pairs.json")
    trusted = read_json("benchmark_cr1_upper_vs_fda_head_process_tbb_500pairs.json")

    paired = lattice["cases"] + sobol["cases"]
    refs = reference["cases"]
    old_boundary = closure["cases"]
    trusted_row = next(r for r in trusted["cases"] if r["case_id"] == "trusted_physical_cr0_control")
    current_source = {
        name.rsplit("/", 1)[-1]: sha256((ROOT / name).read_bytes()) for name in SOURCE_FILES
    }
    base_source = {}
    for name in SOURCE_FILES:
        blob = subprocess.run(
            ["git", "show", f"59a23aec2300f08b50124d7f35dec4b228388d37:{name}"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
        ).stdout
        base_source[name.rsplit("/", 1)[-1]] = sha256(blob)

    payload = {
        "status": "validated_finite_sample_cr0_density_extension_3m_to4m",
        "source_base_commit": "59a23aec2300f08b50124d7f35dec4b228388d37",
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "reference_code_modified": False,
        "production_source_changed": False,
        "source_sha256": {
            "current_worktree": current_source,
            "source_base_59a23ae": base_source,
            "reference_code_used": reference["source_sha256"]["reference.py"],
        },
        "region": {
            "Omega_bh2": [3000000.0, 4000000.0], "Omega_ch2": 0.1201075,
            "cr": 0.0, "n_t": [-0.5, 0.0], "r": 0.01,
            "T_re": 2000.0, "kappa10": 0.01,
            "tilt_contract": "CR0 input n_t is independent",
        },
        "paired_screen": {
            "point_count": len(paired),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in paired),
            "candidate_successes": sum(bool(r["candidate_fast"]["success"]) for r in paired),
            "old_failure_reasons": sorted({r["baseline_fast"]["failure_reason"] for r in paired if not r["baseline_fast"]["success"]}),
            "candidate_full_frequency_count": 75,
            "candidate_positive_finite_h2_all": all(r["candidate_fast"]["h2_positive_finite"] for r in paired),
            "candidate_finite_cutoff_all": all(r["candidate_fast"]["derived_N_inf"] is not None for r in paired),
            "candidate_warning_count_max": max(r["candidate_fast"]["warning_count"] for r in paired),
            "sobol_seed": sobol["sampling_design"]["seed"],
            "sobol_point_count": len(sobol["cases"]),
        },
        "reference_checks": {
            "z_tail": 8.0, "rtol": 1e-9,
            "selected_new_region_cases": [
                {key: row[key] for key in ("case_id", "params", "fast_root", "closed_reference_root", "root_relative_error", "map_closure_residual", "spectrum_max_abs_dex", "spectrum_valid_count")}
                for row in refs
            ],
            "max_root_error": max(r["root_relative_error"] for r in refs),
            "max_spectrum_error_dex": max(r["spectrum_max_abs_dex"] for r in refs),
            "max_map_closure": max(r["map_closure_residual"] for r in refs),
            "3m_origin_audit": {
                "old_closure_field_withdrawn": True,
                "current_fast_replays_bitwise": all(
                    all(r[key] for key in ("root_bitwise_equal", "frequency_bitwise_equal", "spectrum_bitwise_equal"))
                    for r in replay["cases"]
                ),
                "corrected_3m_cases": [
                    {
                        "case_id": row["case_id"],
                        "input_dn_eff_origin": row["input_dn_eff_origin"],
                        "fast_root": row["fast_root"],
                        "closed_reference_root": row["map_history"][-1]["mapped_dn_gw"],
                        "fast_root_relative_to_corrected_reference": row["fast_root_relative_to_corrected_reference"],
                        "map_closure_relative": row["map_closure_relative"],
                    }
                    for row in old_boundary
                ],
                "max_corrected_3m_root_error": max(r["fast_root_relative_to_corrected_reference"] for r in old_boundary),
                "max_corrected_3m_map_closure": max(r["map_closure_relative"] for r in old_boundary),
            },
        },
        "performance": {
            "resources": perf["resources"],
            "100_pair_full_solver": [
                {key: row[key] for key in ("case_id", "pairs", "base_p50_seconds", "current_p50_seconds", "p50_ratio", "base_p95_seconds", "current_p95_seconds", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex", "first_call")}
                for row in perf["cases"]
            ],
            "500_pair_trusted_control": {key: trusted_row[key] for key in ("pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex")},
        },
        "artifact_sha256": {name: sha256(artifact_path(name).read_bytes()) for name in ARTIFACT_FILES},
        "manifest_self_hash_included": False,
        "continuum_claim": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
