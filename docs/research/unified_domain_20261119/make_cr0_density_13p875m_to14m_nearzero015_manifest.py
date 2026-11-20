"""Create source and experiment hashes for the 14M CR0 edge slice."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr0_density_13p875m_to14m_nearzero015_manifest.json"
SOURCE_FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "screen_h2_identity_region.py",
    "screen_cr0_density_14m_nearzero015_slice.py",
    "screen_cr0_density_14m_nearzero015_slice.json",
    "screen_cr0_density_14m_nearzero015_slice_sobol.py",
    "screen_cr0_density_14m_nearzero015_slice_sobol.json",
    "screen_cr0_density_14m_upper_nearzero015_slice.py",
    "screen_cr0_density_14m_upper_nearzero015_slice.json",
    "screen_cr0_density_14m_upper_nearzero015_slice_sobol.py",
    "screen_cr0_density_14m_upper_nearzero015_slice_sobol.json",
    "probe_cr0_density_14m_nearzero015_reference_boundary.py",
    "probe_cr0_density_14m_nearzero015_reference_boundary.json",
    "validate_cr0_density_14m_nearzero015_slice_z8.py",
    "validate_cr0_density_14m_nearzero015_slice_z8.json",
    "probe_cr0_density_14m_upper_holdout05_reference_map.py",
    "probe_cr0_density_14m_upper_holdout05_reference_map.json",
    "diagnose_cr0_density_14m_upper_holdout05_reference_closure.py",
    "diagnose_cr0_density_14m_upper_holdout05_reference_closure.json",
    "benchmark_cr0_density_13m_to13_984375m_self_ab_500pairs.py",
    "benchmark_cr0_density_14m_upper_nearzero015_self_ab_500pairs.py",
    "benchmark_cr0_density_14m_upper_nearzero015_self_ab_tbb_500pairs.json",
    "plot_cr0_density_14m_upper_nearzero015_slice.py",
    "cr0_density_14m_upper_nearzero015_slice.png",
    "cr0_density_13p875m_to14m_nearzero015_summary.md",
    "make_cr0_density_13p875m_to14m_nearzero015_manifest.py",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = read_json("screen_cr0_density_14m_upper_nearzero015_slice.json")["cases"]
    sobol = read_json("screen_cr0_density_14m_upper_nearzero015_slice_sobol.json")["cases"]
    wide_lattice = read_json("screen_cr0_density_14m_nearzero015_slice.json")["cases"]
    wide_sobol = read_json("screen_cr0_density_14m_nearzero015_slice_sobol.json")["cases"]
    reference = read_json("validate_cr0_density_14m_nearzero015_slice_z8.json")
    holdout = read_json("probe_cr0_density_14m_upper_holdout05_reference_map.json")["case"]
    closure = read_json("diagnose_cr0_density_14m_upper_holdout05_reference_closure.json")
    timing = read_json("benchmark_cr0_density_14m_upper_nearzero015_self_ab_tbb_500pairs.json")
    accepted = lattice + sobol
    wide = wide_lattice + wide_sobol
    trusted = next(x for x in timing["cases"] if x["case_id"] == "trusted_physical_cr0_control")
    current_source = {name.rsplit("/", 1)[-1]: sha256((ROOT / name).read_bytes()) for name in SOURCE_FILES}
    baseline_source = {}
    for name in SOURCE_FILES:
        blob = subprocess.run(
            ["git", "show", f"59a23aec2300f08b50124d7f35dec4b228388d37:{name}"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
        ).stdout
        baseline_source[name.rsplit("/", 1)[-1]] = sha256(blob)
    ref_cases = reference["cases"]
    all_ref_roots = [r["root_relative_error"] for r in ref_cases] + [holdout["root_relative_error_at_fast_input"]]
    all_ref_spectra = [r["spectrum_max_abs_dex"] for r in ref_cases] + [holdout["spectrum_max_abs_dex"]]
    new_region_timing = [x for x in timing["cases"] if x["case_id"].startswith("new_region_")]
    payload = {
        "status": "validated_finite_sample_cr0_density_14m_nearzero015_edge",
        "source_base_commit": "d4ee1b9c8f9fad401c5e22c2316094a608074727",
        "old_failure_baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "performance_baseline_commit": "59a23aec2300f08b50124d7f35dec4b228388d37",
        "production_source_changed": False,
        "reference_code_modified": False,
        "source_sha256": {
            "current_worktree": current_source,
            "performance_baseline_59a23ae": baseline_source,
            "reference_used": current_source["reference.py"],
        },
        "accepted_region": {
            "Omega_bh2": [13875000.0, 14000000.0], "Omega_ch2": 0.1201075,
            "cr": 0.0, "n_t": [-0.015625, 0.0], "r": 0.01,
            "T_re": 2000.0, "kappa10": 0.01,
            "tilt_contract": "CR0 input n_t is independent",
        },
        "paired_screen": {
            "point_count": len(accepted),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in accepted),
            "candidate_successes": sum(bool(r["candidate_fast"]["success"]) for r in accepted),
            "old_failure_reasons": sorted({r["baseline_fast"]["failure_reason"] for r in accepted if not r["baseline_fast"]["success"]}),
            "candidate_positive_finite_h2_all": all(r["candidate_fast"]["h2_positive_finite"] for r in accepted),
            "candidate_min_sampled_h2": min(r["candidate_fast"]["h2_min"] for r in accepted),
            "candidate_finite_cutoff_all": all(r["candidate_fast"]["derived_N_inf"] is not None for r in accepted),
            "candidate_frequency_count_all_75": all(r["candidate_fast"]["n_freq"] == 75 for r in accepted),
            "candidate_warning_count_max": max(r["candidate_fast"]["warning_count"] for r in accepted),
            "sobol_seed": read_json("screen_cr0_density_14m_upper_nearzero015_slice_sobol.json")["sampling_design"]["seed"],
            "sobol_point_count": len(sobol),
        },
        "wider_exploratory_box": {
            "Omega_bh2": [13750000.0, 14000000.0], "n_t": [-0.015625, 0.0],
            "point_count": len(wide),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in wide),
            "candidate_successes": sum(bool(r["candidate_fast"]["success"]) for r in wide),
            "excluded_reference_boundary": read_json("probe_cr0_density_14m_nearzero015_reference_boundary.json")["cases"][0],
        },
        "reference_checks": {
            "z_tail": 8.0, "rtol": 1.0e-9,
            "selected_cases": [
                {k: r[k] for k in ("case_id", "params", "fast_root", "closed_reference_root", "root_relative_error", "map_closure_residual", "spectrum_max_abs_dex", "spectrum_valid_count")}
                for r in ref_cases
            ],
            "independent_sobol_single_map": {
                "params": holdout["params"],
                "root_relative_error": holdout["root_relative_error_at_fast_input"],
                "fixed_point_residual": holdout["fixed_point_residual_at_fast_input"],
                "spectrum_max_abs_dex": holdout["spectrum_max_abs_dex"],
                "additional_closure_map": closure["closure_step"],
            },
            "max_root_or_residual": max(all_ref_roots),
            "max_spectrum_error_dex": max(all_ref_spectra),
        },
        "performance": {
            "resources": timing["resources"],
            "500_pair_full_solver": [
                {k: r[k] for k in ("case_id", "pairs", "base_p50_seconds", "current_p50_seconds", "p50_ratio", "base_p95_seconds", "current_p95_seconds", "p95_ratio", "paired_ratio_p95", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex", "first_call")}
                for r in timing["cases"]
            ],
            "max_new_region_p50_ratio": max(r["p50_ratio"] for r in new_region_timing),
            "max_new_region_p95_ratio": max(r["p95_ratio"] for r in new_region_timing),
            "trusted_control": {k: trusted[k] for k in ("pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex")},
        },
        "artifact_sha256": {name: sha256((HERE / name).read_bytes()) for name in ARTIFACT_FILES},
        "manifest_self_hash_included": False,
        "continuum_claim": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
