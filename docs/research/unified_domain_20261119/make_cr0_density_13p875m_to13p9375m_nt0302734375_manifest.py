"""Create source and experiment hashes for the density-dependent tilt slice."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "cr0_density_13p875m_to13p9375m_nt0302734375_manifest.json"
SOURCE_FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/freq_adaptive.py",
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "screen_h2_identity_region.py",
    "screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice.py",
    "screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice.json",
    "screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice_sobol.py",
    "screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice_sobol.json",
    "probe_cr0_density_14m_nt0302734375_reference.py",
    "probe_cr0_density_14m_nt0302734375_reference.json",
    "probe_cr0_density_13p875m_to13p9375m_nt0302734375_sobol_reference.py",
    "probe_cr0_density_13p875m_to13p9375m_nt0302734375_sobol_reference.json",
    "probe_cr0_density_13p9375m_nt0302734375_closure_step.py",
    "probe_cr0_density_13p9375m_nt0302734375_closure_step.json",
    "probe_cr0_density_14m_nearzero015_reference_boundary.py",
    "probe_cr0_density_14m_nearzero015_reference_boundary.json",
    "benchmark_cr0_density_13m_to13_984375m_self_ab_500pairs.py",
    "benchmark_cr0_density_13p875m_to13p9375m_nt0302734375_self_ab_500pairs.py",
    "benchmark_cr0_density_13p875m_to13p9375m_nt0302734375_self_ab_tbb_500pairs.json",
    "plot_cr0_density_13p875m_to13p9375m_nt0302734375_slice.py",
    "cr0_density_13p875m_to13p9375m_nt0302734375_slice.png",
    "cr0_density_13p875m_to13p9375m_nt0302734375_summary.md",
    "make_cr0_density_13p875m_to13p9375m_nt0302734375_manifest.py",
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = read_json("screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice.json")["cases"]
    sobol = read_json("screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice_sobol.json")["cases"]
    screened = lattice + sobol
    boundary = read_json("probe_cr0_density_14m_nt0302734375_reference.json")["cases"]
    holdout = read_json("probe_cr0_density_13p875m_to13p9375m_nt0302734375_sobol_reference.json")["case"]
    closure = read_json("probe_cr0_density_13p9375m_nt0302734375_closure_step.json")
    timing = read_json("benchmark_cr0_density_13p875m_to13p9375m_nt0302734375_self_ab_tbb_500pairs.json")
    trusted = next(r for r in timing["cases"] if r["case_id"] == "trusted_physical_cr0_control")
    accepted_refs = [r for r in boundary if r["reference_status"] == "success"]
    current_source = {name.rsplit("/", 1)[-1]: sha256((ROOT / name).read_bytes()) for name in SOURCE_FILES}
    baseline_source = {}
    for name in SOURCE_FILES:
        blob = subprocess.run(
            ["git", "show", f"59a23aec2300f08b50124d7f35dec4b228388d37:{name}"],
            cwd=ROOT, check=True, stdout=subprocess.PIPE,
        ).stdout
        baseline_source[name.rsplit("/", 1)[-1]] = sha256(blob)
    new_timing = [r for r in timing["cases"] if r["case_id"].startswith("new_region_")]
    root_errors = [r["root_relative_error"] for r in accepted_refs] + [holdout["root_relative_error"]]
    spectrum_errors = [r["spectrum_max_abs_dex"] for r in accepted_refs] + [holdout["spectrum_max_abs_dex"]]
    payload = {
        "status": "validated_finite_sample_cr0_density_dependent_tilt_extension_to_minus_0p0302734375",
        "source_base_commit": "ba66e074ab8630829441c520794246d3d8b8db58",
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
            "Omega_bh2": [13875000.0, 13937500.0], "Omega_ch2": 0.1201075,
            "cr": 0.0, "n_t": [-0.0302734375, 0.0], "r": 0.01,
            "T_re": 2000.0, "kappa10": 0.01,
            "tilt_contract": "CR0 input n_t is independent",
        },
        "paired_screen": {
            "point_count": len(screened),
            "old_successes": sum(bool(r["baseline_fast"]["success"]) for r in screened),
            "candidate_successes": sum(bool(r["candidate_fast"]["success"]) for r in screened),
            "old_failure_reasons": sorted({r["baseline_fast"]["failure_reason"] for r in screened if not r["baseline_fast"]["success"]}),
            "candidate_positive_finite_h2_all": all(r["candidate_fast"]["h2_positive_finite"] for r in screened),
            "candidate_min_sampled_h2": min(r["candidate_fast"]["h2_min"] for r in screened),
            "candidate_finite_cutoff_all": all(r["candidate_fast"]["derived_N_inf"] is not None for r in screened),
            "candidate_frequency_count_all_75": all(r["candidate_fast"]["n_freq"] == 75 for r in screened),
            "candidate_warning_count_max": max(r["candidate_fast"]["warning_count"] for r in screened),
            "sobol_seed": read_json("screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice_sobol.json")["sampling_design"]["seed"],
            "sobol_point_count": len(sobol),
        },
        "reference_checks": {
            "z_tail": 8.0, "rtol": 1.0e-9,
            "accepted_density_cases": accepted_refs,
            "excluded_14m_boundary_case": next(r for r in boundary if r["case_id"] == "upper_density_edge"),
            "independent_sobol_holdout": holdout,
            "13p9375m_second_closure_step": closure["closure_step"],
            "max_accepted_root_residual": max(root_errors),
            "max_accepted_spectrum_error_dex": max(spectrum_errors),
        },
        "performance": {
            "resources": timing["resources"],
            "500_pair_full_solver": [
                {k: r[k] for k in (
                    "case_id", "pairs", "base_p50_seconds", "current_p50_seconds", "p50_ratio",
                    "base_p95_seconds", "current_p95_seconds", "p95_ratio", "paired_ratio_p95",
                    "root_relative_delta_max", "frequency_digest_mismatch_pairs",
                    "spectrum_absdiff_max_dex", "first_call",
                )}
                for r in timing["cases"]
            ],
            "max_new_region_p50_ratio": max(r["p50_ratio"] for r in new_timing),
            "max_new_region_p95_ratio": max(r["p95_ratio"] for r in new_timing),
            "trusted_control": {k: trusted[k] for k in ("pairs", "p50_ratio", "p95_ratio", "root_relative_delta_max", "frequency_digest_mismatch_pairs", "spectrum_absdiff_max_dex")},
        },
        "continuum_claim": False,
        "artifact_sha256": {name: sha256((HERE / name).read_bytes()) for name in ARTIFACT_FILES},
        "manifest_self_hash_included": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
