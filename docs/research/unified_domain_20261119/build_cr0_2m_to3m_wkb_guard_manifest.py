"""Build source/artifact hashes for the CR0 2M–3M WKB guard evidence."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "cr0_2m_to3m_wkb_basis_guard_manifest.json"

SOURCE_FILES = (
    "stiffgwpy_fast/_wkb_transfer.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/reference.py",
)
ARTIFACT_FILES = (
    "tests/test_wkb3_fast_transfer.py",
    "tests/test_wkb_singular_mode_completion.py",
    "docs/research/unified_domain_20261119/screen_cr0_density_2m_to3m_wkb_guard.py",
    "docs/research/unified_domain_20261119/screen_cr0_density_2m_to3m_wkb_guard.json",
    "docs/research/unified_domain_20261119/screen_cr0_density_2m_to3m_wkb_guard_sobol.py",
    "docs/research/unified_domain_20261119/screen_cr0_density_2m_to3m_wkb_guard_sobol.json",
    "docs/research/unified_domain_20261119/validate_cr0_3m_wkb_guard_z8.py",
    "docs/research/unified_domain_20261119/validate_cr0_3m_wkb_guard_z8_multi.json",
    "docs/research/unified_domain_20261119/benchmark_wkb_basis_guard_ab.py",
    "docs/research/unified_domain_20261119/benchmark_wkb_basis_guard_ab.json",
    "docs/research/unified_domain_20261119/plot_cr0_2m_to3m_wkb_guard.py",
    "docs/research/unified_domain_20261119/cr0_density_2m_to3m_wkb_guard_domain.png",
    "docs/research/unified_domain_20261119/cr0_2m_to3m_wkb_basis_guard_summary.md",
    "docs/research/unified_domain_20261119/build_cr0_2m_to3m_wkb_guard_manifest.py",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    grid = json.loads(
        (HERE / "screen_cr0_density_2m_to3m_wkb_guard.json").read_text(encoding="utf-8")
    )
    sobol = json.loads(
        (HERE / "screen_cr0_density_2m_to3m_wkb_guard_sobol.json").read_text(
            encoding="utf-8"
        )
    )
    references = json.loads(
        (HERE / "validate_cr0_3m_wkb_guard_z8_multi.json").read_text(encoding="utf-8")
    )
    benchmark = json.loads(
        (HERE / "benchmark_wkb_basis_guard_ab.json").read_text(encoding="utf-8")
    )
    rows = grid["cases"] + sobol["cases"]
    reference_rows = references["cases"]
    source_hashes = {name: _sha256(ROOT / name) for name in SOURCE_FILES}
    artifact_hashes = {name: _sha256(ROOT / name) for name in ARTIFACT_FILES}

    payload = {
        "status": "complete_finite_sample_domain_evidence",
        "candidate_parent_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "historical_failure_baseline": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "helper_ab_baseline": benchmark["base_commit"],
        "region": {
            "cr": 0.0,
            "independent_n_t": True,
            "Omega_bh2": [2.0e6, 3.0e6],
            "n_t": [-0.5, 0.0],
            "fixed_parameters": {
                key: value
                for key, value in grid["cases"][0]["input"]["params"].items()
                if key not in ("Omega_bh2", "n_t")
            },
            "grid_and_interior_old_new_successes": [
                sum(bool(row["baseline_fast"]["success"]) for row in grid["cases"]),
                sum(bool(row["candidate_fast"]["success"]) for row in grid["cases"]),
            ],
            "sobol": {
                "seed": sobol.get("sampling_design", {}).get("seed"),
                "n": len(sobol["cases"]),
                "old_new_successes": [
                    sum(bool(row["baseline_fast"]["success"]) for row in sobol["cases"]),
                    sum(bool(row["candidate_fast"]["success"]) for row in sobol["cases"]),
                ],
            },
            "old_failure_classes": sorted(
                {
                    row["baseline_fast"].get("failure_reason")
                    for row in rows
                    if not row["baseline_fast"].get("success")
                }
            ),
            "current_failures": sum(
                not bool(row["candidate_fast"].get("success")) for row in rows
            ),
            "all_h2_positive_finite": all(
                bool(row["candidate_fast"].get("h2_positive_finite")) for row in rows
            ),
            "minimum_sampled_h2": min(
                float(row["candidate_fast"]["h2_min"]) for row in rows
            ),
            "all_N_inf_finite": all(
                row["candidate_fast"].get("derived_N_inf") is not None for row in rows
            ),
            "all_output_frequency_counts": sorted(
                {
                    row["candidate_fast"].get("n_freq")
                    for row in rows
                    if row["candidate_fast"].get("success")
                }
            ),
            "candidate_warning_count": sum(
                int(row["candidate_fast"].get("warning_count", 0)) for row in rows
            ),
        },
        "reference": {
            "z_tail": 8.0,
            "rtol": 1.0e-9,
            "reference_code_modified": False,
            "sample_count": len(reference_rows),
            "max_root_relative_error": max(
                float(row["root_relative_error"]) for row in reference_rows
            ),
            "max_spectrum_error_dex": max(
                float(row["spectrum_max_abs_dex"]) for row in reference_rows
            ),
            "max_reference_closure_residual": max(
                float(row["reference_closure_residual"]) for row in reference_rows
            ),
            "samples": [
                {
                    "case_id": row["case_id"],
                    "params": row["params"],
                    "root_relative_error": row["root_relative_error"],
                    "spectrum_max_abs_dex": row["spectrum_max_abs_dex"],
                    "reference_closure_residual": row["reference_closure_residual"],
                }
                for row in reference_rows
            ],
        },
        "performance": {
            "baseline_helper_sha256": benchmark["baseline_helper_sha256"],
            "candidate_helper_sha256": benchmark["candidate_helper_sha256"],
            "resources": benchmark["resources"],
            "cases": [
                {
                    "case_id": row["case_id"],
                    "first_call_compilation_cold": row["first_call_compilation_cold"],
                    "first_baseline_seconds": row["first_base"]["seconds"],
                    "first_candidate_seconds": row["first_candidate"]["seconds"],
                    "p50_ratio": row["p50_ratio"],
                    "p95_ratio": row["p95_ratio"],
                    "root_relative_difference_max": row[
                        "root_relative_difference_max"
                    ],
                    "frequency_digest_mismatch_pairs": row[
                        "frequency_digest_mismatch_pairs"
                    ],
                    "spectrum_absdiff_max_dex": row["spectrum_absdiff_max_dex"],
                }
                for row in benchmark["cases"]
            ],
        },
        "numerical_method": {
            "trigger": "WKB3 basis determinant is nonfinite, subnormal, or within its float64 product-subtraction error bound",
            "fallback": "existing _phase_segment_base phase-substep transfer",
            "error_bound": "4u*(abs(a11*a22)+abs(a12*a21)), u=2.220446049250313e-16",
            "literature": {
                "source": "Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., SIAM",
                "chapters_read": ["2.1-2.2", "3.1"],
                "relation": "standard floating-point rounding model and dot-product product-sum error bound",
                "assumptions": "finite binary64 operations; subnormal determinant rejected explicitly",
            },
        },
        "local_validation": {
            "pytest": "438 passed, 9 deselected, 2 warnings",
            "focused_wkb_tests": "6 passed",
            "ruff": "passed",
            "mypy": "passed: 2 configured source files",
            "compileall": "passed",
            "wheel_sdist_build": "passed",
            "installed_wheel_guard_smoke": "passed",
            "github_actions": "pending push",
        },
        "source_sha256": source_hashes,
        "artifact_sha256": artifact_hashes,
        "boundary": "finite sampled support only; not a continuous-domain proof; cr>0 is not covered",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"manifest: {OUTPUT}")
    print(f"source hashes: {json.dumps(source_hashes, sort_keys=True)}")


if __name__ == "__main__":
    main()
