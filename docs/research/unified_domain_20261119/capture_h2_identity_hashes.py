"""Capture source/runner hashes for the accepted Fast candidate artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUT = HERE / "source_hashes_h2_identity_current.json"
FILES = (
    "stiffgwpy_fast/exact_background.py",
    "stiffgwpy_fast/fast_sgwb.py",
    "stiffgwpy_fast/reference.py",
    "tests/test_long_reheating_domain.py",
    "tests/test_physics_limits.py",
    "tests/test_pressure_fast_domain.py",
    "docs/research/unified_domain_20261119/screen_h2_identity_region.py",
    "docs/research/unified_domain_20261119/screen_h2_identity_density_slice.py",
    "docs/research/unified_domain_20261119/screen_h2_identity_cr1_density.py",
    "docs/research/unified_domain_20261119/validate_h2_identity_region_z8.py",
    "docs/research/unified_domain_20261119/validate_h2_identity_density_slice_z8.py",
    "docs/research/unified_domain_20261119/validate_h2_identity_cr1_density_z8.py",
    "docs/research/unified_domain_20261119/validate_h2_identity_cr1_adaptive_corner_z8.py",
    "docs/research/unified_domain_20261119/validate_h2_identity_cr1_corner_z10.py",
    "docs/research/unified_domain_20261119/replay_h2_identity_cr0_reference_points.py",
    "docs/research/unified_domain_20261119/replay_h2_identity_cr1_reference_points.py",
    "docs/research/unified_domain_20261119/benchmark_h2_candidate_ab.py",
    "docs/research/unified_domain_20261119/plot_h2_identity_domain.py",
)


def main() -> None:
    payload = {
        "baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "reference_code_modified": False,
        "sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in FILES
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
