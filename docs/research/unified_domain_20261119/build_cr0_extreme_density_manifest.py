"""Build source/data hashes and the compact acceptance summary."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    artifacts = (
        "cr0_extreme_density_closure_report.md",
        "cr0_extreme_density_expansion.png",
        "screen_cr0_extreme_density_h0.json",
        "screen_cr0_extreme_density_reheat.json",
        "cr0_extreme_density_x_oracle.json",
        "cr0_extreme_density_reheat_x_oracle.json",
        "benchmark_extreme_density_closure_trusted_ab.json",
        "screen_cr0_extreme_density_slices.py",
        "verify_cr0_extreme_density_x_oracle.py",
        "benchmark_extreme_density_closure_trusted_ab.py",
        "plot_cr0_extreme_density_expansion.py",
        "build_cr0_extreme_density_manifest.py",
    )
    h0 = json.loads((HERE / "screen_cr0_extreme_density_h0.json").read_text())
    reheat = json.loads((HERE / "screen_cr0_extreme_density_reheat.json").read_text())
    oracle_h0 = json.loads((HERE / "cr0_extreme_density_x_oracle.json").read_text())
    oracle_reheat = json.loads((HERE / "cr0_extreme_density_reheat_x_oracle.json").read_text())
    paired = {}
    for key, data in (("density_h0", h0), ("density_reheat", reheat)):
        old = [row["baseline_fast"] for row in data["cases"]]
        new = [row["candidate_fast"] for row in data["cases"]]
        paired[key] = {
            "case_count": len(old),
            "head_successes": sum(bool(row["success"]) for row in old),
            "candidate_successes": sum(bool(row["success"]) for row in new),
            "head_failure_reasons": {
                reason: sum(row["failure_reason"] == reason for row in old)
                for reason in sorted({row["failure_reason"] for row in old})
            },
        }
    payload = {
        "status": "validated_local_candidate",
        "baseline_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            name: sha256(ROOT / "stiffgwpy_fast" / name)
            for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
        },
        "test_sha256": sha256(ROOT / "tests" / "test_physics_limits.py"),
        "artifact_sha256": {name: sha256(HERE / name) for name in artifacts},
        "paired_screens": paired,
        "independent_full_grid_checks": {
            "density_h0": oracle_h0["comparison"],
            "density_reheat": oracle_reheat["comparison"],
        },
        "trusted_ab": json.loads(
            (HERE / "benchmark_extreme_density_closure_trusted_ab.json").read_text()
        ),
        "scope_note": "Finite sample expansion evidence only; no continuous accuracy certification.",
    }
    output = HERE / "cr0_extreme_density_expansion_manifest.json"
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
