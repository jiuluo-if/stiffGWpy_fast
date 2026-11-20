"""Plot and hash paired evidence for the CR0 amplitude/blue-tilt rectangle."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIGURE = HERE / "cr0_b175_amplitude_tilt_evidence.png"
MANIFEST = HERE / "cr0_b175_amplitude_tilt_manifest.json"
REPORT = HERE / "cr0_b175_amplitude_tilt_summary.md"


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    lattice = load("screen_cr0_b175_amplitude_tilt_lattice.json")
    sobol = load("screen_cr0_b175_amplitude_tilt_sobol.json")
    reference = load("validate_cr0_b175_amplitude_tilt_z8.json")
    closure = load("close_reference_cr0_b175_amplitude_tilt_z8.json")
    benchmark = load("benchmark_cr0_b175_amplitude_tilt_wkb_ab_100pairs.json")
    close_by_id = {row["case_id"]: row for row in closure["cases"]}

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    ax = axes[0]
    for dataset, marker, label in (
        (lattice, "o", "3x3 lattice + interior"),
        (sobol, "^", "independent Sobol holdout"),
    ):
        for index, case in enumerate(dataset["cases"]):
            params = case["input"]["params"]
            x, y = np.log10(params["r"]), params["n_t"]
            ax.scatter(x, y, marker=marker, color="#238b45", s=56,
                       edgecolor="black", linewidth=0.35,
                       label=label if index == 0 else None)
            if not case["baseline_fast"]["success"]:
                ax.scatter(x, y, marker="x", color="#cb181d", s=25,
                           linewidth=1.1,
                           label="old Fast failure" if index == 0 and dataset is lattice else None)
    for case in reference["cases"]:
        params = case["params"]
        ax.scatter(np.log10(params["r"]), params["n_t"], marker="*",
                   color="#006d2c", s=125, edgecolor="black", linewidth=0.35)
    ax.add_patch(plt.Rectangle((-12.0, 0.05), 4.0, 0.05, fill=False,
                               edgecolor="#006d2c", linewidth=1.6))
    ax.set(xlim=(-12.3, -7.7), ylim=(0.045, 0.105),
           xlabel=r"$\log_{10}(r)$", ylabel=r"$n_t$",
           title="Paired old failure / current success")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[1]
    refs = {case["case_id"]: case for case in reference["cases"]}
    case_order = [
        "low_amplitude_low_tilt", "low_amplitude_high_tilt",
        "high_amplitude_low_tilt", "high_amplitude_high_tilt",
        "independent_sobol_holdout",
    ]
    labels, roots, spectra = [], [], []
    for case_id in case_order:
        case = refs[case_id]
        labels.append(case_id.replace("_amplitude", "").replace("_tilt", ""))
        roots.append(close_by_id[case_id]["fast_root_relative_to_third_map"]
                     if case_id in close_by_id else case["root_relative_error_proxy"])
        spectra.append(case["spectrum_max_abs_dex"])
    xpos = np.arange(len(labels))
    width = 0.36
    ax.bar(xpos - width/2, roots, width, label="root error (map/closed map)")
    ax.bar(xpos + width/2, spectra, width, label="same-root max spectrum error")
    ax.axhline(2e-4, color="#d95f0e", linestyle="--", linewidth=1,
               label="root limit 2e-4")
    ax.axhline(1e-3, color="#756bb1", linestyle=":", linewidth=1,
               label="spectrum limit 1e-3 dex")
    ax.set_xticks(xpos, labels, rotation=24, ha="right")
    ax.set(ylabel="relative error / dex", title="Current-source z8 Reference")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    ax = axes[2]
    names = [case["case_id"].replace("cr0_", "")
             for case in benchmark["cases"]]
    x = np.arange(len(names))
    p50 = [case["p50_ratio"] for case in benchmark["cases"]]
    p95 = [case["p95_ratio"] for case in benchmark["cases"]]
    ax.bar(x - width/2, p50, width, label="P50 candidate / old helper")
    ax.bar(x + width/2, p95, width, label="P95 candidate / old helper")
    ax.axhline(1.05, color="#d95f0e", linestyle="--", linewidth=1,
               label="+5% limit")
    ax.axhline(1.0, color="black", linewidth=0.7)
    ax.set_xticks(x, names, rotation=16, ha="right")
    ax.set(ylabel="runtime ratio", title="100 alternating full-solve pairs")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    fig.suptitle("CR0 small-amplitude, positive-tilt Fast domain (finite samples)")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)

    source_hashes = {
        "exact_background.py": reference["source_sha256"]["exact_background.py"],
        "fast_sgwb.py": reference["source_sha256"]["fast_sgwb.py"],
        "_wkb_transfer.py": benchmark["candidate_helper_sha256"],
        "reference.py": reference["source_sha256"]["reference.py"],
    }
    artifact_paths = {
        name: HERE / name for name in (
            "screen_cr0_b175_amplitude_tilt.py",
            "screen_h2_identity_region.py",
            "screen_cr0_b175_amplitude_tilt_lattice.json",
            "screen_cr0_b175_amplitude_tilt_sobol.json",
            "validate_cr0_b175_amplitude_tilt_z8.py",
            "validate_cr0_b175_amplitude_tilt_z8.json",
            "close_reference_cr0_b175_amplitude_tilt_z8.py",
            "close_reference_cr0_b175_amplitude_tilt_z8.json",
            "benchmark_cr0_b175_amplitude_tilt_wkb_ab.py",
            "benchmark_cr0_b175_amplitude_tilt_wkb_ab_100pairs.json",
            "benchmark_cr1_b175_wkb_guard_ab.py",
            "cr0_b175_amplitude_tilt_summary.md",
            Path(__file__).name,
            FIGURE.name,
        )
    }
    manifest = {
        "status": "finite_sample_cr0_amplitude_tilt_region_evidence",
        "baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "candidate_commit": "c3c6ab3fc34e705ebdea3d16739740a227d1acb7",
        "source_hashes_referenced_from_raw_validation": True,
        "source_sha256": source_hashes,
        "reference_code_modified": False,
        "artifacts_sha256": {
            name: sha256(path) for name, path in artifact_paths.items()
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(FIGURE)
    print(MANIFEST)


if __name__ == "__main__":
    main()
