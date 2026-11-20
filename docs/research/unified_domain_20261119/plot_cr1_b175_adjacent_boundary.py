"""Plot and hash CR1 adjacent-domain paired, precision, and performance evidence."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIGURE = HERE / "cr1_b175_adjacent_boundary_evidence.png"
MANIFEST = HERE / "cr1_b175_adjacent_manifest.json"
REPORT = HERE / "cr1_b175_adjacent_summary.md"


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    adjacent_lattice = load("screen_cr1_b175_adjacent_lattice.json")
    adjacent_sobol = load("screen_cr1_b175_adjacent_sobol.json")
    upper_lattice = load("screen_cr1_b175_upper_lattice.json")
    upper_sobol = load("screen_cr1_b175_upper_sobol.json")
    boundary = load("validate_cr1_b175_boundary_probe_z8.json")
    holdout_ref = load("validate_cr1_b175_adjacent_sobol_z8.json")
    closure = load("close_reference_cr1_b175_adjacent_z8.json")
    upper_ref = load("validate_cr1_b175_upper_z8.json")
    benchmark = load("benchmark_cr1_b175_adjacent_wkb_guard_ab_100pairs.json")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    ax = axes[0]
    for dataset, marker, color, label in (
        (upper_lattice, "o", "#9e9e9e", "upper box: Fast screen success"),
        (upper_sobol, "^", "#9e9e9e", "upper Sobol screen success"),
        (adjacent_lattice, "o", "#238b45", "accepted adjacent lattice + interior"),
        (adjacent_sobol, "^", "#238b45", "accepted independent Sobol"),
    ):
        for index, case in enumerate(dataset["cases"]):
            p = case["input"]["params"]
            ax.scatter(p["Omega_bh2"] / 1e6, p["Omega_ch2"] / 1e6,
                       marker=marker, color=color, s=48, alpha=0.78,
                       edgecolor="black", linewidth=0.3,
                       label=label if index == 0 else None)
    old_failure_labeled = False
    for dataset in (upper_lattice, upper_sobol, adjacent_lattice, adjacent_sobol):
        for case in dataset["cases"]:
            p = case["input"]["params"]
            ax.scatter(p["Omega_bh2"] / 1e6, p["Omega_ch2"] / 1e6,
                       marker="x", color="#cb181d", s=24, linewidth=1.0,
                       label="old Fast failure" if not old_failure_labeled else None)
            old_failure_labeled = True

    accepted_ref = [
        (0.75e6, 0.75e6), (0.8e6, 0.8e6),
        (784905.9135653079, 799701.8547728658),
    ]
    for index, (x, y) in enumerate(accepted_ref):
        ax.scatter(x / 1e6, y / 1e6, marker="*", color="#006d2c", s=140,
                   edgecolor="black", linewidth=0.4,
                   label="Reference accuracy pass" if index == 0 else None)
    rejected_ref = [
        (1.0e6, 1.0e6), (1.3e6, 1.3e6),
        (1.1836529400199652e6, 1.2204854907467961e6),
    ]
    for index, (x, y) in enumerate(rejected_ref):
        ax.scatter(x / 1e6, y / 1e6, marker="X", color="#b2182b", s=82,
                   edgecolor="black", linewidth=0.3,
                   label="Reference accuracy reject" if index == 0 else None)
    ax.add_patch(plt.Rectangle((0.7, 0.7), 0.1, 0.1, fill=False,
                               edgecolor="#006d2c", linewidth=1.6))
    ax.set(xlim=(0.65, 1.35), ylim=(0.65, 1.35),
           xlabel=r"$\Omega_{bh^2}$ ($10^6$ units)",
           ylabel=r"$\Omega_{ch^2}$ ($10^6$ units)",
           title="Paired screen and precision boundary")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[1]
    close_by_id = {row["case_id"]: row for row in closure["cases"]}
    boundary_by_id = {row["case_id"]: row for row in boundary["cases"]}
    holdout_case = holdout_ref["cases"][0]
    upper_by_id = {row["case_id"]: row for row in upper_ref["cases"]}
    values = [
        ("sum 1.5M", close_by_id["boundary_diagonal_1500k"]["fast_root_relative_to_third_map"], boundary_by_id["boundary_diagonal_1500k"]["spectrum_max_abs_dex"], "pass"),
        ("sum 1.6M", close_by_id["boundary_diagonal_1600k"]["fast_root_relative_to_third_map"], boundary_by_id["boundary_diagonal_1600k"]["spectrum_max_abs_dex"] if "spectrum_max_abs_dex" in close_by_id["boundary_diagonal_1600k"] else boundary_by_id["boundary_diagonal_1600k"]["spectrum_max_abs_dex"], "pass"),
        ("Sobol", close_by_id["adjacent_independent_sobol_upper_sum"]["fast_root_relative_to_third_map"], holdout_case["spectrum_max_abs_dex"], "pass"),
        ("sum 2.0M", upper_by_id["upper_square_mid_corner_2000k_sum"]["root_relative_error_proxy"], upper_by_id["upper_square_mid_corner_2000k_sum"]["spectrum_max_abs_dex"], "reject"),
        ("sum 2.6M", upper_by_id["upper_square_high_corner_2600k_sum"]["root_relative_error_proxy"], upper_by_id["upper_square_high_corner_2600k_sum"]["spectrum_max_abs_dex"], "reject"),
        ("upper Sobol", upper_by_id["upper_square_independent_sobol_holdout"]["root_relative_error_proxy"], upper_by_id["upper_square_independent_sobol_holdout"]["spectrum_max_abs_dex"], "reject"),
    ]
    x = np.arange(len(values))
    width = 0.36
    roots = [item[1] for item in values]
    spectra = [item[2] for item in values]
    colors = ["#238b45" if item[3] == "pass" else "#b2182b" for item in values]
    ax.bar(x - width/2, roots, width, color=colors, label="Fast root error")
    ax.bar(x + width/2, spectra, width, color="#f4a582", label="max spectrum error (dex)")
    ax.axhline(2e-4, color="#2166ac", linestyle="--", linewidth=1,
               label="root limit 2e-4")
    ax.axhline(1e-3, color="#762a83", linestyle=":", linewidth=1,
               label="spectrum limit 1e-3 dex")
    ax.set_xticks(x, [item[0] for item in values], rotation=25, ha="right")
    ax.set(ylabel="relative error / dex", title="Current-source z8 Reference")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=7)

    ax = axes[2]
    names = [case["case_id"].replace("adjacent_", "")
             for case in benchmark["cases"]]
    p50 = [case["p50_ratio"] for case in benchmark["cases"]]
    p95 = [case["p95_ratio"] for case in benchmark["cases"]]
    ax.bar(x[:2] - width/2, p50, width, label="P50 candidate / old helper")
    ax.bar(x[:2] + width/2, p95, width, label="P95 candidate / old helper")
    ax.axhline(1.05, color="#b2182b", linestyle="--", linewidth=1,
               label="+5% limit")
    ax.axhline(1.0, color="black", linewidth=0.7)
    ax.set_xticks(x[:2], names, rotation=18, ha="right")
    ax.set(ylabel="runtime ratio", title="100 alternating full-solve pairs")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(fontsize=7)

    fig.suptitle("CR1 adjacent density domain and rejected upper boundary")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)

    sources = {
        name: sha256(HERE.parents[2] / "stiffgwpy_fast" / name)
        for name in ("exact_background.py", "fast_sgwb.py", "_wkb_transfer.py", "reference.py")
    }
    artifact_names = [
        "screen_cr1_b175_square.py", "screen_h2_identity_region.py",
        "screen_cr1_b175_adjacent_lattice.json", "screen_cr1_b175_adjacent_sobol.json",
        "screen_cr1_b175_upper_lattice.json", "screen_cr1_b175_upper_sobol.json",
        "validate_cr1_b175_boundary_probe_z8.py", "validate_cr1_b175_boundary_probe_z8.json",
        "validate_cr1_b175_adjacent_sobol_z8.py", "validate_cr1_b175_adjacent_sobol_z8.json",
        "validate_cr1_b175_upper_z8.py", "validate_cr1_b175_upper_z8.json",
        "close_reference_cr1_b175_adjacent_z8.py", "close_reference_cr1_b175_adjacent_z8.json",
        "benchmark_cr1_b175_adjacent_wkb_guard_ab.py",
        "benchmark_cr1_b175_adjacent_wkb_guard_ab.json",
        "benchmark_cr1_b175_adjacent_wkb_guard_ab_100pairs.json",
        "benchmark_cr1_b175_wkb_guard_ab.py", "cr1_b175_adjacent_summary.md",
        Path(__file__).name, FIGURE.name,
    ]
    manifest = {
        "status": "finite_sample_cr1_adjacent_region_evidence",
        "candidate_commit": "51c90300076278057ab1fc94f4d7ca3605a82cd6",
        "baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "source_sha256": sources,
        "reference_code_modified": False,
        "artifacts_sha256": {name: sha256(HERE / name) for name in artifact_names},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(FIGURE)
    print(MANIFEST)


if __name__ == "__main__":
    main()
