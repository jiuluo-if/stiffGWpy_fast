"""Plot paired CR1 reheat/kappa coverage, Reference accuracy and A/B timing."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIGURE = HERE / "cr1_b175_reheat_kappa_evidence.png"
MANIFEST = HERE / "cr1_b175_reheat_kappa_manifest.json"
REPORT = HERE / "cr1_b175_reheat_kappa_summary.md"


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    lattice, sobol = load("screen_cr1_b175_reheat_kappa_lattice.json"), load("screen_cr1_b175_reheat_kappa_sobol.json")
    reference = load("validate_cr1_b175_reheat_kappa_z8.json")
    closure = load("close_reference_cr1_b175_reheat_kappa_z8.json")
    base_replay = load("replay_cr1_b175_reheat_kappa_base_corner.json")
    old_base_ref = load("validate_cr1_b175_square_z8.json")
    old_base_closure = load("close_reference_cr1_b175_square_z8.json")
    benchmark = load("benchmark_cr1_b175_reheat_kappa_wkb_ab_100pairs.json")
    ref_by_id = {row["case_id"]: row for row in reference["cases"]}
    close_by_id = {row["case_id"]: row for row in closure["cases"]}
    old_ref = next(row for row in old_base_ref["cases"] if row["case_id"] == "high_square_corner")
    old_close = next(row for row in old_base_closure["cases"] if row["case_id"] == "high_square_corner")

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    ax = axes[0]
    for dataset, marker, label in (
        (lattice, "o", "3x3 lattice + interior"),
        (sobol, "^", "independent Sobol holdout"),
    ):
        for index, case in enumerate(dataset["cases"]):
            p = case["input"]["params"]
            ax.scatter(p["T_re"], p["kappa10"], marker=marker,
                       color="#238b45", edgecolor="black", linewidth=0.35,
                       s=56, label=label if index == 0 else None)
            if not case["baseline_fast"]["success"]:
                ax.scatter(p["T_re"], p["kappa10"], marker="x",
                           color="#cb181d", s=24, linewidth=1.0,
                           label="old Fast failure" if index == 0 and dataset is lattice else None)
    ref_points = [
        (2000.0, 0.01), (2000.0, 0.015), (2500.0, 0.01),
        (2500.0, 0.015), (2372.528196312487, 0.014170661885291338),
    ]
    for index, (t_re, kappa) in enumerate(ref_points):
        ax.scatter(t_re, kappa, marker="*", color="#006d2c", s=125,
                   edgecolor="black", linewidth=0.35,
                   label="Reference-checked" if index == 0 else None)
    ax.add_patch(plt.Rectangle((2000.0, 0.01), 500.0, 0.005, fill=False,
                               edgecolor="#006d2c", linewidth=1.6))
    ax.set(xlim=(1970, 2530), ylim=(0.0095, 0.0155),
           xlabel=r"$T_{re}$", ylabel=r"$\kappa_{10}$",
           title="Paired old failure / current success")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8, loc="lower right")

    labels = ["base corner", "low T / high k", "high T / low k",
              "high T / high k", "Sobol holdout"]
    root_errors = [
        old_close["fast_root_relative_to_third_map"],
        close_by_id["low_reheat_high_kappa"]["fast_root_relative_to_third_map"],
        ref_by_id["high_reheat_low_kappa"]["root_relative_error_proxy"],
        ref_by_id["high_reheat_high_kappa"]["root_relative_error_proxy"],
        close_by_id["independent_sobol_holdout"]["fast_root_relative_to_third_map"],
    ]
    spectra = [
        old_ref["spectrum_max_abs_dex"],
        ref_by_id["low_reheat_high_kappa"]["spectrum_max_abs_dex"],
        ref_by_id["high_reheat_low_kappa"]["spectrum_max_abs_dex"],
        ref_by_id["high_reheat_high_kappa"]["spectrum_max_abs_dex"],
        ref_by_id["independent_sobol_holdout"]["spectrum_max_abs_dex"],
    ]
    ax = axes[1]
    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x-width/2, root_errors, width, label="root error")
    ax.bar(x+width/2, spectra, width, label="max spectrum error")
    ax.axhline(2e-4, color="#d95f0e", linestyle="--", linewidth=1,
               label="root limit 2e-4")
    ax.axhline(1e-3, color="#756bb1", linestyle=":", linewidth=1,
               label="spectrum limit 1e-3 dex")
    ax.set_xticks(x, labels, rotation=20, ha="right", fontsize=7)
    ax.set(ylabel="relative error / dex", title="Current-source z8 Reference")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    ax = axes[2]
    names = [case["case_id"].replace("cr1_reheat_kappa_", "")
             for case in benchmark["cases"]]
    p50 = [case["p50_ratio"] for case in benchmark["cases"]]
    p95 = [case["p95_ratio"] for case in benchmark["cases"]]
    positions = np.arange(len(names))
    ax.bar(positions-width/2, p50, width, label="P50 candidate / old helper")
    ax.bar(positions+width/2, p95, width, label="P95 candidate / old helper")
    ax.axhline(1.05, color="#d95f0e", linestyle="--", linewidth=1,
               label="+5% limit")
    ax.axhline(1.0, color="black", linewidth=0.7)
    ax.set_xticks(positions, names, rotation=14, ha="right")
    ax.set(ylabel="runtime ratio", title="100 alternating full-solve pairs")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    fig.suptitle("CR1 reheating / kappa10 Fast domain (finite samples)")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)

    source_hashes = {
        name: reference["source_sha256"][name]
        for name in ("exact_background.py", "fast_sgwb.py", "reference.py")
    }
    source_hashes["_wkb_transfer.py"] = benchmark["candidate_helper_sha256"]
    source_hashes["freq_adaptive.py"] = base_replay["current_source_sha256"]["freq_adaptive.py"]
    artifact_names = [
        "screen_cr1_b175_reheat_kappa.py", "screen_h2_identity_region.py",
        "screen_cr1_b175_reheat_kappa_lattice.json", "screen_cr1_b175_reheat_kappa_sobol.json",
        "validate_cr1_b175_reheat_kappa_z8.py", "validate_cr1_b175_reheat_kappa_z8.json",
        "close_reference_cr1_b175_reheat_kappa_z8.py", "close_reference_cr1_b175_reheat_kappa_z8.json",
        "replay_cr1_b175_reheat_kappa_base_corner.py", "replay_cr1_b175_reheat_kappa_base_corner.json",
        "validate_cr1_b175_square_z8.py", "validate_cr1_b175_square_z8.json",
        "close_reference_cr1_b175_square_z8.py", "close_reference_cr1_b175_square_z8.json",
        "benchmark_cr1_b175_reheat_kappa_wkb_ab.py",
        "benchmark_cr1_b175_reheat_kappa_wkb_ab_100pairs.json",
        "benchmark_cr1_b175_wkb_guard_ab.py", REPORT.name, Path(__file__).name, FIGURE.name,
    ]
    manifest = {
        "status": "finite_sample_cr1_reheat_kappa_evidence",
        "baseline_commit": "872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
        "candidate_commit": "8e71fb74cadc50c1e3493a877cd19975a47e984b",
        "reference_code_modified": False,
        "source_sha256": source_hashes,
        "artifacts_sha256": {name: sha(HERE/name) for name in artifact_names},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(FIGURE)
    print(MANIFEST)


if __name__ == "__main__":
    main()
