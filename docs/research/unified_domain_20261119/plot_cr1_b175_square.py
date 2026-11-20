"""绘制并核验 CR1 b175 参数方域的配对屏幕与误差、计时证据。"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
LATTICE = HERE / "screen_cr1_b175_square_lattice.json"
SOBOL = HERE / "screen_cr1_b175_square_sobol.json"
REFERENCE = HERE / "validate_cr1_b175_square_z8.json"
CLOSURE = HERE / "close_reference_cr1_b175_square_z8.json"
BENCHMARK = HERE / "benchmark_cr1_b175_wkb_guard_ab.json"
FIGURE = HERE / "cr1_b175_square_evidence.png"
MANIFEST = HERE / "cr1_b175_square_manifest.json"
REPORT = HERE / "cr1_b175_square_summary.md"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    lattice, sobol = load(LATTICE), load(SOBOL)
    reference, closure, benchmark = map(
        load, (REFERENCE, CLOSURE, BENCHMARK)
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)

    ax = axes[0]
    for dataset, marker, label in (
        (lattice, "o", "3x3 grid + interior"),
        (sobol, "^", "independent Sobol holdout"),
    ):
        for case in dataset["cases"]:
            params = case["input"]["params"]
            x, y = params["Omega_bh2"] / 1e5, params["Omega_ch2"] / 1e5
            old_ok = case["baseline_fast"]["success"]
            new_ok = case["candidate_fast"]["success"]
            color = "#238b45" if new_ok else "#cb181d"
            ax.scatter(x, y, marker=marker, s=64, color=color,
                       edgecolor="black", linewidth=0.35,
                       label=label if case is dataset["cases"][0] else None)
            if not old_ok:
                ax.scatter(x, y, marker="x", s=25, color="#cb181d", linewidth=1.2)
    ax.set(xlim=(0.7, 7.3), ylim=(0.7, 7.3),
           xlabel=r"$\Omega_{bh^2}$ ($10^5$ units)",
           ylabel=r"$\Omega_{ch^2}$ ($10^5$ units)",
           title="Paired old failure / candidate success")
    ax.grid(alpha=0.22)
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    close_by_id = {row["case_id"]: row for row in closure["cases"]}
    labels, roots, spectra = [], [], []
    for case in reference["cases"]:
        labels.append(case["case_id"].replace("square_", ""))
        roots.append(close_by_id[case["case_id"]]["fast_root_relative_to_third_map"])
        spectra.append(case["spectrum_max_abs_dex"])
    xpos = np.arange(len(labels))
    width = 0.36
    ax.bar(xpos - width/2, roots, width, label="Fast root / 3rd Reference map")
    ax.bar(xpos + width/2, spectra, width, label="Max spectrum error at Fast root")
    ax.axhline(2e-4, color="#d95f0e", linestyle="--", linewidth=1,
               label="root limit 2e-4")
    ax.axhline(1e-3, color="#756bb1", linestyle=":", linewidth=1,
               label="spectrum limit 1e-3 dex")
    ax.set_xticks(xpos, labels, rotation=18, ha="right")
    ax.set(ylabel="relative error / dex", title="Reference accuracy")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    ax = axes[2]
    names = [case["case_id"].replace("cr1_square_", "")
             for case in benchmark["cases"]]
    x = np.arange(len(names))
    p50 = [case["p50_ratio"] for case in benchmark["cases"]]
    p95 = [case["p95_ratio"] for case in benchmark["cases"]]
    ax.bar(x - width/2, p50, width, label="P50 candidate / old helper")
    ax.bar(x + width/2, p95, width, label="P95 candidate / old helper")
    ax.axhline(1.05, color="#d95f0e", linestyle="--", linewidth=1,
               label="+5% regression limit")
    ax.axhline(1.0, color="black", linewidth=0.7)
    ax.set_xticks(x, names, rotation=18, ha="right")
    ax.set(ylabel="runtime ratio", title="32 alternating full-solve pairs")
    ax.grid(axis="y", alpha=0.22)
    ax.legend(fontsize=7)

    fig.suptitle("CR1 [1e5, 7e5]^2 Fast domain screen (finite samples)")
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)

    sources = {
        name: sha256(HERE.parents[2] / "stiffgwpy_fast" / name)
        for name in ("exact_background.py", "fast_sgwb.py", "_wkb_transfer.py", "reference.py")
    }
    artifacts = [
        LATTICE, SOBOL, REFERENCE, CLOSURE, BENCHMARK, REPORT, Path(__file__),
        HERE / "screen_cr1_b175_square.py",
        HERE / "screen_h2_identity_region.py",
        HERE / "validate_cr1_b175_square_z8.py",
        HERE / "close_reference_cr1_b175_square_z8.py",
        HERE / "benchmark_cr1_b175_wkb_guard_ab.py",
        FIGURE,
    ]
    manifest = {
        "status": "finite_sample_cr1_region_evidence",
        "candidate_commit": "b175f0b3ffceccd2d046808921cb2dceffda1611",
        "source_sha256": sources,
        "reference_code_modified": False,
        "artifacts_sha256": {path.name: sha256(path) for path in artifacts},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(FIGURE)
    print(MANIFEST)


if __name__ == "__main__":
    main()
