"""Plot finite-sample coverage, accuracy, Reference limits and runtime."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = load("screen_cr0_density_14m_upper_nearzero015_slice.json")["cases"]
    sobol = load("screen_cr0_density_14m_upper_nearzero015_slice_sobol.json")["cases"]
    reference = load("validate_cr0_density_14m_nearzero015_slice_z8.json")["cases"]
    holdout = load("probe_cr0_density_14m_upper_holdout05_reference_map.json")["case"]
    timing = load("benchmark_cr0_density_14m_upper_nearzero015_self_ab_tbb_500pairs.json")["cases"]
    rows = lattice + sobol
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)
    fig.suptitle("CR0 high-density edge: 13.875M–14M, near-zero tilt (sampled evidence)", fontsize=14)

    ax = axes[0, 0]
    ax.fill_between([13.875, 14.0], -0.015625, 0.0, color="#d9ead3", alpha=0.75,
                    label="screened box")
    for row in rows:
        params = row["input"]["params"]
        old_ok = row["baseline_fast"]["success"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="s" if old_ok else "x",
                   s=43, color="#cb181d" if not old_ok else "#3182bd", zorder=3)
        if row["candidate_fast"]["success"] and not old_ok:
            ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="o", s=22,
                       color="#238b45", zorder=4)
    for row in reference:
        params = row["params"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="*", s=115,
                   color="#ffbf00", edgecolor="black", zorder=5)
    holdout_params = holdout["params"]
    ax.scatter(holdout_params["Omega_bh2"] / 1e6, holdout_params["n_t"], marker="*", s=115,
               color="#6a51a3", edgecolor="black", zorder=5, label="Sobol Reference map")
    ax.scatter(13.75, -0.015625, marker="v", s=70, color="#d7301f", zorder=6,
               label="Reference step failure; excluded")
    ax.set(xlabel=r"$\Omega_{bh2}$ ($10^6$ units)", ylabel=r"$n_t$",
           title="Old Fast failures ×; new Fast success ○; Reference ★")
    ax.set_xlim(13.72, 14.015)
    ax.set_ylim(-0.018, 0.002)
    ax.grid(alpha=.25)
    ax.legend(fontsize=7, loc="lower right")
    ax.text(.02, .96, "Old Fast 0/19; current Fast 19/19\nOld failure: nonfinite_background_primitive",
            transform=ax.transAxes, fontsize=8, va="top")

    samples = reference + [{
        "case_id": "independent Sobol (single converged map)",
        "root_relative_error": holdout["root_relative_error_at_fast_input"],
        "spectrum_max_abs_dex": holdout["spectrum_max_abs_dex"],
    }]
    x = np.arange(len(samples))
    labels = ["low-density edge", "14M edge", "interior", "Sobol"]
    ax = axes[0, 1]
    ax.bar(x, [r["root_relative_error"] for r in samples], color="#3182bd")
    ax.axhline(2e-4, color="#cb181d", linestyle="--", label="root gate 2e−4")
    ax.set_xticks(x, labels, rotation=12, ha="right")
    ax.set(ylabel="relative root error / fixed-point residual", title="Reference root checks")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.bar(x, [r["spectrum_max_abs_dex"] for r in samples], color="#756bb1")
    ax.axhline(1e-3, color="#cb181d", linestyle="--", label="spectrum gate 1e−3 dex")
    ax.set_xticks(x, labels, rotation=12, ha="right")
    ax.set(ylabel="maximum spectrum error (dex)",
           title="75-frequency Reference comparison; extra Sobol closure map failed")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    labels = ["upper edge", "interior", "Sobol", "trusted"]
    x = np.arange(len(labels))
    width = .36
    ax.bar(x-width/2, [r["p50_ratio"] for r in timing], width, label="P50", color="#31a354")
    ax.bar(x+width/2, [r["p95_ratio"] for r in timing], width, label="P95", color="#a1d99b")
    ax.axhline(1.05, color="#cb181d", linestyle="--", label="5% regression limit")
    ax.set_xticks(x, labels)
    ax.set(ylabel="candidate / Fast baseline", title="500 paired full Fast solves")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    output = HERE / "cr0_density_14m_upper_nearzero015_slice.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
