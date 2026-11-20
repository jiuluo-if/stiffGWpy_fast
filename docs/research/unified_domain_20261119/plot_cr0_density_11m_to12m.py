"""Plot finite-sample CR0 11M–12M coverage, accuracy and runtime evidence."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = load("screen_cr0_density_11m_to12m.json")["cases"]
    sobol = load("screen_cr0_density_11m_to12m_sobol.json")["cases"]
    reference = load("validate_cr0_density_11m_to12m_z8.json")["cases"]
    timing = load("benchmark_cr0_density_11m_to12m_self_ab_tbb_500pairs.json")["cases"]
    rows = lattice + sobol
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    fig.suptitle("CR0 density continuation: 11M–12M (finite sampled evidence)", fontsize=15)

    ax = axes[0, 0]
    ax.fill_between([11.0, 12.0], -0.5, 0.0, color="#d9ead3", alpha=0.7, label="screened box")
    for row in rows:
        result = row["input"]["params"]
        old = row["baseline_fast"]["success"]
        ax.scatter(result["Omega_bh2"] / 1e6, result["n_t"], marker="s" if old else "x",
                   s=42, color="#cb181d" if not old else "#3182bd", zorder=3)
        if row["candidate_fast"]["success"] and not old:
            ax.scatter(result["Omega_bh2"] / 1e6, result["n_t"], marker="o", s=23,
                       color="#238b45", zorder=4)
    for row in reference:
        result = row["params"]
        ax.scatter(result["Omega_bh2"] / 1e6, result["n_t"], marker="*", s=125,
                   color="#ffbf00", edgecolor="black", zorder=5)
    ax.set(xlabel=r"$\Omega_{bh2}$ ($10^6$ units)", ylabel=r"$n_t$",
           title="Paired old/current Fast; ★ current-source Reference")
    ax.set_xlim(10.92, 12.08)
    ax.set_ylim(-0.56, 0.06)
    ax.grid(alpha=.25)
    ax.text(.02, .03, "Old Fast 0/19; all failures: nonfinite_background_primitive",
            transform=ax.transAxes, fontsize=8, va="bottom")

    labels = ["12M −.5", "12M 0", "Sobol"]
    x = np.arange(len(labels))
    ax = axes[0, 1]
    ax.bar(x, [row["root_relative_error"] for row in reference], color="#3182bd")
    ax.axhline(2e-4, color="#cb181d", linestyle="--", label="root gate 2e−4")
    ax.set_xticks(x, labels)
    ax.set(ylabel="relative root error", title="Correct-origin z8 Reference closure")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.bar(x, [row["spectrum_max_abs_dex"] for row in reference], color="#756bb1")
    ax.axhline(1e-3, color="#cb181d", linestyle="--", label="spectrum gate 1e−3 dex")
    ax.set_xticks(x, labels)
    ax.set(ylabel="maximum same-root spectrum error (dex)", title="Shared Fast root and frequency grid")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    timing_labels = ["12M corner", "Sobol", "trusted CR0"]
    width = .36
    ax.bar(x-width/2, [row["p50_ratio"] for row in timing], width, label="P50", color="#31a354")
    ax.bar(x+width/2, [row["p95_ratio"] for row in timing], width, label="P95", color="#a1d99b")
    ax.axhline(1.05, color="#cb181d", linestyle="--", label="5% regression limit")
    ax.set_xticks(x, timing_labels, rotation=15, ha="right")
    ax.set(ylabel="current / source baseline", title="500 paired complete Fast solves")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)
    output = HERE / "cr0_density_11m_to12m_domain.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
