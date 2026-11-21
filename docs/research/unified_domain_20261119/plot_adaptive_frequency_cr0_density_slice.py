"""Plot adaptive-vs-current-goal Fast results for the CR0 density slice."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE / "screen_cr0_density_slice_adaptive_frequency.json"
OUTPUT = HERE / "adaptive_frequency_cr0_density_slice.png"


def main() -> None:
    cases = json.loads(INPUT.read_text(encoding="utf-8"))["cases"]
    density = np.asarray([row["params"]["Omega_bh2"] / 1.0e6 for row in cases])
    tilt = np.asarray([row["params"]["n_t"] for row in cases])
    correction = 100.0 * np.asarray(
        [row["adaptive_root_relative_correction_vs_current_goal"] for row in cases]
    )
    root_error = 1.0e4 * np.asarray(
        [row["root_relative_error_vs_dense_fast"] for row in cases]
    )
    spectrum_error = 1.0e6 * np.asarray(
        [row["max_native_spectrum_error_dex"] for row in cases]
    )
    holdout = np.asarray([row["family"] == "holdout" for row in cases])
    labels = [row["case_id"] for row in cases]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
    scatter = axes[0].scatter(
        density[~holdout], tilt[~holdout], c=correction[~holdout], s=100,
        marker="o", cmap="viridis", edgecolor="black", linewidth=0.6,
        label="rectangle samples",
    )
    axes[0].scatter(
        density[holdout], tilt[holdout], c=correction[holdout], s=140,
        marker="^", cmap="viridis", edgecolor="black", linewidth=0.8,
        label="independent holdouts",
    )
    colorbar = fig.colorbar(scatter, ax=axes[0])
    colorbar.set_label("Adaptive root correction vs current goal (%)")
    axes[0].set_xlabel(r"$\Omega_b h^2$ ($10^6$ units)")
    axes[0].set_ylabel(r"$n_t$")
    axes[0].set_title("CR0 density–tilt sample coverage")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="best", frameon=True)

    x = np.arange(len(cases))
    for i, is_holdout in enumerate(holdout):
        axes[1].scatter(
            i, root_error[i], marker="^" if is_holdout else "o",
            color="#2463a6", s=58, zorder=3,
        )
    axes[1].axhline(2.0, color="#2463a6", linestyle="--", linewidth=1,
                    label=r"root gate $2\times10^{-4}$")
    axes[1].set_ylabel(r"Root error vs dense Fast ($10^{-4}$)")
    axes[1].set_xticks(x, labels, rotation=55, ha="right", fontsize=8)
    axes[1].grid(axis="y", alpha=0.25)
    right = axes[1].twinx()
    right.plot(x, spectrum_error, color="#bc4b36", marker="s", markersize=4,
               linewidth=1, label="native spectrum difference")
    right.axhline(1000, color="#bc4b36", linestyle=":", linewidth=1,
                  label=r"spectrum gate $10^{-3}$ dex")
    right.set_ylabel("Max native-node spectrum difference ($10⁻⁶ dex)",
                     color="#bc4b36")
    right.tick_params(axis="y", colors="#bc4b36")
    lines, names = axes[1].get_legend_handles_labels()
    lines2, names2 = right.get_legend_handles_labels()
    axes[1].legend(lines + lines2, names + names2, loc="upper left", fontsize=8)
    axes[1].set_title("Same-input Fast support checks")
    fig.suptitle(
        "Adaptive frequency support: corrected roots, not domain certification",
        fontsize=12,
    )
    fig.savefig(OUTPUT, dpi=180)
    print(OUTPUT)


if __name__ == "__main__":
    main()
