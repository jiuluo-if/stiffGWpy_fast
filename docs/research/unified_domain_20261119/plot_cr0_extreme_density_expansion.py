"""Plot paired Fast status and independent accuracy for the two new slices."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DATA = Path(__file__).resolve().parent


def load(name: str) -> dict:
    return json.loads((DATA / name).read_text(encoding="utf-8"))


def _slice(ax, data: dict, dim_y: str, title: str, oracle: dict) -> None:
    rows = data["cases"]
    for row in rows:
        params = row["input"]["params"]
        x = np.log10(params["Omega_bh2"])
        y = params[dim_y]
        old = row["baseline_fast"]
        new = row["candidate_fast"]
        ax.scatter(x, y, marker="x", s=48, color="#c0392b", linewidth=1.5,
                   label="HEAD failed" if row is rows[0] else None, zorder=2)
        if new["success"]:
            ax.scatter(x, y, marker="o", s=28, color="#238b45", alpha=0.8,
                       label="candidate solved" if row is rows[0] else None, zorder=3)
    representative = oracle["inputs"]
    ax.scatter(np.log10(representative["Omega_bh2"]), representative[dim_y],
               marker="*", s=240, facecolor="#f1c40f", edgecolor="#333333",
               linewidth=0.8, label="full independent oracle", zorder=4)
    ax.set_title(title)
    ax.set_xlabel(r"$\log_{10}(\Omega_b h^2)$")
    ax.set_ylabel("H0" if dim_y == "H0" else "T_re (GeV)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", frameon=True, fontsize=8)


def main() -> None:
    h0_screen = load("screen_cr0_extreme_density_h0.json")
    re_screen = load("screen_cr0_extreme_density_reheat.json")
    h0_oracle = load("cr0_extreme_density_x_oracle.json")
    re_oracle = load("cr0_extreme_density_reheat_x_oracle.json")
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), constrained_layout=True)
    _slice(axes[0], h0_screen, "H0", "CR0 slice A: density × H0", h0_oracle)
    _slice(axes[1], re_screen, "T_re", "CR0 slice B: density × reheating temperature", re_oracle)
    labels = ["density × H0", "density × T_re"]
    roots = [h0_oracle["comparison"]["root_relative_error"] / 2e-4,
             re_oracle["comparison"]["root_relative_error"] / 2e-4]
    spectra = [h0_oracle["comparison"]["spectrum_max_error_dex"] / 1e-3,
               re_oracle["comparison"]["spectrum_max_error_dex"] / 1e-3]
    pos = np.arange(2)
    axes[2].bar(pos - 0.18, roots, width=0.36, color="#2471a3",
                label="root relative error")
    axes[2].bar(pos + 0.18, spectra, width=0.36, color="#8e44ad",
                label="max spectrum error (dex)")
    axes[2].axhline(1.0, color="#333333", linestyle="--", linewidth=1,
                    label="acceptance limit")
    axes[2].set_xticks(pos, labels, rotation=15)
    axes[2].set_ylabel("fraction of acceptance limit")
    axes[2].set_title("Independent full-grid precision")
    axes[2].grid(True, axis="y", alpha=0.25)
    axes[2].legend(fontsize=8, loc="best")
    fig.suptitle("CR0 extreme-density expansion: paired Fast failures and precision checks", fontsize=13)
    output = DATA / "cr0_extreme_density_expansion.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
