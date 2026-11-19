"""Plot paired Fast coverage and direct accuracy for the CR0 tilt slice."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def _read(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = _read("screen_cr0_density_negative_tilt_exploratory.json")["cases"]
    sobol = _read("screen_cr0_density_negative_tilt_sobol.json")["cases"]
    reference = _read("validate_cr0_density_negative_tilt_z8.json")["cases"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    ax = axes[0]
    for row in lattice:
        params = row["input"]["params"]
        old = row["baseline_fast"]
        new = row["candidate_fast"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="s", s=75,
                   facecolor="#d95f02" if not old["success"] else "#1b9e77",
                   edgecolor="#222222", linewidth=0.6, zorder=3)
        if new["success"]:
            ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="o", s=29,
                       facecolor="#1b9e77", edgecolor="white", linewidth=0.5, zorder=4)
    for row in sobol:
        params = row["input"]["params"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="x", s=47,
                   color="#386cb0", linewidth=1.7, zorder=5)
    ax.set(xlabel=r"$\Omega_{bh2}$ ($10^6$ units)", ylabel=r"$n_t$",
           title="Paired Fast coverage")
    ax.grid(alpha=0.25)
    ax.text(0.02, 0.02, "square: lattice old failure\ncircle: candidate success\nx: 8-point scrambled Sobol",
            transform=ax.transAxes, fontsize=8, va="bottom")

    ax = axes[1]
    old_failures: dict[str, int] = {}
    candidate_success = 0
    for row in [*lattice, *sobol]:
        old = row["baseline_fast"]
        new = row["candidate_fast"]
        if not old["success"]:
            reason = old.get("failure_reason") or "unspecified"
            old_failures[reason] = old_failures.get(reason, 0) + 1
        candidate_success += int(new["success"])
    labels = list(old_failures) + ["candidate success"]
    values = list(old_failures.values()) + [candidate_success]
    colors = ["#d95f02"] * len(old_failures) + ["#1b9e77"]
    ax.barh(labels, values, color=colors)
    ax.invert_yaxis()
    ax.set(xlabel="sample count", title="Failure type and recovery")
    ax.grid(axis="x", alpha=0.25)
    for y, value in enumerate(values):
        ax.text(value + 0.15, y, str(value), va="center", fontsize=9)

    ax = axes[2]
    root = np.asarray([row["root_relative_error"] for row in reference])
    spectrum = np.asarray([row["spectrum_max_abs_dex"] for row in reference])
    labels = [row["case_id"].replace("independent_sobol_", "holdout_") for row in reference]
    y = np.arange(len(labels))
    ax.scatter(root / 2e-4, y - 0.12, label="root error / root gate", marker="o", color="#386cb0")
    ax.scatter(spectrum / 1e-3, y + 0.12, label="max spectrum error / spectrum gate", marker="s", color="#7b3294")
    ax.axvline(1.0, color="#222222", linestyle="--", linewidth=1, label="acceptance limit")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set(xlabel="fraction of respective acceptance limit", title="Direct z8 Reference errors")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")

    fig.suptitle(r"CR0 sampled expansion: $\Omega_{bh2}\in[10^5,10^6]$, $n_t\in[-0.5,0]$", fontsize=12)
    fig.savefig(HERE / "cr0_density_negative_tilt_domain.png", dpi=180)


if __name__ == "__main__":
    main()
