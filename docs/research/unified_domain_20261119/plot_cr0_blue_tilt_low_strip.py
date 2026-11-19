"""Plot the accepted CR0 low-positive-tilt sampled expansion."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def _read(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = _read("screen_cr0_blue_tilt_low_strip.json")["cases"]
    sobol = _read("screen_cr0_blue_tilt_low_strip_sobol.json")["cases"]
    refs = _read("validate_cr0_blue_tilt_low_strip_reference.json")["cases"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    ax = axes[0]
    for row in lattice:
        params = row["input"]["params"]
        old, new = row["baseline_fast"], row["candidate_fast"]
        x = np.log10(params["r"])
        y = params["n_t"]
        ax.scatter(x, y, marker="s", s=73,
                   facecolor="#d95f02" if not old["success"] else "#1b9e77",
                   edgecolor="#222222", linewidth=0.6, zorder=3)
        if new["success"]:
            ax.scatter(x, y, marker="o", s=27, facecolor="#1b9e77",
                       edgecolor="white", linewidth=0.5, zorder=4)
    for row in sobol:
        params = row["input"]["params"]
        ax.scatter(np.log10(params["r"]), params["n_t"], marker="x",
                   color="#386cb0", s=46, linewidth=1.6, zorder=5)
    ax.scatter([-10.0, -10.0, -12.0], [0.2, 0.3, 0.5], marker="x",
               color="#b2182b", s=55, linewidth=1.7)
    ax.set(xlabel=r"$\log_{10}(r)$", ylabel=r"$n_t$",
           title="Paired Fast coverage")
    ax.grid(alpha=0.25)
    ax.text(0.98, 0.98,
            "square: old Fast failure\ncircle: candidate success\nblue x: scrambled Sobol\nred x: Reference root-gate failure",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})

    ax = axes[1]
    rows = [*lattice, *sobol]
    old_fail: dict[str, int] = {}
    new_success = 0
    for row in rows:
        old, new = row["baseline_fast"], row["candidate_fast"]
        if not old["success"]:
            reason = old.get("failure_reason") or "unspecified"
            old_fail[reason] = old_fail.get(reason, 0) + 1
        new_success += int(new["success"])
    labels = list(old_fail) + ["candidate success"]
    values = list(old_fail.values()) + [new_success]
    colors = ["#d95f02"] * len(old_fail) + ["#1b9e77"]
    ax.barh(labels, values, color=colors)
    ax.invert_yaxis()
    ax.set(xlabel="sample count", title="Failure type and recovery")
    ax.grid(axis="x", alpha=0.25)
    for y, value in enumerate(values):
        ax.text(value + 0.15, y, str(value), va="center", fontsize=9)

    ax = axes[2]
    root = np.asarray([row["root_relative_error"] / 2e-4 for row in refs])
    spectrum = np.asarray([row["spectrum_max_abs_dex"] / 1e-3 for row in refs])
    labels = [row["case_id"].replace("_", " ") for row in refs]
    y = np.arange(len(labels))
    ax.scatter(root, y - 0.12, label="root error / gate", marker="o", color="#386cb0")
    ax.scatter(spectrum, y + 0.12, label="spectrum error / gate", marker="s", color="#7b3294")
    ax.axvline(1.0, color="#222222", linestyle="--", linewidth=1, label="acceptance limit")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set(xlabel="fraction of respective acceptance limit", title="Direct z8 Reference errors")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")

    fig.suptitle(r"CR0 sampled Fast expansion: $\log_{10}(r)\in[-12,-8]$, $n_t\in[0.05,0.1]$", fontsize=12)
    fig.savefig(HERE / "cr0_blue_tilt_low_strip_domain.png", dpi=180)


if __name__ == "__main__":
    main()
