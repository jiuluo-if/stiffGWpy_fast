"""Plot finite-sample CR1 upper-kappa domain and its accuracy/runtime gates."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


grid = load("screen_cr1_b175_reheat_kappa_upper_lattice.json")["cases"]
sobol = load("screen_cr1_b175_reheat_kappa_upper_sobol.json")["cases"]
checked = load("validate_cr1_b175_reheat_kappa_upper_z8.json")["cases"]
reuse_low_high = load("validate_cr1_b175_reheat_kappa_narrow_z8.json")["cases"]
reuse_high_high = load("validate_cr1_b175_reheat_kappa_z8.json")["cases"]
old_success = load("audit_cr1_upper_old_success_point_z8.json")
benchmark = load("benchmark_cr1_b175_reheat_kappa_upper_wkb_100pairs.json")["cases"]

fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
fig.suptitle("CR1 upper-kappa reheating extension (finite samples)", fontsize=16)

ax = axes[0, 0]
ax.fill_between([2000, 2500], 0.01, 0.015, color="#6baed6", alpha=0.2,
                label="previously checked CR1 box")
ax.fill_between([1750, 2000], 0.0125, 0.015, color="#74c476", alpha=0.27,
                label="candidate extension")
all_cases = grid + sobol
for item in all_cases:
    p = item["input"]["params"]
    old, new = item["baseline_fast"], item["candidate_fast"]
    if old["success"]:
        ax.scatter(p["T_re"], p["kappa10"], marker="s", s=65,
                   facecolor="#3182bd", edgecolor="black", zorder=5,
                   label="old and current success" if item is next(
                       x for x in all_cases if x["baseline_fast"]["success"]
                   ) else None)
    else:
        ax.scatter(p["T_re"], p["kappa10"], marker="x", s=35,
                   color="#cb181d", linewidth=1.6, zorder=4,
                   label="old failure" if item is grid[0] else None)
        if new["success"]:
            ax.scatter(p["T_re"], p["kappa10"], marker="o", s=25,
                       color="#238b45", zorder=3,
                       label="current success" if item is grid[0] else None)
checked_points = [
    (1750.0, 0.0125),
    (2000.0, 0.0125),
    (float(checked[2]["params"]["T_re"]), float(checked[2]["params"]["kappa10"])),
    (1750.0, 0.015),
    (2000.0, 0.015),
    (old_success["params"]["T_re"], old_success["params"]["kappa10"]),
]
for index, (t_re, kappa) in enumerate(checked_points):
    ax.scatter(t_re, kappa, marker="*", s=105, facecolor="#ffbf00",
               edgecolor="black", zorder=6,
               label="Reference checked" if index == 0 else None)
ax.set(xlabel=r"$T_{re}$", ylabel=r"$\kappa_{10}$",
       title="Paired old/current Fast coverage")
ax.grid(alpha=0.25)
ax.legend(fontsize=8, loc="lower right")

root_labels = ["1750/.0125", "2000/.0125", "Sobol", "1750/.015",
               "2000/.015", "old success\n(current)"]
root_errors = [
    checked[0]["root_relative_error_proxy"],
    checked[1]["root_relative_error_proxy"],
    checked[2]["root_relative_error_proxy"],
    reuse_low_high[1]["root_relative_error_proxy"],
    reuse_high_high[0]["root_relative_error_proxy"],
    old_success["candidate_root_relative_error_proxy"],
]
ax = axes[0, 1]
positions = np.arange(len(root_labels))
ax.bar(positions, root_errors, color="#3182bd")
ax.axhline(2e-4, color="#cb181d", linestyle="--", label="2e-4 gate")
ax.set_yscale("log")
ax.set_xticks(positions, root_labels, rotation=20, ha="right", fontsize=8)
ax.set_ylabel("relative root difference")
ax.set_title("Candidate vs current-source Reference")
ax.grid(axis="y", alpha=0.25)
ax.legend(fontsize=8)

spectrum_labels = ["1750/.0125", "2000/.0125", "Sobol", "1750/.015",
                   "2000/.015", "old success\n(current)"]
spectrum_errors = [
    checked[0]["spectrum_max_abs_dex"],
    checked[1]["spectrum_max_abs_dex"],
    checked[2]["spectrum_max_abs_dex"],
    reuse_low_high[1]["spectrum_max_abs_dex"],
    reuse_high_high[0]["spectrum_max_abs_dex"],
    old_success["spectrum_max_abs_dex"],
]
ax = axes[1, 0]
positions = np.arange(len(spectrum_labels))
ax.bar(positions, spectrum_errors, color="#fd8d3c")
ax.axhline(1e-3, color="#cb181d", linestyle="--", label="1e-3 dex gate")
ax.set_xticks(positions, spectrum_labels, rotation=20, ha="right", fontsize=8)
ax.set_ylabel("maximum spectrum difference (dex)")
ax.set_title("Same-root frequency-bin comparison")
ax.grid(axis="y", alpha=0.25)
ax.legend(fontsize=8)

ax = axes[1, 1]
labels = ["low corner", "Sobol holdout"]
positions = np.arange(2)
width = 0.34
p50 = [row["p50_ratio"] for row in benchmark]
p95 = [row["p95_ratio"] for row in benchmark]
ax.bar(positions - width / 2, p50, width, label="P50 ratio", color="#3182bd")
ax.bar(positions + width / 2, p95, width, label="P95 ratio", color="#fd8d3c")
ax.axhline(1.05, color="#cb181d", linestyle="--", label="+5% gate")
ax.set_xticks(positions, labels)
ax.set_ylabel("current helper / previous helper")
ax.set_title("100 alternating complete Fast solves")
ax.grid(axis="y", alpha=0.25)
ax.legend(fontsize=8)

fig.savefig(HERE / "cr1_b175_reheat_kappa_upper_evidence.png", dpi=180)
