"""Plot the CR1 log-wide r domain, paired failures, and raw-linear accuracy."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
lattice = json.loads((HERE / "screen_cr1_r_logwide_lattice.json").read_text(encoding="utf-8"))["cases"]
sobol = json.loads((HERE / "screen_cr1_r_logwide_sobol.json").read_text(encoding="utf-8"))["cases"]
refs = json.loads((HERE / "validate_cr1_r_logwide_z8.json").read_text(encoding="utf-8"))["cases"]
cases = lattice + sobol

fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), constrained_layout=True)
ax = axes[0]
for row in cases:
    params = row["input"]["params"]
    old = row["baseline_fast"]["success"]
    new = row["candidate_fast"]["success"]
    if not old and new:
        ax.scatter(params["T_re"], params["r"], marker="o", facecolors="none",
                   edgecolors="#2ca02c", linewidths=1.4, s=46)
    elif old and new:
        ax.scatter(params["T_re"], params["r"], marker="s", c="#4c78a8", s=42)
    else:
        ax.scatter(params["T_re"], params["r"], marker="x", c="#d62728", s=40)
ax.set_yscale("log")
ax.set(xlabel=r"$T_{re}$", ylabel=r"$r$ (log scale)", title="CR1 log-r slice")
ax.set_xlim(1738, 2012)
ax.grid(alpha=0.25, which="both")
ax.plot([], [], "o", mfc="none", mec="#2ca02c", label="old fail → Fast success")
ax.plot([], [], "s", c="#4c78a8", label="old and new success")
ax.legend(loc="best", fontsize=8)

old_success = sum(row["baseline_fast"]["success"] for row in cases)
old_failure = len(cases)-old_success
new_success = sum(row["candidate_fast"]["success"] for row in cases)
axes[1].bar(["old failure", "old success", "current success"],
            [old_failure, old_success, new_success],
            color=["#d62728", "#4c78a8", "#2ca02c"])
axes[1].set(ylabel="sample count", title="Paired outcomes (19 inputs)")
axes[1].tick_params(axis="x", labelrotation=18)
axes[1].grid(axis="y", alpha=0.25)

ax = axes[2]
names = [
    ("cutoff T=2000" if row["params"]["r"] == 4e-43 and row["params"]["T_re"] == 2000
     else row.get("case_id", f"r={row['params']['r']:.1e}")).replace("_", "\n")
    for row in refs
]
x = np.arange(len(refs))
root = np.array([row["root_relative_error_proxy"] for row in refs])
spectrum = np.array([row["spectrum_max_abs_dex_from_raw_linear"] for row in refs])
width = 0.34
bars_root = ax.bar(x-width/2, root, width, label="root relative error", color="#4c78a8")
bars_spectrum = ax.bar(x+width/2, spectrum, width, label="spectrum max (dex)", color="#f2a541")
ax.axhline(2e-4, color="#4c78a8", linestyle="--", linewidth=1, label="root gate")
ax.axhline(1e-3, color="#f2a541", linestyle=":", linewidth=1, label="spectrum gate")
ax.set_xticks(x, names, fontsize=6)
ax.set(ylabel="error", title="z8 Reference raw-linear checks")
ax.legend(fontsize=6)
ax.grid(axis="y", alpha=0.25)
for bars in (bars_root, bars_spectrum):
    ax.bar_label(bars, fmt="%.2g", fontsize=5, padding=2)

fig.suptitle("Finite-sample CR1 cutoff-to-r=.005 extension; no continuum claim")
output = HERE / "cr1_r_logwide_evidence.png"
fig.savefig(output, dpi=180)
print(output)
