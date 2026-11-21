"""Plot paired CR1 coverage and representative precision checks."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
lattice = json.loads((HERE / "screen_cr1_nearby_reheat_kappa_lattice.json").read_text(encoding="utf-8"))["cases"]
sobol = json.loads((HERE / "screen_cr1_nearby_reheat_kappa_sobol.json").read_text(encoding="utf-8"))["cases"]
refs = json.loads((HERE / "validate_cr1_nearby_reheat_kappa_z8.json").read_text(encoding="utf-8"))["cases"]
cases = lattice + sobol
params = [row["input"]["params"] for row in cases]
old_ok = [row["baseline_fast"]["success"] for row in cases]
new_ok = [row["candidate_fast"]["success"] for row in cases]
reason_counts: dict[str, int] = {}
for row in cases:
    reason = row["baseline_fast"].get("failure_reason") or "success"
    reason_counts[reason] = reason_counts.get(reason, 0) + 1

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
ax = axes[0]
for row, p in zip(cases, params):
    old = row["baseline_fast"]["success"]
    new = row["candidate_fast"]["success"]
    if old and new:
        ax.scatter(p["T_re"], p["kappa10"], marker="s", c="#4c78a8", s=36)
    elif not old and new:
        ax.scatter(p["T_re"], p["kappa10"], marker="o", facecolors="none",
                   edgecolors="#2ca02c", linewidths=1.5, s=48)
    elif old and not new:
        ax.scatter(p["T_re"], p["kappa10"], marker="x", c="#d62728", s=40)
    else:
        ax.scatter(p["T_re"], p["kappa10"], marker="x", c="#9467bd", s=40)
ax.set(xlabel=r"$T_{re}$", ylabel=r"$\kappa_{10}$", title="CR1 parameter slice")
ax.set_xlim(1678, 1759)
ax.set_ylim(0.0121, 0.0154)
ax.grid(alpha=0.25)
ax.plot([], [], "o", mfc="none", mec="#2ca02c", label="old fail → Fast success")
ax.legend(loc="best", fontsize=8)

ax = axes[1]
labels = list(reason_counts)
counts = [reason_counts[label] for label in labels]
labels.append("current Fast success")
counts.append(sum(new_ok))
colors = ["#d62728" if label != "success" else "#4c78a8" for label in labels[:-1]] + ["#2ca02c"]
ax.barh(labels, counts, color=colors)
ax.set(xlabel="sample count", title="Paired outcomes (19 inputs)")
ax.grid(axis="x", alpha=0.25)

ax = axes[2]
names = [row["case_id"].replace("_", "\n") for row in refs]
x = np.arange(len(refs))
root = np.array([row["root_relative_error_proxy"] for row in refs])
spectrum = np.array([row["spectrum_max_abs_dex"] for row in refs])
width = 0.34
bars1 = ax.bar(x - width/2, root, width, label="root relative error", color="#4c78a8")
bars2 = ax.bar(x + width/2, spectrum, width, label="spectrum max (dex)", color="#f2a541")
ax.axhline(2e-4, color="#4c78a8", linestyle="--", linewidth=1, label="root gate")
ax.axhline(1e-3, color="#f2a541", linestyle=":", linewidth=1, label="spectrum gate")
ax.set_xticks(x, names, fontsize=7)
ax.set(ylabel="error", title="z8 Reference checks")
ax.legend(fontsize=7)
ax.grid(axis="y", alpha=0.25)
for bars in (bars1, bars2):
    ax.bar_label(bars, fmt="%.2g", fontsize=6, padding=2)

fig.suptitle("Finite-sample CR1 extension evidence; no continuum claim")
output = HERE / "cr1_nearby_reheat_kappa_evidence.png"
fig.savefig(output, dpi=180)
print(output)
