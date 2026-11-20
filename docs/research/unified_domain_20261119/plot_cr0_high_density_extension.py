from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = ROOT / "docs/research/unified_domain_20261119"
grid = json.loads((HERE / "screen_cr0_high_density_extension.json").read_text(encoding="utf-8"))
sobol = json.loads((HERE / "screen_cr0_high_density_extension_sobol.json").read_text(encoding="utf-8"))
reference = json.loads((HERE / "validate_cr0_high_density_extension_z8.json").read_text(encoding="utf-8"))

rows = grid["cases"] + sobol["cases"]
fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)

ax = axes[0, 0]
families = {"rectangle": ("#2a6fbb", "s"), "holdout": ("#e07a23", "D")}
for row in rows:
    params = row["input"]["params"]
    family = row["input"]["family"]
    color, marker = families[family]
    density = params["Omega_bh2"] / 1.0e6
    tilt = params["n_t"]
    ax.scatter(density, tilt, marker="x", color="#555555", s=42, linewidth=1.5,
               label="Baseline failure" if row is rows[0] else None, zorder=2)
    ax.scatter(density + 0.025, tilt, marker=marker, color=color, s=38,
               label=f"Current Fast {family}" if not any(h.get_label() == f"Current Fast {family}" for h in ax.collections) else None,
               zorder=3)
for row in reference["cases"]:
    params = row["params"]
    ax.scatter(params["Omega_bh2"] / 1.0e6, params["n_t"], marker="*", color="#b51f37",
               edgecolor="white", linewidth=0.5, s=145,
               label="Direct z8 Reference" if row is reference["cases"][0] else None, zorder=4)
ax.set_xlim(0.94, 2.08)
ax.set_ylim(-0.54, 0.04)
ax.set_xlabel(r"$\Omega_{bh} h^2$ (millions; current markers offset +0.025 for visibility)")
ax.set_ylabel(r"Independent $n_t$ (cr=0)")
ax.set_title("Paired coverage: 3×3 grid + interior + Sobol")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(alpha=0.2)

ax = axes[0, 1]
labels = ["Baseline", "Current Fast"]
success = [
    sum(bool(row["baseline_fast"]["success"]) for row in rows),
    sum(bool(row["candidate_fast"]["success"]) for row in rows),
]
nonfinite = [
    sum(row["baseline_fast"]["failure_reason"] == "nonfinite_background_primitive" for row in rows),
    sum(row["candidate_fast"]["failure_reason"] == "nonfinite_background_primitive" for row in rows),
]
xx = np.arange(2)
ax.bar(xx, success, color=["#888888", "#3d8b63"], label="Success")
ax.bar(xx, nonfinite, bottom=success, color=["#b85252", "#dddddd"], label="nonfinite background primitive")
ax.set_xticks(xx, labels)
ax.set_ylabel("Samples (n=19)")
ax.set_title("Outcome and failure class")
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=0.2)

ref_cases = reference["cases"]
ref_labels = ["2M, nt=0", "2M, nt=-.5", "1M, nt=0", "Sobol holdout"]
xx = np.arange(len(ref_cases))

ax = axes[1, 0]
root_errors = [row["root_relative_error"] for row in ref_cases]
ax.bar(xx, root_errors, color="#4778a8")
ax.axhline(2e-4, color="#b33b3b", linestyle="--", linewidth=1, label="Root limit 2e-4")
ax.set_xticks(xx, ref_labels, rotation=15)
ax.set_ylabel("Closed-root relative error")
ax.set_title("Direct z8 Reference: root")
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=0.2)

ax = axes[1, 1]
spectrum_errors = [row["spectrum_max_abs_dex"] for row in ref_cases]
ax.bar(xx, spectrum_errors, color="#7b66a7")
ax.axhline(1e-3, color="#b33b3b", linestyle="--", linewidth=1, label="Spectrum limit 1e-3 dex")
ax.set_xticks(xx, ref_labels, rotation=15)
ax.set_ylabel("Maximum spectrum error (dex)")
ax.set_title("Direct z8 Reference: spectrum")
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=0.2)

fig.suptitle("CR0 high-density extension: finite-sample evidence", fontsize=14)
output = HERE / "cr0_high_density_extension_domain.png"
fig.savefig(output, dpi=180)
print(output)
