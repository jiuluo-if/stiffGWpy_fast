"""Plot paired coverage, outcomes, Reference error, and WKB guard A/B."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = ROOT / "docs" / "research" / "unified_domain_20261119"
grid = json.loads((HERE / "screen_cr0_density_2m_to3m_wkb_guard.json").read_text())
sobol = json.loads(
    (HERE / "screen_cr0_density_2m_to3m_wkb_guard_sobol.json").read_text()
)
reference = json.loads(
    (HERE / "validate_cr0_3m_wkb_guard_z8_multi.json").read_text()
)
benchmark = json.loads((HERE / "benchmark_wkb_basis_guard_ab.json").read_text())

rows = grid["cases"] + sobol["cases"]
fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)

ax = axes[0, 0]
families = {"rectangle": ("#2a6fbb", "s"), "holdout": ("#e07a23", "D")}
for row in rows:
    params = row["input"]["params"]
    family = row["input"]["family"]
    color, marker = families[family]
    x = params["Omega_bh2"] / 1.0e6
    y = params["n_t"]
    ax.scatter(
        x,
        y,
        marker="x",
        color="#555555",
        s=42,
        linewidth=1.4,
        label="Pinned baseline failed" if row is rows[0] else None,
        zorder=2,
    )
    ax.scatter(
        x + 0.018,
        y,
        marker=marker,
        color=color,
        s=38,
        label=f"Current Fast {family}"
        if not any(c.get_label() == f"Current Fast {family}" for c in ax.collections)
        else None,
        zorder=3,
    )
for index, row in enumerate(reference["cases"]):
    params = row["params"]
    ax.scatter(
        params["Omega_bh2"] / 1.0e6,
        params["n_t"],
        marker="*",
        color="#b51f37",
        edgecolor="white",
        linewidth=0.5,
        s=145,
        label="Direct z8 Reference" if index == 0 else None,
        zorder=4,
    )
ax.set_xlim(1.92, 3.08)
ax.set_ylim(-0.54, 0.04)
ax.set_xlabel(r"$\Omega_{bh} h^2$ (millions; Fast markers offset +0.018)")
ax.set_ylabel(r"Independent $n_t$ ($cr=0$)")
ax.set_title("Paired coverage: 3×3 grid, 2 interior, 8 Sobol")
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(alpha=0.2)

ax = axes[0, 1]
success = [
    sum(bool(row["baseline_fast"]["success"]) for row in rows),
    sum(bool(row["candidate_fast"]["success"]) for row in rows),
]
old_nonfinite = sum(
    row["baseline_fast"]["failure_reason"] == "nonfinite_background_primitive"
    for row in rows
)
other_old = len(rows) - success[0] - old_nonfinite
xx = np.arange(2)
ax.bar(xx, success, color=["#888888", "#3d8b63"], label="Fast success")
ax.bar(
    xx,
    [old_nonfinite, 0],
    bottom=success,
    color="#b85252",
    label="nonfinite_background_primitive",
)
if other_old:
    ax.bar(xx[0], other_old, bottom=success[0] + old_nonfinite, color="#d79539", label="Other failure")
ax.set_xticks(xx, ["Pinned old Fast\n0/19 success", "WKB guard Fast\n19/19 success"])
ax.set_ylabel("Samples (n=19)")
ax.set_title("Outcome and failure class")
ax.text(0, 9.4, "19 × nonfinite_background_primitive", color="white", rotation=90,
        ha="center", va="center", fontsize=8)
ax.text(1, 9.4, "19 successful", color="white", rotation=90,
        ha="center", va="center", fontsize=8)
ax.grid(axis="y", alpha=0.2)

ax = axes[1, 0]
ref_cases = reference["cases"]
ref_labels = [row["case_id"].replace("_", "\n") for row in ref_cases]
xx = np.arange(len(ref_cases))
root_errors = [row["root_relative_error"] / 2e-4 for row in ref_cases]
spectrum_errors = [row["spectrum_max_abs_dex"] / 1e-3 for row in ref_cases]
ax.bar(xx - 0.18, root_errors, width=0.36, color="#4778a8", label="Root error / 2e-4")
ax.bar(xx + 0.18, spectrum_errors, width=0.36, color="#7b66a7", label="Spectrum error / 1e-3 dex")
ax.axhline(1.0, color="#b33b3b", linestyle="--", linewidth=1, label="Acceptance limit")
ax.set_xticks(xx, ref_labels, fontsize=8)
ax.set_ylabel("Error / acceptance limit")
ax.set_title("Same-Fast-root Reference checks")
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax.grid(axis="y", alpha=0.2)

ax = axes[1, 1]
bench_cases = benchmark["cases"]
bench_labels = [row["case_id"].replace("_", "\n") for row in bench_cases]
xx = np.arange(len(bench_cases))
p50 = [row["p50_ratio"] for row in bench_cases]
p95 = [row["p95_ratio"] for row in bench_cases]
ax.bar(xx - 0.18, p50, width=0.36, color="#3d8b63", label="P50 candidate/base")
ax.bar(xx + 0.18, p95, width=0.36, color="#70ad91", label="P95 candidate/base")
ax.axhline(1.05, color="#b33b3b", linestyle="--", linewidth=1, label="5% regression limit")
ax.axhline(1.0, color="#555555", linestyle=":", linewidth=0.8)
ax.set_xticks(xx, bench_labels, fontsize=8)
ax.set_ylabel("Runtime ratio")
ax.set_title("Full-solve WKB guard A/B (32 paired runs)")
ax.legend(frameon=False, fontsize=7)
ax.grid(axis="y", alpha=0.2)

fig.suptitle("CR0 density extension after float64 WKB inverse guard", fontsize=14)
output = HERE / "cr0_density_2m_to3m_wkb_guard_domain.png"
fig.savefig(output, dpi=180)
print(output)
