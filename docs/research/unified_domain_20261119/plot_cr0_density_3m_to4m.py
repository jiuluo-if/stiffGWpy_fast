"""Plot paired CR0 3M–4M coverage, precision and full-solver timing."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = load("screen_cr0_density_3m_to4m.json")["cases"]
    sobol = load("screen_cr0_density_3m_to4m_sobol.json")["cases"]
    reference = load("validate_cr0_density_3m_to4m_z8.json")["cases"]
    corrected_boundary = load("close_cr0_3m_reference_origin_corrected.json")["cases"]
    old_boundary_spectrum = load("validate_cr0_3m_wkb_guard_z8_multi.json")["cases"]
    timing = load("benchmark_cr0_density_3m_to4m_self_ab_tbb_100pairs.json")["cases"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    fig.suptitle("CR0 density continuation: 3M–4M (finite sampled evidence)", fontsize=15)

    ax = axes[0, 0]
    ax.fill_between([3.0, 4.0], -0.5, 0.0, color="#d9ead3", alpha=0.7,
                    label="screened box")
    all_rows = lattice + sobol
    for row in all_rows:
        p = row["input"]["params"]
        old = row["baseline_fast"]["success"]
        new = row["candidate_fast"]["success"]
        marker = "x" if not old else "s"
        ax.scatter(p["Omega_bh2"] / 1e6, p["n_t"], marker=marker, s=40,
                   color="#cb181d" if not old else "#3182bd", zorder=3)
        if new and not old:
            ax.scatter(p["Omega_bh2"] / 1e6, p["n_t"], marker="o", s=22,
                       color="#238b45", zorder=4)
    for item in reference:
        p = item["params"]
        ax.scatter(p["Omega_bh2"] / 1e6, p["n_t"], marker="*", s=120,
                   color="#ffbf00", edgecolor="black", zorder=5)
    for item in corrected_boundary:
        p = item["params"]
        ax.scatter(p["Omega_bh2"] / 1e6, p["n_t"], marker="D", s=48,
                   facecolor="white", edgecolor="#54278f", linewidth=1.4, zorder=6)
    ax.set(xlabel=r"$\Omega_{bh2}$ ($10^6$ units)", ylabel=r"$n_t$",
           title="Paired old/current Fast; ★ Reference; ◇ corrected 3M closure")
    ax.set_xlim(2.92, 4.08)
    ax.set_ylim(-0.56, 0.06)
    ax.grid(alpha=.25)
    ax.text(.02, .03, "Old Fast: 0/19 success; all failures = nonfinite_background_primitive",
            transform=ax.transAxes, fontsize=8, va="bottom")

    corrected_by_id = {r["case_id"]: r for r in corrected_boundary}
    legacy_by_id = {r["case_id"]: r for r in old_boundary_spectrum}
    ref_by_id = {r["case_id"]: r for r in reference}
    labels = []
    root_errors = []
    spectrum_errors = []
    for case_id in ("upper_density_zero_tilt", "upper_density_negative_tilt",
                    "independent_sobol_upper_density"):
        row = corrected_by_id[case_id]
        labels.append("3M 0" if case_id.endswith("zero_tilt") else
                      "3M −.5" if case_id.endswith("negative_tilt") else "3M Sobol")
        root_errors.append(row["fast_root_relative_to_corrected_reference"])
        spectrum_errors.append(legacy_by_id[case_id]["spectrum_max_abs_dex"])
    for case_id in ("upper_density_low_tilt", "upper_density_high_tilt",
                    "independent_sobol_holdout"):
        row = ref_by_id[case_id]
        labels.append("4M −.5" if case_id.endswith("low_tilt") else
                      "4M 0" if case_id.endswith("high_tilt") else "4M Sobol")
        root_errors.append(row["root_relative_error"])
        spectrum_errors.append(row["spectrum_max_abs_dex"])

    x = np.arange(len(labels))
    ax = axes[0, 1]
    ax.bar(x, root_errors, color="#3182bd")
    ax.axhline(2e-4, color="#cb181d", linestyle="--", label="root gate 2e−4")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set(ylabel="relative root error", title="Correct-origin Reference closure")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.bar(x, spectrum_errors, color="#756bb1")
    ax.axhline(1e-3, color="#cb181d", linestyle="--", label="spectrum gate 1e−3 dex")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set(ylabel="maximum same-root spectrum error (dex)", title="Same Fast root and frequency grid")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    perf_labels = [r["case_id"].replace("new_region_", "").replace("trusted_physical_cr0_control", "trusted CR0")
                   for r in timing]
    p50 = [r["p50_ratio"] for r in timing]
    p95 = [r["p95_ratio"] for r in timing]
    width = .36
    ax.bar(x[:len(timing)]-width/2, p50, width, label="P50", color="#31a354")
    ax.bar(x[:len(timing)]+width/2, p95, width, label="P95", color="#a1d99b")
    ax.axhline(1.05, color="#cb181d", linestyle="--", label="5% regression limit")
    ax.set_xticks(x[:len(timing)], perf_labels, rotation=20, ha="right")
    ax.set(ylabel="candidate / source baseline", title="100 paired complete Fast solves")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    output = HERE / "cr0_density_3m_to4m_domain.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
