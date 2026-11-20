"""Plot lower-density CR0 tilt support, Reference accuracy and runtime."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def main() -> None:
    lattice = load("screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice.json")["cases"]
    sobol = load("screen_cr0_density_13p875m_to13p9375m_nt0302734375_slice_sobol.json")["cases"]
    reference = load("probe_cr0_density_14m_nt0302734375_reference.json")["cases"]
    holdout = load("probe_cr0_density_13p875m_to13p9375m_nt0302734375_sobol_reference.json")["case"]
    closure = load("probe_cr0_density_13p9375m_nt0302734375_closure_step.json")
    timing = load("benchmark_cr0_density_13p875m_to13p9375m_nt0302734375_self_ab_tbb_500pairs.json")["cases"]
    probes = [r for r in reference if r["reference_status"] == "success"]
    rows = lattice + sobol
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)
    fig.suptitle("CR0 density-dependent tilt edge: 13.875M–13.9375M, −0.0302734375", fontsize=14)

    ax = axes[0, 0]
    ax.fill_between([13.875, 13.9375], -0.0302734375, 0.0, color="#d9ead3", alpha=.75,
                    label="screened box")
    for row in rows:
        params = row["input"]["params"]
        old_ok = row["baseline_fast"]["success"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="s" if old_ok else "x",
                   s=43, color="#cb181d" if not old_ok else "#3182bd", zorder=3)
        if row["candidate_fast"]["success"] and not old_ok:
            ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="o", s=22,
                       color="#238b45", zorder=4)
    for row in probes:
        params = row["params"]
        ax.scatter(params["Omega_bh2"] / 1e6, params["n_t"], marker="*", s=115,
                   color="#ffbf00", edgecolor="black", zorder=5)
    hp = holdout["params"]
    ax.scatter(hp["Omega_bh2"] / 1e6, hp["n_t"], marker="*", s=120,
               color="#6a51a3", edgecolor="black", zorder=6, label="Sobol Reference")
    ax.scatter(14.0, -0.0302734375, marker="v", s=70, color="#d7301f", zorder=6,
               label="14M Reference failure; excluded")
    ax.set(xlabel=r"$\Omega_{bh2}$ ($10^6$ units)", ylabel=r"$n_t$",
           title="Old Fast failures ×; new Fast success ○; Reference ★")
    ax.set_xlim(13.86, 14.015)
    ax.set_ylim(-0.033, 0.002)
    ax.grid(alpha=.25)
    ax.legend(fontsize=7, loc="lower right")
    ax.text(.02, .96, "Old Fast 0/19; current Fast 19/19\nOld failure: nonfinite_background_primitive",
            transform=ax.transAxes, fontsize=8, va="top")

    samples = probes + [{
        "case_id": "independent Sobol",
        "root_relative_error": holdout["root_relative_error"],
        "spectrum_max_abs_dex": holdout["spectrum_max_abs_dex"],
    }]
    labels = ["lower edge", "interior", "Sobol"]
    x = np.arange(len(samples))
    reference_interior = next(r for r in reference if r["case_id"] == "interior")
    closure_residual = abs(closure["closure_step"]["mapped_dn_gw"]-reference_interior["reference_mapped_root"]) / closure["closure_step"]["mapped_dn_gw"]
    ax = axes[0, 1]
    ax.bar(x, [r["root_relative_error"] for r in samples], color="#3182bd")
    ax.axhline(2e-4, color="#cb181d", linestyle="--", label="root gate 2e−4")
    ax.set_xticks(x, labels, rotation=12, ha="right")
    ax.set(ylabel="relative root error / fixed-point residual", title="z8 Reference root checks")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)
    ax.text(.02, .84, f"13.9375M second-map residual {closure_residual:.2e}",
            transform=ax.transAxes, fontsize=8, va="top")

    ax = axes[1, 0]
    ax.bar(x, [r["spectrum_max_abs_dex"] for r in samples], color="#756bb1")
    ax.axhline(1e-3, color="#cb181d", linestyle="--", label="spectrum gate 1e−3 dex")
    ax.set_xticks(x, labels, rotation=12, ha="right")
    ax.set(ylabel="maximum spectrum error (dex)", title="75-frequency Reference comparison")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    labels = ["lower edge", "interior", "Sobol", "trusted"]
    x = np.arange(len(labels))
    width = .36
    ax.bar(x-width/2, [r["p50_ratio"] for r in timing], width, label="P50", color="#31a354")
    ax.bar(x+width/2, [r["p95_ratio"] for r in timing], width, label="P95", color="#a1d99b")
    ax.axhline(1.05, color="#cb181d", linestyle="--", label="5% regression limit")
    ax.set_xticks(x, labels)
    ax.set(ylabel="candidate / Fast baseline", title="500 paired full Fast solves")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=8)

    output = HERE / "cr0_density_13p875m_to13p9375m_nt0302734375_slice.png"
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
