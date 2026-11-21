"""Plot paired classes and Reference accuracy for the CR0 temperature slice."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


accepted = (
    load("screen_cr0_density_reheat_negative_tilt_lattice.json")["cases"]
    + load("screen_cr0_density_reheat_negative_tilt_sobol.json")["cases"]
)
broad = (
    load("screen_cr0_density_reheat_tilt_lattice.json")["cases"]
    + load("screen_cr0_density_reheat_tilt_sobol.json")["cases"]
)
positive_tilt = (
    load("screen_cr0_tilt_temperature_lattice.json")["cases"]
    + load("screen_cr0_tilt_temperature_sobol.json")["cases"]
)
reference = load("validate_cr0_density_reheat_negative_tilt_z8.json")["cases"]
closure = load("close_reference_cr0_density_reheat_negative_tilt_z8.json")["case"]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

ax = axes[0]
for case in accepted:
    params = case["input"]["params"]
    old = case["baseline_fast"]
    if old.get("success"):
        ax.scatter(params["T_re"], params["n_t"], marker="s", c="#4c78a8", s=36, zorder=3)
    else:
        ax.scatter(params["T_re"], params["n_t"], marker="o", facecolors="none",
                   edgecolors="#2ca02c", linewidths=1.1, s=36, zorder=2)
for case in broad:
    params = case["input"]["params"]
    current = case["candidate_fast"]
    if current.get("failure_reason") == "shared_Neff_guard":
        ax.scatter(params["T_re"], params["n_t"], marker="^", c="#d62728", s=44, zorder=4)
for case in positive_tilt:
    params = case["input"]["params"]
    if case["candidate_fast"].get("failure_reason") == "shared_Neff_guard":
        ax.scatter(params["T_re"], params["n_t"], marker="^", c="#d62728", s=34, zorder=4)
ax.add_patch(plt.Rectangle((4e3, -0.5), 1e7 - 4e3, 0.25, fill=False,
                           edgecolor="#2ca02c", linewidth=1.4))
ax.axvline(2000.0, color="#4c78a8", linestyle="--", linewidth=1,
           label="prior CR0 high-density T sample")
ax.set_xscale("log")
ax.set_xlim(1e3, 1.5e7)
ax.set(xlabel=r"$T_{re}$ (GeV)", ylabel=r"$n_t$",
       title="CR0 paired coverage and guard boundary")
ax.set_ylim(-0.54, 0.55)
ax.grid(alpha=.25, which="both")
ax.plot([], [], "o", mfc="none", mec="#2ca02c", label="old primitive failure → current success")
ax.plot([], [], "s", c="#4c78a8", label="old success retained")
ax.plot([], [], "^", c="#d62728", label="shared N_eff guard rejection")
ax.legend(fontsize=7, loc="best")

ax = axes[1]
labels = ["old Fast", "current Fast", "wider box", "positive tilt"]
heights = [17, 19, 3, 19]
colors = ["#d62728", "#2ca02c", "#9467bd", "#9467bd"]
bars = ax.bar(labels, heights, color=colors)
for bar, value, text in zip(bars, heights, ("17 primitive failures\n2 prior successes", "19/19 successes", "3 guard rejects", "19 guard rejects")):
    ax.text(bar.get_x() + bar.get_width()/2, value + .35, text, ha="center", va="bottom", fontsize=8)
ax.set(ylabel="paired sample count", title="Failure categories remain separate")
ax.set_ylim(0, 23)
ax.tick_params(axis="x", labelrotation=12)
ax.grid(axis="y", alpha=.25)

ax = axes[2]
for case in reference:
    root = case["root_relative_error_proxy"]
    if case["case_id"] == closure["case_id"]:
        root = closure["fast_root_relative_to_third_map"]
    ax.scatter(root, case["spectrum_max_abs_dex_from_raw_linear"],
               c="#2ca02c", s=36, alpha=.85)
ax.axvline(2e-4, color="#4c78a8", linestyle="--", linewidth=1, label="root gate")
ax.axhline(1e-3, color="#f2a541", linestyle=":", linewidth=1, label="spectrum gate")
ax.set(xlabel="root relative error", ylabel="spectrum max error (dex)",
       title="Current-source z8 Reference checks")
ax.set_xlim(0, .00021)
ax.set_xticks([0.0, 1e-4, 2e-4], labels=["0", "1e-4", "2e-4"])
ax.set_ylim(0, .00105)
ax.grid(alpha=.25)
ax.legend(fontsize=7, loc="best")

fig.suptitle("CR0 density × reheating × negative-tilt slice; finite samples, not continuum certification")
output = HERE / "cr0_density_reheat_negative_tilt_evidence.png"
fig.savefig(output, dpi=180)
print(output)
