"""Plot paired classifications and Reference accuracy for the CR1 T extension."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def rows(filename: str) -> list[dict]:
    return json.loads((HERE / filename).read_text(encoding="utf-8"))["cases"]


extension = rows("screen_cr1_r_temperature_extension_lattice.json") + rows(
    "screen_cr1_r_temperature_extension_sobol.json"
)
nextband = rows("screen_cr1_r_nextband_lattice.json") + rows("screen_cr1_r_nextband_sobol.json")
narrow = rows("screen_cr1_r_narrowband_lattice.json") + rows("screen_cr1_r_narrowband_sobol.json")
extension_refs = json.loads(
    (HERE / "validate_cr1_r_temperature_extension_z8.json").read_text(encoding="utf-8")
)["cases"]
closure = json.loads(
    (HERE / "close_reference_cr1_r_temperature_extension_sobol_z8.json").read_text(encoding="utf-8")
)["case"]
next_refs = json.loads((HERE / "validate_cr1_r_nextband_z8.json").read_text(encoding="utf-8"))["cases"]
narrow_refs = json.loads((HERE / "validate_cr1_r_narrowband_z8.json").read_text(encoding="utf-8"))["cases"]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

ax = axes[0]
for case in extension:
    params = case["input"]["params"]
    baseline = case["baseline_fast"]
    if baseline.get("success"):
        ax.scatter(params["T_re"], params["r"], marker="s", c="#4c78a8", s=34, zorder=3)
    else:
        ax.scatter(params["T_re"], params["r"], marker="o", facecolors="none",
                   edgecolors="#2ca02c", linewidths=1.1, s=34, zorder=2)
for group in (nextband, narrow):
    for case in group:
        params = case["input"]["params"]
        ax.scatter(params["T_re"], params["r"], marker="x", c="#d62728",
                   linewidths=1.0, s=28, zorder=1)
ax.add_patch(plt.Rectangle((1875, .005), 125, .005, fill=False,
                           edgecolor="#4c78a8", linewidth=1.2, linestyle="--"))
ax.add_patch(plt.Rectangle((2000, .005), 250, .005, fill=False,
                           edgecolor="#2ca02c", linewidth=1.5))
ax.set(xlabel=r"$T_{re}$ (GeV)", ylabel=r"$r$",
       title="CR1 coverage: prior and temperature extension")
ax.set_xlim(1850, 2275)
ax.set_ylim(.0045, .0158)
ax.grid(alpha=.25)
ax.plot([], [], "o", mfc="none", mec="#2ca02c", label="old primitive failure → current success")
ax.plot([], [], "s", c="#4c78a8", label="old Fast success retained")
ax.plot([], [], "x", c="#d62728", label="screened higher-r point")
ax.legend(fontsize=7, loc="best")

ax = axes[1]
labels = ["T extension", "higher-r rejected", "narrow-r rejected"]
counts = [
    (sum(bool(c["baseline_fast"].get("success")) for c in extension), len(extension)),
    (sum(bool(c["baseline_fast"].get("success")) for c in nextband), len(nextband)),
    (sum(bool(c["baseline_fast"].get("success")) for c in narrow), len(narrow)),
]
x = range(len(labels))
ax.bar([i - .18 for i in x], [n - old for old, n in counts], width=.36,
       color="#2ca02c", label="old failure → current success")
ax.bar([i + .18 for i in x], [old for old, _ in counts], width=.36,
       color="#4c78a8", label="old success")
for i, (_, total) in enumerate(counts):
    ax.text(i, total + .35, f"current {total}/{total}", ha="center", fontsize=8)
ax.set_xticks(list(x), labels, rotation=10)
ax.set(ylabel="paired sample count", title="Failure classes kept separate")
ax.set_ylim(0, 25)
ax.grid(axis="y", alpha=.25)
ax.legend(fontsize=7, loc="upper right")

ax = axes[2]
for group, color, label in (
    (extension_refs, "#2ca02c", "temperature extension"),
    (next_refs, "#d62728", "higher-r probes"),
    (narrow_refs, "#f2a541", "narrow-strip probes"),
):
    for case in group:
        root = case["root_relative_error_proxy"]
        if case["case_id"] == closure["case_id"]:
            root = closure["fast_root_relative_to_third_map"]
        spec = case.get("spectrum_max_abs_dex_from_raw_linear")
        if spec is not None:
            ax.scatter(root, spec, color=color, s=34, alpha=.85, label=label)
            label = "_nolegend_"
ax.axvline(2e-4, color="#4c78a8", linestyle="--", linewidth=1, label="root gate")
ax.axhline(1e-3, color="#f2a541", linestyle=":", linewidth=1, label="spectrum gate")
ax.set(xlabel="root relative error", ylabel="spectrum max error (dex)",
       title="Reference checks; rejected probes shown")
ax.set_xlim(0, .00031)
ax.set_ylim(0, .0011)
ax.grid(alpha=.25)
handles, legend_labels = ax.get_legend_handles_labels()
dedup = dict(zip(legend_labels, handles))
ax.legend(dedup.values(), dedup.keys(), fontsize=7, loc="best")

fig.suptitle("CR1 finite-sample expansion evidence (not a continuum certification)")
output = HERE / "cr1_r_temperature_extension_evidence.png"
fig.savefig(output, dpi=180)
print(output)
