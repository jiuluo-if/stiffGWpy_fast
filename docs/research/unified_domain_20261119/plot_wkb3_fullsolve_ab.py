from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[3]
data = json.loads(
    (root / "docs/research/unified_domain_20261119/benchmark_wkb3_fullsolve_screen.json")
    .read_text(encoding="utf-8")
)
cases = data["cases"]
labels = [row["case"].replace("_", "\n") for row in cases]
p50 = np.array([row["p50_ratio"] for row in cases])
p95 = np.array([row["p95_ratio"] for row in cases])
spectral = np.array([row["max_spectrum_absdiff_dex"] for row in cases])

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
x = np.arange(len(cases))
width = 0.34
axes[0].bar(x - width / 2, p50, width, label="P50 candidate / baseline", color="#3973ac")
axes[0].bar(x + width / 2, p95, width, label="P95 candidate / baseline", color="#e28e2c")
axes[0].axhline(1.05, color="#a33", linestyle="--", linewidth=1, label="5% regression limit")
axes[0].axhline(1.0, color="#555", linewidth=0.8)
axes[0].set_xticks(x, labels)
axes[0].set_ylabel("Runtime ratio (below 1 is faster)")
axes[0].set_title("64 alternating full-solve pairs")
axes[0].legend(frameon=False, fontsize=8)
axes[0].grid(axis="y", alpha=0.2)

positive = np.maximum(spectral, 1e-13)
axes[1].bar(x, positive, color=["#5b9a69"] * 2 + ["#8b6bb1"] * 3)
axes[1].axhline(1e-3, color="#a33", linestyle="--", linewidth=1, label="Spectrum limit: 1e-3 dex")
axes[1].set_yscale("log")
axes[1].set_xticks(x, labels)
axes[1].set_ylabel("Maximum candidate-baseline spectrum difference (dex)")
axes[1].set_title("Self-consistency root identical in every pair")
axes[1].legend(frameon=False, fontsize=8)
axes[1].grid(axis="y", which="both", alpha=0.2)

fig.suptitle("Local adaptive WKB3: trusted cases and high-tilt CR0", fontsize=13)
output = root / "docs/research/unified_domain_20261119/wkb3_fullsolve_ab.png"
fig.savefig(output, dpi=180)
print(output)
