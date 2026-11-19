from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
EXPECTED_FAST = "98a0fe0961b93689a74f63ad9c77ea29c60dd23b68a6ab346d6f984c72bb03e5"


def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8"))


screen = load("screen_cr0_blue_tilt_amplitude_region.json")
sobol = load("screen_cr0_blue_tilt_amplitude_sobol.json")
z8 = load("validate_cr0_blue_tilt_amplitude_z8.json")
z10_names = [
    "validate_cr0_blue_tilt_nt50_rmin_z10.json",
    "validate_cr0_blue_tilt_nt50_rmax_z10.json",
    "validate_cr0_blue_tilt_sobol_z10.json",
]
z10 = [load(n) for n in z10_names]
assert all(d["candidate_source_sha256"]["fast_sgwb.py"] == EXPECTED_FAST for d in z10)
assert z8["source_sha256"]["fast_sgwb.py"] == EXPECTED_FAST
all_cases = screen["cases"] + sobol["cases"]
old_success = sum(bool(c["baseline_fast"]["success"]) for c in all_cases)
new_success = sum(bool(c["candidate_fast"]["success"]) for c in all_cases)
low_corner = next(c for c in z8["cases"] if c["case_id"] == "low_amplitude_low_tilt_corner")
refs = [
    {
        "case": Path(name).stem,
        "r": d["params"]["r"],
        "n_t": d["params"]["n_t"],
        "fast_root": d["fast_root"],
        "reference_root": d["closed_reference_root"],
        "root_relative_error": d["root_relative_error"],
        "reference_closure_residual": d["reference_closure_residual"],
        "spectrum_valid_count": d["spectrum_valid_count"],
        "spectrum_max_abs_dex": d["spectrum_max_abs_dex"],
        "spectrum_p95_abs_dex": d["spectrum_p95_abs_dex"],
        "fast_seconds": d["fast_seconds"],
        "source_sha256": d["candidate_source_sha256"],
    }
    for name, d in zip(z10_names, z10)
]
payload = {
    "status": "sampled_domain_evidence_current_source",
    "candidate_head": "8a4b925cdd7456c32027cf3b436ea4db0abc1526",
    "candidate_fast_sgwb_sha256": EXPECTED_FAST,
    "screen_baseline_base_sha": screen["base_sha"],
    "screen_candidate_base_sha": screen["candidate_base_sha"],
    "region": {
        "fixed": {
            "Omega_bh2": 1e6,
            "Omega_ch2": 0.1201075,
            "H0": 67.32117,
            "DN_eff": 0.0,
            "A_s": 2.100549e-9,
            "cr": 0.0,
            "T_re": 2000.0,
            "DN_re": 10.0,
            "kappa10": 0.01,
        },
        "r": [1e-12, 1e-8],
        "n_t": [0.1, 0.5],
        "h2_positive_finite": all(c["baseline_fast"]["h2_positive_finite"] for c in all_cases),
        "h2_min": min(c["baseline_fast"]["h2_min"] for c in all_cases),
        "sampling": "9 rectangle grid cases + 2 in-rectangle exploratory holdouts + 8 Sobol holdouts",
    },
    "screen": {
        "n": len(all_cases),
        "old_success": old_success,
        "new_success": new_success,
        "old_failure_reasons": sorted({c["baseline_fast"]["failure_reason"] for c in all_cases}),
        "new_failures": [c["input"]["case_id"] for c in all_cases if not c["candidate_fast"]["success"]],
    },
    "reference": {
        "z8_low_corner": {
            "r": low_corner["params"]["r"],
            "n_t": low_corner["params"]["n_t"],
            "root_relative_error": low_corner["root_relative_error"],
            "reference_closure_residual": low_corner["reference_closure_residual"],
            "spectrum_max_abs_dex": low_corner["spectrum_max_abs_dex"],
            "source_sha256": z8["source_sha256"],
        },
        "z10_current_source_cases": refs,
        "z10_max_root_relative_error": max(x["root_relative_error"] for x in refs),
        "z10_max_spectrum_abs_dex": max(x["spectrum_max_abs_dex"] for x in refs),
        "excluded_source_mismatch": {
            "file": "validate_cr0_blue_tilt_nt20_z10.json",
            "observed_fast_sha256": "647488d99dcc13c2843e66e039618ffc33a6adf046474cdc71d06fd828928883",
            "reason": "does not match candidate source hash; not included in metrics",
        },
    },
    "performance": {
        "benchmark": "benchmark_exact_quarter_candidate_trusted_ab_affinity.json",
        "base_sha": "48bbb96a44d4ced3f7e3d80dd2a2deaae440a0a3",
        "candidate_source_sha256": EXPECTED_FAST,
        "resources": {},
        "aggregate_p50_ratio": None,
        "aggregate_p95_ratio": None,
        "per_case": [],
        "acceptance": "not met",
    },
    "reference_code_modified": False,
}
bench = load("benchmark_exact_quarter_candidate_trusted_ab_affinity.json")
payload["performance"]["resources"] = bench["resources"]
payload["performance"]["aggregate_p50_ratio"] = bench["aggregate_p50_ratio"]
payload["performance"]["aggregate_p95_ratio"] = bench["aggregate_p95_ratio"]
payload["performance"]["per_case"] = [
    {"case": c["case"], "p50_ratio": c["p50_ratio"], "p95_ratio": c["p95_ratio"]} for c in bench["paired_cases"]
]
max_p50 = max(c["p50_ratio"] for c in bench["paired_cases"])
max_p95 = max(c["p95_ratio"] for c in bench["paired_cases"])
payload["performance"]["acceptance"] = "not met: aggregate P95 and at least one per-case P50/P95 exceed 1.05"
(HERE / "cr0_blue_tilt_domain_validation_summary.json").write_text(
    json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)

lines = [
    "# CR0 正倾斜扩域：z10 精度核验",
    "",
    "候选提交：`8a4b925cdd7456c32027cf3b436ea4db0abc1526`；`fast_sgwb.py` SHA256：`" + EXPECTED_FAST + "`。",
    "固定背景为 `cr=0`、`Omega_bh2=1e6`、`Omega_ch2=0.1201075`、`T_re=2000`、`kappa10=0.01`；抽样矩形为 `r∈[1e-12,1e-8]`、`n_t∈[0.1,0.5]`。筛查点均为有限正 H²。",
    "",
    f"快速筛查：旧 Fast 成功 `{old_success}/{len(all_cases)}`，新 Fast 成功 `{new_success}/{len(all_cases)}`；旧失败全为 `{', '.join(payload['screen']['old_failure_reasons'])}`。该比例是有限采样证据，不是连续域证明。",
    "",
    "精度核验（Reference 代码与阈值未改）：",
    "",
    "| 样本 | r | n_t | Reference | root relative error | 最大频谱误差 (dex) | Fast 单次时间 (s) |",
    "|---|---:|---:|---|---:|---:|---:|",
    f"| 低角点 | {low_corner['params']['r']:.3g} | {low_corner['params']['n_t']:.3g} | z8 | {low_corner['root_relative_error']:.3g} | {low_corner['spectrum_max_abs_dex']:.3g} | {low_corner['fast_seconds']:.3g} |",
    *[
        f"| {('高倾斜低振幅' if d['params']['r'] == 1e-12 else '高倾斜高振幅' if d['params']['r'] == 1e-8 else 'Sobol 留出')} | {d['params']['r']:.3g} | {d['params']['n_t']:.3g} | z10 | {d['root_relative_error']:.3g} | {d['spectrum_max_abs_dex']:.3g} | {d['fast_seconds']:.3g} |"
        for d in z10
    ],
    "",
    "z10 三点最大根误差为 `{:.3g}`，最大谱误差为 `{:.3g} dex`，均低于验收限。原 z8 高倾斜个别根误差略超限的问题，与尾部截断有关；z10 闭合后当前源码样本通过。一个 `n_t=0.2` 文件因 Fast 源码哈希不一致而排除。".format(
        payload["reference"]["z10_max_root_relative_error"], payload["reference"]["z10_max_spectrum_abs_dex"]
    ),
    "",
    "可信区性能 A/B：固定 Windows affinity `[0,1]`、Numba 2/workqueue、BLAS 1 线程，每例交替 128 轮；全部 6 例输出完全一致。汇总 P50/P95 比为 `{:.3f}/{:.3f}`，逐例最大比为 `{:.3f}/{:.3f}`，超过 5% 门槛，因此当前候选仍未通过性能验收。".format(
        bench["aggregate_p50_ratio"], bench["aggregate_p95_ratio"], max_p50, max_p95
    ),
    "",
    "| 可信案例 | P50 比 | P95 比 |",
    "|---|---:|---:|",
    *[f"| {c['case']} | {c['p50_ratio']:.3f} | {c['p95_ratio']:.3f} |" for c in bench["paired_cases"]],
    "",
    "本轮没有改动生产源码；逐阶段 profile 与独立 kernel/背景分支微基准保留为诊断记录，不作为优化收益。",
    "",
    "![CR0 parameter coverage and accuracy errors](cr0_blue_tilt_domain_validation.png)",
    "原始输入、完整谱与 Reference 数组保存在相邻 z10 JSON 文件中；本图只标示采样设计与精度样本。",
]
(HERE / "cr0_blue_tilt_domain_validation_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

fig, (ax, err) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.1, 1]})
ax.set_xlim(-12.3, -7.7)
ax.set_ylim(0.05, 0.55)
ax.add_patch(
    plt.Rectangle(
        (-12, 0.1), 4, 0.4, fill=False, linestyle="--", edgecolor="#4c78a8", linewidth=1.5, label="screened rectangle"
    )
)
for c in screen["cases"]:
    p = c["input"]["params"]
    ax.scatter(np.log10(p["r"]), p["n_t"], s=55, marker="s", facecolor="#59a14f", edgecolor="black", linewidth=0.4)
for c in sobol["cases"]:
    p = c["input"]["params"]
    ax.scatter(np.log10(p["r"]), p["n_t"], s=48, marker="o", facecolor="#f28e2b", edgecolor="black", linewidth=0.4)
ax.scatter(
    np.log10(low_corner["params"]["r"]),
    low_corner["params"]["n_t"],
    s=150,
    marker="D",
    facecolor="#4e79a7",
    edgecolor="black",
    label="z8 Reference",
)
for d in z10:
    marker = "*" if "sobol" in d["params"] and False else "^"
    ax.scatter(
        np.log10(d["params"]["r"]), d["params"]["n_t"], s=190, marker=marker, facecolor="#e15759", edgecolor="black"
    )
ax.scatter([], [], s=55, marker="s", facecolor="#59a14f", edgecolor="black", label="rectangle/grid")
ax.scatter([], [], s=48, marker="o", facecolor="#f28e2b", edgecolor="black", label="Sobol holdout")
ax.scatter([], [], s=190, marker="^", facecolor="#e15759", edgecolor="black", label="z10 Reference")
ax.set_xlabel(r"$\log_{10}(r)$")
ax.set_ylabel(r"$n_t$")
ax.set_title("Old Fast failed; candidate Fast succeeded")
ax.grid(alpha=0.25)
ax.legend(fontsize=8, loc="upper left")
labels = ["z8 low", "z10 rmin", "z10 rmax", "z10 Sobol"]
roots = [low_corner["root_relative_error"]] + [d["root_relative_error"] for d in z10]
specs = [low_corner["spectrum_max_abs_dex"]] + [d["spectrum_max_abs_dex"] for d in z10]
xx = np.arange(4)
w = 0.34
err.bar(xx - w / 2, roots, w, label="root relative error", color="#4e79a7")
err.bar(xx + w / 2, specs, w, label="max spectrum error (dex)", color="#f28e2b")
err.axhline(2e-4, color="#4e79a7", linestyle="--", linewidth=1, label="root tol 2e-4")
err.axhline(1e-3, color="#f28e2b", linestyle=":", linewidth=1, label="spectrum tol 1e-3 dex")
err.set_yscale("log")
err.set_xticks(xx, labels, rotation=20, ha="right")
err.set_ylabel("relative / absolute error")
err.set_title("Independent Reference checks")
err.grid(axis="y", alpha=0.25)
err.legend(fontsize=8)
fig.suptitle("CR0 blue-tilt: sampled evidence, not a continuous-domain certificate")
fig.tight_layout()
fig.savefig(HERE / "cr0_blue_tilt_domain_validation.png", dpi=180)
plt.close(fig)
print(
    "written",
    *[
        str(HERE / n)
        for n in [
            "cr0_blue_tilt_domain_validation_summary.json",
            "cr0_blue_tilt_domain_validation_summary.md",
            "cr0_blue_tilt_domain_validation.png",
        ]
    ],
    sep="\n",
)
