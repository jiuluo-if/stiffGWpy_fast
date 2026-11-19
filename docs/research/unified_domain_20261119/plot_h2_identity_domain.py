"""Plot paired failure classes and selected Reference precision samples."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "h2_identity_domain_coverage.png"


def _load(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def _status(row: dict[str, object]) -> str:
    if row.get("success"):
        return "success"
    return str(row.get("failure_reason") or row.get("exception") or "failed")


def _points(artifact: dict[str, object], x_key: str, y_key: str):
    points = []
    for pair in artifact["cases"]:
        inp = pair["input"]["params"]
        x = float(inp[x_key])
        y = float(inp[y_key])
        points.append((x, y, _status(pair["baseline_fast"]), _status(pair["candidate_fast"])))
    return points


def _read_precision() -> list[tuple[str, float, float, str]]:
    results = []
    files = (
        ("cr0 Ωb/T", "validate_h2_identity_region_z8.json"),
        ("cr0 Ωb/Ωch", "validate_h2_identity_density_slice_z8.json"),
        ("cr1 Ωb/Ωch", "validate_h2_identity_cr1_density_z8.json"),
    )
    for label, name in files:
        artifact = _load(name)
        for row in artifact["cases"]:
            if "root_relative_error_proxy" not in row:
                continue
            results.append((
                f"{label} {row['case_id']}",
                float(row["root_relative_error_proxy"]),
                float(row["spectrum_max_abs_dex"]),
                "cr1" if label.startswith("cr1") else "cr0",
            ))
    low_cr1 = _load("validate_h2_identity_cr1_low_rectangle_z8.json")
    mid_cr1 = _load("validate_h2_identity_cr1_mid_rectangle_z8.json")
    upper_cr1 = _load("validate_h2_identity_cr1_upper_boundary_z8.json")
    square_cr1 = _load("validate_h2_identity_cr1_600_square_z8.json")
    rejected_boundary = _load("validate_h2_identity_cr1_700_square_z8.json")
    for row in (
        low_cr1["cases"] + mid_cr1["cases"]
        + upper_cr1["cases"] + square_cr1["cases"]
        + rejected_boundary["cases"]
    ):
        results.append((
            f"cr1 低密度 {row['case_id']}",
            float(row["root_relative_error_proxy"]),
            float(row["spectrum_max_abs_dex"]),
            "cr1",
        ))
    symmetry_holdout = _load(
        "validate_h2_identity_cr1_asymmetric_holdout_symmetry.json"
    )
    results.append((
        "cr1 非对称矩形独立留出 800k",
        float(symmetry_holdout["root_relative_error_proxy"]),
        float(symmetry_holdout["spectrum_max_abs_dex"]),
        "cr1",
    ))
    square_holdout = _load(
        "validate_h2_identity_cr1_600_square_holdout_symmetry.json"
    )
    results.append((
        "cr1 扩展方形独立留出 1100k",
        float(square_holdout["root_relative_error_proxy"]),
        float(square_holdout["spectrum_max_abs_dex"]),
        "cr1",
    ))
    for label, name in (
        ("cr1 exact-quarter 1300k", "validate_exact_quarter_transfer_1300_candidate.json"),
        ("cr1 exact-quarter 1350k", "validate_exact_quarter_transfer_1350_candidate.json"),
        ("cr1 exact-quarter 1400k", "validate_exact_quarter_transfer_1400_candidate.json"),
    ):
        row = _load(name)
        results.append((
            label,
            float(row["fast_root_relative_to_closed_reference"]),
            float(row["spectrum_max_abs_dex_at_fast_root"]),
            "cr1",
        ))
    return results


def main() -> None:
    for font_path in (
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ):
        if font_path.exists():
            plt.rcParams["font.family"] = font_manager.FontProperties(
                fname=str(font_path)
            ).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False
    slices = (
        ("cr=0: Ωb × Tre", "screen_h2_identity_region.json", "Omega_bh2", "T_re"),
        ("cr=0: Ωb × Ωch", "screen_h2_identity_density_slice.json", "Omega_bh2", "Omega_ch2"),
        ("cr=1: Ωb × Ωch", "screen_h2_identity_cr1_density.json", "Omega_bh2", "Omega_ch2"),
    )
    fig = plt.figure(figsize=(13.2, 14.0), constrained_layout=True)
    grid = fig.add_gridspec(5, 2, height_ratios=(1.0, 1.0, 1.0, 0.8, 1.35))
    colors = {
        "success": "#16804a",
        "nonfinite_background_primitive": "#c6423a",
        "resource_timeout": "#e09b22",
        "nonpositive_h2": "#7f4ab2",
        "invalid_cutoff": "#586a80",
    }
    markers = {"baseline": "x", "candidate": "o"}
    for row_index, (title, file, x_key, y_key) in enumerate(slices):
        artifact = _load(file)
        points = _points(artifact, x_key, y_key)
        for col_index, version in enumerate(("baseline", "candidate")):
            ax = fig.add_subplot(grid[row_index, col_index])
            for x, y, old_status, new_status in points:
                status = old_status if version == "baseline" else new_status
                ax.scatter(
                    np.log10(x), np.log10(y), s=74 if version == "baseline" else 65,
                    marker=markers[version], c=colors.get(status, "#777777"),
                    linewidths=2 if version == "baseline" else 0.7,
                    edgecolors="black" if version == "candidate" else None,
                    zorder=3,
                )
            ax.grid(True, alpha=0.25)
            if row_index == 2 and version == "candidate":
                ax.scatter(
                    np.log10(7.0e5), np.log10(7.0e5), marker="D", s=105,
                    facecolors="white", edgecolors="#b23b32", linewidths=2,
                    zorder=5,
                )
                ax.annotate(
                    "Fast成功，根门槛失败", (np.log10(7.0e5), np.log10(7.0e5)),
                    xytext=(-92, 10), textcoords="offset points", fontsize=8,
                    color="#9f2f28",
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
                )
                lower = 5.0
                upper = np.log10(6.0e5)
                ax.add_patch(Rectangle(
                    (lower, lower), upper - lower, upper - lower,
                    fill=False, edgecolor="#16804a", linestyle="--",
                    linewidth=1.8, zorder=4,
                ))
                ax.text(
                    lower + 0.015, upper - 0.02, "精度核验子域",
                    color="#12643b", fontsize=8, va="top",
                    bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
                )
                ax.add_patch(Rectangle(
                    (lower, lower), np.log10(7.0e5) - lower,
                    np.log10(3.0e5) - lower, fill=False,
                    edgecolor="#0d7d87", linestyle=":", linewidth=1.8,
                    zorder=4,
                ))
                ax.add_patch(Rectangle(
                    (lower, lower), np.log10(7.0e5) - lower,
                    np.log10(6.0e5) - lower, fill=False,
                    edgecolor="#2476a6", linestyle="-.", linewidth=1.7,
                    zorder=4,
                ))
                ax.text(
                    5.02, np.log10(3.0e5) - 0.02, "非对称子域",
                    color="#0a6972", fontsize=8, va="top",
                    bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
                )
                ax.text(
                    np.log10(6.0e5) - 0.2, np.log10(6.0e5) - 0.03,
                    "当前接收矩形", color="#1f668f", fontsize=8, va="top",
                    bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
                )
            if row_index == 2:
                lower_patch = np.log10(6.5e5)
                upper_patch = np.log10(7.0e5)
                ax.add_patch(Rectangle(
                    (lower_patch, lower_patch),
                    upper_patch - lower_patch,
                    upper_patch - lower_patch,
                    fill=False,
                    edgecolor="#8a3fa0",
                    linestyle="-",
                    linewidth=2.2,
                    zorder=4,
                ))
                patch = _load("screen_exact_quarter_expanded_patch.json")
                patch_marker = "x" if version == "baseline" else "o"
                patch_color = "#c6423a" if version == "baseline" else "#16804a"
                for pair in patch["paired_rows"]:
                    params = pair["input"]["params"]
                    ax.scatter(
                        np.log10(float(params["Omega_bh2"])),
                        np.log10(float(params["Omega_ch2"])),
                        marker=patch_marker,
                        s=75,
                        color=patch_color,
                        edgecolors="black" if version == "candidate" else None,
                        linewidths=0.7,
                        zorder=6,
                    )
            ax.set_title(f"{title} — {'旧版 Fast' if version == 'baseline' else '候选 Fast'}")
            ax.set_xlabel(f"log10({x_key})")
            ax.set_ylabel(f"log10({y_key})")
            statuses = sorted({_status(r[f"{version}_fast"]) for r in artifact["cases"]})
            counts = {status: sum(_status(r[f"{version}_fast"]) == status for r in artifact["cases"])
                      for status in statuses}
            ax.text(
                0.02, 0.98,
                "\n".join(f"{status}: {count}" for status, count in counts.items()),
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
            )

    patch = _load("screen_exact_quarter_expanded_patch.json")
    low_patch = np.log10(6.5e5)
    high_patch = np.log10(7.0e5)
    for col_index, version in enumerate(("baseline", "candidate")):
        ax = fig.add_subplot(grid[3, col_index])
        ax.add_patch(Rectangle(
            (low_patch, low_patch), high_patch-low_patch, high_patch-low_patch,
            fill=False, edgecolor="#8a3fa0", linewidth=2.0, zorder=3,
        ))
        for pair in patch["paired_rows"]:
            params = pair["input"]["params"]
            is_holdout = pair["input"]["case_id"] == "patch_04"
            ax.scatter(
                np.log10(float(params["Omega_bh2"])),
                np.log10(float(params["Omega_ch2"])),
                marker=("D" if is_holdout else ("x" if version == "baseline" else "o")),
                s=72 if is_holdout else 54,
                color=("#8a3fa0" if is_holdout else ("#c6423a" if version == "baseline" else "#16804a")),
                edgecolors=("black" if version == "candidate" and not is_holdout else None),
                linewidths=0.8,
                zorder=4,
            )
        ax.set_xlim(low_patch - 0.012, high_patch + 0.012)
        ax.set_ylim(low_patch - 0.012, high_patch + 0.012)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.25)
        ax.set_title(f"cr=1 扩展区局部：{'旧版 Fast' if version == 'baseline' else '候选 Fast'}")
        ax.set_xlabel("log10(Omega_bh2)")
        ax.set_ylabel("log10(Omega_ch2)")
        ax.text(
            0.02, 0.98,
            f"紫框：验证矩形；紫菱形：中心留出\n"
            f"成功点：{sum(pair[f'{version}_fast']['success'] for pair in patch['paired_rows'])}/5",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
        )

    precision = _read_precision()
    ax_root = fig.add_subplot(grid[4, 0])
    ax_spec = fig.add_subplot(grid[4, 1])
    labels = [item[0] for item in precision]
    positions = np.arange(len(labels))
    point_colors = ["#b23b32" if item[3] == "cr1" else "#1f6d9c" for item in precision]
    root_errors = [item[1] for item in precision]
    spec_errors = [item[2] for item in precision]
    ax_root.scatter(root_errors, positions, c=point_colors, s=30)
    ax_root.axvline(2e-4, color="#b23b32", linestyle="--", label="根误差门槛 2e-4")
    ax_root.set_xscale("log")
    ax_root.set_yticks(positions, labels, fontsize=7)
    ax_root.invert_yaxis()
    ax_root.set_xlabel("Reference 映射根相对差")
    ax_root.set_title("自洽根精度")
    ax_root.grid(True, axis="x", alpha=0.25)
    ax_root.legend(fontsize=8, loc="lower right")
    ax_spec.scatter(spec_errors, positions, c=point_colors, s=30)
    ax_spec.axvline(1e-3, color="#b23b32", linestyle="--", label="谱误差门槛 1e-3 dex")
    ax_spec.set_xscale("log")
    ax_spec.set_yticks(positions, labels, fontsize=7)
    ax_spec.invert_yaxis()
    ax_spec.set_xlabel("频谱最大绝对差 (dex)")
    ax_spec.set_title("全频谱最大误差")
    ax_spec.grid(True, axis="x", alpha=0.25)
    ax_spec.legend(fontsize=8, loc="lower right")
    fig.suptitle("H² 恒等背景候选：旧失败类型、求解覆盖与 Reference 精度", fontsize=15)
    fig.savefig(OUTPUT, dpi=180)
    print(OUTPUT)


if __name__ == "__main__":
    main()
