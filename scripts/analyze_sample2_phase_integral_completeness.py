"""Audit the 787/800 validated Phase 6c modes and the unresolved integral mass.

The oracle-filled hybrid integral is sensitivity analysis only; Reference
values are never counted as Phase candidate outputs.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from scripts._resource_budget import apply_environment, telemetry
except ImportError:  # pragma: no cover
    from _resource_budget import apply_environment, telemetry

apply_environment()

import numpy as np  # noqa: E402

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast import reference as R  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PHASE = ROOT / "docs" / "sample2_phase_reference_root_res8_tail5_gamma0_20261105.json"
OFFSETS = ROOT / "docs" / "sample2_phase_start_offset_probe_20261105.json"
DUAL = ROOT / "docs" / "sample2_phase_residual_orders_dual_reference_20261105.json"
REF_JSON = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.json"
REF_NPZ = ROOT / "docs" / "sample2_reference_root_res8_tail5_gamma0_20261105.npz"
OUTPUT = ROOT / "docs" / "sample2_phase_integral_completeness_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _write(record: dict) -> None:
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    phase_data = json.loads(PHASE.read_text(encoding="utf-8"))
    offset_data = json.loads(OFFSETS.read_text(encoding="utf-8"))
    dual_data = json.loads(DUAL.read_text(encoding="utf-8"))
    ref_data = json.loads(REF_JSON.read_text(encoding="utf-8"))
    with np.load(REF_NPZ) as archive:
        freqs = np.asarray(archive["freqs"], dtype=np.float64)
        ref_logs = np.asarray(archive["reference_log10OmegaGW"], dtype=np.float64)
        ref_used_tail = np.asarray(archive["reference_used_tail"], dtype=bool)
        ref_omega = np.asarray(archive["reference_Ogw"] - archive["reference_Oj"], dtype=np.float64)
    n = freqs.size
    if n != 800 or ref_logs.shape != (n,) or ref_omega.shape != (n,):
        raise RuntimeError("sample2 Reference archive has an unexpected shape")
    expected_freq_sha = hashlib.sha256(freqs.tobytes()).hexdigest()
    if expected_freq_sha != phase_data["frequency_sha256"]:
        raise RuntimeError("phase and Reference frequency arrays do not match")
    if hashlib.sha256(REF_NPZ.read_bytes()).hexdigest() != phase_data["reference_npz_sha256"]:
        raise RuntimeError("phase artifact points to a different Reference NPZ")

    selected: dict[int, dict] = {}
    sources: dict[str, set[int]] = {"gamma0_fullgrid_order32": set(), "residual_dual_reference": set(), "phase_start_offsets": set()}

    def add(index: int, log_omega: float, order: int, source: str, tail_used: bool) -> None:
        index = int(index)
        if not 0 <= index < n:
            raise RuntimeError(f"out-of-range mode index {index}")
        if not np.isfinite(log_omega):
            raise FloatingPointError(f"non-finite candidate log spectrum at mode {index}")
        if bool(tail_used) != bool(ref_used_tail[index]):
            raise RuntimeError(f"candidate/Reference tail semantics differ at mode {index}")
        sources[source].add(index)
        if index not in selected:
            selected[index] = {
                "log10_omega": float(log_omega),
                "order": int(order),
                "source": source,
                "tail_used": bool(tail_used),
            }

    for row in phase_data["rows"]:
        if row["status"] != "completed":
            continue
        order32 = next(item for item in row["orders"] if int(item["order"]) == 32)
        add(row["index"], order32["log10_omega"], 32, "gamma0_fullgrid_order32", order32["tail_used"])

    for row in dual_data["rows"]:
        if row.get("validated") is not True:
            continue
        accepted = [
            item for item in row["candidate_orders"]
            if item.get("validation", {}).get("accepted") is True
        ]
        if not accepted:
            continue
        item = min(accepted, key=lambda value: int(value["order"]))
        candidate = item["candidate"]
        add(row["index"], candidate["log10_omega"], item["order"], "residual_dual_reference", candidate["tail_used"])

    for row in offset_data["rows"]:
        if row.get("validated") is not True:
            continue
        candidate = row["selected_candidate"]["candidate"]
        add(row["index"], candidate["log10_omega"], candidate["order"], "phase_start_offsets", candidate["tail_used"])

    missing = sorted(set(range(n)) - set(selected))
    phase_rows = {int(row["index"]): row for row in phase_data["rows"]}
    baseline_complete_rows = [row for row in phase_data["rows"] if row["status"] == "completed"]
    missing_classification = []
    for index in missing:
        row = phase_rows[index]
        messages = [str(item.get("message", "")) for item in row.get("attempts", [])]
        is_nonpositive_q = any("normal-form Q" in message or "Q<=0" in message for message in messages)
        missing_classification.append({
            "index": index,
            "frequency_log10_hz": float(freqs[index]),
            "reference_log10_omega": float(ref_logs[index]),
            "reference_used_tail": bool(ref_used_tail[index]),
            "phase_failure_class": "nonpositive_or_nonfinite_Q" if is_nonpositive_q else "residual_or_other_unresolved",
        })

    if len(selected) == 787 and len(missing) != 13:
        raise RuntimeError(f"expected 787 validated modes / 13 unresolved, got {len(selected)}/{len(missing)}")
    candidate_logs = np.full(n, np.nan, dtype=np.float64)
    for index, row in selected.items():
        candidate_logs[index] = row["log10_omega"]
    supported_error = np.abs(candidate_logs[np.isfinite(candidate_logs)] - ref_logs[np.isfinite(candidate_logs)])

    # The Reference-filled values in this hybrid are strictly a sensitivity
    # analysis of the 787 candidate modes, never a Phase candidate map.
    hybrid_logs = np.where(np.isfinite(candidate_logs), candidate_logs, ref_logs)
    hybrid_omega = np.power(10.0, hybrid_logs)
    hybrid_g2, hybrid_qerr, hybrid_ierr = R.integrate_spectrum(freqs, hybrid_omega, np.zeros(n))
    omega_nu = gp.Omega_nh2 / (LCDM_SG(**phase_data["params"]).derived_param["h"] ** 2)
    hybrid_dn = float(gp.Neff0 * hybrid_g2 / omega_nu)

    order = np.argsort(freqs)
    x = freqs[order]
    y = ref_omega[order]
    miss_sorted = ~np.isfinite(candidate_logs[order])
    widths = np.diff(x)
    trapezoid_weights = np.empty(n, dtype=np.float64)
    trapezoid_weights[0] = 0.5 * widths[0]
    trapezoid_weights[-1] = 0.5 * widths[-1]
    trapezoid_weights[1:-1] = 0.5 * (widths[:-1] + widths[1:])
    unresolved_mass_fraction_trapezoid = float(
        np.sum((trapezoid_weights * y)[miss_sorted]) / np.sum(trapezoid_weights * y)
    )

    # Measure the current production Fast cost on the same cosmology. This is
    # a public Fast solve at its formal grid/tolerance, not a test or oracle run.
    model = LCDM_SG(**phase_data["params"])
    fast_started = time.perf_counter()
    with warnings.catch_warnings(record=True) as fast_warnings:
        warnings.simplefilter("always")
        fast_result = model.SGWB_iter(
            engine="fast",
            fallback=False,
            tol=1e-6,
            z_tail=5.0,
            freq_res=1.0,
            threads=1,
            accuracy_mode="fast",
            physical_guard=False,
        )
    fast_runtime = time.perf_counter() - fast_started
    fast_success = fast_result is model
    fast_dn = float(model.DN_gw[-1]) if fast_success else None
    fast_frequency_count = int(model.f.size) if fast_success else None
    order32_supported_work_estimate = float(sum(
        float(row["low_seconds"])
        + float(next(item for item in row["orders"] if int(item["order"]) == 32)["total_seconds"])
        for row in baseline_complete_rows
    ))

    record = {
        "schema_version": 1,
        "experiment": "sample2 phase mode-union, unresolved-integral and Fast-cost audit",
        "status": "complete",
        "source_commit": _git("rev-parse", "HEAD"),
        "branch": _git("branch", "--show-current"),
        "reference_root": float(ref_data["params"]["DN_eff"]),
        "frequency_count": int(n),
        "frequency_sha256": expected_freq_sha,
        "coverage": {
            "gamma0_fullgrid_strict_modes": len(sources["gamma0_fullgrid_order32"]),
            "residual_dual_reference_modes": len(sources["residual_dual_reference"]),
            "phase_start_offset_validated_modes": len(sources["phase_start_offsets"]),
            "unique_validated_modes": len(selected),
            "unresolved_mode_count": len(missing),
            "unresolved_modes": missing_classification,
            "source_mode_indices": {key: sorted(value) for key, value in sources.items()},
        },
        "supported_candidate_spectrum": {
            "max_abs_dex_error_vs_cartesian": float(np.max(supported_error)),
            "p95_abs_dex_error_vs_cartesian": float(np.percentile(supported_error, 95)),
            "candidate_order_counts": {
                str(k): int(sum(row["order"] == k for row in selected.values()))
                for k in sorted({row["order"] for row in selected.values()})
            },
        },
        "unresolved_reference_weight": {
            "trapezoid_nodal_mass_fraction": unresolved_mass_fraction_trapezoid,
            "interpretation": "Observed Reference nodal weight for missing modes; not an a-priori error bound for PCHIP interpolation.",
        },
        "oracle_filled_hybrid_sensitivity_only": {
            "DN_gw": hybrid_dn,
            "relative_delta_vs_reference_root": abs(hybrid_dn - float(ref_data["reference_DN_gw"])) / max(abs(float(ref_data["reference_DN_gw"])), np.finfo(float).tiny),
            "quadrature_error": float(hybrid_qerr),
            "interpolation_error": float(hybrid_ierr),
            "uses_reference_values_for_unresolved_modes": len(missing),
            "counts_as_candidate_map": False,
        },
        "cost": {
            "phase_gamma0_fullgrid_runtime_s": float(phase_data["phase_runtime_s"]),
            "phase_completed_modes": int(phase_data["mode_count_completed"]),
            "single_order32_supported_mode_work_estimate_s": order32_supported_work_estimate,
            "single_order32_estimate_limit": "Estimated low-frequency DOP853 + order32 phase work on the 737 gamma0-completed modes only; excludes 63 failed-mode attempts, the 50 separately screened modes, and any additional recovery work.",
            "cartesian_reference_runtime_s": float(ref_data["runtime_s"]),
            "phase_to_cartesian_reference_runtime_ratio": float(phase_data["phase_runtime_s"] / ref_data["runtime_s"]),
            "fast_public_solve_success": bool(fast_success),
            "fast_public_runtime_s": fast_runtime,
            "fast_frequency_count": fast_frequency_count,
            "fast_DN_gw": fast_dn,
            "fast_warnings": [str(item.message) for item in fast_warnings],
            "phase_to_fast_runtime_ratio": float(phase_data["phase_runtime_s"] / fast_runtime),
            "single_order32_supported_work_to_fast_runtime_ratio": float(order32_supported_work_estimate / fast_runtime),
            "fast_configuration": {"engine": "fast", "accuracy_mode": "fast", "freq_res": 1.0, "z_tail": 5.0, "tol": 1e-6, "threads": 1},
        },
        "source_hashes": {
            "phase_artifact": _sha(PHASE),
            "offset_artifact": _sha(OFFSETS),
            "dual_artifact": _sha(DUAL),
            "reference_npz": _sha(REF_NPZ),
            "reference.py": _sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fast_sgwb.py": _sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "driver.py": _sha(Path(__file__).resolve()),
        },
        "limitations": [
            "The 787-mode union combines three-order grid results with separately dual-validated residual and phase-start candidates; it is a support audit, not one uniformly configured solver sweep.",
            "The 13 missing modes prevent a complete Phase candidate integral map. The oracle-filled hybrid replaces them with Reference values and is sensitivity analysis only.",
            "The measured gamma-zero phase run is slower than the 800-mode Cartesian Reference and is orders slower than the formal Fast public solve; it does not pass the Fast cost gate.",
            "No production or Reference source was changed.",
        ],
        "resources": telemetry(workers=1, threads=1),
    }
    _write(record)
    print(json.dumps({
        "coverage": record["coverage"]["unique_validated_modes"],
        "unresolved": record["coverage"]["unresolved_mode_count"],
        "oracle_filled_hybrid_relative_delta": record["oracle_filled_hybrid_sensitivity_only"]["relative_delta_vs_reference_root"],
        "phase_to_cartesian_ratio": record["cost"]["phase_to_cartesian_reference_runtime_ratio"],
        "phase_to_fast_ratio": record["cost"]["phase_to_fast_runtime_ratio"],
        "fast_runtime_s": fast_runtime,
        "output": str(OUTPUT),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
