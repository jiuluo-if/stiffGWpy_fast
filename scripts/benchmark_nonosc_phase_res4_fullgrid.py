"""Screen the phase-transfer prototype against the archived 260-mode res4 grid.

中文：逐模态记录相位原型的精度、成本和明确失败分类；不替代 Cartesian Reference。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import prototype_nonosc_phase_res4_mode as phase  # noqa: E402
from scripts._resource_budget import telemetry  # noqa: E402

ARCHIVE = ROOT / "docs" / "reference_prufer_cartesian_paired_fullgrid_20261104.json"
OUTPUT = ROOT / "docs" / "nonosc_phase_res4_fullgrid_20261104.json"
PHASE_STEP = 0.5
PHASE_START_MIN = 0.0
PHASE_START_MAX = 10.0


def run_mode(row, params):
    freq = float(row["freq_log10_hz"])
    result = {
        "index": int(row["index"]),
        "freq_log10_hz": freq,
        "reference_log10_omega": float(row["cartesian_log10_omega"]),
        "reference_used_tail": bool(row["used_tail"]),
        "attempts": [],
    }
    phase_start = 5.0
    tried_starts = set()
    while phase_start not in tried_starts and PHASE_START_MIN <= phase_start <= PHASE_START_MAX:
        tried_starts.add(phase_start)
        phase.configure_case(params["DN_eff"], freq, phase_start)
        low_started = time.perf_counter()
        try:
            lo = phase.solve_ivp(
                phase.low_rhs,
                (phase.NSTART, phase.bg.N_inf),
                (phase.Z0, 1.0, 0.0),
                method="DOP853",
                rtol=1e-12,
                atol=[1e-12, 1e-13, 1e-13],
                events=[phase.low_event],
                dense_output=False,
            )
            if not lo.success or not lo.t_events[0].size:
                result["attempts"].append(
                    {"phase_start": phase_start, "status": "handoff_not_reached", "message": lo.message}
                )
                phase_start -= PHASE_STEP
                continue
            low_s = time.perf_counter() - low_started
            n_phase = float(lo.t_events[0][0])
            low_state = lo.y_events[0][0].tolist()
            orders = []
            for order in (16, 24, 32):
                t0 = time.perf_counter()
                state, bounds, segments, phase_s, tail_used, zend = phase.candidate(
                    order, low_state, n_phase, phase_start
                )
                end_n = bounds[-1][0]
                power = (
                    phase.spectrum_from_tail(state, end_n, zend)
                    if tail_used
                    else phase.spectrum_today(state, zend)
                )
                orders.append(
                    {
                        "order": order,
                        "tail_used": tail_used,
                        "end_N": end_n,
                        "end_z": zend,
                        "log10_omega": math.log10(power["OmegaGW"]),
                        "phase_seconds": phase_s,
                        "total_seconds": time.perf_counter() - t0 + low_s,
                        "segments": segments,
                    }
                )
            for item in orders:
                item["abs_dex_error_vs_cartesian"] = abs(
                    item["log10_omega"] - result["reference_log10_omega"]
                )
            result.update(
                status="completed",
                phase_start=phase_start,
                tail_semantics_match=bool(row["used_tail"]) == bool(orders[0]["tail_used"]),
                low_rhs_nfev=int(lo.nfev),
                low_seconds=low_s,
                orders=orders,
                max_dense_residual_scaled=max(
                    (seg["dense_residual_scaled"] for item in orders for seg in item["segments"]),
                    default=None,
                ),
            )
            return result
        except ValueError as exc:
            message = str(exc)
            if "high-frequency screen" in message:
                result["attempts"].append(
                    {"phase_start": phase_start, "status": "screen_rejected", "message": message}
                )
                phase_start += PHASE_STEP
                continue
            result["attempts"].append(
                {"phase_start": phase_start, "status": "value_error", "message": message}
            )
            break
        except ArithmeticError as exc:
            message = str(exc)
            if "Riccati residual" in message:
                result["attempts"].append(
                    {"phase_start": phase_start, "status": "residual_rejected", "message": message}
                )
                phase_start += PHASE_STEP
                continue
            result["attempts"].append(
                {
                    "phase_start": phase_start,
                    "status": "numerical_exception",
                    "exception_type": type(exc).__name__,
                    "message": message,
                }
            )
            break
        except Exception as exc:  # persist unexpected numerical failures instead of hiding them
            result["attempts"].append(
                {
                    "phase_start": phase_start,
                    "status": "numerical_exception",
                    "exception_type": type(exc).__name__,
                    "message": str(exc),
                }
            )
            break
    result.setdefault("status", "unsupported_or_failed")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--indices",
        help="comma-separated archive indices for a reproducible subset; default runs all modes",
    )
    parser.add_argument(
        "--minimum-gamma",
        type=float,
        default=phase.PHASE_GAMMA_MIN,
        help="minimum phase-span screen; use 0 only for a separately recorded residual-validation experiment",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not math.isfinite(args.minimum_gamma) or args.minimum_gamma < 0.0:
        raise ValueError("minimum-gamma must be finite and nonnegative")
    phase.PHASE_GAMMA_MIN = args.minimum_gamma
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    ref_module = ROOT / "stiffgwpy_fast" / "reference.py"
    module_sha = hashlib.sha256(ref_module.read_bytes()).hexdigest()
    if module_sha != archive["reference_module_sha256"]:
        raise RuntimeError("Cartesian archive and current Reference module hashes differ")
    params = dict(phase.PARAMS)
    params.update(archive["parameters"])
    for name, value in phase.PARAMS.items():
        if params.get(name) != value:
            raise RuntimeError("phase prototype and Cartesian archive parameters differ at %s" % name)
    if float(archive["z_tail"]) != phase.ZTAIL:
        raise RuntimeError("phase prototype and Cartesian archive tail depths differ")
    rows = archive["rows"]
    if args.indices:
        selected = {int(value) for value in args.indices.split(",")}
        rows = [row for row in rows if int(row["index"]) in selected]
        if {int(row["index"]) for row in rows} != selected:
            raise ValueError("requested archive index is absent")
    rows = sorted(rows, key=lambda row: int(row["index"]))
    started = time.perf_counter()
    results = [run_mode(row, params) for row in rows]
    completed = [row for row in results if row["status"] == "completed"]
    rows_by_index = {int(row["index"]): row for row in rows}
    order32_seconds = [row["orders"][2]["total_seconds"] for row in completed]
    paired_cartesian_seconds = [float(rows_by_index[row["index"]]["cartesian_runtime_s"]) for row in completed]
    payload = {
        "schema_version": 1,
        "experiment": "newton_kantorovich_phase_transfer_vs_cartesian_res4_grid",
        "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "phase_script_sha256": hashlib.sha256(Path(phase.__file__).read_bytes()).hexdigest(),
        "reference_module_sha256": module_sha,
        "reference_archive": str(ARCHIVE.relative_to(ROOT)),
        "reference_archive_sha256": hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
        "parameters": params,
        "z_tail": archive["z_tail"],
        "solver_settings": {
            "low_frequency_method": "DOP853",
            "low_frequency_rtol": 1e-12,
            "low_frequency_atol": [1e-12, 1e-13, 1e-13],
            "phase_orders": [16, 24, 32],
            "newton_update_tolerance": 2e-14,
            "maximum_newton_steps": 20,
            "collocation_residual_limit": 1e-10,
            "dense_residual_limit": 1e-9,
            "minimum_phase_interval_gamma": phase.PHASE_GAMMA_MIN,
            "phase_start_initial": 5.0,
            "phase_start_step": PHASE_STEP,
            "phase_start_limits": [PHASE_START_MIN, PHASE_START_MAX],
        },
        "phase_start_initial": 5.0,
        "phase_start_step": PHASE_STEP,
        "phase_start_limits": [PHASE_START_MIN, PHASE_START_MAX],
        "mode_count_requested": len(rows),
        "mode_count_completed": len(completed),
        "mode_count_unsupported_or_failed": len(results) - len(completed),
        "max_order32_abs_dex_error": max(
            (row["orders"][2]["abs_dex_error_vs_cartesian"] for row in completed), default=None
        ),
        "p95_order32_abs_dex_error": (
            float(np.percentile([row["orders"][2]["abs_dex_error_vs_cartesian"] for row in completed], 95))
            if completed
            else None
        ),
        "order32_candidate_seconds_p50": float(np.percentile(order32_seconds, 50)) if completed else None,
        "order32_candidate_seconds_p95": float(np.percentile(order32_seconds, 95)) if completed else None,
        "order32_candidate_seconds_sum": float(sum(order32_seconds)),
        "paired_cartesian_seconds_sum": float(sum(paired_cartesian_seconds)),
        "aggregate_speedup_vs_archived_cartesian": (
            float(sum(paired_cartesian_seconds) / sum(order32_seconds)) if order32_seconds else None
        ),
        "tail_semantics_mismatch_count": sum(not row["tail_semantics_match"] for row in completed),
        "total_wall_seconds": time.perf_counter() - started,
        "resources": telemetry(workers=1, threads=1),
        "rows": results,
        "limitations": [
            "The Cartesian comparison is an archived exact-grid baseline; its module hash is checked before execution.",
            "This is a single-parameter res4 grid and does not establish parameter-domain reliability.",
            "Unsupported or failed modes remain explicitly classified and are not counted as successful accuracy results.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
