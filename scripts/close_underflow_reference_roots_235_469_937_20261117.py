"""Close Reference roots for the current-source direct 235/469/937 grid screen."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts import underflow_grid_convergence_20261115 as grid_runner  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_current_20261117.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_reference_root_closure_235_469_937_current_20261117.json"
ROOT_RESIDUAL_TOLERANCE = 2e-4
MAX_ROOT_MAPS = 5


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes(order="C")).hexdigest()


def spectrum_delta(new: np.ndarray, old: np.ndarray) -> dict[str, float | int]:
    delta = np.abs(np.asarray(new, dtype=np.float64) - np.asarray(old, dtype=np.float64))
    if delta.size == 0 or np.any(~np.isfinite(delta)):
        raise FloatingPointError("Reference root-iteration spectrum contains invalid log values")
    return {
        "count": int(delta.size),
        "max_abs_dex": float(np.max(delta)),
        "p95_abs_dex": float(np.quantile(delta, 0.95)),
    }


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source.get("status") != "complete_fixed_background_direct_nested_grid_screen":
        raise RuntimeError("Input grid screen is not complete")
    for rel, expected in source["source_sha256"].items():
        if sha256(ROOT / rel) != expected:
            raise RuntimeError(f"Input source hash differs for {rel}")

    model = grid_runner._reference_unit_model()
    physical = grid_runner.LCDM_SG(**grid_runner.PARAMS)
    physical.gen_expansion()
    log_amplitude = grid_runner._log_tensor_amplitude(
        physical.cosmo_param, physical.derived_param
    )
    partial: dict[str, object] = {
        "schema_version": 1,
        "status": "running_reference_root_closure_235_469_937_current_source",
        "scope": "Continue direct Cartesian Reference maps at each current-source Fast self-consistent root on the exact nested 235/469/937 grids until root-relative residual closure; each frequency node is solved directly.",
        "input_artifact": INPUT.relative_to(ROOT).as_posix(),
        "parameters": grid_runner.PARAMS,
        "background": source["background"],
        "root_residual_tolerance": ROOT_RESIDUAL_TOLERANCE,
        "max_root_maps": MAX_ROOT_MAPS,
        "source_sha256": {
            **source["source_sha256"],
            "input_artifact": sha256(INPUT),
            "closure_runner": sha256(Path(__file__)),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": grid_runner.THREADS,
        },
        "grids": {},
    }
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    for count in (235, 469, 937):
        key = str(count)
        grid_record = source["frequency_grids"][key]
        freqs = np.asarray(grid_record["log10_hz"], dtype=np.float64)
        if freqs.size != count or grid_sha256(freqs) != grid_record["sha256_f64le"]:
            raise RuntimeError(f"Direct frequency grid {count} failed its hash check")
        fast_root = float(source["fast_self_consistent"][key]["root_DN_eff"])
        initial = source["reference_map_at_fast_root"][key]
        if not np.isclose(float(initial["input_DN_eff"]), fast_root, rtol=2e-14, atol=0.0):
            raise RuntimeError(f"Reference starting map does not use the Fast root for {count}")
        initial_logs = np.asarray(initial["log10_OmegaGW"], dtype=np.float64)
        x = float(initial["mapped_DN_gw"])
        records: list[dict[str, object]] = []
        for iteration in range(1, MAX_ROOT_MAPS + 1):
            started = time.perf_counter()
            mapped = grid_runner._reference_map(model, freqs, x, log_amplitude)
            elapsed = time.perf_counter() - started
            y = float(mapped["mapped_DN_gw"])
            residual = y / x - 1.0
            logs = np.asarray(mapped["log10_OmegaGW"], dtype=np.float64)
            record: dict[str, object] = {
                "iteration": iteration,
                "input_DN_eff": x,
                "mapped_DN_gw": y,
                "relative_map_residual": residual,
                "runtime_s": elapsed,
                "raw_spectrum_positive_finite_nodes": count,
                "runtime_warnings": mapped["runtime_warnings"],
                "normalized_nodes_underflowed_to_zero": mapped["normalized_nodes_underflowed_to_zero"],
                "quadrature_error_relative": mapped["quadrature_error_relative"],
                "interpolation_error_relative": mapped["interpolation_error_relative"],
                "spectrum_change_vs_reference_at_fast_root": spectrum_delta(logs, initial_logs),
                "log10_OmegaGW": logs.tolist(),
            }
            records.append(record)
            partial["grids"][key] = {
                "frequency_count": count,
                "frequency_sha256_f64le": grid_record["sha256_f64le"],
                "fast_root_DN_eff": fast_root,
                "reference_map_at_fast_root": initial,
                "root_DN_eff": y,
                "root_relative_difference_fast_vs_reference": fast_root / y - 1.0,
                "root_iterations": records,
                "closed": abs(residual) <= ROOT_RESIDUAL_TOLERANCE,
            }
            partial["status"] = "running_reference_root_closure_235_469_937_current_source"
            OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            if abs(residual) <= ROOT_RESIDUAL_TOLERANCE:
                break
            x = y
        else:
            partial["status"] = "reference_root_closure_iteration_limit"
            OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            raise RuntimeError(f"Reference root did not close on {count} direct nodes")

    partial["status"] = "complete_reference_root_closure_235_469_937_current_source"
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "roots": {
            k: {
                "root_DN_eff": v["root_DN_eff"],
                "root_relative_difference_fast_vs_reference": v["root_relative_difference_fast_vs_reference"],
                "residual": v["root_iterations"][-1]["relative_map_residual"],
                "maps": len(v["root_iterations"]),
            }
            for k, v in partial["grids"].items()
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
