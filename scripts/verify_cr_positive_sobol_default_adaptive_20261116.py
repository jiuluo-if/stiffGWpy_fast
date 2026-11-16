"""Recheck the Sobol Fast screen with the current default adaptive handoff."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

PREFLIGHT_INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_holdout_20261116.json"
SOBOL1_REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index1_reference_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_default_adaptive_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def reference_record(index: int, preflight: dict[str, object], sobol1: dict[str, object]) -> tuple[dict[str, object] | None, str | None]:
    if index == 0:
        ref = preflight["samples"][index]["direct_reference"]
        if ref["source_sha256"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
            raise RuntimeError("Sobol point 0 Reference SHA differs from current Reference")
        return ref, "docs/research/results/cr_positive_sobol_cutoff_holdout_20261116.json"
    if index == 1:
        ref = sobol1["reference"]
        if ref["source_sha256"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
            raise RuntimeError("Sobol point 1 Reference SHA differs from current Reference")
        return ref, "docs/research/results/cr_positive_sobol_index1_reference_20261116.json"
    return None, None


def summarize_reference(
    ref: dict[str, object], fast_logs: np.ndarray, freq: np.ndarray, current_root: float
) -> dict[str, object]:
    ref_logs = np.asarray([
        np.nan if value is None else float(value)
        for value in ref["z_tails"]["12.0"]["raw_log10_OmegaGW"]
    ], dtype=np.float64)
    valid = np.isfinite(ref_logs) & np.isfinite(fast_logs)
    if not valid.any():
        raise RuntimeError("direct Reference z12 has no comparable positive finite spectrum nodes")
    delta = fast_logs[valid] - ref_logs[valid]
    worst = int(np.argmax(np.abs(delta)))
    indices = np.flatnonzero(valid)
    return {
        "reference_artifact_source": None,
        "reference_root_DN_eff": float(ref["fixed_DN_eff"]),
        "frequency_sha256_f64le": ref["frequency_sha256_f64le"],
        "reference_z8_to_z10_map_relative_change": float(ref["z8_to_z10_map_relative_change"]),
        "reference_z10_to_z12_map_relative_change": float(ref["z10_to_z12_map_relative_change"]),
        "reference_z10_to_z12_max_abs_dex": float(ref["z10_to_z12_max_abs_dex"]),
        "fast_root_relative_error_vs_reference_z12_map": float(
            current_root / float(ref["z_tails"]["12.0"]["DN_gw"]) - 1.0
        ),
        "fast_default_adaptive_vs_reference_z12": {
            "compared_nodes": int(np.count_nonzero(valid)),
            "uncompared_nodes": int(freq.size - np.count_nonzero(valid)),
            "max_abs_dex": float(np.max(np.abs(delta))),
            "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
            "worst_log10f": float(freq[indices[worst]]),
            "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
        },
        "reference_warnings": {
            z: row["runtime_warnings"]
            for z, row in ref["z_tails"].items()
        },
    }


def main() -> int:
    preflight = json.loads(PREFLIGHT_INPUT.read_text(encoding="utf-8"))
    sobol1 = json.loads(SOBOL1_REFERENCE_INPUT.read_text(encoding="utf-8"))
    for name, source in (
        ("Sobol preflight", preflight),
        ("Sobol index 1 Reference", sobol1),
    ):
        source_sha = source["source_sha256"]
        if source_sha["fast_sgwb.py"] != sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"):
            raise RuntimeError(f"{name} Fast source SHA differs from current Fast")
        if source_sha["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
            raise RuntimeError(f"{name} Reference source SHA differs from current Reference")

    rows = []
    for original in preflight["samples"]:
        index = int(original["sobol_index"])
        params = {key: float(value) for key, value in original["parameters"].items()}
        model = LCDM_SG(**params)
        original_goal = FA.goal_oriented_freqs
        original_construct = FA.construct_fast_frequency_grid
        start = time.perf_counter()
        caught: list[warnings.WarningMessage]
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", RuntimeWarning)
                result = model.SGWB_iter(
                    engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                    threads=1, accuracy_mode="fast", freq_grid="goal",
                    physical_guard=True,
                )
        finally:
            FA.goal_oriented_freqs = original_goal
            FA.construct_fast_frequency_grid = original_construct
        if result is not model:
            raise RuntimeError(f"Sobol point {index} current Fast failed: {model.fast_failure_reason}")
        runtime = time.perf_counter() - start
        freq = np.asarray(model.f, dtype=np.float64)
        logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
        grid_hash = digest(freq)
        prior = original["fast"]
        if grid_hash != prior["frequency_sha256_f64le"]:
            raise RuntimeError(f"Sobol point {index} current default Fast grid changed from its preflight grid")
        root = float(model.cosmo_param["DN_eff"])
        if not np.isclose(root, float(prior["root_DN_eff"]), rtol=2e-12, atol=0.0):
            raise RuntimeError(f"Sobol point {index} current default Fast root changed from its preflight root")
        selected = int(getattr(model, "_fast_adaptive_handoff_selected", 0))
        fraction = float(getattr(model, "_fast_adaptive_handoff_fraction", 0.0))
        ref, ref_path = reference_record(index, preflight, sobol1)
        ref_metrics = None
        if ref is not None:
            ref_metrics = summarize_reference(ref, logs, freq, root)
            ref_metrics["reference_artifact_source"] = ref_path
            ref_hash = digest(freq)
            if ref_hash != ref["frequency_sha256_f64le"]:
                raise RuntimeError(f"Sobol point {index} direct Reference grid changed")
        rows.append({
            "sobol_index": index,
            "unit_point": original["unit_point"],
            "parameters": params,
            "cutoff": original["cutoff"],
            "supplied_n_t": params["n_t"],
            "derived_n_t": float(model.derived_param["nt"]),
            "stored_A_t": float(model.derived_param["A_t"]),
            "fast_default_adaptive": {
                "runtime_s": runtime,
                "root_DN_eff": root,
                "preflight_root_DN_eff": float(prior["root_DN_eff"]),
                "root_change_vs_preflight": root / float(prior["root_DN_eff"]) - 1.0,
                "mapped_DN_gw": float(model.DN_gw[-1]),
                "fixed_point_residual_relative": float(model.DN_gw[-1] / root - 1.0),
                "frequency_count": int(freq.size),
                "frequency_sha256_f64le": grid_hash,
                "adaptive_selected_count": selected,
                "adaptive_selected_fraction": fraction,
                "public_spectrum_finite": bool(np.isfinite(logs).all()),
                "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught],
                "log10_OmegaGW": logs.tolist(),
            },
            "direct_reference_comparison": ref_metrics,
            "status": "fast_default_adaptive_reference_compared" if ref_metrics else "fast_default_adaptive_reference_pending",
        })

    output = {
        "schema_version": 1,
        "status": "complete_current_default_adaptive_sobol_preflight_and_two_reference_holdouts",
        "scope": "Recheck all eight fixed Sobol cr>0 cutoff samples with the current default adaptive Fast path. Confirm exact root/grid replay against the first screen, and compare the two archived direct Reference z12 holdouts without interpolating spectra between nodes. The legacy preflight field named fast_z5_spectrum was generated with the default adaptive handoff active; this corrected artifact labels the current Fast output explicitly.",
        "design": preflight["design"],
        "preflight_input_sha256": sha(PREFLIGHT_INPUT),
        "sobol_index1_reference_input_sha256": sha(SOBOL1_REFERENCE_INPUT),
        "samples": rows,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "script": sha(Path(__file__)),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "sample_count": len(rows),
        "reference_compared_indices": [row["sobol_index"] for row in rows if row["direct_reference_comparison"]],
        "samples": [{
            "sobol_index": row["sobol_index"],
            "selected": row["fast_default_adaptive"]["adaptive_selected_count"],
            "root_change": row["fast_default_adaptive"]["root_change_vs_preflight"],
            "reference_error": None if not row["direct_reference_comparison"] else row["direct_reference_comparison"]["fast_default_adaptive_vs_reference_z12"]["max_abs_dex"],
        } for row in rows],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
