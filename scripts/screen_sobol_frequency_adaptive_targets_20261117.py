"""Screen existing Fast adaptive frequency targets against direct z10 root anchors."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

SAMPLES = {
    3: ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index3_frequency_refinement_z10_305_20261117.json",
    7: ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index7_frequency_refinement_z10_305_20261117.json",
}
OUTPUT = ROOT / "docs" / "research" / "results" / "sobol_fast_adaptive_target_root_screen_20261117.json"
TARGETS = (3e-4, 2e-4, 1e-4, 5e-5, 1e-5)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def root_anchor(data: dict[str, object]) -> float:
    return float(data["root_iteration"]["iterations"][-1]["mapped_DN_gw"])


def main() -> int:
    partial: dict[str, object] = {
        "schema_version": 1,
        "status": "running_sobol_fast_adaptive_target_root_screen",
        "scope": "For fixed Sobol indices 3 and 7, compare existing Fast freq_grid='adaptive' target settings against independently closed direct Reference z10 root anchors on nested 305-node grids. This screen compares roots only; it does not compare spectra on identical grids or measure performance A/B.",
        "targets": list(TARGETS),
        "samples": {},
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "threads": 1,
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "stiff_SGWB.py": sha(ROOT / "stiffgwpy_fast" / "stiff_SGWB.py"),
            "freq_adaptive.py": sha(ROOT / "stiffgwpy_fast" / "freq_adaptive.py"),
            "scripts": {},
            "inputs": {},
        },
    }
    for index, path in SAMPLES.items():
        data = json.loads(path.read_text(encoding="utf-8"))
        for key, relative in (("fast_sgwb.py", "stiffgwpy_fast/fast_sgwb.py"),
                              ("reference.py", "stiffgwpy_fast/reference.py")):
            if data["source_sha256"][key] != sha(ROOT / relative):
                raise RuntimeError(f"Sobol index {index} source SHA differs for {key}")
        params = {key: float(value) for key, value in data["parameters"].items()}
        anchor = root_anchor(data)
        partial["source_sha256"]["inputs"][str(index)] = sha(path)
        partial["samples"][str(index)] = {
            "parameters": params,
            "reference_z10_root_anchor_305": anchor,
            "reference_grid_count": 305,
            "reference_grid_sha256_f64le": data["frequency_grids"]["refined_305_sha256_f64le"],
            "runs": [],
        }
        for target in TARGETS:
            model = LCDM_SG(**params)
            caught: list[warnings.WarningMessage]
            start = time.perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", RuntimeWarning)
                result = model.SGWB_iter(
                    engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                    threads=1, accuracy_mode="fast", physical_guard=True,
                    freq_grid="adaptive", freq_grid_target=target,
                    freq_grid_max_points=1500,
                )
            runtime = time.perf_counter() - start
            row = {
                "target": target,
                "returned_model": result is model,
                "failure_reason": model.fast_failure_reason,
                "runtime_s": runtime,
                "root_DN_eff": (None if result is not model else float(model.cosmo_param["DN_eff"])),
                "mapped_DN_gw": (None if result is not model else float(model.DN_gw[-1])),
                "fixed_point_residual_relative": (
                    None if result is not model
                    else float(model.DN_gw[-1] / model.cosmo_param["DN_eff"] - 1.0)
                ),
                "root_relative_error_vs_reference_z10_305": (
                    None if result is not model
                    else float(model.cosmo_param["DN_eff"] / anchor - 1.0)
                ),
                "frequency_count": (None if result is not model else int(model.f.size)),
                "frequency_sha256_f64le": (
                    None if result is not model else grid_digest(np.asarray(model.f, dtype=np.float64))
                ),
                "freq_grid_error": getattr(model, "freq_grid_error", None),
                "public_spectrum_finite": (
                    False if result is not model else bool(np.all(np.isfinite(model.log10OmegaGW)))
                ),
                "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught],
            }
            partial["samples"][str(index)]["runs"].append(row)
            partial["source_sha256"]["script"] = sha(Path(__file__))
            OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    partial["status"] = "complete_sobol_fast_adaptive_target_root_screen"
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        index: [
            {"target": row["target"], "count": row["frequency_count"],
             "root_error": row["root_relative_error_vs_reference_z10_305"],
             "grid_error": row["freq_grid_error"]}
            for row in sample["runs"]
        ]
        for index, sample in partial["samples"].items()
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
