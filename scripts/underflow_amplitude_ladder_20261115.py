"""Probe cr=0 Fast scaling across represented and underflowed amplitudes."""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.freq_adaptive import _log_tensor_amplitude  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1e-316,
    "n_t": 17.38,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
AMPLITUDES = (1e-316, 1e-317, 1e-318, 1e-320, 5e-324, 0.0)
BASE_RESULT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "underflow_amplitude_ladder_20261115.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_case(
    r: float, freqs: np.ndarray, *, frozen_dn: float | None = None
) -> dict[str, object]:
    params = {**BASE_PARAMS, "r": r}
    model = LCDM_SG(**params)
    if frozen_dn is not None:
        model.cosmo_param["DN_eff"] = frozen_dn
        model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = frozen_dn
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    if frozen_dn is not None:
        model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    error = None
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False,
            tol=0.99 if frozen_dn is not None else 1e-7, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=False,
        )
    except Exception as exc:
        result = None
        error = f"{type(exc).__name__}: {exc}"
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    elapsed = time.perf_counter() - started
    log_amp = float(_log_tensor_amplitude(model.cosmo_param, model.derived_param))
    spectrum = getattr(model, "log10OmegaGW", None)
    mapped = getattr(model, "DN_gw", None)
    root = float(model.cosmo_param["DN_eff"])
    evaluation_background = root if frozen_dn is None else frozen_dn
    map_end = None if mapped is None else float(np.asarray(mapped)[-1])
    finite_spectrum = None if spectrum is None else np.isfinite(spectrum)
    return {
        "solve_mode": "fixed_background" if frozen_dn is not None else "self_consistent",
        "frozen_DN_eff": frozen_dn,
        "r": r,
        "A_t_stored": float(model.derived_param["A_t"]),
        "log_A_t": log_amp if math.isfinite(log_amp) else None,
        "exact_zero": r == 0.0,
        "fast_returned_model": result is not None,
        "fast_failure_reason": model.fast_failure_reason,
        "error": error,
        "frequency_count": int(freqs.size),
        "root_DN_eff": root if result is not None and frozen_dn is None else None,
        "DN_eff_after_call": root,
        "evaluation_background_DN_eff": evaluation_background,
        "mapped_DN_gw": map_end,
        "map_residual_vs_evaluation_background_relative": (
            None if map_end is None or evaluation_background == 0.0
            else (map_end - evaluation_background) / abs(evaluation_background)
        ),
        "root_map_residual_relative": (
            None if frozen_dn is not None or map_end is None or root == 0.0
            else (map_end - root) / abs(root)
        ),
        "runtime_s": elapsed,
        "finite_log_spectrum_count": (
            None if finite_spectrum is None
            else int(np.count_nonzero(finite_spectrum))
        ),
        "log10_OmegaGW": (
            None if spectrum is None else [
                float(value) if np.isfinite(value) else None
                for value in np.asarray(spectrum, dtype=np.float64)
            ]
        ),
    }


def main() -> int:
    baseline = json.loads(BASE_RESULT.read_text(encoding="utf-8"))
    freqs = np.asarray(baseline["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)
    self_consistent = [run_case(r, freqs) for r in AMPLITUDES]
    frozen_dn = float(self_consistent[0]["root_DN_eff"])
    fixed_background = [
        run_case(r, freqs, frozen_dn=frozen_dn) for r in AMPLITUDES
    ]
    rows = self_consistent + fixed_background
    for group in (self_consistent, fixed_background):
        reference = next(row for row in group if row["r"] == BASE_PARAMS["r"])
        reference_spectrum = np.asarray(reference["log10_OmegaGW"], dtype=np.float64)
        reference_map = float(reference["mapped_DN_gw"])
        log10_ref_amp = math.log10(BASE_PARAMS["r"]) + math.log10(BASE_PARAMS["A_s"])
        for row in group:
            spectrum = row["log10_OmegaGW"]
            if (
                spectrum is None
                or row["exact_zero"]
                or any(value is None for value in spectrum)
            ):
                row["normalized_spectrum_max_delta_dex"] = None
                row["map_amplitude_scale_ratio"] = None
                row["root_amplitude_scale_ratio"] = None
                continue
            log10_amp = math.log10(float(row["r"])) + math.log10(BASE_PARAMS["A_s"])
            normalized = np.asarray(spectrum, dtype=np.float64) - log10_amp
            reference_normalized = reference_spectrum - log10_ref_amp
            row["normalized_spectrum_max_delta_dex"] = float(
                np.max(np.abs(normalized - reference_normalized))
            )
            row["map_amplitude_scale_ratio"] = (
                float(row["mapped_DN_gw"]) / reference_map
                / (float(row["r"]) / BASE_PARAMS["r"])
            )
            reference_root = reference.get("root_DN_eff")
            row["root_amplitude_scale_ratio"] = (
                None
                if row["root_DN_eff"] is None or reference_root is None
                else float(row["root_DN_eff"]) / float(reference_root)
                / (float(row["r"]) / BASE_PARAMS["r"])
            )
    output = {
        "schema_version": 1,
        "status": "complete_cr0_amplitude_ladder_fixed_235_grid",
        "scope": "One discrete cr=0, n_t=17.38 point; direct Fast self-consistent roots and fixed-background maps on the exact native 235-node grid from the prior convergence experiment.",
        "parameters_except_r": {key: value for key, value in BASE_PARAMS.items() if key != "r"},
        "amplitude_ladder_r": list(AMPLITUDES),
        "frequency_log10_hz": freqs.tolist(),
        "frequency_sha256_f64le": hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest(),
        "cases": rows,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "threads": 1,
        },
        "source_sha256": {
            name: sha(ROOT / path)
            for name, path in {
                "fast_sgwb.py": "stiffgwpy_fast/fast_sgwb.py",
                "freq_adaptive.py": "stiffgwpy_fast/freq_adaptive.py",
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "diagnostic_script": "scripts/underflow_amplitude_ladder_20261115.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "cases": [{
            "solve_mode": row["solve_mode"],
            "r": row["r"], "A_t_stored": row["A_t_stored"],
            "root_DN_eff": row["root_DN_eff"],
            "root_map_residual_relative": row["root_map_residual_relative"],
            "normalized_spectrum_max_delta_dex": row["normalized_spectrum_max_delta_dex"],
            "map_amplitude_scale_ratio": row["map_amplitude_scale_ratio"],
            "root_amplitude_scale_ratio": row["root_amplitude_scale_ratio"],
            "failure": row["fast_failure_reason"],
        } for row in rows],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
