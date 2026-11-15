"""Test integral-first scaling on the low-temperature blue-tilt holdout."""

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

from stiffgwpy_fast import freq_adaptive as FA  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 5e-324,
    "n_t": 17.45,
    "cr": 0.0,
    "T_re": 10.0,
    "DN_re": 10.0,
    "kappa10": 1e-3,
}
HOLDOUT_INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_holdout_20261115.json"
REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "adaptive_handoff_background_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "integral_first_low_temperature_holdout_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(freqs: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(freqs, dtype="<f8").tobytes()).hexdigest()


def run_fast(freqs: np.ndarray, *, frozen_dn: float | None, z_tail: float) -> dict[str, object]:
    model = LCDM_SG(**BASE)
    if frozen_dn is not None:
        model.cosmo_param["DN_eff"] = float(frozen_dn)
        model.gen_expansion()
    original_gen = model.gen_expansion
    original_goal = FA.goal_oriented_freqs
    original_construct = FA.construct_fast_frequency_grid

    def fixed_background_gen(*args, **kwargs):
        model.cosmo_param["DN_eff"] = float(frozen_dn)
        return original_gen(*args, **kwargs)

    FA.goal_oriented_freqs = lambda *args, **kwargs: freqs.copy()
    FA.construct_fast_frequency_grid = lambda *args, **kwargs: freqs.copy()
    if frozen_dn is not None:
        model.gen_expansion = fixed_background_gen
    started = time.perf_counter()
    try:
        result = model.SGWB_iter(
            engine="fast", fallback=False,
            tol=0.99 if frozen_dn is not None else 1e-7,
            z_tail=z_tail, threads=1, accuracy_mode="fast",
            freq_grid="goal", physical_guard=True,
        )
    finally:
        FA.goal_oriented_freqs = original_goal
        FA.construct_fast_frequency_grid = original_construct
    spectrum = getattr(model, "log10OmegaGW", None)
    mapped = None if not hasattr(model, "DN_gw") else float(model.DN_gw[-1])
    return {
        "returned": result is model,
        "solve_mode": "self_consistent" if frozen_dn is None else "fixed_background",
        "frozen_DN_eff": frozen_dn,
        "root_DN_eff": float(model.cosmo_param["DN_eff"]) if result is model and frozen_dn is None else None,
        "DN_eff_after_call": float(model.cosmo_param["DN_eff"]),
        "mapped_DN_gw": mapped,
        "fixed_point_residual_relative": (
            None if mapped is None or frozen_dn is None or frozen_dn == 0.0
            else mapped / frozen_dn - 1.0
        ),
        "failure_reason": model.fast_failure_reason,
        "failure_detail": getattr(model, "fast_failure_detail", None),
        "runtime_s": time.perf_counter() - started,
        "log10_OmegaGW": None if spectrum is None else np.asarray(spectrum, dtype=np.float64).tolist(),
    }


def main() -> int:
    holdout = json.loads(HOLDOUT_INPUT.read_text(encoding="utf-8"))
    reference = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    if holdout["parameters"] != BASE or reference["parameters"] != BASE:
        raise RuntimeError("holdout/reference parameters do not match this sample")
    if reference["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("Reference source has changed since its direct z10 computation")
    freqs = np.asarray(holdout["frequency_log10_hz"], dtype=np.float64)
    grid_sha = digest(freqs)
    if freqs.size != 221 or grid_sha != holdout["frequency_sha256_f64le"]:
        raise RuntimeError("saved direct 221-node grid hash does not match")
    if grid_sha != reference["frequency_sha256_f64le"]:
        raise RuntimeError("direct Reference z10 was not computed on this exact grid")
    frozen_dn = float(reference["frozen_DN_eff"])
    ref_z8 = np.asarray(reference["reference"]["8.0"]["log10_OmegaGW"], dtype=np.float64)
    ref_z10 = np.asarray(reference["reference"]["10.0"]["log10_OmegaGW"], dtype=np.float64)
    if ref_z8.shape != freqs.shape or ref_z10.shape != freqs.shape:
        raise RuntimeError("direct Reference tail spectra do not match the saved native grid")

    root_run = run_fast(freqs, frozen_dn=None, z_tail=5.0)
    if not root_run["returned"]:
        raise RuntimeError("current Fast self-consistent low-temperature holdout failed")
    current_root = float(root_run["root_DN_eff"])
    fixed_at_current_root = run_fast(freqs, frozen_dn=current_root, z_tail=5.0)
    fixed_at_reference_root = run_fast(freqs, frozen_dn=frozen_dn, z_tail=5.0)
    fast_z8 = run_fast(freqs, frozen_dn=frozen_dn, z_tail=8.0)
    for label, row in (
        ("fixed_at_current_root", fixed_at_current_root),
        ("fixed_at_reference_root", fixed_at_reference_root),
        ("fast_z8", fast_z8),
    ):
        if not row["returned"]:
            raise RuntimeError(f"current Fast {label} run failed")

    fast_spectrum = np.asarray(fast_z8["log10_OmegaGW"], dtype=np.float64)
    delta_z10 = fast_spectrum - ref_z10
    delta_z8 = fast_spectrum - ref_z8
    worst = int(np.argmax(np.abs(delta_z10)))
    output = {
        "schema_version": 1,
        "status": "complete_current_integral_first_low_temperature_holdout",
        "scope": "Current Fast direct self-consistent/fixed-root maps and z8 spectrum for the independent n_t=17.45,T_re=10,kappa10=.001,r=minsub holdout; compare same exact 221-node frozen background to direct Reference z8/z10 spectra.",
        "parameters": BASE,
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_sha,
        "reference_artifacts": {
            "holdout_grid": str(HOLDOUT_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "direct_z10": str(REFERENCE_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "source_sha256_reference_matches_current": True,
            "frozen_DN_eff": frozen_dn,
            "z8_to_z10_max_abs_dex": float(np.max(np.abs(ref_z8 - ref_z10))),
            "z8_to_z10_p95_abs_dex": float(np.percentile(np.abs(ref_z8 - ref_z10), 95.0)),
            "map_z8_to_z10_relative_change": (
                float(reference["reference"]["8.0"]["mapped_DN_gw"])
                / float(reference["reference"]["10.0"]["mapped_DN_gw"]) - 1.0
            ),
            "z10_interpolation_error_relative": reference["reference"]["10.0"]["interpolation_error_relative"],
            "z10_quadrature_error_relative": reference["reference"]["10.0"]["quadrature_error_relative"],
            "z10_warnings": reference["reference"]["10.0"]["warnings"],
            "z10_normalized_nodes_underflowed_to_zero": reference["reference"]["10.0"]["normalized_nodes_underflowed_to_zero"],
        },
        "current_fast": {
            "z5_self_consistent": root_run,
            "z5_fixed_at_current_root": fixed_at_current_root,
            "z5_fixed_at_reference_root": fixed_at_reference_root,
            "z8_fixed_at_reference_root": fast_z8,
        },
        "holdout_metrics": {
            "current_root_relative_to_reference_frozen_root": current_root / frozen_dn - 1.0,
            "fixed_at_current_root_residual_relative": fixed_at_current_root["mapped_DN_gw"] / current_root - 1.0,
            "fixed_at_reference_root_map_relative_error_vs_reference_z10": (
                fixed_at_reference_root["mapped_DN_gw"]
                / reference["reference"]["10.0"]["mapped_DN_gw"] - 1.0
            ),
            "z8_map_relative_error_vs_reference_z10": (
                fast_z8["mapped_DN_gw"] / reference["reference"]["10.0"]["mapped_DN_gw"] - 1.0
            ),
            "z8_vs_reference_z10_max_abs_dex": float(np.max(np.abs(delta_z10))),
            "z8_vs_reference_z10_p95_abs_dex": float(np.percentile(np.abs(delta_z10), 95.0)),
            "z8_vs_reference_z10_worst_log10f": float(freqs[worst]),
            "z8_vs_reference_z10_passes_1e-3_dex": bool(np.max(np.abs(delta_z10)) <= 1e-3),
            "z8_vs_reference_z8_max_abs_dex": float(np.max(np.abs(delta_z8))),
        },
        "spectra": {
            "fast_z8_log10_OmegaGW": fast_spectrum.tolist(),
            "reference_z8_log10_OmegaGW": ref_z8.tolist(),
            "reference_z10_log10_OmegaGW": ref_z10.tolist(),
        },
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
                "exact_background.py": "stiffgwpy_fast/exact_background.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "holdout_grid_input": "docs/research/results/adaptive_handoff_background_holdout_20261115.json",
                "reference_z10_input": "docs/research/results/adaptive_handoff_background_z10_20261115.json",
                "script": "scripts/screen_integral_first_low_temperature_holdout_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "current_root": current_root,
        "fixed_root_residual": output["holdout_metrics"]["fixed_at_current_root_residual_relative"],
        "map_error": output["holdout_metrics"]["fixed_at_reference_root_map_relative_error_vs_reference_z10"],
        "spectrum_max": output["holdout_metrics"]["z8_vs_reference_z10_max_abs_dex"],
        "spectrum_p95": output["holdout_metrics"]["z8_vs_reference_z10_p95_abs_dex"],
        "worst_log10f": output["holdout_metrics"]["z8_vs_reference_z10_worst_log10f"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
