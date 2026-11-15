"""Test integral-first scaling on the independent blue-tilt minimum-r holdout."""

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
    "n_t": 17.38,
    "cr": 0.0,
    "T_re": 2000.0,
    "DN_re": 10.0,
    "kappa10": 1e-2,
}
GRID_INPUT = ROOT / "docs" / "research" / "results" / "underflow_grid_235_469_937_20261115.json"
REFERENCE_INPUT = ROOT / "docs" / "research" / "results" / "underflow_min_subnormal_reference_z10_20261115.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "integral_first_blue_tilt_holdout_20261116.json"
COUNTS = (235, 469, 937)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def nested_midpoints(descending: np.ndarray) -> np.ndarray:
    ascending = np.asarray(descending, dtype=np.float64)[::-1]
    refined = np.empty(2 * ascending.size - 1, dtype=np.float64)
    refined[::2] = ascending
    refined[1::2] = 0.5 * (ascending[:-1] + ascending[1:])
    return refined[::-1]


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
    mapped = None if not hasattr(model, "DN_gw") else float(model.DN_gw[-1])
    spectrum = getattr(model, "log10OmegaGW", None)
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
        "fast_failure_reason": model.fast_failure_reason,
        "fast_failure_detail": getattr(model, "fast_failure_detail", None),
        "runtime_s": time.perf_counter() - started,
        "log10_OmegaGW": None if spectrum is None else np.asarray(spectrum, dtype=np.float64).tolist(),
    }


def main() -> int:
    grid_data = json.loads(GRID_INPUT.read_text(encoding="utf-8"))
    ref_data = json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))
    ref_sha = ref_data["source_sha256"]["reference.py"]
    if ref_sha != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("Reference source changed since the direct z10 artifact")
    if ref_data["parameters"] != BASE:
        raise RuntimeError("Direct Reference z10 parameters do not match this holdout")
    if float(ref_data["frozen_DN_eff"]) <= 0.0:
        raise RuntimeError("Direct Reference z10 is missing its frozen positive root")

    grids = {235: np.asarray(grid_data["frequency_grids"]["235"]["log10_hz"], dtype=np.float64)}
    grids[469] = nested_midpoints(grids[235])
    grids[937] = nested_midpoints(grids[469])
    if tuple(grids[n].size for n in COUNTS) != COUNTS:
        raise AssertionError("the recorded nested 235/469/937 grid sizes changed")
    for count in COUNTS:
        expected = grid_data["frequency_grids"][str(count)]["sha256_f64le"]
        if digest(grids[count]) != expected:
            raise RuntimeError(f"the reconstructed {count}-node direct grid hash changed")
    if digest(grids[937]) != ref_data["frequency_sha256_f64le"]:
        raise RuntimeError("the direct 937-node Reference z10 grid does not match")
    ref_spectrum = np.asarray(ref_data["reference"]["log10_OmegaGW"], dtype=np.float64)
    ref_map = float(ref_data["reference"]["mapped_DN_gw"])
    ref8_spectrum = np.asarray(ref_data["reference_z8"]["log10_OmegaGW"], dtype=np.float64)
    reference_root = float(ref_data["frozen_DN_eff"])
    ref_tail_max = float(ref_data["tail_comparison"]["reference_z8_to_z10_max_abs_dex"])
    ref_tail_p95 = float(ref_data["tail_comparison"]["reference_z8_to_z10_p95_abs_dex"])

    current_grid_hashes = {
        str(count): {"count": int(grids[count].size), "sha256_f64le": digest(grids[count])}
        for count in COUNTS
    }
    results = {}
    root235 = None
    for count in COUNTS:
        freq = grids[count]
        self_consistent = run_fast(freq, frozen_dn=None, z_tail=5.0)
        if not self_consistent["returned"]:
            raise RuntimeError(f"current Fast self-consistent solve failed on {count} nodes")
        if count == 235:
            root235 = float(self_consistent["root_DN_eff"])
        root = float(self_consistent["root_DN_eff"])
        fixed_at_root = run_fast(freq, frozen_dn=root, z_tail=5.0)
        fixed_at_reference_root = run_fast(freq, frozen_dn=reference_root, z_tail=5.0)
        if not fixed_at_root["returned"] or not fixed_at_reference_root["returned"]:
            raise RuntimeError(f"current Fast fixed-root map failed on {count} nodes")
        results[str(count)] = {
            "self_consistent": self_consistent,
            "fixed_at_own_root": fixed_at_root,
            "fixed_at_reference_root": fixed_at_reference_root,
        }
    assert root235 is not None

    fast_z8 = run_fast(grids[937], frozen_dn=reference_root, z_tail=8.0)
    if not fast_z8["returned"]:
        raise RuntimeError("current Fast z8 spectrum comparison failed")
    fast_spectrum = np.asarray(fast_z8["log10_OmegaGW"], dtype=np.float64)
    spectrum_delta = fast_spectrum - ref_spectrum
    spectrum_idx = int(np.argmax(np.abs(spectrum_delta)))
    fixed_map = float(results["937"]["fixed_at_reference_root"]["mapped_DN_gw"])
    map_error = fixed_map / ref_map - 1.0
    root_changes = {
        count: float(results[count]["self_consistent"]["root_DN_eff"]) / root235 - 1.0
        for count in ("469", "937")
    }
    fixed_map_changes = {
        count: float(results[count]["fixed_at_reference_root"]["mapped_DN_gw"])
        / float(results["235"]["fixed_at_reference_root"]["mapped_DN_gw"]) - 1.0
        for count in ("469", "937")
    }
    output = {
        "schema_version": 1,
        "status": "complete_current_integral_first_blue_tilt_holdout",
        "scope": "Current Fast direct self-consistent roots on nested 235/469/937 grids for the independent minimum-subnormal blue-tilt point; 937-node z8 spectrum/map compared on the same frozen 235-root and direct grid against the existing source-hash-verified Reference z10 artifact.",
        "parameters": BASE,
        "fast_options": {"tol": 1e-7, "threads": 1, "physical_guard": True},
        "reference_artifact": {
            "path": str(REFERENCE_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "z_tail": 10.0,
            "rtol": ref_data["reference"]["rtol"],
            "frozen_DN_eff": float(ref_data["frozen_DN_eff"]),
            "mapped_DN_gw": ref_map,
            "z8_to_z10_max_abs_dex": ref_tail_max,
            "z8_to_z10_p95_abs_dex": ref_tail_p95,
            "reference_source_sha256_matches_current": ref_sha == sha(ROOT / "stiffgwpy_fast" / "reference.py"),
        },
        "frequency_grids": current_grid_hashes,
        "current_fast_roots_and_fixed_maps": results,
        "holdout_metrics": {
            "root_change_469_vs_235": root_changes["469"],
            "root_change_937_vs_235": root_changes["937"],
            "fixed_map_change_469_vs_235": fixed_map_changes["469"],
            "fixed_map_change_937_vs_235": fixed_map_changes["937"],
            "937_map_at_reference_root_relative_error_vs_direct_reference_z10": map_error,
            "reference_z10_map_residual_relative_to_235_root": ref_map / reference_root - 1.0,
            "937_fast_z8_vs_reference_z10_max_abs_dex": float(np.max(np.abs(spectrum_delta))),
            "937_fast_z8_vs_reference_z10_p95_abs_dex": float(np.percentile(np.abs(spectrum_delta), 95.0)),
            "937_fast_z8_vs_reference_z10_worst_log10f": float(grids[937][spectrum_idx]),
            "937_fast_z8_vs_reference_z10_passes_1e-3_dex": bool(np.max(np.abs(spectrum_delta)) <= 1e-3),
            "reference_z10_interpolation_error_relative": ref_data["reference"]["interpolation_error_relative"],
            "reference_z10_quadrature_error_relative": ref_data["reference"]["quadrature_error_relative"],
        },
        "spectra": {
            "fast_z8_log10_OmegaGW_937": fast_spectrum.tolist(),
            "reference_z10_log10_OmegaGW_937": ref_spectrum.tolist(),
            "reference_z8_log10_OmegaGW_937": ref8_spectrum.tolist(),
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
                "canonical_grid_input": "docs/research/results/underflow_grid_235_469_937_20261115.json",
                "reference_z10_input": "docs/research/results/underflow_min_subnormal_reference_z10_20261115.json",
                "script": "scripts/screen_integral_first_blue_tilt_holdout_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "roots": [(n, results[str(n)]["self_consistent"]["root_DN_eff"]) for n in COUNTS],
        "map_changes": fixed_map_changes,
        "ref_map_error": map_error,
        "spectrum_max": output["holdout_metrics"]["937_fast_z8_vs_reference_z10_max_abs_dex"],
        "spectrum_p95": output["holdout_metrics"]["937_fast_z8_vs_reference_z10_p95_abs_dex"],
        "spectrum_worst_log10f": output["holdout_metrics"]["937_fast_z8_vs_reference_z10_worst_log10f"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
