"""Compare Fast panel-budget roots with read-only Prüfer Reference maps.

For the two high-blue target cases, evaluate Reference fixed-background maps
on the original and 10%-curvature-refined grids at each corresponding Fast
candidate root and a +1e-3 relative perturbation. The resulting secants are
screening estimates, not certified root enclosures.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

from scipy.interpolate import PchipInterpolator  # noqa: E402

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

NESTED = ROOT / "docs" / "adaptive_nested_reference_midpoints_20261105.json"
BUDGET = ROOT / "docs" / "adaptive_panel_fraction_budget_20261105.json"
OUTPUT = ROOT / "docs" / "adaptive_panel_candidate_root_screen_20261105.json"
TARGETS = (
    "adaptive_embedded_fast_high_blue_target1e5_20261105.json",
    "adaptive_embedded_fast_high_blue_target1e6_20261105.json",
)
RTOL = 1e-9
DELTA = 1e-3


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _dn_gw(freqs: np.ndarray, spectrum: np.ndarray, model: LCDM_SG) -> float:
    integral = float(PchipInterpolator(freqs, spectrum).integrate(freqs[0], freqs[-1]))
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    return float(gp.Neff0 * integral * np.log(10.0) / omega_nu)


def _map(params: dict, dn_eff: float, z_tail: float, freqs: np.ndarray) -> dict:
    model = LCDM_SG(**params)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=freqs,
            z_tail=z_tail,
            rtol=RTOL,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    elapsed = time.perf_counter() - started
    if caught:
        raise RuntimeError(f"Reference warnings: {[str(w.message) for w in caught]}")
    returned = np.asarray(result["freqs"], dtype=np.float64)
    order = np.argsort(returned)
    returned = returned[order]
    if not np.array_equal(returned, freqs):
        raise RuntimeError("Reference did not return the exact requested frequency grid")
    spectrum = np.asarray(result["Ogw"], dtype=np.float64)[order] - np.asarray(
        result["Oj"], dtype=np.float64
    )[order]
    if not np.isfinite(spectrum).all() or np.any(spectrum <= 0.0):
        raise FloatingPointError("Reference returned a nonpositive or nonfinite spectrum")
    dn_gw = _dn_gw(freqs, spectrum, model)
    return {
        "DN_eff": dn_eff,
        "DN_gw_reference": dn_gw,
        "root_residual": dn_eff - float(params.get("DN_eff", 0.0)) - dn_gw,
        "runtime_s": elapsed,
        "spectrum": spectrum.tolist(),
    }


def _secant(a: dict, b: dict) -> dict:
    xa, xb = float(a["DN_eff"]), float(b["DN_eff"])
    fa, fb = float(a["root_residual"]), float(b["root_residual"])
    if fb == fa:
        raise ZeroDivisionError("Reference secant denominator is zero")
    root = xb - fb * (xb - xa) / (fb - fa)
    return {
        "secant_root_estimate": root,
        "residual_bracketed": bool(fa == 0.0 or fb == 0.0 or np.signbit(fa) != np.signbit(fb)),
        "secant_slope": (fb - fa) / (xb - xa),
        "fast_root_relative_difference": abs(xa - root) / abs(root),
        "perturbed_fast_root_relative_difference": abs(xb - root) / abs(root),
    }


def _case(target_name: str, nested: dict, budget: dict) -> dict:
    if not any(c["probe_artifact"] == target_name for c in nested["cases"]):
        raise KeyError(f"nested reference artifact has no case {target_name!r}")
    probe_path = ROOT / "docs" / target_name
    probe_data = json.loads(probe_path.read_text(encoding="utf-8"))
    source_path = ROOT / "docs" / probe_data["source_artifact"]
    source = json.loads(source_path.read_text(encoding="utf-8"))
    params = source.get("parameters", source.get("parameters_input"))
    z_tail = float(source["reference"]["z_tail"])
    dn_eff_input = float(params.get("DN_eff", 0.0))
    base_fast = next(c for c in budget["cases"] if c["probe_artifact"] == target_name)
    candidate = base_fast["end_to_end_10pct_midpoint_candidate"]
    x_base = float(candidate["baseline_DN_eff_final"][0])
    x_candidate = float(candidate["candidate_DN_eff_final"][0])
    base_freqs_raw = np.asarray(source["fast"]["frequencies_log10_hz"], dtype=np.float64)
    base_freqs = np.sort(base_freqs_raw)
    mids_raw = np.asarray(probe_data["actual_fast_midpoint_samples"]["frequencies_log10_hz"], dtype=np.float64)
    mids = np.sort(mids_raw)
    base_spectrum = np.asarray(source["fast"]["omega_gw"], dtype=np.float64)[np.argsort(base_freqs_raw)]
    curvature = np.abs(
        np.asarray(PchipInterpolator(base_freqs, np.log10(base_spectrum)).derivative(2)(base_freqs))
    )
    widths = np.diff(base_freqs)
    weights = 0.5 * (base_spectrum[:-1] + base_spectrum[1:]) * widths
    scores = curvature[:-1] * widths * widths / 8.0 * weights
    count = max(1, int(np.ceil(0.10 * mids.size)))
    selected = np.sort(np.argsort(scores)[::-1][:count])
    refined_freqs = np.sort(np.concatenate((base_freqs, mids[selected])))

    maps = {}
    for name, grid, x in (
        ("base", base_freqs, x_base),
        ("candidate", refined_freqs, x_candidate),
    ):
        left = _map(params, x, z_tail, grid)
        right = _map(params, x * (1.0 + DELTA), z_tail, grid)
        maps[name] = {
            "frequency_count": int(grid.size),
            "fast_root": x,
            "fast_root_map": left,
            "perturbed_map": right,
            "secant": _secant(left, right),
        }
    base_root = maps["base"]["secant"]["secant_root_estimate"]
    candidate_root = maps["candidate"]["secant"]["secant_root_estimate"]
    return {
        "probe_artifact": target_name,
        "source_artifact": source_path.name,
        "source_sha256": _sha(source_path),
        "parameters": params,
        "DN_eff_input": dn_eff_input,
        "fixed_point_perturbation_relative": DELTA,
        "reference_rtol": RTOL,
        "selected_fraction": 0.10,
        "selected_midpoint_count": int(count),
        "base": maps["base"],
        "candidate": maps["candidate"],
        "candidate_reference_secant_root_vs_base_reference_secant_root": abs(candidate_root - base_root) / abs(base_root),
        "limitation": "Finite two-map secant estimates on each grid; not a rigorous root enclosure or a full Reference root solve.",
    }


def main() -> None:
    nested = json.loads(NESTED.read_text(encoding="utf-8"))
    budget = json.loads(BUDGET.read_text(encoding="utf-8"))
    cases = [_case(name, nested, budget) for name in TARGETS]
    record = {
        "schema_version": 1,
        "experiment": "adaptive_panel_candidate_reference_root_secant_screen",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "target_artifacts": list(TARGETS),
        "nested_artifact_sha256": _sha(NESTED),
        "budget_artifact_sha256": _sha(BUDGET),
        "runtime_environment": platform.platform(),
        "resources": telemetry(threads=1),
        "cases": cases,
        "script_sha256": _sha(Path(__file__)),
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": OUTPUT.name,
        "cases": [
            {
                "probe": case["probe_artifact"],
                "base": case["base"]["secant"],
                "candidate": case["candidate"]["secant"],
                "reference_root_shift": case["candidate_reference_secant_root_vs_base_reference_secant_root"],
            }
            for case in cases
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
