"""Validate adaptive Fast midpoint refinement on a nested Prüfer grid.

This driver solves the exact Fast-grid interval midpoints with the unchanged
Reference/Prüfer backend, combines them with the archived same-grid Reference
map, and compares base, full-midpoint, and curvature-ranked nested PCHIP
integrals at one fixed background. It does not modify or accelerate Reference.
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
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment, telemetry  # noqa: E402

apply_environment()

from stiffgwpy_fast import global_param as gp  # noqa: E402
from stiffgwpy_fast.fast_sgwb import integrate_frequency_pchip, ln10  # noqa: E402
from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES = (
    "adaptive_embedded_integral_fast_midpoints_20261105.json",
    "adaptive_embedded_fast_high_blue_target1e6_20261105.json",
    "adaptive_embedded_fast_high_blue_target1e5_20261105.json",
    "adaptive_embedded_fast_high_blue_kappa100_20261105.json",
    "adaptive_cr1_prufer_embedded_midpoints_20261105.json",
    "adaptive_embedded_fast_cr1_lowT_20261105.json",
    "adaptive_embedded_fast_cr1_red_extreme_20261105.json",
    "adaptive_embedded_fast_ninf222_20261105.json",
)
RTOL = 1e-9
OUTPUT = ROOT / "docs" / "adaptive_nested_reference_midpoints_20261105.json"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _dn_pchip(freqs: np.ndarray, spectrum: np.ndarray, model: LCDM_SG) -> float:
    integral = integrate_frequency_pchip(freqs, spectrum) * ln10
    omega_nu = gp.Omega_nh2 / model.derived_param["h"] ** 2
    return float(gp.Neff0 * integral / omega_nu)


def _case(probe_path: Path) -> dict:
    probe = json.loads(probe_path.read_text(encoding="utf-8"))
    source_path = ROOT / "docs" / probe["source_artifact"]
    source = json.loads(source_path.read_text(encoding="utf-8"))
    params = source.get("parameters", source.get("parameters_input"))
    if not isinstance(params, dict):
        raise RuntimeError(f"missing parameter mapping in {source_path.name}")
    dn_eff = float(source.get("fixed_DN_eff", source["fast"]["DN_eff_final"]))
    z_tail = float(source["reference"]["z_tail"])

    base_freqs = np.asarray(source["fast"]["frequencies_log10_hz"], dtype=np.float64)
    base_fast = np.asarray(source["fast"]["omega_gw"], dtype=np.float64)
    base_ref_key = probe.get("adaptive_reference_key", "adaptive_exact_fast_grid")
    base_ref_record = source["reference_maps"][base_ref_key]
    ref_base_freqs = np.asarray(base_ref_record["frequencies_log10_hz"], dtype=np.float64)
    ref_base = np.asarray(base_ref_record["omega_gw"], dtype=np.float64)
    if not np.array_equal(base_freqs, ref_base_freqs):
        raise RuntimeError("archived Fast and same-grid Reference frequencies differ")

    mids = np.asarray(
        probe["actual_fast_midpoint_samples"]["frequencies_log10_hz"], dtype=np.float64
    )
    fast_mid = np.asarray(probe["actual_fast_midpoint_samples"]["omega_gw"], dtype=np.float64)
    if mids.size != base_freqs.size - 1 or fast_mid.size != mids.size:
        raise RuntimeError("Fast probe does not provide one actual midpoint per base panel")

    model = LCDM_SG(**params)
    started = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ref_mid_result = run_reference(
            model,
            dn_eff=dn_eff,
            freq_subset=mids,
            z_tail=z_tail,
            rtol=RTOL,
            self_consistent=False,
            workers=1,
            physical_guard=False,
            backend="prufer",
        )
    ref_mid_runtime = time.perf_counter() - started
    if caught:
        raise RuntimeError(f"Reference warnings: {[str(item.message) for item in caught]}")
    ref_mid_freqs = np.asarray(ref_mid_result["freqs"], dtype=np.float64)
    ref_mid = np.asarray(ref_mid_result["Ogw"], dtype=np.float64) - np.asarray(
        ref_mid_result["Oj"], dtype=np.float64
    )
    order = np.argsort(ref_mid_freqs)
    ref_mid_freqs = ref_mid_freqs[order]
    ref_mid = ref_mid[order]
    mids_ordered = np.sort(mids)
    fast_mid_ordered = fast_mid[np.argsort(mids)]
    if not np.array_equal(ref_mid_freqs, mids_ordered):
        raise RuntimeError("Reference did not return the requested exact midpoint frequencies")
    if not np.all(np.isfinite(ref_mid)) or not np.all(ref_mid > 0.0):
        raise FloatingPointError("Reference midpoint spectrum is not finite and positive")

    base_x = np.sort(base_freqs)
    base_fast_sorted = base_fast[np.argsort(base_freqs)]
    base_ref_sorted = ref_base[np.argsort(ref_base_freqs)]
    nested_freqs = np.sort(np.concatenate((base_x, mids_ordered)))
    nested_fast = np.empty_like(nested_freqs)
    nested_ref = np.empty_like(nested_freqs)
    nested_fast[::2] = base_fast_sorted
    nested_fast[1::2] = fast_mid_ordered
    nested_ref[::2] = base_ref_sorted
    nested_ref[1::2] = ref_mid

    base_fast_dn = _dn_pchip(base_x, base_fast_sorted, model)
    base_ref_dn = _dn_pchip(base_x, base_ref_sorted, model)
    nested_fast_dn = _dn_pchip(nested_freqs, nested_fast, model)
    nested_ref_dn = _dn_pchip(nested_freqs, nested_ref, model)
    dense_key = probe["dense_reference_key"]
    selected_dense_ref = source["reference_maps"].get(dense_key, source.get(dense_key))
    if not isinstance(selected_dense_ref, dict):
        raise RuntimeError(f"missing archived dense Reference map {dense_key!r}")
    selected_grid_ref_dn = float(selected_dense_ref["DN_gw"])

    # Recompute the production curvature×panel-weight ranking on the base Fast spectrum.
    widths = np.diff(base_x)
    curvature = np.abs(
        np.asarray(PchipInterpolator(base_x, np.log10(base_fast_sorted)).derivative(2)(base_x))
    )
    omega_weights = 0.5 * (base_fast_sorted[:-1] + base_fast_sorted[1:]) * widths
    scores = curvature[:-1] * widths * widths / 8.0 * omega_weights
    selected = np.sort(np.argsort(scores)[::-1][: max(1, int(np.ceil(0.25 * mids.size)))])
    selected_unsorted_freqs = np.concatenate((base_x, mids_ordered[selected]))
    selected_order = np.argsort(selected_unsorted_freqs)
    selected_freqs = selected_unsorted_freqs[selected_order]
    selected_fast = np.concatenate((base_fast_sorted, fast_mid_ordered[selected]))[selected_order]
    selected_ref = np.concatenate((base_ref_sorted, ref_mid[selected]))[selected_order]
    selected_fast_dn = _dn_pchip(selected_freqs, selected_fast, model)
    selected_ref_dn = _dn_pchip(selected_freqs, selected_ref, model)

    fast_vs_ref_base = np.abs(np.log10(base_fast_sorted / base_ref_sorted))
    fast_vs_ref_nested = np.abs(np.log10(nested_fast / nested_ref))
    fast_vs_ref_selected = np.abs(np.log10(selected_fast / selected_ref))
    return {
        "probe_artifact": probe_path.name,
        "probe_sha256": _sha(probe_path),
        "source_artifact": source_path.name,
        "source_sha256": _sha(source_path),
        "base_reference_map": base_ref_key,
        "parameters": params,
        "fixed_DN_eff": dn_eff,
        "z_tail": z_tail,
        "frequency_counts": {
            "base": int(base_x.size),
            "all_exact_midpoints": int(mids_ordered.size),
            "nested": int(nested_freqs.size),
            "curvature_selected_midpoints": int(selected.size),
            "curvature_selected_total": int(selected_freqs.size),
        },
        "reference_midpoint_runtime_s": ref_mid_runtime,
        "pchip_DN_gw": {
            "base_fast": base_fast_dn,
            "base_reference": base_ref_dn,
            "nested_fast": nested_fast_dn,
            "nested_reference": nested_ref_dn,
            "curvature_selected_fast": selected_fast_dn,
            "curvature_selected_reference": selected_ref_dn,
            "archived_reference_grid_DN_gw": selected_grid_ref_dn,
        },
        "relative_changes": {
            "fast_nested_vs_base": abs(nested_fast_dn - base_fast_dn) / max(abs(nested_fast_dn), abs(base_fast_dn)),
            "reference_nested_vs_base": abs(nested_ref_dn - base_ref_dn) / max(abs(nested_ref_dn), abs(base_ref_dn)),
            "fast_base_vs_reference_base": abs(base_fast_dn - base_ref_dn) / max(abs(base_fast_dn), abs(base_ref_dn)),
            "fast_nested_vs_reference_nested": abs(nested_fast_dn - nested_ref_dn) / max(abs(nested_fast_dn), abs(nested_ref_dn)),
            "curvature_selected_fast_vs_reference_selected": abs(selected_fast_dn - selected_ref_dn) / max(abs(selected_fast_dn), abs(selected_ref_dn)),
            "nested_reference_vs_archived_dense_grid": abs(nested_ref_dn - selected_grid_ref_dn) / max(abs(nested_ref_dn), abs(selected_grid_ref_dn)),
        },
        "spectrum_abs_dex": {
            "fast_vs_reference_base_max": float(np.max(fast_vs_ref_base)),
            "fast_vs_reference_base_p95": float(np.quantile(fast_vs_ref_base, 0.95)),
            "fast_vs_reference_nested_max": float(np.max(fast_vs_ref_nested)),
            "fast_vs_reference_nested_p95": float(np.quantile(fast_vs_ref_nested, 0.95)),
            "fast_vs_reference_curvature_selected_max": float(np.max(fast_vs_ref_selected)),
            "fast_vs_reference_curvature_selected_p95": float(np.quantile(fast_vs_ref_selected, 0.95)),
        },
        "nested_reference_spectrum": {
            "frequencies_log10_hz": nested_freqs.tolist(),
            "omega_gw": nested_ref.tolist(),
        },
        "limitation": "One fixed background per case; nested-grid PCHIP changes are convergence evidence, not a rigorous error bound or self-consistent root certificate.",
    }


def main() -> None:
    cases = [_case(ROOT / "docs" / name) for name in CASES]
    record = {
        "schema_version": 1,
        "experiment": "nested_exact_prufer_reference_midpoint_validation",
        "created_at": datetime.now().astimezone().isoformat(),
        "source_commit": _git("rev-parse", "HEAD"),
        "reference_rtol": RTOL,
        "cases": cases,
        "resources": telemetry(threads=1),
        "runtime_environment": platform.platform(),
        "script_sha256": _sha(Path(__file__)),
    }
    OUTPUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": OUTPUT.name,
        "cases": [
            {
                "probe": case["probe_artifact"],
                "reference_midpoint_runtime_s": case["reference_midpoint_runtime_s"],
                "relative_changes": case["relative_changes"],
                "spectrum_abs_dex": case["spectrum_abs_dex"],
            }
            for case in cases
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
