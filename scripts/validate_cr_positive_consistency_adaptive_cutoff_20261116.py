"""Validate the current default adaptive handoff at a cr>0 cutoff point."""

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

from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_consistency_cutoff_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_consistency_adaptive_cutoff_20261116.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def main() -> int:
    baseline = json.loads(INPUT.read_text(encoding="utf-8"))
    params = {key: float(value) for key, value in baseline["parameters"].items()}
    direct = baseline["direct_reference"]
    if direct["source_sha256"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("direct Reference source hash differs from current Reference")
    reference_spectrum = np.asarray(
        direct["results"]["z12"]["raw_log10_OmegaGW"], dtype=np.float64
    )
    model = LCDM_SG(**params)
    started = time.perf_counter()
    result = model.SGWB_iter(
        engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
        threads=1, accuracy_mode="fast", freq_grid="goal", physical_guard=True,
    )
    runtime = time.perf_counter() - started
    if result is not model:
        raise RuntimeError(f"current default Fast failed: {model.fast_failure_reason}")

    freqs = np.asarray(model.f, dtype=np.float64)
    if freqs.size != baseline["frequency_grid"]["count"]:
        raise RuntimeError("current default Fast grid has a different node count")
    if digest(freqs) != baseline["frequency_grid"]["log10f_sha256_f64le"]:
        raise RuntimeError("current default Fast grid differs from the direct Reference grid")
    fast_spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
    if reference_spectrum.shape != fast_spectrum.shape or not (
        np.isfinite(reference_spectrum).all() and np.isfinite(fast_spectrum).all()
    ):
        raise FloatingPointError("current Fast or direct Reference spectrum is non-finite/incomplete")
    delta = fast_spectrum - reference_spectrum
    worst = int(np.argmax(np.abs(delta)))
    derived = model.derived_param
    output = {
        "schema_version": 1,
        "status": "complete_current_default_cr_positive_consistency_adaptive_handoff",
        "scope": "Current production default Fast at a physical-amplitude cr>0 point just above the analytic cutoff. Compare its self-consistent result and full spectrum with direct same-grid Reference z12; the original baseline artifact contains an earlier pre-adaptive z5 Fast result.",
        "baseline_reference_input": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "parameters": params,
        "tensor_tilt": {
            "input_n_t": float(params["n_t"]),
            "derived_nt_used_by_solver": float(derived["nt"]),
            "input_n_t_controls_spectrum": False,
        },
        "frequency_grid": {
            "count": int(freqs.size),
            "log10f_sha256_f64le": digest(freqs),
            "log10f": freqs.tolist(),
        },
        "current_fast": {
            "self_consistent_returned": True,
            "root_DN_eff": float(model.cosmo_param["DN_eff"]),
            "mapped_DN_gw": float(model.DN_gw[-1]),
            "fixed_point_residual_relative": float(model.DN_gw[-1] / model.cosmo_param["DN_eff"] - 1.0),
            "iteration_delta": float(model.dn_converged_delta),
            "runtime_s": runtime,
            "z_tail_used": float(model.z_tail_used),
            "adaptive_handoff_selected_count": int(getattr(model, "_fast_adaptive_handoff_selected", 0)),
            "adaptive_handoff_selected_fraction": float(getattr(model, "_fast_adaptive_handoff_fraction", 0.0)),
            "A_t": float(derived["A_t"]),
            "A_t_normal": bool(derived["A_t"] >= np.finfo(float).tiny),
            "N_inf": float(derived["N_inf"]),
            "N_re": float(derived["N_re"]),
            "handoff_eps_at_prior_worst_node": float(model.handoff_eps[68]),
        },
        "direct_reference": {
            "source_sha256_matches_current": True,
            "z12_fixed_DN_eff": float(direct["fixed_DN_eff"]),
            "z12_DN_gw": float(direct["results"]["z12"]["DN_gw"]),
            "z10_to_z12_map_relative_change": float(direct["z10_to_z12_map_relative_change"]),
            "z10_to_z12_max_abs_dex": float(direct["z10_to_z12_max_abs_dex_on_positive_nodes"]),
            "warnings": direct["warnings"],
            "public_log_floor_nodes": int(direct["results"]["z12"]["public_log_floor_nodes"]),
        },
        "comparison": {
            "root_relative_error_vs_reference_z12": float(
                model.cosmo_param["DN_eff"] / direct["results"]["z12"]["DN_gw"] - 1.0
            ),
            "max_abs_spectrum_error_vs_reference_z12_dex": float(np.max(np.abs(delta))),
            "p95_abs_spectrum_error_vs_reference_z12_dex": float(np.percentile(np.abs(delta), 95.0)),
            "worst_log10f": float(freqs[worst]),
            "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3),
            "fast_log10_OmegaGW": fast_spectrum.tolist(),
            "reference_z12_raw_log10_OmegaGW": reference_spectrum.tolist(),
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
                "stiff_SGWB.py": "stiffgwpy_fast/stiff_SGWB.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "baseline_input": "docs/research/results/cr_positive_consistency_cutoff_20261116.json",
                "diagnostic_script": "scripts/validate_cr_positive_consistency_adaptive_cutoff_20261116.py",
            }.items()
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "root_relative_error_vs_reference_z12": output["comparison"]["root_relative_error_vs_reference_z12"],
        "max_spectrum_error_dex": output["comparison"]["max_abs_spectrum_error_vs_reference_z12_dex"],
        "p95_spectrum_error_dex": output["comparison"]["p95_abs_spectrum_error_vs_reference_z12_dex"],
        "handoff_eps_at_prior_worst_node": output["current_fast"]["handoff_eps_at_prior_worst_node"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
