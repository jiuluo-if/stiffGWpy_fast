"""Preflight a cr>0 consistency-branch sample at the valid cutoff boundary."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.reference import run_reference  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

BOUNDARY_INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_cutoff_amplitude_boundary_20261116.json"
DEFAULT_OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_consistency_cutoff_20261116.json"
RATIO_ABOVE_CUTOFF = 1.1
BASE = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "A_s": 2.100549e-9,
    "n_t": 17.3,
    "cr": 1.0,
    "DN_re": 20.0,
    "DN_eff": 0.0,
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-t", type=float, default=17.3)
    parser.add_argument("--t-re", type=float, default=2000.0)
    parser.add_argument("--kappa10", type=float, default=0.01)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    boundary = json.loads(BOUNDARY_INPUT.read_text(encoding="utf-8"))
    selected = next(
        row for row in boundary["cases"]
        if row["T_re"] == args.t_re and row["kappa10"] == args.kappa10
    )
    r_cutoff = float(selected["cutoff_root_r"])
    params = {**BASE, "n_t": args.n_t, "T_re": args.t_re, "kappa10": args.kappa10,
              "r": r_cutoff * RATIO_ABOVE_CUTOFF}
    model = LCDM_SG(**params)
    derived = model.derived_param
    initial_delta_n = (
        None if derived["N_inf"] is None
        else float(derived["N_inf"] - derived["N_re"])
    )
    stored_at = float(derived["A_t"])
    derived_nt = float(derived["nt"])
    output: dict[str, object] = {
        "schema_version": 1,
        "status": "preflight_only",
        "scope": f"cr>0 consistency-branch cutoff-valid preflight with input n_t={args.n_t}, T_re={args.t_re}, kappa10={args.kappa10}. The actual tilt is derived by the consistency relation; input n_t is recorded but does not control the spectrum. Preserve physical background construction and use Fast physical_guard; Reference runs only after Fast returns a finite admitted result.",
        "parameters": params,
        "tensor_tilt": {
            "input_n_t": float(args.n_t),
            "derived_nt_used_by_solver": derived_nt,
            "input_n_t_controls_spectrum": False,
        },
        "cutoff": {
            "boundary_input": str(BOUNDARY_INPUT.relative_to(ROOT)).replace("\\", "/"),
            "r_cutoff": r_cutoff,
            "r_over_r_cutoff": params["r"] / r_cutoff,
            "initial_Delta_N": initial_delta_n,
            "initial_N_inf": None if derived["N_inf"] is None else float(derived["N_inf"]),
            "initial_N_re": float(derived["N_re"]),
            "derived_nt": derived_nt,
            "stored_A_t": stored_at,
            "A_t_class": "rounded_zero" if stored_at == 0.0 else "normal" if stored_at >= np.finfo(float).tiny else "subnormal",
            "log_A_t_finite": math.isfinite(math.log(stored_at)) if stored_at > 0.0 else False,
        },
        "current_fast": {},
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
                "LCDM_stiff_Neff.py": "stiffgwpy_fast/LCDM_stiff_Neff.py",
                "reference.py": "stiffgwpy_fast/reference.py",
                "boundary_input": "docs/research/results/cr_positive_cutoff_amplitude_boundary_20261116.json",
                "script": "scripts/screen_cr_positive_consistency_cutoff_20261116.py",
            }.items()
        },
    }

    if initial_delta_n is None or initial_delta_n <= 0.0:
        output["status"] = "invalid_cutoff_preflight"
        output["current_fast"] = {"called": False, "reason": "initial_cutoff_invalid"}
    elif stored_at <= 0.0 or not math.isfinite(stored_at):
        output["status"] = "invalid_tensor_amplitude_preflight"
        output["current_fast"] = {"called": False, "reason": "stored_A_t_not_positive_finite"}
    else:
        started = time.perf_counter()
        try:
            result = model.SGWB_iter(
                engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
                threads=1, accuracy_mode="fast", freq_grid="goal",
                physical_guard=True,
            )
            runtime = time.perf_counter() - started
            succeeded = result is model
            freqs = np.asarray(getattr(model, "f", []), dtype=np.float64)
            output["status"] = "fast_admitted_reference_pending" if succeeded else "fast_guard_or_solver_rejected"
            output["current_fast"] = {
                "called": True,
                "returned_model": succeeded,
                "failure_reason": model.fast_failure_reason,
                "failure_detail": getattr(model, "fast_failure_detail", None),
                "runtime_s": runtime,
                "root_DN_eff": float(model.cosmo_param["DN_eff"]) if succeeded else None,
                "derived_nt_used_by_solver": float(model.derived_param["nt"]) if succeeded else None,
                "mapped_DN_gw": float(model.DN_gw[-1]) if succeeded and len(model.DN_gw) else None,
                "root_delta": float(model.dn_converged_delta) if succeeded else None,
                "physical_guard": True,
                "frequency_count": int(freqs.size),
                "frequency_sha256_f64le": digest(freqs) if freqs.size else None,
                "public_spectrum_finite": bool(np.all(np.isfinite(model.log10OmegaGW))) if succeeded else None,
                "reference_deferred_until_fast_admitted": not succeeded,
            }
            if succeeded:
                output["status"] = "fast_admitted_reference_required"
                root = float(model.cosmo_param["DN_eff"])
                fast_spectrum = np.asarray(model.log10OmegaGW, dtype=np.float64)
                ref_model = LCDM_SG(**{**params, "DN_eff": root})
                refs: dict[float, dict[str, object]] = {}
                ref_runtime: dict[float, float] = {}
                warnings_seen: list[str] = []
                for z_tail in (8.0, 10.0, 12.0):
                    tick = time.perf_counter()
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always")
                        refs[z_tail] = run_reference(
                            ref_model, dn_eff=root, freq_subset=freqs,
                            z_tail=z_tail, rtol=1e-9,
                            self_consistent=False, physical_guard=True,
                        )
                    ref_runtime[z_tail] = time.perf_counter() - tick
                    warnings_seen.extend(
                        f"z_tail={z_tail}: {w.category.__name__}: {w.message}"
                        for w in caught
                    )

                raw_logs: dict[float, np.ndarray] = {}
                valid_masks: dict[float, np.ndarray] = {}
                ref_summary: dict[str, object] = {}
                for z_tail, ref in refs.items():
                    raw = np.asarray(ref["Ogw"] - ref["Oj"], dtype=np.float64)
                    valid = np.isfinite(raw) & (raw > 0.0)
                    logs = np.full(raw.shape, np.nan, dtype=np.float64)
                    logs[valid] = np.log10(raw[valid])
                    raw_logs[z_tail] = logs
                    valid_masks[z_tail] = valid
                    ref_summary[f"z{int(z_tail)}"] = {
                        "DN_gw": float(ref["DN_gw"]),
                        "quadrature_error": float(ref["quadrature_error"]),
                        "interpolation_error": float(ref["interpolation_error"]),
                        "used_tail_fraction": float(np.mean(np.asarray(ref["used_tail"], dtype=bool))),
                        "runtime_s": ref_runtime[z_tail],
                        "raw_signal_positive_finite_nodes": int(np.count_nonzero(valid)),
                        "raw_signal_invalid_nodes": int(raw.size - np.count_nonzero(valid)),
                        "public_log_floor_nodes": int(np.count_nonzero(np.asarray(ref["log10OmegaGW"]) == -40.0)),
                        "raw_log10_OmegaGW": [float(v) if ok else None for v, ok in zip(logs, valid)],
                    }
                z12_mask = valid_masks[12.0] & np.isfinite(fast_spectrum)
                z10_z12_mask = valid_masks[10.0] & valid_masks[12.0]
                if np.any(z12_mask):
                    delta = fast_spectrum[z12_mask] - raw_logs[12.0][z12_mask]
                    worst_local = int(np.argmax(np.abs(delta)))
                    common_indices = np.flatnonzero(z12_mask)
                    worst_index = int(common_indices[worst_local])
                    spectrum_comparison = {
                        "compared_nodes": int(np.count_nonzero(z12_mask)),
                        "uncompared_nodes": int(freqs.size - np.count_nonzero(z12_mask)),
                        "max_abs_dex": float(np.max(np.abs(delta))),
                        "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)),
                        "worst_log10f": float(freqs[worst_index]),
                        "passes_1e-3_dex_on_compared_nodes": bool(np.max(np.abs(delta)) <= 1e-3),
                    }
                else:
                    spectrum_comparison = {
                        "compared_nodes": 0,
                        "uncompared_nodes": int(freqs.size),
                        "max_abs_dex": None,
                        "p95_abs_dex": None,
                        "worst_log10f": None,
                        "passes_1e-3_dex_on_compared_nodes": False,
                    }
                tail_delta = raw_logs[12.0][z10_z12_mask] - raw_logs[10.0][z10_z12_mask]
                output["status"] = "fast_admitted_reference_complete"
                output["frequency_grid"] = {
                    "count": int(freqs.size),
                    "log10f_sha256_f64le": digest(freqs),
                    "log10f": freqs.tolist(),
                }
                output["direct_reference"] = {
                    "source_sha256": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
                    "fixed_DN_eff": root,
                    "frequency_count": int(freqs.size),
                    "results": ref_summary,
                    "z8_to_z10_map_relative_change": float(refs[10.0]["DN_gw"] / refs[8.0]["DN_gw"] - 1.0),
                    "z10_to_z12_map_relative_change": float(refs[12.0]["DN_gw"] / refs[10.0]["DN_gw"] - 1.0),
                    "z10_to_z12_max_abs_dex_on_positive_nodes": float(np.max(np.abs(tail_delta))) if tail_delta.size else None,
                    "warnings": warnings_seen,
                }
                output["comparison"] = {
                    "fast_spectrum_vs_reference_z12_raw": spectrum_comparison,
                    "fast_map_relative_error_vs_reference_z12": float(model.DN_gw[-1] / refs[12.0]["DN_gw"] - 1.0),
                    "reference_z10_to_z12_map_change": output["direct_reference"]["z10_to_z12_map_relative_change"],
                }
        except Exception as exc:  # Record the failure stage without hiding the exception.
            output["status"] = "fast_exception_reference_not_run"
            output["current_fast"] = {
                "called": True,
                "exception_type": type(exc).__name__,
                "exception": str(exc),
                "runtime_s": time.perf_counter() - started,
                "reference_deferred_until_fast_admitted": True,
            }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "status": output["status"],
        "initial_Delta_N": initial_delta_n,
        "A_t_class": output["cutoff"]["A_t_class"],
        "current_fast": output["current_fast"],
        "comparison": output.get("comparison"),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
