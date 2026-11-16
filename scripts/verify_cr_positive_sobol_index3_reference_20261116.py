"""Add a direct Reference z8/z10/z12 comparison for Sobol point 3."""

from __future__ import annotationsimport hashlibimport jsonimport mathimport platformimport sysimport timeimport warningsfrom pathlib import Pathimport numpy as npimport scipyROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stiffgwpy_fast.reference import run_reference  # noqa: E402from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402INPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_cutoff_holdout_20261116.json"
OUTPUT = ROOT / "docs" / "research" / "results" / "cr_positive_sobol_index3_reference_20261116.json"
REFERENCE_TAILS = (8.0, 10.0, 12.0)
REFERENCE_RTOL = 1e-9


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def raw_log_spectrum(result: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    signal = np.asarray(result["Ogw"], dtype=np.float64) - np.asarray(result["Oj"], dtype=np.float64)
    valid = np.isfinite(signal) & (signal > 0.0)
    logs = np.full(signal.shape, np.nan, dtype=np.float64)
    logs[valid] = np.log10(signal[valid])
    return logs, valid


def main() -> int:
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    if source["source_sha256"]["reference.py"] != sha(ROOT / "stiffgwpy_fast" / "reference.py"):
        raise RuntimeError("Sobol preflight Reference source SHA differs from current Reference")
    sample = next(row for row in source["samples"] if row["sobol_index"] == 3)
    params = {key: float(value) for key, value in sample["parameters"].items()}
    expected_root = float(sample["fast"]["root_DN_eff"])

    model = LCDM_SG(**params)
    start_fast = time.perf_counter()
    fast_warnings: list[warnings.WarningMessage]
    with warnings.catch_warnings(record=True) as fast_warnings:
        warnings.simplefilter("always", RuntimeWarning)
        result = model.SGWB_iter(
            engine="fast", fallback=False, tol=1e-7, z_tail=5.0,
            threads=1, accuracy_mode="fast", freq_grid="goal",
            physical_guard=True,
        )
    fast_runtime = time.perf_counter() - start_fast
    if result is not model:
        raise RuntimeError(f"Sobol point 3 current Fast failed: {model.fast_failure_reason}")
    root = float(model.cosmo_param["DN_eff"])
    freqs = np.asarray(model.f, dtype=np.float64)
    fast_logs = np.asarray(model.log10OmegaGW, dtype=np.float64)
    grid_hash = digest(freqs)
    if grid_hash != sample["fast"]["frequency_sha256_f64le"]:
        raise RuntimeError("Sobol point 3 current Fast grid differs from its preflight grid")
    if not math.isclose(root, expected_root, rel_tol=2e-12, abs_tol=0.0):
        raise RuntimeError("Sobol point 3 current Fast root differs from its preflight root")

    fast_record = {
        "runtime_s": fast_runtime,
        "root_DN_eff": root,
        "preflight_root_DN_eff": expected_root,
        "root_relative_change_vs_preflight": root / expected_root - 1.0,
        "mapped_DN_gw": float(model.DN_gw[-1]),
        "fixed_point_residual_relative": float(model.DN_gw[-1] / root - 1.0),
        "frequency_count": int(freqs.size),
        "frequency_sha256_f64le": grid_hash,
        "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in fast_warnings],
    }
    script_hash = sha(Path(__file__))
    partial = {
        "schema_version": 1,
        "status": "running_sobol_index3_direct_reference",
        "scope": "Follow-up Sobol index 3 from the fixed 8-point cr>0 cutoff design. Rerun current Fast, verify exact root/grid against preflight, then directly solve Reference z8/z10/z12 at that same fixed root and exact frequency grid. Reference source remains unchanged.",
        "sobol_preflight_input": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "sobol_index": 3,
        "parameters": params,
        "current_fast": fast_record,
        "reference": {
            "source_sha256": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fixed_DN_eff": root,
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": grid_hash,
            "rtol": REFERENCE_RTOL,
            "z_tails": {},
        },
        "source_sha256": {
            "fast_sgwb.py": sha(ROOT / "stiffgwpy_fast" / "fast_sgwb.py"),
            "reference.py": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "preflight_input": sha(INPUT),
            "script": script_hash,
        },
    }
    OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    ref_model = LCDM_SG(**{**params, "DN_eff": root})
    refs: dict[str, object] = {}
    for z_tail in REFERENCE_TAILS:
        start = time.perf_counter()
        caught: list[warnings.WarningMessage]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ref = run_reference(
                ref_model, dn_eff=root, freq_subset=freqs,
                z_tail=z_tail, rtol=REFERENCE_RTOL,
                self_consistent=False, physical_guard=True,
            )
        logs, valid = raw_log_spectrum(ref)
        refs[str(z_tail)] = {
            "DN_gw": float(ref["DN_gw"]),
            "quadrature_error": float(ref["quadrature_error"]),
            "interpolation_error": float(ref["interpolation_error"]),
            "used_tail_fraction": float(np.mean(np.asarray(ref["used_tail"], dtype=bool))),
            "runtime_s": time.perf_counter() - start,
            "raw_signal_positive_finite_nodes": int(np.count_nonzero(valid)),
            "raw_signal_invalid_nodes": int(valid.size - np.count_nonzero(valid)),
            "public_log_floor_nodes": int(np.count_nonzero(np.asarray(ref["log10OmegaGW"]) == -40.0)),
            "runtime_warnings": [f"{item.category.__name__}: {item.message}" for item in caught],
            "raw_log10_OmegaGW": [float(value) if okay else None for value, okay in zip(logs, valid)],
        }
        partial["reference"]["z_tails"] = refs
        OUTPUT.write_text(json.dumps(partial, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    ref10 = np.asarray([np.nan if v is None else float(v) for v in refs["10.0"]["raw_log10_OmegaGW"]])
    ref12 = np.asarray([np.nan if v is None else float(v) for v in refs["12.0"]["raw_log10_OmegaGW"]])
    tail_valid = np.isfinite(ref10) & np.isfinite(ref12)
    tail_delta = ref12[tail_valid] - ref10[tail_valid]
    compare_valid = np.isfinite(fast_logs) & np.isfinite(ref12)
    delta = fast_logs[compare_valid] - ref12[compare_valid]
    ref_map_z8 = float(refs["8.0"]["DN_gw"])
    ref_map_z10 = float(refs["10.0"]["DN_gw"])
    ref_map_z12 = float(refs["12.0"]["DN_gw"])
    output = {
        "schema_version": 1,
        "status": "complete_sobol_index3_direct_reference_holdout",
        "scope": "Follow-up Sobol index 3 from the fixed 8-point cr>0 cutoff design. Rerun current Fast, verify exact root/grid against preflight, then directly solve Reference z8/z10/z12 at that same fixed root and exact frequency grid. Reference source remains unchanged.",
        "sobol_preflight_input": str(INPUT.relative_to(ROOT)).replace("\\", "/"),
        "sobol_index": 3,
        "unit_point": sample["unit_point"],
        "parameters": params,
        "supplied_n_t": params["n_t"],
        "derived_n_t": float(model.derived_param["nt"]),
        "cutoff": sample["cutoff"],
        "stored_A_t": float(model.derived_param["A_t"]),
        "current_fast": {
            **fast_record,
        },
        "reference": {
            "source_sha256": sha(ROOT / "stiffgwpy_fast" / "reference.py"),
            "fixed_DN_eff": root,
            "frequency_count": int(freqs.size),
            "frequency_sha256_f64le": grid_hash,
            "rtol": REFERENCE_RTOL,
            "z_tails": refs,
            "z8_to_z10_map_relative_change": ref_map_z10 / ref_map_z8 - 1.0,
            "z10_to_z12_map_relative_change": ref_map_z12 / ref_map_z10 - 1.0,
            "z10_to_z12_max_abs_dex": float(np.max(np.abs(tail_delta))) if tail_delta.size else None,
            "z10_to_z12_compared_nodes": int(tail_delta.size),
        },
        "comparison": {
            "fast_root_relative_error_vs_reference_z12_map": root / ref_map_z12 - 1.0,
            "fast_default_adaptive_vs_reference_z12": {
                "compared_nodes": int(delta.size),
                "uncompared_nodes": int(freqs.size - delta.size),
                "max_abs_dex": float(np.max(np.abs(delta))) if delta.size else None,
                "p95_abs_dex": float(np.percentile(np.abs(delta), 95.0)) if delta.size else None,
                "passes_1e-3_dex": bool(np.max(np.abs(delta)) <= 1e-3) if delta.size else False,
            },
        },
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
            "preflight_input": sha(INPUT),
            "script": script_hash,
        },
        "spectra": {
            "fast_default_adaptive_log10_OmegaGW": fast_logs.tolist(),
        },
    }
    OUTPUT.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "parameters": {"T_re": params["T_re"], "kappa10": params["kappa10"], "r": params["r"]},
        "root_relative_error_vs_reference_z12": output["comparison"]["fast_root_relative_error_vs_reference_z12_map"],
        "spectrum_max_abs_dex": output["comparison"]["fast_default_adaptive_vs_reference_z12"]["max_abs_dex"],
        "z10_to_z12_root_change": output["reference"]["z10_to_z12_map_relative_change"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

