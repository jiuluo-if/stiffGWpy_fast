"""Compare the same frozen Fast cases across two source revisions."""

from __future__ import annotations

import argparse
import contextlib
import cProfile
import hashlib
import importlib.metadata
import io
import json
import math
import os
import platform
import pstats
import queue
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path

BASE_PARAMS = {
    "Omega_bh2": 0.0223828,
    "Omega_ch2": 0.1201075,
    "H0": 67.32117,
    "DN_eff": 0.0,
    "A_s": 2.100549e-9,
    "r": 1.0e-2,
    "n_t": 0.5,
    "cr": 0.0,
    "T_re": 2.0e3,
    "DN_re": 10.0,
    "kappa10": 1.0e-2,
}


def classify_pair(old_status: str, new_status: str) -> str:
    """Classify a pair without conflating numerical failure and input validity."""
    statuses = (old_status, new_status)
    if "resource_exceeded" in statuses:
        return "resource_exceeded"
    if statuses == ("mathematically_invalid", "mathematically_invalid"):
        return "mathematically_invalid"
    if "mathematically_invalid" in statuses:
        return "classification_disagreement"
    if statuses == ("success", "success"):
        return "persistent_success"
    if statuses == ("success", "failed"):
        return "old_success_to_new_failure"
    if statuses == ("failed", "success"):
        return "old_failure_to_new_success"
    if statuses == ("failed", "failed"):
        return "persistent_failure"
    return "classification_disagreement"


def compare_result_pair(old: dict, new: dict) -> dict:
    """Compare two solver records only after enforcing exact input identity."""
    if old.get("parameter_sha256") != new.get("parameter_sha256"):
        raise ValueError("paired results must have identical parameter SHA-256")
    comparison = {
        "parameter_sha256": old["parameter_sha256"],
        "classification": classify_pair(old["status"], new["status"]),
        "DN_gw_relative_difference": None,
        "spectrum_comparable": False,
        "spectrum_max_abs_difference_dex": None,
        "spectrum_comparison_bins": 0,
    }
    old_root = old.get("DN_gw")
    new_root = new.get("DN_gw")
    if old_root is not None and new_root is not None and float(old_root) != 0.0:
        comparison["DN_gw_relative_difference"] = (
            abs(float(new_root) - float(old_root)) / abs(float(old_root))
        )

    old_frequency = old.get("log_frequency")
    new_frequency = new.get("log_frequency")
    old_spectrum = old.get("log10OmegaGW")
    new_spectrum = new.get("log10OmegaGW")
    if (
        isinstance(old_frequency, list)
        and old_frequency == new_frequency
        and isinstance(old_spectrum, list)
        and isinstance(new_spectrum, list)
        and len(old_spectrum) == len(new_spectrum) == len(old_frequency)
    ):
        old_valid = [value is not None for value in old_spectrum]
        new_valid = [value is not None for value in new_spectrum]
        if old_valid == new_valid:
            differences = [
                abs(float(left) - float(right))
                for left, right, valid in zip(old_spectrum, new_spectrum, old_valid)
                if valid
            ]
            comparison["spectrum_comparable"] = True
            comparison["spectrum_comparison_bins"] = len(differences)
            if differences:
                comparison["spectrum_max_abs_difference_dex"] = max(differences)
    return comparison


def build_comparison_config(config_type, threads: int):
    """Build the explicit public fast profile shared by both source workers."""
    return config_type(
        h=0.005,
        col_step=8,
        z_tail=5.0,
        phase_max=0.25,
        freq_grid="goal",
        threads=threads,
        kink_split=True,
        transition_refine_radius=2.0,
        transition_refine_factor=8,
        parallel_chunksize=None,
    )


def environment_snapshot(threads: int) -> dict:
    versions = {}
    for package in ("numpy", "scipy", "numba"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    try:
        import psutil

        affinity = psutil.Process().cpu_affinity()
    except (ImportError, AttributeError, OSError):
        affinity = None
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "cpu_affinity": affinity,
        "threads_argument": threads,
        "library_versions": versions,
        "blas_thread_environment": {
            key: os.environ.get(key)
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }


def profile_callable(function, *args, limit: int = 30, **kwargs):
    """Run one diagnostic call and return bounded cumulative Python profile text."""
    profiler = cProfile.Profile()
    result = profiler.runcall(function, *args, **kwargs)
    output = io.StringIO()
    pstats.Stats(profiler, stream=output).strip_dirs().sort_stats("cumtime").print_stats(limit)
    return result, output.getvalue()


def profile_stage_calls(functions: dict, call):
    """Time selected real function calls during one diagnostic invocation."""
    timings = {name: 0.0 for name in functions}
    originals = []
    for name, (module, attribute) in functions.items():
        original = getattr(module, attribute)

        def timed(*args, _name=name, _original=original, **kwargs):
            started = time.perf_counter()
            try:
                return _original(*args, **kwargs)
            finally:
                timings[_name] += time.perf_counter() - started

        originals.append((module, attribute, original))
        setattr(module, attribute, timed)
    try:
        return call(), timings
    finally:
        for module, attribute, original in originals:
            setattr(module, attribute, original)


def reference_from_artifact(path: Path, params: dict[str, float]) -> dict:
    """Extract a source-hashed, same-input Reference map and its exact support."""
    path = Path(path)
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("params") != params:
        raise ValueError("Reference artifact parameters must exactly match the case input")
    if artifact.get("reference_production_modified") is not False:
        raise ValueError("Reference artifact does not confirm an unchanged Reference solver")
    source_hashes = artifact.get("source_sha256")
    if not isinstance(source_hashes, dict) or not source_hashes:
        raise ValueError("Reference artifact must include source SHA-256 values")
    saved_reference = artifact.get("reference", {})
    saved_fast = artifact.get("fast", {})
    frequencies = saved_fast.get("freqs_log10_hz")
    ogw = saved_reference.get("Ogw")
    oj = saved_reference.get("Oj")
    if not (
        isinstance(frequencies, list)
        and isinstance(ogw, list)
        and isinstance(oj, list)
        and len(frequencies) == len(ogw) == len(oj)
    ):
        raise ValueError("Reference map and saved frequency support are incomplete")
    log_spectrum = []
    for gw, shear in zip(ogw, oj):
        signal = float(gw) - float(shear)
        log_spectrum.append(
            _safe_float(math.log10(signal)) if math.isfinite(signal) and signal > 0.0 else None
        )
    return {
        "root": _safe_float(saved_reference.get("root")),
        "log_frequency": frequencies,
        "log10OmegaGW": log_spectrum,
        "source_sha256": source_hashes,
        "artifact_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "artifact_status": artifact.get("status"),
        "map_count": saved_reference.get("map_count"),
        "rtol": saved_reference.get("rtol"),
        "z_tail": saved_reference.get("z_tail"),
        "root_closed": False,
    }


def compare_reference_result(result: dict, reference: dict) -> dict:
    """Check a result against a saved Reference map without overclaiming closure."""
    comparison = {
        "certification": "anchor_screen_only",
        "precision_status": "not_comparable",
        "reference_root_relative_error": None,
        "reference_root_pass": None,
        "reference_spectrum_max_abs_error_dex": None,
        "reference_spectrum_pass": None,
        "reference_spectrum_comparable": False,
    }
    result_root = result.get("DN_gw")
    reference_root = reference.get("root")
    if result_root is not None and reference_root not in (None, 0.0):
        root_error = abs(float(result_root) - float(reference_root)) / abs(float(reference_root))
        comparison["reference_root_relative_error"] = root_error
        comparison["reference_root_pass"] = root_error <= 2.0e-4

    frequencies = result.get("log_frequency")
    spectrum = result.get("log10OmegaGW")
    ref_frequencies = reference.get("log_frequency")
    ref_spectrum = reference.get("log10OmegaGW")
    if (
        result.get("status") == "success"
        and frequencies == ref_frequencies
        and isinstance(spectrum, list)
        and isinstance(ref_spectrum, list)
        and len(spectrum) == len(ref_spectrum)
    ):
        valid = [value is not None for value in ref_spectrum]
        if all(spectrum[index] is not None for index, use in enumerate(valid) if use):
            differences = [
                abs(float(spectrum[index]) - float(ref_spectrum[index]))
                for index, use in enumerate(valid)
                if use
            ]
            comparison["reference_spectrum_comparable"] = True
            if differences:
                maximum = max(differences)
                comparison["reference_spectrum_max_abs_error_dex"] = maximum
                comparison["reference_spectrum_pass"] = maximum <= 1.0e-3

    checks = (comparison["reference_root_pass"], comparison["reference_spectrum_pass"])
    if all(value is True for value in checks):
        comparison["precision_status"] = "anchor_pass"
    elif any(value is False for value in checks):
        comparison["precision_status"] = "anchor_fail"
    return comparison


def parameter_fingerprint(params: dict[str, float]) -> str:
    """Return a stable digest for exact same-input pairing."""
    encoded = json.dumps(
        params, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def historical_record_to_case(record: dict, source: dict) -> dict:
    """Adapt one archived scan row without upgrading finite output to trust."""
    params = record.get("params")
    if not isinstance(params, dict):
        raise ValueError("historical record must contain a params mapping")
    digest = parameter_fingerprint(params)
    if record.get("parameter_sha256") != digest:
        raise ValueError("historical record parameter SHA-256 mismatch")
    status = record.get("status")
    if status == "success":
        sample_kind = "historical_output_unverified"
    elif status in ("failed", "exception", "returned_nonfinite_output"):
        sample_kind = "historical_failure"
    elif status in ("preflight_reject", "mathematically_invalid"):
        sample_kind = "mathematically_invalid"
    elif status in ("preflight_skip", "resource_budget_interrupted", "resource_exceeded"):
        sample_kind = "resource_case"
    else:
        sample_kind = "historical_unclassified"
    index = int(record["sample_index"])
    cr = float(params.get("cr", 0.0))
    return {
        "label": f"history_{index}_{record.get('label', 'case')}",
        "sample_kind": sample_kind,
        "region_id": f"historical_{source.get('suite', 'scan')}",
        "sample_index": index,
        "seed": source.get("suite_seed"),
        "params": params,
        "parameter_sha256": digest,
        "n_t_effective": (
            -float(params["r"]) / 8.0 if cr > 0.0 else _safe_float(params.get("n_t"))
        ),
        "historical_status": status,
        "historical_classification": record.get("classification"),
        "historical_reason": record.get("reason"),
        "history_source": source,
    }


def parse_historical_jsonl(
    raw: str | bytes,
    *,
    artifact_ref: str,
    artifact_path: str,
    sample_indices: list[int],
) -> list[dict]:
    """Select exact rows from an archived scan and bind their source digest."""
    raw_bytes = raw if isinstance(raw, bytes) else raw.encode("utf-8")
    records = [
        json.loads(line)
        for line in raw_bytes.decode("utf-8-sig").splitlines()
        if line.strip()
    ]
    if not records or records[0].get("record_type") != "manifest":
        raise ValueError("historical JSONL must start with a manifest record")
    manifest = records[0]
    source = {
        "artifact_ref": artifact_ref,
        "artifact_path": artifact_path,
        "artifact_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "scan_source_commit": manifest.get("commit"),
        "suite": manifest.get("suite"),
        "suite_seed": manifest.get("seed"),
    }
    by_index = {}
    for record in records[1:]:
        if "sample_index" in record:
            index = int(record["sample_index"])
            if index in by_index:
                raise ValueError(f"duplicate historical sample index: {index}")
            by_index[index] = record
    if len(sample_indices) != len(set(sample_indices)):
        raise ValueError("historical sample indices must be unique")
    cases = []
    for index in sample_indices:
        record = by_index.get(int(index))
        if record is None:
            raise ValueError(f"historical sample index not found: {index}")
        cases.append(historical_record_to_case(record, source))
    return cases


def load_historical_cases(
    git_ref: str,
    artifact_path: str,
    sample_indices: list[int],
) -> list[dict]:
    """Read selected JSONL records from a committed historical scan artifact."""
    relative_path = Path(artifact_path)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ValueError("historical artifact path must stay within the repository")
    repo_root = Path(__file__).resolve().parents[1]
    commit = _resolve_git_commit(repo_root, git_ref)
    object_name = f"{commit}:{relative_path.as_posix()}"
    raw = subprocess.check_output(["git", "show", object_name], cwd=repo_root)
    return parse_historical_jsonl(
        raw,
        artifact_ref=commit,
        artifact_path=relative_path.as_posix(),
        sample_indices=sample_indices,
    )


def expand_historical_cases(
    cases: list[dict],
    neighborhood_axes: dict,
    *,
    power: int,
    seed: int,
) -> list[dict]:
    """Add declared local Sobol neighborhoods around archived failure rows."""
    expanded = []
    for case in cases:
        if case.get("sample_kind") != "historical_failure" or not neighborhood_axes:
            expanded.append(case)
            continue
        anchor = {
            "label": case["label"],
            "sample_kind": "historical_failure",
            "region_id": case["region_id"],
            "params": case["params"],
            "axes": neighborhood_axes,
        }
        related = generate_anchor_cases(anchor, power=power, seed=seed)
        for row in related:
            row["historical_sample_index"] = case["sample_index"]
            row["historical_status"] = case["historical_status"]
            row["historical_classification"] = case["historical_classification"]
            row["historical_reason"] = case["historical_reason"]
            row["history_source"] = case["history_source"]
        expanded.extend(related)
    return expanded


def _cases_fingerprint(cases: list[dict]) -> str:
    encoded = "\n".join(
        json.dumps(case, sort_keys=True, separators=(",", ":"), allow_nan=False)
        for case in cases
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_case_manifest(path: Path, cases: list[dict], metadata: dict) -> None:
    """Write an immutable-input JSONL manifest with per-case and set hashes."""
    if not isinstance(metadata, dict):
        raise TypeError("metadata must be a mapping")
    labels = set()
    for case in cases:
        if not isinstance(case, dict) or not isinstance(case.get("params"), dict):
            raise ValueError("each case must contain a params mapping")
        if case.get("parameter_sha256") != parameter_fingerprint(case["params"]):
            raise ValueError(f"parameter SHA-256 mismatch for {case.get('label')!r}")
        label = case.get("label")
        if not isinstance(label, str) or not label or label in labels:
            raise ValueError("case labels must be non-empty and unique")
        labels.add(label)

    manifest = {
        "record_type": "manifest",
        "schema_version": 1,
        "case_count": len(cases),
        "cases_sha256": _cases_fingerprint(cases),
        "metadata": metadata,
    }
    lines = [json.dumps(manifest, ensure_ascii=False, sort_keys=True, allow_nan=False)]
    lines.extend(
        json.dumps(case, ensure_ascii=False, sort_keys=True, allow_nan=False)
        for case in cases
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")
    temporary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary_path.replace(path)


def read_case_manifest(path: Path) -> tuple[dict, list[dict]]:
    """Load a frozen manifest and verify every parameter and set digest."""
    records = [
        json.loads(line)
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not records or records[0].get("record_type") != "manifest":
        raise ValueError("case manifest must start with a manifest record")
    manifest = records[0]
    cases = records[1:]
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported case manifest schema")
    if manifest.get("case_count") != len(cases):
        raise ValueError("case count does not match manifest")
    for case in cases:
        if case.get("parameter_sha256") != parameter_fingerprint(case.get("params", {})):
            raise ValueError(f"parameter SHA-256 mismatch for {case.get('label')!r}")
    if manifest.get("cases_sha256") != _cases_fingerprint(cases):
        raise ValueError("case set SHA-256 mismatch")
    return manifest, cases


def generate_region_cases(
    region: dict, *, power: int, seed: int, sample_kind: str = "sobol_interior"
) -> list[dict]:
    """Generate a reproducible Sobol interior sample for one declared region.

    Axes are supplied as ``[low, high]`` for linear sampling or as
    ``{"bounds": [low, high], "scale": "log"}`` for logarithmic sampling.
    A positive-``cr`` region derives ``n_t`` from ``r`` and never samples it.
    """
    from scipy.stats import qmc

    region_id = str(region.get("region_id", "")).strip()
    if not region_id:
        raise ValueError("region_id must be a non-empty string")
    if isinstance(power, bool) or not isinstance(power, int) or not 0 <= power <= 16:
        raise ValueError("power must be an integer in [0, 16]")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if sample_kind not in ("sobol_interior", "sobol_holdout"):
        raise ValueError("sample_kind must be sobol_interior or sobol_holdout")

    fixed = dict(BASE_PARAMS)
    fixed.update(region.get("fixed", {}))
    cr = float(region.get("cr", fixed.get("cr", 0.0)))
    if not math.isfinite(cr):
        raise ValueError("cr must be finite")
    fixed["cr"] = cr

    axes = region.get("axes")
    if not isinstance(axes, dict) or not axes:
        raise ValueError("axes must be a non-empty mapping")
    unknown = set(fixed) - set(BASE_PARAMS)
    unknown.update(set(axes) - set(BASE_PARAMS))
    if unknown:
        raise ValueError(f"unknown parameter axes or fixed values: {sorted(unknown)}")
    if "cr" in axes:
        raise ValueError("cr selects a physical branch and cannot be an axis")
    if cr > 0.0 and "n_t" in axes:
        raise ValueError("n_t is derived from r when cr > 0")

    axis_specs = []
    for key, spec in sorted(axes.items()):
        if isinstance(spec, dict):
            bounds = spec.get("bounds")
            scale = spec.get("scale", "linear")
        else:
            bounds = spec
            scale = "linear"
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError(f"axis {key} must have two bounds")
        low, high = map(float, bounds)
        if not math.isfinite(low) or not math.isfinite(high) or high <= low:
            raise ValueError(f"axis {key} bounds must be finite and increasing")
        if scale not in ("linear", "log"):
            raise ValueError(f"axis {key} scale must be 'linear' or 'log'")
        if scale == "log" and low <= 0.0:
            raise ValueError(f"log axis {key} bounds must be positive")
        axis_specs.append((key, low, high, scale))

    points = qmc.Sobol(
        d=len(axis_specs), scramble=True, seed=seed
    ).random_base2(power)
    cases = []

    def add_case(point, sample_kind: str, label_suffix: str) -> None:
        params = fixed.copy()
        for coordinate, (key, low, high, scale) in enumerate(axis_specs):
            unit = float(point[coordinate])
            if scale == "log":
                params[key] = math.exp(math.log(low) + unit * (math.log(high) - math.log(low)))
            else:
                params[key] = low + unit * (high - low)
        if cr > 0.0:
            if float(params["r"]) <= 0.0:
                raise ValueError("r must be positive when cr > 0")
            params["n_t"] = -float(params["r"]) / 8.0
        cases.append({
            "label": f"{region_id}_{label_suffix}",
            "sample_kind": sample_kind,
            "region_id": region_id,
            "sample_index": len(cases),
            "seed": seed,
            "params": params,
            "parameter_sha256": parameter_fingerprint(params),
        })

    for index, point in enumerate(points):
        add_case(point, sample_kind, f"{sample_kind}_{seed}_{index:04d}")

    center = [0.5] * len(axis_specs)
    for coordinate, (key, _low, _high, _scale) in enumerate(axis_specs):
        for side, bound_name in ((0.0, "low"), (1.0, "high")):
            point = center.copy()
            point[coordinate] = side
            add_case(
                point,
                "region_boundary",
                f"boundary_{key}_{bound_name}",
            )
    return cases


def generate_anchor_cases(anchor: dict, *, power: int, seed: int) -> list[dict]:
    """Keep an audited anchor and generate its declared Sobol neighborhood."""
    label = str(anchor.get("label", "")).strip()
    kind = anchor.get("sample_kind")
    if not label:
        raise ValueError("anchor label must be a non-empty string")
    if kind not in ("historical_failure", "trusted_point"):
        raise ValueError("unsupported anchor sample_kind")
    raw_params = anchor.get("params")
    if not isinstance(raw_params, dict):
        raise ValueError("anchor params must be a mapping")

    params = dict(raw_params)
    effective_n_t = _safe_float(params.get("n_t"))
    if float(params.get("cr", 0.0)) > 0.0:
        r = float(params["r"])
        if not math.isfinite(r) or r <= 0.0:
            raise ValueError("r must be positive and finite when cr > 0")
        effective_n_t = -r / 8.0
    parameter_fingerprint(params)

    axes = {}
    for key, spec in (anchor.get("axes") or {}).items():
        if not isinstance(spec, dict) or "bounds" in spec:
            axes[key] = spec
            continue
        scale = spec.get("scale", "linear")
        center_value = float(params[key])
        if scale == "log" and "factor" in spec:
            factor = float(spec["factor"])
            if center_value <= 0.0 or not math.isfinite(factor) or factor <= 1.0:
                raise ValueError(f"log neighborhood factor for {key} must be > 1 at a positive anchor")
            axes[key] = {
                "bounds": [center_value / factor, center_value * factor],
                "scale": "log",
            }
        elif scale == "linear" and "delta" in spec:
            delta = float(spec["delta"])
            if not math.isfinite(delta) or delta <= 0.0:
                raise ValueError(f"linear neighborhood delta for {key} must be positive")
            axes[key] = {
                "bounds": [center_value - delta, center_value + delta],
                "scale": "linear",
            }
        else:
            raise ValueError(f"axis {key} must specify bounds or a relative factor/delta")

    region_id = str(anchor.get("region_id", label)).strip()
    center = {
        "label": f"{label}_anchor",
        "sample_kind": kind,
        "region_id": region_id,
        "sample_index": 0,
        "seed": seed,
        "input_params": raw_params,
        "params": params,
        "n_t_effective": effective_n_t,
        "parameter_sha256": parameter_fingerprint(params),
    }
    reference_path = anchor.get("reference_artifact")
    if reference_path is not None:
        center["reference"] = reference_from_artifact(reference_path, params)
    region = {
        "region_id": f"{region_id}_{label}",
        "cr": params.get("cr", 0.0),
        "fixed": params,
        "axes": axes,
    }
    neighbors = generate_region_cases(region, power=power, seed=seed)
    for case in neighbors:
        location = case["sample_kind"]
        case["label"] = f"{label}_{case['label']}"
        case["sample_kind"] = f"{kind}_neighborhood"
        case["sample_location"] = location
        case["anchor_label"] = label
        case["input_params"] = dict(case["params"])
        case["sample_index"] += 1
    return [center, *neighbors]


def generate_case_set(
    regions: list[dict],
    anchors: list[dict],
    *,
    power: int,
    seed: int,
) -> list[dict]:
    """Assemble region, independent holdout, history, and trusted samples."""
    if not isinstance(regions, list) or not isinstance(anchors, list):
        raise TypeError("regions and anchors must be lists")
    cases = []
    region_ids = set()
    for region in regions:
        region_id = str(region.get("region_id", "")).strip()
        if not region_id or region_id in region_ids:
            raise ValueError("region_id values must be non-empty and unique")
        region_ids.add(region_id)
        sample_kind = region.get("sample_kind", "sobol_interior")
        region_seed = region.get("seed", seed)
        if sample_kind == "sobol_holdout" and region_seed == seed:
            raise ValueError("Sobol holdout seed must differ from calibration seed")
        cases.extend(generate_region_cases(
            region,
            power=power,
            seed=region_seed,
            sample_kind=sample_kind,
        ))
    for anchor in anchors:
        cases.extend(generate_anchor_cases(
            anchor,
            power=power,
            seed=anchor.get("seed", seed),
        ))

    labels = [case["label"] for case in cases]
    if len(labels) != len(set(labels)):
        raise ValueError("generated case labels must be unique")
    return cases


def _safe_float(value):
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def serialize_log_frequency(values) -> list[float | None]:
    """Serialize Fast's existing log10(Hz) grid without transforming it again."""
    return [_safe_float(value) for value in values]


def success_runtime_ratio(old: dict, new: dict) -> float | None:
    """Return a timing ratio only when both solvers completed the same case."""
    if old.get("status") != "success" or new.get("status") != "success":
        return None
    old_time = old.get("elapsed_s")
    new_time = new.get("elapsed_s")
    if old_time is None or new_time is None or float(old_time) <= 0.0:
        return None
    return float(new_time) / float(old_time)


def _json_safe(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return _safe_float(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "item"):
        return _json_safe(value.item())
    return str(value)


def _parameter_domain_reason(params: dict) -> str | None:
    """Reject inputs outside the model's physical parameter definitions."""
    required = (
        "Omega_bh2", "Omega_ch2", "H0", "DN_eff", "A_s", "r", "cr",
        "T_re", "kappa10",
    )
    if any(name not in params for name in required):
        return None
    try:
        values = {name: float(params[name]) for name in required}
    except (TypeError, ValueError, OverflowError):
        return None
    if not all(math.isfinite(value) for value in values.values()):
        return "cosmological parameters must be finite"
    if values["H0"] <= 0.0:
        return "H0 must be positive"
    if values["T_re"] <= 0.0:
        return "T_re must be positive"
    for name in ("Omega_bh2", "Omega_ch2"):
        if values[name] < 0.0:
            return f"{name} is a physical density and cannot be negative"
    if values["A_s"] < 0.0 or values["r"] < 0.0:
        return "A_s and r must be non-negative power parameters"
    if values["kappa10"] < 0.0:
        return "kappa10 is a density ratio and cannot be negative"
    cr = values["cr"]
    if cr <= 0.0:
        try:
            dn_re = float(params["DN_re"])
            n_t = float(params["n_t"])
        except (KeyError, TypeError, ValueError, OverflowError):
            return None
        if not math.isfinite(n_t):
            return "n_t must be finite when cr <= 0"
        if not math.isfinite(dn_re):
            return "DN_re must be finite when cr <= 0"
        if dn_re < 0.0:
            return "DN_re is a duration and cannot be negative"
    elif values["r"] <= 0.0:
        return "r must be positive when cr > 0"
    from stiffgwpy_fast.global_param import Neff0

    if Neff0 + values["DN_eff"] < 0.0:
        return "DN_eff makes the total effective neutrino count negative"
    return None


def _worker_evaluate(case: dict, max_n_inf: float, threads: int, fast, LCDM_SG) -> dict:
    params = case["params"]
    digest = parameter_fingerprint(params)
    if digest != case.get("parameter_sha256"):
        raise ValueError(f"case parameter SHA-256 mismatch: {case.get('label')}")

    result = {
        "label": case["label"],
        "sample_kind": case["sample_kind"],
        "region_id": case["region_id"],
        "sample_index": case["sample_index"],
        "seed": case["seed"],
        "parameter_sha256": digest,
        "status": "failed",
        "classification": "solver_failure",
        "reason": None,
        "DN_gw": None,
        "log_frequency": None,
        "log10OmegaGW": None,
        "elapsed_s": None,
        "rss_mb": None,
    }
    started = time.perf_counter()
    try:
        domain_reason = _parameter_domain_reason(params)
        if domain_reason is not None:
            result.update(
                status="mathematically_invalid",
                classification="mathematically_invalid",
                reason=domain_reason,
            )
            return result
        model = LCDM_SG(**params)
        derived = model.derived_param
        n_inf = _safe_float(derived.get("N_inf"))
        result["derived"] = {
            "N_inf": n_inf,
            "n_t_effective": _safe_float(derived.get("nt")),
            "DN_re_effective": _safe_float(model.cosmo_param.get("DN_re")),
        }
        if n_inf is None or n_inf <= 0.0:
            result.update(
                status="mathematically_invalid",
                classification="mathematically_invalid",
                reason="derived N_inf is undefined, non-finite, or non-positive",
            )
            return result
        if n_inf > max_n_inf:
            result.update(
                status="resource_exceeded",
                classification="resource_exceeded",
                reason=f"N_inf={n_inf:.8g} exceeds preflight limit {max_n_inf:g}",
            )
            return result

        config = build_comparison_config(fast.FastSolverConfig, threads)
        result["solver_config"] = asdict(config)
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            solve_args = {
                "tol": 1.0e-6,
                "freq_res": 1.0,
                "config": config,
                "physical_guard": True,
            }
            if case.get("profile_stages"):
                from stiffgwpy_fast import exact_background, freq_adaptive

                stage_functions = {
                    "gen_fast": (fast, "gen_fast"),
                    "prep_frequency_only": (fast, "prep_frequency_only"),
                    "prep_frequency_kernel": (fast, "prep_frequency_kernel"),
                    "solve_kernel": (fast, "solve_kernel"),
                    "int_SGWB_W": (fast, "int_SGWB_W"),
                    "fast_phi_s2_split": (exact_background, "fast_phi_s2_split"),
                    "goal_oriented_freqs": (freq_adaptive, "goal_oriented_freqs"),
                    "grid_independent_freqs": (freq_adaptive, "grid_independent_freqs"),
                }
                output, result["stage_timings_s"] = profile_stage_calls(
                    stage_functions,
                    lambda: fast.SGWB_iter_fast(model, **solve_args),
                )
            elif case.get("profile_python"):
                output, result["python_profile"] = profile_callable(
                    fast.SGWB_iter_fast, model, **solve_args
                )
            else:
                output = fast.SGWB_iter_fast(model, **solve_args)
        result["console_output"] = captured.getvalue().strip() or None
        result["reason"] = getattr(model, "fast_failure_reason", None)
        result["failure_detail"] = getattr(model, "fast_failure_detail", None)
        if output is None:
            if result["reason"] in (
                "invalid_r", "invalid_cutoff", "invalid_parameter_domain",
                "shared_Neff_guard",
            ):
                result["status"] = "mathematically_invalid"
                result["classification"] = "mathematically_invalid"
            else:
                result["classification"] = result["reason"] or "solver_failure"
            return result

        import numpy as np

        log_spectrum = np.asarray(model.log10OmegaGW, dtype=float)
        frequency = np.asarray(model.f, dtype=float)
        if log_spectrum.shape != frequency.shape:
            result["reason"] = "frequency and spectrum shapes differ"
            return result
        physical_spectrum = np.asarray(model.Ogw_today) - np.asarray(model.Oj_today)
        dn_gw = float(model.DN_gw[-1])
        exact_zero = (
            dn_gw == 0.0
            and np.all(physical_spectrum == 0.0)
            and np.isneginf(log_spectrum).all()
        )
        if exact_zero or (np.isfinite(log_spectrum).all() and np.all(physical_spectrum > 0.0)):
            result["status"] = "success"
            result["classification"] = (
                "exact_zero_tensor_spectrum" if exact_zero else "finite_positive_spectrum"
            )
        else:
            result["reason"] = "returned non-finite or degenerate spectrum"
            result["classification"] = "returned_nonfinite_output"
        result["DN_gw"] = _safe_float(dn_gw)
        result["log_frequency"] = serialize_log_frequency(frequency)
        result["log10OmegaGW"] = [_safe_float(value) for value in log_spectrum]
        result["frequency_count"] = int(frequency.size)
        result["finite_spectrum_bins"] = int(np.count_nonzero(np.isfinite(log_spectrum)))
        result["local_error_budget"] = _json_safe(
            getattr(model, "local_error_budget", None)
        )
    except Exception as exc:  # A solver exception remains a numerical failure record.
        result["reason"] = repr(exc)
        result["classification"] = "solver_exception"
    finally:
        result["elapsed_s"] = time.perf_counter() - started
        try:
            import psutil

            result["rss_mb"] = psutil.Process().memory_info().rss / (1024.0**2)
        except ImportError:
            result["rss_mb"] = None
    return result


def _worker_main(source_root: Path, threads: int, max_n_inf: float) -> int:
    source_root = source_root.resolve()
    sys.path.insert(0, str(source_root))
    os.chdir(source_root)
    from stiffgwpy_fast import fast_sgwb
    from stiffgwpy_fast.stiff_SGWB import LCDM_SG

    for line in sys.stdin:
        captured = io.StringIO()
        try:
            case = json.loads(line)
            with contextlib.redirect_stdout(captured):
                result = _worker_evaluate(case, max_n_inf, threads, fast_sgwb, LCDM_SG)
        except Exception as exc:
            result = {
                "label": None,
                "status": "failed",
                "classification": "worker_exception",
                "reason": repr(exc),
            }
        if captured.getvalue().strip():
            result["worker_console_output"] = captured.getvalue().strip()
        print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
    return 0


def _resolve_git_commit(repo_root: Path, ref: str) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "--verify", f"{ref}^{{commit}}"],
        cwd=repo_root,
        text=True,
    ).strip()


def _tree_fingerprint(source_root: Path) -> str:
    digest = hashlib.sha256()
    package = source_root / "stiffgwpy_fast"
    for path in sorted(package.rglob("*")):
        if path.is_file():
            relative = path.relative_to(source_root).as_posix().encode("utf-8")
            digest.update(relative + b"\0")
            with path.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
    return digest.hexdigest()


def _extract_revision(repo_root: Path, commit: str, destination: Path) -> Path:
    archive = subprocess.run(
        ["git", "archive", "--format=tar", commit, "stiffgwpy_fast"],
        cwd=repo_root,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
        members = tar.getmembers()
        for member in members:
            target = (destination / member.name).resolve()
            if root not in target.parents and target != root:
                raise ValueError(f"source archive path escapes destination: {member.name}")
            if member.issym() or member.islnk():
                raise ValueError(f"source archive contains a link: {member.name}")
        tar.extractall(destination, members=members)
    return destination


class _PersistentWorker:
    def __init__(
        self,
        source_root: Path,
        threads: int,
        max_n_inf: float,
        timeout_s: float,
    ):
        self.source_root = source_root
        self.threads = threads
        self.max_n_inf = max_n_inf
        self.timeout_s = timeout_s
        self.process = None
        self.responses = None
        self.reader = None
        self.generation = 0
        self._start()

    def _start(self) -> None:
        self.responses = queue.Queue()
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker-source-root",
            str(self.source_root),
            "--worker-threads",
            str(self.threads),
            "--worker-max-n-inf",
            str(self.max_n_inf),
        ]
        self.started_at = time.perf_counter()
        self.generation += 1
        self.process = subprocess.Popen(
            command,
            cwd=self.source_root,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=None,
            text=True,
            encoding="utf-8",
            bufsize=1,
        )

        process = self.process
        responses = self.responses

        def read_responses():
            for line in process.stdout:
                responses.put(line)
            responses.put(None)

        self.reader = threading.Thread(target=read_responses, daemon=True)
        self.reader.start()
        self.first_response = True

    def request(self, case: dict) -> dict:
        if self.process is None or self.process.poll() is not None:
            self._start()
        started = time.perf_counter()
        try:
            self.process.stdin.write(
                json.dumps(case, ensure_ascii=False, allow_nan=False) + "\n"
            )
            self.process.stdin.flush()
            line = self.responses.get(timeout=self.timeout_s)
        except queue.Empty:
            self.close()
            return {
                "label": case["label"],
                "parameter_sha256": case["parameter_sha256"],
                "status": "resource_exceeded",
                "classification": "resource_timeout",
                "reason": f"worker exceeded {self.timeout_s:g} second timeout",
                "elapsed_s": None,
                "rss_mb": None,
                "worker_cold_start": False,
                "worker_generation": self.generation,
                "request_wall_s": time.perf_counter() - started,
            }
        if line is None:
            exit_code = self.process.poll()
            self.close()
            return {
                "label": case["label"],
                "parameter_sha256": case["parameter_sha256"],
                "status": "failed",
                "classification": "worker_terminated",
                "reason": f"worker exited with status {exit_code}",
                "elapsed_s": None,
                "rss_mb": None,
                "worker_cold_start": False,
                "worker_generation": self.generation,
                "request_wall_s": time.perf_counter() - started,
            }
        result = json.loads(line)
        result["worker_cold_start"] = self.first_response
        result["worker_generation"] = self.generation
        result["request_wall_s"] = time.perf_counter() - started
        if self.first_response:
            result["cold_start_wall_s"] = time.perf_counter() - self.started_at
            self.first_response = False
        return result

    def close(self) -> None:
        if self.process is None:
            return
        if self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        self.process = None


def _percentiles(values: list[float]) -> dict:
    if not values:
        return {"p50": None, "p95": None}
    import numpy as np

    return {
        "p50": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
    }


def _run_paired_scan(
    *,
    cases_path: Path,
    baseline_ref: str,
    candidate_ref: str | None,
    output_path: Path,
    threads: int,
    max_n_inf: float,
    repeats: int,
    timeout_s: float,
) -> dict:
    repo_root = Path(__file__).resolve().parents[1]
    cases_manifest, cases = read_case_manifest(cases_path)
    baseline_commit = _resolve_git_commit(repo_root, baseline_ref)
    current_commit = _resolve_git_commit(repo_root, candidate_ref or "HEAD")
    warm_times = {
        side: {"success": [], "not_success": []}
        for side in ("baseline", "candidate")
    }
    paired_ratios = []
    all_records = []
    classifications = Counter()
    rss_max = {"baseline": None, "candidate": None}
    cold_starts = {"baseline": [], "candidate": []}
    batch_started = time.perf_counter()

    with tempfile.TemporaryDirectory(prefix="fast-domain-pair-") as temporary:
        temporary_root = Path(temporary)
        baseline_root = _extract_revision(
            repo_root, baseline_commit, temporary_root / "baseline"
        )
        candidate_root = temporary_root / "candidate"
        if candidate_ref is None:
            candidate_root.mkdir()
            shutil.copytree(repo_root / "stiffgwpy_fast", candidate_root / "stiffgwpy_fast")
            candidate_source = "current worktree package snapshot"
        else:
            candidate_root = _extract_revision(
                repo_root, current_commit, candidate_root
            )
            candidate_source = "Git revision package snapshot"
        source_hashes = {
            "baseline_package_sha256": _tree_fingerprint(baseline_root),
            "candidate_package_sha256": _tree_fingerprint(candidate_root),
        }
        workers = {
            "baseline": _PersistentWorker(
                baseline_root, threads, max_n_inf, timeout_s
            ),
            "candidate": _PersistentWorker(
                candidate_root, threads, max_n_inf, timeout_s
            ),
        }
        try:
            for case_index, case in enumerate(cases):
                results = {"baseline": [], "candidate": []}
                for repeat in range(repeats):
                    order = (
                        ("baseline", "candidate")
                        if (case_index + repeat) % 2 == 0
                        else ("candidate", "baseline")
                    )
                    for side in order:
                        result = workers[side].request(case)
                        if result.get("parameter_sha256") != case["parameter_sha256"]:
                            raise ValueError(
                                f"{side} worker changed input identity for {case['label']}"
                            )
                        results[side].append(result)
                        if result.get("elapsed_s") is not None:
                            if result.get("worker_cold_start"):
                                cold_starts[side].append(
                                    result.get("cold_start_wall_s")
                                )
                            else:
                                outcome = (
                                    "success" if result.get("status") == "success"
                                    else "not_success"
                                )
                                warm_times[side][outcome].append(
                                    float(result["elapsed_s"])
                                )
                        rss = result.get("rss_mb")
                        if rss is not None:
                            rss_max[side] = max(rss_max[side] or 0.0, float(rss))
                    baseline_result = results["baseline"][-1]
                    candidate_result = results["candidate"][-1]
                    ratio = success_runtime_ratio(baseline_result, candidate_result)
                    if (
                        repeat > 0
                        and ratio is not None
                        and not baseline_result.get("worker_cold_start")
                        and not candidate_result.get("worker_cold_start")
                    ):
                        paired_ratios.append(ratio)

                baseline_result = results["baseline"][-1]
                candidate_result = results["candidate"][-1]
                comparison = compare_result_pair(baseline_result, candidate_result)
                reference_validation = None
                if case.get("reference") is not None:
                    reference_validation = {
                        "baseline": compare_reference_result(
                            baseline_result, case["reference"]
                        ),
                        "candidate": compare_reference_result(
                            candidate_result, case["reference"]
                        ),
                        "reference_artifact_sha256": case["reference"][
                            "artifact_sha256"
                        ],
                        "source_sha256": case["reference"]["source_sha256"],
                    }
                classifications[comparison["classification"]] += 1
                all_records.append({
                    "record_type": "paired_case",
                    "label": case["label"],
                    "sample_kind": case["sample_kind"],
                    "region_id": case["region_id"],
                    "sample_index": case["sample_index"],
                    "seed": case["seed"],
                    "parameter_sha256": case["parameter_sha256"],
                    "classification": comparison["classification"],
                    "comparison": comparison,
                    "reference_validation": reference_validation,
                    "baseline": baseline_result,
                    "candidate": candidate_result,
                    "baseline_elapsed_s_samples": [
                        row.get("elapsed_s") for row in results["baseline"]
                    ],
                    "candidate_elapsed_s_samples": [
                        row.get("elapsed_s") for row in results["candidate"]
                    ],
                    "baseline_status_samples": [
                        row.get("status") for row in results["baseline"]
                    ],
                    "candidate_status_samples": [
                        row.get("status") for row in results["candidate"]
                    ],
                })
        finally:
            for worker in workers.values():
                worker.close()

    elapsed_wall = time.perf_counter() - batch_started
    metadata = {
        "record_type": "manifest",
        "schema_version": 1,
        "baseline_ref": baseline_ref,
        "baseline_commit": baseline_commit,
        "candidate_commit": current_commit,
        "candidate_ref": candidate_ref,
        "candidate_source": candidate_source,
        **source_hashes,
        "case_manifest": str(cases_path.resolve()),
        "case_manifest_sha256": hashlib.sha256(cases_path.read_bytes()).hexdigest(),
        "case_set_sha256": cases_manifest["cases_sha256"],
        "case_count": len(cases),
        "threads": threads,
        "environment": environment_snapshot(threads),
        "max_n_inf_preflight": max_n_inf,
        "repeats_per_case": repeats,
        "worker_timeout_s": timeout_s,
        "solver_config": {
            "entrypoint": "SGWB_iter_fast",
            "tol": 1.0e-6,
            "freq_res": 1.0,
            "physical_guard": True,
            "FastSolverConfig": {
                "h": 0.005,
                "col_step": 8,
                "z_tail": 5.0,
                "phase_max": 0.25,
                "freq_grid": "goal",
                "threads": threads,
                "kink_split": True,
                "transition_refine_radius": 2.0,
                "transition_refine_factor": 8,
                "parallel_chunksize": None,
            },
        },
        "precision_note": "Fast-to-Fast consistency is not independent accuracy certification.",
    }
    summary = {
        "warm_runtime_s_by_status": {
            side: {status: _percentiles(values) for status, values in outcomes.items()}
            for side, outcomes in warm_times.items()
        },
        "paired_success_candidate_to_baseline_runtime_ratio": _percentiles(paired_ratios),
        "cold_start_wall_s": cold_starts,
        "peak_worker_rss_mb": rss_max,
        "classification_counts": dict(classifications),
        "warm_measurement_count": {
            side: {
                status: len(values)
                for status, values in outcomes.items()
            }
            for side, outcomes in warm_times.items()
        },
        "paired_success_ratio_count": {
            "paired_ratio": len(paired_ratios),
        },
        "batch_wall_s": elapsed_wall,
        "batch_throughput_cases_per_s": len(cases) / elapsed_wall if elapsed_wall else None,
        "accuracy_certification": "not established by paired Fast outputs",
    }
    output_path = output_path if output_path.is_absolute() else repo_root / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as stream:
        stream.write(json.dumps(metadata, ensure_ascii=False, allow_nan=False) + "\n")
        for record in all_records:
            stream.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
    summary_path = output_path.with_suffix(".summary.json")
    summary_path.write_text(
        json.dumps({"manifest": metadata, **summary}, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    return {"manifest": metadata, **summary}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, help="frozen JSONL case manifest")
    parser.add_argument("--case-spec", type=Path, help="JSON region/anchor case specification")
    parser.add_argument(
        "--baseline-ref",
        default="872a70dc5ab3bb2ff9d3c0d78f05850341b42a45",
    )
    parser.add_argument(
        "--candidate-ref",
        help="optional candidate Git ref; default snapshots the current worktree package",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-n-inf", type=float, default=120.0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--timeout-s", type=float, default=1800.0)
    parser.add_argument("--worker-source-root", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-threads", type=int, default=2, help=argparse.SUPPRESS)
    parser.add_argument("--worker-max-n-inf", type=float, default=120.0, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.worker_source_root is not None:
        return _worker_main(
            args.worker_source_root,
            args.worker_threads,
            args.worker_max_n_inf,
        )
    if (args.cases is None) == (args.case_spec is None):
        parser.error("provide exactly one of --cases or --case-spec")
    if args.output is None:
        parser.error("--output is required")
    if args.threads < 1 or args.max_n_inf <= 0.0 or args.timeout_s <= 0.0:
        parser.error("threads, max-n-inf, and timeout-s must be positive")
    if args.repeats < 2:
        parser.error("--repeats must be at least 2 to separate cold and warm calls")

    cases_path = args.cases
    if args.case_spec is not None:
        specification = json.loads(args.case_spec.read_text(encoding="utf-8"))
        cases = generate_case_set(
            specification.get("regions", []),
            specification.get("anchors", []),
            power=int(specification["power"]),
            seed=int(specification["seed"]),
        )
        for history in specification.get("history_sources", []):
            historical_cases = load_historical_cases(
                history["git_ref"],
                history["artifact_path"],
                history["sample_indices"],
            )
            cases.extend(expand_historical_cases(
                historical_cases,
                history.get("neighborhood_axes", {}),
                power=int(history.get("neighborhood_power", specification["power"])),
                seed=int(history.get("neighborhood_seed", specification["seed"])),
            ))
        cases_path = args.output.with_name(args.output.stem + ".cases.jsonl")
        write_case_manifest(
            cases_path,
            cases,
            metadata={
                "case_spec_sha256": hashlib.sha256(
                    args.case_spec.read_bytes()
                ).hexdigest(),
                "seed": specification["seed"],
                "power": specification["power"],
            },
        )

    summary = _run_paired_scan(
        cases_path=cases_path,
        baseline_ref=args.baseline_ref,
        candidate_ref=args.candidate_ref,
        output_path=args.output,
        threads=args.threads,
        max_n_inf=args.max_n_inf,
        repeats=args.repeats,
        timeout_s=args.timeout_s,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
