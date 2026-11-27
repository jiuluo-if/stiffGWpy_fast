from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.fast_domain_compare import (
    build_comparison_config,
    classify_pair,
    compare_reference_result,
    compare_result_pair,
    environment_snapshot,
    expand_historical_cases,
    generate_anchor_cases,
    generate_case_set,
    generate_region_cases,
    historical_record_to_case,
    parameter_fingerprint,
    parse_historical_jsonl,
    profile_callable,
    profile_stage_calls,
    read_case_manifest,
    reference_from_artifact,
    serialize_log_frequency,
    success_runtime_ratio,
    write_case_manifest,
)
from stiffgwpy_fast import FastSolverConfig


def test_classify_pair_separates_solver_and_domain_outcomes():
    cases = [
        ("failed", "success", "old_failure_to_new_success"),
        ("success", "failed", "old_success_to_new_failure"),
        ("success", "success", "persistent_success"),
        ("failed", "failed", "persistent_failure"),
        ("mathematically_invalid", "mathematically_invalid", "mathematically_invalid"),
        ("resource_exceeded", "success", "resource_exceeded"),
        ("mathematically_invalid", "success", "classification_disagreement"),
    ]

    for old_status, new_status, expected in cases:
        assert classify_pair(old_status, new_status) == expected


def test_comparison_config_freezes_the_same_fast_numerics_per_worker():
    from dataclasses import asdict

    config = build_comparison_config(FastSolverConfig, threads=2)

    assert asdict(config) == {
        "h": 0.005,
        "col_step": 8,
        "z_tail": 5.0,
        "phase_max": 0.25,
        "freq_grid": "goal",
        "threads": 2,
        "kink_split": True,
        "transition_refine_radius": 2.0,
        "transition_refine_factor": 8,
        "parallel_chunksize": None,
    }


def test_environment_snapshot_records_runtime_and_thread_controls():
    snapshot = environment_snapshot(threads=2)

    assert snapshot["python"]
    assert snapshot["platform"]
    assert snapshot["logical_cpu_count"] >= 1
    assert snapshot["threads_argument"] == 2
    assert set(snapshot["library_versions"]) == {"numpy", "scipy", "numba"}


def test_profile_callable_returns_result_and_bounded_python_profile():
    value, profile = profile_callable(sum, [1, 2, 3])

    assert value == 6
    assert "function calls" in profile
    assert "builtins.sum" in profile


def test_profile_stage_calls_times_real_wrapped_functions():
    module = SimpleNamespace(add=lambda left, right: left + right)

    result, timings = profile_stage_calls(
        {"add": (module, "add")},
        lambda: module.add(2, 3),
    )

    assert result == 5
    assert set(timings) == {"add"}
    assert timings["add"] >= 0.0


def test_compare_result_pair_requires_exact_input_and_reports_spectrum_delta():
    params_sha = parameter_fingerprint({"cr": 0.0, "r": 0.01})
    old = {
        "parameter_sha256": params_sha,
        "status": "success",
        "DN_gw": 2.0,
        "log_frequency": [-2.0, -1.0],
        "log10OmegaGW": [-6.0, -5.0],
    }
    new = {
        **old,
        "DN_gw": 2.2,
        "log10OmegaGW": [-6.01, -4.98],
    }

    comparison = compare_result_pair(old, new)

    assert comparison["classification"] == "persistent_success"
    assert comparison["DN_gw_relative_difference"] == pytest.approx(0.1)
    assert comparison["spectrum_max_abs_difference_dex"] == pytest.approx(0.02)

    new["parameter_sha256"] = "different-input"
    with pytest.raises(ValueError, match="identical parameter SHA-256"):
        compare_result_pair(old, new)


def test_log_frequency_serializer_preserves_existing_log10_coordinates():
    assert serialize_log_frequency([-14.5, -1.0, 0.0, 2.0]) == [
        -14.5, -1.0, 0.0, 2.0
    ]


def test_failed_old_path_has_no_paired_success_runtime_ratio():
    assert success_runtime_ratio(
        {"status": "failed", "elapsed_s": 0.001},
        {"status": "success", "elapsed_s": 0.02},
    ) is None
    assert success_runtime_ratio(
        {"status": "success", "elapsed_s": 0.01},
        {"status": "success", "elapsed_s": 0.0125},
    ) == pytest.approx(1.25)


def test_historical_record_keeps_failure_class_and_does_not_trust_finite_output():
    params = {"cr": 0.0, "r": 0.01, "n_t": 0.5}
    source = {
        "commit": "a" * 40,
        "artifact_path": "docs/fast_domain_pressure_boundary_20261104.jsonl",
        "artifact_sha256": "b" * 64,
        "suite_seed": 20261030,
    }
    failed_record = {
        "sample_index": 118,
        "label": "boundary_A_s_cr0_0",
        "params": params,
        "parameter_sha256": parameter_fingerprint(params),
        "status": "failed",
        "classification": "nonfinite_numerical_state",
        "reason": "nonfinite",
    }
    finite_record = {
        **failed_record,
        "sample_index": 2,
        "status": "success",
        "classification": "finite_positive_spectrum",
    }

    failure_case = historical_record_to_case(failed_record, source)
    finite_case = historical_record_to_case(finite_record, source)

    assert failure_case["sample_kind"] == "historical_failure"
    assert failure_case["historical_classification"] == "nonfinite_numerical_state"
    assert finite_case["sample_kind"] == "historical_output_unverified"
    assert finite_case["historical_classification"] == "finite_positive_spectrum"


def test_historical_jsonl_parser_selects_exact_indices_and_records_source_hash():
    params = {"cr": 0.0, "r": 0.01, "n_t": 0.5}
    record = {
        "sample_index": 118,
        "label": "boundary_A_s_cr0_0",
        "params": params,
        "parameter_sha256": parameter_fingerprint(params),
        "status": "failed",
        "classification": "nonfinite_numerical_state",
        "reason": "nonfinite",
    }
    raw = "\n".join([
        json.dumps({"record_type": "manifest", "commit": "f" * 40, "suite": "boundaries", "seed": 33}),
        json.dumps(record),
    ])

    cases = parse_historical_jsonl(
        raw,
        artifact_ref="872a70dc",
        artifact_path="docs/fast_domain_pressure_boundary_20261104.jsonl",
        sample_indices=[118],
    )

    assert len(cases) == 1
    assert cases[0]["historical_classification"] == "nonfinite_numerical_state"
    assert cases[0]["history_source"]["scan_source_commit"] == "f" * 40
    assert len(cases[0]["history_source"]["artifact_sha256"]) == 64


def test_historical_failures_expand_to_source_linked_neighborhoods():
    params = {"cr": 0.0, "r": 0.01, "n_t": 0.5}
    record = {
        "sample_index": 118,
        "label": "boundary_A_s_cr0_0",
        "params": params,
        "parameter_sha256": parameter_fingerprint(params),
        "status": "failed",
        "classification": "nonfinite_numerical_state",
        "reason": "nonfinite",
    }
    source = {"commit": "a" * 40, "suite": "boundaries", "suite_seed": 33}
    original = historical_record_to_case(record, source)

    cases = expand_historical_cases(
        [original],
        {"r": {"scale": "log", "factor": 1.1}},
        power=1,
        seed=99,
    )

    assert len(cases) == 5
    assert cases[0]["sample_kind"] == "historical_failure"
    assert all(case["history_source"] == source for case in cases)
    assert all(
        case["historical_classification"] == "nonfinite_numerical_state"
        for case in cases
    )
    assert all(case["sample_kind"] == "historical_failure_neighborhood" for case in cases[1:])


def test_parameter_fingerprint_ignores_mapping_insertion_order():
    first = {"r": 1.0e-8, "cr": 0.0, "n_t": 0.5}
    second = {"n_t": 0.5, "r": 1.0e-8, "cr": 0.0}

    assert parameter_fingerprint(first) == parameter_fingerprint(second)


def test_cr_positive_sobol_cases_use_the_consistency_relation():
    region = {
        "region_id": "cr1_probe",
        "cr": 1.0,
        "fixed": {"r": 0.01, "n_t": 0.5, "T_re": 2000.0},
        "axes": {"Omega_bh2": [0.02, 0.03]},
    }

    cases = generate_region_cases(region, power=3, seed=1741)

    interior = [case for case in cases if case["sample_kind"] == "sobol_interior"]
    boundary = [case for case in cases if case["sample_kind"] == "region_boundary"]
    assert len(interior) == 8
    assert len(boundary) == 2
    assert all(case["params"]["n_t"] == -case["params"]["r"] / 8.0 for case in cases)
    assert [case["params"] for case in cases] == [
        case["params"] for case in generate_region_cases(region, power=3, seed=1741)
    ]


def test_region_boundary_cases_hold_other_axes_at_the_midpoint():
    region = {
        "region_id": "cr0_boundary_probe",
        "cr": 0.0,
        "fixed": {"r": 0.01, "T_re": 2000.0},
        "axes": {
            "n_t": [-0.2, 0.4],
            "kappa10": {"bounds": [1.0e-4, 1.0e-2], "scale": "log"},
        },
    }

    cases = generate_region_cases(region, power=1, seed=1741)
    boundary = [case for case in cases if case["sample_kind"] == "region_boundary"]

    assert len(boundary) == 4
    assert sorted(case["params"]["n_t"] for case in boundary) == pytest.approx(
        [-0.2, 0.1, 0.1, 0.4]
    )
    assert sorted(case["params"]["kappa10"] for case in boundary) == pytest.approx(
        [1.0e-4, 1.0e-3, 1.0e-3, 1.0e-2]
    )


def test_sobol_mapping_is_independent_of_axis_mapping_order():
    region = {
        "region_id": "stable_axis_order",
        "cr": 0.0,
        "fixed": {"r": 0.01, "T_re": 2000.0},
        "axes": {"n_t": [-0.2, 0.4], "Omega_bh2": [0.01, 0.04]},
    }
    reversed_region = {**region, "axes": dict(reversed(list(region["axes"].items())))}

    first = generate_region_cases(region, power=2, seed=1741)
    second = generate_region_cases(reversed_region, power=2, seed=1741)

    assert [case["params"] for case in first] == [case["params"] for case in second]


def test_historical_anchor_neighborhood_uses_derived_cr_positive_tilt():
    anchor = {
        "label": "legacy_failure_129",
        "sample_kind": "historical_failure",
        "region_id": "blue_tilt_boundary",
        "params": {"cr": 1.0, "r": 0.01, "n_t": 0.5, "T_re": 2000.0},
        "axes": {"r": {"bounds": [0.008, 0.012], "scale": "log"}},
    }

    cases = generate_anchor_cases(anchor, power=2, seed=1741)
    center = cases[0]

    assert center["sample_kind"] == "historical_failure"
    assert center["input_params"]["n_t"] == 0.5
    assert center["params"]["n_t"] == 0.5
    assert center["n_t_effective"] == -0.01 / 8.0
    assert center["parameter_sha256"] == parameter_fingerprint(center["params"])
    assert all(case["sample_kind"] == "historical_failure_neighborhood" for case in cases[1:])
    assert all(case["params"]["n_t"] == -case["params"]["r"] / 8.0 for case in cases[1:])


def test_anchor_neighborhood_requires_a_supported_sample_kind():
    anchor = {
        "label": "unknown_anchor",
        "sample_kind": "random_point",
        "params": {"cr": 0.0, "r": 0.01, "n_t": 0.2},
        "axes": {"n_t": [0.1, 0.3]},
    }

    with pytest.raises(ValueError, match="sample_kind"):
        generate_anchor_cases(anchor, power=1, seed=1741)


def test_anchor_neighborhood_accepts_relative_log_and_linear_radii():
    anchor = {
        "label": "failure_nearby",
        "sample_kind": "historical_failure",
        "params": {"cr": 0.0, "r": 0.01, "n_t": 0.2, "T_re": 2000.0},
        "axes": {
            "r": {"scale": "log", "factor": 1.1},
            "T_re": {"scale": "linear", "delta": 50.0},
        },
    }

    cases = generate_anchor_cases(anchor, power=1, seed=1741)
    neighbors = cases[1:]

    assert len(neighbors) == 2 + 2 * 2
    assert all(
        0.01 / 1.1 * (1.0 - 1.0e-14)
        <= case["params"]["r"]
        <= 0.01 * 1.1 * (1.0 + 1.0e-14)
        for case in neighbors
    )
    assert all(1950.0 <= case["params"]["T_re"] <= 2050.0 for case in neighbors)


def test_saved_reference_artifact_produces_same_grid_anchor_errors():
    root = Path(__file__).resolve().parents[1]
    path = root / (
        "tests/data/"
        "validate_cr1_t0p004_0p005_r0000000000010641088010743262_"
        "0000000000010675_unseen_holdout_T00475_rlo_z10.json"
    )
    artifact = json.loads(path.read_text(encoding="utf-8"))
    params = artifact["params"]

    reference = reference_from_artifact(path, params)
    saved_fast = artifact["fast"]
    result = {
        "status": "success",
        "DN_gw": saved_fast["root"],
        "log_frequency": saved_fast["freqs_log10_hz"],
        "log10OmegaGW": saved_fast["log10_spectrum"],
    }

    comparison = compare_reference_result(result, reference)

    assert comparison["reference_root_relative_error"] == pytest.approx(
        artifact["comparison"]["root_relative_error"]
    )
    assert comparison["reference_spectrum_max_abs_error_dex"] == pytest.approx(
        artifact["comparison"]["spectrum_max_abs_error_dex"]
    )
    assert comparison["precision_status"] == "anchor_fail"
    assert comparison["certification"] == "anchor_screen_only"


def test_case_set_keeps_cr_branches_and_independent_holdout_seed():
    regions = [
        {
            "region_id": "cr0_calibration",
            "cr": 0.0,
            "sample_kind": "sobol_interior",
            "fixed": {"r": 0.01, "T_re": 2000.0},
            "axes": {"n_t": [-0.2, 0.4]},
        },
        {
            "region_id": "cr1_holdout",
            "cr": 1.0,
            "sample_kind": "sobol_holdout",
            "seed": 9907,
            "fixed": {"r": 0.01, "T_re": 2000.0},
            "axes": {"T_re": {"bounds": [1000.0, 3000.0], "scale": "log"}},
        },
    ]

    cases = generate_case_set(regions, [], power=2, seed=1741)
    calibration = [case for case in cases if case["sample_kind"] == "sobol_interior"]
    holdout = [case for case in cases if case["sample_kind"] == "sobol_holdout"]

    assert len(calibration) == len(holdout) == 4
    assert all(case["seed"] == 1741 for case in calibration)
    assert all(case["seed"] == 9907 for case in holdout)
    assert all(case["params"]["n_t"] == -case["params"]["r"] / 8.0 for case in holdout)


def test_cr_positive_region_rejects_nt_as_an_independent_axis():
    region = {
        "region_id": "invalid_cr1_axes",
        "cr": 1.0,
        "fixed": {"r": 0.01, "T_re": 2000.0},
        "axes": {"n_t": [-0.1, 0.0]},
    }

    with pytest.raises(ValueError, match="n_t is derived"):
        generate_region_cases(region, power=2, seed=1741)


def test_case_manifest_round_trip_and_rejects_modified_parameters(tmp_path):
    region = {
        "region_id": "manifest_probe",
        "cr": 0.0,
        "fixed": {"r": 0.01, "T_re": 2000.0},
        "axes": {"n_t": [-0.2, 0.4]},
    }
    cases = generate_region_cases(region, power=1, seed=1741)
    path = tmp_path / "cases.jsonl"

    write_case_manifest(path, cases, metadata={"seed": 1741})
    manifest, loaded = read_case_manifest(path)

    assert manifest["case_count"] == len(cases)
    assert manifest["metadata"] == {"seed": 1741}
    assert loaded == cases

    lines = path.read_text(encoding="utf-8").splitlines()
    tampered = json.loads(lines[1])
    tampered["params"]["n_t"] += 0.01
    lines[1] = json.dumps(tampered)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="parameter SHA-256 mismatch"):
        read_case_manifest(path)
