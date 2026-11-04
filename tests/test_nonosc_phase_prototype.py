"""Regression checks for the standalone phase-transfer research prototype."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from scripts import benchmark_nonosc_phase_res4_fullgrid as phase_grid
from scripts import prototype_nonosc_phase_res4_mode as phase


def test_phase_handoff_at_tail_threshold_uses_reference_tail_even_if_today_z_falls_below(
    monkeypatch,
):
    for name in ("model", "bg", "NSTART", "Z0", "NCUT", "NCUT100"):
        monkeypatch.setattr(phase, name, getattr(phase, name))
    params = dict(phase.PARAMS)
    params.update(
        r=0.19649188544593907,
        cr=0.0,
        n_t=1.0622501865308032,
        T_re=1761.314039059193,
        kappa10=0.011598483151630803,
        DN_re=10.0,
    )
    monkeypatch.setattr(phase, "PARAMS", params)
    monkeypatch.setattr(phase, "ZTAIL", 5.0)
    phase.configure_case(2.291411867022825, -16.318140930557433, 5.0)

    low = phase.solve_ivp(
        phase.low_rhs,
        (phase.NSTART, phase.bg.N_inf),
        (phase.Z0, 1.0, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=[1e-12, 1e-13, 1e-13],
        events=[phase.low_event],
        dense_output=False,
    )
    assert low.success and low.t_events[0].size
    n_phase = float(low.t_events[0][0])
    low_state = low.y_events[0][0].tolist()

    _, bounds, _, _, tail_used, z_end = phase.candidate(
        16, low_state, n_phase, phase.ZPHASE
    )

    assert tail_used is True
    assert bounds[-1][0] == pytest.approx(n_phase, abs=1e-12)
    assert z_end == pytest.approx(phase.ZTAIL, abs=1e-12)


def test_fullgrid_driver_accepts_a_valid_zero_segment_tail_transfer(monkeypatch):
    root_result = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "docs"
            / "sample2_phase_reference_root_res8_tail5_20261104.json"
        ).read_text(encoding="utf-8")
    )
    reference_row = next(row for row in root_result["rows"] if row["index"] == 646)
    phase_module = phase_grid.phase
    for name in ("model", "bg", "NSTART", "Z0", "NCUT", "NCUT100"):
        monkeypatch.setattr(phase_module, name, getattr(phase_module, name))
    monkeypatch.setattr(
        phase_module,
        "PARAMS",
        {key: value for key, value in root_result["params"].items() if key != "DN_eff"},
    )
    monkeypatch.setattr(phase_module, "ZTAIL", 5.0)
    monkeypatch.setattr(phase_module, "PHASE_GAMMA_MIN", 10.0)

    result = phase_grid.run_mode(
        {
            "index": reference_row["index"],
            "freq_log10_hz": reference_row["freq_log10_hz"],
            "cartesian_log10_omega": reference_row["reference_log10_omega"],
            "used_tail": reference_row["reference_used_tail"],
        },
        root_result["params"],
    )

    assert result["status"] == "completed"
    assert result["tail_semantics_match"] is True
    assert result["orders"][2]["segments"] == []
    assert result["max_dense_residual_scaled"] is None
    assert result["orders"][2]["abs_dex_error_vs_cartesian"] < 1e-3


@pytest.mark.parametrize(("index", "order"), [(656, 24), (687, 16)])
def test_phase_collocation_stops_once_existing_residual_gate_is_met(monkeypatch, index, order):
    root_grid = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "docs"
            / "sample2_phase_reference_root_res8_tail5_gamma0_20261105.json"
        ).read_text(encoding="utf-8")
    )
    probe = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "docs"
            / "sample2_phase_root_residual_order_probe_20261105.json"
        ).read_text(encoding="utf-8")
    )
    row = next(item for item in root_grid["rows"] if item["index"] == index)
    probe_row = next(item for item in probe["rows"] if item["index"] == index)
    module = phase
    params = {key: value for key, value in root_grid["params"].items() if key != "DN_eff"}
    monkeypatch.setattr(module, "PARAMS", params)
    monkeypatch.setattr(module, "ZTAIL", root_grid["z_tail"])
    monkeypatch.setattr(module, "PHASE_GAMMA_MIN", 0.0)
    phase_start = probe_row["phase_start"]
    module.configure_case(root_grid["params"]["DN_eff"], row["freq_log10_hz"], phase_start)
    low = module.solve_ivp(
        module.low_rhs,
        (module.NSTART, module.bg.N_inf),
        (module.Z0, 1.0, 0.0),
        method="DOP853",
        rtol=1e-12,
        atol=[1e-12, 1e-13, 1e-13],
        events=[module.low_event],
        dense_output=False,
    )
    assert low.success and low.t_events[0].size

    state, bounds, segments, _, tail_used, z_end = module.candidate(
        order,
        low.y_events[0][0].tolist(),
        float(low.t_events[0][0]),
        phase_start,
    )
    end_n = bounds[-1][0]
    power = (
        module.spectrum_from_tail(state, end_n, z_end)
        if tail_used
        else module.spectrum_today(state, z_end)
    )
    assert max(item["collocation_residual_scaled"] for item in segments) <= 1e-10
    assert max(item["dense_residual_scaled"] for item in segments) <= 1e-9
    assert any(item["newton_residual_converged"] for item in segments)
    assert abs(math.log10(power["OmegaGW"]) - row["reference_log10_omega"]) < 1e-3
