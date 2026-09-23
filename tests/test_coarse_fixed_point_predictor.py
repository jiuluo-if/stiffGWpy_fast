"""Contract tests for the standalone coarse fixed-point predictor."""

import numpy as np

from scripts.benchmark_coarse_fixed_point_predictor import _evaluate_once, _run_case
from scripts.benchmark_phase_recurrence import CASES
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import global_param as gp
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_coarse_predictor_rejection_is_detected(fast_settings):
    row = _run_case("default", coarse_count=32, repeats=1, threads=2)
    assert row["status_equal"] is True
    assert row["deterministic_equal"] is True
    assert row["candidate_spectrum_max_dex"] > 1e-3


def test_predictor_map_uses_production_pchip_log_frequency_measure(fast_settings):
    FS.apply_accuracy_mode("fast")
    FS.set_threads(2)
    model = LCDM_SG(**CASES["default"])
    result = _evaluate_once(model, 0.0, coarse_count=None)
    integrand = result["Ogw"] - result["Oj"]
    expected_g2 = float(
        np.sum(FS._pchip_integrals_vectorized(result["freqs"], integrand))
        * FS.ln10
    )
    expected_dn = gp.Neff0 * expected_g2 / (
        gp.Omega_nh2 / model.derived_param["h"] ** 2
    )
    assert result["dn_gw"] == expected_dn
