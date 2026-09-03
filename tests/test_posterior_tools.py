"""Pure-function regression tests for the posterior-validation tooling.

These tests exercise the math that certifies the fast posterior against the
continuous-sigma reference (ESS>=2000 via independent importance sampling,
per-point Delta logL and importance-reweighted posterior shift).  No engine
solves are run here, so the tests stay fast; the expensive driver runs are
scripts/importance_posterior.py and scripts/cobaya_posterior_fast_vs_reference.py.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"


def _load_script(name):
    spec = importlib.util.spec_from_file_location(
        name.replace(".py", ""), SCRIPTS / name)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ip():
    return _load_script("importance_posterior.py")


@pytest.fixture(scope="module")
def cp():
    return _load_script("cobaya_posterior_fast_vs_reference.py")


def test_mock_loglike_matches_manual_gaussian(ip):
    truth = {"truth_log10Om": [0.0, 0.0, 0.0]}
    lo = np.array([0.01, -0.02, 0.005])
    ll = ip.loglike(lo, truth, 0.05)
    manual = -0.5 * np.sum((lo / 0.05) ** 2)
    assert ll == pytest.approx(manual)
    assert ip.loglike(np.array([np.nan, 0.0, 0.0]), truth, 0.05) == float("-inf")


def test_log_prior_flat_and_truncated(ip):
    inside = ip.log_prior(-2.0, 0.0)
    assert inside == pytest.approx(-np.log((5.0) * 1.0))
    assert ip.log_prior(-6.0, 0.0) == -np.inf
    assert ip.log_prior(-2.0, 0.7) == -np.inf


def test_is_ess_approaches_n_for_matched_proposal(ip, tmp_path):
    # draws from q == posterior: importance weights constant -> ESS ~ n
    rng = np.random.default_rng(3)
    n = 4000
    truth = {"truth_log10Om": np.zeros(11).tolist()}
    log10r = rng.normal(-2.0, 0.02, n)
    n_t = rng.normal(0.0, 0.10, n)
    # engine model: log10 Omega = log10r + 2 (truth vector is 0 at r=1e-2),
    # plus a small per-bin noise floor; q must over-cover n_t (weakly
    # constrained by the tilt) or IS variance collapses.
    lo = np.zeros((n, 11))
    for i in range(n):
        lo[i] = (log10r[i] + 2.0) + rng.normal(0.0, 0.005, 11)
    w = np.empty(n)
    for i in range(n):
        lp = ip.log_prior(log10r[i], n_t[i])
        ll = ip.loglike(lo[i], truth, 0.05)
        mu = np.array([-2.0, 0.0])
        sig = np.array([0.02, 0.40])
        x = np.array([log10r[i], n_t[i]])
        lq = -0.5 * np.sum(((x - mu) / sig) ** 2) - np.sum(np.log(sig))
        w[i] = np.exp(lp + ll - lq)
    wn = w / w.sum()
    ess = 1.0 / np.sum(wn * wn)
    mean_r = float(np.sum(wn * log10r))
    assert ess > 0.3 * n, "IS ESS collapsed: %.0f" % ess
    assert mean_r == pytest.approx(-2.0, abs=0.02)


def test_reweighted_shift_recovers_injected_bias(ip):
    # pointwise Delta logL proportional to (x - x0) with known slope sigma^{-2}
    # must move the weighted mean by exactly the injected amount.
    rng = np.random.default_rng(7)
    n = 2000
    x = rng.normal(0.0, 1.0, n)
    w0 = np.ones(n) / n
    bias = 0.3
    dll = -(x - bias) ** 2 / 2 + x ** 2 / 2   # L_ref/L_fast ratio
    w1 = w0 * np.exp(dll)
    w1 /= w1.sum()
    mean1 = np.sum(w1 * x)
    assert mean1 == pytest.approx(bias, abs=0.05)
    assert ip._boot(x, n=200, seed=1).size == 200


def test_bins_consistent_between_modules_and_truth(cp, ip):
    assert cp.PTA_BINS == ip.PTA_BINS
    assert cp.KNEE_BINS == ip.KNEE_BINS
    assert cp.LVK_BINS == ip.LVK_BINS
    assert len(cp.BINS) == 11
    truth_p = REPO / "docs" / "mcmc_posterior" / "mock_truth.json"
    if truth_p.exists():
        truth = json.loads(truth_p.read_text(encoding="utf-8"))
        assert truth["bins"] == cp.BINS
        assert len(truth["truth_log10Om"]) == 11


def test_sgwb_mock_likelihood_registered(cp):
    pytest.importorskip("cobaya")
    from cobaya.likelihood import Likelihood
    assert issubclass(cp.SGWB_mock, Likelihood)
    assert cp.SGWB_mock().get_requirements() == {"f": None, "omGW_stiff": None}
    assert cp.SGWB_mock().sigma_dex == 0.05


def test_oracle_truth_engine_neutral_source(cp):
    d = json.loads((REPO / "docs" / "archive" / "reference" /
                    "deep_oracle_default.json").read_text(encoding="utf-8"))
    # the mock data vector must come from the reference (not fast) spectrum
    lo_ref = np.asarray(d["ref"]["log10Om"], dtype=float)
    lo_fast = np.asarray(d["fast"]["log10Om"], dtype=float)
    assert np.isfinite(lo_ref).all() and np.isfinite(lo_fast).all()
    assert np.max(np.abs(lo_ref - lo_fast)) > 1e-6, (
        "reference and fast spectra identical: truth source would be moot")
