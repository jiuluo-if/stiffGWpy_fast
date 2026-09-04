"""Cobaya integration smoke tests (kept separate from the fast unit suite)."""

import importlib
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


def test_adapter_exports_canonical_derived_names():
    pytest.importorskip("cobaya")
    from stiffgwpy.cobaya.stiffGW import stiffGW

    expected = {
        "Delta_Neff_GW", "Delta_Neff_total",
        "log10hc_prim_fyr", "f_end", "Delta_Neff_GW_error",
    }
    adapter = stiffGW()
    assert set(adapter.get_can_provide_params()) == expected
    assert set(adapter.DERIVED_PARAMS) == expected
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] /
                          "stiffgwpy" / "cobaya" / "stiffGW.yaml").read_text(
                              encoding="utf-8"))
    assert tuple(cfg["params"]) == adapter.DERIVED_PARAMS


def test_adapter_exposes_error_telemetry():
    pytest.importorskip("cobaya")
    from stiffgwpy.cobaya.stiffGW import stiffGW

    class Model:
        fast_evals = 5
        fast_failures = 0
        fast_guard_rejections = 0
        fast_physical_rejections = 0
        lsoda_evals = 0
        lsoda_fallbacks = 0
        reference_evals = 1
        escalations = 1
        escalated_from = "production"
        DN_gw_error = 1.0e-2
        spectrum_error = 0.07

    adapter = stiffGW()
    adapter.stiffGW_model = Model()
    stats = adapter.engine_stats
    assert stats['reference_evals'] == 1
    assert stats['escalations'] == 1
    assert stats['escalated_from'] == 'production'
    assert stats['DN_gw_error'] == 1.0e-2
    assert stats['spectrum_error'] == 0.07


def test_reference_engine_sets_state(monkeypatch):
    pytest.importorskip("cobaya")
    import numpy as np

    from stiffgwpy.cobaya.stiffGW import stiffGW

    freqs = np.array([6.0, 4.0, 2.0, 0.0, -2.0, -4.0])
    logO = np.array([-15.6, -15.6, -11.6, -7.9, -11.7, -13.5])
    fake = dict(
        freqs=freqs, log10OmegaGW=logO, DN_eff=2.27e-3, DN_gw=2.27e-3,
        kappa_r=1.94e-3, g2=2.81e-8, n_iter=2)

    def fake_run(m, **kw):
        return fake

    monkeypatch.setattr("stiffgwpy.reference.run_reference", fake_run)
    adapter = stiffGW()
    adapter.engine = "reference"
    adapter.freq_res = 1.0
    m = adapter.stiffGW_model
    params = {k: 0.0 for k in ("A_s", "r", "n_t", "cr", "DN_re", "kappa10")}
    params.update(Omega_bh2=0.022, Omega_ch2=0.12, H0=67.0, DN_eff=0.0,
                  T_re=2e3)
    state = {}
    assert adapter.calculate(state, want_derived=True, **params)
    assert np.allclose(state['f'], freqs)
    assert np.allclose(state['omGW_stiff'], logO)
    assert state['derived']['Delta_Neff_GW'] == pytest.approx(2.27e-3, rel=1e-9)
    assert m.reference_evals == 1


def test_qualified_module_import_resolves_class():
    pytest.importorskip("cobaya")
    module = importlib.import_module("stiffgwpy.cobaya.stiffGW")
    cls = getattr(module, "stiffGW")
    assert cls.__module__ + "." + cls.__name__ == \
        "stiffgwpy.cobaya.stiffGW.stiffGW"


def test_yaml_modes_supply_science_defaults_via_sentinel():
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] /
                          "stiffgwpy" / "cobaya" / "stiffGW.yaml").read_text(
                              encoding="utf-8"))
    assert cfg["accuracy_mode"] == "production"
    # The science defaults now come from the named accuracy_mode preset, not
    # from the yaml: a 0 sentinel means "use the preset" so the yaml can never
    # silently override accuracy_mode.
    assert cfg["h"] == 0
    assert cfg["col_step"] == 0
    assert cfg["z_tail"] == 0
    assert cfg["freq_res"] == 0
    assert cfg["fast_threads"] == 8


def test_close_warns_on_high_fallback_fraction(monkeypatch):
    pytest.importorskip("cobaya")
    from stiffgwpy.cobaya.stiffGW import stiffGW

    class Log:
        def __init__(self):
            self.warning_calls = []

        def info(self, *args, **kwargs):
            pass

        def warning(self, *args, **kwargs):
            self.warning_calls.append((args, kwargs))

    class Model:
        fast_evals = 10
        fast_failures = 1
        lsoda_evals = 1
        lsoda_fallbacks = 1
        last_engine = "lsoda"
        last_fast_error = "guard"

    adapter = stiffGW()
    adapter.log = Log()
    adapter.stiffGW_model = Model()
    adapter.close()
    assert adapter.log.warning_calls


def test_engine_stats_exposes_failure_fraction():
    pytest.importorskip("cobaya")
    from stiffgwpy.cobaya.stiffGW import stiffGW

    class Model:
        fast_evals = 20
        fast_failures = 3
        fast_guard_rejections = 2
        fast_physical_rejections = 1
        lsoda_evals = 3
        lsoda_fallbacks = 3

    adapter = stiffGW()
    adapter.stiffGW_model = Model()
    stats = adapter.engine_stats
    assert stats['fast_failure_fraction'] == pytest.approx(0.15)
    assert stats['fallback_fraction'] == pytest.approx(0.15)


@pytest.mark.cobaya
def test_qualified_name_real_cobaya_smoke_from_other_cwd(tmp_path):
    """Run Cobaya 3.6.x through the installed qualified class name.

    The probe executes from a temporary working directory and performs a real
    Cobaya run, exercising import resolution, theory calculation, derived
    parameters, and the fast/fallback path.
    """
    cobaya = pytest.importorskip("cobaya")
    version = tuple(int(x) for x in cobaya.__version__.split('.')[:2])
    assert version == (3, 6), f"smoke contract targets Cobaya 3.6.x, got {cobaya.__version__}"
    repo = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, str(repo / "scripts" / "_smoke_probe.py")],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        timeout=180,
    )
    assert proc.returncode == 0, proc.stdout + "\n" + proc.stderr
    assert "stiffgwpy.cobaya.stiffGW.stiffGW" in proc.stdout
    assert "OK" in proc.stdout

def test_engine_stats_includes_eval_status():
    pytest.importorskip("cobaya")
    from stiffgwpy.cobaya.stiffGW import stiffGW

    class Model:
        fast_evals = 3
        fast_failures = 0
        fast_guard_rejections = 0
        fast_physical_rejections = 0
        lsoda_evals = 0
        lsoda_fallbacks = 0
        reference_evals = 0
        escalations = 1
        escalated_from = "production"
        last_eval_status = "FAST_ESCALATED"
        eval_status_counts = {"FAST": 2, "FAST_ESCALATED": 1}
        dlogl_estimated = 3.0e-3
        Delta_Neff_abs_error = 1.0e-7
        DN_gw_error = 1.0e-2
        spectrum_error = 0.07

    adapter = stiffGW()
    adapter.stiffGW_model = Model()
    stats = adapter.engine_stats
    assert stats['last_eval_status'] == 'FAST_ESCALATED'
    assert stats['eval_status_counts'] == {"FAST": 2, "FAST_ESCALATED": 1}
    assert stats['escalation_fraction'] == pytest.approx(1.0 / 3.0)
    assert stats['dlogl_estimated'] == 3.0e-3
    assert stats['Delta_Neff_abs_error'] == 1.0e-7


def _fake_fast_solve(m, **kw):
    """Fill a model with minimal solve + telemetry so the wrapper's local
    budget and likelihood-aware gate can run without a real solve."""
    import numpy as np
    m.cosmo_param['DN_eff'] = 0.00227
    m.DN_gw = [0.00227]
    m.SGWB_converge = True
    m.f = np.array([-7.0, -6.0, -5.0, -4.0])
    m.log10OmegaGW = np.array([-15.0, -11.0, -8.0, -12.0])
    m.kappa_r = 0.0019
    m.handoff_eps = np.array([5.0e-4, 5.0e-4, 5.0e-4])
    m.z_tail_used = 8.0
    m.phase_max_used = 0.5
    m.freq_grid_used = 'adaptive'
    m.freq_grid_error = 1.0e-5
    m.quadrature_error_local = 1.0e-5
    m.floating_point_error = 1.0e-12
    m.dn_converged_delta = 2.0e-9
    m.freq_res_used = 1.0
    m.freq_grid_n = 247
    m.transition_refine_used = True
    m.sigma_exact_used = False
    return m


def test_auto_escalate_likelihood_gate_escalates(monkeypatch):
    """With a tight likelihood sigma the estimated |Delta logL| must exceed
    dlogl_tol and the evaluation must be marked FAST_ESCALATED (never a
    silent fallback)."""
    pytest.importorskip("cobaya")
    from stiffgwpy import fast_sgwb
    from stiffgwpy.stiff_SGWB import LCDM_SG

    monkeypatch.setattr(fast_sgwb, 'SGWB_iter_fast', _fake_fast_solve)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m.SGWB_iter(engine='fast', accuracy_mode='production', auto_escalate=True,
                likelihood_sigma=1.0e-6, dlogl_tol=1.0e-3)
    assert m.last_eval_status == 'FAST_ESCALATED'
    assert m.escalations == 1
    assert m.escalated_from == 'production'
    assert m.dlogl_estimated is not None
    assert m.dlogl_estimated > 1.0e-3


def test_auto_escalate_likelihood_gate_stays_fast(monkeypatch):
    """With a loose likelihood sigma the estimated |Delta logL| stays below
    the budget and the evaluation remains a plain FAST."""
    pytest.importorskip("cobaya")
    from stiffgwpy import fast_sgwb
    from stiffgwpy.stiff_SGWB import LCDM_SG

    monkeypatch.setattr(fast_sgwb, 'SGWB_iter_fast', _fake_fast_solve)
    m = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    m.SGWB_iter(engine='fast', accuracy_mode='production', auto_escalate=True,
                likelihood_sigma=1.0e-2, dlogl_tol=1.0e-3)
    assert m.last_eval_status == 'FAST'
    assert m.escalations == 0
    assert m.dlogl_estimated is not None
    assert m.dlogl_estimated < 1.0e-3
    assert m.eval_status_counts['FAST'] == 1


def test_adapter_forwards_likelihood_kwargs(monkeypatch):
    pytest.importorskip("cobaya")
    import numpy as np

    from stiffgwpy import fast_sgwb
    from stiffgwpy.cobaya.stiffGW import stiffGW

    captured = {}

    def fake_sgwb_iter(self, **kw):
        captured.update(kw)
        self.SGWB_converge = True
        self.f = np.array([-5.0, -6.0, -7.0])
        self.log10OmegaGW = np.array([-8.0, -11.0, -15.0])
        self.DN_gw = [0.00227]
        self.kappa_r = 0.0019
        self.cosmo_param['DN_eff'] = 0.00227
        return self

    from stiffgwpy.stiff_SGWB import LCDM_SG
    monkeypatch.setattr(LCDM_SG, 'SGWB_iter', fake_sgwb_iter)
    adapter = stiffGW()
    adapter.engine = 'fast'
    adapter.accuracy_mode = 'production'
    adapter.auto_escalate = True
    adapter.escalate_to_reference = True
    adapter.error_tol = 1.0e-3
    adapter.likelihood_sigma = 0.25
    adapter.dlogl_tol = 0.01
    monkeypatch.setattr(fast_sgwb, 'max_threads', lambda: 4)
    adapter.fast_threads = 8
    params = {k: 0.0 for k in ("A_s", "r", "n_t", "cr", "DN_re", "kappa10")}
    params.update(Omega_bh2=0.022, Omega_ch2=0.12, H0=67.0, DN_eff=0.0,
                  T_re=2e3)
    state = {}
    assert adapter.calculate(state, want_derived=True, **params)
    assert captured.get('likelihood_sigma') == 0.25
    assert captured.get('dlogl_tol') == 0.01
    assert captured.get('error_tol') == 1.0e-3
    assert captured.get('auto_escalate') is True
    assert captured.get('escalate_to_reference') is True
    assert captured.get('threads') == 4




