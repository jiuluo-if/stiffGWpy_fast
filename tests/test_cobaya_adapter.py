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
        "log10hc_prim_fyr", "f_end",
    }
    adapter = stiffGW()
    assert set(adapter.get_can_provide_params()) == expected
    assert set(adapter.DERIVED_PARAMS) == expected
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] /
                          "stiffgwpy" / "cobaya" / "stiffGW.yaml").read_text(
                              encoding="utf-8"))
    assert tuple(cfg["params"]) == adapter.DERIVED_PARAMS


def test_qualified_module_import_resolves_class():
    pytest.importorskip("cobaya")
    module = importlib.import_module("stiffgwpy.cobaya.stiffGW")
    cls = getattr(module, "stiffGW")
    assert cls.__module__ + "." + cls.__name__ == \
        "stiffgwpy.cobaya.stiffGW.stiffGW"


def test_yaml_keeps_production_science_defaults():
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] /
                          "stiffgwpy" / "cobaya" / "stiffGW.yaml").read_text(
                              encoding="utf-8"))
    assert cfg["accuracy_mode"] == "production"
    assert cfg["h"] == 0.01
    assert cfg["col_step"] == 4
    assert cfg["z_tail"] == 7.0
    assert cfg["freq_res"] == 1.0
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
