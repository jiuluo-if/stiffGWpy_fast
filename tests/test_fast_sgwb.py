# -*- coding: utf-8 -*-
import numpy as np
import pytest
from numba import get_num_threads

from stiffgwpy_fast import LCDM_SG
from stiffgwpy_fast import fast_sgwb as FS


def _make_model(**kw):
    defaults = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    defaults.update(kw)
    return LCDM_SG(**defaults)


def test_set_col_step_validates():
    FS.set_col_step(4)
    for bad in (0, 9, -1, 'x'):
        with pytest.raises(ValueError):
            FS.set_col_step(bad)
    FS.set_col_step(1)
    FS.set_col_step(8)


def test_set_threads_validates():
    FS.set_threads(1)
    for bad in (0, -2, 'x'):
        with pytest.raises(ValueError):
            FS.set_threads(bad)
    with pytest.raises(ValueError):
        FS.set_threads(FS._MAX_THREADS + 1)
    FS.set_threads(FS._MAX_THREADS)


def test_set_h_validates():
    FS.set_h(0.005)
    for bad in (0, -1, 0.2, 'x'):
        with pytest.raises(ValueError):
            FS.set_h(bad)
    FS.set_h(0.01)


def test_set_z_tail_validates():
    FS.set_z_tail(7.0)
    for bad in (1.5, 20.0, -5, 'x'):
        with pytest.raises(ValueError):
            FS.set_z_tail(bad)
    FS.set_z_tail(5.0)


def test_pchip_frequency_integral_handles_descending_native_grid():
    """PCHIP candidate must integrate a linear log-frequency spectrum exactly."""
    freqs = np.array([1.0, 0.5, 0.0])
    integrand = 2.0 * freqs + 1.0
    assert FS.integrate_frequency_pchip(freqs, integrand) == pytest.approx(2.0)


def test_frequency_quadrature_methods_cover_q1_candidates():
    """Q1 methods integrate a quadratic spectrum on descending nodes."""
    freqs = np.array([1.0, 0.75, 0.5, 0.25, 0.0])
    integrand = freqs**2 + 2.0 * freqs + 1.0
    expected = 7.0 / 3.0
    methods = (
        'simpson', 'pchip', 'log_pchip', 'gauss2', 'gauss3', 'gauss5',
        'natural_cubic', 'chebyshev',
    )
    values = {
        method: FS.integrate_frequency_quadrature(freqs, integrand, method)
        for method in methods
    }
    assert values['simpson'] == pytest.approx(expected)
    assert values['pchip'] == pytest.approx(expected)
    assert values['gauss5'] == pytest.approx(expected, rel=1e-12)
    assert values['natural_cubic'] == pytest.approx(expected, rel=1e-12)
    assert values['chebyshev'] == pytest.approx(expected, rel=1e-12)
    assert all(np.isfinite(value) for value in values.values())


def test_frequency_quadrature_rejects_unknown_method():
    with pytest.raises(ValueError, match='unknown frequency quadrature'):
        FS.integrate_frequency_quadrature(
            np.array([0.0, 1.0]), np.array([1.0, 2.0]), 'unknown')


def test_local_frequency_quadrature_estimator_returns_interval_budget():
    freqs = np.array([0.0, 0.2, 0.55, 1.0, 1.6, 2.0])
    integrand = np.exp(-freqs) * (1.0 + 0.2 * np.sin(3.0 * freqs))
    errors, candidate, baseline = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'pchip')
    assert errors.shape == (freqs.size - 1,)
    assert candidate.shape == errors.shape
    assert baseline.shape == errors.shape
    assert np.all(errors >= 0.0)
    assert np.sum(candidate) == pytest.approx(
        FS.integrate_frequency_quadrature(freqs, integrand, 'pchip'))
    assert np.all(np.isfinite(baseline))
    assert np.any(errors > 0.0)


def test_local_frequency_quadrature_supports_chebyshev_candidate():
    freqs = np.linspace(0.0, 1.0, 7)
    integrand = 1.0 + freqs**2
    errors, candidate, baseline = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'chebyshev')
    assert np.all(np.isfinite(errors))
    assert np.sum(candidate) == pytest.approx(
        FS.integrate_frequency_quadrature(freqs, integrand, 'chebyshev'))
    assert np.sum(baseline) == pytest.approx(
        FS.integrate_frequency_quadrature(freqs, integrand, 'simpson'))


def test_local_frequency_quadrature_full_panel_allocation_is_conservative():
    freqs = np.array([0.0, 0.2, 0.55, 1.0, 1.6, 2.0])
    integrand = np.exp(-freqs) * (1.0 + 0.2 * np.sin(3.0 * freqs))
    errors, candidate, baseline = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'pchip', allocation='full_panel')
    for start in range(0, errors.size - 1, 2):
        expected = abs(np.sum(candidate[start:start + 2]) -
                       np.sum(baseline[start:start + 2]))
        assert errors[start] == pytest.approx(expected)
        assert errors[start + 1] == pytest.approx(expected)


def test_local_frequency_quadrature_panel_envelope_dominates_full_panel():
    freqs = np.linspace(0.0, 1.0, 9)
    integrand = 1.0 + freqs**4
    full, _, _ = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'pchip', allocation='full_panel')
    envelope, _, _ = FS.estimate_frequency_quadrature_local(
        freqs, integrand, 'pchip', allocation='panel_envelope')
    assert np.all(envelope >= full)


def test_estimator_coverage_summary_reports_quantile_calibration():
    from scripts.benchmark_quadrature_estimator_coverage import (
        summarize_estimator_coverage)

    rows = [
        {'predicted_rel': 2.0, 'actual_rel': 1.0},
        {'predicted_rel': 1.0, 'actual_rel': 2.0},
        {'predicted_rel': 4.0, 'actual_rel': 1.0},
    ]
    summary = summarize_estimator_coverage(rows, 'test')
    assert summary['n'] == 3
    assert summary['false_safe'] == 1
    assert summary['coverage'] == pytest.approx(2.0 / 3.0)
    assert summary['actual_over_prediction_p95'] > 1.0


def test_fixed_spectrum_dense_reference_is_interval_local():
    from scripts.benchmark_quadrature_estimator_coverage import (
        fixed_spectrum_dense_reference)

    freqs = np.array([0.0, 0.25, 0.8, 1.5])
    integrand = 1.0 + 2.0 * freqs
    local = fixed_spectrum_dense_reference(freqs, integrand, subdivisions=32)
    expected = ((freqs[1:] - freqs[:-1]) +
                (freqs[1:]**2 - freqs[:-1]**2))
    assert local.shape == (freqs.size - 1,)
    assert np.allclose(local, expected, rtol=1e-10, atol=1e-12)


def test_quadrature_coverage_includes_gauss_estimators():
    from scripts.benchmark_quadrature_estimator_coverage import (
        ESTIMATOR_METHODS)

    assert ESTIMATOR_METHODS == ('pchip', 'gauss2', 'gauss3', 'gauss5')


def test_ensemble_local_error_is_pointwise_maximum():
    from scripts.benchmark_same_grid_reference import ensemble_local_error

    values = [np.array([1.0, 4.0, 2.0]), np.array([3.0, 1.0, 5.0])]
    assert np.array_equal(ensemble_local_error(values),
                          np.array([3.0, 4.0, 5.0]))


def test_positive_tilt_benchmark_disables_consistency_relation():
    from scripts.benchmark_same_grid_reference import CASES

    positive = CASES['positive_tilt']
    assert positive['cr'] == 0
    assert positive['n_t'] > 0.0
    assert 'DN_re' in positive


def test_pchip_frequency_quadrature_is_opt_in():
    """PCHIP DN quadrature is explicit and leaves the default path unchanged."""
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='goal', threads=1,
                              kink_split=True)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg,
                             frequency_quadrature='pchip') is m
    assert m.frequency_quadrature_used == 'pchip'
    omega_nu = FS.gp.Omega_nh2 / m.derived_param['h']**2
    expected = (FS.gp.Neff0 * FS.ln10 *
                FS.integrate_frequency_pchip(m.f, m.Ogw_today - m.Oj_today) /
                omega_nu)
    assert m.DN_gw[-1] == pytest.approx(expected, rel=1e-12)
    assert m.estimated_DN_quadrature_error > 0.0
    assert m.estimated_DN_quadrature_error_rel > 0.0
    assert m.quadrature_error_estimator_method == 'pchip'
    assert m.quadrature_error_local_by_interval.shape == (m.f.size - 1,)
    assert np.all(m.quadrature_error_local_by_interval >= 0.0)
    assert m.quadrature_error_local == pytest.approx(
        m.estimated_DN_quadrature_error_rel)
    assert m.quadrature_error_estimator_max_rel > 0.0


def test_gauss_frequency_quadrature_is_available_in_solver():
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='goal', threads=1,
                              kink_split=True)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg,
                             frequency_quadrature='gauss3') is m
    assert m.frequency_quadrature_used == 'gauss3'
    assert np.isfinite(m.DN_gw[-1])


def test_r_le_zero_returns_none(model):
    m = _make_model(r=0.0)
    assert FS.SGWB_iter_fast(m) is None


def test_nonfinite_dn_gw_aborts_and_restores(monkeypatch, model):
    orig_dn = model.cosmo_param['DN_eff']

    def bad_solve(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
                  fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
                  Ogw, Oj, Opgw, h_arr=None, Sv=None, phase_max=0.0,
                  handoff_eps=None):
        Ogw[...] = np.nan
        Oj[...] = 0.0
        Opgw[...] = 0.0

    monkeypatch.setattr(FS, 'solve_kernel', bad_solve)
    assert FS.SGWB_iter_fast(model) is None
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.DN_eff_orig is None
    assert model.SGWB_converge is False


def test_outer_probe_skips_column_assembly_until_final_solve(monkeypatch):
    """The first outer probe only needs the final column for DN_gw."""
    calls = []
    original = FS.solve_kernel

    def wrapped(*args, **kwargs):
        calls.append(args[11])
        return original(*args, **kwargs)

    monkeypatch.setattr(FS, 'solve_kernel', wrapped)
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='construct', threads=1)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg, kink_split=True) is m
    assert len(calls) >= 2
    assert calls[0] == 0
    assert calls[-1] == 1


def test_goal_path_reuses_stable_first_full_solve(monkeypatch):
    """Formal goal fast path may skip a redundant second kernel solve."""
    calls = []
    original = FS.solve_kernel

    def wrapped(*args, **kwargs):
        calls.append(args[11])
        return original(*args, **kwargs)

    monkeypatch.setattr(FS, 'solve_kernel', wrapped)
    monkeypatch.setattr(FS, '_OUTER_FULL_REUSE_ENABLED', True)
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='goal', threads=1,
                              kink_split=True)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg) is m
    assert calls == [1]
    assert m.outer_full_reuse_used is True


def test_kink_path_uses_frequency_only_preparation(monkeypatch):
    """The formal kink path must skip unused uniform primitive preparation."""
    calls = []
    original = FS.prep_frequency_kernel

    def wrapped(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(FS, 'prep_frequency_kernel', wrapped)
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='construct', threads=1)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg, kink_split=True) is m
    assert len(calls) == 2


def test_kink_path_reuses_stable_exact_primitive(monkeypatch):
    """Small background updates may reuse the exact primitive safely."""
    from stiffgwpy_fast import exact_background as EB

    calls = []
    original = EB.fast_phi_s2_split

    def wrapped(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(EB, 'fast_phi_s2_split', wrapped)
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='construct', threads=1)
    assert FS.SGWB_iter_fast(m=_make_model(), tol=1e-6,
                             config=cfg, kink_split=True) is not None
    assert len(calls) == 1

    calls.clear()
    high = _make_model(T_re=2e4)
    assert FS.SGWB_iter_fast(high, tol=1e-6, config=cfg,
                             kink_split=True) is high
    assert len(calls) == 2


def test_per_call_thread_config_is_restored_on_failure(monkeypatch, model):
    """A per-call thread override must not leak through a failed solve."""
    before = get_num_threads()
    target = 1 if before != 1 else min(2, FS._MAX_THREADS)
    config = FS.FastSolverConfig(threads=target)

    def bad_solve(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
                  fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
                  Ogw, Oj, Opgw, h_arr=None, Sv=None, phase_max=0.0,
                  handoff_eps=None):
        Ogw[...] = np.nan
        Oj[...] = 0.0
        Opgw[...] = 0.0

    monkeypatch.setattr(FS, 'solve_kernel', bad_solve)
    assert FS.SGWB_iter_fast(model, config=config) is None
    assert get_num_threads() == before


def test_max_iterations_abort_and_restore(monkeypatch, model):
    orig_dn = model.cosmo_param['DN_eff']
    rng = np.random.default_rng(42)
    monkeypatch.setattr(FS, 'MAX_ITER', 6)

    def fake_solve(Nv, Phi_grid, Phi_mid, S2, S2inv, j0s, z0s, P_t, ev_minus,
                   fp_minus, fp_freq, assemble, n_coarse, col_step, h, z_tail,
                   Ogw, Oj, Opgw, h_arr=None, Sv=None, phase_max=0.0,
                   handoff_eps=None):
        # Produce a random DN_gw_new in (0.05, 4.9) that never converges.
        Omega_nu = FS.gp.Omega_nh2 / model.derived_param['h'] ** 2
        span = np.log10(fp_freq.max()) - np.log10(fp_freq.min())
        val = rng.uniform(0.05, 4.9)
        const = val * Omega_nu / (FS.gp.Neff0 * FS.ln10 * span)
        Ogw[...] = 0.0
        Oj[...] = 0.0
        Opgw[...] = 0.0
        Ogw[:, -1] = const

    monkeypatch.setattr(FS, 'solve_kernel', fake_solve)
    assert FS.SGWB_iter_fast(model) is None
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.DN_eff_orig is None
    assert model.SGWB_converge is False


def test_stale_full_evolution_attrs_cleared(model):
    for name in ('N_hc', 'Th', 'Oj', 'Ogw', 'Opgw'):
        setattr(model, name, [np.array([1.0])])
    assert FS.SGWB_iter_fast(model) is model
    assert model.SGWB_converge
    for name in ('N_hc', 'Th', 'Oj', 'Ogw', 'Opgw'):
        assert not hasattr(model, name), name


def test_success_leaves_dn_eff_orig_for_reset(model):
    orig_dn = model.cosmo_param['DN_eff']
    FS.SGWB_iter_fast(model)
    assert model.SGWB_converge
    assert model.DN_eff_orig == orig_dn
    final_dn = model.cosmo_param['DN_eff']
    assert final_dn > orig_dn
    model.reset()
    assert model.cosmo_param['DN_eff'] == orig_dn
    assert model.SGWB_converge is False


def test_no_nan_and_deterministic_under_memory_churn():
    """Regression for the stale-heap NaN: solve_kernel starts each channel at
    j0 (horizon crossing + 3 decades), so early coarse columns are never
    written for channels with ``j0 % col_step != 0``.  int_SGWB_W still reads
    those cells, which must be zero-filled (np.zeros), not np.empty, otherwise
    reused heap contents (including NaN) leak into the g2/w2 curves."""
    FS.set_col_step(4)
    FS.set_h(0.04)
    FS.set_z_tail(5.0)
    try:
        first_dn = None
        for _ in range(3):
            # Scramble heap reuse so np.empty would return stale NaN pages.
            junk = [np.empty((2000, 2000)) for _ in range(3)]
            for arr in junk:
                arr.fill(np.nan)
            del junk
            m = _make_model()
            FS.SGWB_iter_fast(m, tol=1e-7, freq_res=1.0)
            for name in ('g2', 'w2', 'DN_gw', 'log10OmegaGW',
                         'Ogw_today', 'Opgw_today'):
                a = np.asarray(getattr(m, name), float)
                assert np.isfinite(a).all(), 'non-finite %s' % name
            dn = float(m.DN_gw[-1])
            if first_dn is None:
                first_dn = dn
            assert dn == first_dn, 'non-deterministic DN_gw across runs'
    finally:
        FS.set_col_step(1)
        FS.set_h(0.01)
        FS.set_z_tail(5.0)

# ============ physics-first regression / convergence / error budget ============

_LOCAL_BUDGET_KEYS = (
    'background_model', 'sigma_transition', 'ode_integration',
    'horizon_crossing', 'wkb_handoff', 'interpolation', 'frequency_grid',
    'quadrature', 'tail_approximation', 'floating_point', 'self_consistency')


def _solve_dn(**kw):
    """Run one fast solve on the default point and return (DN_eff, n_freq).

    ``phase_max_set`` / ``z_tail_set`` are applied as module settings (they
    are not SGWB_iter_fast kwargs); everything else is forwarded."""
    pm = kw.pop('phase_max_set', None)
    zt = kw.pop('z_tail_set', None)
    if pm is not None:
        FS.set_phase_max(pm)
    if zt is not None:
        FS.set_z_tail(zt)
    opts = dict(tol=1e-7, transition_refine=True, freq_grid='construct',
                freq_res=1.0)
    opts.update(kw)
    m = _make_model()
    FS.SGWB_iter_fast(m, **opts)
    return float(m.cosmo_param['DN_eff']), len(m.f)


def test_solve_kernel_legacy_signature_regression():
    """The pre-horizon-crossing-adaptivity 20-arg solve_kernel call must still
    be accepted: the four new trailing parameters are keyword-only with
    defaults, so old positional callers keep working unchanged."""
    import inspect
    sig = inspect.signature(FS.solve_kernel)
    params = list(sig.parameters)
    legacy = ('Nv', 'Phi_grid', 'Phi_mid', 'S2', 'S2inv', 'j0s', 'z0s', 'P_t',
              'ev_minus', 'fp_minus', 'fp_freq', 'assemble', 'n_coarse',
              'col_step', 'h', 'z_tail', 'Ogw', 'Oj', 'Opgw')
    assert params[:len(legacy)] == list(legacy)
    for name, default in (('h_arr', None), ('Sv', None), ('phase_max', 0.0),
                          ('handoff_eps', None)):
        assert sig.parameters[name].default == default


def test_phase_max_convergence(fast_settings):
    """Horizon-crossing phase sub-stepping converges: halving phase_max must
    change Delta N_eff by < 1e-5 relative and the change must not grow."""
    FS.apply_accuracy_mode('production')
    FS.set_z_tail(5.0)
    d1, _ = _solve_dn(phase_max_set=0.5)
    d2, _ = _solve_dn(phase_max_set=0.25)
    d3, _ = _solve_dn(phase_max_set=0.125)
    chg12 = abs(d2 - d1) / d1
    chg23 = abs(d3 - d2) / d2
    assert chg12 < 1e-5, chg12
    assert chg23 < 1e-5, chg23
    assert chg23 <= chg12 + 1e-9


def test_z_tail_convergence(fast_settings):
    """The analytic frozen-tail handoff must be flat in z_tail: moving the
    threshold from 7 to 8 changes Delta N_eff by < 2e-4 relative."""
    FS.apply_accuracy_mode('production')
    d7, _ = _solve_dn(z_tail_set=7.0)
    d8, _ = _solve_dn(z_tail_set=8.0)
    rel = abs(d8 - d7) / d7
    assert rel < 2e-4, rel


def test_tail_match_uses_physical_oscillation_amplitude():
    """The leading tail match uses d(ln|T|)/dN=-1, independent of sigma.

    For a sub-horizon tensor mode the physical solution has T proportional to
    a^-1 at leading order.  The amplitude quadrature is therefore x + y/omega,
    not the damping-removal coefficient p/2 used by the older match.
    """
    for sigma in (1.0, 4.0 / 3.0, 2.0):
        assert FS._tail_match_gamma(sigma) == pytest.approx(1.0)


@pytest.mark.slow
def test_freq_grid_convergence(fast_settings):
    """The adaptive frequency grid converges with freq_res: fr1 -> fr2 -> fr4
    must shrink the Delta N_eff change monotonically."""
    FS.apply_accuracy_mode('production')
    d1, n1 = _solve_dn(freq_grid='adaptive', freq_res=1.0)
    d2, n2 = _solve_dn(freq_grid='adaptive', freq_res=2.0)
    d4, n4 = _solve_dn(freq_grid='adaptive', freq_res=4.0)
    assert n1 < n2 < n4
    chg12 = abs(d2 - d1) / d1
    chg24 = abs(d4 - d2) / d2
    assert chg12 < 5e-4, chg12
    assert chg24 < chg12


def test_estimate_local_error_categories(fast_settings):
    """A real production solve must produce an 11-category local a-posteriori
    error budget with sane (positive, classified) entries."""
    FS.apply_accuracy_mode('production')
    FS.set_z_tail(8.0)
    m = _make_model()
    FS.SGWB_iter_fast(m, tol=1e-7, transition_refine=True,
                      freq_grid='adaptive', freq_res=2.0)
    b = FS.estimate_local_error(m)
    assert tuple(b['categories']) == _LOCAL_BUDGET_KEYS
    for name, cat in b['categories'].items():
        assert cat['value'] >= 0.0, name
        assert cat['kind'] in ('local', 'calibrated'), name
    assert b['DN_gw_error'] > 0.0
    assert b['Delta_Neff_abs_error'] > 0.0
    assert b['systematic_error'] >= 0.0
    assert b['random_rss'] >= 0.0
    assert b['handoff_eps_max'] >= 0.0
    assert b['z_tail_used'] == 8.0
    assert b['freq_grid_used'] == 'adaptive'


def test_estimate_local_error_missing_telemetry():
    """estimate_local_error must not raise on a model with no solve telemetry
    and must still return a conservative budget (defaults only)."""
    class Bare:
        pass
    b = FS.estimate_local_error(Bare())
    assert tuple(b['categories']) == _LOCAL_BUDGET_KEYS
    assert b['DN_gw_error'] > 0.0
    assert b['Delta_Neff_abs_error'] == 0.0
    for cat in b['categories'].values():
        assert cat['value'] >= 0.0


def test_wrapper_wires_production_preset_and_local_budget():
    """The SGWB_iter wrapper must forward the fast preset (kink split
    transition + goal frequency grid) and attach the local error budget,
    eval status and telemetry."""
    m = _make_model()
    m.SGWB_iter(engine='fast', accuracy_mode='fast', tol=1e-7)
    assert m.SGWB_converge
    assert m.transition_refine_used is False
    assert m.kink_split_used is True
    assert m.freq_grid_used == 'goal'
    assert m.last_eval_status == 'FAST'
    assert m.eval_status_counts['FAST'] == 1
    assert m.local_error_budget is not None
    assert m.DN_gw_error > 0.0
    assert m.Delta_Neff_abs_error > 0.0
    assert m.error_estimates is not None
    assert m.quadrature_error >= 0.0


def test_kink_split_variant_inserts_only_reheating_breakpoint():
    """Phase A uses one exact N_re breakpoint without global refinement."""
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.0, freq_grid='construct', threads=1)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg, kink_split=True) is m

    plain = _make_model()
    FS.gen_fast(plain, cfg.h)
    assert m.Nv.size == plain.Nv.size
    assert m.kink_split_index >= 0
    assert 0.0 < m.kink_split_fraction < 1.0
    assert m.sigma[m.kink_split_index] == pytest.approx(1.0)
    assert m.sigma[m.kink_split_index + 1] > 1.0
    assert m.kink_split_used is True
    assert m.transition_refine_used is False


def test_phase_envelope_caps_each_split_transfer_segment():
    """Phase B applies the phase cap independently on both sides of ``N_re``."""
    assert FS._phase_substeps(0.02, 4.0, 0.0) == 1
    assert FS._phase_substeps(0.02, 4.0, 0.5) >= 2
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.5, freq_grid='construct', threads=1)
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg, kink_split=True) is m
    assert m.phase_max_used == pytest.approx(0.5)
    assert np.all(np.isfinite(m.log10OmegaGW))


def test_goal_frequency_grid_solves_native_eval_nodes():
    """The fast solver accepts the sparse goal grid without post interpolation."""
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.25, freq_grid='goal', threads=1)
    ev = np.array([-2.0, -1.0, 0.0, 1.0])
    m = _make_model()
    assert FS.SGWB_iter_fast(m, tol=1e-6, config=cfg, kink_split=True,
                              eval_freqs=ev) is m
    assert 64 <= m.f.size <= 120
    assert np.all(np.min(np.abs(m.f[:, None] - ev[None, :]), axis=0) == 0.0)
    assert m.freq_grid_used == 'goal'




def test_eval_freqs_are_native_grid_nodes():
    """eval_freqs force-adds log10(f/Hz) evaluation points into the solve grid
    as native nodes: downstream point evaluations (likelihood bins) must not
    inherit spline interpolation error on steep spectral features.  Points
    must land exactly on a node (nearest-node distance 0) and each requested
    one must increase the node count."""
    ev = [-2.0, -1.0, 0.0, 1.0, 1.40]
    m0 = _make_model()
    assert FS.SGWB_iter_fast(m0, tol=1e-7, transition_refine=True,
                             freq_grid='grid_independent', freq_res=1.0) is m0
    fs0 = np.sort(np.asarray(m0.f, dtype=float))
    m1 = _make_model()
    assert FS.SGWB_iter_fast(m1, tol=1e-7, transition_refine=True,
                             freq_grid='grid_independent', freq_res=1.0,
                             eval_freqs=ev) is m1
    fs1 = np.sort(np.asarray(m1.f, dtype=float))
    assert fs1.size == fs0.size + len(ev), (fs1.size, fs0.size)
    dmin = np.min(np.abs(fs1 - np.asarray(ev)[:, None]), axis=1)
    assert np.allclose(dmin, 0.0), dmin


def test_eval_freqs_do_not_change_self_consistent_dn():
    """Likelihood-only native nodes must not perturb the bolometric DN integral."""
    cfg = FS.FastSolverConfig(h=0.02, col_step=8, z_tail=5.0,
                              phase_max=0.25, freq_grid='goal', threads=1)
    base = _make_model()
    assert FS.SGWB_iter_fast(base, tol=1e-7, config=cfg, kink_split=True) is base
    with_eval = _make_model()
    assert FS.SGWB_iter_fast(
        with_eval, tol=1e-7, config=cfg, kink_split=True,
        eval_freqs=np.array([-1.73, -0.42, 0.37, 1.21])) is with_eval
    assert with_eval.f.size > base.f.size
    assert with_eval.DN_gw[-1] == pytest.approx(base.DN_gw[-1], rel=1e-12, abs=1e-15)
    assert base.estimated_DN_quadrature_error > 0.0
    assert base.estimated_DN_quadrature_error_rel > 0.0
