# -*- coding: utf-8 -*-
"""Cobaya theory adapter for :class:`stiffgwpy.stiff_SGWB.LCDM_SG`.

The theory is registered under the fully-qualified class name
``stiffgwpy.cobaya.stiffGW.stiffGW`` so it can be used from any working
directory after ``pip install .`` (no ``python_path`` needed):

.. code-block:: yaml

   theory:
     stiffgwpy.cobaya.stiffGW.stiffGW:
       engine: fast
       fallback: True
       accuracy_mode: production
       fast_threads: 8

The class exposes these derived parameters (must stay in sync with
:meth:`get_can_provide_params` and the ``params`` section of ``stiffGW.yaml``):

* ``Delta_Neff_GW``     -- SGWB contribution to Delta N_eff today
* ``Delta_Neff_total``  -- total Delta N_eff after the self-consistent loop
* ``log10hc_prim_fyr``  -- log10 h_c of the primordial SGWB at f = 1/yr
* ``f_end``             -- UV cutoff frequency (Hz) of the spectrum

Per-run observability: the underlying ``LCDM_SG`` instance counts fast
evaluations / fast failures / LSODA fallbacks and records the engine used on
the last solve, so a run can report the fallback fraction (see
:attr:`engine_stats`).
"""

import math

import astropy.units as u
import numpy as np
from cobaya.theory import Theory
from scipy import interpolate


class stiffGW(Theory):
    # Type annotations only: the default *values* live in stiffGW.yaml (the
    # same-name file Cobaya loads next to this module).  Declaring values here
    # AND in the yaml raises a duplicate-option error in Cobaya 3.6.
    speed = 0.1
    engine: str
    fallback: bool
    fast_threads: int
    h: float
    col_step: int
    z_tail: float
    freq_res: float
    accuracy_mode: str
    auto_escalate: bool
    error_tol: float
    # Likelihood-aware escalation: |Delta logL| = 0.5*(Delta_Neff_abs_error /
    # likelihood_sigma)^2 is compared against dlogl_tol (natural log-units).
    # The effective sigma is the likelihood's sensitivity to Delta N_eff
    # (e.g. ~0.25 for CMB/BAO-scale Neff constraints; set per likelihood).
    likelihood_sigma: float
    dlogl_tol: float
    escalate_to_reference: bool
    reference_rtol: float
    reference_z_tail: float
    # Canonical public names.  In particular, do not use the ambiguous
    # historical ``Delta_Neff`` alias: GW-only and total contributions are
    # distinct quantities and are kept explicit throughout the adapter.
    DERIVED_PARAMS = ('Delta_Neff_GW', 'Delta_Neff_total',
                      'log10hc_prim_fyr', 'f_end', 'Delta_Neff_GW_error')

    def initialize(self):
        """Create the model object used for every parameter evaluation."""
        try:
            from ..stiff_SGWB import LCDM_SG
        except ImportError:  # imported as a top-level module by cobaya
            from stiffgwpy.stiff_SGWB import LCDM_SG
        self.stiffGW_model = LCDM_SG()
        self.log.info("Initialized!")

    def initialize_with_provider(self, provider):
        self.provider = provider

    def close(self):
        # Cobaya calls ``close`` when a run finishes.  Emit the accumulated
        # telemetry so numerical fallbacks are visible in MCMC logs rather
        # than silently disappearing behind a successful LSODA retry.
        stats = self.engine_stats
        if stats is not None:
            counts = stats['eval_status_counts']
            self.log.info(
                "SGWB engine summary: fast_evals=%d fast_failures=%d "
                "fast_guard_rejections=%d fast_physical_rejections=%d "
                "lsoda_evals=%d lsoda_fallbacks=%d "
                "fast_failure_fraction=%.6g fallback_fraction=%.6g "
                "escalation_fraction=%.6g eval_status=%s",
                stats['fast_evals'], stats['fast_failures'],
                stats['fast_guard_rejections'],
                stats['fast_physical_rejections'], stats['lsoda_evals'],
                stats['lsoda_fallbacks'],
                stats['fast_failure_fraction'],
                stats['fallback_fraction'],
                stats['escalation_fraction'],
                {k: v for k, v in counts.items() if v})
            if (stats['fallback_fraction'] > 0.05
                    or stats['escalation_fraction'] > 0.05):
                self.log.warning(
                    "SGWB fast fallback fraction %.3f exceeds 5%%; "
                    "inspect fast numerical failures before using this chain.",
                    stats['fallback_fraction'])

    def get_requirements(self):
        return {'Omega_bh2': None, 'Omega_ch2': None, 'H0': None, 'DN_eff': None,
                'A_s': None, 'r': None, 'n_t': None, 'cr': None,
                'T_re': None, 'DN_re': None, 'kappa10': None}

    def get_can_provide(self):
        return ['f', 'omGW_stiff', 'hubble', 'kappa_s', 'kappa_r']

    def get_can_provide_params(self):
        # Must be exactly the keys written into state['derived'] in calculate().
        return list(self.DERIVED_PARAMS)

    @property
    def engine_stats(self):
        """Fast/fallback counters accumulated on the shared model instance."""
        m = getattr(self, 'stiffGW_model', None)
        if m is None:
            return None
        return dict(fast_evals=getattr(m, 'fast_evals', 0),
                    fast_failures=getattr(m, 'fast_failures', 0),
                    lsoda_evals=getattr(m, 'lsoda_evals', 0),
                    lsoda_fallbacks=getattr(m, 'lsoda_fallbacks', 0),
                    last_engine=getattr(m, 'last_engine', None),
                    last_fast_error=getattr(m, 'last_fast_error', None),
                    fast_failure_reason=getattr(m, 'fast_failure_reason', None),
                    fast_guard_rejections=getattr(m, 'fast_guard_rejections', 0),
                    fast_physical_rejections=getattr(
                        m, 'fast_physical_rejections', 0),
                    reference_evals=getattr(m, 'reference_evals', 0),
                    escalations=getattr(m, 'escalations', 0),
                    escalated_from=getattr(m, 'escalated_from', None),
                    DN_gw_error=getattr(m, 'DN_gw_error', None),
                    spectrum_error=getattr(m, 'spectrum_error', None),
                    fast_failure_fraction=_failure_fraction(m),
                    fallback_fraction=_fallback_fraction(m),
                    last_eval_status=getattr(m, 'last_eval_status', None),
                    eval_status_counts=dict(getattr(m, 'eval_status_counts', {})),
                    escalation_fraction=_escalation_fraction(m),
                    dlogl_estimated=getattr(m, 'dlogl_estimated', None),
                    Delta_Neff_abs_error=getattr(m, 'Delta_Neff_abs_error', None))

    def calculate(self, state, want_derived=True, **params_values_dict):
        """Evaluate the model; on failure return False (logpost -> -inf)."""
        m = self.stiffGW_model
        m.reset()
        for key in m.cosmo_param:
            if key in params_values_dict:
                m.cosmo_param[key] = params_values_dict[key]

        sgwb_kwargs = {}
        if self.accuracy_mode:
            sgwb_kwargs['accuracy_mode'] = self.accuracy_mode
        if self.h:
            sgwb_kwargs['h'] = self.h
        if self.col_step:
            sgwb_kwargs['col_step'] = self.col_step
        if self.z_tail:
            sgwb_kwargs['z_tail'] = self.z_tail
        if self.freq_res and self.freq_res != 1.0:
            sgwb_kwargs['freq_res'] = self.freq_res
        if self.auto_escalate:
            sgwb_kwargs['auto_escalate'] = True
            if getattr(self, 'escalate_to_reference', False):
                sgwb_kwargs['escalate_to_reference'] = True
        if self.error_tol:
            sgwb_kwargs['error_tol'] = self.error_tol
        if getattr(self, 'likelihood_sigma', None):
            sgwb_kwargs['likelihood_sigma'] = self.likelihood_sigma
        if getattr(self, 'dlogl_tol', None):
            sgwb_kwargs['dlogl_tol'] = self.dlogl_tol
        if self.engine == 'fast' and self.fast_threads:
            try:
                from .. import fast_sgwb
            except ImportError:
                from stiffgwpy import fast_sgwb
            # Clamp to the machine's thread budget so an over-large yaml
            # default cannot hard-fail on small cores.
            fast_sgwb.set_threads(min(int(self.fast_threads),
                                      fast_sgwb._MAX_THREADS))
        if self.engine == 'reference':
            # Independent continuous-sigma high-accuracy pipeline (slow; for
            # certification/benchmark points, not for MCMC thermal path).
            from ..reference import apply_reference_to_model
            apply_reference_to_model(
                m, freq_res=getattr(self, 'freq_res', 1.0) or 1.0,
                z_tail=getattr(self, 'reference_z_tail', 5.0) or 5.0,
                rtol=getattr(self, 'reference_rtol', 1e-11) or 1e-11)
        else:
            m.SGWB_iter(engine=self.engine, fallback=self.fallback, **sgwb_kwargs)

        if not m.SGWB_converge:
            # logpost will be -inf; make sure any previous derived values are
            # not reused for this point.
            state.pop('derived', None)
            return False
        self.last_eval_status = getattr(m, 'last_eval_status', 'UNKNOWN')

        state['f'] = m.f                                     # log10(f/Hz)
        state['omGW_stiff'] = m.log10OmegaGW                 # log10 Omega_GW(f)
        state['hubble'] = m.derived_param['H_0']             # s^-1
        state['kappa_s'] = m.derived_param['kappa_s']
        state['kappa_r'] = m.kappa_r
        if want_derived:
            yr = u.yr.to(u.s)
            log10f_yr = -math.log10(yr)
            if m.f[0] >= log10f_yr:
                f_t = np.flip(state['f'])
                Ogw_t = np.flip(state['omGW_stiff'])
                spec_prim = interpolate.CubicSpline(f_t[f_t > -13],
                                                    Ogw_t[f_t > -13])
                omGW_stiff_fyr = spec_prim(log10f_yr)
            else:
                omGW_stiff_fyr = -100.0
            derived = {
                'Delta_Neff_GW': m.DN_gw[-1],
                'Delta_Neff_total': m.cosmo_param['DN_eff'],
                'log10hc_prim_fyr': (omGW_stiff_fyr / 2
                                     + math.log10(math.sqrt(1.5) *
                                                  state['hubble'] / math.pi)
                                     - log10f_yr),
                'f_end': np.power(10., m.f[0]),
                'Delta_Neff_GW_error': getattr(m, 'DN_gw_error', 0.0),
            }
            if tuple(derived) != self.DERIVED_PARAMS:
                raise RuntimeError('derived parameter contract drift: %s != %s'
                                   % (tuple(derived), self.DERIVED_PARAMS))
            state['derived'] = derived
        else:
            # Cobaya may reuse the state dictionary between calls; never leave
            # derived values from a previous parameter point visible when the
            # caller did not request them.
            state.pop('derived', None)
        return True


def _fallback_fraction(m):
    fast_total = getattr(m, 'fast_evals', 0)
    if fast_total <= 0:
        return 0.0
    return getattr(m, 'lsoda_fallbacks', 0) / float(fast_total)


def _failure_fraction(m):
    fast_total = getattr(m, 'fast_evals', 0)
    if fast_total <= 0:
        return 0.0
    return getattr(m, 'fast_failures', 0) / float(fast_total)


def _escalation_fraction(m):
    total = getattr(m, 'fast_evals', 0)
    if total <= 0:
        return 0.0
    return getattr(m, 'escalations', 0) / float(total)
