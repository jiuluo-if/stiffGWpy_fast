import math

import astropy.units as u
import numpy as np
from cobaya.theory import Theory
from scipy import interpolate

try:
    from mpi4py import MPI
except ImportError:  # pragma: no cover - mpi4py is an optional dependency
    MPI = None

class stiffGW(Theory):
    speed = 0.1
    engine: str = 'lsoda'    # solver engine: 'lsoda' (default) or 'fast';
                             # can be overridden per-run through the theory yaml
    fallback: bool = True    # engine='fast': rerun with LSODA on failure
    fast_threads: int = 0    # engine='fast': OpenMP threads (0 = module default)
    h: float = 0.0           # engine='fast': step size (0 = module default)
    col_step: int = 0        # engine='fast': column stride (0 = module default)
    z_tail: float = 0.0      # analytic-tail threshold (0 = module default)
    freq_res: float = 1.0    # frequency-grid density (1.0 = default grid)
    accuracy_mode: str = ''  # engine='fast': 'reference'|'production'|'ultra-fast'
                             # ('' = no preset; explicit knobs still apply)
    params = {'Delta_Neff_GW': {'derived': True, 'latex': r'\Delta N_\mathrm{eff,GW}'},
              'Delta_Neff_total': {'derived': True, 'latex': r'\Delta N_\mathrm{eff,total}'},
              'log10hc_prim_fyr': {'derived': True, 'latex': r'\log_{10}h_{c,\mathrm{prim}}'},
              'f_end': {'derived': True, 'latex': r'f_\mathrm{end}'},
             }

    def initialize(self):
        """called from __init__ to initialize"""
        from ..stiff_SGWB import LCDM_SG
        self.stiffGW_model = LCDM_SG()
        #self.comm = MPI.COMM_WORLD
        #self.rank = self.comm.Get_rank()
        self.log.info("Initialized!")

    def initialize_with_provider(self, provider):
        """
        Initialization after other components initialized, using theory.Provider class instance.
        It is used to return any dependencies (requirements of this theory)
        via methods like "provider.get_X()" and "provider.get_param(‘Y’)".
        """
        self.provider = provider

    def close(self):
        pass


    def get_requirements(self):
        """
        Return dictionary of quantities that are always needed by this component
        and should be calculated by another component or provided by input parameters.
        """
        return {'Omega_bh2': None, 'Omega_ch2': None, 'H0': None, 'DN_eff': None,
                'A_s': None, 'r': None, 'n_t': None, 'cr': None,
                'T_re': None, 'DN_re': None, 'kappa10': None}

#    def must_provide(self, **requirements):
#        if 'A' in requirements:
#            # e.g. calculating A requires B computed using same kmax (default 10)
#            return {'B': {'kmax': requirements['A'].get('kmax', 10)}}

    def get_can_provide(self):
        return ['f', 'omGW_stiff', 'hubble', 'kappa_s', 'kappa_r',]

    def get_can_provide_params(self):
        return ['Delta_Neff_GW', 'Delta_Neff_total', 'log10hc_prim_fyr', 'f_end',]


    def calculate(self, state, want_derived=True, **params_values_dict):
        """
        The 'Theory.calculate()' method takes a dictionary 'params_values_dict'
        of the parameter values as keyword arguments and saves all needed results
        in the 'state' dictionary (which is cached and reused as needed).
        """

        # Set parameters
        self.stiffGW_model.reset()
        #args = {p: v for p, v in params_values_dict.items()}
        #self.log.debug("Setting parameters: %r", args)
        #print(self.rank, ": ", params_values_dict)
        for key in self.stiffGW_model.cosmo_param:
            if key in params_values_dict:
                self.stiffGW_model.cosmo_param[key] = params_values_dict[key]

        # Compute!  The engine is configurable through the theory yaml
        # ('lsoda' by default; 'fast' opts into the experimental accelerated
        # solver with automatic LSODA fallback).
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
        if self.engine == 'fast' and self.fast_threads:
            from .. import fast_sgwb
            fast_sgwb.set_threads(self.fast_threads)
        self.stiffGW_model.SGWB_iter(engine=self.engine,
                                     fallback=self.fallback, **sgwb_kwargs)


        if self.stiffGW_model.SGWB_converge:
            state['f'] = self.stiffGW_model.f                                    # Output frequency in log10(f/Hz)
            state['omGW_stiff'] = self.stiffGW_model.log10OmegaGW                # log10(Omega_GW(f))
            state['hubble'] = self.stiffGW_model.derived_param['H_0']            # H_0 in units of s^-1
            state['kappa_s'] = self.stiffGW_model.derived_param['kappa_s']       # kappa_stiff(T_i) for AlterBBN
            state['kappa_r'] = self.stiffGW_model.kappa_r                        # kappa_rad(T_i) for AlterBBN, related to Delta_Neff

            if want_derived:
                yr = u.yr.to(u.s)
                log10f_yr = -math.log10(yr)
                if self.stiffGW_model.f[0] >= log10f_yr:
                    f_t = np.flip(state['f'])
                    Ogw_t = np.flip(state['omGW_stiff'])
                    spec_prim = interpolate.CubicSpline(f_t[f_t>-13], Ogw_t[f_t>-13])
                    omGW_stiff_fyr = spec_prim(log10f_yr)    # log10(Omega_GW(f_yr))
                else:
                    omGW_stiff_fyr = -100.

                state['derived'] = {'Delta_Neff_GW': self.stiffGW_model.DN_gw[-1],                  # Delta N_eff due to the primordial SGWB today
                                    'Delta_Neff_total': self.stiffGW_model.cosmo_param['DN_eff'],   # Total Delta N_eff after GW calculation
                                    'log10hc_prim_fyr': omGW_stiff_fyr/2 + math.log10(math.sqrt(1.5)*state['hubble']/math.pi)-log10f_yr,
                                                                                                    # log10(h_c(f_yr)) of the primordial SGWB
                                    'f_end': np.power(10., self.stiffGW_model.f[0]),                # Hz, UV cutoff frequency
                                   }
        else:
            #self.log.debug("SGWB calculation not converged, mostly due to total N_eff too large. Assigning 0 likelihood and going on.")
            return False


#    def get_A(self, normalization=1):
#        return self.current_state['A'] * normalization
