"""Cobaya NANOGrav likelihood using a packaged density grid.

中文：读取随包提供的频率、``log10(rho)`` 网格和似然密度数组，计算 NANOGrav PTA 似然。
"""

import math

import astropy.units as u
import numpy as np
from cobaya.likelihood import Likelihood
from scipy import interpolate


class NANOGrav(Likelihood):
    
    def initialize(self):
        """Load the density grid and prefit its fixed per-frequency splines.

        中文：读取 `NANOGrav.yaml` 指定的资源，并预拟合各频点的似然密度样条。
        """
        self.data = np.load(self.kde_file)                   # log-likelihood for each frequency bin
        self.log10rhogrid = np.load(self.log10rhogrid_file)  # log10(delay/s) = log10(sqrt(rho/s^2))
        self.freqs = np.load(self.freq_file)                 # Hz, frequency bins of the PTA experiment
        # These splines depend only on packaged likelihood data, not on theory parameters.
        # 这些样条只依赖随包分发的似然数据，与 theory 参数无关。
        self._density_splines = tuple(
            interpolate.CubicSpline(self.log10rhogrid, self.data[0, i])
            for i in range(self.Nfreqs))
        self._log10rhogrid_min = float(np.min(self.log10rhogrid))
        self._log10rhogrid_max = float(np.max(self.log10rhogrid))
    
    def close(self):
        pass
    
            
    def get_requirements(self):
        """Request frequency, SGWB spectrum, and Hubble value from the theory.

        中文：声明此似然需要 theory 提供 `f`、`omGW_stiff` 和 `hubble`。
        """
        return {'f': None, 'omGW_stiff': None, 'hubble': None}
    
    
    def logp(self, _derived=None, **params_values):
        """Return the likelihood for the current theory and nuisance parameters.

        中文：取回 theory 频率、SGWB 谱和 Hubble 值，转为递增频率顺序后计算 NANOGrav 似然。
        """
        f_theory = np.flip(self.provider.get_result('f'))  # log10(f/Hz)
        Ogw_theory = np.flip(self.provider.get_result('omGW_stiff'))  # log10(Omega_GW(f))
        H_0 = self.provider.get_result('hubble')  # s^{-1}
        
        return self.log_likelihood(f_theory, Ogw_theory, H_0, **params_values)

    
    def log_likelihood(self, f_theory, Ogw_theory, H_0, **data_params):
        """Evaluate the NANOGrav density grid; ``f_theory`` must be increasing.

        中文：组合原初 SGWB 与 SMBHB 前景；超出数据网格的模型值按边界密度计算。
        """
        yr = u.yr.to(u.s)               # s, one Julian year
        T_base = 1/self.freqs[0]        # s, baseline of the PTA data
        
        Nf = self.Nfreqs                # Number of work frequencies, i = 1 - Nf
        work_freqs = self.freqs[:Nf]    # work frequencies in Hz
        # Limiting frequency bins can retain stiff-amplified models that higher-frequency data would exclude.
        # 限制工作频点可能保留本会被更高频数据排除的 stiff 放大模型（i=Nf+1-Ntot）。
        
        # Primordial SGWB
        rho_prim = np.zeros_like(work_freqs)
        if f_theory[-1] >= math.log10(work_freqs[0]):
            cond = (f_theory > -13) & (f_theory < -5)
            f_t = f_theory[cond]; Ogw_t = Ogw_theory[cond]
            spec_prim = interpolate.CubicSpline(f_t, Ogw_t)

            cond = (np.log10(work_freqs)<=f_t[-1])
            # Calculate theoretical Omega_GW ONLY for PTA frequency bins in its range, i.e., <= f_end
            Ogw_prim = spec_prim(np.log10(work_freqs[cond]))        
            rho_prim[cond] = np.divide(np.power(10., Ogw_prim) * H_0**2, 8*np.pi**4 * work_freqs[cond]**5 * T_base)  # s^2

        # SGWB from SMBHBs
        log10hc_BH = data_params['A_BBH']      # log10(h_c) at f_yr
        gamma_BH = data_params['gamma_BBH']
        rho_BH = np.power(10., log10hc_BH*2) * np.power(work_freqs*yr, -gamma_BH) * yr**3 / (12*np.pi**2 * T_base)   # s^2
        
        log10rho_Model = np.log10(np.sqrt(rho_prim + rho_BH))    # log10(delay/s) = log10(sqrt(rho/s^2))

        logL = 0
        for i in range(Nf):
            density_spline = self._density_splines[i]
            if log10rho_Model[i]>self._log10rhogrid_max:
                logL += density_spline(self._log10rhogrid_max)
            elif log10rho_Model[i]<self._log10rhogrid_min:
                logL += density_spline(self._log10rhogrid_min)
            else:
                logL += density_spline(log10rho_Model[i])

        return logL
