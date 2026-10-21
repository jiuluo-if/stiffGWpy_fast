"""Cobaya IPTA likelihood using configured sample and frequency resources.

中文：根据 PTA 频点上的原初 SGWB 谱和 SMBHB 背景，使用每个频点的样本分布计算 IPTA 似然。
"""

import math

import astropy.units as u
import numpy as np
from cobaya.likelihood import Likelihood
from scipy import interpolate
from scipy.stats import gaussian_kde as kde


class IPTA(Likelihood):
    
    def initialize(self):
        """Load IPTA samples and prefit one KDE for each configured frequency.

        中文：读取 `IPTA.yaml` 指定的样本和频率，并预拟合每个频点的 KDE。
        """
        self.data = np.loadtxt(self.sample_file)
        self.freqs = np.loadtxt(self.freq_file)
        # Per-bin bounds and KDEs depend only on these fixed samples.
        # 每个频点的边界和 KDE 只依赖固定样本，不依赖 theory 输出。
        self._log10rho_max = np.max(self.data, axis=0)
        self._log10rho_min = np.min(self.data, axis=0)
        self._kde_by_frequency = tuple(
            kde(self.data[:, i]) for i in range(self.Nfreqs))
    
    def close(self):
        pass
    
            
    def get_requirements(self):
        """Request frequency, SGWB spectrum, and Hubble value from the theory.

        中文：声明此似然需要 theory 提供 `f`、`omGW_stiff` 和 `hubble`。
        """
        return {'f': None, 'omGW_stiff': None, 'hubble': None}
    
    
    def logp(self, _derived=None, **params_values):
        """Return the IPTA log-likelihood for the current theory and nuisance parameters.

        中文：取回频谱和 Hubble 值，转为递增频率顺序后传给独立的似然计算函数。
        """
        f_theory = np.flip(self.provider.get_result('f'))  # log10(f/Hz)
        Ogw_theory = np.flip(self.provider.get_result('omGW_stiff'))  # log10(Omega_GW(f))
        H_0 = self.provider.get_result('hubble')  # s^{-1}
        
        return self.log_likelihood(f_theory, Ogw_theory, H_0, **params_values)

    
    def log_likelihood(self, f_theory, Ogw_theory, H_0, **data_params):
        """Evaluate the per-frequency PTA likelihood; ``f_theory`` must be increasing.

        The prediction combines the primordial spectrum with the SMBHB power-law
        nuisance component. 中文：在频率支持范围内构造总残差功率，并对各频点 KDE 的 log-PDF 求和。
        """
        yr = u.yr.to(u.s)                # s, one Julian year
        T_base = 1/self.freqs[0]         # s, baseline of the PTA data
        
        Nf = self.Nfreqs                 # Number of work frequencies, i = 1 - Nf
        work_freqs = self.freqs[:Nf]     # work frequencies in Hz
        # Limiting work frequencies can retain stiff models that higher-frequency data would exclude (i=Nf+1-Ntot).
        # 限制工作频点可能保留本会被高频数据排除的 stiff 模型（i=Nf+1-Ntot）。

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
        log10rho_max = self._log10rho_max
        log10rho_min = self._log10rho_min
        kde_by_frequency = self._kde_by_frequency
        
        logL = 0
        for i in range(Nf):
            if log10rho_Model[i]>log10rho_max[i]:
                logL += kde_by_frequency[i].logpdf(log10rho_max[i])[0]
            elif log10rho_Model[i]<log10rho_min[i]:
                logL += kde_by_frequency[i].logpdf(log10rho_min[i])[0]
            else:
                logL += kde_by_frequency[i].logpdf(log10rho_Model[i])[0]
        
        return logL
