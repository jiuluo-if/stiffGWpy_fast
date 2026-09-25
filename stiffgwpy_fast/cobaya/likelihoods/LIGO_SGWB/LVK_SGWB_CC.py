"""Cobaya likelihood for LVK stochastic-background cross-correlation data.

中文：读取 YAML 指定的交叉相关数据，使用 theory 提供的频率和 SGWB 谱计算高斯似然。
"""

from cobaya.likelihood import Likelihood
from cobaya.log import LoggedError, get_logger
import numpy as np
import os
from scipy import interpolate


class LVK_SGWB_CC(Likelihood):
    
    def initialize(self):
        """Load the cross-correlation table selected by ``CC_file``.

        中文：读取配置中的 `CC_file`；该选项默认指向随包提供的 LVK 数据表。
        """
        self.data = np.loadtxt(self.CC_file)
    
    def close(self):
        pass
        
    def get_requirements(self):
        """Request the theory frequency grid and SGWB spectrum.

        中文：声明似然需要 theory 提供的频率 `f` 和 `omGW_stiff`。
        """
        return {'f': None, 'omGW_stiff': None,}

    
    def logp(self, _derived=None, **params_values):
        """Return the LVK log-likelihood for the current theory spectrum.

        Cobaya supplies nuisance parameters through ``params_values``.
        中文：取回 theory 频率和谱，并转为递增频率顺序后交给 `log_likelihood`。
        """
        f_theory = self.provider.get_result('f'); f_theory = np.flip(f_theory)
        Ogw_theory = self.provider.get_result('omGW_stiff'); Ogw_theory = np.flip(Ogw_theory)
        
        #if _derived is not None:
        #    _derived['N_eff'] = self.provider.get_param('Delta_Neff_GW') + 3.044
        
        return self.log_likelihood(f_theory, Ogw_theory, **params_values)

    
    def log_likelihood(self, f_theory, Ogw_theory, **data_params):
        """Interpolate the model on supported bins and return ``-chi2/2``.

        ``f_theory`` must be increasing in ``log10(Hz)`` order. 中文：将模型谱转换为线性
        ``Omega_GW``，按数据不确定度归一化残差并求和；不在模型频率支持内的 bin 保持零模型值。
        """
        f_LVK = np.log10(self.data[:,0])
        Cf_LVK = self.data[:,1]
        sigma_LVK = self.data[:,2]
        
        Ogw_Model = np.zeros_like(Cf_LVK)
        if f_theory[-1]>=f_LVK[0]:   # Calculate theoretical Omega_GW ONLY for LIGO frequency bins in its range, i.e., <= f_end
            f_t = f_theory[(f_theory >= -5)]; Ogw_t = Ogw_theory[(f_theory >= -5)]
            #print(f_t, Ogw_t)
            spec = interpolate.interp1d(f_t, Ogw_t, kind='cubic')
            
            cond = (f_LVK<=f_theory[-1])
            Ogw_Model[cond] = np.power(10., spec(f_LVK[cond]))  
            
        #print(Ogw_Model)
             
        chi2_array = np.square(np.divide((Cf_LVK-Ogw_Model), sigma_LVK))
        chi2 = sum(chi2_array)
        
        return -chi2 / 2
