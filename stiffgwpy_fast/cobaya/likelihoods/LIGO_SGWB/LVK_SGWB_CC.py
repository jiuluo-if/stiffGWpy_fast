"""Cobaya likelihood for LVK stochastic-background cross-correlation data.

中文：读取 YAML 指定的交叉相关数据，使用 theory 提供的频率和 SGWB 谱计算高斯似然。
"""

import numpy as np
from cobaya.likelihood import Likelihood
from scipy import interpolate


class LVK_SGWB_CC(Likelihood):

    def initialize(self):
        """Load the selected cross-correlation table and prepare fixed columns.

        中文：读取配置中的 `CC_file`，并缓存不随理论参数变化的对数频率、观测值和不确定度列。
        该选项默认指向随包提供的 LVK 数据表。
        """
        self.data = np.loadtxt(self.CC_file)
        self.f_LVK = np.log10(self.data[:, 0])
        self.Cf_LVK = self.data[:, 1]
        self.sigma_LVK = self.data[:, 2]

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
        f_theory = np.flip(self.provider.get_result('f'))
        Ogw_theory = np.flip(self.provider.get_result('omGW_stiff'))
        return self.log_likelihood(f_theory, Ogw_theory, **params_values)


    def log_likelihood(self, f_theory, Ogw_theory, **data_params):
        """Interpolate the model on supported bins and return ``-chi2/2``.

        ``f_theory`` must be increasing in ``log10(Hz)`` order.

        中文：将模型谱转换为线性 ``Omega_GW``，按数据不确定度归一化残差并求和；
        不在模型频率支持内的 bin 保持零模型值。
        """
        f_LVK = self.f_LVK
        Cf_LVK = self.Cf_LVK
        sigma_LVK = self.sigma_LVK

        Ogw_Model = np.zeros_like(Cf_LVK)
        # Calculate theoretical Omega_GW only within the LIGO frequency range (<= f_end).
        # 仅对理论频谱覆盖范围内的 LIGO 频点计算 Omega_GW。
        if f_theory[-1]>=f_LVK[0]:
            theory_support = f_theory >= -5
            f_t = f_theory[theory_support]
            Ogw_t = Ogw_theory[theory_support]
            spec = interpolate.interp1d(f_t, Ogw_t, kind='cubic')

            cond = (f_LVK >= f_t[0]) & (f_LVK <= f_t[-1])
            Ogw_Model[cond] = np.power(10., spec(f_LVK[cond]))

        chi2 = np.sum(
            np.square((Cf_LVK - Ogw_Model) / sigma_LVK), dtype=np.float64)

        return -chi2 / 2
