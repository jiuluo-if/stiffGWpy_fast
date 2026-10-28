"""LCDM plus stiff-matter cosmology and derived-parameter calculations.

中文：定义基础 `LCDM_SN` 模型、宇宙学输入参数及其派生量，供 SGWB solver 和 Cobaya adapter
共享。单位和参数含义也见根目录双语 README。
"""

import math
from collections.abc import Mapping
from os import PathLike
from pathlib import Path

import numpy as np
import yaml

from .functions import int_FD
from .global_param import (
    N_10,
    M_Pl,
    N_fin,
    N_i,
    N_max,
    Omega_mnuh2,
    Omega_nh2,
    Omega_orh2,
    Omega_ph2,
    T_i,
    T_max,
    TCMB_GeV,
    a_10,
    ln10,
    nu_today,
    parsec,
    rho_th,
    spl_rho,
    spl_rhop,
    spl_T_N,
)

_DEFAULT_COSMO_PARAMS = {
    'Omega_bh2': 0.0223828,
    'Omega_ch2': 0.1201075,
    'H0': 67.32117,
    'DN_eff': 0.0,
    'A_s': 2.100549e-9,
    'r': 0.0,
    'n_t': 0.0,
    'cr': 0,
    'T_re': 1e12,
    'DN_re': 0,
    'kappa10': 0.0,
}


class LCDM_SN:
    """
    Cosmological model: LCDM + stiff component + constant Delta N_eff

    The free/base parameters of the model are:
    'Omega_bh2', 'Omega_ch2', 'H0', 'DN_eff', 'A_s', 'r', 'n_t', 'cr', 'T_re', 'DN_re', 'kappa10'.
    - [H0] = km s^-1 Mpc^-1, [T_re] = GeV.
    - Set cr > 0 if the consistency relation is assumed, otherwise set cr <= 0 and provide (n_t, DN_re).
    - DN_re is the number of e-folds from the end of inflation to the end of reheating,
      assuming a^{-3} matter-like evolution.
    - kappa10 := rho_stiff/rho_photon at 10 MeV.

    There are three ways to instantiate a model with desired base parameters:
    1. input a yaml file,  2. input a dictionary,
    3. specify the parameters by keyword arguments.

    They will be stored in the 'obj_name.cosmo_param' attribute as a dictionary
    and can be modified later on.

    中文：模型参数保存在 `cosmo_param` 字典中。`derived_param` 根据当前参数计算共享的
    膨胀率、物质组分和张量初始条件；改变物理参数后应重置上层求解状态。

    """

    def __init__(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError("LCDM_SN accepts at most one positional configuration")

        if not args:
            params = _DEFAULT_COSMO_PARAMS.copy()
        elif isinstance(args[0], Mapping):
            params = dict(args[0])
        elif isinstance(args[0], (str, PathLike)):
            config_path = Path(args[0])
            if config_path.suffix.lower() not in ('.yml', '.yaml'):
                raise ValueError("configuration file must have a .yml or .yaml suffix")
            if not config_path.is_file():
                raise FileNotFoundError("configuration file not found: %s" % config_path)
            with config_path.open('r', encoding='utf-8') as stream:
                loaded = yaml.safe_load(stream)
            if not isinstance(loaded, Mapping):
                raise ValueError("configuration file must contain a YAML mapping")
            params = dict(loaded)
        else:
            raise TypeError("configuration must be a mapping or a YAML path")

        allowed = set(_DEFAULT_COSMO_PARAMS)
        unknown = (set(params) | set(kwargs)) - allowed
        if unknown:
            raise TypeError("unknown cosmological parameter(s): %s" % ', '.join(sorted(unknown)))
        params.update(kwargs)
        if args:
            missing = allowed - set(params)
            if missing:
                raise ValueError("configuration is missing required parameter(s): %s"
                                 % ', '.join(sorted(missing)))
        self.cosmo_param = params

        self.DN_eff_GW = 0          # 再加热结束前引力波贡献的相对论自由度。


    @property
    def derived_param(self):
        fast_cache = (bool(getattr(self, '_fast_derived_cache_enabled', False))
                      and not bool(getattr(
                          self, '_fast_derived_cache_disabled', False)))
        if fast_cache:
            # During one fast solve all cosmological inputs except DN_eff are
            # immutable; the outer loop changes only this self-consistency
            # variable.  The fast driver clears this cache at solve start.
            # 单次 fast 求解期间，除 DN_eff 外的宇宙学输入不变；外层迭代只更新该自洽变量。
            # 求解器会在每轮开始时清除此缓存。
            cache_key = self.cosmo_param['DN_eff']
            if (getattr(self, '_fast_derived_cache_key', None)
                    == cache_key):
                return self._fast_derived_cache
        derived_dict = {
            'h': self.cosmo_param['H0']/100,
            'H_0': self.cosmo_param['H0']/(10*parsec),     # 哈勃常数，单位 s^-1。

            'Omega_mh2': self.cosmo_param['Omega_bh2'] + self.cosmo_param['Omega_ch2'],                  # 重子与冷暗物质的总密度参数。
            'Omega_sh2': Omega_ph2 * self.cosmo_param['kappa10'] * (1e-2*a_10/TCMB_GeV)**4 * a_10**2,    # 当前时刻的 Omega_stiff*h^2。

            'A_t': self.cosmo_param['r'] * self.cosmo_param['A_s'],
            'nt': input_nt(self.cosmo_param),

            'kappa_s': self.cosmo_param['kappa10'] * (1e-2/T_i)**4 * math.exp(6*(N_i-N_10))              # AlterBBN 初始温度 T_i 处的 kappa_stiff。
        }

        if self.cosmo_param['T_re'] >= T_max:
            derived_dict['N_re'] = N_max + math.log(self.cosmo_param['T_re']/T_max)
            derived_dict['rho_re'] = rho_th[-1]
        else:
            derived_dict['N_re'] = spl_T_N(math.log10(self.cosmo_param['T_re']))
            derived_dict['rho_re'] = spl_rho(derived_dict['N_re'])
        self.rhorad_re = (derived_dict['rho_re'] + 7/8*(4/11)**(4/3)*self.cosmo_param['DN_eff']) * (TCMB_GeV*math.exp(N_10))**4
        # 再加热温度处的辐射密度系数，含额外辐射，单位为 GeV^4。
        self.rhostiff_re = self.cosmo_param['kappa10'] * 1e-8 * math.exp(2*(derived_dict['N_re']-N_10))
        # 再加热温度处的刚性物质密度系数，单位为 GeV^4。

        derived_dict['N_inf'] = None
        if self.cosmo_param['cr'] > 0:
            if self.cosmo_param['r'] > 0:
                derived_dict['V_inf'] = (1.5*derived_dict['A_t'])**.25 * math.pi**.5 * M_Pl
                # 单场慢滚暴胀的特征能标 V^(1/4)，单位 GeV。

                Delta_N = (N_10-derived_dict['N_re'])*4/3 + math.log(M_Pl)*4/3 + math.log(45/2*derived_dict['A_t']/(self.rhorad_re+self.rhostiff_re))/3
                # 假设暴胀结束至再加热结束期间按 a^{-3} 物质型演化，计算这段回望 e 折数。

                if Delta_N >= 0:
                    derived_dict['N_inf'] = derived_dict['N_re'] + Delta_N
                    self.cosmo_param['DN_re'] = Delta_N
                else:
                    print('Derived inflation-to-reheating interval is negative (Delta_N < 0). '
                          'Adjust relevant input parameters: r, T_re, kappa10.')
            else:
                print('r cannot be set to zero in a single-field, slow-roll inflation. Use positive r!')
        else:
            derived_dict['N_inf'] = derived_dict['N_re'] + self.cosmo_param['DN_re']

        if fast_cache:
            # ``cr > 0`` may update DN_re while constructing N_inf; retain the
            # post-computation key so the next unchanged lookup is reusable.
            # 构造 N_inf 时 `cr > 0` 可能更新 DN_re；缓存使用计算后的键值，供后续相同参数查询复用。
            cache_key = self.cosmo_param['DN_eff']
            self._fast_derived_cache_key = cache_key
            self._fast_derived_cache = derived_dict
        return derived_dict




    def gen_expansion(self):
        """
        Generate the expansion history of an extended Lambda-CDM model
        with an extra stiff matter and a constant Delta_Neff (which can
        mimic the effect of the primordial SGWB)

        N_re and N_inf are LOOKBACK numbers of e-folds
        at the end of reheating and inflation, respectively
        N_re = ln (T_re/T_CMB), N_BBN < N_re < N_inf

        Run this function only when self.derived_param['N_inf'] is not None.

        中文：生成含刚性物质和常数 ``Delta_Neff`` 的扩展 Lambda-CDM 膨胀历史；后者可模拟原初 SGWB 的影响。
        ``N_re`` 和 ``N_inf`` 分别是再加热结束和暴胀结束时的回溯 e-fold 数，且
        ``N_re = ln(T_re/T_CMB)``、``N_BBN < N_re < N_inf``。仅在
        ``self.derived_param['N_inf']`` 已有值时调用。

        """

        # 读取一次派生参数，供本次膨胀历史计算复用。
        derived = self.derived_param

        Omh2 = derived['Omega_mh2']
        Osh2 = derived['Omega_sh2']

        Oerh2 = Omega_ph2 * 7/8 * (4/11)**(4/3) * self.cosmo_param['DN_eff']         # 额外辐射的 Omega*h^2。
        Otrh2 = Omega_orh2 + Oerh2                                                   # 电子-正电子湮灭后的总辐射密度参数。
        Otreh2 = Omega_ph2 * rho_th[-1] + Oerh2                                 # 标准模型相变前的总辐射密度参数。

        OLh2 = derived['h']**2 - Omh2 - Omega_mnuh2 - Omega_nh2*2/3 - Omega_ph2 - Oerh2 - Osh2     # 暗能量密度参数 Omega_Lambda*h^2。


        # Construct output arrays / 构造输出数组。

        len_inf = math.floor(derived['N_inf']*100)+1
        Nv = np.arange(0, len_inf)*.01
        if Nv[-1] != derived['N_inf']:
            Nv = np.append(Nv, derived['N_inf'])
        # Nv is equivalent to ln(a), its present-day value is set at Nv[-1] = N_inf.
        # Nv 对应 ln(a)，并将今天的端点设为 Nv[-1]=N_inf。
        Sv = np.zeros_like(Nv)
        fv = np.zeros_like(Nv)  # fv 与视界频率 f_H=aH/(2π) 成比例。


        # Main loop: calculate the expansion history / 主循环：计算膨胀历史。

        index_re = int(np.argmin(np.abs(Nv - (derived['N_inf'] - derived['N_re']))))
        # 再加热边界落在网格上最近的节点；连续边界位置由 N_inf-N_re 给出。

        for i in range(index_re, len(Nv), 1):
            eN = math.exp(Nv[-1]-Nv[i])
            e3N = math.exp(3.0*(Nv[-1]-Nv[i]))  # 分别对应 a^-1 和 a^-3。
            nu = nu_today / eN          # 中微子质量与热能之比 (m_nu*c^2)/(kB*T_nu)。
            if (nu > 100):              # 有质量中微子处于深度非相对论区。
                H2 = Omh2 + Omega_mnuh2 + (Omega_ph2+2/3*Omega_nh2+Oerh2)*eN + Osh2*e3N + OLh2/e3N    # h^2*a^3 项。
                Sv[i] = (Omh2 + Omega_mnuh2 + 4/3*(Omega_ph2+2/3*Omega_nh2+Oerh2)*eN + 2*Osh2*e3N)/H2
            elif (nu <= 100) and (nu >= 0.1):
                [rho_nu, p_nu] = int_FD(nu)
                H2 = Omh2 + (Omega_ph2+(2/3+rho_nu/3)*Omega_nh2+Oerh2)*eN + Osh2*e3N + OLh2/e3N
                Sv[i] = (Omh2 + 4/3*(Omega_ph2+2/3*Omega_nh2+Oerh2)*eN + (rho_nu+p_nu)*Omega_nh2/3*eN + 2*Osh2*e3N)/H2
            elif (nu < 0.1) and (Nv[i] > Nv[-1]-N_fin):     # 有质量中微子处于深度相对论区。
                H2 = Omh2 + Otrh2*eN + Osh2*e3N + OLh2/e3N
                Sv[i] = (Omh2 + 4/3*Otrh2*eN + 2*Osh2*e3N)/H2
            elif (Nv[i] <= Nv[-1]-N_fin) and (Nv[i] >= Nv[-1]-N_max):   # 在 20 keV 至 10^6 GeV 区间使用热史查找表。
                rho_i = spl_rho(Nv[-1]-Nv[i])
                rhop_i = spl_rhop(Nv[-1]-Nv[i])
                H2 = Omh2 + (Omega_ph2*rho_i + Oerh2)*eN + Osh2*e3N + OLh2/e3N
                Sv[i] = (Omh2 + (Omega_ph2*rhop_i + 4/3*Oerh2)*eN + 2*Osh2*e3N)/H2
            else:
                H2 = Omh2 + Otreh2*eN + Osh2*e3N + OLh2/e3N
                Sv[i] = (Omh2 + 4/3*Otreh2*eN + 2*Osh2*e3N)/H2

            fv[i] = -0.5*Nv[i] + 0.5*math.log(H2)    # 计算 N 与 ln(H) 的组合量。

        Sv[0:index_re] = 1
        fv[0:index_re] = fv[index_re] - 0.5*(Nv[0:index_re] - Nv[index_re])
        # reheating is assumed to be MD since the end of inflation
        # 暴胀结束至再加热结束按 MD（物质主导）分支延拓，故 sigma=1。

        # Convert fv into physical units: here \tilde f_H = ln (f_H/Hz)
        # 将 fv 转换为 ln(f_H/Hz)，最终再除以 ln(10) 输出 log10(f_H/Hz)。
        f0 = fv[-1]
        Delta_f = math.log(2*math.pi/derived['H_0'])
        fv = fv - f0 - Delta_f

        # Output the expansion history and frequencies / 写出膨胀历史与频率数组。

        self.Nv = Nv
        self.N = Nv - Nv[-1]
        self.sigma = Sv
        self.f_hor = fv / ln10       # 以 Hz 为单位的视界频率 log10(f_H/Hz)。
        self.f_re = self.f_hor[index_re]


def input_nt(params):
    """
    Input the tensor spectral index.
    Use the inflationary consistency relation or the provided n_t.
    中文：读取张量谱指数 ``n_t``；可按暴胀一致性关系从模型参数推导，或使用用户显式提供的值。
    """
    if params['cr'] > 0:
        return -params['r']/8        # 暴胀一致性关系。
    else:
        return params['n_t']
