"""Physical constants and baseline cosmological values shared by the model.

中文：集中定义公式使用的常数和基准宇宙学量；这些模块级值属于数值模型，不是逐次求解配置。
需要调参时应通过模型参数或 `FastSolverConfig`，不要临时改写本模块常量。
"""

import math

import astropy.units as u
import numpy as np
from astropy import constants as const
from scipy import interpolate

from ._resources import package_path
from .functions import int_FD

#######   Shared constants; do not change during a solve.   #######
#######   共享常数；求解过程中不要修改。   #######

# Base physical constants
# 基本物理常数

kB = const.k_B.to(u.eV/u.K).value          # eV K^-1
hbar = const.hbar.to(u.eV*u.s).value       # eV s
c = const.c.to(u.cm/u.s).value             # cm s^-1
eV = (u.eV/(const.c**2)).to(u.g).value     # g*c^2
GN = const.G.value*1e3                     # cm^3 g^-1 s^-2
parsec = u.parsec.to(u.cm)                 # cm
yr = u.yr.to(u.s)                          # 儒略年的秒数。
ln10 = math.log(10)

# Derived physical constants
# 派生物理常数

rc = math.pi**2/(15*(hbar*c)**3)           # 辐射常数，单位 eV^-3 cm^-3。
m_Pl = math.sqrt(hbar*c**3/(GN*eV))*1e-9   # 普朗克质量，单位 GeV/c^2。
M_Pl = m_Pl/math.sqrt(8*math.pi)           # 约化普朗克质量，单位 GeV/c^2。
f_piv = 0.05 *c/(2*math.pi*1e6*parsec)     # CMB 枢轴波数 0.05 Mpc^-1 对应的频率，单位 s^-1。


#######   Cosmological parameters   #######
#######   宇宙学参数   #######

# Base parameters
# 基础参数

TCMB = 2.7255  # Planck 2018 的 CMB 温度，单位 K；固定该值以免导入 astropy.cosmology。
TCMB_GeV = 1e-9*kB*TCMB   # GeV

m_nu = 60                 # 单个有质量中微子本征态的质量，单位 meV。

Neff0 = 3.044
# Derived
# 派生参数

Omega_ph2 = rc * (kB*TCMB)**4 * (8*math.pi*GN*eV/3) *(0.1*parsec)**2         # 光子密度参数。
Omega_nh2 = Omega_ph2 * 7/8 * (4/11)**(4/3) * Neff0                          # 相对论标准模型中微子密度参数。
Omega_orh2 = Omega_ph2 + Omega_nh2                                           # 电子-正电子湮灭后的标准模型辐射密度参数。

Tnu = TCMB * (4/11)**(1/3)                # 今天所有中微子共用的有效热温度。
nu_today = m_nu*1e-3 / (kB*Tnu)
[rho_nu0, p_nu0] = int_FD(nu_today)
Omega_mnuh2 = Omega_nh2/3 * rho_nu0       # 有质量中微子本征态对 Omega_nu h^2 的贡献。


######      Thermal history from 20 keV to 10^6 GeV     ######
######      20 keV 至 10^6 GeV 的热史     ######

with package_path('stiffgwpy_fast', 'th.dat') as _th_path:
    thdata = np.loadtxt(_th_path)
T_th = thdata[:, 0]
T_max = T_th[-1]
N_th = thdata[:, 2]
N_max = N_th[-1]
N_fin = N_th[0]
rho_th = thdata[:, 5]
rhop_th = thdata[:, 6]

spl_rho = interpolate.CubicSpline(N_th, rho_th)
spl_rhop = interpolate.CubicSpline(N_th, rhop_th)

spl_T_N = interpolate.CubicSpline(np.log10(T_th), N_th)
N_10 = spl_T_N(-2)
a_10 = np.exp(-N_10)

#####       AlterBBN initial conditions      #########
#####       AlterBBN 初始条件      #########

T_i = 27. * 8.617330637338339e-5   # AlterBBN 演化的初始温度，对应 27×10^9 K。
N_i = spl_T_N(math.log10(T_i))
z_ratio = TCMB_GeV/T_i*np.exp(N_i)  # 温度比 T_CMB/(T_i*a_i)。


#####    Gravitational-wave reference frequencies    #####
#####    引力波参考频率    #####

f_yr = 1/yr     # PTA 使用的参考频率，单位 Hz。
f_LIGO = 25     # LIGO-Virgo 使用的参考频率，单位 Hz。


#####     Obsolete observational limits    #####
#####     已弃用的观测上限    #####

Neff_A = (2.86 + 0.57)*Neff0/3      # BBN 95% 上限，Aver 等（2015）与 Cooke 等（2014）。
Neff_I = 3.58 + 0.40                # BBN 95% 上限，Izotov 等（2014）与 Cooke 等（2014）。
# 以下 Neff 观测约束不参与当前计算，保留为历史参数。

Neff_l = 2.99 + 0.43                # CMB+BAO+Y_P 的 95% 上限（Aver 等，2015；含低红移声学尺度和 Silk 阻尼尺度）。
Neff_lnp = 2.97 + 0.58              # 同一分析但不含 Y_P 的上限。

T_np = 1.293              # 中子/质子比冻结温度，单位 MeV。
N_np = math.log(T_np*1e6/kB/TCMB)
T_D = 1/14                # 氘核合成温度，单位 MeV。
N_D = math.log(T_D*1e6/kB/TCMB)
N_BBN = [18.2, 23.3]      # 约覆盖 0.02–3 MeV，以包含完整的 BBN 过程。
