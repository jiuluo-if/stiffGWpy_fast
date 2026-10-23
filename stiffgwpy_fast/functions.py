"""Massive-neutrino thermodynamics and legacy LSODA tensor-mode helpers.

The high-level `lsoda` engine uses these routines; the formal fast engine has
its own Numba stepping implementation in `fast_sgwb.py`.

中文：这里提供费米–狄拉克中微子热力学积分和旧版 LSODA 单频率张量积分辅助函数。正式 fast
档位使用 `fast_sgwb.py` 中独立的 Numba 步进路径。
"""

import math

import numpy as np
from scipy import interpolate
from scipy.integrate import quad, solve_ivp


def int_FD(y):
    """
    Evaluate energy density and pressure for a massive neutrino in Fermi-Dirac dist.
    x = pc/(kB*T), y = mc^2/(kB*T) = mc^2/(kB*T0)*a
    中文：计算有质量中微子的费米-狄拉克能量密度和压强积分，并按无质量费米子的能量密度归一化。
    其中 ``x = pc/(kB*T)``，``y = mc^2/(kB*T) = mc^2/(kB*T0)*a``。

    """

    I_rho = quad(integrand_rho, 0, 30, args=(y))
    I_p = quad(integrand_p, 0, 30, args=(y))

    coeff_massless = 7*math.pi**4/120

    rho = I_rho[0]/coeff_massless
    p = I_p[0]/coeff_massless

    return np.array((rho, p))


def integrand_rho(x,y):
    return x**2 * math.sqrt(x**2+y**2)/(math.exp(x)+1)

def integrand_p(x,y):
    return x**4 / (3 * math.sqrt(x**2+y**2) * (math.exp(x)+1))



def _make_subhorizon_event(z_tail):
    """Factory for the deep-subhorizon terminal event (default z_tail = 5).

    中文：创建深次视界终止事件；当 ``z`` 向上穿过 ``z_tail`` 时停止积分，默认阈值为 5。
    """
    def subhorizon(N, state, spline):
        return state[0] - z_tail
    subhorizon.terminal = True
    subhorizon.direction = 1
    return subhorizon


def solve_SGWB(Nv, Sv, j0, z0, z_tail=5.0, rtol=1e-6, atol=None):
    """Integrate one tensor mode with LSODA on a precomputed background grid.

    ``Nv`` and ``Sv`` hold the expansion grid and sampled ``sigma(N)``;
    ``j0`` and ``z0`` define the mode's initial grid point and
    ``z = ln(k/aH)``. A terminal event stops integration as ``z`` crosses
    ``z_tail`` upward. ``rtol`` and ``atol`` are passed to ``solve_ivp``.
    The returned result includes a dense solution in ``result.sol``.

    中文：在预先生成的背景网格上用 LSODA 积分单个张量模态。``Nv`` 和 ``Sv`` 分别是膨胀网格与
    采样的 ``sigma(N)``；``j0`` 和 ``z0`` 指定模态的起始网格点及 ``z = ln(k/aH)`` 初值。
    当 ``z`` 向上穿过 ``z_tail`` 时终止；``rtol`` 和 ``atol`` 传给 ``solve_ivp``，返回结果的
    ``result.sol`` 提供稠密解。
    """

    N_span = (Nv[j0], Nv[-1])
    ini_state = [z0, 0, math.exp(z0)]
    param = (interpolate.InterpolatedUnivariateSpline(Nv, Sv),)
    if atol is None:
        atol = [1e-10, 1e-20, 1e-20]

    result = solve_ivp(tensor, N_span, ini_state,
                       method='LSODA',
                       dense_output=True,
                       events=[_make_subhorizon_event(z_tail)],
                       args=param,
                       rtol=rtol, atol=atol,
                       jac=jacobian,
                      )

    return result



def tensor(N, state, spline):
    r"""
    Dynamical system for tensor modes:

    中文：以 ``N`` 为自变量，将张量模态方程写成 ``(z, x, y)`` 的一阶系统；``sigma`` 由背景样条提供。

    z = ln(2*pi*f/aH)
    x = \dot T_h / H
    y = (2*pi*f/aH) * T_h
    sigma = -2 \dot H / 3H^2

    z' = 1.5*sigma - 1
    x' = -3*x + 1.5*sigma*x - exp(z)*y
    y' = -y + 1.5*sigma*y + exp(z)*x

    """

    z, x, y = state
    sigma = float(spline(N))

    dz = 1.5*sigma - 1
    dx = -3*x + 1.5*sigma*x - math.exp(z)*y
    dy = -y + 1.5*sigma*y + math.exp(z)*x

    return [dz, dx, dy]


def jacobian(N, state, spline):

    z, x, y = state
    sigma = float(spline(N))

    jac = np.array(
        (
            (0, 0, 0),
            (-math.exp(z)*y, -3+1.5*sigma, -math.exp(z)),
            (math.exp(z)*x, math.exp(z), -1+1.5*sigma),
        ),
        dtype = np.float64,
    )

    return jac


def subhorizon(N, state, spline):         # 深亚视界阈值：k/aH = exp(5)，约为 150
    return state[0]-5
