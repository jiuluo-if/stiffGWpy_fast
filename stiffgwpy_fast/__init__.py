# -*- coding: utf-8 -*-
"""stiffgwpy_fast: LCDM + stiff matter + primordial stochastic GW background (SGWB).

中文：提供 `LCDM_SG` 模型、fast/reference/LSODA 引擎入口和逐次调用的快求解器配置。
"""
from .config import FastSolverConfig
from .LCDM_stiff_Neff import LCDM_SN
from .stiff_SGWB import LCDM_SG

__version__ = "0.2.1"
__all__ = ["LCDM_SG", "LCDM_SN", "FastSolverConfig"]
