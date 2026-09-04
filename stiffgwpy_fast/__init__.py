# -*- coding: utf-8 -*-
"""stiffgwpy_fast: LCDM + stiff matter + primordial stochastic GW background (SGWB)."""
from .config import FastSolverConfig
from .LCDM_stiff_Neff import LCDM_SN
from .stiff_SGWB import LCDM_SG

__version__ = "0.2.1"
__all__ = ["LCDM_SG", "LCDM_SN", "FastSolverConfig"]
