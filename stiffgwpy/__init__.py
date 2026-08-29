# -*- coding: utf-8 -*-
"""stiffgwpy: LCDM + stiff matter + primordial stochastic GW background (SGWB).

Exposes the canonical model class ``LCDM_SG`` (equivalent to the legacy
``from stiff_SGWB import LCDM_SG``) and the accelerated solver
``stiffgwpy.fast_sgwb.SGWB_iter_fast``.
"""
from .stiff_SGWB import LCDM_SG
from .LCDM_stiff_Neff import LCDM_SN

__version__ = "0.2.0"
__all__ = ["LCDM_SG", "LCDM_SN"]
