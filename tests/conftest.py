# -*- coding: utf-8 -*-
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from stiffgwpy import LCDM_SG


@pytest.fixture
def model():
    """Baseline model (case A of the 12-case grid)."""
    return LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
