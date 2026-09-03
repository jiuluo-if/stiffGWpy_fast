# -*- coding: utf-8 -*-
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from stiffgwpy import LCDM_SG
from stiffgwpy import fast_sgwb as FS


@pytest.fixture
def model():
    """Baseline model (case A of the 12-case grid)."""
    return LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)


@pytest.fixture
def fast_settings():
    """Snapshot and restore process-global fast-solver settings."""
    saved = FS.get_settings()
    freq_grid = FS._FREQ_GRID
    yield saved
    FS.set_threads(saved['threads'])
    FS.set_col_step(saved['col_step'])
    FS.set_h(saved['h'])
    FS.set_z_tail(saved['z_tail'])
    FS.set_phase_max(saved['phase_max'])
    FS.set_freq_grid(freq_grid)
