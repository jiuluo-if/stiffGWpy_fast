from types import SimpleNamespace

import numpy as np
import pytest

from stiffgwpy_fast import LCDM_SG, reference
from stiffgwpy_fast import functions as legacy


def test_lsoda_wrapper_rejects_failed_partial_solution(monkeypatch):
    monkeypatch.setattr(
        legacy,
        "solve_ivp",
        lambda *args, **kwargs: SimpleNamespace(
            success=False, status=-1, message="forced incomplete integration"
        ),
    )

    with pytest.raises(RuntimeError, match="LSODA.*incomplete integration"):
        legacy.solve_SGWB(np.arange(4.0), np.ones(4), 0, -1.0)


def test_reference_rejects_failed_partial_solution(monkeypatch):
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    failed_result = SimpleNamespace(
        success=False,
        status=-1,
        message="forced incomplete integration",
        t_events=[np.array([])],
        nfev=3,
        y=np.array([[-10.0, -9.0], [0.0, 1e-20], [1e-5, 2e-5]]),
    )
    monkeypatch.setattr(reference, "solve_ivp", lambda *args, **kwargs: failed_result)

    with pytest.raises(RuntimeError, match="reference.*incomplete integration"):
        reference.solve_reference_mode(model, -6.0, 0.0)
