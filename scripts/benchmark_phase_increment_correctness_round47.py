"""Named/edge/Sobol full-outer correctness gate for the Round 47 twin."""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts._resource_budget import apply_environment  # noqa: E402

apply_environment()

import numpy as np  # noqa: E402

from scripts.benchmark_phase_increment_kernel_twin import (  # noqa: E402
    solve_kernel_phase_increment,
)
from scripts.benchmark_phase_recurrence import CASES  # noqa: E402
from scripts.benchmark_same_grid_reference import CASES as ORACLE_CASES  # noqa: E402
from scripts.validate_edges_vs_reference import EDGE_POINTS  # noqa: E402
from stiffgwpy_fast import fast_sgwb as FS  # noqa: E402
from stiffgwpy_fast.stiff_SGWB import LCDM_SG  # noqa: E402

CASES_TO_RUN = {
    **CASES,
    **{name: EDGE_POINTS[name] for name in (
        'edge_r_lo', 'edge_tre_hi', 'edge_kap_lo', 'edge_nt_blue')},
    **{name: ORACLE_CASES[name] for name in (
        'sobol_000', 'sobol_006', 'sobol_010', 'sobol_015')},
}


def digest(value):
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def run_one(name, kwargs, candidate):
    original = FS.solve_kernel
    if candidate:
        FS.solve_kernel = lambda *args, **kw: solve_kernel_phase_increment(
            *args, 32)
    try:
        model = LCDM_SG(**kwargs)
        result = FS.SGWB_iter_fast(
            model, kink_split=True, freq_grid='goal',
            frequency_quadrature='pchip')
        return model, result
    finally:
        FS.solve_kernel = original


def compare(name, kwargs):
    base, base_result = run_one(name, kwargs, False)
    cand, cand_result = run_one(name, kwargs, True)
    if base_result is None or cand_result is None:
        return {
            'name': name,
            'baseline_result_none': base_result is None,
            'candidate_result_none': cand_result is None,
            'baseline_converged': bool(getattr(base, 'SGWB_converge', False)),
            'candidate_converged': bool(getattr(cand, 'SGWB_converge', False)),
            'baseline_failure': getattr(base, 'fast_failure_reason', None),
            'candidate_failure': getattr(cand, 'fast_failure_reason', None),
        }
    base_spec = np.asarray(base.log10OmegaGW, dtype=np.float64)
    cand_spec = np.asarray(cand.log10OmegaGW, dtype=np.float64)
    base_dn = float(np.asarray(base.DN_gw)[-1])
    cand_dn = float(np.asarray(cand.DN_gw)[-1])
    return {
        'name': name,
        'baseline_result_none': base_result is None,
        'candidate_result_none': cand_result is None,
        'baseline_converged': bool(getattr(base, 'SGWB_converge', False)),
        'candidate_converged': bool(getattr(cand, 'SGWB_converge', False)),
        'baseline_failure': getattr(base, 'fast_failure_reason', None),
        'candidate_failure': getattr(cand, 'fast_failure_reason', None),
        'baseline_iterations': int(len(np.asarray(base.DN_gw))),
        'candidate_iterations': int(len(np.asarray(cand.DN_gw))),
        'spectrum_max_abs_dex': float(np.max(np.abs(cand_spec - base_spec))),
        'spectrum_max_rel_linear': float(np.max(
            np.abs(10.0 ** cand_spec - 10.0 ** base_spec) /
            np.maximum(np.abs(10.0 ** base_spec), 1e-300))),
        'DN_gw_relative': abs(cand_dn - base_dn) / max(abs(base_dn), 1e-300),
        'finite_candidate': bool(
            np.isfinite(cand_spec).all() and np.isfinite(cand_dn) and
            np.isfinite(np.asarray(cand.g2)).all() and
            np.isfinite(np.asarray(cand.w2)).all()),
        'digests': {
            'f_equal': digest(base.f) == digest(cand.f),
            'spectrum_equal': digest(base.log10OmegaGW) == digest(cand.log10OmegaGW),
            'DN_gw_equal': digest(base.DN_gw) == digest(cand.DN_gw),
            'g2_equal': digest(base.g2) == digest(cand.g2),
            'w2_equal': digest(base.w2) == digest(cand.w2),
        },
    }


def main():
    FS.apply_accuracy_mode('fast')
    FS.set_threads(2)
    records = [compare(name, kwargs) for name, kwargs in CASES_TO_RUN.items()]
    payload = {
        'experiment': 'round47_phase_increment_full_outer_correctness',
        'production_path_changed': False,
        'commit': os.popen('git rev-parse HEAD').read().strip(),
        'settings': {'threads': 2, 'kink_split': True,
                     'freq_grid': 'goal', 'frequency_quadrature': 'pchip',
                     'reanchor': 32},
        'records': records,
    }
    out = ROOT / 'docs/phase_increment_correctness_round47_20260923.json'
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n',
                   encoding='utf-8')
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
