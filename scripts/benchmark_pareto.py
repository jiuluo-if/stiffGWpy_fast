# -*- coding: utf-8 -*-
"""benchmark_pareto.py -- runtime vs physical-error Pareto for the SGWB engines.

Measures, on the default point, each engine's wall-clock time and its integrated
``Delta N_eff`` relative error against the *independent continuous-sigma
reference* (``stiffgwpy.reference``) rather than the old LSODA path.  This is the
physics-first accuracy-vs-speed trade-off requested by the audit, not a
fast-vs-LSODA speedup report.

Engines:
  * fast (plain grid, fastest; the historical -0.2% grid path),
  * fast transition-refine (production; kink-aware, ~-0.04%),
  * lsoda (original adaptive LSODA path),
  * reference (continuous-sigma DOP853; the physics anchor, error 0).

The reference and LSODA runs are slow (tens of seconds to minutes each), so this
benchmark is meant for a single point / certification, not for CI.
"""

import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy import fast_sgwb as FS
from stiffgwpy.stiff_SGWB import LCDM_SG

KW = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)


def env_meta():
    meta = {'numpy': np.__version__}
    for mod in ('numba', 'scipy'):
        try:
            meta[mod] = __import__(mod).__version__
        except Exception:
            pass
    try:
        out = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                             text=True, timeout=10)
        if out.returncode == 0:
            meta['commit'] = out.stdout.strip()
    except Exception:
        pass
    return meta


def meas_fast(tr):
    """Measure a fast config; returns (runtime_s, DN_eff)."""
    FS.apply_accuracy_mode('production')
    FS.set_z_tail(5.0)
    m = LCDM_SG(**KW)
    t0 = time.perf_counter()
    m = FS.SGWB_iter_fast(m, tol=1e-7, transition_refine=tr)
    dt = time.perf_counter() - t0
    return dt, float(m.cosmo_param['DN_eff'])


def meas_lsoda():
    m = LCDM_SG(**KW)
    t0 = time.perf_counter()
    m.SGWB_iter(engine='lsoda')
    dt = time.perf_counter() - t0
    return dt, float(m.cosmo_param['DN_eff'])


def meas_reference():
    from stiffgwpy import reference as REF
    m = LCDM_SG(**KW)
    t0 = time.perf_counter()
    ref = REF.run_reference(m, freq_res=1.0, z_tail=5.0, rtol=1e-11,
                            self_consistent=True)
    dt = time.perf_counter() - t0
    return dt, float(ref['DN_gw'])


def main():
    out_dir = os.path.join('docs', 'reference')
    os.makedirs(out_dir, exist_ok=True)
    records = []

    def rec(r):
        records.append(r)
        print(json.dumps(r, ensure_ascii=False), flush=True)

    ref_t, ref_dn = meas_reference()
    rec({'engine': 'reference', 'runtime_s': ref_t, 'DN_eff': ref_dn,
         'rel_err': 0.0})
    for name, tr in (('fast_transition', True), ('fast_grid', False)):
        t, dn = meas_fast(tr)
        rec({'engine': name, 'runtime_s': t, 'DN_eff': dn,
             'rel_err': (dn - ref_dn) / ref_dn})
    t, dn = meas_lsoda()
    rec({'engine': 'lsoda', 'runtime_s': t, 'DN_eff': dn,
         'rel_err': (dn - ref_dn) / ref_dn})
    rec({'meta': env_meta(), 'kw': KW})
    path = os.path.join(out_dir, 'pareto_default.json')
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(records, fh, indent=2)
    print('wrote %s' % path)


if __name__ == '__main__':
    main()
