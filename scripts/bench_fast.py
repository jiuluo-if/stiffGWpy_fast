# -*- coding: utf-8 -*-
"""
bench_fast.py -- reproducible performance comparison for stiffGWpy.

Runs the original LSODA-based SGWB_iter() and the accelerated
fast_sgwb.SGWB_iter_fast() on the same parameter grid, printing wall-clock
times, speedups and physics agreement for each case.

Usage:  python bench_fast.py
"""
import sys
import time
import statistics

import numpy as np

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root

from stiffgwpy import fast_sgwb          # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG    # noqa: E402

CASES = {
    "A baseline r=1e-2 cr=1 T_re=2e3 k10=1e-2": dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "B low r=1e-3":                           dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    "C r=3.6e-2":                             dict(r=3.6e-2, cr=1, T_re=2e3, kappa10=1e-2),
    "D high r=0.1":                           dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    "E low T_re=10 GeV":                      dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    "F high T_re=1e4 GeV":                    dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    "G low kappa10=1e-3":                     dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-3),
    "H high kappa10=1":                       dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    "I cr=0 kappa10=1 T_re=1e3":              dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    "J cr=0 baseline":                        dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    "K extreme r=0.1 T_re=10 k10=1":          dict(r=1e-1, cr=1, T_re=1e1, kappa10=1.0),
    "L cr=0 T_re=1e4 k10=0.1":                dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
}

NREP = 3  # fast warm repeats per case


def main():
    print("stiffGWpy benchmark: original SGWB_iter vs fast_sgwb.SGWB_iter_fast")
    print("=" * 100)
    print("%-32s %12s %14s %10s %18s" % ("case", "orig (s)", "fast med (ms)", "speedup", "DN_gw[-1] rel"))
    print("-" * 100)
    for name, kw in CASES.items():
        m = LCDM_SG(**kw)
        if m.derived_param['N_inf'] is None:
            print("%-32s %12s" % (name, "invalid combo"))
            continue
        t0 = time.perf_counter()
        m.SGWB_iter()
        t_orig = time.perf_counter() - t0

        fast_sgwb.SGWB_iter_fast(LCDM_SG(**kw))  # warm-up (JIT cache)
        times = []
        last = None
        for _ in range(NREP):
            mf = LCDM_SG(**kw)
            t0 = time.perf_counter()
            fast_sgwb.SGWB_iter_fast(mf)
            times.append(time.perf_counter() - t0)
            last = mf
        med = statistics.median(times)
        rel_dn = 0.0
        if m.SGWB_converge and last.SGWB_converge and m.DN_gw[-1] != 0:
            rel_dn = abs(m.DN_gw[-1] - last.DN_gw[-1]) / abs(m.DN_gw[-1])
        print("%-32s %12.3f %14.3f %9.0fx %18.2e" % (name, t_orig, med * 1e3, t_orig / med, rel_dn))
    print("-" * 100)
    print("Note: original runs are the dominant cost; total runtime ~2-4 min.")


if __name__ == "__main__":
    main()
