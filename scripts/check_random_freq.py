# -*- coding: utf-8 -*-
"""
check_random_freq.py -- random 10-frequency spot check between the original
LCDM_SG.SGWB_iter() (old) and fast_sgwb.SGWB_iter_fast() (new).

Both solvers are run to self-consistency for a fixed parameter set, then 10
frequency channels are drawn at random (seeded, reproducible) from the shared
frequency grid. The two energy-spectrum curves log10OmegaGW(f) are plotted on
the same axes, the 10 sampled frequencies are marked, and the agreement
(|dlog10OmegaGW| and relative difference in OmegaGW) is quantified both at the
10 random frequencies and over the whole spectrum.

Usage:
    python check_random_freq.py [seed] [output.png]
"""

import os
import sys

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root
from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

PARAMS = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)   # baseline case A
TOL = 1e-4                                             # algorithm's own convergence tolerance


def main(seed=20260829, out_png=None):
    rng = np.random.default_rng(seed)

    # ---- old method (original LSODA solver) ----
    m_old = LCDM_SG(**PARAMS)
    m_old.SGWB_iter(engine='lsoda')
    if not m_old.SGWB_converge:
        raise RuntimeError('old method did not converge')

    # ---- new method (fast numba/Magnus solver) ----
    m_new = LCDM_SG(**PARAMS)
    FS.SGWB_iter_fast(m_new)
    if not m_new.SGWB_converge:
        raise RuntimeError('new method did not converge')

    f_old = np.asarray(m_old.f, float)
    f_new = np.asarray(m_new.f, float)
    if f_old.shape != f_new.shape or not np.allclose(f_old, f_new, rtol=1e-12, atol=1e-14):
        raise RuntimeError('frequency grids of the two runs differ')
    f = f_old

    log10Om_old = np.asarray(m_old.log10OmegaGW, float)
    log10Om_new = np.asarray(m_new.log10OmegaGW, float)
    Om_old = np.power(10.0, log10Om_old)
    Om_new = np.power(10.0, log10Om_new)
    rel_diff = np.abs(Om_old - Om_new) / np.maximum(Om_new, 1e-300)
    dex_diff = np.abs(log10Om_old - log10Om_new)

    # ---- randomly pick 10 frequency channels ----
    idx = np.sort(rng.choice(len(f), size=10, replace=False))

    print('=' * 72)
    print('Random 10-frequency spot check  (params: %s)' % PARAMS)
    print('seed=%d   N_freq_channels=%d   converged: old=%s new=%s'
          % (seed, len(f), m_old.SGWB_converge, m_new.SGWB_converge))
    print('-' * 72)
    print('%-6s %-14s %-16s %-16s %-12s %-12s' % ('#', 'log10(f/Hz)',
          'OmGW_old', 'OmGW_new', '|dlog10Om|', 'rel diff'))
    for k, i in enumerate(idx):
        print('%-6d %-14.4f %-16.4e %-16.4e %-12.3e %-12.3e'
              % (k + 1, f[i], Om_old[i], Om_new[i], dex_diff[i], rel_diff[i]))
    print('-' * 72)

    def stats(d, name):
        print('%-34s  at 10 random: max=%.3e  mean=%.3e   |  full grid: max=%.3e  mean=%.3e'
              % (name, d[idx].max(), d[idx].mean(), d.max(), d.mean()))

    stats(dex_diff, '|dlog10OmegaGW| (dex)')
    stats(rel_diff, 'relative diff in OmegaGW')
    frac10 = float((rel_diff[idx] < TOL).mean())
    frac_all = float((rel_diff < TOL).mean())
    print('%-34s  at 10 random: %d/10 < 1e-4  |  full grid: %d/%d < 1e-4'
          % ('channels within 1e-4 tolerance', int(round(frac10 * 10)),
             int(frac_all * len(f)), len(f)))
    print('=' * 72)

    # ---- figure: spectra overlay + residual ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [2.2, 1.4]})
    ax1.plot(f, log10Om_old, 'o-', ms=3, lw=1.2, color='#1f77b4', label='old: SGWB_iter (LSODA)')
    ax1.plot(f, log10Om_new, 'x--', ms=4, lw=1.0, color='#d62728', label='new: SGWB_iter_fast (Magnus)')
    for i in idx:
        ax1.axvline(f[i], color='gray', ls=':', lw=0.8, alpha=0.7)
    ax1.plot(f[idx], log10Om_old[idx], 'o', ms=9, mfc='none', mec='#1f77b4', mew=1.5)
    ax1.plot(f[idx], log10Om_new[idx], 'x', ms=10, mec='#d62728', mew=1.8)
    ax1.set_ylabel(r'$\log_{10}\Omega_{\mathrm{GW}}(f)$')
    ax1.set_title('Energy spectrum: old vs new solver (10 random frequency channels marked)')
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(alpha=0.3)

    ax2.semilogy(f, rel_diff, '-', color='#7f7f7f', lw=0.9, label='|dOm|/Om (full grid)')
    ax2.semilogy(f[idx], rel_diff[idx], 'o', ms=7, mfc='none', mec='#d62728', mew=1.6,
                 label='10 random frequencies')
    ax2.axhline(TOL, color='#2ca02c', ls='--', lw=1.2,
                label='algorithm tolerance 1e-4')
    ax2.set_ylabel('rel. diff  $|\\Omega_{\\mathrm{old}}-\\Omega_{\\mathrm{new}}|/\\Omega_{\\mathrm{new}}$')
    ax2.set_xlabel(r'$\log_{10}(f/\mathrm{Hz})$')
    ax2.legend(loc='best', fontsize=9)
    ax2.grid(alpha=0.3, which='both')

    fig.tight_layout()
    if out_png is None:
        out_png = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'random10_freq_spectra.png')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print('figure saved to: %s' % out_png)


if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20260829
    out = sys.argv[2] if len(sys.argv) > 2 else None
    main(seed=seed, out_png=out)
