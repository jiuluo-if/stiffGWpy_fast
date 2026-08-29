# -*- coding: utf-8 -*-
"""
plot_convergence.py -- audit phase 2: error-vs-resolution figures + summary.

Reads the JSONL written by ``convergence_study.py`` and produces:

* docs/convergence/conv_h_dn.png       DN_gw[-1] error vs expansion step h
* docs/convergence/conv_h_curve.png    curve relative-error stats vs h
* docs/convergence/conv_h_spectrum.png spectrum (dex) stats vs h
* docs/convergence/conv_colstep.png    col_step sensitivity (fast only)
* docs/convergence/conv_ztail.png      analytic-tail z_tail sensitivity
* docs/convergence/conv_freq.png       frequency-grid density sensitivity
* docs/convergence/summary.json        headline metrics + convergence orders

Convergence order is the slope of log10(error) vs log10(h) over the finest
3 points of the shared-grid path (engine vs LSODA@h_finest).  The same-grid
fast-vs-LSODA error is engine-only and is reported as a flat floor.
"""
import argparse
import json
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

POINT = 'default'


def load(path):
    rows = [json.loads(line) for line in open(path, encoding='utf-8') if line.strip()]
    return rows


def order(hs, errs, window=3):
    hs = np.asarray(hs, float)
    errs = np.asarray([max(float(e), 1e-300) for e in errs], float)
    if len(hs) < window + 1:
        return None
    m = slice(len(hs) - window, None)
    if np.allclose(hs[m], hs[m][0]):
        return None
    return float(np.polyfit(np.log10(hs[m]), np.log10(errs[m]), 1)[0])


def plot_h(rows, outdir):
    ref = next(r for r in rows if r['kind'] == 'lsoda_ref')
    hs = sorted({r['h'] for r in rows if r['kind'] == 'fast'})
    fast_dn = {r['h']: r['DN_gw_last'] for r in rows if r['kind'] == 'fast'}
    lsoda_dn = {r['h']: r['DN_gw_last'] for r in rows if r['kind'] == 'lsoda'}
    samegrid = {r['h']: r for r in rows if r['kind'] == 'samegrid_fast_vs_lsoda'}
    vs_ref = {}
    for r in rows:
        if r['kind'] == 'vs_ref':
            vs_ref.setdefault(r['engine'], {})[r['h']] = r
    ref_dn = ref['DN_gw_last']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    # --- DN_gw[-1] ---
    ax = axes[0]
    x = np.array(hs, float)
    ax.loglog(x, [abs(fast_dn[h]-ref_dn)/abs(ref_dn) for h in hs],
              'o-', label='fast vs LSODA finest', color='C0')
    xl = np.array([h for h in hs if h in lsoda_dn], float)
    ax.loglog(xl, [abs(lsoda_dn[h]-ref_dn)/abs(ref_dn) for h in xl],
              's--', label='LSODA vs LSODA finest', color='C1')
    xs = np.array([h for h in hs if h in samegrid], float)
    ax.loglog(xs, [samegrid[h]['m_DN_gw_last_rel'] for h in xs],
              '^:', label='fast vs LSODA same grid', color='C2')
    ax.set_xlabel('expansion step h (e-folds)')
    ax.set_ylabel('|rel err| of final $\\Delta N_{eff}$')
    ax.set_title('final $\\Delta N_{eff}$ vs h')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)

    # --- curve stats (signal region) ---
    ax = axes[1]
    for name, fmt in (('DN_gw', 'o-'), ('g2', 's--'), ('w2', '^:')):
        for engine, c in (('fast', 'C0'), ('lsoda', 'C1')):
            hh = [h for h in hs if h in vs_ref.get(engine, {})]
            errs = [vs_ref[engine][h]['m_%s' % name]['median'] for h in hh]
            ax.loglog(hh, errs, fmt, color=c, alpha=0.55,
                      label='%s %s median' % (engine, name))
    ax.set_xlabel('expansion step h')
    ax.set_ylabel('median rel err (signal region)')
    ax.set_title('DN_gw/g2/w2 curves vs LSODA finest')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7)

    # --- spectrum ---
    ax = axes[2]
    for engine, c in (('fast', 'C0'), ('lsoda', 'C1')):
        hh = [h for h in hs if h in vs_ref.get(engine, {})]
        dex = [vs_ref[engine][h]['m_spectrum']['dex_max'] for h in hh]
        dex50 = [vs_ref[engine][h]['m_spectrum']['dex_p50'] for h in hh]
        ax.loglog(hh, dex, 'o-', color=c, label='%s dex_max' % engine)
        ax.loglog(hh, dex50, 's--', color=c, alpha=0.55, label='%s dex_p50' % engine)
    hh = [h for h in hs if h in samegrid]
    same = [samegrid[h]['m_spectrum']['dex_max'] for h in hh]
    ax.loglog(hh, same, '^:', color='C2', label='same-grid dex_max')
    ax.set_xlabel('expansion step h')
    ax.set_ylabel('log10 Omega abs err (dex)')
    ax.set_title('spectrum vs LSODA finest')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'conv_h_dn.png'), dpi=130)
    plt.close(fig)

    # convergence order (shared-grid path, finest 3 points)
    hhs = np.array(hs, float)
    o_fast = order(hhs, [abs(fast_dn[h]-ref_dn) for h in hs])
    o_lsoda = order([h for h in hs if h in lsoda_dn],
                    [abs(lsoda_dn[h]-ref_dn) for h in hs if h in lsoda_dn])
    return {'h_fast_order': o_fast, 'h_lsoda_order': o_lsoda,
            'h_samegrid_dn_rel': {str(h): samegrid[h]['m_DN_gw_last_rel']
                                  for h in sorted(samegrid)},
            'h_fast_dn': {str(h): fast_dn[h] for h in hs},
            'h_lsoda_dn': {str(h): lsoda_dn[h] for h in hs if h in lsoda_dn},
            'h_ref_dn': ref_dn}


def plot_knob(rows, outdir, kind, xname, title, fname):
    recs = {r[xname]: r for r in rows if r['kind'] == kind}
    xs = sorted(recs)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    dn = [recs[x]['m_DN_gw_last_rel'] for x in xs]
    ax.semilogx(xs, dn, 'o-')
    ax.set_xlabel(xname)
    ax.set_ylabel('|rel err| final $\\Delta N_{eff}$')
    ax.set_title(title + ': final $\\Delta N_{eff}$')
    ax.grid(True, which='both', alpha=0.3)
    ax = axes[1]
    for name in ('DN_gw', 'g2', 'w2'):
        med = [recs[x]['m_%s' % name]['median'] for x in xs]
        mx = [recs[x]['m_%s' % name]['max'] for x in xs]
        ax.semilogx(xs, med, 'o-', label='%s median' % name)
        ax.semilogx(xs, mx, 's--', alpha=0.5, label='%s max' % name)
    ax.set_xlabel(xname)
    ax.set_ylabel('curve rel err (signal region)')
    ax.set_title(title + ': curves')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, fname), dpi=130)
    plt.close(fig)
    return {str(x): {'DN_gw_last_rel': recs[x]['m_DN_gw_last_rel'],
                     'DN_gw_median': recs[x]['m_DN_gw']['median'],
                     'g2_median': recs[x]['m_g2']['median'],
                     'w2_median': recs[x]['m_w2']['median']} for x in xs}


def plot_freq(rows, outdir):
    same = {r['freq_res']: r for r in rows if r['kind'] == 'freq_samegrid'}
    fast = {r['freq_res']: r for r in rows if r['kind'] == 'freq_fast'}
    xs = sorted(same)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    ax.semilogx(xs, [same[x]['m_DN_gw_last_rel'] for x in xs], 'o-',
                label='fast vs LSODA same grid')
    ax.semilogx(xs, [fast[x]['m_DN_gw_last_rel'] for x in xs], 's--',
                label='fast vs LSODA freq_res=2 ref')
    ax.set_xlabel('freq_res (density multiplier)')
    ax.set_ylabel('|rel err| final $\\Delta N_{eff}$')
    ax.set_title('frequency grid: final $\\Delta N_{eff}$')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.semilogx(xs, [same[x]['m_spectrum']['dex_max'] for x in xs], 'o-',
                label='same-grid dex_max')
    ax.semilogx(xs, [fast[x]['m_spectrum']['dex_max'] for x in xs], 's--',
                label='vs freq_res=2 dex_max')
    ax.semilogx(xs, [fast[x]['m_spectrum']['dex_p50'] for x in xs], '^:',
                label='vs freq_res=2 dex_p50')
    ax.set_xlabel('freq_res')
    ax.set_ylabel('dex error')
    ax.set_title('frequency grid: spectrum')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'conv_freq.png'), dpi=130)
    plt.close(fig)
    return {str(x): {'samegrid_dn_rel': same[x]['m_DN_gw_last_rel'],
                     'vs2_dn_rel': fast[x]['m_DN_gw_last_rel'],
                     'samegrid_dex_max': same[x]['m_spectrum']['dex_max'],
                     'vs2_dex_max': fast[x]['m_spectrum']['dex_max'],
                     'vs2_dex_p50': fast[x]['m_spectrum']['dex_p50']} for x in xs}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--point', default=POINT)
    ap.add_argument('--in', dest='inpath', default=None)
    ap.add_argument('--out', default=os.path.join('docs', 'convergence'))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    inpath = args.inpath or os.path.join(args.out, 'convergence_%s.jsonl' % args.point)
    rows = load(inpath)
    summary = {'meta': next(r for r in rows if r['kind'] in ('lsoda_ref', 'colstep_ref',
                                                             'ztail_ref', 'freq_ref')).get('meta', {}),
               'point': args.point}
    summary.update(plot_h(rows, args.out))
    if any(r['kind'] == 'colstep' for r in rows):
        summary['colstep'] = plot_knob(rows, args.out, 'colstep', 'col_step',
                                       'coarse-column step', 'conv_colstep.png')
    if any(r['kind'] == 'ztail' for r in rows):
        summary['ztail'] = plot_knob(rows, args.out, 'ztail', 'z_tail',
                                     'analytic-tail cutoff', 'conv_ztail.png')
    if any(r['kind'] == 'freq_samegrid' for r in rows):
        summary['freq'] = plot_freq(rows, args.out)
    with open(os.path.join(args.out, 'summary.json'), 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print('wrote %s' % os.path.join(args.out, 'summary.json'))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
