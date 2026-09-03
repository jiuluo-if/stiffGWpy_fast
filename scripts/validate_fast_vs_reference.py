# -*- coding: utf-8 -*-
"""validate_fast_vs_reference.py -- physics-first, LSODA-free validation.

Fast (production / matched-grid) vs the independent continuous-sigma DOP853
reference in ``stiffgwpy.reference``.  LSODA is never used as a truth anchor.

Phases (checkpointed to JSONL; reruns skip completed point ids):

1. ``--phase fast-sweep``:  Sobol/LHS draws over the physical parameter box,
   solved by the production fast engine (z8 adaptive grid).  Records the
   converged DN_eff, validity status and the local error budget.
2. ``--phase reference``:  named single points (default/stiff/lowT/highT/
   rad_dominant/tiny_r/transition/cr0_blue/extreme) plus stratified Sobol
   draws; each point solved by fast and by the continuous-sigma reference on
   the SAME grid-independent frequency grid and the SAME z_tail (matched
   tails, so the comparison isolates engine error).  Per point: per-mode dex
   and linear-Omega relative-error stats in the signal and transition bands,
   and the integrated Delta N_eff relative error.
3. ``--phase summary``:  max / p95 / p99 statistics over every reference
   point and the acceptance-gate verdicts.

Usage::

    python scripts/validate_fast_vs_reference.py --phase fast-sweep --n-sweep 240
    python scripts/validate_fast_vs_reference.py --phase reference --n-ref 16
    python scripts/validate_fast_vs_reference.py --phase summary
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scipy.stats import qmc  # noqa: E402

from stiffgwpy import fast_sgwb as FS  # noqa: E402
from stiffgwpy import global_param as gp  # noqa: E402
from stiffgwpy import reference as REF  # noqa: E402
from stiffgwpy._metrics import dex_abs, rel_linear_omega  # noqa: E402
from stiffgwpy.freq_adaptive import grid_independent_freqs  # noqa: E402
from stiffgwpy.stiff_SGWB import LCDM_SG  # noqa: E402

SEED = 20260831
ln10 = math.log(10.0)

BASE = dict(Omega_bh2=0.0223828, Omega_ch2=0.1201075, H0=67.32117,
            DN_eff=0.0, A_s=2.100549e-9)

# The GW physics is set by (r, n_t, cr, T_re, DN_re, kappa10); the LCDM
# anchors stay at BASE so Sobol coverage concentrates on Omega_GW physics.
SAMPLED = [
    ('r', 1e-6, 1e-1, 'log'),
    ('n_t', -0.5, 0.5, 'lin'),
    ('cr', 0.0, 1.0, 'disc'),
    ('T_re', 1e1, 1e6, 'log'),
    ('DN_re', 0.0, 30.0, 'lin'),
    ('kappa10', 1e-6, 1.0, 'log'),
]

SINGLE_POINTS = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'stiff': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    'highT': dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    'rad_dominant': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-6),
    'tiny_r': dict(r=1e-6, cr=1, T_re=2e3, kappa10=1e-2),
    'transition': dict(r=5e-3, cr=1, T_re=5e2, kappa10=1e-1),
    # cr=0 blue-tilt branch: fully radiative scaling, zero stiff plateau.
    'cr0_blue': dict(r=1e-2, cr=0, n_t=0.2, DN_re=5.0, T_re=1e3,
                     kappa10=1e-3),
    # near the shared-Neff guard but physical (extreme stiff-enhanced r).
    'extreme': dict(r=3e-2, cr=1, T_re=1e4, kappa10=1e-1),
}

SIGNAL_BAND = (-6.0, 1.0)
TRANSITION_BAND = (-2.0, 0.0)


def env_meta():
    meta = {'python': sys.version.split()[0], 'platform': sys.platform,
            'numpy': np.__version__}
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


def to_physical(u, lo, hi, mode):
    if mode == 'log':
        return lo * (hi / lo) ** u
    if mode == 'disc':
        return 1.0 if u >= 0.5 else 0.0
    return lo + (hi - lo) * u


def physical_params(uvec):
    kw = dict(BASE)
    for (name, lo, hi, mode) in SAMPLED:
        kw[name] = to_physical(float(uvec[name]), lo, hi, mode)
    return kw


def sobol_design(n, seed=SEED):
    d = len(SAMPLED)
    sampler = qmc.Sobol(d=d, scramble=True, seed=seed)
    u = sampler.random(n)
    names = [s[0] for s in SAMPLED]
    return [dict(zip(names, u[i])) for i in range(n)]


def band_mask(freqs, band):
    return (freqs >= band[0]) & (freqs <= band[1])


def _stats_of(arr):
    arr = np.asarray(arr, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    return {'n': int(arr.size), 'max': float(arr.max()),
            'p95': float(np.percentile(arr, 95)),
            'p99': float(np.percentile(arr, 99)),
            'mean': float(arr.mean())}


def band_stats(lo_ref, lo_test, freqs, band):
    """dex (log10) and linear-Omega relative error stats inside a band."""
    mask = band_mask(np.asarray(freqs), band)
    dex = dex_abs(np.asarray(lo_ref)[mask], np.asarray(lo_test)[mask])
    rel = rel_linear_omega(np.asarray(lo_ref)[mask], np.asarray(lo_test)[mask])
    return {'dex': _stats_of(dex), 'rel': _stats_of(rel)}


def fast_point_sweep(kw, z_tail=8.0):
    cfg = FS.apply_accuracy_mode('production')
    FS.set_z_tail(z_tail)
    m = LCDM_SG(**kw)
    t0 = time.perf_counter()
    try:
        res = FS.SGWB_iter_fast(m, tol=cfg['tol'],
                                transition_refine=cfg.get('transition_refine', False))
    except Exception as exc:
        return {'status': 'exception',
                'error': '%s: %s' % (type(exc).__name__, exc),
                'dt': time.perf_counter() - t0}
    dt = time.perf_counter() - t0
    if res is None:
        return {'status': 'rejected',
                'reason': getattr(m, 'fast_failure_reason', 'unknown'), 'dt': dt}
    if not getattr(m, 'SGWB_converge', False):
        return {'status': 'rejected', 'reason': 'not_converged', 'dt': dt}
    b = FS.estimate_local_error(m)
    return {
        'status': 'ok', 'dt': dt,
        'DN_eff': float(m.cosmo_param['DN_eff']),
        'DN_gw': float(np.asarray(m.DN_gw)[-1]),
        'DN_gw_error': float(b['DN_gw_error']),
        'Delta_Neff_abs_error': float(b['Delta_Neff_abs_error']),
        'n_freq': int(len(m.f)),
        'transition_refine_used': bool(getattr(m, 'transition_refine_used', False)),
        'handoff_eps_max': float(b['handoff_eps_max']),
    }


def reference_point(kw, z_tail=7.0, rtol=1e-9, freq_res=1.0,
                    label=None, kind='single'):
    """Fast + continuous-sigma reference, SAME grid and SAME z_tail."""
    rec = {'label': label or 'sobol', 'kind': kind, 'params': dict(kw),
           'z_tail': z_tail, 'rtol': rtol, 'freq_res': freq_res,
           'status': 'ok'}
    cfg = FS.apply_accuracy_mode('production')
    FS.set_z_tail(z_tail)
    m = LCDM_SG(**kw)
    t0 = time.perf_counter()
    try:
        res = FS.SGWB_iter_fast(m, tol=cfg['tol'],
                                transition_refine=cfg.get('transition_refine', False),
                                freq_grid='grid_independent', freq_res=freq_res)
    except Exception as exc:
        rec['status'] = 'fast_exception'
        rec['error'] = '%s: %s' % (type(exc).__name__, exc)
        return rec
    rec['fast_dt'] = time.perf_counter() - t0
    if res is None:
        rec['status'] = 'fast_rejected'
        rec['reason'] = getattr(m, 'fast_failure_reason', 'unknown')
        return rec
    dn_total = float(m.cosmo_param['DN_eff'])
    dn_gw_fast = float(np.asarray(m.DN_gw)[-1])
    f_fast = np.sort(np.asarray(m.f, dtype=float))
    lo_fast = np.asarray(m.log10OmegaGW, dtype=float)
    order = np.argsort(np.asarray(m.f, dtype=float))
    f_fast = np.asarray(m.f, dtype=float)[order]
    lo_fast = lo_fast[order]
    G = grid_independent_freqs(m, freq_res)[0].astype(np.float64)
    order_g = np.argsort(G)
    G = G[order_g]
    t0 = time.perf_counter()
    try:
        Ogw, Oj, Opgw, used = REF.spectrum_reference(m, G, dn_total,
                                                     z_tail=z_tail, rtol=rtol,
                                                     workers=8)
    except Exception as exc:
        rec['status'] = 'ref_exception'
        rec['error'] = '%s: %s' % (type(exc).__name__, exc)
        return rec
    rec['ref_dt'] = time.perf_counter() - t0
    g2, qerr, ierr = REF.integrate_spectrum(G, Ogw, Oj)
    Omega_nu = gp.Omega_nh2 / m.derived_param['h'] ** 2
    dn_gw_ref = float(gp.Neff0 * g2 / Omega_nu)
    lo_ref = np.log10(np.maximum(Ogw - Oj, 1e-300))
    # Node alignment: fast grid == reference grid up to the 1e-13 background
    # drift of the self-consistency loop; nearest-node match is exact.
    j = np.argmin(np.abs(f_fast[:, None] - G[None, :]), axis=1)
    f_match = f_fast
    lo_ref_m = lo_ref[j]
    grid_dev = float(np.max(np.abs(f_match - G)))
    rec['grid_dev'] = grid_dev
    rec['DN_gw_fast'] = dn_gw_fast
    rec['DN_gw_ref'] = dn_gw_ref
    rec['DN_gw_rel'] = (dn_gw_fast - dn_gw_ref) / dn_gw_ref
    rec['quadrature_error'] = float(qerr)
    rec['interpolation_error'] = float(ierr)
    rec['used_tail_frac'] = float(np.mean(used))
    rec['n_freq'] = int(G.size)
    for bname, band in (('all', (-30.0, 10.0)),
                        ('signal', SIGNAL_BAND),
                        ('transition', TRANSITION_BAND)):
        rec[bname] = band_stats(lo_ref_m, lo_fast, f_match, band)
    return rec


def load_lines(path):
    out = []
    if os.path.exists(path):
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                line = line.strip()
                if line:
                    out.append(json.loads(line))
    return out


def append_line(path, rec):
    with open(path, 'a', encoding='utf-8') as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + '\n')


def done_ids(path):
    seen = set()
    for rec in load_lines(path):
        seen.add(rec.get('label'))
    return seen


def run_fast_sweep(args):
    path = os.path.join(args.out, 'fast_sweep.jsonl')
    os.makedirs(args.out, exist_ok=True)
    design = sobol_design(args.n_sweep)
    seen = done_ids(path)
    n_ok = 0
    t_start = time.perf_counter()
    for i, u in enumerate(design):
        kw = physical_params(u)
        label = 'fast_%04d' % i
        if label in seen:
            recs = [r for r in load_lines(path) if r.get('label') == label]
            if recs and recs[-1].get('status') == 'ok':
                n_ok += 1
            continue
        rec = fast_point_sweep(kw, z_tail=args.fast_z_tail)
        rec['label'] = label
        rec['u'] = {k: float(v) for k, v in u.items()}
        rec['params'] = kw
        append_line(path, rec)
        if rec['status'] == 'ok':
            n_ok += 1
        if (i + 1) % 20 == 0:
            print('fast sweep %d/%d ok=%d dt=%.0fs' % (
                i + 1, args.n_sweep, n_ok, time.perf_counter() - t_start),
                flush=True)
    print('fast sweep done: ok=%d / %d (%.1f%%)' % (
        n_ok, args.n_sweep, 100.0 * n_ok / max(1, args.n_sweep)), flush=True)


def run_reference_phase(args):
    path = os.path.join(args.out, 'reference_points.jsonl')
    os.makedirs(args.out, exist_ok=True)
    seen = done_ids(path)
    cases = []
    for name, kw in SINGLE_POINTS.items():
        cases.append((name, kw, 'single'))
    design = sobol_design(args.n_sweep)
    n_sobol = max(0, args.n_ref - len(cases))
    stride = max(1, args.n_sweep // max(1, n_sobol))
    for i in range(0, args.n_sweep, stride):
        if n_sobol <= 0:
            break
        cases.append(('sobol_%04d' % i, physical_params(design[i]), 'sobol'))
        n_sobol -= 1
    if args.points:
        wanted = [p.strip() for p in args.points.split(',') if p.strip()]
        cases = [c for c in cases if any(
            c[0] == w or c[0].startswith(w) for w in wanted)]
        if not cases:
            sys.exit('--points matched no cases; choose from %s'
                     % ', '.join(list(SINGLE_POINTS) + ['sobol']))
    for label, kw, kind in cases:
        if label in seen:
            continue
        rec = reference_point(kw, z_tail=args.ref_z_tail, rtol=args.ref_rtol,
                              freq_res=args.ref_freq_res, label=label,
                              kind=kind)
        append_line(path, rec)
        seen.add(label)
        print('%-16s %-9s DN_gw_fast=%.9f DN_gw_ref=%.9f rel=%+.3e '
              'dex_sig_max=%.2e rel_sig_max=%.2e dt_f=%.1fs dt_r=%.1fs'
              % (label, kind, rec.get('DN_gw_fast', float('nan')),
                 rec.get('DN_gw_ref', float('nan')),
                 rec.get('DN_gw_rel', float('nan')),
                 (rec.get('signal') or {}).get('dex', {}).get('max', float('nan')),
                 (rec.get('signal') or {}).get('rel', {}).get('max', float('nan')),
                 rec.get('fast_dt', 0.0), rec.get('ref_dt', 0.0)), flush=True)
    print('reference phase done.', flush=True)


def run_summary(args):
    path = os.path.join(args.out, 'reference_points.jsonl')
    recs = [r for r in load_lines(path) if r.get('status') == 'ok'
            and r.get('DN_gw_ref')]
    if not recs:
        sys.exit('no completed reference points under %s' % path)
    rels = np.abs([r['DN_gw_rel'] for r in recs])
    summary = {'n_points': len(recs), 'meta': env_meta(),
               'settings': {'z_tail': recs[0].get('z_tail'),
                            'rtol': recs[0].get('rtol'),
                            'freq_res': recs[0].get('freq_res')},
               'DN_gw_rel_abs': {
                   'max': float(rels.max()),
                   'p95': float(np.percentile(rels, 95)),
                   'p99': float(np.percentile(rels, 99)),
                   'median': float(np.median(rels))},
               'acceptance': {}}
    for band in ('signal', 'transition'):
        dex = [r[band]['dex']['max'] for r in recs
               if (r.get(band) or {}).get('dex')]
        rel = [r[band]['rel']['max'] for r in recs
               if (r.get(band) or {}).get('rel')]
        summary[band] = {'dex_max': max(dex), 'dex_p95': float(
            np.percentile(dex, 95)), 'rel_max': max(rel),
            'rel_p95': float(np.percentile(rel, 95))}
        summary['acceptance']['%s_rel_lt_1e-3' % band] = bool(max(rel) < 1e-3)
    summary['acceptance']['DN_gw_rel_lt_1e-4'] = bool(rels.max() < 1e-4)
    out = os.path.join(args.out, 'validation_summary.json')
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1), flush=True)
    print('saved', out, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phase', default='summary',
                    choices=('fast-sweep', 'reference', 'summary'))
    ap.add_argument('--out', default='docs/paramsweep_ref')
    ap.add_argument('--n-sweep', type=int, default=240)
    ap.add_argument('--n-ref', type=int, default=16)
    ap.add_argument('--points', default=None,
                    help='comma list to restrict the reference phase '
                         '(single names or sobol prefix)')
    ap.add_argument('--fast-z-tail', type=float, default=8.0)
    ap.add_argument('--ref-z-tail', type=float, default=7.0)
    ap.add_argument('--ref-rtol', type=float, default=1e-9)
    ap.add_argument('--ref-freq-res', type=float, default=2.0)
    args = ap.parse_args(argv)
    if args.phase == 'fast-sweep':
        run_fast_sweep(args)
    elif args.phase == 'reference':
        run_reference_phase(args)
    else:
        run_summary(args)


if __name__ == '__main__':
    main()

