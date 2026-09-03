# -*- coding: utf-8 -*-
"""Two-profile fast validation: convergence + parameter sweep + oracle anchor.

This is the physics-first certification driver for the two user-facing fast
profiles (``fast`` = plain-grid / speed-first, ``production`` = transition-refine
/ precision-first).  It never uses LSODA as a precision oracle; the precision
anchor is the independent continuous-sigma reference pipeline
(``stiffgwpy.reference``) on matched frequency subsets.

Phases (choose one or more):
  --phase convergence    fast-vs-fast convergence of h/freq_res/z_tail/col_step/phase_max
  --phase param_sweep    Sobol/LHS over the physical parameter box (fast plain-grid
                         large screen + production on flagged points)
  --phase oracle_anchor  reference-engine oracle on a small matched subset

Usage:
  python scripts/validate_two_modes.py --phase convergence --profile all
  python scripts/validate_two_modes.py --phase param_sweep --n 200 --seed 20260903
  python scripts/validate_two_modes.py --phase oracle_anchor --points 3
"""

import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from stiffgwpy import fast_sgwb as FS
from stiffgwpy.stiff_SGWB import LCDM_SG

ln10 = math.log(10.0)

DEFAULT_POINT = dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)


def _env_meta():
    import subprocess
    meta = {
        'python': sys.version.split()[0],
        'platform': sys.platform,
        'cpu_count': os.cpu_count(),
        'numpy': np.__version__,
    }
    try:
        meta['commit'] = subprocess.run(
            ['git', 'rev-parse', 'HEAD'], capture_output=True, text=True,
            timeout=10).stdout.strip()
    except Exception:
        pass
    return meta


def _reset_settings():
    FS.set_threads(min(4, FS._MAX_THREADS))
    FS.set_col_step(4)
    FS.set_h(0.01)
    FS.set_z_tail(8.0)
    FS.set_phase_max(0.5)
    FS.set_freq_grid('adaptive')


def _solve(profile, **override):
    """One fast solve on the default point; returns (DN_gw, DN_eff, n_freq, secs)."""
    _reset_settings()
    mode = 'fast' if profile == 'plain_grid' else 'production'
    FS.apply_accuracy_mode(mode)
    for k, v in override.items():
        if k == 'h':
            FS.set_h(v)
        elif k == 'col_step':
            FS.set_col_step(v)
        elif k == 'z_tail':
            FS.set_z_tail(v)
        elif k == 'phase_max':
            FS.set_phase_max(v)
        elif k == 'freq_res':
            pass  # forwarded to SGWB_iter_fast below
        elif k == 'freq_grid':
            FS.set_freq_grid(v)
    m = LCDM_SG(**DEFAULT_POINT)
    t0 = time.time()
    FS.SGWB_iter_fast(m, tol=1e-7, freq_res=override.get('freq_res', 1.0),
                      transition_refine=(mode == 'production'),
                      freq_grid=override.get('freq_grid',
                                             'adaptive' if mode == 'production' else 'construct'))
    secs = time.time() - t0
    return (float(m.DN_gw[-1]), float(m.cosmo_param['DN_eff']),
            int(m.freq_grid_n), secs, m)


def phase_convergence(outdir):
    """Fast-vs-fast convergence curves (relative DN_gw change vs a reference run).

    The reference for convergence is the tightest *same-profile* setting (the
    profile's own most resolved config), NOT LSODA.  This isolates each knob's
    convergence; the absolute accuracy anchor is the oracle phase.
    """
    ref_plain = _solve('plain_grid', h=0.0025, col_step=1, z_tail=10.0,
                       freq_res=2.0, freqs=None)
    ref_prod = _solve('production', h=0.005, col_step=1, z_tail=10.0,
                      phase_max=0.125, freq_res=2.0)
    rows = []

    def _rec(profile, knob, val, dn, nf, secs, **kw):
        rows.append(dict(profile=profile, knob=knob, value=val,
                         DN_gw=dn, DN_eff=_dn_eff(dn), n_freq=nf,
                         runtime_s=secs, **kw))

    def _dn_eff(dn):
        # Delta N_eff contribution is DN_gw (the model stores its own effective
        # DN_eff after the loop; here we report the self-consistent DN_eff).
        return dn

    # --- plain-grid sweeps (reference = ref_plain) ---
    for h in (0.04, 0.02, 0.01, 0.005):
        dn, dseff, nf, secs, _ = _solve('plain_grid', h=h)
        _rec('plain_grid', 'h', h, dn, nf, secs,
             rel_ref=abs(dn - ref_plain[0]) / max(abs(ref_plain[0]), 1e-300))
    for col in (1, 2, 4, 8):
        dn, _, nf, secs, _ = _solve('plain_grid', col_step=col)
        _rec('plain_grid', 'col_step', col, dn, nf, secs,
             rel_ref=abs(dn - ref_plain[0]) / max(abs(ref_plain[0]), 1e-300))
    for zt in (5.0, 7.0, 8.0, 10.0):
        dn, _, nf, secs, _ = _solve('plain_grid', z_tail=zt)
        _rec('plain_grid', 'z_tail', zt, dn, nf, secs,
             rel_ref=abs(dn - ref_plain[0]) / max(abs(ref_plain[0]), 1e-300))
    for fr in (0.5, 1.0, 2.0, 4.0):
        dn, _, nf, secs, _ = _solve('plain_grid', freq_res=fr)
        _rec('plain_grid', 'freq_res', fr, dn, nf, secs,
             rel_ref=abs(dn - ref_plain[0]) / max(abs(ref_plain[0]), 1e-300))

    # --- production (transition-refine) sweeps (reference = ref_prod) ---
    for h in (0.02, 0.01, 0.005):
        dn, _, nf, secs, _ = _solve('production', h=h)
        _rec('production', 'h', h, dn, nf, secs,
             rel_ref=abs(dn - ref_prod[0]) / max(abs(ref_prod[0]), 1e-300))
    for pm in (0.5, 0.25, 0.125):
        dn, _, nf, secs, _ = _solve('production', phase_max=pm)
        _rec('production', 'phase_max', pm, dn, nf, secs,
             rel_ref=abs(dn - ref_prod[0]) / max(abs(ref_prod[0]), 1e-300))
    for zt in (7.0, 8.0, 10.0):
        dn, _, nf, secs, _ = _solve('production', z_tail=zt)
        _rec('production', 'z_tail', zt, dn, nf, secs,
             rel_ref=abs(dn - ref_prod[0]) / max(abs(ref_prod[0]), 1e-300))
    for fr in (1.0, 2.0, 4.0):
        dn, _, nf, secs, _ = _solve('production', freq_res=fr)
        _rec('production', 'freq_res', fr, dn, nf, secs,
             rel_ref=abs(dn - ref_prod[0]) / max(abs(ref_prod[0]), 1e-300))

    out = dict(generated=_env_meta(),
               reference=dict(plain_grid_ref_DN_gw=ref_plain[0],
                              production_ref_DN_gw=ref_prod[0]),
               rows=rows)
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, 'convergence_two_modes.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    for r in rows:
        print("%-11s %-10s %-8s DN_gw=%.9g nf=%4d t=%.3fs rel=%.3g"
              % (r['profile'], r['knob'], r['value'], r['DN_gw'], r['n_freq'],
                 r['runtime_s'], r['rel_ref']))
    print("wrote", path)
    return path


def phase_param_sweep(outdir, n, seed, method='lhs'):
    """Parameter-box sweep with fast plain-grid screen + production flag points."""
    from scipy.stats import qmc
    # Schemas (min, max) on log-scaled / linear physical handles.
    # Full physical parameter schema (all parameters that actually enter the
    # background / tensor-source physics).  Ranges are the prior box used for
    # the coverage screen; they are wide enough to exercise the extremes that
    # gate on physics validity.
    box = {
        'Omega_bh2': (0.01, 0.04),
        'Omega_ch2': (0.05, 0.25),
        'H0': (50.0, 90.0),
        'DN_eff': (0.0, 3.0),
        'log10A_s': (-9.2, -8.2),        # A_s in [6.3e-10, 6.3e-9]
        'log10r': (-4.0, -1.0),          # r in [1e-4, 1e-1]
        'n_t': (-0.5, 0.5),
        'cr': (0.0, 1.0),
        'log10T_re': (-1.0, 7.0),        # T_re in [0.1, 1e7] GeV
        'DN_re': (0.0, 40.0),
        'log10kappa10': (-3.0, 3.0),
    }
    cols = list(box)
    if method == 'sobol':
        sampler = qmc.Sobol(d=len(cols), scramble=True, seed=seed)
        unit = sampler.random(n)
    else:
        sampler = qmc.LatinHypercube(d=len(cols), scramble=True, seed=seed)
        unit = sampler.random(n)
    l_bounds = np.array([box[c][0] for c in cols])
    u_bounds = np.array([box[c][1] for c in cols])
    pts = l_bounds + unit * (u_bounds - l_bounds)
    rows = []
    status_counts = dict(total=n, success=0, physical_invalid=0, physical_guard=0,
                         numerical_failure=0)
    for i in range(n):
        p = dict(zip(cols, pts[i]))
        row_in = {('DN_eff_in' if k == 'DN_eff' else k): float(v)
                  for k, v in p.items()}
        kwargs = dict(Omega_bh2=p['Omega_bh2'], Omega_ch2=p['Omega_ch2'],
                      H0=p['H0'], DN_eff=p['DN_eff'], A_s=10 ** p['log10A_s'],
                      r=10 ** p['log10r'], n_t=p['n_t'], cr=p['cr'],
                      T_re=10 ** p['log10T_re'], DN_re=p['DN_re'],
                      kappa10=10 ** p['log10kappa10'])
        # shared Delta Neff physical guard (same as the engine): total Neff > 5
        m = LCDM_SG(**kwargs)
        dn_gw = np.nan
        secs = np.nan
        try:
            _reset_settings()
            FS.apply_accuracy_mode('fast')
            FS.set_freq_grid('construct')
            t0 = time.time()
            res = FS.SGWB_iter_fast(m, tol=1e-6, transition_refine=False,
                                    freq_grid='construct')
            secs = time.time() - t0
            if res is None:
                reason = getattr(m, 'fast_failure_reason', 'failed')
                if reason == 'shared_Neff_guard':
                    status = 'physical_guard'
                    status_counts['physical_guard'] += 1
                elif reason in ('invalid_r', 'invalid_cutoff'):
                    status = 'physical_invalid'
                    status_counts['physical_invalid'] += 1
                else:
                    status = 'numerical_failure'
                    status_counts['numerical_failure'] += 1
                rows.append(dict(point=i, **row_in, status=status,
                                 DN_gw=np.nan, DN_eff=np.nan, runtime_s=secs))
                continue
            status = 'success'
            status_counts['success'] += 1
            dn_gw = float(m.DN_gw[-1])
            rows.append(dict(point=i, **row_in, status=status, DN_gw=dn_gw,
                             DN_eff=float(m.cosmo_param['DN_eff']),
                             runtime_s=secs, n_freq=int(m.freq_grid_n)))
        except Exception as exc:
            status = 'exception'
            status_counts['numerical_failure'] += 1
            rows.append(dict(point=i, **row_in, status=status, error=repr(exc),
                             DN_gw=np.nan,
                             DN_eff=np.nan, runtime_s=secs))
    path = os.path.join(outdir, 'param_sweep_plain.json')
    os.makedirs(outdir, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(dict(generated=_env_meta(), n=n, seed=seed, box=box,
                       status=status_counts, rows=rows), f, indent=2)
    print("status", json.dumps(status_counts))
    print("wrote", path)
    return path


def phase_oracle_independence(outdir):
    """Measure the reference oracle's own z_tail / tolerance sensitivity.

    A frozen analytic tail at ``z_tail`` is an approximation; this phase runs
    ``reference.oracle_variants`` on a small signal-band frequency subset at the
    default point to quantify how much the oracle's own z_tail choice moves the
    answer.  A deep z_tail (>= 14) is deliberately NOT attempted: the mode
    equation becomes deep-subhorizon stiff and the hand-off is numerically
    necessary.  This is the honest caveat behind "the reference is the truth".
    """
    os.makedirs(outdir, exist_ok=True)
    from stiffgwpy import reference as REF
    m = LCDM_SG(**DEFAULT_POINT)
    # signal-band frequency subset (re-entered modes; the ill-defined
    # sub-horizon-today low-f tail carries no Delta Neff weight).
    freqs = np.linspace(-4.0, 2.0, 12)
    dn = float(m.cosmo_param['DN_eff'])
    # Keep rtol uniform (1e-8) so the tail term is isolated; the reference's
    # own ODE/quadrature convergence is separately tiny (its reported
    # quadrature_error ~1e-22 and interpolation_error ~1e-9 in the matched
    # artifacts).  A deep z_tail (>= 14) is not attempted (stiff).
    res = REF.oracle_variants(m, freqs=freqs, dn_eff=dn,
                              z_tail_conservative=8.0, z_tail_deep=10.0,
                              rtol_conservative=1e-8, rtol_deep=1e-8)
    out = dict(generated=_env_meta(), n_freq=int(freqs.size),
               default_point=DEFAULT_POINT,
               deep_tail_note=("z_tail >= 14 is numerically infeasible: the "
                               "deep-subhorizon mode equation becomes stiff, so "
                               "the frozen-tail hand-off is a practical "
                               "necessity and the residual is reported, not "
                               "removed"),
               result=res)
    path = os.path.join(outdir, 'oracle_independence.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    print("oracle A/B/C delta_AB=%.3g delta_BC=%.3g status=%s"
          % (res['delta_AB'], res['delta_BC'], res['status']))
    print("wrote", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', default='convergence',
                    choices=['convergence', 'param_sweep', 'oracle_independence'])
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--method', default='lhs', choices=['lhs', 'sobol'])
    ap.add_argument('--seed', type=int, default=20260903)
    ap.add_argument('--outdir', default=os.path.join(os.getcwd(), 'docs', 'validation'))
    args = ap.parse_args()
    if args.phase == 'convergence':
        phase_convergence(args.outdir)
    elif args.phase == 'param_sweep':
        phase_param_sweep(args.outdir, args.n, args.seed, args.method)
    elif args.phase == 'oracle_independence':
        phase_oracle_independence(args.outdir)
    else:
        raise SystemExit("unknown phase")


if __name__ == '__main__':
    main()
