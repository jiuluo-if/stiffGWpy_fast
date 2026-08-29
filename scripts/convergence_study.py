# -*- coding: utf-8 -*-
"""
convergence_study.py -- audit phase 2: numerical convergence of the fast solver.

Sweeps the fast-solver resolution knobs and cross-validates against the LSODA
reference **on the same sigma grid** (both engines consume the same (Nv, Sv,
f_hor) grid produced by ``fast_sgwb.gen_fast`` at spacing ``h``), so that:

* fast(h) vs LSODA(h)  -> the fast solver's own error (Magnus stepping +
  tail anchor + assembly), independent of the shared expansion-grid error;
* engine(h) vs LSODA(h_finest) -> the shared expansion-grid convergence
  (sigma-spline error at the instantaneous-reheating kink, frequency-grid
  error, outer Delta N_eff self-consistency).

Sweeps:
  --h-sweep     h in {0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625}
  --colstep     col_step in {1, 2, 4, 8} (fast only; coarse-column/PCHIP)
  --ztail       z_tail in {3, 5, 7, 10, 15} (analytic-tail threshold)
  --freq        freq_res in {0.5, 1.0, 2.0} (frequency-grid density)

Metrics per the audit: Delta N_eff_final, DN_gw[-1], kappa_r, log10OmegaGW
(dex abs + linear-Omega rel on the signal region), and the DN_gw/g2/w2/hubble
curves (max / median / p95 / p99 / p99.9 / L2 / L-inf of the relative error).

Usage:
  python scripts/convergence_study.py --point default --h-sweep --colstep --ztail --freq
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

from stiffgwpy import fast_sgwb as FS
from stiffgwpy import global_param as gp
from stiffgwpy._metrics import dex_abs, rel_abs, rel_linear_omega, signal_mask
from stiffgwpy.stiff_SGWB import LCDM_SG

ln10 = math.log(10.0)

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'cr0': dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}
H_LIST = [0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625]
COLSTEP_LIST = [1, 2, 4, 8]
ZTAIL_LIST = [3.0, 5.0, 7.0, 10.0, 15.0]
FREQRES_LIST = [0.5, 1.0, 2.0]
PERCENTILES = (50, 95, 99, 99.9)


def env_meta():
    meta = {
        'python': sys.version.split()[0],
        'platform': sys.platform,
        'cpu_count': os.cpu_count(),
        'numpy': np.__version__,
    }
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


# ---------------- engines on a common grid ----------------

def _outer_finalize(m, DN_gw_new, Omega_nu, DN_eff_orig, ln10v=ln10):
    m.cosmo_param['DN_eff'] = DN_eff_orig + DN_gw_new
    m.hubble = math.log10(2*math.pi) + m.f_hor + (m.Nv[-1]-m.Nv)/ln10v
    m.DN_gw = (gp.Neff0 * np.multiply(m.g2,
               np.exp(2*(m.f_hor-m.f_hor[-1])*ln10v + 2*(m.Nv-m.Nv[-1])))
               / Omega_nu)
    m.Ogw_today = np.array([m.Ogw[i][-1] for i in range(len(m.f))])
    m.Opgw_today = np.array([m.Opgw[i][-1] for i in range(len(m.f))])
    m.Oj_today = np.array([m.Oj[i][-1] for i in range(len(m.f))])
    m.log10OmegaGW = np.log10(m.Ogw_today - m.Oj_today)
    m.kappa_r = m.cosmo_param['DN_eff']*7/8*(4/11)**(4/3)*gp.z_ratio**4
    m.SGWB_converge = True


def run_lsoda(kw, h, z_tail=5.0, freq_res=1.0, rtol=1e-8, tol=1e-7):
    """LSODA integration on the same (Nv, Sv, f_hor) grid the fast solver uses."""
    m = LCDM_SG(**kw)
    FS.gen_fast(m, h)          # grid at spacing h (same for both engines)
    FS.set_z_tail(z_tail)      # module state only affects gen_fast; grid already built
    Omega_nu = gp.Omega_nh2/m.derived_param['h']**2
    DN_eff_orig = m.cosmo_param['DN_eff']
    DN_gw_list = [0.0]
    DN_gw_new = 0.0
    DN_gw_min = 0.0
    DN_gw_max = 10.0
    iters = 0
    for iters in range(60):
        m.construct_f(freq_res)
        m.run_SGWB(z_tail=z_tail, rtol=rtol, atol=[1e-12, 1e-22, 1e-22])
        m.int_SGWB()
        DN_gw_new = gp.Neff0 * m.g2[-1] / Omega_nu
        if abs((gp.Neff0+DN_eff_orig+DN_gw_new)/(gp.Neff0+DN_eff_orig+DN_gw_list[-1]) - 1) < tol:
            break
        if DN_gw_new > DN_gw_list[-1] > DN_gw_min and DN_gw_max >= DN_gw_list[-1]:
            DN_gw_min = DN_gw_list[-1]
        elif DN_gw_new < DN_gw_list[-1] < DN_gw_max and DN_gw_min <= DN_gw_list[-1]:
            DN_gw_max = DN_gw_list[-1]
        if 0 < DN_gw_min <= DN_gw_max < 10:
            DN_gw_new = (DN_gw_min + DN_gw_max)/2
        m.cosmo_param['DN_eff'] = DN_eff_orig + DN_gw_new
        DN_gw_list.append(DN_gw_new)
    _outer_finalize(m, DN_gw_new, Omega_nu, DN_eff_orig)
    m._iters = iters + 1
    return m


def run_fast(kw, h, col_step, z_tail, freq_res, tol=1e-7):
    FS.set_col_step(col_step)
    FS.set_z_tail(z_tail)
    FS.set_h(h)
    m = LCDM_SG(**kw)
    FS.SGWB_iter_fast(m, tol=tol, freq_res=freq_res)
    m._iters = getattr(m, '_iters', 0)
    return m


# ---------------- metrics ----------------

def _curve_stats(a, b):
    """Curve relative-error stats on the *signal region* only.

    ``a`` is the reference array, ``b`` the test array (same N grid after
    resampling).  Near-zero channels make the plain relative error meaningless
    (e.g. g2/w2 in the deep-subhorizon start), so stats are computed only
    where ``|a|`` exceeds ``max(1e-12, 1e-8 * peak)``; the complementary region
    is reported separately with its maximum *absolute* difference.  Any NaN
    points are counted and excluded from the stats.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    absdiff = np.abs(a - b)
    finite = np.isfinite(a) & np.isfinite(b)
    n_nan = int((~finite).sum())
    a_f = np.where(finite, a, 0.0)
    b_f = np.where(finite, b, 0.0)
    rel = rel_abs(a_f, b_f)
    peak = float(np.max(np.abs(a_f))) if a_f.size else 0.0
    mask = finite & (np.abs(a_f) > max(1e-12, 1e-8 * peak))
    rel_sig = rel[mask] if mask.any() else np.array([0.0])
    nearzero = absdiff[~mask & finite]
    out = {'n_pts': int(a.size), 'n_nan': n_nan,
           'n_signal': int(mask.sum()),
           'max': float(rel_sig.max()),
           'L2': float(np.sqrt(np.mean(rel_sig ** 2))),
           'Linf': float(rel_sig.max()),
           'median': float(np.median(rel_sig)),
           'nearzero_abs_max': float(nearzero.max()) if nearzero.size else 0.0}
    for p in PERCENTILES:
        out['p%d' % p] = float(np.percentile(rel_sig, p))
    return out


def _spectrum_stats(lo_ref, lo_test):
    dex = dex_abs(lo_ref, lo_test)
    mask = signal_mask(lo_test)
    lin = rel_linear_omega(lo_ref[mask], lo_test[mask]) if mask.any() else np.array([0.0])
    out = {'dex_max': float(dex.max())}
    for p in PERCENTILES:
        out['dex_p%d' % p] = float(np.percentile(dex, p))
    out['lin_max'] = float(lin.max())
    out['lin_median'] = float(np.median(lin))
    return out


def metrics(ref, test):
    out = {}
    out['DN_gw_last_ref'] = float(ref.DN_gw[-1])
    out['DN_gw_last_test'] = float(test.DN_gw[-1])
    out['DN_gw_last_rel'] = (float(abs(ref.DN_gw[-1]-test.DN_gw[-1])/abs(ref.DN_gw[-1]))
                             if ref.DN_gw[-1] != 0 else None)
    out['DN_eff_ref'] = float(ref.cosmo_param['DN_eff'])
    out['DN_eff_test'] = float(test.cosmo_param['DN_eff'])
    out['kappa_r_rel'] = (float(abs(ref.kappa_r-test.kappa_r)/abs(ref.kappa_r))
                          if ref.kappa_r != 0 else None)
    for name in ('DN_gw', 'g2', 'w2', 'hubble'):
        a = np.asarray(getattr(ref, name), float)
        b = np.asarray(getattr(test, name), float)
        if a.shape != b.shape:
            # different expansion grids (h sweep): resample the test curve onto
            # the reference N grid before comparing.
            if not a.size or not b.size:
                out[name] = None
                continue
            Nref = np.asarray(getattr(ref, 'Nv'), float)
            Ntst = np.asarray(getattr(test, 'Nv'), float)
            b = np.interp(Nref, Ntst, b)
        if a.size:
            out[name] = _curve_stats(a, b)
        else:
            out[name] = None
    lo_r = np.asarray(ref.log10OmegaGW, float)
    lo_t = np.asarray(test.log10OmegaGW, float)
    if lo_r.shape != lo_t.shape:
        fr = np.asarray(ref.f, float)
        ft = np.asarray(test.f, float)
        # np.interp requires an *ascending* xp; the f grid is built descending.
        if ft.size > 1 and ft[0] > ft[-1]:
            ft = ft[::-1]
            lo_t = lo_t[::-1]
        lo_t = np.interp(fr, ft, lo_t)
    out['spectrum'] = _spectrum_stats(lo_r, lo_t)
    return out


def estimate_order(hs, errs, window=3):
    """Convergence order from the last `window` points (finest h)."""
    if len(hs) < window + 1:
        return None
    pairs = [(math.log(hs[i]), math.log(max(errs[i], 1e-300)))
             for i in range(len(hs)-window, len(hs))]
    # linear fit slope of log(err) vs log(h)
    x = np.array([p[0] for p in pairs])
    y = np.array([p[1] for p in pairs])
    if np.allclose(x, x[0]):
        return None
    return float(np.polyfit(x, y, 1)[0])


# ---------------- CLI ----------------

def main(argv=None):
    ap = argparse.ArgumentParser(description='fast-solver convergence study')
    ap.add_argument('--point', default='default', choices=sorted(CASES))
    ap.add_argument('--out', default=None, help='output directory')
    ap.add_argument('--h-sweep', action='store_true')
    ap.add_argument('--colstep', action='store_true')
    ap.add_argument('--ztail', action='store_true')
    ap.add_argument('--freq', action='store_true')
    ap.add_argument('--rtol', type=float, default=1e-8)
    ap.add_argument('--h-min-lsoda', type=float, default=0.00125,
                    help='finest h at which LSODA is run (cost ~ 1/h^0.5+)')
    args = ap.parse_args(argv)

    kw = CASES[args.point]
    out_dir = args.out or os.path.join('docs', 'convergence')
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, 'convergence_%s' % args.point)
    records = []

    def rec(record):
        records.append(record)
        print(json.dumps(record, ensure_ascii=False), flush=True)

    FS.set_col_step(1)
    FS.set_z_tail(5.0)
    FS.set_h(0.01)

    # ---- h sweep ----
    if args.h_sweep:
        h_lsoda = [h for h in H_LIST if h >= args.h_min_lsoda]
        t0 = time.perf_counter()
        m_ref = run_lsoda(kw, h=h_lsoda[-1], z_tail=5.0, rtol=args.rtol)
        rec({'kind': 'lsoda_ref', 'h': h_lsoda[-1], 'DN_gw_last': float(m_ref.DN_gw[-1]),
             't_s': time.perf_counter()-t0})
        for h in H_LIST:
            mf = run_fast(kw, h=h, col_step=4, z_tail=5.0, freq_res=1.0)
            rec({'kind': 'fast', 'h': h, 'DN_gw_last': float(mf.DN_gw[-1]),
                 'DN_eff': float(mf.cosmo_param['DN_eff']),
                 'conv': bool(mf.SGWB_converge)})
        for h in h_lsoda:
            t0 = time.perf_counter()
            ml = run_lsoda(kw, h=h, z_tail=5.0, rtol=args.rtol)
            rec({'kind': 'lsoda', 'h': h, 'DN_gw_last': float(ml.DN_gw[-1]),
                 'DN_eff': float(ml.cosmo_param['DN_eff']),
                 't_s': time.perf_counter()-t0})
            mf = run_fast(kw, h=h, col_step=4, z_tail=5.0, freq_res=1.0)
            rec({'kind': 'samegrid_fast_vs_lsoda', 'h': h,
                 **{('m_' + k): v for k, v in metrics(ml, mf).items()}})
            rec({'kind': 'vs_ref', 'engine': 'lsoda', 'h': h,
                 **{('m_' + k): v for k, v in metrics(m_ref, ml).items()}})
            rec({'kind': 'vs_ref', 'engine': 'fast', 'h': h,
                 **{('m_' + k): v for k, v in metrics(m_ref, mf).items()}})

    # ---- col_step sweep (fast only) ----
    if args.colstep:
        h0 = 0.01
        m0 = run_fast(kw, h=h0, col_step=1, z_tail=5.0, freq_res=1.0)
        rec({'kind': 'colstep_ref', 'h': h0, 'DN_gw_last': float(m0.DN_gw[-1])})
        for cs in COLSTEP_LIST:
            mf = run_fast(kw, h=h0, col_step=cs, z_tail=5.0, freq_res=1.0)
            rec({'kind': 'colstep', 'col_step': cs,
                 **{('m_' + k): v for k, v in metrics(m0, mf).items()}})

    # ---- z_tail sweep (tail formula is shared with LSODA; fast-only) ----
    if args.ztail:
        h0 = 0.01
        zdeep = ZTAIL_LIST[-1]
        mdeep = run_fast(kw, h=h0, col_step=1, z_tail=zdeep, freq_res=1.0)
        rec({'kind': 'ztail_ref', 'z_tail': zdeep, 'DN_gw_last': float(mdeep.DN_gw[-1])})
        for zt in ZTAIL_LIST:
            mf = run_fast(kw, h=h0, col_step=1, z_tail=zt, freq_res=1.0)
            rec({'kind': 'ztail', 'z_tail': zt,
                 **{('m_' + k): v for k, v in metrics(mdeep, mf).items()}})

    # ---- freq_res sweep ----
    if args.freq:
        h0 = 0.01
        mref = run_lsoda(kw, h=h0, z_tail=5.0, freq_res=2.0, rtol=args.rtol)
        rec({'kind': 'freq_ref', 'freq_res': 2.0, 'DN_gw_last': float(mref.DN_gw[-1])})
        for fr in FREQRES_LIST:
            mf = run_fast(kw, h=h0, col_step=1, z_tail=5.0, freq_res=fr)
            rec({'kind': 'freq_fast', 'freq_res': fr, 'DN_gw_last': float(mf.DN_gw[-1]),
                 **{('m_' + k): v for k, v in metrics(mref, mf).items()}})
            ml = run_lsoda(kw, h=h0, z_tail=5.0, freq_res=fr, rtol=args.rtol)
            rec({'kind': 'freq_samegrid', 'freq_res': fr,
                 **{('m_' + k): v for k, v in metrics(ml, mf).items()}})

    json_path = prefix + '.jsonl'
    with open(json_path, 'w', encoding='utf-8') as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
    print('wrote %s' % json_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
