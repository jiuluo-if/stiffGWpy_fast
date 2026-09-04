# -*- coding: utf-8 -*-
"""benchmark_reference.py -- physics-first benchmark: fast vs production vs reference.

Unlike ``scripts/bench_fast.py`` (which reports speedup and fast-vs-LSODA), this
script benchmarks against the *independent high-accuracy reference* in
``stiffgwpy_fast.reference``.  That reference uses a continuous ``sigma(N)``
evaluator (no fixed-step grid through the reheating kink) and a high-order
adaptive ODE (DOP853), so it exposes the two *shared* errors the old
fast-vs-LSODA framing hid.

On a representative log-frequency subset it reports:

* Omega_GW dex error for fast and production vs reference;
* Delta N_eff relative difference vs the reference bolometric integral;
* an ODE error estimate (reference sensitivity to rtol);
* a tail error estimate (reference sensitivity to z_tail);
* a quadrature / interpolation error estimate from the reference integral;
* runtime for each configuration.

The reference is intentionally slow (no frequency parallelism, tight
tolerances), so it runs on a frequency subset by default.  ``--freq-full``
runs it on the full model grid (much slower).  ``--with-lsoda`` additionally
re-runs the original LSODA path on the full grid for the fast-vs-LSODA column.
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

from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast import reference as REF
from stiffgwpy_fast._metrics import dex_abs, rel_linear_omega, signal_mask
from stiffgwpy_fast.stiff_SGWB import LCDM_SG

ln10 = math.log(10.0)

CASES = {
    'default': dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    'cr0': dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    'lowT': dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
}

# Representative log10(f/Hz) covering the low-f tail, the spectrum knee, the
# stiff enhancement and the high-f cutoff.
FREQ_SUBSET = [-18.4, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, -15.0,
               -14.5, -14.0, -13.0, -12.0, -11.0, -10.0, -9.0, -8.0,
               -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0,
               3.0, 4.0, 5.0, 6.0]


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


def run_fast(kw, z_tail=None):
    cfg = FS.apply_accuracy_mode('production')
    if z_tail is not None:
        FS.set_z_tail(z_tail)
    m = LCDM_SG(**kw)
    t0 = time.perf_counter()
    # Forward the preset's continuous-sigma (kink breakpoint) path: production
    # mode is defined with transition_refine=True (see ACCURACY_MODES), and
    # running the plain spline-sigma path silently loses ~0.2% of Delta N_eff.
    FS.SGWB_iter_fast(m, tol=cfg['tol'], transition_refine=cfg.get('transition_refine', False))
    dt = time.perf_counter() - t0
    return m, dt


def run_lsoda(kw):
    m = LCDM_SG(**kw)
    t0 = time.perf_counter()
    m.SGWB_iter(engine='lsoda')
    dt = time.perf_counter() - t0
    return m, dt


def _spectrum_on_subset(m, subset):
    """Interpolate a model's spectrum (f, log10OmegaGW) onto a log-f subset."""
    fs = np.asarray(m.f, dtype=float)
    lo = np.asarray(m.log10OmegaGW, dtype=float)
    order = np.argsort(fs)
    fs = fs[order]
    lo = lo[order]
    return np.interp(subset, fs, lo)


def _ref_model(kw):
    return LCDM_SG(**kw)


def _rex_stats(lo_ref, lo_test):
    dex = dex_abs(lo_ref, lo_test)
    mask = signal_mask(lo_ref)
    lin = rel_linear_omega(lo_ref[mask], lo_test[mask]) if mask.any() else np.array([0.0])
    return {
        'dex_max': float(dex.max()),
        'dex_p50': float(np.percentile(dex, 50)),
        'dex_p95': float(np.percentile(dex, 95)),
        'lin_max': float(lin.max()),
        'lin_median': float(np.median(lin)),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description='physics-first fast/production/reference benchmark')
    ap.add_argument('--point', default='default', choices=sorted(CASES))
    ap.add_argument('--out', default=None)
    ap.add_argument('--freq-full', action='store_true')
    ap.add_argument('--with-lsoda', action='store_true')
    ap.add_argument('--no-ode-error', action='store_true',
                    help='skip the extra reference rtol=1e-10 run (saves ~1/3 time)')
    ap.add_argument('--no-tail-error', action='store_true',
                    help='skip the extra reference tail-sensitivity run')
    ap.add_argument('--z-tail', type=float, default=7.0,
                    help='reference z_tail for the primary run (default 7.0 = '
                         'matched to fast production mode)')
    args = ap.parse_args(argv)

    kw = CASES[args.point]
    out_dir = args.out or 'docs/archive/reference'
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, 'benchmark_%s' % args.point)
    records = []

    def rec(record):
        records.append(record)
        print(json.dumps(record, ensure_ascii=False), flush=True)

    subset = None if args.freq_full else np.array(FREQ_SUBSET)

    # ---- fast production ----
    m_fast, dt_fast = run_fast(kw, z_tail=args.z_tail)
    dn_fast = float(m_fast.cosmo_param['DN_eff'])
    dn_gw_fast = float(m_fast.DN_gw[-1])
    kappa_r_fast = float(m_fast.kappa_r)
    rec({'kind': 'fast', 'runtime_s': dt_fast, 'DN_eff': dn_fast,
         'DN_gw': dn_gw_fast, 'kappa_r': kappa_r_fast, 'n_freq': len(m_fast.f)})

    # ---- reference (physics-first truth) at the fast-converged DN_eff ----
    t0 = time.perf_counter()
    ref = REF.run_reference(_ref_model(kw),
                            dn_eff=dn_fast, freq_res=1.0, z_tail=args.z_tail, rtol=1e-11,
                            freq_subset=subset, self_consistent=False)
    dt_ref = time.perf_counter() - t0
    rec({'kind': 'reference', 'runtime_s': dt_ref, 'DN_gw': float(ref['DN_gw']),
         'kappa_r': float(ref['kappa_r']), 'n_freq': ref['n_freq'],
         'quadrature_error': ref['quadrature_error'],
         'interpolation_error': ref['interpolation_error'],
         'used_tail_frac': float(np.mean(ref['used_tail']))})

    # Compare fast / production spectrum to the reference on the same freqs.
    if subset is not None:
        lo_ref = ref['log10OmegaGW']
        lo_fast = _spectrum_on_subset(m_fast, subset)
        rec({'kind': 'spec_err', 'engine': 'fast',
             **{('m_' + k): v for k, v in _rex_stats(lo_ref, lo_fast).items()}})

    # ---- reference ODE sensitivity (rtol) ----
    if not args.no_ode_error:
        ref_lo = REF.run_reference(_ref_model(kw),
                                   dn_eff=dn_fast, freq_res=1.0, z_tail=args.z_tail,
                                   rtol=1e-10, freq_subset=subset, self_consistent=False)
        ode_err = dex_abs(ref['log10OmegaGW'], ref_lo['log10OmegaGW'])
        rec({'kind': 'ode_error', 'rtol_lo': 1e-10, 'rtol_hi': 1e-11,
             'dex_max': float(ode_err.max()),
             'dex_p95': float(np.percentile(ode_err, 95)),
             'DN_gw_rel': (float(abs(ref['DN_gw'] - ref_lo['DN_gw']) /
                                 abs(ref['DN_gw'])) if ref['DN_gw'] != 0 else None)})

    # ---- reference tail sensitivity (z_tail) ----
    if not args.no_tail_error:
        ref_tail = REF.run_reference(_ref_model(kw),
                                     dn_eff=dn_fast, freq_res=1.0, z_tail=5.0,
                                     rtol=1e-11, freq_subset=subset, self_consistent=False)
        tail_err = dex_abs(ref['log10OmegaGW'], ref_tail['log10OmegaGW'])
        rec({'kind': 'tail_error', 'zt_lo': 5.0, 'zt_hi': args.z_tail,
             'dex_max': float(tail_err.max()),
             'dex_p95': float(np.percentile(tail_err, 95)),
             'DN_gw_rel': (float(abs(ref['DN_gw'] - ref_tail['DN_gw']) /
                                 abs(ref['DN_gw'])) if ref['DN_gw'] != 0 else None)})

    # ---- optional LSODA ----
    if args.with_lsoda:
        m_ls, dt_ls = run_lsoda(kw)
        rec({'kind': 'lsoda', 'runtime_s': dt_ls, 'DN_eff': float(m_ls.cosmo_param['DN_eff']),
             'DN_gw': float(m_ls.DN_gw[-1]), 'kappa_r': float(m_ls.kappa_r),
             'n_freq': len(m_ls.f)})
        if subset is not None:
            lo_ls = _spectrum_on_subset(m_ls, subset)
            rec({'kind': 'spec_err', 'engine': 'lsoda',
                 **{('m_' + k): v for k, v in _rex_stats(lo_ref, lo_ls).items()}})

    rec({'meta': env_meta(), 'point': args.point, 'kw': kw,
         'subset': None if subset is None else subset.tolist(),
         'freq_full': args.freq_full})

    json_path = prefix + '.jsonl'
    with open(json_path, 'w', encoding='utf-8') as fh:
        for r in records:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
    print('wrote %s' % json_path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
