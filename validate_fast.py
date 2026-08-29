import sys, time, json
import numpy as np
sys.path.insert(0, r"F:\codex\stiffGWpy")
import fast_sgwb as FS
from stiff_SGWB import LCDM_SG

CASES = [
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-3, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=3.6e-2, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-1, cr=1, T_re=2e3, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=1e1, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=1e4, kappa10=1e-2),
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1e-3),
    dict(r=1e-2, cr=1, T_re=2e3, kappa10=1.0),
    dict(r=3.6e-2, cr=0, T_re=1e3, kappa10=1.0),
    dict(r=1e-2, cr=0, T_re=2e3, kappa10=1e-2),
    dict(r=1e-1, cr=1, T_re=1e1, kappa10=1.0),
    dict(r=1e-2, cr=0, T_re=1e4, kappa10=1e-1),
]

def rel(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return np.abs(a-b)/np.maximum(np.abs(b), 1e-300)

def run_case(idx, kw):
    rec = dict(idx=idx, kw={k: repr(v) for k, v in kw.items()})
    m = LCDM_SG(**kw)
    if m.derived_param['N_inf'] is None:
        rec['skip'] = 'N_inf None (invalid combo)'
        return rec
    mo = LCDM_SG(**kw)
    t0 = time.perf_counter(); mo.SGWB_iter(); to = time.perf_counter()-t0
    rec['t_orig'] = to
    mf = LCDM_SG(**kw)
    t0 = time.perf_counter(); FS.SGWB_iter_fast(mf); tf1 = time.perf_counter()-t0
    mf2 = LCDM_SG(**kw)
    t0 = time.perf_counter(); FS.SGWB_iter_fast(mf2); tf2 = time.perf_counter()-t0
    rec['t_fast_first_ms'] = tf1*1e3; rec['t_fast_warm_ms'] = tf2*1e3
    rec['speedup'] = to/tf2
    rec['conv_orig'] = bool(mo.SGWB_converge)
    rec['conv_fast'] = bool(mf2.SGWB_converge)
    rec['status_match'] = bool(mo.SGWB_converge == mf2.SGWB_converge)
    if not (mo.SGWB_converge and mf2.SGWB_converge):
        rec['note'] = 'not both converged'
        return rec
    for name in ['f', 'log10OmegaGW', 'hubble']:
        a = np.asarray(getattr(mo, name), float); b = np.asarray(getattr(mf2, name), float)
        rec['len_'+name] = (int(a.size), int(b.size))
        if a.size == b.size:
            d = np.abs(a-b)
            rec[name+'_maxabs'] = float(d.max())
            rec[name+'_maxrel'] = float(rel(a, b).max()) if b.size else 0.0
        else:
            rec[name+'_maxabs'] = None; rec[name+'_maxrel'] = None
    a = np.asarray(mo.DN_gw, float); b = np.asarray(mf2.DN_gw, float)
    rec['len_DN_gw'] = (int(a.size), int(b.size))
    fin = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 0)
    rec['DN_gw_fin_maxrel'] = float(rel(a, b)[fin].max()) if fin.any() else None
    rec['DN_gw_last'] = (float(a[-1]), float(b[-1]))
    rec['DN_gw_last_rel'] = float(abs(a[-1]-b[-1])/abs(a[-1]))
    rec['kappa_r'] = (float(mo.kappa_r), float(mf2.kappa_r))
    rec['kappa_r_rel'] = float(abs(mo.kappa_r-mf2.kappa_r)/abs(mo.kappa_r))
    rec['DN_eff_final'] = (float(mo.cosmo_param['DN_eff']), float(mf2.cosmo_param['DN_eff']))
    g = rel(mo.g2, mf2.g2); g2f = np.isfinite(mo.g2) & np.isfinite(mf2.g2) & (np.abs(mf2.g2) > 0)
    rec['g2_fin_maxrel'] = float(g[g2f].max())
    return rec

if __name__ == '__main__':
    which = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else range(len(CASES))
    out = []
    for i in which:
        t0 = time.perf_counter()
        r = run_case(i, CASES[i])
        r['wall_s'] = time.perf_counter()-t0
        print(json.dumps(r, ensure_ascii=False), flush=True)
        out.append(r)
    fn = r"C:\Users\联想\AppData\Local\Temp\stiffgw_bench\validate_repo_out.jsonl"
    with open(fn, 'w', encoding='utf-8') as fh:
        for r in out:
            fh.write(json.dumps(r, ensure_ascii=False) + '\n')
