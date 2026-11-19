from pathlib import Path

def edit(name, old, new, count=1):
    path=Path("stiffgwpy_fast")/name
    s=path.read_text(encoding="utf8")
    assert s.count(old)>=count, (name,old)
    path.write_text(s.replace(old,new,count),encoding="utf8")

freq="freq_adaptive.py"
edit(freq,"""        H2 = h2_at(N)
        raw = -0.5 * N + 0.5 * _m.log(H2)
        if N < N_re_abs:
            raw = raw_re - 0.5 * (N - N_re_abs)""","""        if N < N_re_abs:
            raw = raw_re - 0.5 * (N - N_re_abs)
        else:
            raw = -0.5 * N + 0.5 * _m.log(h2_at(N))""")
bg="exact_background.py"
# All public primitive builders keep the old convention by default.
for name, args in [('exact_phi_s2','m, Nv, DN_eff, h'),('exact_phi_s2_grid','m, Nv, DN_eff'),('exact_phi_s2_split','m, Nv, DN_eff, sigma_nodes=None'),('fast_phi_s2_split','m, Nv, DN_eff, sigma_nodes=None')]:
    edit(bg,f'def {name}({args}):',f'def {name}({args}, *, relative_scale=False):')
for name in ['exact_phi_s2','exact_phi_s2_grid','exact_phi_s2_split']:
    p=Path("stiffgwpy_fast")/bg
    s=p.read_text(encoding="utf8");start=s.index('def '+name+'(');end=s.index('\ndef ',start+1)
    part=s[start:end].replace('S2 = np.exp(Psi)\n    S2inv = np.exp(-0.5 * Psi)','S2, S2inv = _primitive_scales(Psi, relative_scale)')
    p.write_text(s[:start]+part+s[end:],encoding="utf8")
edit(bg,'def exact_phi_s2(','''def _primitive_scales(psi, relative_scale):
    if relative_scale:
        units = np.ones_like(psi)
        return units, units
    return np.exp(psi), np.exp(-0.5 * psi)


def exact_phi_s2(''')
edit(bg,'Phi_grid, Phi_mid, S2, S2inv):\n    n = len(Nv)','Phi_grid, Phi_mid, S2, S2inv, relative_scale=False):\n    n = len(Nv)')
edit(bg,'S2[i] = np.exp(Psi)\n        S2inv[i] = np.exp(-0.5 * Psi)','S2[i] = 1.0 if relative_scale else np.exp(Psi)\n        S2inv[i] = 1.0 if relative_scale else np.exp(-0.5 * Psi)')
edit(bg,"workspace['S2'], workspace['S2inv'])","workspace['S2'], workspace['S2inv'], relative_scale)")
# Build variable grids on the actual physical MD/post-MD branches.
old="""    H2 = H2_vec(Nv, m, dn)
    raw = np.empty_like(Nv)
    raw_re = -0.5 * N_re_abs + 0.5 * math.log(float(H2_vec(np.array([N_re_abs]), m, dn)[0]))
    for i in range(Nv.size):
        if Nv[i] < N_re_abs:
            raw[i] = raw_re - 0.5 * (Nv[i] - N_re_abs)
        else:
            raw[i] = -0.5 * Nv[i] + 0.5 * math.log(float(H2[i]))"""
edit(bg,old,"""    post = Nv >= N_re_abs
    raw = np.empty_like(Nv)
    raw[post] = -0.5 * Nv[post] + 0.5 * np.log(H2_vec(Nv[post], m, dn))
    raw_re = raw[np.searchsorted(Nv, N_re_abs)]
    raw[~post] = raw_re - 0.5 * (Nv[~post] - N_re_abs)""",2)

fs='fast_sgwb.py'
edit(fs,'def prep_kernel(Nv, Sv, f_hor, freqs, h, ln10v, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus):','def prep_kernel(Nv, Sv, f_hor, freqs, h, ln10v, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus, relative_scale=False):')
edit(fs,'S2[i] = math.exp(Psi[i])\n        S2inv[i] = math.exp(-0.5 * Psi[i])','S2[i] = 1.0 if relative_scale else math.exp(Psi[i])\n        S2inv[i] = 1.0 if relative_scale else math.exp(-0.5 * Psi[i])')
edit(fs,'def prep_fast(m, Nv, freqs, h, variable_grid=False):','def prep_fast(m, Nv, freqs, h, variable_grid=False, relative_scale=False):')
edit(fs,'prep_kernel(Nv, Sv, f_hor, freqs, h, ln10, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus)','prep_kernel(Nv, Sv, f_hor, freqs, h, ln10, Phi_grid, Phi_mid, Psi, S2, S2inv, j0s, z0s, fp_minus, relative_scale)')
edit(fs,'            h_final = 0.0\n            if transition_refine:', '''            h_final = 0.0
            # WKB forms exp(N_h) and exp(-N); its absolute factors cease
            # to be normal before a translation-invariant product does.
            relative_scale = Nv[-1] >= _LOG_MAX_FLOAT - z_tail
            if transition_refine:''')
edit(fs,'def _uniform_step_terminal_width(', '_LOG_MAX_FLOAT = math.log(np.finfo(np.float64).max)\n\n\n def _uniform_step_terminal_width('.replace('\n def','\ndef'))
edit(fs,'m, Nv, freqs, h, variable_grid=False\n','m, Nv, freqs, h, variable_grid=False, relative_scale=relative_scale\n',2)
edit(fs,'m, Nv, m.cosmo_param["DN_eff"]\n','m, Nv, m.cosmo_param["DN_eff"], relative_scale=relative_scale\n')
edit(fs,'m, Nv, m.cosmo_param["DN_eff"], sigma_nodes=m.sigma\n','m, Nv, m.cosmo_param["DN_eff"], sigma_nodes=m.sigma, relative_scale=relative_scale\n',3)
edit(fs,'exact_phi_s2(m, Nv, m.cosmo_param["DN_eff"], h)','exact_phi_s2(m, Nv, m.cosmo_param["DN_eff"], h, relative_scale=relative_scale)')
edit(fs,'    z_tail_by_mode=None,\n):','    z_tail_by_mode=None,\n    relative_background=None,\n):')
edit(fs,'@njit(parallel=True, cache=True)\ndef solve_kernel(','''@njit(cache=True, inline="always")
def _mode_state(x, y, index, start, log_gain, relative_background):
    if relative_background is None:
        return x, y
    log_scale = relative_background[0]
    # A common exact gauge at the mode start cancels in every observable.
    scale = math.exp((log_scale[index] - log_scale[start]) + log_gain)
    return scale * x, scale * y


@njit(parallel=True, cache=True)
def solve_kernel(''')
edit(fs,'        k = j0\n        # Track','        log_gain = 0.0\n        k = j0\n        # Track')
p=Path('stiffgwpy_fast')/fs;s=p.read_text(encoding='utf8')
import re
s,n=re.subn(r'(?m)^( +)(assemble_main\(Ogw, Oj, Opgw, m, [^\n]+, S2\[k\], )xh, yh(, zz, Pt\))',r'\1ax, ay = _mode_state(xh, yh, k, j0, log_gain, relative_background)\n\1\2ax, ay\3',s);assert n==4
p.write_text(s,encoding='utf8')
edit(fs,'                k += 1\n                zz =','''                if relative_background is not None:
                    size = max(abs(xh), abs(yh))
                    if size > 2.0**450 or 0.0 < size < 2.0**-450:
                        xh /= size
                        yh /= size
                        log_gain += math.log(size)
                k += 1
                zz =''')
edit(fs,'        lyh = yh\n        last_z','        lyh = yh\n        last_gain = log_gain\n        last_z')
edit(fs,'                    lyh = yh\n                    last_z','                    lyh = yh\n                    last_gain = log_gain\n                    last_z')
edit(fs,'        if tail_state is not None and phase_weight[m]','        lxh, lyh = _mode_state(lxh, lyh, kend, j0, last_gain, relative_background)\n        if tail_state is not None and phase_weight[m]')
edit(fs,'            eNz = math.exp(Nv[kend] - last_z)','            eNz = math.exp(Nv[kend] - last_z) if relative_background is None else 0.0')
edit(fs,'                assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus)','''                if relative_background is None:
                    assemble_tail(Ogw, Oj, Opgw, m, slot, kk2, coeff, eNz, fp_i, Pt, ev_minus, fp_minus)
                else:
                    horizon = relative_background[1]
                    decay = (Nv[kend] - Nv[kk2]) - last_z
                    log_q = z0 + (horizon[j0] - horizon[kk2]) * ln10
                    th = coeff * math.exp(decay)
                    xf = coeff * math.exp(decay + log_q)
                    oj = -th * th / 3.0 * Pt
                    op = xf * xf / 36.0 * Pt
                    Oj[m, slot] = oj
                    Opgw[m, slot] = op
                    Ogw[m, slot] = 3.0 * op + oj''')
edit(fs,'                if tail_z_by_mode is None:\n                    solve_kernel','''                if relative_scale:
                    solve_kernel(*solve_args, z_tail_by_mode=tail_z_by_mode,
                                 relative_background=(Phi_grid - Nv, f_hor))
                elif tail_z_by_mode is None:
                    solve_kernel''')

