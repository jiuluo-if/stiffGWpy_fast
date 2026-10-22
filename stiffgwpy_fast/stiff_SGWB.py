"""High-level LCDM + stiff-matter + primordial-SGWB model and engine dispatch.

中文：本模块定义 `LCDM_SG` 和 `SGWB_iter()` 公共入口，在 fast、reference 与 LSODA 路径间
分派，并记录每次求解的引擎和回退遥测。数值实现分别位于 `fast_sgwb.py` 与 `reference.py`。
"""

import math
import multiprocessing as mp
import os
import warnings
from functools import partial

import numpy as np
from numpy import concatenate as cat
from scipy import integrate

from .functions import solve_SGWB
from .global_param import *
from .LCDM_stiff_Neff import LCDM_SN

MAX_ITER = 60            # 外层自洽二分迭代上限。
_DEFAULT_ACCURACY_MODE = object()
# The sentinel distinguishes an omitted argument from explicit ``None``.
# 对象哨兵用于区分“未传 accuracy_mode”与“显式传入 None”的兼容调用。


def _sgwb_pool_size():
    """Number of frequency-parallel worker processes for ``run_SGWB``.

    ``SGWB_POOL_SIZE`` (default 1) overrides the process count.  Under MPI
    (world size > 1) the default is 1 unless the variable is set explicitly:
    the reference path already parallelizes across ranks, and nested
    ``mp.Pool`` workers inside MPI ranks deadlock/oversubscribe.  This is a
    pure process-management knob; it does not change any numerical result.

    中文：该变量只控制旧 LSODA 路径的频率进程数，不改变数值结果。MPI 下默认不启用嵌套
    进程池，以免过度占用资源或发生死锁。
    """
    env = os.environ.get('SGWB_POOL_SIZE')
    if env is not None:
        n = int(env)
        if n < 1:
            raise ValueError('SGWB_POOL_SIZE must be >= 1, got %r' % env)
        return n
    return 1


def _mark_eval_status(m, label):
    """Record the engine status of one evaluation for telemetry.

    中文：记录本次求解的引擎状态，并更新该状态的累计计数。
    """
    counts = getattr(m, 'eval_status_counts', None)
    if counts is None:
        counts = {'FAST': 0, 'FAST_ESCALATED': 0, 'REFERENCE': 0,
                  'LSODA': 0, 'LSODA_FALLBACK': 0}
        m.eval_status_counts = counts
    counts[label] = counts.get(label, 0) + 1
    m.last_eval_status = label


class LCDM_SG(LCDM_SN):
    """
    Cosmological model: LCDM + stiff component + constant Delta N_eff

    The free/base parameters of the model are: 
    'Omega_bh2', 'Omega_ch2', 'H0', 'DN_eff', 'A_s', 'r', 'n_t', 'cr', 'T_re', 'DN_re', 'kappa10'.
    - [H0] = km s^-1 Mpc^-1, [T_re] = GeV.
    - Set cr > 0 if the consistency relation is assumed, otherwise set cr <= 0 and provide (n_t, DN_re).
    - DN_re is the number of e-folds from the end of inflation to the end of reheating, 
      assuming a^{-3} matter-like evolution.
    - kappa10 := rho_stiff/rho_photon at 10 MeV.

    There are three ways to instantiate a model with desired base parameters: 
    1. input a yaml file,  2. input a dictionary, 
    3. specify the parameters by keyword arguments.
    
    They will be stored in the 'obj_name.cosmo_param' attribute as a dictionary
    and can be modified later on.
    
    Once a successful 'SGWB_iter' calculation has been performed for an instance, 
    if for some reason you would like to reuse this instance with a new set of 
    parameters, you MUST run 'obj_name.reset()' first to reset the status, and 
    then modify the 'obj_name.cosmo_param' dictionary with desired values.

    中文：`LCDM_SG` 在基础宇宙学模型上提供 SGWB 迭代和引擎分派。构造参数可以来自 YAML、
    字典或关键字；复用实例并更改参数前，先调用 `reset()`。
    
    """
    def __init__(self, *args,
                 **kwargs):
        LCDM_SN.__init__(self, *args, **kwargs)
        # Evaluation telemetry is intentionally persistent across ``reset``
        # calls.  Cobaya reuses one theory instance and resets the cosmology
        # for every point; clearing counters there would make the reported
        # fallback fraction meaningless.
        # 中文：Cobaya 会复用 theory 实例；重置宇宙学参数时保留累计计数，避免回退比例被清零。
        self.fast_evals = 0
        self.fast_failures = 0
        self.lsoda_evals = 0
        self.lsoda_fallbacks = 0
        self.last_engine = None
        self.last_fast_error = None
        self.fast_failure_reason = None
        self.fast_guard_rejections = 0
        self.fast_physical_rejections = 0
        self.reference_evals = 0
        self.escalations = 0
        self.escalated_from = None
        self.reset()


    def reset(self):
        if hasattr(self, 'DN_eff_orig') and (self.DN_eff_orig is not None):
            self.cosmo_param['DN_eff'] = self.DN_eff_orig
        self.DN_eff_orig = None
        self.SGWB_converge = False     # 标记 SGWB 求解是否成功收敛。
        self.local_error_budget = None
        self.DN_gw_error = None
        self.Delta_Neff_abs_error = None
        self.spectrum_error = None
        self.quadrature_error = None
        self.error_estimates = None


    def run_SGWB_single(self, freq, z_tail=5.0, rtol=1e-6, atol=None):
        """
        Integrate the equation of motion of primordial tensor fluctuations 
        for a single frequency and given expansion history

        中文：在给定膨胀历史下，对单个频率的原初张量扰动方程进行积分。
        """

        P_t = self.Tensor_power(freq)            # 原初张量功率谱。

        j = self.find_index_hc(freq+3)           # solver begins at kc/aH <= 10^{-3}
        # 起始节点选在 kc/aH<=10^-3，即模态仍处于超视界状态时。
        z0 = (freq-self.f_hor[j]) * ln10         # convert to ln(kc/aH)

        result = solve_SGWB(self.Nv, self.sigma, j, z0,
                            z_tail=z_tail, rtol=rtol, atol=atol)

        N_this = self.Nv[(self.Nv >= result.sol.t_min) & (self.Nv <= result.sol.t_max)]
        [zf_this, xf_this, yf_this] = result.sol(N_this)
        Th_this = np.divide(yf_this, np.exp(zf_this))
        Oj_this = np.multiply(xf_this, Th_this)/3 * P_t
        Ogw_this = (np.power(xf_this,2) + np.power(yf_this,2))/24 * P_t + Oj_this
        Opgw_this = (-5*np.power(xf_this,2) + 7*np.power(yf_this,2))/72 * P_t

        # high-frequency, post-horizon-reentry regime
        # 模态重新入视界后的高频区间改用时间平均的振幅延拓。
        if (N_this[-1]<self.Nv[-1]):
            N_hf = self.Nv[self.Nv > N_this[-1]]
            f_hf = self.f_hor[self.Nv > N_this[-1]]

            coeff = math.sqrt(0.5*(xf_this[-1]**2+yf_this[-1]**2))
            Th_hf = coeff * np.exp(-zf_this[-1] + N_this[-1] - N_hf)      # the rms T_h, time-averaged
            # Th_hf 是时间平均后的 T_h 均方根振幅。
            xf_hf = np.multiply(Th_hf, np.power(10, freq-f_hf))           # the rms x_f, time-averaged
            # xf_hf 是对应的时间平均 x_f 均方根振幅。
            Oj_hf = -np.power(Th_hf,2)/3 * P_t
            Opgw_hf = np.power(xf_hf,2)/36 * P_t
            Ogw_hf = Opgw_hf * 3 + Oj_hf

            N_this = cat((N_this, N_hf), axis=None)
            Th_this = cat((Th_this, Th_hf), axis=None)
            Oj_this = cat((Oj_this, Oj_hf), axis=None)
            Ogw_this = cat((Ogw_this, Ogw_hf), axis=None)
            Opgw_this = cat((Opgw_this, Opgw_hf), axis=None)

        return [N_this, Th_this, Oj_this, Ogw_this, Opgw_this]


    def run_SGWB(self, z_tail=5.0, rtol=1e-6, atol=None):
        """
        Solve for selected frequencies and given expansion history
        using multiprocessing parallelism

        中文：在给定膨胀历史下求解选定频率；可通过多进程并行处理不同频率。
        """

        with mp.Pool(processes=_sgwb_pool_size()) as pool:
            worker = partial(self.run_SGWB_single, z_tail=z_tail, rtol=rtol, atol=atol)
            res_it = pool.imap(worker, self.f, chunksize=3)
            res = [y for y in res_it]

        # Each list entry below corresponds to one frequency channel.
        # 以下各列表的同一索引对应一个频率通道。
        self.N_hc = [single_sol[0] for single_sol in res]  # Number of e-folds, starting at the horizon crossing for each frequency
        self.Th = [single_sol[1] for single_sol in res]    # Tensor transfer function
        self.Oj = [single_sol[2] for single_sol in res]    # \dot T_h * T_h / 3H * P_t
        self.Ogw = [single_sol[3] for single_sol in res]   # dOmega_GW / dlnf
        self.Opgw = [single_sol[4] for single_sol in res]  # dOmega_pGW / dlnf


    def SGWB_iter(self, engine='fast', fallback=False, tol=1e-4,
                  z_tail=None, rtol=1e-6, atol=None, freq_res=None,
                  h=None, col_step=None, threads=None,
                  accuracy_mode=_DEFAULT_ACCURACY_MODE,
                  auto_escalate=False, error_tol=None, sigma_exact=False,
                  escalate_to_reference=False, transition_refine=None,
                  likelihood_sigma=None, dlogl_tol=None, eval_freqs=None):
        """
        Main numerical scheme:
        Iteration method that yields self-consistent cosmology including the stiff-amplified primordial SGWB,
        for which the extra radiation due to the SGWB is mimicked by a constant Delta N_eff.

        ``engine`` selects the solver: 'fast' (default; the goal-kink-hybrid fast
        solver) or 'lsoda' (the original adaptive LSODA path) or 'reference'
        (the independent continuous-sigma precision path).  When
        ``engine='fast'`` and the fast solver
        fails (returns None or raises), ``fallback=True`` automatically reruns
        the LSODA path. ``tol`` is the outer Delta N_eff self-consistency
        stopping criterion. The named fast preset uses ``1e-6``; an explicit
        non-default value overrides it. The ``1e-4`` signature default is
        retained for LSODA and legacy manual-settings paths. ``z_tail`` can
        override the selected fast preset and sets the LSODA/reference tail
        threshold; ``rtol`` and ``atol`` tune reference/LSODA ODE tolerances.
        ``freq_res`` scales frequency-grid density (audit-only; 1.0 = default).

        中文说明：默认 ``fast`` 预设的 outer tolerance 为 ``1e-6``；显式传入
        非默认 ``tol`` 会覆盖它。函数签名中的 ``1e-4`` 保留给 LSODA 和旧手动配置路径。

        Fast-path tuning (ignored by the LSODA path): when ``engine='fast'`` and
        ``accuracy_mode`` is omitted, the high-level API uses the combined ``fast``
        goal-kink-hybrid preset. Pass ``accuracy_mode=None`` explicitly for legacy manual module
        settings. ``accuracy_mode='fast'`` is the only formal user preset;
        historical ``production``/``transition_refine``/``ultra-fast`` names
        are deprecated aliases mapped to it. Explicit ``h``, ``col_step`` and
        ``threads`` plus non-default ``z_tail``/``freq_res``/``tol`` override
        the fast preset.
        With ``accuracy_mode=None``, ``h``/``col_step``/``threads`` snapshot
        the legacy settings (env FAST_H/FAST_COL_STEP/FAST_THREADS) and
        ``z_tail``/``freq_res``/``tol`` are applied as passed. The preferred
        named-mode path is per-call and does not mutate that legacy state.
        On success the model object is returned.

        """
        # 快速高层调用默认使用组合后的 fast preset；显式 None 仍保留给需要手动管理
        # 旧模块设置的兼容调用方。
        if engine == 'fast' and accuracy_mode is _DEFAULT_ACCURACY_MODE:
            accuracy_mode = 'fast'
        # None sentinel means "engine default" and is resolved to a concrete
        # value before any engine branch consumes it.  We capture the
        # explicit-override flag first so an explicitly-passed z_tail/freq_res
        # wins over a named accuracy preset (and so a default-valued 5.0/1.0 is
        # not misread as an override, see test_modes).
        # None 哨兵表示“引擎默认值”，会在进入引擎分支前解析为具体参数。先记录调用方是否
        # 显式设置，确保显式传入的 z_tail/freq_res 优先于命名精度档位；默认的 5.0/1.0
        # 不会被误判为显式覆盖（见 test_modes）。
        z_tail_explicit = z_tail is not None
        freq_res_explicit = freq_res is not None
        if z_tail is None:
            z_tail = 5.0
        if freq_res is None:
            freq_res = 1.0
        if engine == 'reference':
            # Independent continuous-sigma high-accuracy pipeline (slow; for
            # certification / benchmark points, not the MCMC thermal path).
            # 独立的连续 sigma 高精度参考流程运行较慢，用于认证/基准点，不用于 MCMC 热路径。
            from .reference import apply_reference_to_model
            apply_reference_to_model(self, freq_res=freq_res, z_tail=z_tail,
                                     rtol=rtol)
            self.last_engine = 'reference'
            _mark_eval_status(self, 'REFERENCE')
            return self
        if engine == 'fast':
            from . import fast_sgwb
            self.fast_evals = getattr(self, 'fast_evals', 0) + 1
            self.last_engine = 'fast'
            self.last_fast_error = None
            self.fast_failure_reason = None
            fallback_attempted = False
            try:
                config_overrides = {}
                if accuracy_mode is not None:
                    canonical = fast_sgwb.normalize_accuracy_mode(accuracy_mode)
                    if canonical in ('production', 'ultra-fast'):
                        warnings.warn(
                            "accuracy_mode=%r is deprecated; use the single user-facing fast mode"
                            % accuracy_mode,
                            DeprecationWarning, stacklevel=2)
                        canonical = 'fast'
                        accuracy_mode = canonical
                    cfg = fast_sgwb.ACCURACY_MODES[canonical]
                    self.accuracy_mode_used = canonical
                    self.accuracy_profile = (
                        'goal-kink-hybrid' if canonical == 'fast'
                        else 'validation')
                    f_tol = cfg['tol']
                    f_freq = cfg['freq_res']
                    tr_mode = cfg.get('transition_refine', False)
                    kink_mode = cfg.get('kink_split', False)
                    if h is not None:
                        config_overrides['h'] = h
                    if col_step is not None:
                        config_overrides['col_step'] = col_step
                    if threads is not None:
                        config_overrides['threads'] = threads
                    if z_tail_explicit:
                        config_overrides['z_tail'] = z_tail
                    if tol != 1e-4:
                        f_tol = tol
                    if freq_res_explicit:
                        f_freq = freq_res
                    config = fast_sgwb.resolve_config(canonical, **config_overrides)
                else:
                    self.accuracy_mode_used = 'none'
                    self.accuracy_profile = 'manual'
                    config_overrides = dict(h=h, col_step=col_step,
                                            threads=threads,
                                            z_tail=z_tail if z_tail_explicit else None)
                    config = fast_sgwb.resolve_config(None, **config_overrides)
                    f_tol = tol
                    f_freq = freq_res
                    tr_mode = bool(transition_refine)
                    kink_mode = False
                if transition_refine is not None:
                    tr_mode = bool(transition_refine)
                    if tr_mode:
                        kink_mode = False
                result = fast_sgwb.SGWB_iter_fast(self, tol=f_tol,
                                                  freq_res=f_freq,
                                                  sigma_exact=sigma_exact,
                                                  transition_refine=tr_mode,
                                                  kink_split=kink_mode,
                                                  config=config,
                                                  eval_freqs=eval_freqs)
            except Exception as exc:
                self.fast_failures = getattr(self, 'fast_failures', 0) + 1
                self.last_fast_error = repr(exc)
                self.fast_failure_reason = 'exception'
                if not fallback:
                    _mark_eval_status(self, 'FAST_FAILED')
                    raise
                print('SGWB_iter: fast engine failed (%s); falling back to LSODA.' % exc)
                self.lsoda_fallbacks = getattr(self, 'lsoda_fallbacks', 0) + 1
                fallback_attempted = True
                self.reset()
                result = self.SGWB_iter(engine='lsoda')
                if result is not None:
                    _mark_eval_status(self, 'LSODA_FALLBACK')
            reason = getattr(self, 'fast_failure_reason', None)
            if result is None and reason == 'shared_Neff_guard':
                # Both engines reject this deterministic non-physical point;
                # retrying LSODA would only waste time and obscure telemetry.
                # 两种引擎都会拒绝这个确定性非物理点；重试 LSODA 不会改变结果，只会增加耗时并掩盖遥测。
                self.fast_guard_rejections = getattr(self, 'fast_guard_rejections', 0) + 1
            elif result is None and reason in ('invalid_r', 'invalid_cutoff'):
                # Input/cutoff validation failures are also deterministic
                # physical rejections, not recoverable numerical failures.
                # 输入/截断验证失败同样是确定性的物理拒绝，不属于可恢复的数值失败。
                self.fast_physical_rejections = getattr(
                    self, 'fast_physical_rejections', 0) + 1
            elif result is None and fallback and not fallback_attempted:
                self.fast_failures = getattr(self, 'fast_failures', 0) + 1
                self.last_fast_error = self.last_fast_error or 'fast solver returned None'
                self.fast_failure_reason = reason or 'failed'
                print('SGWB_iter: fast engine returned no solution; falling back to LSODA.')
                self.lsoda_fallbacks = getattr(self, 'lsoda_fallbacks', 0) + 1
                fallback_attempted = True
                self.reset()
                result = self.SGWB_iter(engine='lsoda')
                if result is not None:
                    _mark_eval_status(self, 'LSODA_FALLBACK')
            elif (result is None and reason != 'shared_Neff_guard'
                  and not fallback_attempted):
                self.fast_failures = getattr(self, 'fast_failures', 0) + 1
                self.last_fast_error = self.last_fast_error or 'fast solver returned None'
            if result is None:
                if reason == 'shared_Neff_guard':
                    status = 'FAST_GUARD_REJECTED'
                elif reason in ('invalid_r', 'invalid_cutoff'):
                    status = 'FAST_PHYSICAL_REJECTED'
                elif fallback_attempted:
                    status = 'LSODA_FALLBACK_FAILED'
                else:
                    status = 'FAST_FAILED'
                _mark_eval_status(self, status)
            if result is not None:
                # Local a-posteriori error budget (physics-first, per point):
                # WKB handoff / frequency grid / quadrature / phase truncation
                # / tail / cancellation / self-consistency, computed from this
                # solve's telemetry (see fast_sgwb.estimate_local_error).
                # 根据本次求解遥测估算 WKB 交接、频率网格、求积、相位截断、尾部、相消和自洽误差。
                try:
                    self.local_error_budget = fast_sgwb.estimate_local_error(self)
                    self.DN_gw_error = self.local_error_budget['DN_gw_error']
                    self.Delta_Neff_abs_error = self.local_error_budget['Delta_Neff_abs_error']
                    self.quadrature_error = self.local_error_budget['categories']['quadrature']['value']
                except Exception:
                    self.local_error_budget = None
                if accuracy_mode is not None:
                    self.error_estimates = fast_sgwb.estimate_error(accuracy_mode)
                    self.spectrum_error = self.error_estimates['spectrum_error']
                _mark_eval_status(self, 'FAST')
                if auto_escalate:
                    # Likelihood-aware accuracy gate: escalate when the
                    # estimated |Delta logL| from the local Delta N_eff error
                    # exceeds the requested budget (sigma = effective sigma of
                    # the likelihood on Delta N_eff; dlogl_tol = budget in
                    # 似然感知门限：当 Delta_Neff 误差导致的估计 |Delta logL| 超过预算时升级精度；
                    # sigma 是似然对 Delta_Neff 的有效误差，dlogl_tol 是允许的对数似然变化。
                    budget_status = (getattr(self, 'local_error_budget', None) or {}).get(
                        'certification_status', 'uncertified')
                    do_escalate = budget_status.startswith('uncertified')
                    if likelihood_sigma is not None:
                        abs_err = float(getattr(self, 'Delta_Neff_abs_error', 0.0) or 0.0)
                        dlogl = 0.5*(abs_err/float(likelihood_sigma))**2
                        self.dlogl_estimated = dlogl
                        do_escalate = do_escalate or (
                            dlogl > (dlogl_tol if dlogl_tol is not None else 1e-3))
                    else:
                        tol_err = error_tol if error_tol is not None else 5e-3
                        do_escalate = do_escalate or (
                            float(getattr(self, 'DN_gw_error', 0.0) or 0.0) > tol_err)
                    if do_escalate:
                        self.escalated_from = accuracy_mode
                        self.escalations = getattr(self, 'escalations', 0) + 1
                        self.reset()
                        if escalate_to_reference:
                            result = self.SGWB_iter(engine='reference')
                        else:
                            result = self.SGWB_iter(engine='fast',
                                                    accuracy_mode='reference')
                        status = ('FAST_ESCALATED' if result is not None
                                  else 'FAST_ESCALATION_FAILED')
                        _mark_eval_status(self, status)
            return result
        if engine != 'lsoda':
            raise ValueError("engine must be 'lsoda' or 'fast', got %r" % (engine,))

        self.lsoda_evals = getattr(self, 'lsoda_evals', 0) + 1
        self.last_engine = 'lsoda'
        _mark_eval_status(self, 'LSODA')

        if self.cosmo_param['r'] <= 0:
            print('Must set a positive r to calculate the inflationary GWs!')
            return None

        if self.derived_param['N_inf'] is None:
            print('High-end cutoff frequency has not been set properly.')
            return None

        if self.SGWB_converge == False:
            Omega_nu = Omega_nh2/self.derived_param['h']**2

            self.DN_eff_orig = self.cosmo_param['DN_eff']
            # Record the original input value of DN_eff before entering the iterations
            # 保存输入的 DN_eff，未收敛时用于恢复原始模型状态。

            self.DN_gw = [0]; DN_gw_new = 0
            DN_gw_min = 0; DN_gw_max = 10
            converged = False
            try:
                for _iter in range(MAX_ITER):    # 外层自洽迭代。
                    self.gen_expansion()
                    self.construct_f(freq_res)
                    self.run_SGWB(z_tail=z_tail, rtol=rtol, atol=atol)
                    self.int_SGWB()

                    DN_gw_new = Neff0 * self.g2[-1] / Omega_nu

                    if not math.isfinite(DN_gw_new):
                        print('SGWB_iter: non-finite DN_gw_new, aborting.')
                        break
                    if self.DN_eff_orig + DN_gw_new > 5:
                        break

                    # Stop once the requested relative self-consistency tolerance is met.
                    # 达到指定的相对自洽容差后停止迭代。
                    if abs((Neff0+self.DN_eff_orig+DN_gw_new)/(Neff0+self.DN_eff_orig+self.DN_gw[-1]) - 1) < tol:
                        self.cosmo_param['DN_eff'] = self.DN_eff_orig + DN_gw_new
                        converged = True
                        break

                    # Use bisection to choose the next trial value when a bracket is available.
                    # 已形成有效区间时用二分法选择下一次试算值。
                    if DN_gw_new > self.DN_gw[-1] > DN_gw_min and DN_gw_max >= self.DN_gw[-1]:
                        DN_gw_min = self.DN_gw[-1]
                    elif DN_gw_new < self.DN_gw[-1] < DN_gw_max and DN_gw_min <= self.DN_gw[-1]:
                        DN_gw_max = self.DN_gw[-1]

                    if 0 < DN_gw_min <= DN_gw_max < 10:
                        DN_gw_new = (DN_gw_min + DN_gw_max)/2

                    self.cosmo_param['DN_eff'] = self.DN_eff_orig + DN_gw_new
                    self.DN_gw.append(DN_gw_new)
                else:
                    print('SGWB_iter: did not converge within %d iterations.' % MAX_ITER)
            finally:
                if not converged:
                    self.cosmo_param['DN_eff'] = self.DN_eff_orig
                    self.DN_eff_orig = None
            if not converged:
                return None


            self.SGWB_converge = True
            self.hubble = math.log10(2*math.pi) + self.f_hor + (self.Nv[-1]-self.Nv)/ln10       # log10(H/s^-1)，H = 2π f_hor/a。
            self.DN_gw = Neff0 * np.multiply(self.g2, np.exp(2*(self.f_hor-self.f_hor[-1])*ln10+2*(self.Nv-self.Nv[-1]))) / Omega_nu
            # 记录 SGWB 对 DN_eff 的完整演化，即 rho_GW(N)/(rho_gamma,0*7/8*(4/11)^(4/3))。
            # 末端满足 self.DN_gw[-1] + self.DN_eff_orig = self.cosmo_param['DN_eff']。

            self.Ogw_today = np.array([self.Ogw[i][-1] for i in range(len(self.f))])
            self.Opgw_today = np.array([self.Opgw[i][-1] for i in range(len(self.f))])
            self.Oj_today = np.array([self.Oj[i][-1] for i in range(len(self.f))])
            self.log10OmegaGW = np.log10(self.Ogw_today - self.Oj_today)  # log10(Omega_GW(f))
            # 暂不扣除超视界阶段 Omega_j 的负贡献。

            # Extra radiation (e.g. SGWB) parameterized as kappa_rad(T_i) for AlterBBN.
            # 供 AlterBBN 使用的额外辐射参数化形式，例如 SGWB 的 kappa_rad(T_i)。

            self.kappa_r = self.cosmo_param['DN_eff']* 7/8*(4/11)**(4/3) * z_ratio**4
            # Using the final asymptotic value of Delta N_eff, since for all reasonable T_re (>~ 1 MeV),
            # Delta N_eff,GW has already (or almost) reached its asymptotic value by T_i=27e9 K for AlterBBN.
            # 对合理的再加热温度（约高于 1 MeV），到 AlterBBN 的 T_i=27e9 K 时，
            # SGWB 引起的 Delta N_eff 已基本达到渐近值，因此这里使用最终值。

        return self if self.SGWB_converge else None

    def Tensor_power(self, freq):
        return self.derived_param['A_t'] * np.power((10**freq)/f_piv, self.derived_param['nt'])


    def find_index_hc(self, freq):
        """
        Find the index in the number of e-folds array at which 
        the horizon size is closest to freq, i.e., freq ~ aH/2pi

        中文：在 e-fold 数组中查找视界频率最接近 ``freq``（即 ``freq ≈ aH/(2π)``）的节点索引。
        """
        j = 0
        while freq <= self.f_hor[j] and j<len(self.N)-1:
            j = j+1
        if (j >= 1): j = j-1

        return j


    def construct_f(self, freq_res=1.0):
        """
        Construct the array of sampled frequencies to calculate tensor transfer functions,
        chosen empirically -- more points around transition!

        ``freq_res`` scales the density of the sampled frequency grid
        (0.25 .. 8.0; 1.0 reproduces the original grid).  It is an
        audit-only knob used to measure the frequency-grid error of the
        Delta N_eff integral; both the LSODA and the fast path consume it.

        中文：构造张量传递函数所用的采样频率，并在转变附近加密。``freq_res`` 将网格密度缩放到
        0.25–8.0 倍（1.0 对应原始网格）；它是审计参数，用于测量 ``Delta N_eff`` 的频率网格误差，
        LSODA 与 fast 路径都会使用。
        """
        r = float(freq_res)
        if not 0.25 <= r <= 8.0:
            raise ValueError('freq_res must be in [0.25, 8.0], got %r' % freq_res)
        fmax = self.f_hor[0]; fmin = min(self.f_hor); fcmb = math.log10(f_piv)
        if self.derived_param['nt']>0:            # 只保留超视界功率小于 1 的模式。
            fmax = min(fmax, (-math.log10(self.derived_param['A_t']))/self.derived_param['nt']+math.log10(f_piv))

        f = fmax+np.zeros(1)
        f = cat((f, f[-1]-np.linspace(1,24,int(24*r))*1e-3/r), axis=None)
        f = cat((f, f[-1]-np.linspace(1,13,int(13*r))*2e-3/r), axis=None)
        f = cat((f, f[-1]+5e-2-np.logspace(-1.28,-1, num=max(2,int(17*r)))), axis=None)
        f = cat((f, f[-1]+1e-1-np.logspace(-1,0, num=max(2,int(38*r)))[1:]), axis=None)
        f = cat((f, f[-1]-np.linspace(1,10,int(10*r))*.2/r), axis=None)

        if fmax >= self.f_re:
            f = f[np.logical_or(f>=self.f_re+.2, f>=fmax)]; f = f[f>=fmax-.4]
            f = cat((f, np.arange(f[-1]-.5, self.f_re+1.5, -.5/r)), axis=None)          # 再加热期间，覆盖 T_re 之前的频率段。
            f = cat((f, np.arange(f[-1]-.1, self.f_re+.8, -.1/r)), axis=None)          # 接近 f_re 的频率段。
            f = cat((f, np.arange(f[-1]-.05, self.f_re+.4, -.05/r)), axis=None)
            f = cat((f, np.arange(f[-1]-.01, self.f_re-.22, -.01/r)), axis=None)        # 穿过 f_re 的频率段。
            f = cat((f, np.arange(f[-1]-.022, self.f_re-.5, -.02/r)), axis=None)
            f = cat((f, np.arange(f[-1]-.05, self.f_re-1, -.05/r)), axis=None)         # 远离 f_re 后的频率段。
            f = cat((f, np.arange(f[-1]-.2, self.f_re-3, -.2/r)), axis=None)

        if self.rhostiff_re > self.rhorad_re:
            f = cat((f, np.arange(f[-1]-1, self.f_re-math.log10(self.rhostiff_re/self.rhorad_re)+1, -1/r)), axis=None)     # 再加热后、刚性物质主导结束前的频率段。
            f = cat((f, np.arange(f[-1]-.2, self.f_re-math.log10(self.rhostiff_re/self.rhorad_re)-2, -.2/r)), axis=None)   # 穿过 T_sr 对应的谱踝部。

        f = f[f>=fcmb]
        f = cat((f, np.arange(f[-1]-1, fcmb, -1/r)), axis=None)                         # 辐射主导期，从 z_eq 之前延伸至枢轴频率。
        f = cat((f, np.arange(f[-1]-.2, fmin-1, -.2/r)), axis=None); f = f[f>=fmin]     # 穿过 z_eq 并延伸至当前最低频率。

        self.f = f


    def int_SGWB(self):
        """
        Integrate the SGWB spectrum to obtain the bolometric energy fractions at each moment.
        This step is important since it calculates the Delta N_eff due to the primordial SGWB.

        中文：沿频率积分 SGWB 谱，计算每个膨胀时刻的玻尔能量分数，并据此得到原初 SGWB 对
        ``Delta N_eff`` 的贡献。
        """
        self.g2 = np.zeros_like(self.Nv)      # 总 Omega_gw。
        self.w2 = np.zeros_like(self.Nv)      # 总 Omega_pgw。

        # Patch the per-frequency solutions into matrices for integration over frequency.
        # 将逐频率求解结果写入矩阵，以便沿频率方向积分。
        M_N_hc = -np.ones((len(self.f),len(self.Nv)))
        M_Ogw = np.zeros_like(M_N_hc); M_Opgw = np.zeros_like(M_N_hc); M_Oj = np.zeros_like(M_N_hc)
        for i in range(len(self.f)):
            # Only count GW energy density when each mode is at least 3 decades inside the horizon
            # 每个模式在视界内深入至少 3 个数量级后，才计入引力波能量密度。
            cond = self.N_hc[i]>=self.N_hc[i][0]
            M_N_hc[i, (mask := np.isin(self.Nv, self.N_hc[i][cond]))] = self.N_hc[i][cond]
            M_Ogw[i,mask] = self.Ogw[i][cond]; M_Opgw[i,mask] = self.Opgw[i][cond]; M_Oj[i,mask] = self.Oj[i][cond]

        M_N_hc = np.transpose(M_N_hc); M_Ogw = np.transpose(M_Ogw); M_Opgw = np.transpose(M_Opgw); M_Oj = np.transpose(M_Oj)

        # Integrate the spectrum at each expansion-history node.
        # 对膨胀历史的每个节点积分频谱。
        for i in range(len(self.Nv)):
            ind_int = (M_N_hc[i]>=0)
            f_int = np.flip(self.f[ind_int])                           # log10(f/Hz)，按频率递增排列。
            Ogw_int = np.flip(M_Ogw[i,ind_int]); Opgw_int = np.flip(M_Opgw[i,ind_int]); Oj_int = np.flip(M_Oj[i,ind_int])

            if len(f_int)>=2:
                # Exclude the negative super-horizon Omega_j contribution from this integral for now.
                # 暂不将超视界区域 Omega_j 的负贡献计入该积分。
                self.g2[i] = integrate.simpson(Ogw_int-Oj_int, f_int) * ln10
                self.w2[i] = integrate.simpson(Opgw_int, f_int) * ln10
            else:
                self.g2[i] = 0; self.w2[i] = 0

def add_boundary_knots(spline):
    """
    Add knots infinitesimally to the left and right.

    Additional intervals are added to have zero 2nd and 3rd derivatives,
    and to maintain the first derivative from whatever boundary condition
    was selected. The spline is modified in place.

    中文：在边界两侧加入极近的节点，并补充区间，使边界处二、三阶导数为零，同时保持所选边界条件
    给出的一阶导数。该函数会原地修改样条对象。
    """
    # Determine the left-edge slope used for the extrapolation interval.
    # 求出左端斜率，供外推区间使用。
    leftx = spline.x[0]
    lefty = spline(leftx)
    leftslope = spline(leftx, nu=1)

    # Add a breakpoint immediately to the left and use the known slope
    # to construct the PPoly coefficients.
    # 在左侧紧邻处添加断点，并用已知斜率构造 PPoly 系数。
    leftxnext = np.nextafter(leftx, leftx - 1)
    leftynext = lefty + leftslope*(leftxnext - leftx)
    leftcoeffs = np.array([0, 0, leftslope, leftynext])
    spline.extend(leftcoeffs[..., None], np.r_[leftxnext])

    # Repeat the construction at the right boundary.
    # 在右端重复相同的构造过程。
    rightx = spline.x[-1]
    righty = spline(rightx)
    rightslope = spline(rightx,nu=1)
    rightxnext = np.nextafter(rightx, rightx + 1)
    rightynext = righty + rightslope * (rightxnext - rightx)
    rightcoeffs = np.array([0, 0, rightslope, rightynext])
    spline.extend(rightcoeffs[..., None], np.r_[rightxnext])
