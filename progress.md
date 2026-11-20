# Progress Log: stiffgwpy_fast 单一 fast 生产求解器

> 本日志按 session 日期记录进展；历史条目中的“当前 HEAD”不代表今天的工作分支状态。当前状态以根 README、docs 索引和最新证据 artifact 为准。

## Session: 2026-09-12 (nested native-frequency quadrature + oracle tail attribution)

### Actions Taken

- 第九原则要求的真实嵌套频率求值已落地：新增
  `scripts/benchmark_nested_frequency_quadrature.py`，流程为 base goal 网格
  -> 局部估计器 top-N 区间 -> 把区间中点并入**积分** support grid（受控替换
  `goal_oriented_freqs`，`finally` 恢复）-> 重新求解 ->
  `E_nested = |DN_refined - DN_base|`，再与同网格独立 reference 的
  `actual_error` 比较，并统计 coverage / false-safe / monotonicity。
- 实现陷阱：`eval_freqs` 只把节点加入 solve 网格、不改变积分 support grid，
  第一版实现因此读到 `E_nested = 0`；必须替换积分 support grid 才有信号。
- 阳性对照：同一加密流程在 `simpson` 下给出 `E_nested = 1.362e-09`
  （default），比 PCHIP 的 `9.07e-14` 大 4-5 个量级，证明测量对积分方法敏感。
- 归属实验（移动 oracle 而不是 fast）：`scripts/benchmark_same_grid_reference.py
  --z-tail 10` 把 Oracle A 的 frozen handoff 从 `z=8` 加深到 `z=10`，default 点
  的 fast-vs-oracle DN 相对差从 `2.944e-4` 降到 `4.956e-6`，与 reference 自报
  的 handoff 缺陷 `|1.5*sigma-1|/e^z`（`3.4e-4 -> 4.5e-5`）同阶下降；z=14
  版本在 `>20 min` 后仍未收敛而被终止（更多模式不再触发解析 handoff）。

### Test Results

- PCHIP 六点（N=4/8/12、2 线程）：`E_nested` 绝对 `1.7e-14`-`6.3e-10`、相对
  `7.4e-12`-`2.8e-9`（lowT 例外 `9.3e-5`）；同点 `actual_rel` 为
  `1.83e-4`-`2.97e-4`，比 `E_nested_rel` 大 5-7 个量级，即 native 频率网格
  已收敛，DN 残差不由网格离散/插值主导。
- 三方差值一致性（default、同网格、同一 DN_eff、z=10）：Oracle A
  `2.2636586329e-3`、Oracle C（Prüfer）`2.2636592187e-3`、fast（z_tail=5）
  `2.2636474151e-3`；两个独立 oracle 互差 `2.4e-7`，fast 距两者
  `4.96e-6`/`5.21e-6`。
- 参考运行成本：z=8 `94.2 s`、z=10 `578.8 s`（6.1 倍）、z=14 `>20 min` 未收敛。

### Artifacts

- `docs/nested_frequency_quadrature_head.json`（六点 PCHIP，新增 `commit` 字段）
- `docs/nested_frequency_quadrature_simpson_control.json`（default 阳性对照）
- `docs/oracle_a_same_grid_z10_default.json`（Oracle A `z_tail=10` 归属证据）
- `scripts/benchmark_nested_frequency_quadrature.py`
- `findings.md` 的「Nested native-frequency quadrature」与
  「Oracle A same-grid `z_tail` attribution」小节

## Session: 2026-09-12

### Actions Taken

- 从最新 HEAD `bc00a35` 重新判断：把 `scripts/profile_fast_breakdown.py` 扩展到
  六个正式代表点后，20 线程阶段分解显示 warm runtime 由 Python/背景准备层
  主导（default 点 `tensor_solve_kernel` 约 `0.81 ms`，而
  `expansion_background` `1.04 ms`、`fast_phi_s2_split` `0.63 ms`、
  `goal_frequency_construction` `0.48 ms`），其中含纯冗余计算。临时分段插桩
  定位到 `goal_oriented_freqs` 内重复的 `H2_vec`/pchip 拟合/
  `_correct_kink_background`，以及每次访问都重算的 `derived_param` property
  （约 `11 us/次`，单次求解 29 次访问）；插桩已 `git checkout` 完全回退。
- 接受四项等价冗余消除：(a) `gen_fast` 目标频率网格由
  `np.unique(np.sort(...))`（nv=13764 时约 `158 us`）改为 `concatenate` +
  稳定排序 + `not_equal` 归并去重（约 `25 us`），并缓存 `n_re_abs`；
  (b) `grid_independent_freqs.f_hor_cont` 中与 `N` 无关的
  `H2_vec(N_inf)`/`H2_vec(N_re_abs)`/`raw_last`/`raw_re`/`Delta_f`/`ln10`
  提到闭包外只算一次；(c) `Ogw`/`Oj`/`Opgw` 在 `np.zeros` 之后立即
  `fill(0.0)` 的冗余写入跳过（`_fresh_buffers` 标记）；
  (d) `_sigma_node_limits`/`fast_phi_s2_split` 内重复的 `m.derived_param`
  求值绑定为单次 `d`。
- 新增 `scripts/benchmark_fast_runtime_ab.py`：六点 warm median/p95/min 与
  `f`/`log10OmegaGW`/`DN_gw`/`g2`/`w2` 的 SHA256 digest，支持
  `--repeats`/`--threads`/`--affinity`/`--json`。

### Test Results

- 等价性：20 线程 4 轮、2 线程 3 轮 A/B 配对，六点 × 五字段 digest 全部逐位
  一致，`converged`/`n_freq`/`fast_failure_reason` 不变（648 + 288 项比较，
  另含 A/B 交叉比较，0 mismatch）。
- 速度：A,B,B,A 对称序配对（20 线程 affinity=20、15 repeats；2 线程
  affinity=20、25 repeats）default ratio `0.9456`/`0.9177`，即改善 `5.4%`/
  `8.2%`；六点中位 ratio `0.816-1.011`，无点退化超 `2%` 噪声门限。
- full pytest `159 passed, 6 deselected`；`pytest -m cobaya` `1 passed,
  164 deselected`；Ruff、mypy（2 source files）、compileall、中文注释门禁、
  `validate_manifest.py`、`git diff --check`、wheel/sdist 构建 +
  `verify_distribution` + installed-wheel smoke 全部通过。

### Artifacts

- `docs/fast_pyoverhead_ab_symmetric.json`（对称序配对汇总与 digest 校验计数）
- `docs/fast_pyoverhead_before_t20.json`、`docs/fast_pyoverhead_after_t20.json`
- `docs/fast_pyoverhead_before_t2.json`
- `scripts/benchmark_fast_runtime_ab.py`；`scripts/profile_fast_breakdown.py`
  扩展为六点 CASE 并新增 `steps_per_channel`
- `findings.md` 的「Fast Python preparation-layer de-duplication (2026-09-12)」

## Session: 2026-09-11

### Actions Taken

- Oracle C（高阶 WKB tail oracle）已完成：新增
  `scripts/benchmark_oracle_c_wkb.py`，推导并验证解析 tail 修正
  `transfer^2_wkb = transfer^2_frozen * (1 + sin(2*theta_f)/omega_f)`。
  在 `default/lowT/highT/stiff` 全 76 频率 native 网格上对比
  `frozen(z=5)`、`wkb(z=5)`、`deep(z=10)`：frozen-vs-deep DN 相对差
  `1.34e-3..3.65e-3`，wkb-vs-deep 降到 `5.46e-6..1.25e-5`（改善
  216-665 倍）；per-mode `Ogw` 残差 median 从 `2.7e-3..4.7e-3` 降到
  `2.1e-5..2.9e-5`。frozen/deep 的 DN 与 Stage B 完全一致
  （default `0.0022718753`/`0.0022636593`），验证脚本正确性。
- 结论：Stage B 的 `3.63e-3` 非单调 tail systematic 是 frozen-amplitude
  的 O(eps) 绝热缺陷（eps_med~6.7e-3..1.1e-2），不是真实物理非单调；
  解析修正把它压到与 eps² 同量级（~1e-5）。Oracle C 接受为解析 tail
  oracle。残余上界来自 deep(z=10) 自身 O(eps(10))~9e-5 的 frozen tail。

- 确认 HEAD `91598e7` 与 `fast/fast_v0.2` 一致、工作树干净、远端 CI 全绿。
- 按文档记录的 Prüfer next experiment，新增 `scripts/benchmark_prufer_fullgrid.py`
  并完成 **完整 native 网格（76 频率）认证**：`default/lowT/highT/stiff` +
  10 个参数轴 edge + 5 个 Sobol 点，`z_tail=5/7`，固定 `DN_eff` 频谱比较 +
  outer 自洽重放。
- Full-grid `DN_gw` 相对差 median `6.44e-10`、max `2.44e-9`；`Ogw/Oj/Opgw`
  p95 分别 `<=4.9e-7/1.9e-6/1.1e-5`；分量 max 异常全部位于未入尾低频模式
  （零点附近相对放大，绝对差 `1e-19..1e-22`，物理无关）；`used_tail`
  逐行一致；26 accepted outer 比较迭代次数全部一致、4 个显式
  `shared_Neff_guard`（`edge_tre_hi`/`edge_nt_blue`）双实现一致。
- 确定性重放：default 全网格 `DN_gw` 与全部频谱数组在 Prüfer/Cartesian
  两条路径均 bitwise 一致。
- 决策：PASS / VERIFIED（相对 Cartesian DOP853 reference）。Prüfer 仍为
  reference-only 原型，不晋升正式 kernel：速度信号在 Python/DOP853 栈上
  测得（8 频 0.53-0.55x，edge 含 outer 后最高 0.93x），不代表 Numba
  kernel 收益；其价值转为独立状态变量 oracle 交叉校验。

- Fast 真实 DN 误差标定（2026-09-11）：新增
  `scripts/benchmark_fast_true_error.py`，以 Oracle C 的 WKB 锚点测量 fast
  真实误差（fast 与 frozen reference 共享 `z_tail` 尾约定，fast-vs-reference
  会抵消共享 defect）。结果（vs WKB(z5)）：default `4.31e-4`、highT `7.50e-4`、
  stiff `1.29e-3`、lowT `1.24e-2`；default/highT/stiff 比此前对外报告的同约定
  PCHIP `2.94e-4` 大 1.5-4 倍，仍未进入 `2e-4` gate；lowT 为独立非 tail 误差
  源（`DN_gw=5.2e-8`），待定位。决策：ACCEPTED as calibration finding，不改
  正式 kernel；artifacts 为 `docs/fast_true_error.json` 与
  `docs/fast_true_error_assessment.md`。

- Fast 真实 DN 误差的按频率分解（2026-09-11）：新增
  `scripts/benchmark_fast_residual_decomposition.py`，用 `build_Wmat` 权重
  复现 fast 的 Simpson DN（与自报值逐位一致），并在同一积分器下比较 fast 与
  Oracle C WKB 锚点的 per-node 值：default `9.49e-06`、highT `9.32e-06`、
  stiff `1.10e-05`、lowT `1.41e-04`；而 Simpson-vs-PCHIP 积分器差为
  `4.20e-04/7.40e-04/1.28e-03/1.26e-02`（1-2 个数量级更大）。结论：真实 DN
  误差由默认 Simpson 频率积分器主导，fast 传播/tail kernel 已 ~1e-5；lowT
  残差 89% 在非 tail 且最低频单点占 52.7%。决策：ACCEPTED as diagnosis
  finding，不改正式 kernel；下一实验为默认 PCHIP 化（含 runtime 与参数空间
  验收）。

- 频率积分器 A/B（Simpson vs PCHIP，2026-09-11）：新增
  `scripts/benchmark_fast_quadrature_ab.py`（交替测量、warmup=3、repeats=25、
  Numba=2）。PCHIP 把真实 DN 误差降到 default `1.07e-05`、highT `9.88e-06`、
  stiff `1.15e-05`、lowT `1.66e-04`（vs Oracle C WKB 锚点，1-2 个数量级），
  四点全部 <2e-4 gate；但 scipy PCHIP 路径 warm runtime 比值
  `1.204..1.307`，超过 `<10%` 预算，故未切换默认。微基准定位成本为
  `PchipInterpolator.integrate` `0.141 ms/次` 与局部 estimator `1.39 ms/次`。
  决策：ACCEPTED as measurement；下一实验为预计算 PCHIP 积分权重（全局向量 +
  逐区间矩阵）并向量化 estimator，使默认切换落在预算内。

- PCHIP 单次拟合共享（2026-09-11）：先前“PCHIP 积分可预计算为固定权重向量”
  的假设被**否证**（Fritsch-Carlson 斜率是节点值的非线性函数：单位矩阵探针
  给出的 `weights @ y` 与 `PchipInterpolator(x, y).integrate()` 相差约 20 倍）。
  改为 `pchip_integral_breakdown` 一次构造 `PchipInterpolator` +
  `antiderivative()`，同时返回全程与逐区间积分，并让
  `estimate_frequency_quadrature_local(..., candidate_intervals=...)` 复用；
  非重叠三点 Simpson 基线改为等价向量化（实测逐位一致，`0.674 -> 0.027 ms`）。
  同会话配对 A/B（`--repeats 50`，前后各一次、NUMBA=2、BLAS=1）：PCHIP 专属
  开销 `2.20/1.92/2.01/2.29 ms -> 0.98/0.59/0.72/0.84 ms`，比值
  `1.366/1.307/1.228/1.240 -> 1.165/1.095/1.073/1.089`
  （default/lowT/highT/stiff）；四点 `DN_gw` 与改前逐位一致（rel `0.0`），
  vs Oracle C WKB 不变。决策：ACCEPTED（opt-in PCHIP 路径严格 Pareto 改善，
  同等观测值、约 60% 更少 PCHIP 专属开销）；**REJECTED：本阶段不切换默认
  `frequency_quadrature`**（预登记 `<10%` 未稳健满足：同会话 default `1.165`，
  16 线程生产刻度探针 `1.101-1.154`）。下一实验：NumPy 向量化 PCHIP 斜率 +
  分段解析积分核 + estimator 分摊循环向量化，再评估默认切换。

- 向量化 PCHIP 积分核与 estimator 分摊（2026-09-11）：新增
  `_pchip_integrals_vectorized`（`_pchip_slopes` / `_pchip_edge_slope` 复刻
  scipy 的 Fritsch-Carlson 斜率与形状保持 guard，分段用 Hermite 闭式积分
  `h/2*(y_i+y_{i+1}) + h^2*(m_i-m_{i+1})/12`），`pchip_integral_breakdown`
  保持 scipy 参考实现；局部 estimator 的逐面板分摊改为等价向量化（三种
  allocation 逐位一致，含零候选和退化为半误差均分）。微基准（76 点 native
  网格、median of 2000）：scipy breakdown `0.143 ms` -> 向量化 `0.024 ms`；
  estimator（复用候选）`0.290 ms -> 0.061 ms`。配对 A/B（`--repeats 50`、
  NUMBA=2 与 16、前后各一次）：2 线程比值
  `1.1622/1.1261/1.0495/1.0901 -> 1.0105/1.0551/1.0397/1.0411`，16 线程
  `1.1987/1.1137/1.1270/1.1203 -> 1.0428/1.0342/1.0180/1.0413`；PCHIP 专属
  开销 `0.06-0.42 ms`（2 线程）/`0.15-0.33 ms`（16 线程）；同会话第二次配对
  复现同一模式。实测 `DN_gw` 相对 scipy PCHIP 路径变化 `3.83e-16/4.98e-16/0/0`
  （1-2 ulp），默认 Simpson 路径逐位不变。决策：ACCEPTED（预算 `<1.10` 在两个
  线程刻度均满足）。下一步（已写定 acceptance）：默认切换 PCHIP + manifest/
  README/ERROR_BUDGET/coverage 同步刷新。

- PCHIP 默认切换落地（2026-09-11）：`_SGWB_iter_fast_impl` 与
  `SGWB_iter_fast` 的默认 `frequency_quadrature` 改为 `pchip`（`simpson` 仍可
  显式选择）；PCHIP 热路径遇到非有限 integrand 时改为把 NaN 交给统一的
  `math.isfinite(DN_gw_new)` guard，保持既有 `fast_failure_reason='nonfinite'`
  的 abort/restore 语义，`_pchip_integrals_vectorized` 继续 fail-loud。权威
  配对 A/B（50 repeats、workers=1、BLAS=1，`docs/fast_quadrature_default_ab_t{2,16}.json`）
  的 `pchip/simpson` warm median 比值：2 线程 `1.019/1.035/1.034/1.024`、
  16 线程 `1.069/1.028/0.987/1.017`（default/lowT/highT/stiff），全 `<1.10`。
  20 线程 25-repeat 矩阵 `docs/benchmark_head_matrix.json`：default warm
  median/p95 `4.932/5.467 ms`、cold `0.222 s`；highT/stiff/high-kappa 仍是
  最慢工况 `7.79/7.79/7.27 ms`。stability screen `guard_count=3`、
  `failure_count=0`；LHS 参数空间 400 点 `254 success / 146 physical_guard /
  0 numerical_failure`。`ERROR_BUDGET` 的频率积分项由 `1.0e-3` 重新标定为
  `2.0e-4`，README/README_zh/CHANGELOG/fast_v02_audit_report 同步刷新。
  决策：ACCEPTED，完整评测见 `docs/fast_quadrature_default_switch_assessment.md`。
  下一步：high-T/stiff/high-kappa runtime 分解与 nested native-frequency
  quadrature（真实新增频率 solve）以区分 same-grid 残差是积分还是传递主导。

### Artifacts

- `docs/oracle_prufer_fullgrid_{default,lowT,highT,stiff}.json`
- `docs/oracle_prufer_fullgrid_edge_{r,tre,dnre,kap,nt}_*.json`
- `docs/oracle_prufer_fullgrid_sobol_{000,002,006,010,015}.json`
- `scripts/benchmark_prufer_fullgrid.py`
- `docs/fast_true_error.json`
- `docs/fast_true_error_assessment.md`
- `scripts/benchmark_fast_true_error.py`
- `docs/fast_residual_decomposition.json`
- `docs/fast_residual_decomposition_assessment.md`
- `scripts/benchmark_fast_residual_decomposition.py`
- `docs/fast_quadrature_ab.json`
- `docs/fast_quadrature_ab_assessment.md`
- `scripts/benchmark_fast_quadrature_ab.py`
- `docs/fast_quadrature_reuse_ab.json`
- `docs/fast_quadrature_reuse_ab_before.json`
- `docs/fast_quadrature_reuse_assessment.md`
- `tests/test_pchip_integral_breakdown.py`
- `docs/fast_quadrature_pchip_kernel_t2_before.json`
- `docs/fast_quadrature_pchip_kernel_t2_after.json`
- `docs/fast_quadrature_pchip_kernel_t16_before.json`
- `docs/fast_quadrature_pchip_kernel_t16_after.json`

## Session: 2026-09-09

### Current Status

- **Phase:** 3 - Implementation and evidence-driven optimization
- **Started:** 2026-09-09

### Actions Taken

- 读取目标文件，确认最终目标为单一 fast 生产 solver，并保留 reference oracle。
- 读取用户指定的 `pasted-text-1.txt`；确认本轮优先级为远端 HEAD fresh profiling、DN_gw 误差分解、嵌入式 frequency quadrature 与可验证 Pareto benchmark。
- 核对当前工作树：分支为 `codex/fast_v0.2`，HEAD 与 `fast/fast_v0.2` 均为 `767056d`；未发现代码工作树改动。
- 读取并启用 `using-superpowers`、`self-evolution`、`planning-with-files`、`brainstorming` 与 `writing-plans` 的适用流程。
- 确认分支 `codex/fast_v0.2`、远端目标 `fast/fast_v0.2`；重大修改按阶段提交并推送。
- 完成首轮关键词和入口检索，发现 README、Cobaya、benchmarks、tests 仍广泛暴露双档位。
- 核对远端：`fast` 指向目标仓库；用户要求重大 Phase 提交并推送远端同名分支。
- Phase A 已加入 uniform-grid 单区间 transfer split 与 one-sided reheating primitive；目标测试覆盖断点属性和单侧积分约定。
- Phase A 的 sparse subset 结果已降级为诊断；benchmark 现默认在候选完整 246 点频率网格上运行独立 reference，避免把 subset 频率积分误报成总 `DN_gw` 误差。
- Phase B 已完成 split transfer segment 的 phase-cap 实现；`phase_max=.25` target test 通过，但单点 oracle 收益有限。
- Phase C 已接入 `freq_grid='goal'`：seed + reheating feature reserve + native `eval_freqs`；测试和 solver native-node 验证通过。
- 正式 fast preset 已切换为单一 `goal-kink-hybrid` 组合；旧 production 保留作 validation/compatibility mode。
- CI #53 的根因是新增英文代码注释未满足中文注释门禁；`4436b9d` 已修复并推送，CI run `34304423555` 的 3.9–3.13 与 Cobaya jobs 全部通过。
- seed64 频率预算已完成独立 reference A/B：76 点比 seed80 的 90 点更省节点，默认点 spectrum p95 `2.719e-3`、DN rel `1.368e-3`；保留 seed64 作为当前实验基线。
- 四阶 Magnus 与背景曲率局部子步实验均回退：没有把 `h=.01` 的 spectrum p95 压到 `<3e-3`，且引入额外热路径成本。
- 新增安全热路径优化：外层首轮 `solve_kernel` 只组装最终列，后续收敛轮再组装完整列；新增测试通过，like-for-like plain warm median 约 `6.727 -> 5.051 ms`，正式 kink 路径尚未达到 `<=4 ms/point`。
- 节点 sigma 复用实验虽数值等价，但没有稳定耗时收益，已回退，不纳入主线。
- exact split 的节点 sigma 缓存已通过逐项等价测试并接入 formal kink 路径；profile steady warm median 约 `10.62 -> 8.90 ms`。外层审计显示第二轮 DN 更新很小，但一次 full solve 会让最终 `g2/w2` 在 high-T 点偏离约 `3%`，所以只保留首轮 probe 优化，不跳过第二次 full solve。
- 批量合并 midpoint/quarter sigma 的尝试数值等价但变慢（约 `8.4 -> 10.6 ms`），已回退；线程扫描也未找到稳定的 `<=4 ms` 配置，8/16 线程稳定 profiler 均约 `8.4 ms`。
- formal kink 路径已跳过不会被使用的 uniform-grid primitive preparation，只保留 frequency-only preparation；新增回归测试通过，10 次 profiler 与 full-prep 对照的 steady median 约 `9.68 -> 8.98 ms`，输出 digest 完全一致。
- 在节点 sigma/f_hor/Nv 变化均小于 `1e-4` 时复用上一轮 exact primitive，default 的 spectrum 最大变化约 `2.98e-4 dex`、DN 变化约 `1.1e-12`，高温点自动回退；formal default warm median 约 `5.96 ms`。独立 reference 显示当前 DN rel 仍为约 `4.04e-3`，需要继续解决积分量误差。
- z_tail 从 5 加深到 6 的对照没有改善 DN（约 `4.04e-3`），只带来很小的频谱变化，已回退为 z_tail=5。
- 新增 `fast_phi_s2_split`：平滑背景区复用节点 sigma，只有 reheating kink 区间保留连续-sigma 单侧探针；默认点与原 primitive 的 DN 差 `3.22e-12`、频谱最大差 `2.80e-6 dex`。
- 新 primitive 的默认点深尾 reference（76 点、z_tail=8、rtol=1e-11）为 spectrum p50 `1.3968e-4`、p95 `4.1466e-4`、max `1.4152e-3`，DN rel `7.1404e-4`；完整回归 `116 passed, 6 deselected`。
- 20-thread warm 速度中位数约 `4.49 ms/point`，接近但尚未满足 `<=4 ms/point`，下一步继续分析 kernel/outer probe。
- 本轮新增正式 goal+kink 路径的 outer full-solve 稳定性门控：首轮完整组装；更新后 `sigma`、`f_hor` 均在 `1e-4` 内则复用首轮结果，否则保留第二次完整 kernel。默认点由两次 kernel 降为一次，high-T/stiff 自动回退；完整回归 `117 passed, 6 deselected`，16-thread profiler warm median 约 `5.37 ms`，仍未满足 `<=4 ms/point`。
- 本轮静态验证：ruff、mypy、manifest、中文注释门禁全部通过。临时 A/B 脚本和未验证 tail 扫描脚本均已删除，未纳入分支。
- 2026-09-10 已 fetch 并确认远端 HEAD `767056d`；新增固定 affinity/threading-layer profiler，完成 1/2/4/8/16/20/32 threads 的 fresh warm/cold 分层测量。20 threads reuse A/B 显示 reuse `5.006 ms`/1 kernel，禁用后 `5.705 ms`/2 kernels。
- 2026-09-10 完成 default fresh full reference：reference `DN_gw=0.002262832966946746`，formal fast `0.002263022064729132`，scalar relative error `8.36e-5`。
- 按 TDD 新增 opt-in `frequency_quadrature='pchip'` 候选及 `integrate_frequency_pchip` helper；同一 spectrum 的 default/lowT/highT/stiff quadrature deltas 为 `4.20e-4/1.275e-2/7.40e-4/1.282e-3`，暂不切换默认。
- 继续 fresh 审计：完成六点 25-repeat standardized matrix；完成六种频率积分表示比较；完成 DN-driven midpoint refinement 原型。原型在 default 的 76→80→86→96 nodes 上 PCHIP DN 误差非单调，拒绝升级。

### Test Results

| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| `git branch --show-current` | `codex/fast_v0.2` | `codex/fast_v0.2` | PASS |
| `git diff --stat` | 初始无代码 diff | 无输出 | PASS |
| `git remote -v` | 目标 `fast` 远端存在 | `https://github.com/jiuluo-if/stiffGWpy_fast.git` | PASS |
| `git ls-remote --heads fast` | 远端同名目标存在 | `fast_v0.2`=`4436b9d`，`main`=`b83aa89` | PASS |
| `python -m pytest -q` | 现有回归通过 | `110 passed, 6 deselected` | PASS |
| `python scripts/bench_fast.py --reps 3 --cases 0` | fresh fast runtime 基线 | warm median `2.461 ms`，fallback `0` | PASS |
| `python scripts/benchmark_reference.py --point default --z-tail 5 --no-ode-error --no-tail-error` | fresh oracle 对照 | production fast vs reference 已输出，误差仍超目标 | PASS |
| `python -m pytest -q tests/test_freq_adaptive.py::test_split_primitive_uses_left_limit_at_reheating_kink tests/test_freq_adaptive.py::test_breakpoint_phi_s2_is_accurate_without_dense_subgrid tests/test_fast_sgwb.py::test_kink_split_variant_inserts_only_reheating_breakpoint` | Phase A target tests | 2 passed, 1 passed after endpoint correction; current rerun pending | PASS |
| `python scripts/benchmark_phase_a.py --reps 3 --subset-reference` | Phase A sparse diagnostic | plain `7.389 ms`, kink `8.305 ms`; transition max `1.281% -> 0.776%`; total DN comparison invalid by design | DIAGNOSTIC |
| `python scripts/benchmark_phase_a.py --reps 3` | Phase A full-grid oracle | plain `7.713 ms`, kink `10.020 ms`; DN rel `1.069% -> 0.411%`; spectrum dex p95 `1.004e-2 -> 4.745e-3`; transition max `3.066% -> 1.317%` | PARTIAL |
| `python scripts/profile_fast_breakdown.py --case A --reps 3` and `--kink-split` | Phase A profiler | steady median total `3.701 ms` plain / `7.441 ms` kink at 16 threads; kink exact primitive overhead is measurable | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 3 --freq-grid goal --phase-max 0.25` | Phase B/C at `h=.02` | 76-point grid; kink spectrum dex p95 `4.792e-3`, DN rel `4.344e-3` | PARTIAL |
| `python scripts/benchmark_phase_a.py --reps 2 --h 0.005 --freq-grid goal --phase-max 0.25` | Phase B/C precision probe | seed64 76-point grid; kink spectrum p95 `2.719e-3`, max `3.040e-3`; DN rel `1.368e-3`; thread=1 median `15.030 ms` | PARTIAL |
| `FAST_THREADS=16 python scripts/bench_fast.py --reps 3 --cases 0` | Formal fast runtime probe | warm median `6.727 ms`, cold `1.919 s`, fallback `0`; speed target not met | PARTIAL |
| `python -m pytest -q` after CI fix and experiment rollback | full regression | `111 passed, 6 deselected in 25.15s` | PASS |
| `python scripts/bench_fast.py --reps 3 --cases 0` with `FAST_THREADS=16` | outer-probe hot-path comparison | plain warm median `5.051 ms`; formal kink probe `8.872 ms`; fallback `0` | PARTIAL |
| scoped ruff + mypy + Chinese comment gate + manifest | CI-maintained checks | all PASS; comments base `93645ec` | PASS |
| `python scripts/run_fixed_profile.py --threads {1,2,4,8,16,20,32} --affinity-count {20,32} --reps 7 --case A --kink-split` | HEAD fresh fixed-environment profiling | warm median `9.78/7.28/6.66/5.04/5.88/5.47/4.81 ms`; digest stable; workqueue | PASS |
| same profiler at 20 threads with `--disable-outer-reuse` | reuse A/B | enabled `5.006 ms`, 1 kernel; disabled `5.705 ms`, 2 kernels | PASS |
| `python scripts/benchmark_reference.py --point default --freq-full --z-tail 8 --no-ode-error --no-tail-error` | fresh full independent reference | reference runtime `652.44 s`, `DN_gw=0.002262832966946746`, `n_freq=242` | PASS |
| `python -m pytest -q` after opt-in quadrature candidate | regression safety | `119 passed, 6 deselected`, 2 pre-existing deprecation warnings | PASS |
| `python scripts/benchmark_head_matrix.py` | 六点 25-repeat runtime/状态基线 | workqueue、20 threads、CPU 0-19；0 failure；default warm median/p95 `4.510/5.056 ms` | PASS |
| `python scripts/compare_frequency_quadrature.py` | 同一 76-node spectrum 的积分表示 A/B | Chebyshev 及固定 PCHIP/cubic 在参数空间出现明显非一致偏差 | DIAGNOSTIC |
| `python scripts/benchmark_dn_refinement.py` | DN-driven midpoint estimator | 76→80→86→96 的 PCHIP DN 相对 dense reference 非单调 | REJECTED |
| `python scripts/benchmark_same_grid_reference.py --point default` | exact native-grid independent reference | 76 nodes；reference PCHIP DN `0.0022643136483710326`；Simpson/PCHIP 相对误差 `7.14e-4/2.94e-4`；reference runtime `98.90 s` | PASS / DIAGNOSTIC |
| `python scripts/benchmark_eval_invariant.py` | 六点 eval-only native-node invariant | 4 个 eval nodes 加入后六点 `DN_gw` delta 全为 0；无 numerical failure | PASS |
| `python scripts/benchmark_node_counts.py` | exact 76/80/90/110-node sweep | DN 相对 dense reference `6.01e-5/5.55e-4/1.40e-4/2.29e-4`；收敛非单调；无 numerical failure | PASS / REJECTED FOR PROMOTION |
| focused tests for embedded telemetry and eval invariant | production behavior guard | `2 passed` | PASS |
| final-H HEAD artifact regeneration | release evidence freshness | quadrature/seed/refinement/eval/node artifacts regenerated from `9bb2697`; same-grid reference regenerated with same source HEAD | PASS / DELIVERY PENDING |
| `python scripts/benchmark_fast_stability.py` | fast-only named + Sobol stability screen | 24 points；0 numerical failure；3 explicit `shared_Neff_guard`；embedded estimator telemetry recorded | PASS / DIAGNOSTIC |
| `python scripts/benchmark_same_grid_reference.py --point lowT/highT` | fresh same-grid oracle for parameter dependence | low-T PCHIP `1.826e-4`；high-T PCHIP `2.932e-4`；PCHIP not globally promoted | PASS / DIAGNOSTIC |
| final default telemetry 50-repeat profiler | remove default PCHIP construction | warm median `4.976 -> 4.544 ms` (`8.7%`), digest unchanged | PASS / ACCEPTED |
| `python scripts/build_two_mode_manifest.py; python scripts/validate_manifest.py` | estimator coverage replay | Simpson coverage `0/6`、PCHIP coverage `6/6`；release gate `NOT VERIFIED`；manifest validation OK | PASS / NOT VERIFIED |
| `python scripts/benchmark_candidate_grid.py --seed 78` + default independent oracle | 89-node/PCHIP candidate A/B | 六工况无 numerical failure；default same-grid DN rel `1.142e-3` Simpson / `2.952e-4` PCHIP；candidate rejected | PASS / REJECTED |
| current-HEAD `error_budget_probe.py` default/lowT/highT/stiff | outer reuse safety | false-safe `0/4` under DN `2e-4` and spectrum-max `1e-3`; default spectrum max delta `9.03e-4`，保留但标注接近预算 | PASS / ACCEPTED WITH LIMIT |
| current-HEAD `benchmark_fast_stability.py` | refresh formal-fast SHA-bound stability artifact | 24 points；0 numerical failure；3 explicit `shared_Neff_guard`；artifact commit=`6a32838` | PASS / DIAGNOSTIC |
| `benchmark_same_grid_reference.py` on cr0/tilt/Sobol extras | parameter-space independent oracle | 5 additional points；no numerical failure；PCHIP still `2.76e-4..2.98e-4`；manifest gate `NOT VERIFIED` | PASS / NOT VERIFIED |
| `benchmark_phase_candidate.py --phase-max 0.25 0.35 0.5` | phase-cap speed/accuracy A/B | 25 repeats x 6 cases；no stable >5% gain；larger caps rejected；spectrum/DN perturbation recorded | PASS / REJECTED |
| `solve_kernel` Numba literal-specialization probe | kernel branch-elimination candidate | 42 focused tests pass but warm kernel regresses to ~150 ms；source reverted，candidate rejected | PASS / REJECTED |
| default same-grid reference `z_tail=5` vs `8` | oracle-independence tail A/B | DN `0.0022718753` vs `0.0022643136`，relative delta `3.34e-3`；记录为 tail/transfer caveat | PASS / DIAGNOSTIC |

| Q1 helper TDD red-green | 新增统一频率积分 API 的测试首次 `2 failed`（API 不存在），实现后 `2 passed` | PASS |
| Q1 solver candidate 入口 | `gauss3` opt-in 与既有 `pchip` 测试 `2 passed` | PASS |
| Q1 相关回归 | `python -m pytest -q tests/test_fast_sgwb.py tests/test_modes.py tests/test_validation_manifest.py` 得 `46 passed, 1 deselected` | PASS |
| 统一 Q1 比较脚本 | HEAD `bd648fc`、20 threads/workqueue、76/77 native nodes；6 工况均无 failure；尚未与 dense reference 比较 | PASS / DIAGNOSTIC |
| Q2 local estimator TDD | 初次实现暴露 Chebyshev spline 未初始化及全局权重分摊广播错误；均已通过回归测试修复 | PASS |
| default dense-reference 八方法比较 | 独立 reference `103.6 s`；PCHIP/Gauss 全局 DN rel `2.942e-4`，natural cubic `9.882e-4`，log-PCHIP `1.567e-3`；局部绝对和 estimator 约 `15.3%`，严重过保守，暂不用于 adaptive production | PASS / REJECTED FOR PRODUCTION |
| panel estimator 修正后 default dense-reference | 独立 reference `99.9 s`；PCHIP/Gauss 实际 DN rel `2.942e-4`，local-sum `3.118e-3`，local-max `5.615e-4`；两者均保守但 sum 过宽，max 作为后续 coverage 候选 | PASS / DIAGNOSTIC |
| Q2 coverage replay | 9 个点（default/lowT/highT/stiff/cr0/tilt/3 Sobol）；local-sum 与 local-max 均 `9/9` coverage、`0` false-safe；local-max 中位预测/实际 `2.01x`，p99 actual/predicted `0.524`，low-T 最大过保守 | PASS / NOT PRODUCTION-READY |
| low-T dense-reference 八方法 | Simpson DN rel `1.241e-2`；PCHIP/Gauss/Chebyshev `1.826e-4`；natural cubic `5.901e-3`；log-PCHIP `9.706e-3`；证实过保守 estimator 主要反映 Simpson 基线误差 | PASS / DIAGNOSTIC |
| high-T dense-reference 八方法 | Simpson `1.033e-3`；PCHIP/Gauss/Chebyshev `2.932e-4`；natural cubic `1.248e-3`；log-PCHIP `1.799e-3`；local-max `5.615e-4` | PASS / DIAGNOSTIC |
| stiff dense-reference 八方法 | Simpson `1.576e-3`；PCHIP/Gauss/Chebyshev `2.965e-4`；natural cubic `8.184e-4`；log-PCHIP `1.389e-3`；local-max `5.962e-4` | PASS / DIAGNOSTIC |
| cr0 blue dense-reference 八方法 | Simpson `8.830e-4`；PCHIP/Gauss/Chebyshev `2.765e-4`；natural cubic `5.353e-4`；log-PCHIP `1.380e-3`；local-max `7.727e-4` | PASS / DIAGNOSTIC |
| Sobol-000 dense-reference 八方法 | Simpson `1.048e-3`；PCHIP/Gauss/Chebyshev `2.964e-4`；natural cubic `7.130e-4`；log-PCHIP `1.313e-3`；local-max `7.724e-4` | PASS / DIAGNOSTIC |
| Sobol-002 dense-reference 八方法 | Simpson `9.784e-4`；PCHIP/Gauss/Chebyshev `2.934e-4`；natural cubic `8.234e-4`；log-PCHIP `1.411e-3`；local-max `7.715e-4` | PASS / DIAGNOSTIC |
| Sobol-006 dense-reference 八方法 | Simpson `8.558e-4`；PCHIP/Gauss/Chebyshev `2.978e-4`；natural cubic `5.918e-4`；log-PCHIP `1.118e-3`；local-max `5.753e-4` | PASS / DIAGNOSTIC |
| positive-tilt 配置修正与 dense-reference | 原 `cr=1` 实际覆盖 `n_t`，已改为 `cr=0,n_t=.2,DN_re=5`；修正后 PCHIP/Gauss/Chebyshev `2.927e-4`，natural cubic `7.054e-4`，log-PCHIP `1.548e-3`；local coverage `50%` | PASS / REJECTED ESTIMATOR |
| final Q2 regression | coverage replay refreshed；full pytest `126 passed, 6 deselected, 2 warnings in 100.14s`；compileall/diff check pass | PASS |
| final Q2 regression refresh | positive-tilt 修正后 full pytest `127 passed, 6 deselected, 2 warnings in 98.64s`；compileall、scoped ruff 与 diff check pass | PASS |
| 76/80/90/110 node same-spectrum sweep | PCHIP DN `2.263647e-3/2.263064e-3/2.263324e-3/2.262561e-3`；local-max `5.615e-4/7.715e-4/7.711e-4/5.616e-4`；节点增加与 estimator 均非单调，不放行 blind refinement | PASS / REJECTED FOR PROMOTION |
| default local-reference estimator 对齐 | local Pearson `r=0.592`、coverage `66.7%`、false-safe `33.3%`；predicted/actual local sum `1.353e-3/1.315e-4`；证明当前 panel estimator 不合格，global total coverage 不足以放行 | PASS / REJECTED |
| Q1 warm integration overhead/shape screen | 50 warm repeats x 6 工况：Simpson `43.9 us`、PCHIP `120.9 us`、natural cubic `105.5 us`、Gauss5 `577.6 us`、log-PCHIP `630.4 us`、Chebyshev `4727.2 us`；PCHIP/Gauss/Chebyshev shape-preserving，natural cubic 多工况 overshoot | PASS / DIAGNOSTIC |
| Sobol-002/006 dense-reference 八方法 | Sobol-002 PCHIP/Gauss/Chebyshev `2.934e-4`；Sobol-006 `2.978e-4`；均未达到 `<2e-4` | PASS / DIAGNOSTIC |
| full regression after Q2 changes | `126 passed, 6 deselected, 2 warnings in 104.43s` | PASS |
| final static verification | `compileall`、scoped ruff (`E501,F401,F821,F841`) 与 `git diff --check` 全部通过 | PASS |
| quadrature-only dense local reference | 固定同一 native spectrum、PCHIP dense Simpson 作为 quadrature-only 真值；9 点上 PCHIP-Simpson interval-max coverage `9/9`，sum coverage `8/9`；max 中位 prediction/actual `1.263x`，仍仅作诊断 | PASS / DIAGNOSTIC |
| Gauss-Simpson quadrature-only replay | `gauss2/3/5` interval-max coverage 均 `9/9`，sum coverage 均 `7/9`；与 PCHIP 相同的 panel 参考下没有校准或速度优势，不晋升 production | PASS / REJECTED FOR PROMOTION |
| Q2 independent-local actual 对齐修正 | 将 actual 改为 formal Simpson panel 对 independent reference 的误差；default coverage `36%`、positive-tilt `28.4%`，虽 Pearson `0.965/0.980`，99% 所需 safety factor `707/1331`，确认当前 interval 分配不合格 | PASS / REJECTED |
| Q2 panel-envelope allocation | 将完整 panel 误差及相邻一个 panel 的最大值传播到 interval；default coverage `92%`、positive-tilt `86.5%`，99% safety factor `2.62/6.63`，仍未达到 95/99% gate，不进入 production | PASS / REJECTED |
| panel-envelope regression | 新 allocation、independent-reference benchmark 与 coverage replay 完成；full pytest `131 passed, 6 deselected, 2 warnings in 104.80s`；compileall、scoped ruff、diff check pass | PASS |
| Q2 ensemble independent replay | PCHIP/Gauss/log-PCHIP/natural-cubic/Chebyshev 的 panel-envelope ensemble：default `100%`、lowT/highT `97.3%`、stiff `94.7%`、cr0-blue `89.0%`、positive-tilt `94.6%`、Sobol-000 `92%`、Sobol-002 `100%`、Sobol-006 `94.6%`；仍未满足全参数 95/99% gate | PASS / REJECTED |
| Phase 1 CI/package audit | 修复仓库级 Ruff `I001`；本机镜像 403 与旧 setuptools 环境问题已区分，使用官方 PyPI 完成 isolated sdist/wheel、distribution boundary、installed-wheel smoke；full pytest `132 passed, 6 deselected, 2 warnings`、Cobaya `1 passed, 137 deselected`、Ruff/mypy/manifest/comment gate 全部通过 | PASS |

| Prüfer reheating/Sobol outer follow-up | `edge_tre_lo`、`sobol_000/002/006` 在 z=5/7 均为 matching outer iterations，最大 `DN_gw` relative difference `8.33e-10/2.12e-9/2.06e-9/2.50e-9`；`edge_tre_hi` 两个深度均为 Prüfer/Cartesian 一致的 physical guard；新增 artifacts 记录资源 telemetry 与 guard status | PASS / STANDALONE ONLY |
| Prüfer complete edge/Sobol replay | 10 个参数轴 edge + `sobol_000/002/006/010/015`，共 15 点、26 个 accepted outer comparisons、4 个 physical guards、0 numerical failure；accepted 最大 outer DN 差 `6.44e-9`、power 差 `3.26e-7`、runtime ratio `0.929` | PASS / STANDALONE ONLY |

| Oracle tail convergence scaffolding | 新增 `reference.summarize_tail_convergence`，显式报告 deepest-tail central value、observed systematic bound、相邻变化和仅描述性的指数衰减率；`run_reference` 新增可选 `workers` 传递；新增 `scripts/benchmark_oracle_tail_convergence.py` 固定 fast `DN_eff` 扫描 `z_tail=5/6/7/8/10`。完整 native-grid default 扫描因 reference 计算过慢停止，未生成不完整 artifact，尚无新的 tail 数值结论 | IMPLEMENTED / NUMERICAL SWEEP PENDING |
| Oracle scaffolding verification | focused reference tests `6 passed, 1 deselected`；full pytest `133 passed, 6 deselected, 2 warnings`；Cobaya `1 passed, 138 deselected`；maintained Ruff、mypy、manifest、compileall/diff check、官方 PyPI isolated build、distribution boundary、installed-wheel smoke 全部通过 | PASS |
| Oracle tail representative-subset sweep | 固定 fast `DN_eff`、reference `z_tail=5/6/7/8/10`；default 12 点子集 observed systematic rel `2.520e-3` 且非单调；`[-4,2]` 信号区间 4 点子集：low-T `6.345e-3`、stiff `2.396e-3`，均非单调；只作为 oracle-floor/tail-sensitive 诊断，不能替代完整 native-grid 认证 | PASS / ORACLE-SENSITIVE / FULL GRID PENDING |
| Tail diagnostic scope simplification | 后续研究测试先采用明确标注的代表频率子集和 focused regression，完整 pytest/Cobaya/CI 仅在代码提交或正式门禁时运行；不得把子集结果写成全频率精度结论 | ACCEPTED PROCESS CHANGE |
| Phase A/B/C test and resource audit | 新增 `docs/test_coverage_matrix.md` 与 3-case compatibility smoke；CI 改为 compatibility 3.9–3.13、canonical 3.11、static、package、Cobaya 五职责 job；slow workflow 与本地 oracle 默认 Numba=2、BLAS=1、reference workers=1；高并行仅显式指定 | IMPLEMENTED / LOCAL VERIFIED |
| Resource-capped canonical verification | 在 `NUMBA_NUM_THREADS=2`、`FAST_THREADS=2`、BLAS=1 下：compatibility `3 passed`；canonical `136 passed, 6 deselected, 2 warnings`；Cobaya `1 passed, 141 deselected`；Ruff/mypy/manifest/compileall/diff/build/distribution/wheel smoke 全部通过 | PASS |
| Oracle checkpoint/resume hardening | tail sweep 绑定当前 commit、reference module SHA、参数、频率、`DN_eff`、`z_tail`、`rtol` 的 immutable cache key；每个 z_tail 原子写 checkpoint，schema/commit/version 不匹配时拒绝恢复；新增纯逻辑 checkpoint 测试 | PASS / STAGE B READY |
| Oracle Stage B default full native grid | 当前 HEAD `886c64f`、76 native frequencies、`z_tail=5/6/7/8/10`、workers=1、Numba=2；runtime `8.34/16.52/40.02/101.82/708.23 s`；central `DN_gw=0.0022636593`；observed systematic rel `3.6295e-3`，tail sequence 非单调；resume 命中 `5/5`，结果保存于 `docs/oracle_tail_convergence_default_full.json` | PASS / ORACLE-SENSITIVE / STAGE C PENDING |
| Oracle Stage C low-T/high-T/stiff minimal grid | 当前 HEAD `095b1a7`、每点 8 个 native representative frequencies、`z_tail=5/6/7/8/10`、workers=1、Numba=2、BLAS=1；low-T/high-T/stiff observed systematic rel `4.9874e-4/4.1532e-3/1.7799e-3`，分别 `non-monotone/non-monotone/monotone`；结果保存于 `docs/oracle_tail_convergence_stageC_min8.json`，无 tail correction promotion | PASS / ORACLE-SENSITIVE / STAGE D PENDING |
| Nested-parallel resource audit | 为 `param_sweep`、`validate_plain_grid_vs_reference`、`validate_edges_vs_reference`、`benchmark_same_grid_reference` 增加 `--threads`；外层 workers/pool 大于 1 时内层自动降为 1，并新增 `nested_thread_budget` 单元测试；默认仍为 workers/pool=1、Numba=2 | PASS / RESOURCE-CAPPED |
| Per-frequency tail/adiabaticity diagnostic | HEAD `f939c07`；default/lowT/highT/stiff 各 8 个 native frequencies、`z_tail=5/6/7/8/10`，共 160 条模式记录；相对 z_tail=10 的单模式最大差异 `0.590%–0.672%`，z_tail=5 最大 `|omega'/omega^2|≈1.35e-2`，z_tail=10 降至约 `9.1e-5`；phase、amplitude、omega、二阶 adiabaticity 与 DN weight 已记录，resume 命中 160/160 | PASS / ORACLE-DIAGNOSTIC / NO PRODUCTION CHANGE |
| Oracle B phase-averaged assessment | 证明当前 reference tail 已用 `sqrt((x_f^2+y_f^2)/2)` 做 phase-averaged envelope；同一 handoff 上的 phase-averaged observable 与 Oracle A 代数等价，不是独立 oracle，不能缩小 systematic；报告 `docs/oracle_b_phase_averaged_assessment.md` | PASS / REJECTED FOR PROMOTION |
| Oracle B finite phase-window prototype | default/low-T/high-T/stiff 各 4 个代表 native frequencies；z=5 后继续 DOP853 到 z=5.5/6/7；最大 z=5 对 z=7 相对变化 `5.321e-3/5.552e-3/3.927e-3/6.897e-3`；低频未入尾部显式保留；未改变正式路径 | PASS / DIAGNOSTIC / NOT PROMOTED |

### Errors

| Error | Resolution |
|-------|------------|
| 计划文件首次补丁未匹配 init-session 生成的简化模板 | 重新读取实际文件后以完整替换方式更新 |
| 直接查询 `numba.threading_layer()` 时线程层尚未初始化并抛出 `ValueError` | 将线程层初始化纳入基准 warmup，再记录实际 layer；不把该探测异常误判为 solver failure |
| PowerShell 内联 Python 装饰器/多行函数探测命令发生 `SyntaxError` | 不再用复杂 `python -c`；改用已有脚本或通过 `apply_patch` 创建可复现诊断脚本 |
| 固定 profiler 首次用 `from scripts import profile_fast_breakdown` 导入失败 | `scripts` 无包初始化文件；改为把脚本目录加入 `sys.path` 后直接导入 |
| 全仓 `ruff check stiffgwpy_fast scripts tests` 暴露 200+ 个既有 lint 错误 | 未对无关旧文件做大范围格式化；改为对本轮新增/修改 profiling scripts 做 scoped ruff，结果通过，并保留全仓门禁未通过事实 |
| quadrature 比较脚本先导入 SciPy/NumPy 后才设置 `NUMBA_THREADING_LAYER`，实际 layer 为 `tbb` | 该次输出作废；已将环境设置移到数值库导入前，重新运行时强制验证 `workqueue` |

## Session: 2026-09-16 (current-HEAD Prüfer kernel feasibility)

### Actions Taken

- fetch `fast/fast_v0.2` 并复核当前 HEAD、`progress.md`、`findings.md`、fresh
  runtime matrix、outer-reuse A/B、最近 rejected experiments；确认无新的远端代码。
- 以 2 threads、workqueue、BLAS=1、reference workers=1 重新运行完整 native-grid
  Prüfer/Cartesian DOP853 对照，覆盖 default/high-T/stiff/high-kappa。
- 创建并删除临时 Numba phase-kernel spike；先测试 midpoint，发现 crossing handoff
  比较口径问题后改为 production 同一 `kend`，再测试 RK4，并补充同初值 DOP853
  与同一 Cartesian exact-step 对照。

### Test Results

| Test | Result | Status |
|---|---|---|
| current-HEAD full-grid Prüfer vs Cartesian, z_tail=5 | 4 points, 76/77 native nodes；DN relative `3.02e-10/2.53e-10/2.66e-10/4.66e-10`；0 numerical failure；outer difference约 `2–5e-10` | PASS / ORACLE VERIFIED |
| Numba phase RK4 propagation spike | 25-repeat warm median `2.41/2.31/2.41/2.33 ms`；handoff amplitude error `1.74%/0.97%/1.48%/2.37%`；phase `0.0136/0.0181/0.0082/0.0103 rad` | REJECTED |
| clean canonical pytest | `159 passed, 6 deselected, 2 warnings` | PASS |
| recent remote CI | run `34628009527`, HEAD `a1d701a5`, 9/9 success | PASS |

### Decision

- 固定步 midpoint/RK4 Prüfer kernel 不进入 production；其百分比 amplitude 偏差
  与无稳定 speedup 不满足 acceptance criteria。当前正式 fast 行为未改变，未产生
  可推送的代码突破。
- 当前最高价值后续实验是 exact Cartesian transfer map 的 phase-coordinate
  重写或 nonoscillatory Riccati phase-function，而不是继续调 fixed-step RK 阶数。

### Errors

| Error | Resolution |
|---|---|
| package gate 首次使用默认镜像时出现非 UTF-8 解码错误，构建工具未能完成 | 切换到官方 PyPI index 后重新执行 sdist、wheel、distribution verification 和 installed-wheel smoke，全部通过 |

## Session: 2026-09-16 (current-HEAD native tail A/B)

### Actions Taken

- fetch 并确认远端 `fast_v0.2` 与本地 HEAD 均为 `16ea185c84f50f081a6a1c17e9c47e2380d9d67a`。
- 重新运行 current-SHA full native-grid Prüfer/Cartesian tail nested A/B，覆盖
  default、low-T、high-T、stiff、high-kappa 的 z=5/6/7；所有 artifact 记录
  Numba=2、workqueue、BLAS=1、workers=1 与代码 SHA。
- 重新运行 25-repeat production warm matrix；同时以 per-process 显式 `set_z_tail`
  运行 production handoff A/B，避免 benchmark preset 覆盖环境变量造成假实验。

### Test Results

| Test | Result | Status |
|---|---|---|
| Prüfer/Cartesian native tail A/B | production-like oracle 的 z=5 到 z=7 DN 漂移约 low-T `1.4e-3`、其他点 `2.8–2.9e-3`；pairwise cross-oracle 差约 `1e-9` | PASS / TAIL SYSTEMATIC CONFIRMED |
| Production z_tail 5/6/7 | DN handoff 变化仅 `4.1e-6/1.0e-4/5.1e-6/5.2e-6/4.5e-6`；z=7 warm median 约 `5.74/5.20/8.22/8.59/8.68 ms` | REJECTED FOR PROMOTION |
| Fresh production matrix, z=5 | 25 repeats；default median/p95 `5.09/5.49 ms`；high-T/stiff/high-kappa median `6.90/6.57/6.09 ms`；0 failure | PASS / PERFORMANCE BASELINE |

### Decision

- handoff 深度不是当前 production DN 主要瓶颈；不改变 production z_tail。
- 下一阶段转向 propagation/grid/model alignment，避免重复 fixed-step phase 与简单
  z_tail 加深方向。

### Fresh breakdown profile

- 当前 SHA 四点 7-repeat profile 的 tensor-solve kernel median 为
  default/high-T/stiff/high-kappa `4.57/4.36/4.75/4.52 ms`；kernel prepare
  `1.37/1.12/1.39/1.08 ms`；PCHIP `0.13/0.12/0.13/0.13 ms`。
- 结论：speed hotspot 已定位到真实 tensor solve/传播路径；frequency quadrature
  不是可接受的主要优化对象。
- profiling artifact 现包含 `commit` 字段；四份当前 SHA profile 已重生成并待随本轮
  provenance tooling 记录一起提交。
- 后续复核确认正式 profile 必须显式 `kink_split=true`；已在 release HEAD 重生成四份
  正式路径 artifact，替代此前 false-kink profile 的速度结论。
## Session: 2026-09-16 (current-HEAD step-size Pareto)

### Actions Taken

- 刷新当前 SHA 的 Oracle C WKB anchors，并以 h=.005 baseline、h=.00625、h=.0075、
  h=.01 做正式 fast runtime A/B。
- 对 h=.0075 与 h=.005 运行相同 10 点 matched continuous-sigma reference 矩阵，覆盖
  named edge/tilt/regime 与 Sobol；另检查当前 quadrature estimator 的 actual/predicted、
  p95/p99、false-safe。

### Test Results

| Test | Result | Status |
|---|---|---|
| step-size runtime A/B | h=.0075 相对 h=.005 快约 `15–34%`；h=.00625 收益不稳定 | PASS / CANDIDATE |
| Oracle C WKB/deep | h=.0075 多数点误差增加；h=.01 更差 | REJECTED FOR PRECISION REGRESSION |
| matched parameter matrix | h=.0075 signal max rel `6.7833e-3`、DN max `3.6052e-3`；baseline 分别 `6.7832e-3`、`3.6014e-3` | REJECTED / MODEL ALIGNMENT STILL OPEN |
| estimator diagnostic | 9 points；coverage `1.0`；false-safe `0`；p95/p99 actual/pred recorded | PASS / DIAGNOSTIC ONLY |

### Decision

- 不改变 fast 默认 h；本轮没有 production algorithm change。
- 最高价值方向仍是 exact Cartesian transfer-map phase representation 或 Riccati
  reduction，且必须直接验证 continuous-sigma model alignment。

## Session: 2026-09-16 (linear-z Magnus transfer-map prototype)

- 预注册并实现 standalone 二阶 Magnus map：每个 native step 对线性插值的 `z(N)`
  构造无迹矩阵指数；当前 production midpoint map 保留作 baseline。
- current SHA `d6aef07` artifact `docs/transfer_map_phase_round.json` 覆盖
  5 个命名 regime、35 个 tail modes、z_tail=5、Cartesian DOP853 rtol=1e-9，
  以及每个 regime 25-repeat runtime A/B；资源固定 Numba=2/workqueue/BLAS=1/
  workers=1。
- focused 对照：最大振幅变化 `6.80e-6`、最大相位变化 `4.23e-4 rad`，没有
  10x 精度改善；baseline -> candidate median 分别为
  `2.56->45.52`、`2.58->49.61`、`3.02->53.14`、`2.65->47.31`、
  `2.75->45.93 ms`（default/low-T/high-T/stiff/high-kappa）。
- 决策：`REJECTED FOR PRODUCTION`；未做 full-grid promotion，生产 fast 默认与
  输出契约不变。下一实验必须转向能减少真实 tensor solve 的解析/Numba transfer
  map 或 Riccati reduction，并同时测 model alignment。

## Session: 2026-09-16 (closed-form Numba Magnus follow-up)

- 将上一轮二阶线性-z Magnus map 从 `scipy.linalg.expm` 换为无迹 2x2 矩阵指数
  解析闭式并用 Numba 编译；artifact `docs/transfer_map_phase_numba_round.json`
  绑定 current SHA `19869ef`。
- 五个命名 regime、35 个 tail modes：最大振幅变化 `1.486e-5`、最大相位变化
  `7.109e-4 rad`，35/35 模式仅有微小方向改善，未达 10x 精度标准。
- 25-repeat propagation median baseline -> candidate 为
  `3.254->5.076`、`2.667->4.950`、`2.602->4.821`、`2.570->4.708`、
  `2.600->4.751 ms`（default/low-T/high-T/stiff/high-kappa），全部超过 5%
  成本预算。
- 决策：线性-z 二阶 Magnus 假设整体 `REJECTED FOR PRODUCTION`；无 full-grid
  promotion，production fast 默认不变。下一项转向 nonoscillatory Riccati 或
  能真实减少 tensor solves 的解析分支。

## Session: 2026-09-16 (fresh four-regime breakdown)

- 在 HEAD `65578a5` 重新以正式 `kink_split=true`、25 repeats、Numba=2、workqueue、
  BLAS=1、workers=1 profile default/high-T/stiff/high-kappa。
- total median 为 `8.384/12.379/12.577/12.449 ms`，tensor solve median 为
  `3.174/5.329/5.575/5.630 ms`，mean steps/channel 为 `8187.4/8143.9/8229.6/
  8043.2`（按 default/high-T/stiff/high-kappa）。
- tensor solve/propagation 是最大真实热点；prepare、background、PCHIP 均非首要
  瓶颈。下一项转向减少真实 channel propagation work 的 Riccati/解析分支。

## Session: 2026-09-16 (assembly lower-bound profile)

- 以 `PROFILE_ASSEMBLE=0` 做理论下界 profile，覆盖 default/high-T/stiff/high-kappa，
  25 repeats、Numba=2/workqueue/BLAS=1；仅作速度诊断，不作数值输出。
- total median 正式路径 -> no-assembly：`8.384->8.090`、`12.379->12.066`、
  `12.577->12.291`、`12.449->11.988 ms`；tensor median：
  `3.174->2.497`、`5.329->4.361`、`5.575->5.117`、`5.630->4.794 ms`。
- total 收益仅 `2–4%`，排除 assembly 作为 10% 级优化方向；下一阶段仍需减少真实
  propagation arithmetic/steps，转向有极点规避的 Riccati phase-function 原型。

## Session: 2026-09-17 (raw Riccati pole audit)

### Actions Taken

- fetch `fast_v0.2` 并确认本地/远端 HEAD 均为
  `d5870311bfcdc814220181d911f3d938a4ebfaf0`。
- 阅读当前 progress/findings、验证清单、既有 profile 与最近 phase/Magnus rejected
  artifacts；未使用旧 profile 推断本轮热点。
- 在不改变 production path 的 standalone 脚本中，从 `reference._tensor_orig` 推导
  并积分 `r=x/y` Riccati；极点作为 terminal fail-closed 状态记录。
- 以 8 个参数点（4 named regimes、2 edge、high-T_re edge、固定 Sobol）和 8 个
  频率完成 64-mode audit，并以独立 Cartesian DOP853 重放每个极点位置。
- 以当前 HEAD、正式 `kink_split=true`、2 threads/workqueue、BLAS=1、25 warm
  repeats fresh profile default/high-T/stiff/high-kappa。

### Test Results

| Test | Result | Status |
|---|---|---|
| raw Riccati pole audit | 64/64 modes pole before z=5；0 numerical failure；0 tail reached | REJECTED |
| Cartesian cross-check at pole | 64/64 solve success；max `abs(y)/abs(x)=1.053e-6` | PASS / POLE CONFIRMED |
| fresh 2-thread four-regime profile | total median `6.442/9.427/9.750/9.719 ms`；tensor median `2.633/4.513/4.931/4.578 ms` | PASS / HOTSPOT CONFIRMED |

### Decision

- raw `x/y` Riccati `REJECTED FOR PRODUCTION`；不进入 runtime A/B、API 或 fallback。
- 当前 propagation/tensor solve 仍是速度首攻；tail、PCHIP 和 assembly 暂无理由优先。
- 下一项 phase 研究必须改为无极点的 complex-log-derivative/nonoscillatory
carrier，并先建立数学残差与 Cartesian/Prüfer 独立交叉验证。

## Session: 2026-09-17 (phase/WKB and closed-form branch consolidation)

- 以当前 HEAD `85cbfc50653e8d11dd3a1398a4496fc9d0a345d0` 重新 fetch，确认远端
  `fast_v0.2` 同 SHA；整理本轮 standalone scripts 与 JSON artifacts，生产
  `fast_sgwb` 路径未改变。
- WKB carrier 在 `z_match=3.625/3.75/4.0` 的四 regime 对照中，`z=3.625`
  的 DN proxy 为 `1.50e-4..1.84e-4`，但仅减少后段传播，按 fresh total profile
  预计总收益不足 10%；其余 match 点同样未形成跨 regime 总体突破。
- nested native-frequency 复核确认已有局部 estimator 在多 regime 上
  `false_safe=true`，不能作为独立 error bound；同插值器内的 PCHIP/Gauss 对照
  只保留为诊断。
- Oracle C WKB 当前 HEAD fresh true-error：default/high-T/stiff/high-kappa
  对深锚点分别约 `5.21e-6/2.24e-5/3.35e-6/2.45e-5`，low-T 约 `1.72e-4`；
  precision 首要剩余风险仍是 low-T/transition 与 tail oracle uncertainty。
- 新假设“辐射主导段常系数闭式 transfer map”在 `z_match=3.75` 四 regime
  的最大振幅误差为 `6.19e-3..1.86e-2`，传播 median 反而增加约 `11.0–17.5%`；
  `REJECTED FOR PRODUCTION`。

## Session: 2026-09-17 (phase-substep exponential recurrence prototype)

### Actions Taken

- 重新 fetch `fast/fast_v0.2`，确认本地与远端均为 `58af509ba5e15c28017ded9c4abceebc3594c5c2`；重新读取 progress/findings、validation manifest、round7 fresh profiles 与 rejected experiments。
- 新增 standalone `scripts/benchmark_phase_recurrence.py`。候选只替换 phase
  substep 内 `exp(z)` 的计算方式，保留完整 assembly、kink split、tail matching
  和 production `phase_max`；正式 `fast_sgwb.py` 未修改。
- 完成 2-thread/16-thread 各 30-repeat kernel A/B、五工况 full outer solve、
  named+edge+Sobol 14 点矩阵，以及当前 HEAD 的 fresh Oracle C focused 对照。

### Test Results

| Test | Result |
|---|---|
| 2-thread kernel median ratio default/high-T/stiff/high-kappa/low-T | `0.873/0.809/0.899/0.898/0.904` |
| 16-thread kernel median ratio default/high-T/stiff/high-kappa/low-T | `0.802/0.851/0.926/0.932/0.916` |
| Full outer solve | 五工况 candidate/production 均 converged |
| Named + edge + Sobol matrix | `14/14` status match，`false_safe_count=0`，最大 spectrum delta `6.99e-6 dex`，最大 DN delta `6.44e-6` |
| Fresh Oracle C | candidate 未新增 oracle systematic；low-T 独立保留，candidate-vs-WKB `1.659e-4` |

### Decision

`ACCEPTED AS STANDALONE PROTOTYPE / NOT PRODUCTION`。recurrence 在目标 regime
达到速度门，但 digest 不逐位一致，且尚未完成 full-grid Oracle A/Prüfer/WKB 与
更大参数空间认证；不修改 production，不更新当前 stale validation manifest。

### Artifacts

- `scripts/benchmark_phase_recurrence.py`
- `docs/phase_recurrence_round_20260917.json`
- `docs/phase_recurrence_round_16t_20260917.json`
- `docs/phase_recurrence_fullsolve_20260917.json`
- `docs/phase_recurrence_parameter_matrix_20260917.json`
- `docs/phase_recurrence_oracle_20260917/`
- 完整整理的 artifacts 包括 WKB、block-transfer、radiation closed-form、
  true-error、outer-reuse、residual decomposition、validation 与 estimator
  复核文件；资源记录为 Numba/workqueue=2、BLAS=1、reference workers=1（涉及
  25-repeat standalone 对照），无 silent fallback 或 API 变化。

## Session: 2026-09-17 (current HEAD fresh outer-reuse headroom)

### Actions Taken

- 从 `fast` 远端 fetch 后确认本地/远端 HEAD 均为
  `f411bbb252e20b7626b33d007ea73d40bf34f46c`；重读现有 progress/findings、manifest、
  benchmark 与 rejected records。
- 在正式 `kink_split=true`、Numba=2/workqueue、BLAS=1、workers=1 下完成六点 25-repeat
  fresh runtime；另完成 default/high-T/stiff/high-kappa 的 25-repeat 阶段 breakdown。
- 对现有 outer reuse 判据与 forced reuse 完成 7 个 accepted 点的 25-repeat A/B；显式
  `shared_Neff_guard` 点单独记录为 guard，不计 numerical failure。

### Test Results

| Test | Result | Status |
|---|---|---|
| fresh six-point runtime | median/p95 default `6.99/8.03 ms`；high-T `11.76/12.23`；stiff `12.38/13.44`；high-kappa `12.20/12.92`；六点 converged、0 failure | PASS / BASELINE |
| fresh stage breakdown | high-T/stiff/high-kappa tensor median `7.07/7.62/7.20 ms`；准备层各项均低于约 `1.6 ms` | PASS / HOTSPOT CONFIRMED |
| forced outer reuse | high-T/stiff/high-kappa runtime ratio `0.735/0.747/0.757`，spectrum max `3.944e-3/1.051e-3/1.547e-2 dex`；DN relative `3.387e-11/3.053e-10/3.283e-11` | REJECTED |
| low-T / independent points | low-T ratio `0.976`、spectrum `0`；positive-tilt ratio `0.705`、spectrum `2.075e-1 dex`；sobol_000 ratio `0.988`、spectrum `0` | REJECTED FOR SIMPLE THRESHOLD WIDENING |
| manifest refresh | `build_two_mode_manifest.py` 重建并绑定当前 HEAD；fast plain-grid `NOT VERIFIED`、transition-refine `PARTIALLY VERIFIED` | PASS / PROVENANCE |

### Decision

- 不放宽固定 outer reuse threshold，也不修改 production solver。DN_gw 单值几乎不变
  不能作为安全判据；频谱误差在 high-T/stiff/high-kappa 与 positive-tilt 显示明显
  observable false-safe。
- 下一候选只能是 standalone observable-aware criterion，组合 `delta log Omega`、
  `delta Phi`、`delta S2`、horizon shift 与 frequency-weighted DN sensitivity；本轮
  尚未具备接受该 criterion 的证据。

### Errors

| Error | Resolution |
|---|---|
| `run_fixed_profile.py` 外层只允许 `A/B` | 改用支持六点的 runtime benchmark，并以独立进程直接调用 profiler 获取四个目标 regime；未据失败输出形成结论 |
| outer-reuse A/B 包含 guard 点后访问缺失的 `DN_gw` | 排除显式 guard 点完成 accepted-point A/B；guard 仍保留为物理边界问题 |

## Session: 2026-09-18 (observable-aware reuse proxy prototype)

### Actions Taken

- fetch `fast_v0.2` 后确认本地/远端 HEAD 为 `083fdf5c44d18e2d417ed9561571110aaa9ee29a`，
  并重读当前 progress/findings、manifest、benchmark 与 rejected records。
- 新增 standalone `scripts/benchmark_outer_observable_proxy.py`，强制关闭 shortcut
  以捕获 first→second outer update，记录 `delta log Omega`、`delta Phi`、`delta S2`、
  horizon shift 与 frequency-weighted DN sensitivity proxy；不改变 production code。
- 修正两处原型数据流错误：`prep_frequency_only` 不含 Phi，改从真实 kernel 参数捕获；
  频率数组降序导致积分符号错误，改为按 log-frequency 升序积分。

### Test Results

| Test | Result | Status |
|---|---|---|
| standalone proxy | 当前 HEAD、Numba=2、BLAS=1、workers=1，5 个需 second solve 的点成功捕获 2 kernels；low-T 正常 1 kernel 收敛 | PASS / DIAGNOSTIC |
| proxy values | high-T `dlogOmega=3.944e-3`、`dPhi=7.694e-5`、`dS2=6.560e-5`、`dH=2.008e-3`；high-kappa `1.547e-2/2.964e-4/1.800e-4/7.728e-3`；stiff `1.051e-3/2.092e-5/1.300e-8/5.349e-4` | OBSERVED |
| focused validation | `py_compile` 与 `git diff --check` 通过 | PASS |

### Decision

- 当前代理只能事后复现 first→second 的频谱差；`delta Phi`、`delta S2`、horizon shift
  与频率加权 proxy 尚未形成可提前保证 spectrum budget 的上界。尤其 high-T/high-kappa
  的小代理值对应明显频谱偏差，暂不进入 runtime A/B 或 production。
- 原型 artifact 为 `docs/outer_observable_proxy_round_20260918.json`。下一步应先把
  kernel 输出对背景扰动的局部线性响应做成独立 sensitivity estimate，再谈
  coverage/p95/p99/false-safe；禁止把本轮诊断称为 estimator。

## Session: 2026-09-19 (local observable response sensitivity)

### Actions Taken

- fetch `fast_v0.2` 后确认本地/远端 HEAD 为 `dca9bd85bd9e0aae19bcbd871077bad06601dc0d`，
  并重读当前研究记录、manifest、benchmark 与 rejected artifacts。
- 扩展 standalone proxy：保存首轮 kernel 输入，分别对 `Phi/Phi_mid` 与 `S2/S2inv`
  注入真实 first→second 背景差分，运行额外 reference-only kernel，测量一阶响应预测
  的 `dlogOmega`；正式 solver 未改变。

### Test Results

| point | Phi predicted | S2 predicted | actual dlogOmega |
|---|---:|---:|---:|
| high-T | `7.443e-3` | `4.017e-3` | `3.944e-3` |
| stiff | `2.840e-3` | `1.063e-3` | `1.051e-3` |
| high-kappa | `2.795e-2` | `1.544e-2` | `1.547e-2` |

`default`/`positive-tilt` 的 Phi/S2 背景差为零，但实际频谱差分别为
`3.923e-4`/`6.226e-4`，说明仍有未捕获的 kernel/尾部响应项；low-T 首轮已收敛，
未形成 second-solve 样本。

### Decision

- 结果支持继续推导“背景输入到 observable 的局部响应”，尤其 S2 响应在 stiff 与
  high-kappa 上接近实际偏差；但当前预测依赖额外 kernel 调用，不是可部署 criterion。
- 状态：`CANDIDATE FOR FURTHER DERIVATION`，不进入 runtime A/B、production 或
  fallback。下一步要把 response 化为不调用第二次 kernel 的解析/低成本估计，并在
  named+Sobol 上验证 coverage、p95/p99、false-safe=0。

## Session: 2026-09-20 (analytic observable-response boundary)

- 当前 HEAD `0253550ee4ac5e920125444ba138775ffdfbefed` 上，standalone proxy 已覆盖
  全部 24 个 named/edge/Sobol 输入；显式 physical guard 标记为 ineligible。
- 新增 S2 端点缩放与固定 j0 horizon-start 的解析诊断项；`14/24` eligible，安全因子
  1.1 下 coverage `4/14`、false-safe `0/4`，eligible actual spectrum-delta
  p95/max 为 `8.217e-2/2.075e-1 dex`，high-T/stiff/high-kappa 为 `0/3` safe。
- Decision: `REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`；下一轮补齐
  `j0/z0/tail/phase-path` 响应并用 Cartesian/Prüfer/WKB oracle 复核，不做 runtime A/B。

原始 artifact：`docs/outer_observable_proxy_full_round_20260920.json`。

## Session: 2026-09-17 (adiabaticity-triggered carrier prototype)

- 在最新 HEAD `a65eaa73d3e0ca0adec31b6da26e2f154578f45f` 上，将 standalone WKB-carrier
  spike 扩展为按 mode 的触发器：连续 3 个 native node 满足
  `|1.5*sigma-1|*exp(-z) <= eps_trigger` 才切换；正式 solver 未改变。
- `eps=3e-4` 的 amplitude p50/p95/max 为 `2.429e-4/5.993e-4/5.993e-4`，但
  runtime ratio 为 `0.857..1.033`；`eps=1e-3` 的 amplitude max 为 `2.376e-3`，
  runtime ratio 为 `0.970..1.252`。`eps=3e-3/5e-3` 的 amplitude max 达
  `6.654e-2/7.375e-2`。
- 低阈值触发点晚于生产 `z_tail=5`，已加入显式限制；未触发 mode 不进入 carrier。

Decision: `REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`；若继续 hybrid，必须加入
二阶 adiabaticity、独立 full-grid oracle 与 25--50 repeat A/B。

## Session: 2026-09-17 (outer goal-grid invariant A/B)

- 当前远端/本地 HEAD 为 `f9b4e9828bad9eab75e09cd7212658211040978b`；fresh 25-repeat
  baseline（Numba=2、affinity=2、BLAS=1、workers=1）为 default/highT/stiff/
  high-kappa `7.907/12.982/14.263/13.584 ms` median，p95 为
  `8.605/13.932/15.351/16.852 ms`。16-thread 分层 profile 的目标点 total median
  为 `6.557/6.556/6.355 ms`，tensor 为 `1.701/1.693/1.533 ms`；准备层与 outer
  重复调用仍是主要可研究部分。
- 新增 standalone `benchmark_outer_static_invariants.py`，确认 high-T/stiff/
  high-kappa 第二轮 `Nv` 仅有 `1e-12` 量级漂移，goal grid 仅有 `1.6e-12` 量级漂移，
  但 sigma/f_hor/primitive 仍随背景更新，不能凭近似相等直接复用。
- 新增 standalone `benchmark_outer_grid_reuse.py`，在 13 个 named/Sobol/edge 点以
  25-repeat A/B 测试复用首轮 goal grid。`10/13` 点 final `f` digest 不一致；最大
  频率差 `1.56e-12`，而 `log10OmegaGW/DN_gw/g2/w2` 的数值差为 0。speed ratio
  范围 `0.960..1.003`，没有稳定 >5% 收益。

Decision: `REJECTED FOR PRODUCTION`。即使科学量数值相同，输出契约要求 digest 完全
一致；该候选不能作为准备层去重。下一候选应转向 background/primitive 的可证明缓存，
并继续保留 low-T 独立检查。

## Session: 2026-09-17 (local derived-param cache A/B)

- 当前 HEAD `a9564eab3f0203a3d57dac46127047142ee30059` 的 fresh profile 显示，
  `gen_expansion` 一次调用内重复读取 `derived_param`；旧实现测试计数为 9 次。
- standalone `scripts/benchmark_derived_param_local_cache.py` 以 25-repeat A/B 覆盖
  low-T、high-T、stiff、high-kappa；Nv、sigma、f_hor、spectrum、DN_gw digest
  全部逐位一致。
- candidate/baseline warm median ratio 为 `0.9963/0.9813/0.9912/0.9979`，最大稳定
  收益约 1.9%，未达到 >5% runtime acceptance gate。

Decision: `REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`。正式实现已撤回，原始
A/B 结果保留在 `docs/derived_param_local_cache_round_20260917.json`。

## Session: 2026-09-17 (H2 endpoint cache acceptance)

- 当前 HEAD `6dff78ef440f9d4f8e363a7ee7f3016ebd307137` fresh 16-thread profile：
  high-T/stiff/high-kappa total median `6.36/6.94/6.95 ms`，对应 background
  `0.76/0.68/0.74 ms`、goal preparation `0.40/0.40/0.47 ms`。
- 在 goal-grid 的 `H2_vec` 端点增加函数内缓存；测试先复现旧实现 4 次端点调用，
  候选降为不超过 3 次。50-repeat、13 个 named/Sobol/edge 点 A/B 的 digest mismatch
  为 `0/13`，spectrum 与 DN 的 candidate-baseline p50/p95/max 均为 `0`。
- 目标 regime candidate/baseline warm median ratio：high-T `0.9838`、stiff `1.0110`、
  high-kappa `0.9346`，即 high-kappa 稳定改善约 `6.54%`，满足 >5% 速度门槛。
- 独立 Cartesian/Prüfer oracle：三点 DN 相对差 `8.88e-10–1.996e-9`，频谱最大相对差
  `1.49e-7–2.73e-7`；未观察到新的 error-budget 项。13/13 输出契约比较通过，
  因而参数空间 false-safe 计数为 `0`。

Decision: `ACCEPTED FOR PRODUCTION`。该改动只缓存同一次 grid build 的纯 H2 端点值，
不跨 outer iteration；正式实现与回归测试保留，完整 A/B 见
`docs/h2_endpoint_cache_round_20260917.json`，oracle 见三个 `docs/h2_endpoint_oracle_*`。

## Session: 2026-09-17 (remove unrelated comment-language CI gate)

- 按用户要求删除 `.github/workflows/ci.yml` 中的中文注释门禁及其专用脚本，
  同时更新测试覆盖矩阵；保留 Ruff、mypy、compileall、manifest、pytest、package、
  Cobaya 和多 Python 版本兼容性测试。
- 编码声明不再作为 CI 失败来源；代码质量检查仍由 Ruff 和类型/编译门禁负责。

## Session: 2026-09-17 (analytic branch eligibility boundary)

- 当前 HEAD `9b979c36b1cce46e3beb268892a09df7e4a7d9a5` 完成 fresh 16-thread profile：
  high-T/stiff/high-kappa total median `6.16/6.74/6.68 ms`，tensor kernel
  `1.62/1.94/1.84 ms`。
- 新增 standalone `benchmark_analytic_branch_eligibility.py`，按实际 mode crossing
  与 `sigma` 统计 z_tail=5 前可进入常系数 radiation/stiff branch 的 segment；目标点
  eligible fraction 为 `16.47%/14.31%/17.66%`，low-T 为 `6.62%`，严格 radiation
  segment 为 `0`。

Decision: `REJECTED FOR PRODUCTION / PROTOTYPE DEFERRED`。当前 eligible work 明显低于
>30% tensor-work 优先门槛，且 low-T 独立边界更差；不重复旧 fixed-z radiation 原型。
原始 artifact：`docs/analytic_branch_eligibility_round_20260917.json`。

## Session: 2026-09-17 (unused Psi preparation buffer)

- 当前 HEAD `f3e5f6b131daa639a036abb9fe04ca1e863171ca` fresh 16-thread profile 的稳定
  preparation/kernel median 为 high-T `1.19/1.66 ms`、stiff `1.27/2.11 ms`、
  high-kappa `1.14/1.74 ms`；继续检查 `prep_kernel` 的 allocation/copy 热点。
- standalone 原型 `scripts/benchmark_prep_kernel_no_psi.py` 去掉未被正式 `solve_kernel`
  消费的 `Psi` 写入和 buffer，仅比较其余七个输出数组。default/high-T/stiff/
  high-kappa/low-T 各 50 次，所有比较数组逐位一致；准备层 candidate/baseline ratio
  分别为 `0.9767/1.0230/1.0239/0.9556/1.0237`。

Decision: `REJECTED FOR PRODUCTION`；最好 high-kappa 也只有约 `4.44%`，未达到稳定
`>5%` runtime gate，且其余目标区间变慢。证据保留于
`docs/prep_kernel_no_psi_round_20260917.json` 与四个 regime artifact；不进入正式代码。

## Session: 2026-09-17 (constant tridiagonal buffers)

- fresh round-4 profile 后检查 `prep_kernel` 的三对角 spline solve；`aa`/`cc` 全程为常数
  1，建立候选以去除两 个只读常数 buffer，同时保留其余计算路径。
- focused tests（进程启动即固定 `NUMBA_NUM_THREADS=2`）`54 passed, 3 deselected`；
  high-T/stiff/high-kappa 25-repeat candidate total median 为 `5.995/6.606/6.043 ms`，
  与 fresh baseline `6.241/6.787/5.818 ms` 的 ratio 为 `0.9605/0.9733/1.0385`。
  三个 regime 的 `f/log10OmegaGW/DN_gw/g2/w2` digest 均一致。

Decision: `REJECTED FOR PRODUCTION`；high-T 约 `3.95%` 改善、stiff 约 `2.67%` 改善，
high-kappa 反而变慢约 `3.85%`，未达到稳定 `>5%`。候选 patch 已回退，candidate profile
artifact 保留作 rejected evidence。

本轮另发现一次本地测试命令环境错误：测试进程未在启动前固定 Numba 线程，随后测试 fixture
改变环境导致 Numba 已初始化线程与环境值冲突。以后测试与 formal benchmark 使用独立进程，
并在启动前分别固定 `NUMBA=2/FAST=2` 或 `NUMBA=16/FAST=16`；远端 CI 当前 run 已验证通过。

## CI resource-environment isolation (2026-09-17)

- [x] 复核 `ci.yml` 与 `slow.yml`，中文注释/编码语言门禁及其脚本已删除，未发现同类
  无关检查；兼容性、数值回归、manifest、构建、Cobaya 和静态检查均保留。
- [x] 将资源环境契约改为干净子进程、数值库导入前验证，避免在已初始化 Numba 的
  pytest 进程中修改 `NUMBA_*` 环境变量。
- [x] 本地 CI 同口径 gate 全部通过：Ruff、mypy、compileall、diff check、manifest，
  `160 passed, 6 deselected`。

## Fresh round-5 profile and phase/primitive candidates (2026-09-17)

- [x] 在 HEAD `759d5cc`、16 threads/workqueue、BLAS=1、`kink_split=true`、25 repeats
  重新 profile：default/high-T/stiff/high-kappa/low-T total median 为
  `5.35/6.86/7.49/7.20/4.32 ms`，对应 p95 为 `5.74/7.89/8.26/8.18/5.67 ms`。
  tensor kernel median 为 `1.08/1.81/1.92/1.85/1.15 ms`，exact phase primitive 为
  `0.90/1.35/1.49/1.46/0.89 ms`。原始 artifacts 为
  `docs/profile_fast_breakdown_round5_*_20260917.json`。
- [x] exact primitive allocation/fusion standalone prototype：五个 regime 各 50 次，
  digest 完全一致、max abs/rel 均为 0，但 fused/baseline median ratio 为
  `1.0355/0.9874/0.9947/1.0024/1.0092`；REJECTED，未进入 production 或 oracle gate。
- [x] phase envelope exponential reuse candidate：focused tests `54 passed, 3 deselected`，
  spectrum/DN digest 与 baseline 一致；kernel median ratio（default/high-T/stiff/
  high-kappa/low-T）约 `1.00/1.01/1.01/0.99/0.91`，没有目标区间稳定 >5% 改善，
  REJECTED 并回退。候选 artifacts 为 `docs/profile_phase_exp_reuse_candidate_*.json`。
- [x] `PROFILE_ASSEMBLE=0` 诊断显示 assembly 可归属约 `0.19/0.17/0.03/0.16 ms`
  （high-T/stiff/high-kappa/low-T），不是足以单独达到目标的主热点；不据此做近似。

## Fresh round-6 profile and kernel specialization candidates (2026-09-17)

- [x] 在 HEAD `42b5526`、16 threads/workqueue、BLAS=1、`kink_split=true`、25 repeats
  fresh profile：default/high-T/stiff/high-kappa/low-T total median 为
  `4.60/6.36/7.27/6.44/4.45 ms`；tensor median 为 `1.03/1.78/2.04/1.75/1.11 ms`。
  原始 artifacts 为 `docs/profile_fast_breakdown_round6_*_20260917.json`。
- [x] static `gen_fast` grid reuse standalone：目标 regime 第二次 `Nv` 仅最后 anchor 改变；
  50-repeat total ratio high-T/stiff/high-kappa 为 `0.96/0.99/0.98`，digest/DN 完全一致，
  REJECTED。
- [x] no-assembly first-probe kernel standalone：25-repeat ratio 为
  default/high-T/stiff/high-kappa/low-T `0.97/0.96/1.00/0.95/1.01`，digest/DN 完全一致，
  REJECTED。
- [x] 组合两项后 high-T/stiff/high-kappa/low-T ratio 为 `0.956/1.015/0.961/1.000`，
  stiff 退化，整体 REJECTED；组合 artifact 为 `docs/probe_grid_combo_20260917.json`。

## Current session update: 2026-09-17 (full native-grid recurrence audit)

- [x] fetch 后确认 local/remote HEAD 均为 `a1b0393358c1a6fccf5704217670b8d5e00ed659`，
  production fast 路径未改；重新读取 required state、manifest、benchmark 和 rejected
  records，旧 manifest 未作为当前认证证据。
- [x] 完成 default/low-T/high-T/stiff/high-kappa 的完整 native goal grid
  production-versus-recurrence-versus-Prüfer/DOP853 对照，全部 `status=ok`；最高
  recurrence spectrum p95/max 为 stiff `4.36e-6/7.60e-6`，最高 full-grid DN
  relative 为 low-T `3.12e-7`。
- [ ] 下一步继续 recurrence 的 reheating/kink edge、扩大 Sobol、16/20-thread total
  runtime 与 Oracle A/WKB 认证；未完成前不得接入正式 fast。

## Current session update: 2026-09-17 (recurrence edge/Sobol and total-runtime gate)

- [x] fetch 后确认 remote/local HEAD，重新读取 progress/findings、manifest、benchmark
  artifacts 与 rejected records；当前 validation manifest 仍绑定旧 commit，未冒充本轮证据。
- [x] recurrence full-grid 扩展到全部 24 个 named/edge/Sobol 点：22 accepted、2
  同步 `shared_Neff_guard`，无 numerical failure。
- [x] full outer matrix：`status_mismatch=0`、`false_safe_count=0`，最大 spectrum
  delta `6.988e-6 dex`，最大 DN relative `6.441e-6`；low-T 独立误差仍受控。
- [x] 16/20-thread 各 30-repeat total runtime A/B 完成；high-T/stiff/high-kappa
  未形成稳定超过 5% 的 total-runtime 改善，recurrence `REJECTED FOR PRODUCTION`，
  保留 standalone 证据。
- [ ] 下一候选回到 fresh profile 的 background/primitive/grid allocation 或 tensor
  kernel 结构；不得重复 recurrence、固定 threshold widening 或已列 rejected 路径。

## Current session update: 2026-09-17 (inline tensor transfer-map spike)

- [x] 重新 fetch 并复核当前 HEAD `10f539b`；发现 proposed 的
  `fast_phi_s2_split` allocation/fusion 已在历史记录中完成且已 rejected，未重复。
- [x] standalone prototype 只将 `scaled_step` 的 2x2 transfer-map 表达式内联到
  phase loop；phase subdivision、kink split、tail、assembly 和 production source
  均未改变。脚本提交于 HEAD `0671e1f`。
- [x] 当前 HEAD 下完成 2-thread kernel/outer 25--30 repeat，以及 16/20-thread
  formal outer A/B；所有 artifact 均带 `commit=0671e1fc549bf45428496935c364956b17b64405`。
- [x] 2-thread outer ratio high-T/stiff/high-kappa 为 `0.9351/0.9488/0.9137`，
  formal 16-thread 为 `0.9856/0.9974/0.9127`，20-thread 为 `1.0690/0.9876/0.9986`；
  不能证明稳定 >5% total-runtime 改善。
- [x] 2-thread outer spectrum delta p50/p95/max（dex）目标点为 high-T
  `8.28e-12/1.38e-6/2.66e-6`、stiff `1.89e-11/1.86e-6/3.30e-6`、
  high-kappa `5.92e-13/8.68e-7/2.42e-6`；DN relative 分别
  `3.01e-13/1.18e-11/7.70e-14`。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。不进入 Oracle、
Sobol 或 production candidate 认证；下一轮继续寻找未重复的 tensor/background
结构候选。

Artifacts：`docs/kernel_inline_spike_round8_20260917.json`、
`docs/kernel_inline_outer_spike_round8_20260917.json`、
`docs/kernel_inline_outer_spike_16t_round8_20260917.json`、
`docs/kernel_inline_outer_spike_20t_round8_20260917.json`。

## Current session update: 2026-09-17 (phase-exp hoist spike)

- [x] 当前 HEAD `0dae703` fresh round9 profile 已重新生成并固化；2-thread total
  median default/high-T/stiff/high-kappa/low-T 为 `6.284/9.257/9.660/9.501/6.587 ms`，
  tensor median 为 `2.505/4.409/4.895/4.663/3.077 ms`。
- [x] standalone phase-exp hoist 只复用 `n_sub=1` 分支中重复计算的同一
  `exp(z_mid)`；production source 未修改。脚本提交于 `06b7df1`。
- [x] 2-thread kernel 30-repeat ratio high-T/stiff/high-kappa 为
  `0.9238/0.9285/0.9911`；full outer 25-repeat ratio 为
  `0.9597/0.9610/0.9865`，未形成稳定 >5% total-runtime 改善。
- [x] full outer spectrum relative p95/max 为 high-T `3.18e-6/6.13e-6`、
  stiff `4.27e-6/7.60e-6`、high-kappa `2.00e-6/5.56e-6`；DN relative
  分别 `3.01e-13/1.18e-11/7.70e-14`。digest 不逐位相同但变化为舍入级。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。速度门在
2-thread full outer 已失败，不进入 formal 16/20-thread、Oracle 或 Sobol 认证。

Artifacts：`docs/phase_exp_hoist_spike_round9_20260917.json`、
`docs/phase_exp_hoist_outer_round9_20260917.json`。

## Current session update: 2026-09-17 (round9 assembly attribution)

- [x] 在当前 HEAD `df6880a`、Numba=2、BLAS=1、workers=1 下重跑目标三点
  `PROFILE_ASSEMBLE=0` 诊断；tensor median 从正式的
  high-T/stiff/high-kappa `4.409/4.895/4.663 ms` 降到
  `3.971/4.284/4.312 ms`。
- [x] 这是删除绝大多数中间 assembly 写入的理论上限诊断，不是可接受的 solver
  候选；完整 total median 没有形成稳定同幅收益，且会破坏中间列输出契约。
- [x] 结论：assembly-only 不能作为下一步 >5% total 优化方向；不修改 production，
  下一轮转向真正减少 propagation 工作量的 hybrid/phase 结构研究。

Artifact：`docs/profile_fast_breakdown_round9_noassemble_highT_20260917.json`、
`docs/profile_fast_breakdown_round9_noassemble_stiff_20260917.json`、
`docs/profile_fast_breakdown_round9_noassemble_high_kappa_20260917.json`。

## Current session update: 2026-09-17 (phase fastmath standalone)

- [x] 在 fresh HEAD `4f8acc3` 上只对 standalone phase-exp hoist twin 添加 Numba
  `fastmath=True`；不改变步长、handoff、公式、tail 或 production source。脚本随后
  以 `164bd84` 单独提交并推送，artifact 均绑定该 SHA。
- [x] 2-thread kernel 30-repeat ratio 为 default/high-T/stiff/high-kappa/low-T
  `0.8692/0.8468/0.9018/0.8695/0.9434`；full-outer ratio 为
  `0.9527/0.9620/0.9061/0.9418/0.9957`。
- [x] 正式 16-thread full-outer ratio 为
  `0.9354/0.9992/0.9840/0.9680/1.0343`，20-thread 为
  `1.0038/1.0050/0.9649/0.9771/0.9921`（顺序均为 default/high-T/stiff/high-kappa/low-T）；
  目标 regime 未形成稳定 >5% 收益。
- [x] 2-thread full-outer 目标点 spectrum delta p50/p95/max（dex）为 high-T
  `8.28e-12/1.38e-6/2.66e-6`、stiff `1.89e-11/1.86e-6/3.30e-6`、
  high-kappa `5.92e-13/8.68e-7/2.42e-6`；DN relative 为
  `3.01e-13/1.18e-11/7.80e-14`。low-T spectrum max `2.12e-6`、DN relative
  `3.12e-7`。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。formal runtime
gate 失败，不进入独立 Oracle、Sobol、coverage 或 false-safe 认证；production source
未修改。下一轮继续真正减少 propagation 工作量的 hybrid phase/Riccati/WKB standalone
prototype。

Artifacts：`docs/phase_exp_fastmath_spike_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_16t_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_20t_round10_20260917.json`。

## Current session update: 2026-09-17 (analytic branch eligibility boundary)

- [x] fresh HEAD `ba95ecc` 重跑实际 segment eligibility；严格 radiation segment 为
  `0`，严格 stiff segment 为 high-T/stiff/high-kappa/low-T
  `40202/28140/48294/0`。
- [x] `|Δsigma|<1e-4` 的 branch fraction 为 high-T/stiff/high-kappa/low-T
  `16.47%/14.31%/17.66%/6.62%`；不能据此直接声称可减少 30% tensor work。
- [x] 固定 `z_match=3.75` radiation closed-form 30-repeat candidate/exact ratio 为
  `1.149/1.156/1.104/1.137`，amplitude max relative 为
  `1.69%/0.62%/1.86%/1.81%`，因此固定 handoff 诊断拒绝。

下一步：只做按实际 sigma 连续区间识别的局部 radiation/stiff 跳跃 prototype，
先验证 output error，再决定是否值得 runtime A/B；production source 保持不变。

Artifacts：`docs/analytic_branch_eligibility_round11_20260917.json`、
`docs/radiation_closed_form_round11_20260917.json`。

## Current session update: 2026-09-17 (radiation local-branch jump)

- [x] standalone prototype 只在实际 `sigma≈4/3` 连续 run 上合并 radiation exact map，
  非候选段保持 production `_phase_segment`。
- [x] tolerance `1e-4` 的 branch segment 总数约 `7723–8914`，每 mode median run
  length 为 0；candidate/baseline runtime ratio default/high-T/stiff/high-kappa/low-T
  为 `1.022/1.070/1.020/1.045/1.195`。
- [x] amplitude p95/max（目标 high-T/stiff/high-kappa）为
  `1.51e-4/3.58e-3`、`3.35e-4/2.09e-3`、`1.51e-4/3.02e-3`；low-T max
  `5.67e-3`，不满足 error budget。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。下一步转向
`|omega'/omega^2|` adiabaticity-trigger prototype，不再调整该 branch threshold。

Artifact：`docs/radiation_local_branch_round11_20260917.json`。

## Current session update: 2026-09-17 (adiabaticity-trigger handoff)

- [x] standalone prototype 用 `|omega'/omega^2|=|1.5*sigma-1|exp(-z)` 触发
  现有 analytic tail，production source 未修改。
- [x] threshold `0.02/0.01/0.005/0.001` 均未通过 accuracy/runtime gate；以
  `0.001` 为例目标点 runtime ratio high-T/stiff/high-kappa 为
  `1.489/1.425/1.514`，spectrum max 为 `0.463/0.400/0.463`，low-T max
  `0.345`、DN relative `8.45e-3`。
- [x] 定位到 `sigma≈2/3` 使一阶 epsilon 因 `omega'≈0` 虚假变小，导致 z<5
  的 premature handoff；criterion 不能单独使用。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。下一步若继续
hybrid，必须组合 phase/curvature 或 mode-local amplitude criterion，不再单独调
epsilon threshold。

Artifact：`docs/adiabaticity_handoff_round11_20260917.json`。

## Current session update: 2026-09-17 (round11 fresh stage profile)

- [x] 当前 HEAD `c5566a8` 重新完成 high-T/stiff/high-kappa 独立 fresh profile，
  资源为 Numba=2、BLAS=1、workers=1、kink split。
- [x] warm median total 为 `9.696/10.052/10.410 ms`，tensor kernel 为
  `4.514/5.090/4.660 ms`；`fast_phi_s2_split` 为 `1.378/1.420/1.496 ms`，
  background 为 `0.670/0.691/0.783 ms`。
- [x] tensor propagation 仍为最大单项，下一候选必须围绕真实 propagation 或带
  高阶安全约束的 hybrid handoff；不再重复已拒绝 micro-optimization。

Artifacts：`docs/profile_fast_breakdown_round11_highT_20260917.json`、
`docs/profile_fast_breakdown_round11_stiff_20260917.json`、
`docs/profile_fast_breakdown_round11_high_kappa_20260917.json`。

## Current session update: 2026-09-17 (second-order adiabaticity trigger)

- [x] 在 standalone WKB carrier 中加入 `epsilon_2=|q^2+q'|/omega^2` 联合
  trigger，要求连续三个 native node；production source 未修改。
- [x] `eps=3e-4` 目标 high-T/stiff/high-kappa amplitude max 为
  `3.83e-6/4.58e-5/3.04e-5`，DN proxy 为
  `2.54e-13/5.71e-12/2.00e-12`，但 runtime ratio 为 `1.071/1.026/0.998`。
- [x] `eps=1e-3/3e-3` 仍无跨目标区稳定收益；`3e-3` high-T/high-kappa
  spectrum max 为 `1.79e-3/1.06e-3`，不满足精度预算。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。不进入 full-grid
Oracle、false-safe 或 production candidate 认证；下一步转向真正减少 propagation
work 的 nonoscillatory phase/Riccati 结构。

Artifacts：`docs/wkb_carrier_second_order_3e-4_round11_20260917.json`、
`docs/wkb_carrier_second_order_1e-3_round11_20260917.json`、
`docs/wkb_carrier_second_order_3e-3_round11_20260917.json`。

## Current session update: 2026-09-17 (kink sigma probe cache)

- [x] fresh HEAD `b8e71bf` 重新完成固定资源阶段 profile：2-thread warm median
  total default/highT/stiff/high-kappa/lowT 为 `6.057/8.910/9.994/9.678/6.388 ms`；
  tensor kernel 为 `2.442/4.131/4.682/4.375/2.909 ms`。
- [x] standalone 缓存同一 solver 调用内重复的三 probe `sigma_vec`，50-repeat
  A/B 五点 digest 与状态逐位一致；runtime ratio 为 `0.996/1.004/1.010/1.008/0.916`。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。缓存对目标
慢工况无稳定 >5% 收益，且 stiff/high-kappa 退化；production 未修改。

Artifacts：`docs/profile_fast_breakdown_round12_default_20260917.json`、
`docs/profile_fast_breakdown_round12_highT_20260917.json`、
`docs/profile_fast_breakdown_round12_stiff_20260917.json`、
`docs/profile_fast_breakdown_round12_high_kappa_20260917.json`、
`docs/profile_fast_breakdown_round12_lowT_20260917.json`、
`docs/kink_sigma_cache_round12_20260917.json`。

## Current session update: 2026-09-17 (observable-aware outer reuse)

- [x] 当前 HEAD `8d8102d` fresh 运行 observable proxy 与 15-repeat 交替 headroom
  A/B；记录 `delta log Omega`、`delta Phi`、`delta S2`、horizon shift、DN proxy。
- [x] forced reuse 对 high-T/stiff/high-kappa 的 runtime ratio 为
  `0.748/0.766/0.729`，但 spectrum max 偏差为 `3.944e-3/1.051e-3/1.547e-2`
  dex；DN 相对差虽为 `3.39e-11/3.05e-10/3.28e-11`，不能替代 spectrum gate。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS OBSERVABLE-AWARE DIAGNOSTIC`。
不放宽 outer reuse threshold，不继续调单一 observable 阈值。

Artifacts：`docs/outer_observable_proxy_round12_20260917.json`、
`docs/outer_reuse_headroom_round12_20260917.json`。

## Current session update: 2026-09-17 (WKB early-tail handoff)

- [x] fresh HEAD `b0cdf54` 上完成 standalone z=4 early-tail kernel twin；修正
  endpoint assembly 后重新运行。
- [x] candidate kernel median ratio `0.708`，但 default 对 Prüfer z=10 oracle 的
  spectrum p50/p95/max 为 `1.292e-2/1.833e-2/2.407e-2`，DN relative `7.405e-3`；
  当前 z=5 baseline DN relative 为 `5.214e-6`。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。一阶 WKB tail
correction 不足以支撑 z=4 handoff，production 未修改。

Artifact：`docs/wkb_tail_handoff_round13_20260917.json`。

## Current session update: 2026-09-17 (WKB carrier preintegration)

- [x] 在 fresh HEAD `cf183c4` 上完成 `z_match=4 -> z_tail=5` 的 standalone
  carrier 预积分 A/B；25 repeats、Numba=2、BLAS=1、workers=1。
- [x] candidate/exact kernel median ratio 为
  `0.910/0.943/1.051/0.928`（default/high-T/stiff/high-kappa），stiff 退化约
  `5.1%`；amplitude max 相对差为 `1.35e-4/1.28e-4/1.07e-4/2.07e-4`，无
  numerical failure。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。收益不跨目标
regime 稳定，且尚未值得进入独立 oracle 认证；production 未修改。

Artifact：`docs/wkb_carrier_round14_20260917.json`。

## Current session update: 2026-09-17 (primitive sigma-node no-copy)

- [x] fresh HEAD `b23190a` 上完成 formal/exact-kink TDD contract，primitive 输出
  逐位一致（`2 passed`）。
- [x] 五工况 outer 25-repeat A/B ratio 为
  `0.997/0.982/0.987/0.997/0.989`（default/high-T/stiff/high-kappa/lowT），
  digest、spectrum、DN、convergence/failure 全部一致。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。最高总耗时
改善约 `1.81%`，不满足速度门槛；production 未修改。

Artifacts：`docs/fast_phi_nocopy_round15_*.json`。

## Current session update: 2026-09-17 (FD interpolator lookup cache)

- [x] fresh HEAD `2dd7547` 上完成 FD lookup identity contract；对象身份与输出
  determinism 均保持一致。
- [x] 25-repeat 的局部 default 改善在扩大到 50 repeats 后消失；五工况 ratio 为
  `1.011/1.004/0.999/1.016/1.012`，无稳定收益。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。production 未
修改，不再重复 Python import/cache 微优化。

Artifacts：`docs/fd_lookup_cache_round16_*.json`。

## Current session update: 2026-09-17 (cross-channel overlap diagnostic)

- [x] fresh round17 profile 重新确认 tensor kernel 为主要热点。
- [x] 统计五工况 `j0`、tail endpoint 与 `(j0, tail_index)` overlap；pair reuse 为
  `0`，不具备简单 exact transfer 共享条件。

Decision: `REJECTED AS EXACT REUSE CANDIDATE / DIAGNOSTIC RETAINED`。production 未
修改，下一步转向 outer allocation/assembly 或新的 mode-local 数学路径。

Artifacts：`docs/profile_fast_breakdown_round17_*.json`、
`docs/channel_overlap_round17_20260917.json`。
## Session: 2026-09-17 (round18 outer assembly attribution)

### Actions Taken

- 先尝试 fetch `fast_v0.2`；远端 TLS handshake 两次失败，但本地 tracking ref
  与远端已知 HEAD 均为 `9173f0d41e2264926d0a5f4c458530f872a17b13`，未以旧 SHA
  替代当前状态。
- 阅读最近 commits、`progress.md`、`findings.md`、既有 profile/validation
  artifacts，并确认 raw Riccati pole audit 与 Prüfer full-grid oracle 已有结论，
  不重复已否决实验。
- 新增只读 stage attribution：warm 后每个 regime 重复 5 次，记录 preparation、
  primitive、每次 `solve_kernel` 的 `assemble` 标志、积分与插值阶段。

### Test Results

| regime | total median | primitive | solve per call | full solves/run |
|---|---:|---:|---:|---:|
| default | 6.182 ms | 0.847 ms | 2.450 ms | 1 |
| high-T | 9.094 ms | 0.588 ms | 2.201 ms | 2 |
| stiff | 9.921 ms | 0.690 ms | 2.452 ms | 2 |
| high-kappa | 9.214 ms | 0.636 ms | 2.115 ms | 2 |
| low-T | 6.323 ms | 0.787 ms | 3.123 ms | 1 |

Artifact: `docs/outer_attribution_round18_20260917.json`；script:
`scripts/benchmark_outer_attribution.py`。

### Decision

`DIAGNOSTIC CONFIRMED / NO PRODUCTION CHANGE`。high-T/stiff/high-kappa 的额外
完整 solve 已被直接观测，但 assembly shortcut 已有 2--4% 上限证据，不能满足
候选门槛；下一轮继续寻找减少 propagation arithmetic/steps 的 standalone 方向。

## Session: 2026-09-17 (round19 tail-factor cache spike)

### Actions Taken

- fresh fetch 后确认本地与 `fast/fast_v0.2` 同为 `3921c8c`，重新读取最近
  commits、profile、validation 和 rejection artifacts。
- 先写 TDD identity contract，再实现 standalone tail-factor cache twin；只将
  `ev_minus * fp_minus` 的预计算结果用于 `xf`，保持 `Th` 的原始表达式。
- 固定 Numba=2、BLAS=1、workers=1、workqueue，完成五工况 25-repeat kernel A/B。

### Test Results

| regime | candidate/base | digest | spectrum max relative | DN relative |
|---|---:|---|---:|---:|
| default | 0.9633 | mismatch | `6.70e-6` | `7.51e-12` |
| high-T | 0.9763 | mismatch | `6.13e-6` | `2.99e-13` |
| stiff | 0.8911 | mismatch | `7.60e-6` | `1.18e-11` |
| high-kappa | 0.9610 | mismatch | `5.56e-6` | `7.79e-14` |
| low-T | 0.8892 | mismatch | `4.88e-6` | `3.12e-7` |

TDD contract：`1 passed`；artifact：
`docs/tail_factor_cache_round19_20260917.json`。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS FAILED SPIKE`。候选局部 kernel 速度有改善，
但所有 regime 均未保持 bitwise digest，且本轮没有理由进入独立 oracle promotion；
production path 未修改。下一候选必须避免仅改变浮点运算顺序。

## Session: 2026-09-17 (round20 phase loop-state reuse)

### Actions Taken

- fresh fetch 后确认 HEAD `6c35c1e`；重新读取 round19 rejection、round17/20
  profiles、outer static-invariant 与 validation artifacts。
- fresh 25-repeat stage profile：default/high-T/stiff/high-kappa/lowT total
  `6.620/9.390/9.860/9.930/6.260 ms`；kernel `2.440/4.670/5.010/4.700/3.000 ms`。
- standalone twin 只复用已计算的 `zz` loop state；修正 kink endpoint 后，kernel
  和 full outer 均完成 bitwise gate。

### Test Results

- kernel 25-repeat：五工况 digest、spectrum p50/p95/max、DN 全部零差；ratio
  `0.880/0.849/0.845/0.895/0.861`。
- full outer 50-repeat、2 threads：ratio
  `0.926/0.997/0.976/0.914/0.956`；digest equal、failure 全为 null。
- full outer 25-repeat、16 threads：`0.984/1.016/0.977/0.985/0.940`；20 threads：
  `0.988/0.978/0.996/0.952/1.007`；全部 digest equal、无 failure。
- TDD contract 与 scoped regression：`56 passed, 1 deselected`（首次未固定线程
  的 pytest 环境冲突已按既有资源契约在干净 2-thread 进程重跑通过）。

### Decision

`ACCEPTED / PRODUCTION PATCH READY`。将 `while` 条件和 `z_node` 改为复用现有
`zz`，不改变任何数值表达式；production source 仅作最小修改，随后绑定新 HEAD
重跑 post-commit verification。

### Post-commit verification

- production commit：`94144c14c0a79869b2043845b83f13c18ab7e4f1`，已 push 到
  `fast_v0.2`。
- 全量 pytest：`165 passed, 6 deselected, 2 warnings`。
- post-commit fresh stage profile：default/high-T/stiff/high-kappa/lowT total
  median `6.04/9.83/10.15/9.98/6.42 ms`，均 converged、无 failure。
- artifacts：`docs/profile_fast_breakdown_round20_postcommit_*.json`，每个 JSON
  的 `commit` 字段均为 `94144c1`。
## Round 21 — H2 cross-layer endpoint reuse triage (2026-09-17)

- Fresh HEAD verified: `3e5bf7ee82d50a5901d5c66aad886ef909819d80`.
- Profiled the current `goal` path from the postcommit 25-repeat artifacts with fixed study resources.
- Candidate: pass H2 endpoint values from grid generation into `_correct_kink_background`.
- Result: rejected before code change. The full stage is only `0.10–0.18 ms` warm median across named cases; removing it cannot reach the >5% runtime gate.
- Artifact: `docs/h2_endpoint_cross_layer_round21_20260917.json`.
- Next: fresh high-T/stiff/high-kappa attribution at tensor/outer-loop granularity; do not alter the H2 return contract for this low-headroom candidate.

## Round 22 — Observable-aware outer reuse screen (2026-09-17)

- Fresh HEAD: `7742f770c1ebd948612351a66ab47d430a4fd56c`.
- Fresh 25-repeat stage profiles: high-T/stiff/high-kappa total `10.42/10.50/10.28 ms`, tensor kernel `4.48/4.74/4.68 ms`; each has two outer iterations.
- Forced-reuse A/B on named, tilt, and Sobol points: second-kernel removal saves `21–28%`, but spectrum deviations range from zero/small-update to `0.208 dex`; DN relative differences alone are not sufficient evidence of spectrum safety.
- Observable proxy screen (`delta S2`, horizon shift, `delta Phi`) separates the observed small-update rows, but remains diagnostic and has no certified DN error bound.
- Decision: **REJECTED_FOR_PRODUCTION**. No formal threshold or production behavior changed.
- Artifacts: `docs/outer_observable_gate_round22_20260917.json` plus the four named/extended profile and outer-reuse/proxy artifacts listed in the findings entry.
- Next: seek a certified predictor only if it can satisfy non-bitwise oracle/Sobol gates; otherwise continue strict-equivalence hotspot attribution.

## Round 24 — Standalone Numba Phi/S2 primitive (2026-09-17)

- Fresh HEAD: `72231d76c47732b35d484c42ce07a0c9a6fa532e`.
- Added only a standalone prototype and contract test; production source was not changed.
- Primitive arithmetic was bitwise identical and `13–21%` faster at 2 threads.
- Full outer 50-repeat ratios at 2 threads: default/high-T/stiff/high-kappa/lowT `0.950/0.968/0.976/1.007/1.003`.
- Formal 20-thread 25-repeat ratios: `1.042/0.875/0.957/0.986/0.967`; extended Sobol/tilt outputs remained digest-equal, with `sobol_006` at `1.016`.
- Decision: **REJECTED_FOR_PRODUCTION** because end-to-end improvement is not stable across the required regimes and formal default regresses.
- Artifact: `docs/phi_s2_numba_spike_round24_summary_20260917.json`.
- Next: continue strict-equivalence attribution; do not promote this helper or revisit the same primitive allocation variants.

## Round 25 — Phi/S2 workspace reuse (2026-09-17)

- Fresh HEAD before the candidate: `b22042329f8da8f6aabcc986c89fa49231b43cca`.
- TDD first: standalone workspace identity test and production buffer-lifetime
  test were red before implementation; both are now green.
- Candidate reuses only model-local arrays whose prior values are dead before the
  next `fast_phi_s2_split` call. Formulas, kink probes, node inputs, and output
  ordering are unchanged.
- Standalone full-outer A/B, 2 threads, 50 repeats, ratios
  default/high-T/stiff/high-kappa/low-T: `0.948/0.927/0.925/0.916/0.917`;
  all digest-equal and DN relative difference `0`.
- Formal 16-thread ratios: `0.877/0.907/0.896/0.915/0.886`; formal 20-thread
  ratios: `0.918/0.919/0.889/0.925/0.894`; all digest-equal.
- Extended tilt/Sobol 2-thread ratios range `0.918–0.977`; all digest-equal and
  DN relative difference `0`.

Decision: **ACCEPTED / PRODUCTION PATCH PENDING POST-COMMIT GATE**. The candidate
passes strict-equivalence and the multi-regime speed screen; production integration
is limited to a Numba fill helper plus model-local workspace allocation/reuse.

### Post-commit verification

- Production commit `bd2cff602339b386efe9cfdcb3c19056d113e6b8` was pushed to
  `fast_v0.2`.
- Fresh 25-repeat fixed-resource profiles at that SHA: total medians in ms
  default/high-T/stiff/high-kappa/low-T = `5.64/7.68/8.55/7.87/5.52`;
  all five cases converged, had no numerical failure, and each repeated digest
  was stable. The default remains above the final `<4 ms` target.
- Fresh artifacts: `docs/profile_fast_breakdown_round25_postcommit_*.json`.

## Round 26 — Outer snapshot reference reuse (2026-09-17)

- Fresh HEAD: `91696a3b53ad6668b63dc7c5713f2f17b52e1f06`.
- Hypothesis: `outer_sigma_prev`, `outer_f_hor_prev`, `exact_Nv`,
  `exact_sigma`, and `exact_f_hor` are only compared until the next
  `gen_fast`, which replaces the model-owned arrays; retaining references could
  remove redundant copies without changing the reuse test.
- TDD identity contract passed, and all A/B outputs were digest-equal with
  zero DN relative difference and no numerical failure.
- 2-thread 50-repeat ratios default/high-T/stiff/high-kappa/low-T:
  `0.991/0.970/0.987/0.971/0.954`.
- Formal 16-thread ratios: `0.989/0.931/0.974/0.957/1.006`; formal 20-thread
  ratios: `0.954/0.935/0.988/0.983/0.986`.

Decision: **REJECTED FOR PRODUCTION**. Although safe and beneficial in several
slow regimes, it does not provide a stable >5% multi-regime improvement and
slightly regresses low-T at 16 threads. Production source was restored; only
the diagnostic artifacts are retained.
## 2026-09-23 — tensor propagation optimization

- Read the pasted specification, relevant memory, existing plan, core fast solver, profiler, resource budget, benchmark docs, and prior 2026-09-17 artifacts.
- Confirmed branch/remote/base: `codex/fast_v0.2` -> `fast/fast_v0.2`, both at `465196c82b1c938af7bfce125933a4dfda22f4a2`; Git email is `2966684515@qq.com`.
- Initial 25-repeat stage profile was intentionally invalidated after detecting `kink_split=False`; no production conclusion was taken from it. Corrected canonical profile is next.
- Corrected canonical five-case profile completed with `kink_split=True`, fixed affinity/resources, and baseline SHA binding. High-T/stiff/high-kappa each require two propagation calls; P0 remains the next candidate because tensor propagation dominates.
- P0 strict repeated-exp twin passed named kernel bitwise gates at 2/16/20 threads but failed stable full-outer runtime gate; artifacts were saved and production code remained unchanged. Next: P1 counted assembly-state twin.
- P1 counted assembly-state twin passed output/trace bitwise gates and showed strong 2-thread kernel/outer savings, but formal 20-thread 50-repeat full-outer regressed highT/stiff/high-kappa. Candidate rejected; production source remains unchanged. Next: LLVM/ASM audit.
- LLVM/ASM audit initially hit cache-disabled inspection; `Dispatcher.recompile()` fixed the audit path. It confirmed duplicate phase exp sites and remaining integer division/remainder operations in generic `solve_kernel`; next is an isolated guarded canonical-fast specialization twin.
- Canonical-fast specialization passed named bitwise gates and improved most formal points, but canonical 20-thread 50-repeat high-kappa regressed 1.29%; no production dispatch was added. Remaining work is fresh verification, review, commit and push of the diagnostic artifacts.
- Verification found all three new candidate tests passing. Existing scoped fast tests have 3 baseline failures caused by stale 19–23-argument solve-kernel monkeypatches and an unchanged shared-N_eff guard in the eval-frequency transition-refine case; these remain unmodified and will be disclosed.
- Autonomous continuation request received. Fresh fetch at the pushed Round 28 commit `1750772a1c70cb3f99edc054fd7841220bdf3601`; no tracked source drift. The next loop must re-profile canonical `kink_split=True`, search distinct numerical/HPC hypotheses, and keep the next concrete experiment in this ledger before execution.
- Round 29 fresh canonical profile completed at the fetched HEAD. Artifacts `docs/profile_fast_breakdown_round29_20260923_{default,lowT,highT,stiff,high_kappa}.json` are bound to `1750772a1c70cb3f99edc054fd7841220bdf3601`; the tensor solve remains the only stage with enough Amdahl headroom for the next experiment.
- Round 29 H2 grouped SoA/AoSoA spike completed as standalone. TDD red/green passed; all five named kernels were bitwise equal but 31–35% slower, and full outer was 18.6–24.3% slower with p95 regression. Candidate LLVM showed no vector double lanes, so it was rejected without formal-thread expansion. Next direction must be materially different from grouped lockstep, not a bucket-width retune.
- Round 30 coarse fixed-point predictor audit: the first prototype had a harness-only PCHIP measure defect (missing `ln(10)`); TDD caught it and the production measure was restored. The corrected gain-1 first iterate still failed the `1e-3 dex` spectrum gate in all five named cases (`3.342e-3`–`5.905e-3 dex`) because production may return a spectrum computed at the original `DN_eff` while only updating the final scalar `DN_eff`. Predictor rejected; production unchanged. Next: residual-controlled local phase/adiabatic transfer prototype with exact fallback.
- Round 31 residual-controlled two-step transfer: threshold-0 exact fallback passed bitwise; threshold-1e-3 passed named/edge/Sobol reference screens and gave 2-thread outer ratios `0.838–0.898`, but formal 16/20-thread full-outer ratios `0.921–0.991` / `0.969–1.021` failed stable >5% performance, including high-kappa regression. Rejected for production; next selected experiment is fourth-order two-exponential commutator-free Magnus over four native intervals with exact fallback.
- Round 32 CF4 Gaussian-node two-exponential Magnus screen: threshold-zero full output remained bitwise equal; independent oracle residuals stayed near baseline, but the full-assembly kernel ratios at threshold 1e-3 were `0.790–1.069` across named regimes with high-T regression. Rejected before formal outer expansion. Next: residual-certified phase-function Chebyshev representation feasibility spike, distinct from raw Riccati and WKB handoff paths.
- Round 33 fresh profile at `00e2bfd` completed. A read-only Kummer positivity screen was finite and positive for sampled smooth post-horizon modes in all five named regimes (`4/4` usable each); this only justifies building the next standalone phase-function residual prototype, not production promotion.
- Round 33 phase-function Chebyshev prototype completed: independent-solution Kummer amplitude stayed positive, but degree-12 fits over useful z-windows had residual `12.5–57.6`; only a too-short window passed. Rejected before runtime; next is boundary-conditioned Kummer phase solve, not a degree/threshold retune.
- Round 34 started after fresh fetch at `e40a2d0f5373843bfe3fc66a0173741be2346882`. The selected hypothesis is a right-boundary-WKB-conditioned complex solution integrated backward on the same smooth post-horizon windows. It is distinct from the rejected arbitrary-basis Chebyshev fit; production remains unchanged.
- Round 34 fresh profile medians were `5.528/5.557/10.170/9.533/10.319 ms` and tensor medians `2.336/2.754/4.767/5.010/5.003 ms` for default/lowT/highT/stiff/high-kappa under fixed 2-thread resources; hard cases retained two propagation calls. The boundary-conditioned WKB phase screen failed with Kummer residual `1.0–41.9` and Cartesian transfer error `0.25–1.60`, so it was rejected before timing. Next concrete experiment is the certified small-angle polynomial transfer twin, not further phase-function tuning.
- Round 35 fresh profile at `27b9fc5` measured total/tensor medians `5.637/2.330`, `5.479/2.887`, `9.053/4.250`, `8.478/4.633`, `9.721/4.462` ms for default/lowT/highT/stiff/high-kappa. The small-angle `sin/cos` twin passed local/reference/guard gates but failed stable full end-to-end performance across 2/16/20-thread resources; it is rejected. Next concrete experiment is range-reduced per-transfer `exp(z)` with exact production subdivision decisions.
- Round 36 fresh profile at `f94785e` measured total/tensor medians `6.209/2.750`, `5.543/3.125`, `9.157/4.486`, `9.082/4.908`, `9.244/4.649` ms. Range-reduced per-transfer exp passed local precision but ran `1.986x` slower in the isolated kernel and was rejected before oracle/full outer. Next concrete experiment is tangent/sensitivity-assisted outer correction with explicit final-spectrum semantics.
- Round 37 fresh profile at `c2ff97a` measured total/tensor medians `5.604/2.462`, `5.479/2.913`, `8.237/4.355`, `8.246/4.515`, `8.321/4.341` ms. Full-state tangent matched one-step finite differences but became non-finite over long propagation and cost `41.1x` even without assembly/tail; rejected. Next concrete experiment is bounded observable/adjoint sensitivity, not another full tangent state.
- Round 38 fresh profile at `85cec7a` measured total/tensor medians `5.419/2.429`, `5.247/2.902`, `8.280/4.321`, `8.194/4.508`, `8.329/4.343` ms. Scalar adjoint transpose was locally exact but non-finite and `24.3x` slower without assembly/tail; it also cannot satisfy full-spectrum output. Rejected. Next concrete experiment is combined transfer-map lookup/interpolation.
- Round 39 fresh profile at `0d0d77c` measured total/tensor medians `5.337/2.489`, `5.127/2.757`, `7.823/4.222`, `8.608/4.723`, `8.086/4.486` ms. The 4097-node/32-row cubic combined transfer table failed local error (`9.73e-10` vs `1e-11`) before kernel timing. Rejected; next concrete experiment is Levin-collocation nonoscillatory phase feasibility.
- Round 40 fresh profile at `1e88d31` measured total/tensor medians `5.274/2.277`, `5.192/2.760`, `7.952/4.262`, `8.362/4.698`, `8.106/4.330` ms. Levin interaction-envelope compression failed (`1.36e10–1.97e10` fit error; `0.0058–0.0562` transfer error) before runtime. Rejected; next concrete experiment is explicit four-mode lockstep LLVM/vectorization feasibility.
- Round 41 fresh profile at `633d3670f022b40cc7f6e2ca4910fa6a1c4babc4` measured total/tensor medians `5.470/2.338`, `5.041/2.696`, `7.795/4.040`, `8.087/4.543`, `8.030/4.460` ms for default/lowT/highT/stiff/high-kappa with `kink_split=true`. The explicit four-lane, `fastmath=false` LLVM probe was deterministic but exposed no `<N x double>` vector lane. Rejected before full kernel/outer timing; artifact `docs/fixed_lane_probe_round41_20260923.json`. Next concrete experiment must be a distinct real-work reduction path, not a grouped/unrolled lane retune.
- Round 42 fresh profile at `7a39df3d2cc2b27d96749989521ad0bb57e06f75` measured total/tensor medians `5.419/2.359`, `5.131/2.786`, `8.225/4.203`, `8.193/4.534`, `7.797/4.255` ms for default/lowT/highT/stiff/high-kappa. Sparse propagation of 39/76-77 modes reduced fixed-DN candidate medians to `0.722/0.744/0.644/0.661/0.635` of baseline, but scalar PCHIP reconstruction failed fixed-DN spectrum/DN gates (`0.282–0.907 dex`, `2.19e-3–5.18e-3` DN relative). Rejected before outer/oracle; artifact `docs/sparse_frequency_round42_20260923.json`. Next concrete experiment is phase-aware complex-state reconstruction, not scalar-spectrum stride tuning.
- Round 43 fresh profile at `ae28e6cf71df2ab3e2f500804b1d739530820b2a` measured total/tensor medians `5.239/2.372`, `5.014/2.807`, `7.656/4.041`, `7.849/4.478`, `8.081/4.102` ms. Phase-aware Cartesian-state reconstruction at the same 39/76-77 modes remained finite but failed fixed-DN spectrum/DN gates (`1.734–2.478 dex`, `1.62e-2–4.04e-2` DN relative). Rejected before oracle/outer; artifact `docs/phase_aware_sparse_frequency_round43_20260923.json`. Next concrete experiment is tail-assembly factor attribution/hoisting, gated by Amdahl.
- Round 44 fresh profile at `873b02483d05cc0fa52180428cc1d74a6340da5d` measured total/tensor medians `5.374/2.330`, `5.218/2.724`, `9.320/4.357`, `8.045/4.402`, `8.182/4.434` ms and total p95 `6.387/6.518/11.243/9.636/9.709` ms. Tail assembly attribution found only `5.92%–8.04%` full-kernel share; common-factor reuse was non-bitwise and its estimated end-to-end ceiling was `1.6%–3.7%`. Rejected by Amdahl before full twin/outer; artifact `docs/tail_assembly_attribution_round44_20260923.json`. Next concrete experiment must target propagation-level work, not tail/layout micro-optimization.
- Round 45 fresh profile at `724b765881662acd4faee13b714f746cbe6ba9fa` measured total/tensor medians `5.608/2.523`, `5.432/2.769`, `7.906/4.048`, `8.564/4.631`, `8.469/4.566` ms and total p95 `6.868/6.388/9.327/10.494/9.655` ms. Cayley-regularized projective Riccati stayed finite with `1.2e-15–6.9e-15` power error but cost `2.22–2.33x` per transfer and reduced no work. Rejected before full kernel; artifact `docs/cayley_riccati_round45_20260923.json`. Next concrete experiment must have an explicit transfer-work reduction mechanism.
- Round 46 fresh profile at `5e33e5c375fb869d8d9b48b4b3973ac36ba31b54` completed under the canonical fixed resources. Phase-increment recurrence feasibility passed synthetic and actual z-path screens: 20 named-mode records were finite with power relative error `8.53e-5–5.59e-4`, while synthetic transfer cost was `0.45–0.47x`. Not production yet; full kernel twin with exact fallback is next. Artifact `docs/phase_increment_recurrence_round46_20260923.json`.
- Round 47 completed the full phase-increment kernel twin. Fixed-DN and named/edge/Sobol/oracle correctness gates passed at experimental tolerance, but alternating full-outer ratios were unstable: 2-thread `0.948–0.990`, 16-thread `0.921–1.013`, and 20-thread `0.983–1.021`, with formal high-T/high-kappa regressions. Rejected; production unchanged. Next selected pure-math experiment: exact algebraic state-basis transfer rewrite, not AI, learned prediction, phase-recurrence retuning, or exp-hoist reopening.
- Round 47 correction: the first timing run had an experimental re-anchor age-counter bug and is superseded. Corrected full-outer ratios are 2-thread `0.981–0.994`, 16-thread `0.972–1.016`, and 20-thread `0.981–1.022`; the candidate remains rejected with production unchanged. Fixed-DN/oracle values were also rerun after the fix.
- Round 48 completed the pure algebraic state-basis screen at new HEAD `0b4b208`. Long-chain power error was `1.16e-15–7.42e-15`, but loop ratios were only `0.983/0.991/1.000`; Amdahl bounds the full-path gain below 1%, so it was rejected before full-kernel construction. Next selected pure-math direction: native-grid Prüfer/WKB hybrid feasibility.
- Round 49 completed native-grid Prüfer transfer feasibility at `2485e99`. Power error was `4.81e-16–3.05e-15`, but loop ratios were `1.656–2.129` due to per-step `sin/cos/atan2/sqrt`; rejected before full kernel. Next selected strict-equivalence arithmetic: verify and test repeated `w*si` coefficient reuse only if LLVM/ASM leaves it unmerged.
- Round 50 completed at fresh HEAD `d230d93`. The canonical profile was noisy (`total/tensor` medians default `6.230/2.406`, lowT `7.927/3.944`, highT `11.476/2.591`, stiff `8.229/2.029`, high-kappa `7.956/1.635` ms). LLVM already CSEs the repeated `w*si` coefficient in `scaled_step`, so the strict-equivalence candidate was rejected as already optimized without timing. Next concrete experiment: a tiny residual-certified two-interval uniform-asymptotic defect screen, with production unchanged and no AI/learned model.
- Round 51 completed at fresh HEAD `b929cb7`. Fresh profile medians were total `4.46/3.93/6.26/6.38/6.57 ms` and tensor `1.78/2.03/1.41/1.78/1.43 ms` for default/lowT/highT/stiff/high-kappa; hard cases still used two propagations. The log-mean two-interval block passed local TDD and had ideal transfer ratio `0.768`, but the embedded matrix defect guard cost `3.76x` baseline (`0.1865/0.0496 ms`) and accepted no timing blocks. Rejected by Amdahl before full kernel/outer; production unchanged. Artifact `docs/uniform_asymptotic_defect_round51_20260923.json`. Next concrete experiment: cheap conservative endpoint defect bound precomputation, only if it retains the same mathematical safety gate. No AI/learned model.
- CI audit: GitHub Actions run #219 for `b929cb7` had canonical failure while static/package/cobaya and compatibility 3.9-3.13 passed. Local workflow-equivalent canonical rerun had 197 passed, 4 pre-existing failures, 6 deselected; tracked production/tests had no diff, so no unrelated CI fix was made.
- Round 52 completed at fresh HEAD `8e8b710`. Fresh profile medians were total/tensor `5.11/1.86`, `4.86/2.16`, `7.04/1.46`, `6.98/1.88`, `7.30/1.57` ms for default/lowT/highT/stiff/high-kappa. The conservative endpoint defect bound covered the matrix defect (actual/bound <= `0.273`) but its guard chain was `1.565x` baseline at `1e-8` and `1.544x` at `1e-6`, with little block acceptance. Rejected by Amdahl before full kernel/outer; production unchanged. Artifact `docs/uniform_defect_bound_round52_20260923.json`. Next concrete experiment: scaled Riccati/WKB adaptive representation with projective normalization, tiny finite/defect screen only; no AI/learned model.
- Round 53 completed at fresh HEAD `2969196`. Fresh profile was noisy (`total/tensor` default/lowT/highT/stiff/high-kappa `8.62/3.44`, `8.22/4.32`, `12.74/2.65`, `7.38/2.21`, `9.72/2.07` ms). The scaled Riccati/log-amplitude prototype was finite on sampled paths with max power relative `5.32e-6`, but same-subdivision timing was `2.005x` baseline (`0.0842/0.0420 ms`, p95 `0.0888 ms`); rejected by Amdahl before full kernel/outer. Artifact `docs/scaled_riccati_round53_20260923.json`. Next concrete experiment: exact two-by-two block action with a precomputed schedule only if it removes real phase substeps. No AI/learned model.
- Round 54 completed at fresh HEAD `6f14d18`. Fresh profile medians were total/tensor `5.04/1.90`, `6.67/2.70`, `9.71/1.70`, `8.62/2.06`, `8.30/1.70` ms for default/lowT/highT/stiff/high-kappa. Exact block composition retained all `4096` phase maps and added `2048` matrix products in sampled schedules, so it was rejected by Amdahl without kernel timing. Artifact `docs/exact_block_schedule_round54_20260923.json`. Next concrete experiment: fresh LLVM/ASM audit for vector lanes or invariant loads; no AI/learned model.
- Round 55 completed at fresh HEAD `7e005f1`. Fresh profile medians were total/tensor `5.00/1.99`, `4.67/2.23`, `6.41/1.43`, `6.58/1.86`, `7.03/1.51` ms for default/lowT/highT/stiff/high-kappa. LLVM/ASM again showed no vector lanes or new strict-equivalence cost; generic solve-kernel div/mod remains already covered by the rejected counted-state path. Rejected as already covered without timing; production unchanged. Artifacts `docs/profile_fast_breakdown_round55_20260923.json` and `docs/audit_fast_kernel_llvm_round55_20260923.json`. Next concrete experiment: bounded analytic outer derivative/predictor screen, exact final-spectrum correction, no AI/learned model.
- Round 56 completed at fresh HEAD `4b483ce`. Fresh profile medians were total/tensor `4.66/1.83`, `4.19/2.09`, `6.37/1.42`, `6.73/1.89`, `11.43/2.05` ms for default/lowT/highT/stiff/high-kappa. Safeguarded exact-map secant residuals were near zero, but final `DN` differed by `0.00518–0.00831` and spectrum by `0.00334–0.00590 dex`; each candidate used three exact maps. Rejected by correctness before production; artifact `docs/safeguarded_outer_secant_round56_20260923.json`. Next concrete experiment: representation that actually removes phase substeps, no AI/learned model.

- Round 57 completed at fresh HEAD `493add8`. Fresh 25-repeat canonical totals were `6.711/7.318/11.246/11.731/11.309 ms` median for default/lowT/highT/stiff/high-kappa, with tensor shares `0.41/0.53/0.48/0.52/0.48`; `kink_split=true` and the fixed resource contract were confirmed.
- The Bessel variable-coefficient transfer passed standalone DOP853 accuracy (`1.4e-15–1.0e-14` power error) and finite named-path gates, but was `166.7–188.7x` slower than the native midpoint transfer. Rejected by Amdahl before full kernel/full outer; production unchanged.
- CI run 225: static/package/cobaya and Python 3.9–3.13 compatibility passed; canonical failed with the known baseline failures.
- Next: tiny asymptotic Bessel/WKB coefficient-recurrence operation-count screen, only if it can remove the special-function cost by >=100x. No AI/learned model.

- Round 58 completed at fresh HEAD `6cb3bc7`. The 25-repeat full-path profile
  was host-noisy (`10.066/7.771/14.792/17.094/13.325 ms` median for
  default/lowT/highT/stiff/high-kappa), but digests matched Round 57.
- Full-native WKB phase-integral prototype failed correctness on all five
  actual regimes: the first 32-interval mode probe crossed a turning-point /
  forbidden interval and became non-finite. Rejected before timing and full
  kernel; production unchanged.
- Next: tiny uniform turning-point Airy/Bessel connection feasibility screen,
  only if it can prove finite horizon crossing. No AI/learned model.

- Round 59 completed at `8fd4975`: the simple-turning-point Airy connection
  was finite on all five actual canonical crossing intervals and had local
  DOP853 power error `8.10e-11–8.12e-11` versus midpoint `1.61e-10–1.62e-10`.
- Only one native interval per mode was replaced, so the building block passed
  local feasibility but was rejected for production by Amdahl; no full-kernel
  promotion was attempted and production remains unchanged.
- Next: operation-count screen for a composite forbidden-to-oscillatory
  Airy/WKB representation. No AI/learned model.

- Round 60 completed at `c02ed31`: the 4000-interval composite contained
  `2852–2892` forbidden, `1–3` Airy, and `1106–1147` WKB intervals.
- The candidate stayed finite but accumulated catastrophic power error versus
  DOP853 (`2.66e4–1.18e9`), while midpoint production remained `<=1.31e-4`.
  Rejected by correctness before production/full user-facing benchmark.
- Next: full-chain residual-certified phase-function feasibility, only if an
  independent residual bound is established first. No AI/learned model.

- Round 61 completed at `7090dbd`: fresh canonical total medians/p95 were
  `6.723/8.355`, `7.364/8.128`, `10.712/12.554`, `13.646/26.556`, and
  `11.164/12.535 ms` for default/lowT/highT/stiff/high-kappa; digests were
  unchanged.
- The exp coefficient recurrence was finite and nearly bitwise on actual
  paths (`0--1.71e-15` power error), but its transfer ratios were `1.25--1.78x`
  baseline, so it was rejected by Amdahl before full-kernel/full-outer timing.
- CI run 229 for the pushed Round 60 HEAD failed only in canonical;
  static/package/cobaya and compatibility 3.9--3.13 passed. Production source
  remains unchanged.
- [ ] Next: direct residual-controlled phase-function integration over a full
  oscillatory chain, distinct from the rejected Chebyshev and boundary-WKB
  representations. No AI/learned model.

- Round 62 completed at `40a07c9`: fresh totals/p95 were
  `5.888/6.759`, `5.395/6.437`, `12.983/16.719`, `12.667/14.439`, and
  `9.754/11.598 ms`; tensor medians were `2.580/3.062/6.068/6.622/4.676 ms`.
- Direct Ermakov--Kummer phase integration stayed finite and step-stable, but
  transfer error was `2.61--2.63` and cost `12.6--14.7x`; rejected by
  correctness and Amdahl before production integration.
- Production source remains unchanged. Round 61 CI run 229 had only
  canonical failure; the Round 62 push will create the next CI run.
- [ ] Next: local polynomial-coefficient macrostep transfer with an
  independent defect bound and no global WKB amplitude transport. No AI.

- Round 63 completed at `e63fc9e`: fresh canonical total medians/p95 were
  `5.601/6.562`, `6.150/8.198`, `9.389/11.472`, `12.872/14.588`, and
  `12.704/15.391 ms`; tensor medians were `2.619/3.029/4.587/6.756/5.737 ms`.
- Pairwise commutator Magnus halved the map count but had component error
  `1.23e-3--5.34e-3`, power error `5.48e-6--1.46e-5`, and warmed ratios
  `1.69--1.91x`; rejected by correctness and Amdahl. Production unchanged.
- Round 62 push is confirmed at remote SHA `e63fc9e`; the new CI run was not
  yet visible, while the latest visible run #229 had canonical failure only.
- [ ] Next: defect-canceling local polynomial transfer reusing the production
  trigonometric map, with an independent bound before timing. No AI.

- Local CI-equivalent full pytest at Round 63: `217 passed, 4 failed,
  6 deselected`; all four failures match the documented baseline and all
  Round 61--63 tests passed. No unrelated fix was applied.

- Round 64 completed at `c219a24`: fresh totals/p95 were
  `5.592/7.602`, `5.998/7.514`, `12.464/14.307`, `13.216/14.735`, and
  `12.130/15.461 ms`; tensor medians were `2.351/2.985/5.881/6.673/5.981 ms`.
- Pairwise Duhamel correction halved map count but produced power error
  `1.40e-3--7.60e-3`, component error `8.46e-4--1.78e-2`, and ratios
  `1.07--4.73x`; rejected by correctness/performance. Production unchanged.
- [ ] Next: bounded deterministic outer predictor/corrector, only if it can
  eliminate the second propagation while preserving the full spectrum.
  No AI.

- Round 65 completed at `2f6cf08`: fresh five-case profiles used fixed
  resources and `kink_split=true`; medians/p95 were
  `5.851/177.853`, `6.860/219.624`, `8.644/193.029`, `9.815/202.410`, and
  `10.495/212.629 ms` for default/lowT/highT/stiff/high-kappa. The tensor
  medians were `2.889/3.010/4.367/5.232/4.889 ms`.
- Frozen-transfer envelope rescaling passed its two TDD tests but was
  ineligible on four cases because discrete horizon-start indices changed
  between outer maps; lowT lacked a comparable second snapshot. Rejected by
  correctness before formal E2E timing; production unchanged.
- CI remains externally stale: current remote SHA is `2f6cf08`, but GitHub
  still only shows run #229 with canonical failure and all other jobs passed.
- [ ] Next: pre-kernel interval certificate that handles discrete `j0`
  changes, or reject the outer-predictor family and move to a new certified
  transfer representation. No AI/learned model.

- CI fix pushed as `6efc945`: stale kink-extended `solve_kernel` test mocks
  were updated, and the existing fast-settings fixture is now autouse so
  standalone spikes cannot leak `_KINK_SPLIT` into later tests.
- Full local canonical equivalent now passes: `224 passed, 6 deselected, 2
  warnings`; static, manifest, package-boundary, and wheel smoke checks also
  pass. Remote Actions is not yet showing the new run, so remote green remains
  unverified.

## CI follow-up — Round 62 cross-platform diagnostic (2026-09-23)

- [x] Retrieved the protected failure logs with the repository Actions
  account: run #235 (`6efc945`) and run #236 (`7b373d2`) each had exactly one
  canonical failure, `test_direct_phase_function_is_finite_and_stable_in_an_isolated_process`.
- [x] Distinguished the failure from production behavior: Ubuntu rejected the
  matrix-reconstruction norm stability metric, while the direct Ermakov state
  remained finite and stable. The matrix norm is sensitive to endpoint matrix
  inversion and BLAS platform differences.
- [x] Changed only the standalone diagnostic gate to compare the endpoint
  `(rho, rho_prime, theta)` state at full versus half max step; retained the
  old matrix metric in the JSON artifact. No production source or numerical
  contract changed.
- [x] Focused test, full resource-capped canonical pytest (`224 passed, 6
  deselected, 2 warnings`), Ruff, compileall, and diff check pass locally.
- [x] Pushed as `e4d2a1f`; Actions run #237 is fully green for canonical,
  static, package, Cobaya, and compatibility 3.9--3.13. No AI/learned model.
- [ ] Resume the next standalone mathematical optimization experiment only
  after recording a fresh HEAD/profile and the ranked candidate pool.

## Round 66 — event-aware outer tangent endpoint screen (2026-09-23)

- [x] Reconfirmed local/remote HEAD `8b56cb7` and refreshed the five-case
  fixed-resource profile; formal `kink_split=true` and 2-thread/workqueue
  resources were recorded.
- [x] Generated multiple hypotheses after Amdahl review. Selected an
  event-aware tangent correction because it was the only remaining candidate
  that could remove a hard-case second propagation without changing the
  transfer equation.
- [x] Standalone prototype and child-process TDD gate passed. Default mode 20
  endpoint relative error was `1.84e-15`; tangent/direct runtime was `1.74x`.
  Low-T had no comparable second snapshot; high-T/stiff/high-kappa changed
  the discrete kink representation and failed closed.
- [x] Rejected before full-outer timing on correctness coverage and Amdahl;
  production source unchanged. Artifact:
  `docs/outer_sensitivity_round66_20260923.json`.
- [ ] Next: implement only a tiny pre-kernel monotone transfer-perturbation
  certificate with explicit event guards; reject immediately if its bound is
  not informative. No AI/learned model.

## CI follow-up — Round 66 optional dependency (2026-09-23)

- [x] Remote run #239 exposed one new issue: canonical imported `psutil`,
  absent from declared Ubuntu dev dependencies; all other eight jobs passed.
- [x] Removed the hard dependency in `39a31ec`; POSIX uses
  `os.sched_setaffinity`, Windows uses an optional fallback. Focused test,
  Ruff, compileall, and diff checks pass locally.
- [x] Current HEAD `ef8f293` run #241 is fully green across all nine jobs;
  the optional-dependency failure is closed.
- [ ] Continue with the pre-kernel certificate experiment.

## Round 67 — outer transfer perturbation certificate (2026-09-23)

- [x] Reconfirmed HEAD `f064599` against `fast/fast_v0.2`; latest Actions run
  #242 is fully green. Reused the fresh Round66 2-thread/workqueue,
  affinity-`[0,1]`, BLAS-one profile and its Amdahl headroom assessment.
- [x] Built a standalone norm-bound prototype using the production midpoint
  transfer matrices. The bound explicitly checks frequency start index, kink
  representation, phase subdivision, and finite values; no production source
  was changed.
- [x] TDD and Ruff passed (`1 passed`; all checks passed). The five-case run
  produced `default` 74/76 event-stable modes but an infinite Frobenius bound
  after thousands of transfers; low-T had no second snapshot and the other
  hard cases changed kink representation. Direct default endpoint perturbation
  was around `1e-10` for the first modes, so the infinity is certificate
  looseness, not evidence of a large physical difference.
- [x] Rejected before full-outer timing as
  `REJECTED_BY_VACUOUS_BOUND_OR_EVENT_CHANGE`; artifact:
  `docs/outer_certificate_round67_20260923.json`.
- [ ] Next concrete experiment: a standalone residual-certified transfer
  representation with an invariant residual/energy screen, not an outer
  reuse predictor. No AI/learned model.

## CI follow-up — no-copy standalone cross-backend assertion (2026-09-23)

- [x] Read run #243 logs. Only canonical failed, in the two existing
  `test_fast_phi_nocopy_spike.py` bitwise assertions; all other jobs passed.
- [x] Confirmed the failure is isolated to NumPy `cumsum` versus production
  Numba sequential accumulation on Ubuntu, and is unrelated to Round67 or
  production fast-path code.
- [x] Changed only the standalone spike tests to use a documented bounded
  forward-error gate (`rtol=5e-12`, `atol=1e-15`, `equal_nan=False`) while
  retaining exact kink metadata checks. Focused tests: `3 passed`; Ruff clean.
- [x] Pushed the cross-backend no-copy assertion repair at `c7d4766`; GitHub Actions run `35839494161` passed all nine jobs.

## Documentation freshness pass — 2026-09-24

- Revalidated current fast mode from source: `USER_FAST_PROFILES=('fast',)`;
  only the goal-kink-hybrid profile is user-facing. Fetched `fast_v0.2` and
  confirmed current HEAD `3c33668` is the remote target; source code is
  unchanged from the Round 68 solver profile commit `c7d4766`.
- Refreshed English/Chinese root READMEs, package description, docs indexes,
  experiment catalogue, benchmark summary, numerical-method/physics/accuracy/
  validation/Cobaya/reproducibility descriptions, and archival status labels.
- Corrected dated parameter-screen counts to the source JSON (254/146) and
  MCMC per-step comparison range to the report table (2.9–4.9x). Preserved
  original experiment measurements and commit provenance.
- Updated historical validation generator descriptions and the manifest's
  evidence-scope wording. The validation numbers themselves were not rerun.
- Refreshed the local-only upstream comparison snapshot to HEAD `3c33668`,
  corrected its ahead count to 287, and added `README_zh.md`; it remains outside
  the shared commit scope.
- Remaining: audit all tracked Markdown/generator wording, verify README pairs,
  relative links and JSON, inspect the complete diff, then commit/push only the
  shared reviewed paths to `fast_v0.2`.

## CI follow-up — MCMC benchmark Ruff repair (2026-09-24)

- Latest failed Actions runs `36006258849`, `36006594124`, and `36006595813`
  had only the `static` job fail. The reported issues were all in
  `scripts/benchmark_mcmc_sagenet_compare.py`: import order, one unused import,
  one unused local, a local import order, and an unnecessary f-string prefix.
- Reordered imports, removed only unused bindings, and kept the report content
  and benchmark logic unchanged. Local Ruff and Python compilation pass.
- Pushed as `4c59d80`; Actions run `36007199461` passed canonical, static,
  package, Cobaya, and Python 3.9–3.13 compatibility jobs.

## Repository documentation editorial audit — 2026-09-25

- Refined 18 tracked guides and reports and pushed `5165ec76c6a63bb2af67764c7b39493da2e8b7c3` to `fast/fast_v0.2`.
- Actions run `36084129180` completed successfully: canonical, static, package, Cobaya, and all Python 3.9–3.13 compatibility jobs passed.
- At this first-pass checkpoint, the full prose review was still underway. The local audit matrix tracks all 46 tracked Markdown files and the local-only comparison dossier; raw JSON/JSONL and numeric `.txt` evidence remains unchanged.
- Existing uncommitted `findings.md` entries and local performance artifacts were not staged. The upstream comparison dossier was refreshed to HEAD `5165ec7` and remains local-only.

## Session continuation — 2026-11-05

- Read the user-referenced research specification and current skill/planning instructions.
- Confirmed repository branch `fast_v0.2`, clean HEAD `16787b798902ebfcad54e85701a573ea27da31c4`, remote branch equality, and final CI run `41969511573` success (9/9 jobs).
- Applied latest user steering: do not optimize Reference; use it only as an independent validator while continuing Fast-domain work.
- Restored cumulative plan/findings/progress and added a Fast-only continuation checkpoint. Next: inspect pressure artifacts and current Fast implementation, then isolate one unresolved failure.

- Fast-only root cause narrowed: adaptive grid termination uses max node-local curvature contribution; high-blue same-background nested Prüfer data show global PCHIP integral changes exceeding both the requested target and the sum of the local indicators.
- No source change yet. Next action: use systematic debugging and TDD to test a proposed embedded global integral estimator before modifying solver behavior.

## 2026-11-05 — adaptive integral stopping-rule continuation

- Re-read the current Fast adaptive production branch, generic midpoint helper, pressure report, and archived positive-cr full-grid validation artifact.
- Clarified that formal Fast adaptation has no actual midpoint check; the generic helper does, but its previous Reference experiments were expensive and remained uncertified. The candidate under consideration is a global integral delta, not pointwise interpolation-error refinement.
- No source or test changes made. No tests run. Reference implementation remains untouched. Existing modifications in `task_plan.md`, `findings.md`, and `progress.md` were preserved and extended.
- Next: test the global-integral estimator contract and measure its expected propagation work before any production implementation.

- Ran a read-only embedded global-integral probe on the archived positive-cr Prüfer spectra. Using the denser nonnested map to interpolate 255 midpoint values gave a `1.7701e-4` relative PCHIP integral change; the archived independent two-grid `DN_gw` change is `9.9973e-5`. The single-point estimate is conservative here, but the midpoint values are interpolated, not direct, and the grids are nonnested; no accuracy guarantee follows.
- Next action: run a small actual-Fast midpoint propagation diagnostic and measure its cost before implementing the estimator. No production source changed and no test was run.

## 2026-11-05 — actual Fast midpoint probe

- Added standalone diagnostic `scripts/probe_adaptive_embedded_integral.py`; it replays the public adaptive Fast solve, checks the archived grid/root, freezes that converged model background, then propagates true interval midpoints through the production Fast kernel.
- Reproduced both existing artifacts. Positive-cr: 8.5425e-5 embedded global relative delta vs 9.9973e-5 archived Prüfer adaptive-to-res2 change; 256→511 frequencies raised kernel median from 145.5 ms to 300.6 ms (2.07x). High-blue cr=0: 7.2628e-5 vs 7.9403e-5; 73→145 frequencies raised median 8.55 ms to 16.49 ms (1.93x). Base spectra match the stored Fast results within 5.30e-4 dex and 1.60e-5 dex, respectively.
- Both probes show slight underestimation and roughly doubled kernel work. This is candidate evidence only, not a bound or production acceptance. Reference code was not changed; no production solver code changed.
- `ruff check scripts/probe_adaptive_embedded_integral.py` and `python -m py_compile scripts/probe_adaptive_embedded_integral.py` passed. The script itself was run on both fixed-background cases.
- Harness errors (wrong module-default config; missing optional artifact metadata) were corrected; final runs verify archived frequency grid and `DN_eff` replay before measuring.
- Next: test an integral-contribution panel selector for reducing midpoint propagation cost while preserving the global delta signal; then evaluate additional held-out cases and decide accept/reject.

## 2026-11-05 — panel-ranked midpoint feasibility

- Extended the fixed-background diagnostic with per-panel PCHIP integral changes and a comparison between oracle ranking (actual midpoint panel deltas) and the existing curvature×panel-weight ranking (available before midpoint solves).
- Signed cancellation is material: absolute panel-change mass is 1.65× the net delta on high-blue and 9.20× on positive-cr. Ranking by current curvature×weight and sampling the top quarter captured 80.36% and 99.989% of the oracle absolute mass, respectively.
- In those diagnostic subsets, `DN_gw` discrepancy vs archived dense Prüfer fell to `2.56e-5` (high-blue) and `2.38e-5` (positive-cr); kernel medians were about 1.17× and 1.43× base. Oracle ranking is explicitly offline/optimistic because it uses all midpoint values.
- The run's one-off Fast map checks passed; Ruff and `py_compile` passed. No production source changed.
- CI for commit `67b73910ddaa3fe2049412b78af8af8f634141d3` passed all nine jobs in run `41972973367`.
- Next: test cumulative curvature-budget selection on more held-out points, including neighboring and boundary cases, and track spectrum error as well as `DN_gw`. If it stays promising, implement behind the opt-in adaptive path with a focused regression first.

## 2026-11-05 — expanded panel-selection screens

- Replayed the diagnostic across seven archived cases: high-blue at three grid targets, high-blue `kappa10=100`, standard positive-cr, low-T positive-cr, extreme-red positive-cr, and high-`N_inf` blue. Each run matched its archived Fast grid and fixed-point background.
- Top-quarter curvature-ranked midpoints improved dense-Prüfer `DN_gw` discrepancies substantially in the high-blue target sweep and standard positive-cr/kappa corner (e.g. `7.56e-4→6.88e-5` at high-blue target `1e-5`; `1.10e-4→2.38e-5` standard positive-cr).
- It did not help uniformly: low-T positive-cr `2.04e-4→2.17e-4`, extreme-red positive-cr `6.89e-5→6.98e-5`, and high-`N_inf` only `1.15e-4→1.03e-4`. Embedded deltas were very small relative to the archived reference-grid changes on the first two cases, suggesting a different dominant error source or reference-grid behavior.
- No production change. Ruff, `py_compile`, and `git diff --check` pass. Commit `39c0343` is pushed; CI run `41973818558` is still in progress at the last query.
- Next: trace Fast-versus-Reference residuals and background/tail/propagation contributions for low-T and extreme-red, and prototype an online cumulative selector using only accessible coarse-grid signals. Require `DN_gw`, spectrum, cost, and held-out checks before production work.
# 2026-11-05 midday continuation

- Replaced the misleading nonnested-grid comparison with exact Prüfer evaluations at every Fast adaptive midpoint for eight fixed-background cases.
- Confirmed the previous low-T/extreme-red apparent discrepancies were interpolation artifacts; true nested integral changes are only `3.39e-6` and `1.36e-6`.
- Confirmed nonzero high-blue nested-grid changes across three target settings (`6.5255e-4`, `2.4525e-4`, `8.1061e-5`). Current curvature-ranked subset remains a research candidate only.
- Read four phase-function papers and recorded methods, assumptions, limitations, project mapping, citations, and local PDF hashes in `docs/research/phase_function_papers_batch1_20261105.md`.
- Production solver unchanged. Remaining: resolve selector cost and false-safe behavior, update phase status from complete evidence, then commit/push the verified research artifacts.

## 2026-11-05 small panel-budget continuation

- Screened curvature-ranked 5%, 10%, and 15% midpoint budgets against all eight exact nested Prüfer fixed-background maps. 10% passed the `2e-4 DN_gw` and `1e-3 dex` pointwise screens; 5% missed on two maps.
- Replayed 10% candidate grids through the public Fast outer solve with the archived coarse grid fixed. Warm P50 ratio was `1.014–1.125`, P95 ratio up to `1.129`: rejected against the `5%` speed budget.
- A coarse summed-curvature trigger would select every screened case, including two with nested integral changes below `4e-6`, so it is not a selective general trigger.
- Read-only Prüfer root secants on base and 10% grids passed `2e-4` at high-blue targets `1e-5`/`1e-6`; target `1e-5` root discrepancy worsened from `4.91e-5` to `7.80e-5`, target `1e-6` improved from `3.63e-5` to `9.28e-6`. No rigorous root enclosure.
- Scripts and artifacts: `screen_adaptive_panel_fraction_budget.py` / `adaptive_panel_fraction_budget_20261105.json`, `screen_adaptive_panel_candidate_roots.py` / `adaptive_panel_candidate_root_screen_20261105.json`. Production sources unchanged; fresh Ruff, compile, JSON validation, and diff checks remain before commit.

## 2026-11-05 literature/project audit continuation

- Read full Agocs–Barnett (2024) and Piron–Tacu (2025) algorithm sections and downloaded/validated both source PDFs; updated `docs/research/phase_function_papers_batch1_20261105.md` with equations, assumptions, failure modes, and prior-experiment deduplication.
- Source-audited `oscode` at `4dd643b`; it uses adaptive WKB/RK local steps, C++17/Eigen, and interpolated coefficient arrays. Local execution was unavailable because `pyoscode` is absent and no C++ compiler is installed; a pip build-dependency attempt was stopped, with no package or project source changed.
- Previous CI run `41979029318` passed. Remaining: commit/push this research record and continue the original broad literature/open-source and domain algorithm work.

## 2026-11-05 official `riccati` reproduction

- Read the upstream Python implementation and executed its official Airy and high-frequency burst tests: 2 passed, 8 invalid-escape warnings, 0.94 s.
- A separate analytic constant-frequency oscillator (`ω=40`, `[0,10]`) used one Riccati step and had zero reported endpoint-state and derivative relative error.
- This is only an ARDC implementation reproduction, not project-domain validation; next continue the remaining named papers and repositories and complete the current Phase 6c gates.

## 2026-11-05 WKB micro-window paper

- Read and archived the full Rosemeier–Klein 2026 preprint; the method uses short micro-simulations to fit slow modal amplitude/phase and extrapolate larger steps under a fixed skew-symmetric fast operator plus weak nonlinear remainder.
- The project’s time-varying frequency, horizon turn, and reheating kink do not yet meet that model contract. No code was imported or run; next derive whether a local constant-operator split with bounded residual exists before any prototype.

## 2026-11-05 semiclassical WKB paper

- Read the complete 2010 author manuscript for the required 2011 SIAM WKB-based scheme; logged its transformed system, second-order error bounds, positive-coefficient assumption, and phase-quadrature floor.
- A deep-subhorizon mapping might exist after removing damping, but it is not uniform at horizon entry and has not been shown to beat the exact phase transfer. Keep this as an unproven inference; no production edits.

## 2026-11-05 Neumann–Filon paper

- Read and recorded Perczyński–Madejski’s Neumann-series/Filon scheme and local error bounds. It treats a finite Fourier time forcing multiplying an elliptic spatial generator; this does not directly match the homogeneous, time-varying tensor-mode matrix equation.
- Existing Duhamel defect correction is already rejected by project accuracy/runtime evidence. Continue literature review without relabeling that method.

## 2026-11-05 Magnus review and deduplication

- Downloaded and validated the arXiv PDF for Blanes–Casas–Oteo–Ros, *The Magnus expansion and some of its applications*; read the Magnus theorem/convergence sections, integrator construction, matrix-exponential cost, commutator-free methods, and oscillatory Schrödinger application. Verified journal metadata and DOI against the publisher record.
- Recorded the sufficient `∫||A||<π` convergence condition, local-series limitation, structure preservation, step splitting, and exponential costs in `docs/research/phase_function_papers_batch1_20261105.md`.
- Cross-checked against project experiments: linear-z closed 2×2 Magnus took 1.8–2.0× baseline time; CF4 took 1.42–1.63× baseline time in the unfiltered kernel screen and showed no meaningful independent-oracle accuracy gain. No new Magnus or CF4 production candidate is justified without a substantive error-control/cost change.
- No production source changed. Six named first-batch papers and six of eight named open-source projects remain unaudited; broader domain algorithm and region gates remain open.

## 2026-11-05 reduced-basis literature screen

- Downloaded and checked Ballani–Kressner's complete author-hosted EPFL PDF; read the reduced-basis error estimator, affine operator representation, rbTensor algorithm, numerical setup/results, and stated limitations. Verified journal metadata and DOI with SIAM.
- Recorded that the proposed many-query reduced basis / low-rank tensor method depends on parameter separability, greedy training, and a posteriori residual bounds; the paper's experiments use finite parameter grids and note limitations in black-box reliability and continuous-domain control.
- Project fit remains a research hypothesis only: no certified conversion of nonnormal tensor-mode propagation, spectrum integration, and self-consistent-root errors follows from the paper's elliptic coercivity bounds. No code or production source changed; no reproduction was claimed.
- Five of the sixteen named first-batch papers and six of eight named open-source projects remain unaudited; broad Fast-domain and performance gates remain open.

## 2026-11-05 greedy reduced-basis theory screen

- Retrieved the complete author manuscript from the RWTH Aachen IGPM publication archive and cross-checked the final SIAM metadata. Read the weak-greedy definition, polynomial/exponential width-rate theorems, noisy-snapshot robustness result, and exact-greedy sharpness example.
- The main useful condition is explicit: the greedy method inherits polynomial decay only when the solution manifold already has matching Kolmogorov-width decay and the practical residual selection is uniformly `γ`-weak; bounded snapshot errors create an error floor. This does not certify the project’s spectrum or self-consistent root.
- Added the source hash and project mapping to `docs/research/phase_function_papers_batch1_20261105.md`. No code or production source changed; no numerical reproduction was attempted.
- Four of sixteen named first-batch papers and six of eight named open-source projects remain unaudited; continuous-region algorithm and accuracy/speed gates remain open.

## 2026-11-05 adjoint sparse-grid error literature screen

- Downloaded and parsed the full Sandia/OSTI author manuscript for Jakeman–Wildey and read the weighted-residual identity, sparse-grid sample estimator, dimension/local refinement indicators, estimator surrogate, cost discussion, and numerical examples. Cross-checked DOI and final JCP metadata.
- The estimator targets a chosen QoI of a parameterized PDE/ODE and can reduce redundant sparse-grid model solves, but requires forward/adjoint state surrogates and nontrivial residual work; its own error is approximate and tied to smoothness and forward/adjoint accuracy.
- Compared the proposal with the existing project scalar-adjoint screen: 65.1 ms versus 2.68 ms baseline (24.3×), nonfinite output, and no full-spectrum reconstruction. Thus no direct production candidate is justified without a fundamentally cheaper and validated discrete-adjoint path through frequency integration and the fixed-point root.
- Added source PDF/hash and applicability/deduplication notes. No production source changed and no new numerical reproduction was claimed. Three of sixteen named first-batch papers and six of eight named open-source projects remain unaudited; broad Fast-domain gates remain open.

## 2026-11-05 empirical-interpolation error-bound screen

- Downloaded and parsed the authors’ final manuscript for Eftang–Grepl–Patera; checked the DOI and version against MIT DSpace. Read the full EIM construction, uniform derivative bound, Proposition 2.1, proof, and scalar/two-parameter examples.
- Recorded the rigorous uniform error bound’s exact ingredients: derivative suprema over the parameter domain, finite-set covering radius, derivative interpolation errors, and the EIM Lebesgue constant. The bound is valid across the domain only if these inputs are actually available; numerical examples show a plateau as interpolation order grows.
- Candidate relevance: a possible certificate for an EIM approximation of Fast’s nonaffine mode generator, but coefficient error is not yet output/root error. Need compact smooth coordinates, analytical derivative bounds that respect parameter branches and reheating, and stability propagation into the oscillatory transfer before any implementation.
- Search of current docs/scripts/source found no existing EIM/DEIM prototype. No production code changed and no reproduction was attempted. Two of sixteen specified papers and six of eight named open-source projects remain unaudited.

## 2026-11-05 EIM transfer smoothness audit

- Traced `cr`, `T_re`, neutrino-temperature, thermal-history, reheating-kink, `sigma/H²`, `f_hor`, `Phi/S2`, and tensor-kernel dependencies in current source.
- Concrete fixed-coordinate counterexample on `cr=0`: at fixed `N=10`, sweeping `DN_re=9.999999, 10, 10.000001` moves `N_re_abs` across that coordinate and changes `sigma(10)` from `1.9999999990742416` to `1.9999999990742432` to `1.0`. The branch surface therefore creates a finite jump as a function of the input parameter, so the global smooth EIM premise fails in native `(N, μ)` coordinates. On the usual `cr=1` point, the same physical kink has `N_re_abs=31.041498487849346` and one-sided values `1.0` / `1.9999999990742432`.
- `cr` changes the actual parameter map (`n_t`, `N_inf`, and consistency-derived `DN_re`); `sigma/H²` also switch thermal formula regimes and spline pieces. A branch-aligned EIM with exact event matching is not ruled out, but its derivative bounds and coefficient-to-root error propagation remain unproved. The detailed evidence is in `findings.md`.
- No Fast or Reference code changed. Next EIM work should quantify whether smooth branch-local coordinates give useful derivative bounds/cost at all; do not treat the source-paper coefficient certificate as a Fast output certificate.

## 2026-11-05 normalized background-coordinate probe

- Derived an exact normalized coordinate for the `N >= N_re_abs` background branch: `s=(N-N_re_abs)/N_re`, so `N_inf-N=(1-s)N_re`. This removes the independent inflation-length offset from the `sigma_vec`/`H2_vec` expressions within that segment while preserving the event boundary explicitly.
- Replayed `sigma_vec` and `H2_vec` at five common `s` values for otherwise identical `cr=0` models with `DN_re=8` and `12`: same `N_re`, bitwise-identical `sigma` and `H²` arrays. The exact kink still causes a finite jump when represented in native fixed-`N` coordinates, so this suggests branch alignment rather than smoothing.
- This is background-level evidence only; it says nothing yet about `Phi/S2`, the transfer state, frequency support, `DN_gw`, self-consistent root or runtime. Next determine whether the normalized invariance survives the complete kernel and whether it reduces EIM derivative bounds/cost enough to matter.

## 2026-11-05 complete Fast normalized-coordinate probe

- [x] Added a reproducible diagnostic script and raw JSON output covering the aligned background, primitive shifts, frequency support, physical-guarded Fast roots, and eight common-frequency spectrum samples for `DN_re=8` and `12`.
- [x] Corrected resource-budget import order so Numba imports only after environment defaults are installed; reran both cases successfully with one Fast thread and captured the restored thread budget.
- [x] Recorded that post-event `sigma`, `H²`, and `f_hor` match, but full spectra differ by up to `0.006073 dex` at common frequencies and the upper support changes. Cross-parameter differences are physical sensitivity, not a numerical accuracy estimate.
- [x] `ruff check`, `py_compile`, `git diff --check`, and the standalone probe pass. No production source or Reference implementation changed.
- [ ] Commit/push the reviewed research artifact, then test whether a branch-aligned coordinate offers a useful general Fast mechanism; current evidence is only two parameter points and does not meet regional accuracy/speed gates.

## 2026-11-15 underflow mode-state attribution

- Root cause of the earlier large fixed-background spectrum difference: three Fast experiment scripts seeded `DN_eff` after the initial expansion. Corrected the initialization order, regenerated the 235/469/937 JSON, and refreshed the dependent tail, step, handoff, PCHIP, and state diagnostics. Fast/Reference production sources remain unchanged.
- With matched backgrounds, z5 max errors are `0.000480`, `0.000648`, and `0.002385 dex` at 235/469/937. Z8 reduces them to `0.000390`, `0.000390`, and `0.000468 dex`. The 937 z5 outlier is insensitive to h refinement (`0.002385→0.002377 dex` from h=.005 to .000625), while tail z8 reduces that mode to `3.54e-6 dex`.
- Added direct 15-node local frequency refinement around the 937 outlier. At z8, max/P95 Fast-to-Reference z10 error is `1.44e-5/1.22e-5 dex`; the original 937-node Reference PCHIP misses these nodes by up to `0.00882 dex`.
- Replayed the formal Fast settings and reconstructed the z8 `assemble_tail` `Opgw` result from the saved state; prediction/production ratio is `0.99999999999994`. This validates one analytic-tail path; phase-resolved blending remains unproven.
- Completed direct Reference z10 on all 235/469/937 nodes. Reference z8→z10 max spectral shift is `2.069e-4 dex` on each grid (P95 about `1.42e-4`); Fast z8 vs Reference z10 max errors are `0.000353`, `0.000380`, and `0.000391 dex`.
- Next: run alternating repeated z5/z8 cost measurements at the same fixed background/root/grid/thread, then test one independent non-underflow parameter point before considering any profile change.
## 2026-11-15 tail-depth cost study

- Ran 12 alternating z5/z8 pairs at the same frozen root `0.027412228089632874`, native 937-node frequency grid, Fast profile, and one thread. Warm P50/P95 were z5 `0.1299/0.1351 s` and z8 `0.2034/0.2128 s`, a 1.566 P50 ratio. Both replay spectra were bitwise repeatable (`max spread=0`); peak RSS was about 173.7 MB.
- First z5 call took `0.446 s`; the first z8 call was `0.208 s` after Numba had already been warmed by z5, so this is not an independent z8 cold-start measurement. Fixed-map residuals were `6.31e-8` (z5) and `5.38e-8` (z8).
- z8 reduces the observed tail-sensitive maximum spectrum error below the Reference z10 threshold on these three grids but costs about 56.6% more steady runtime. Keep it validation-only pending a second parameter regime and profile-level performance evidence.
## 2026-11-15 non-underflow tail-depth continuation

- Tested a normal-amplitude `cr=0`, `r=.01`, `n_t=0` case. Fast self-consistent root converged, and fixed-root z8 versus Reference z10 passed both the `1e-3 dex` spectrum and `2e-4` root residual gates at this one grid.
- The Reference z8 fixed map missed the root residual threshold (`2.94e-4` relative); Reference z10 residual was `4.02e-6`, confirming the need to check root/tail convergence separately from pointwise spectra.
- An `n_t=.5` attempted sample was rejected by `shared_Neff_guard` before spectral comparison.
- Next: evaluate the `cr>0` consistency branch with physical amplitude and direct Reference; then determine whether the increased z8 cost is acceptable across more than two `cr=0` points.

## 2026-11-15 `cr>0` consistency-branch continuation

- Added a physical-amplitude Reference comparison at `cr=1`; no unit-amplitude scaling bypassed the amplitude-dependent background. Fast and Reference derived the same nonzero `A_t` exactly.
- Fast root closed; z5 max spectrum error vs Reference z10 exceeded `1e-3 dex`, while z8 passed on the 76-node grid. Reference z8 map missed the root threshold but z10 closed to `4.96e-6` relative.
- Next: screen the `cr>0` amplitude-underflow boundary and record whether rejection is correctly classified; do not infer it from the successful normal-amplitude point.

## 2026-11-15 consistency-background underflow repair

- Added a failing regression test for nonzero `r,A_s` whose product underflows in the `cr>0` consistency background; observed the original `math domain error` before implementation.
- Rewrote only `V_inf` and `Delta_N` evaluation in log space when the stored `A_t` is zero but both inputs are positive. The physical cutoff remains enforced: this sample computes `Delta_N=-209.962` and Fast rejects it as `invalid_cutoff`; no spectrum is reported.
- Preserved distinct exact-zero behavior: `cr>0,r=0` returns `invalid_r`, while exact `A_s=0` returns `invalid_cutoff`; the `cr=0` log-power path still solves the nonzero underflow sample.
- Three targeted consistency tests pass, followed by the complete suite/static/Mypy/build gates. Normal-amplitude `cr>0` z5/z8 outputs were rechecked unchanged.
- 2026-11-15 Round 34 follow-up: the 3×3 `cr>0`, `r=1e-316` matrix across `T_re={10,2000,10000}` and `kappa10={1e-3,1e-2,1}` yielded finite log-domain background values but `Delta_N<0` at all nine discrete points, all rejected before primitive construction as `invalid_cutoff`. Re-running the existing ordinary-amplitude stress control reached `S2` and reproduced `nonfinite_background_primitive` at index 13763. Script/JSON are `scripts/scan_consistency_underflow_boundary_20261115.py` and `docs/research/results/consistency_underflow_boundary_20261115.json`; this separates cutoff rejection from `S2` overflow without broadening the `cr>0` underflow domain.
- 2026-11-15 Round 35: ran a six-value `cr=0`, `n_t=17.38` amplitude ladder through the saved native 235-node grid. Fixed-background scaled spectra and integral maps agree to near float precision; all five positive underflow amplitudes remain distinct from exact zero and return finite log spectra. Self-consistent root/map scale with amplitude to ~`1.2e-11`, while normalized self-consistent spectra shift ~`1.92e-3 dex` from background changes. Raw JSON and script are in `underflow_amplitude_ladder_20261115.*`; independent Reference and 469/937 convergence remain required.
- Round 35 Reference follow-up: one direct 235-node unit-amplitude Cartesian Reference at the common frozen root, analytically log-scaled across the five nonzero amplitudes, gives max Fast-vs-Reference spectrum error `0.000480448 dex` and map relative error `8.35e-7` for each. This confirms the fixed-background underflow scaling against the independent pipeline at this discrete point; self-consistent roots other than the base still lack Reference/grid-convergence checks.
- 2026-11-15 Round 36: separated Fast frequency-discretization effects from root movement over three nonzero underflow amplitudes. 235→469 roots shift `−2.2160e-6`; 235→937 shift `−2.3031e-6`, with fixed-root maps consistent across amplitudes. Independent unit Reference at each 235 root yields max spectral differences below `4.81e-4 dex`. The Reference PCHIP warnings persist (4 per sample); higher-grid Reference and z8→z10 checks for the new roots remain open. See the round 36 JSON/script and research note.
- 2026-11-15 Round 37: swept root tolerances and independently remapped each of three roots. Roots stayed identical across `1e-5..1e-11`; independent residuals were at roundoff. Background evaluation counts varied for the minimum subnormal at `tol=1e-11`, confirming a tighter stopping decision without a root shift. This remains scoped to the one 235-node parameter slice.
- 2026-11-15 Round 38: retained original Reference and compared its 235-node underflow integral with stable positive linear PCHIP arithmetic. Warnings reproduced (four per call), stable route emitted none, maps agreed within `7.11e-15` relative; recorded 13 zero and 2 subnormal normalized nodes. Corrected and reran the amplitude-ladder JSON after finding its Reference payload had been overwritten by a local baseline variable. The stable route remains diagnostic-only; no Reference production change.
- 2026-11-15 Round 39: minimum positive `r=5e-324` direct grid study found a z5 937-point max spectral failure (`0.002441 dex`). Fast z8 vs direct Reference z10 passed (`0.0003722 max`, `0.0000233 P95`); Reference z8→z10 changed `0.0001649 max`, `0.0001421 P95`. Reference z10 cost 1673 s and kept the known PCHIP warnings. z8 remains validation-only; no default change.
- Round 39 mode attribution: the failing `log10f=-16.291944809` mode has Fast z5 handoff epsilon `2.496e-3`; with z8, available evolution ends at `z=4.99072<8`, so Fast and Reference both propagate numerically to today. Point error contracts from `−2.441e-3` to `−4.34e-6 dex`. This identifies early WKB switching as the mechanism for this point.
- 2026-11-15 Round 40: `z_today<5` delayed-handoff mask passed same-root spectral/map screens at three discrete cr0 points, including an independent tilt plus low-T/low-kappa 221-node sample. It selects ~5% of modes but has no real per-mode runtime measurement. High-T/high-kappa exploratory root crossed physical guard and was excluded. Candidate remains diagnostic-only.
- 2026-11-15 Round 41: extended the mask screen to thresholds 5–8 and one normal-amplitude `cr=1` point with direct Reference z10. Threshold 5 fails that point; threshold 6 passes all five discrete cr0/cr1 screens (`max ≤5.78e-4 dex`, map/root residuals ≤`1.02e-5`) while selecting 5.5–9.2% of rows for z8. No production/kernel changes or speed claims.
- 2026-11-15 Round 42: added a private per-mode kernel handoff path, disabled unless the model explicitly carries experimental attributes. The first full run exposed four default-path regressions because the kernel call was mistakenly nested under `tail_state`; moved it outside the state-reset guard. Targeted and full regression then pass. Full suite: `408 passed, 9 deselected`; Ruff, configured Mypy and wheel/sdist build pass. Same-root/grid one-thread 12-pair A/B at 937 nodes gives z5/per-mode/z8 P50 `0.2270/0.2283/0.3079s`, P95 `0.2471/0.2341/0.3563s`; per-mode max error is `4.855e-4 dex` against direct Reference z10, with zero repeated-spectrum spread. Other holdouts have single-run accuracy checks but not Reference z10 convergence. Experimental path remains private; defaults unchanged.
- Round 42 Reference follow-up: the low-T/high-tilt 221-node holdout now has direct Reference z8/z10 tails; spectral shift is `1.648e-4 dex` max / `1.598e-4 dex` P95. Actual per-mode vs z10 is `5.692e-4 dex` max, mapped root residual `−2.22e-16`. Two failed script runs were diagnosed as serialization/list-arithmetic errors after completing the expensive solve; checkpoints and NumPy-to-list conversion now preserve the computed evidence. One more z10 background closes, but `n_t=17.30` and higher-grid holdouts remain.
- Round 42 `n_t=17.30` Reference follow-up: z8→z10 shifts `1.550e-4 dex` max / `1.418e-4 dex` P95; actual per-mode vs z10 max/P95 is `4.746e-4/3.119e-5 dex`, map error `−4.48e-6`, fixed-root residual 0. The per-mode candidate now has four direct Reference z10 comparison points; changed-root frequency-grid and multi-resource evidence still required.
- Round 42 resource follow-up: same-grid/root 12-triplet A/B at 1/2/8 threads shows per-mode P50 within 1% of z5 and `1.16–1.55×` faster than z8; max z10 spectrum error stays `4.855e-4 dex`, repeats are bitwise stable. Multi-resource gate for this single sample is complete; parameter and grid generalization remains.
- 2026-11-15 Round 43 start: confirmed detached worktree HEAD equals fetched `fast/fast_v0.2` (`39b23179c411098fc393fe7a3bf3d24ac3b7a850`); tracked tree clean. New boundary screen shows finite `log_A_t` when `r=5e-324` stores `A_t=0`, while exact `r=0` remains a separate zero-limit success. Began a unit-transfer Reference log-integral falsifier on the exact Fast-generated 75-node grid.
- Round 43 diagnostic attempt 1 completed the expensive Reference computation but JSON serialization rejected `-inf`. Backtrace plus data flow identified the field: exact-zero `fast_case` serialized `_log_tensor_amplitude == -inf`, not an invalid integral. Updated only the diagnostic serializer to record non-finite `log_A_t` as JSON `null`; no Reference or production code changed. Re-run needed to persist the already-computed-but-lost measurement.
- 2026-11-16: Verified the real collection failure at base `8dfbee5`: its public-history test imported through the isolated research validator, which depended on the missing `prototype_cr0_event_chart_kernel`. The already-committed `d6ddbb2` moves test helpers to `scripts/cr0_event_chart_public_history.py`; the optional prototype import remains only in the validator. Current GitHub CI for `39b2317` is successful and local full suite passes, so this priority needs no further code edit.
- Round 43 implementation QA: full `pytest -q` passed (`408 passed, 9 deselected`, 2 deprecation warnings); `mypy` and `python -m build` passed. Initial Ruff run caught one import-order issue in the new z12 script; after adding the required import separation, `ruff check .` passed.
- Round 43 tail/fixed-point completion: 457-node direct Reference z12 completed in `5936.14 s`; z10→z12 map change `−6.136e-6`, interpolation estimate `1.92808e-4`, quadrature estimate `1.11e-14`, no warnings. Direct Fast fixed-root reruns on 235 and 457 nodes each had zero map residual and identical grid hashes; `Nv`, `sigma`, `f_hor` match the zero-DN background bitwise. The 457-node Fast root differs from Reference z12 by `−1.168e-6` relative. This closes the discrete `cr=0,n_t=0,r=1e-315` root/tail/frequency check only.
- Round 43 pointwise accuracy completion: direct Reference z10 vs Fast z5 on the same 457 nodes took `800.73 s`; maximum spectrum error `4.67178e-4 dex`, P95 `3.07251e-5 dex`, worst `log10(f/Hz)=-15.84708`; mapped-root difference `−7.303e-6`. Full nodewise spectra are in `subnormal_spectrum_reference_z10_20261116.json`; no broader parameter-domain or performance claim.
- Round 43 final code verification after adding the stable `log_pchip` sub-branch: full `pytest -q` passed (`409 passed, 9 deselected`, two existing deprecation warnings); Ruff, configured Mypy, `python -m build`, and `git diff --check` passed. All new and historical JSON outputs contain the exact source hashes used; the first baseline/boundary and 235/457 z8/z10/z12 outputs intentionally retain their earlier Fast source SHA, while the independent fixed-root and pointwise-spectrum artifacts use the final source SHA.
- Round 43 delivery: commit `65dbd55e3207394d1391e9e5217ea65179ac1f1c` was authored/committed with `2966684515@qq.com`, pushed to `fast_v0.2`, and remote SHA verified by fetch. GitHub Actions run `43738808345` completed successfully for Canonical, static, package, Cobaya, and Python 3.9–3.13 compatibility jobs.
- 2026-11-16 Round 44 independent `n_t=17.38,r=5e-324` holdout: verified current Reference SHA and identical nested grid hashes, reran current Fast self-consistent 235/469/937 grids plus direct fixed-root maps, and compared all 937 Fast z8 spectrum nodes to the archived direct Reference z10 spectrum at the same frozen root. Map error `−3.2534e-6`, spectrum max/P95 `3.7216e-4/2.3253e-5 dex`, roots shift `−2.2160e-6/−2.3031e-6`, fixed-root residual 0. The raw JSON retains the Reference warnings and full spectra.
- 2026-11-16 Round 45 low-T/high-tilt holdout: with the source/hash-verified direct Reference z10 on the exact 221-node grid, current Fast `r=minsub,n_t=17.45,T_re=10,kappa10=.001,cr=0` root differs from frozen Reference root by `−2.44e-15`; own-root remap residual 0, Fast z8 vs Reference z10 max/P95 `5.6918e-4/2.0551e-5 dex`; Reference z8→z10 max/P95 `1.6483e-4/1.5985e-4`, map shift `−1.2447e-5`. Artifacts/runner added; no production source changed.
- Round44/45 delivery: commit `90e968c9ab13538c019d8e2c9b771a98ec71c222` was pushed and fetched at `fast_v0.2`; GitHub Actions run `43740889847` completed successfully for all workflow jobs.
- 2026-11-16 Round46 `kappa10` decade holdout: current Fast at `n_t=17.45,T_re=10,kappa10=.01,r=minsub,cr=0` used a direct 227-node grid and completed self-consistency. Direct Reference z8/z10 on that exact grid yielded Reference tail max/P95 `1.6498e-4/1.5849e-4 dex`, z8→z10 map `+1.2049e-5`; current Fast z8 vs Ref z10 max/P95 `3.1394e-4/2.2217e-5 dex`, map error `−3.294e-6`, fixed-root residual 0. Artifacts include Reference checkpoints, full spectra, input hashes, and current source hashes; no production code changed.
- 2026-11-16 Round46 kappa-decade holdout: current Fast at `n_t=17.45,T_re=10,kappa10=.01,r=minsub,cr=0` used a direct 227-node grid and completed self-consistency. Direct Reference z8/z10 on that exact grid yielded Reference tail max/P95 `1.6498e-4/1.5849e-4 dex`, z8→z10 map `+1.2049e-5`; current Fast z8 vs Ref z10 max/P95 `3.1394e-4/2.2217e-5 dex`, map error `−3.294e-6`, fixed-root residual 0. Artifacts include intermediate checkpoints, source hashes and full spectra; no production code changed.
- Round46 delivery: commit `285ddab02c1c2d07f1cd6e035c99a9468e04c832` pushed/fetched at `fast_v0.2`; GitHub Actions run `43743231707` completed successfully for all jobs.
- Round47 consistency-branch isolation: current `cr=1` normal-amplitude Fast z8 agrees with same-grid direct Reference z10 at `3.5795e-4 dex` max; physical `A_t` is normal and the integral-first helper was not called. Current-source re-run of the 3×3 `r=1e-316,cr=1` background matrix again yielded all `Delta_N<0`, `N_inf=None`, pre-primitive `invalid_cutoff`; the separate S2 pressure control remained `primitive:S2`. Targeted tests for nonzero underflow/exact-zero `cr>0` cases passed.
- 2026-11-16 Round48 `cr>0` valid-cutoff sample: analytically selected `r=1.1*r_cut=3.962216254776857e-43` at `T_re=2000,kappa10=.01`; initial `Delta_N=.03177006>0`. Current Fast self-consistent root is `1.067074488428421e-43`; direct Reference z8/z10/z12 used the exact new 77-node Fast grid and same root. Fast fixed-root map vs Reference z12 differs `−6.086e-6`; Fast z8 vs Reference z12 raw spectrum max/P95 is `3.776e-4/2.567e-4 dex`; Reference z10→z12 map/max-spectrum shifts are `−1.953e-5/3.647e-4 dex`. Reference public log spectrum floors all 77 nodes at `−40`, so comparisons use raw positive finite `Ogw-Oj`; no warnings. z12 required 1290.7 s. No production code change or performance claim.
- Round48 boundary and valid-point scripts/JSON plus the research ledger were committed as `797c82a2530f89d9eab0b8f3b60653add932407a`; required Author/Committer email was verified, push/fetch matched the remote SHA, and GitHub Actions run `43748405989` passed all jobs. Next: independent changed-`T_re`/`kappa10` cutoff-valid `cr>0` holdout.
- 2026-11-16 Round49 low-temperature `cr>0` cutoff-valid holdout: exact analytic boundary from the hashed 3×3 matrix, `r=1.1*r_cut` at `T_re=10,kappa10=.001`; direct Fast/reference z8/z10/z12 on the same new 76-point grid. Fast root `3.287879454649206e-64`; fixed Fast map vs Reference z12 `+2.571e-6`; Fast z8 vs Reference z12 raw log-spectrum error max/P95 `5.693e-4/4.881e-4 dex`; Reference z10→z12 map shift `+3.389e-6`, max spectrum shift `5.410e-4 dex`. All 76 public Reference log nodes hit the `−40` floor; raw signal was finite/positive, with zero warnings. Runner is now parameterized; this remains one discrete low-temperature check. Next test the high-temperature/high-kappa boundary corner.
- 2026-11-16 Round50 high/high `cr>0` cutoff corner: current Fast at `T_re=10000,kappa10=1,r=1.1*r_cut`, direct Reference z8/z10/z12 on exact 77-point grid. Fast root vs z12 Reference `−6.074e-6`, fixed-map comparison `−5.402e-6`; spectrum max/P95 `8.357e-5/2.859e-5 dex`. Reference z10→z12 map and maximum spectrum changes `−3.942e-5` and `8.712e-5 dex`; z8→z10 map was `−2.528e-4` and is explicitly retained as above the root screen. 59 nodes hit Reference's public `−40` floor; raw positive finite `Ogw-Oj` used, no warnings. z12 took 1198 s. Low/mid/high diagonal points now exist but T and kappa co-vary; next hold kappa at `.01` and vary T.
- 2026-11-16 Round51 fixed-`kappa10=.01` temperature control at `T_re=10`: direct z8/z10/z12 Reference on the exact 76-node Fast grid. Root vs Reference z12 `+6.326e-6`; fixed-map difference `−1.235e-6`; Fast z8 vs Reference z12 max/P95 `3.752e-4/2.980e-4 dex`. Reference z10→z12 map/spectrum change `+4.972e-6/3.741e-4 dex`, no warnings. All public Reference log-spectrum nodes are floored; comparisons use positive finite raw `Ogw-Oj`. z12 took 1467 s. Next use `T_re=10000,kappa10=.01` to complete a temperature sweep at fixed kappa.
- 2026-11-16 Round52 fixed-`kappa10=.01` high-temperature control: `T_re=10000`, exact 77-point Fast grid, direct Reference z8/z10/z12. Fast root vs Reference z12 `−9.856e-6`, fixed-map difference `−9.103e-6`; max/P95 spectrum error `3.845e-4/2.394e-4 dex`. z10→z12 map/spectrum `−3.432e-5/3.549e-4 dex`; z8→z10 map `−2.542e-4` is logged above the screen, not hidden. All public Reference logs floor, raw positive finite signal used, no warnings; z12 took 1243 s. Fixed-kappa temperature samples at 10/2000/10000 are now available. Next preflight blue tilt on valid `cr>0` cutoff side with physical guard.
- 2026-11-16 Round53: corrected a branch-semantic mistake: `cr>0` derives `nt` from the consistency relation and ignores the supplied `n_t`; the `17.3/17.45` inputs both yield `nt=−4.952770318471071e-44`, so those are duplicate physical simulations, not high-tilt holdouts. A separate `cr=0` minimum-subnormal sample retains true `nt=17.45`.
- Round53 diagnosed the cutoff-valid `cr>0` z5 spectrum miss at `log10f=−16.0917663`: WKB handoff epsilon `3.3378e-3`, while fixed-point root error vs Reference z12 is only `6.68e-6`. Solved the existing endpoint ripple bound for `1e-3 dex`, obtaining `z_today=6.0748686`; deferring only eligible modes to `z_tail+3` removes that local error.
- Round53 same-root/grid/thread current-source A/B: near-cutoff cr1 z5/adaptive/z8 P50/P95 `4.00/5.51`, `4.45/6.41`, `13.05/16.12 ms`; max spectrum errors `1.886e-3/5.639e-4/3.776e-4 dex`. A normal-amplitude cr1 control z5/adaptive max errors `1.429e-3/3.741e-4 dex`. The cr0 subnormal phase-active control remains `5.695e-4 dex` with zero delayed modes. A/Bs are discrete points, not broad speed claims.
- Round53 production implementation uses the ripple-derived z_today threshold by default, preserves phase-resolved active modes, and avoids per-mode dispatch when no mode is selected. Its original note preceded delivery; commit `5a612e1` was pushed and CI `43815188536` passed, then CI-result ledger commit `ba8c6b4` was pushed and CI `43815931166` passed.
- Round53 937-node independent cr0 minimum-subnormal holdout completed on the exact nested direct grid (no interpolation): Fast z5/adaptive/z8 max errors vs direct Reference z10 are `2.440865e-3/4.855091e-4/3.721593e-4 dex`; root error `−2.3031e-6`; Reference z8→z10 map/spectrum `3.6809e-6/1.6489e-4 dex`. The 44 selected modes include a nonmonotone-background mode with `z_today=4.99072<z_tail`; the final conservative mask therefore delays every phase-inactive mode with `z_today<z_cut`. Adaptive P50/P95 `126.22/132.32 ms`, z5 `125.39/127.89 ms`, z8 `199.15/219.92 ms`; gain is versus z8 only. Reference still reports four PCHIP divide-overflow warnings and 49 normalized zeros. Raw result and runner: `cr0_underflow_937_handoff_20261116.json` and `diagnose_cr0_underflow_937_handoff_20261116.py`. Latest-mask full suite passed (`411 passed, 9 deselected`); Ruff, Mypy, sdist/wheel build, and `git diff --check` passed. Commit `5a612e1dbf2a3de1f703da753764a2229d46e3b4` was authored and committed with the required email, pushed and fetch-verified at `fast_v0.2`; GitHub CI run `43815188536` passed static, package, Cobaya, canonical, and Python 3.9–3.13 compatibility jobs.
- 2026-11-16 Reference PCHIP 937 diagnostic: exact canonical direct Reference z10 nodes yield four divide-overflow warnings and 49 zero / 10 subnormal nodes after global scaling over `362.22 dex`. Log arithmetic of the same positive linear PCHIP gives the same total integral to `1.42e-14` relative with no warnings; zero-adjacent panels sum to only `10^-322.716` of total integral. Midpoint comparison shows 48 original PCHIP values underflow to zero; among 889 representable outputs P95 difference is `2.84e-14 dex`, max `0.3009 dex` at `log10(f/Hz)=-15.51694` in the weak tail. This explains the stable map alongside local interpolation risk; Reference source remains untouched. Diagnostic script/JSON hash checks passed, `pytest tests/test_stable_log_pchip.py` passed (3 tests), Ruff passed. No general replacement claim.
- 2026-11-16 expanded PCHIP direct-data check: current Reference SHA and exact frequency hashes verified for same-background 469 z8, 937 z8/z10, and a different low-temperature 221 z10 background. At the three 2000 K grid/tail variants, global scaling gives 25/5, 49/10, 49/10 zero/positive-subnormal nodes and four warnings; stable same-linear-PCHIP integral differences are `0`, `1.42e-14`, `1.42e-14`, with zero-adjacent panel fractions `10^-321.72`, `10^-322.72`, `10^-322.72`. The distinct low-T case has 0/0 nodes, no warnings, and bit-identical integral. Weak-tail pointwise discrepancies persist up to `0.3009 dex`; no Reference production change. Multi-grid runner generated the raw JSON and passed Ruff.
- 2026-11-16 current-source adaptive WKB at direct 469 nodes: same frozen root/grid/thread z5/adaptive/z8 max spectrum errors `5.5744e-4/4.7462e-4/3.7591e-4 dex`; all pass. Adaptive selected 22 nodes and is `1.33%` slower than z5 at P50 but `1.58x` faster than z8. Current selfconsistent root vs Reference frozen root `−2.2160e-6`, own map residual zero, root vs Reference z8 map `+2.3209e-7`. Reference retains four PCHIP warnings/25 normalized zeros. Runner completed and Ruff passed; artifacts saved. This point already passes under z5 and does not demonstrate broad speedup.
- 2026-11-16 current-source `cr>0` adaptive WKB high-T holdout: `T_re=10000,kappa10=.01,r=1.1*r_cut`, exact 77-node direct Reference z12. Supplied `n_t=0` maps to actual derived `nt=−7.7617e-40`; cutoff margin `.03177006`, `A_t=1.3043e-47` normal. Reference z8→z10 root shift `−2.5422e-4` retained; z10→z12 `−3.4323e-5`, spectrum max `3.5486e-4 dex`. z5/adaptive/z8 spectrum max `1.3817e-3/5.4243e-4/3.8447e-4 dex`; adaptive selects 9 modes, final P50/P95 `4.50/6.15 ms` vs z5 `4.15/5.64`, z8 `12.69/20.15`; root vs Reference z12 `−9.8558e-6`. Commit `a86b1ed` pushed and CI `43819844787` passed all jobs.
- Matched low-T cr>0 cutoff handoff at fixed `kappa10=.01,T_re=10`: 76-node direct Ref z12; actual derived nt `−4.6891e-58`, normal `A_t=7.8797e-66`, valid cutoff margin `.03177006`. Ref z10→z12 map/spectrum changes `+4.9719e-6/3.7411e-4 dex`; root error `+5.8515e-6`. z5/adaptive/z8 max errors `2.5217e-3/4.3519e-4/3.7522e-4 dex`, adaptive delays `12/76`; P50/P95 `4.69/6.11 ms` vs z5 `4.14/5.93`, z8 `13.99/22.31`. Artifacts and reusable input/output runner saved; Ruff and both direct runs passed; commit/CI pending.
- Matched low-T cutoff-boundary kappa decade completed: at `T_re=10,kappa10=.001,r=1.1*r_cut`, z5/adaptive/z8 max errors `1.9847e-3/5.6956e-4/5.6931e-4 dex`; Reference z10→z12 map/spectrum `+3.3891e-6/5.4096e-4 dex`, Fast root error `+7.4970e-7`. Adaptive delays `13/76`, P50/P95 `4.14/4.25 ms` vs z5 `3.75/4.00`, z8 `14.12/14.87`; root and spectrum gates pass. The r value follows each kappa's cutoff, so this is a decade-spaced boundary pair, not fixed-r isolation. Script, artifact, and updated research ledger are ready for commit/CI.
- High-T cr>0 cutoff-edge `kappa10=1` completed at fixed `T_re=10000`: z5/adaptive/z8 max spectrum errors `1.7319e-3/7.0358e-4/8.3570e-5 dex`; adaptive delays 8/77 and root vs Reference z12 `−6.0744e-6`. Reference z8→z10 map `−2.5281e-4` retained; z10→z12 root/spectrum `−3.9420e-5/8.7119e-5 dex`, no warnings. Adaptive P50/P95 `3.99/4.60 ms`, 6.4% slower than z5 and `2.94x` faster than z8. This is an r-follows-cutoff sample, not fixed-r kappa isolation. Artifact saved; script and Ruff pending final commit/CI.
- 2026-11-16 Sobol cutoff-domain screen complete: scrambled 8-point design, fixed seed 20261116 across logT `[1,4]` and log-kappa `[-3,0]`, each with r=`1.1*r_cut`. All eight Fast samples are cutoff-valid, normal `A_t`, guarded/admitted, finite, warning-free, and map-closed. Direct z8/z10/z12 Reference completed for indices 0 and 1 at exact current Fast roots/grids: root errors `−5.246e-6/−7.590e-6`; default adaptive max spectrum errors `5.600e-4/8.133e-4 dex`; Reference z10→z12 root/spectrum `−3.041e-5/2.022e-4` and `+3.847e-6/2.417e-4`. Both pass; z8→z10 shifts above `2e-4` retained. Cutoff-probe `r=minsub` warning classified by separate audit; safe `.01` probe reproduces exact cutoff with bitwise-identical background components. Corrected default-adaptive replay distinguishes legacy label. Six Sobol points remain Fast-only; no continuous claim. Files saved; Ruff/source hash checks passed; commit/CI pending.
- Sobol index 2 direct Reference completed at exact preflight grid/root: low T/high kappa (`16.7307,.75557`), 76 positive finite nodes; Reference z8/z10/z12 no warnings. Adaptive Fast spectrum max/P95 `4.5093e-4/6.4026e-5 dex`, root error `−6.3226e-6`; z10→z12 root/spectrum `−4.3112e-6/1.2089e-4 dex`. z12 runtime `1333.47 s`. Metadata index label corrected with preserved calculation snapshot hash; 3/8 direct Reference holds completed. Source hashes and Ruff pass; commit/CI pending.
- Sobol index 3, the high-temperature/low-kappa corner `(T_re,kappa10)=(2623.915,.001642)`, was solved directly in Reference at the exact current Fast root and 77-node grid. Fast replay matched preflight exactly; all z8/z10/z12 Reference spectra had 77 positive finite nodes and no warnings. Default-adaptive Fast vs z12 max/P95 spectrum error is `7.202e-4/5.089e-4 dex`, and root error is `−8.337e-6`. Keep z8→z10 map shift `−2.990e-4`; z10→z12 map/spectrum changes are `+6.185e-6/4.854e-4 dex`; z12 takes `1297.17 s`. Four of eight fixed Sobol points now have direct Reference checks; four remain Fast-only. Artifact `cr_positive_sobol_index3_reference_20261116.json`; calculation snapshot and corrected replay runner are separately hashed.
- Sobol index 4, `(T_re,kappa10)=(4843.8098,.329228)`, was directly solved with Reference z8/z10/z12 at the exact current Fast root and 77-node grid. All Reference raw spectra had 77 positive finite nodes and no warnings (public logs all floor at −40, so comparisons use raw positive `Ogw-Oj`). Default-adaptive Fast vs z12 max/P95 error `2.140e-4/8.723e-5 dex`; root difference `−7.848e-6`. Retain z8→z10 map shift `−2.530e-4`; z10→z12 map/spectrum changes `−3.831e-5/1.414e-4 dex`; z12 took `1315.55 s`. This is fixed-Fast-root Reference mapping evidence, not an independently closed Reference self-consistent root. Five of eight fixed Sobol points now have direct Reference checks; three remain Fast-only. Artifact `cr_positive_sobol_index4_reference_20261116.json`.
- Sobol index 5, the mid-low-temperature/low-kappa point `(T_re,kappa10)=(45.75163,.00375797)`, was directly solved at the exact Fast root and 76-node grid. Reference z8/z10/z12 have 76 positive finite raw nodes and no warnings; all public logs floor at −40, so comparison uses raw positive `Ogw-Oj`. Default-adaptive Fast vs z12 max/P95 spectrum error `4.731e-4/3.881e-4 dex`; Fast root difference `−3.186e-6`. z8→z10 map shift `−1.380e-5`; z10→z12 map/spectrum changes `−4.979e-6/4.788e-4 dex`; z12 takes `1582.18 s`. This is fixed-Fast-root mapping evidence, not a separately converged Reference self-consistent root. Six of eight fixed Sobol points now have direct Reference checks; two remain Fast-only. Artifact `cr_positive_sobol_index5_reference_20261116.json`.
- Sobol index 6 `(T_re,kappa10)=(68.27371,.1635397)` was directly solved at the exact Fast root and 77-node grid. Reference z8/z10/z12 have 77 positive finite raw nodes and no warnings; public logs all floor at −40, so comparison uses raw positive `Ogw-Oj`. Default-adaptive Fast vs z12 max/P95 spectrum error `2.654e-4/1.438e-4 dex`; Fast root difference `−8.215e-6`. Preserve z8→z10 map shift `−2.745e-4`; z10→z12 map/spectrum changes `+2.795e-6/1.386e-4 dex`; z12 takes `1402.81 s`. This is fixed-Fast-root mapping evidence, not a separately closed Reference self-consistent root. Six of eight points now have direct Reference checks; only index 7 remains Fast-only after this result. Artifact `cr_positive_sobol_index6_reference_20261116.json`.
- Sobol index 7 `(T_re,kappa10)=(642.6992,.0116999)` completed direct Reference z8/z10/z12 at the exact Fast root and 77-node grid. All raw signals are positive finite, no warnings; public logs all floor at −40, so the raw positive `Ogw-Oj` logs are used. Default-adaptive Fast vs z12 max/P95 error `4.984e-4/2.940e-4 dex`; root difference `−7.002e-6`. Preserve z8→z10 map shift `−3.111e-4`; z10→z12 map/spectrum changes `+1.577e-5/3.433e-4 dex`; z12 took `1406.07 s`. This is the eighth direct Reference point in the fixed 8-point design, not continuous-domain certification or a separately closed Reference self-consistent root. Artifact `cr_positive_sobol_index7_reference_20261116.json`.
- Reference root-closure audit: `run_reference(dn_eff=...)` is one fixed map; its built-in self-consistent stop condition measures relative change against total `N_eff`, which is not itself a relative residual for tiny `DN_gw`. At Sobol index 7 on the exact 77-node z12 grid, direct Reference mapping from the Fast root changes it by `+7.0016e-6`; a second direct map at that updated root returns the same `1.655519043161566e-47` exactly, giving residual 0 and zero pointwise raw-spectrum change across 77 nodes, no warnings. z12 quadrature/interpolation error estimates `2.28e-66/2.04e-55`. This closes one point under an explicit `|F(x)/x-1|<=2e-4` criterion, not all eight. Artifact/script `cr_positive_sobol_index7_reference_root_closure_20261116.json` / `verify_sobol_index7_reference_root_closure_20261116.py`.
- Independent Reference root-closure repeat at Sobol index 3: starting from the fixed Fast root, z12 Reference maps `2.526153432897177e-44` to `2.526174494813028e-44` (`+8.3375e-6`); a second direct z12 map at the updated root returns exactly the same root, residual 0, zero pointwise spectrum change across 77 raw-positive nodes, no warnings. The second map takes `1335.97 s`; quadrature/interpolation estimates are `3.48e-63/5.18e-52`. Together with index 7 this provides two independent root-closure samples, not all eight; six still only have fixed-Fast-root map checks. Artifact/script `cr_positive_sobol_index3_reference_root_closure_20261116.json` / `verify_sobol_index3_reference_root_closure_20261116.py`.
- Current-source cr>0 frequency refinement at Sobol index 7: hold Fast root and z10, directly solve a nested 153-node grid by adding all 76 log-frequency midpoints to the 77 native nodes. Shared 77 spectra match exactly. Direct new nodes differ from the original linear-signal 77-node PCHIP by P95/max `0.04556/0.52462 dex` (max at `log10 f=-18.2978`, weak tail). Reference map changes `−3.5165e-4` between 77 and 153 nodes and the refined fixed point is `−3.6042e-4` from Fast, exceeding the `2e-4` root screen at z10. A second map on the refined grid closes at residual 0 with no spectrum change. QUADPACK relative error is `1.38e-19`; same-grid half-density interpolation estimates are `1.23e-8` (77) and `4.36e-9` (153), so these estimates miss the omitted-node grid error. Existing 77-node z10→z12 map change is only `+1.58e-5`; 153-node z12 is still needed to quantify tail/grid interaction. Artifact `cr_positive_sobol_index7_frequency_refinement_z10_20261117.json`; current-source calculation and corrected replay scripts are hashed separately.
- Sobol index 7 frequency refinement was repeated at z12 on the same nested 153-node direct Reference grid. At the fixed Fast root, the 77→153 map shift is `−3.51162e-4` at z12, close to `−3.51646e-4` at z10; 153-node z10→z12 map shift is `+1.62585e-5`, close to the 77-node tail shift `+1.57748e-5`. Updated-root mapping on 153-node z12 closes at residual 0, while that mapped root differs from Fast by `−3.44163e-4`, exceeding the `2e-4` root screen. For this one discrete sample, the 77-node root agreement does not pass; the 305-node z12 direct Reference grid and root closure are now recorded separately for index 7; this remains a discrete comparison, not a production-grid decision. The two fine z12 Reference solves took `2654.98/2618.47 s`; artifact `cr_positive_sobol_index7_frequency_refinement_z12_20261117.json` and runner are separate from the z10 data.
- Sobol index 7 nested z10 refinement continued from 153 to 305 directly solved Reference nodes. The 153→305 map change is `−4.1022e-5`; cumulative 77→305 is `−3.9265e-4`. All 153 shared nodes match exactly, no warnings; the 152 new nodes differ from 153-node linear-signal PCHIP by P95/max `0.03291/0.24721 dex`. A second 305-node map closes at residual 0, but the 305-node mapped root remains `−4.0142e-4` from the Fast root, outside the `2e-4` screen for this sample. z10 map runtimes are `740.05/742.97 s`; z12 refinement currently stops at 153 nodes, so the 305-node tail interaction remains unresolved. Artifact/script `cr_positive_sobol_index7_frequency_refinement_z10_305_20261117.json` / `verify_sobol_index7_frequency_refinement_z10_305_20261117.py`.
- Independent Sobol index 3 nested z10 77/153/305 direct Reference: map shifts are `−2.3371e-4` (77→153), `−2.4965e-5` (153→305), cumulative `−2.5867e-4`; all shared 153 spectra match exactly and 305 root iteration closes at residual 0. The 152 new-node linear-PCHIP discrepancy has P95/max `0.00766/0.24074 dex`. The 305-node mapped root remains `−2.5651e-4` from Fast, slightly beyond the root screen. Across indexes 3 and 7, 153→305 map changes are within `2e-4` while the 305-node roots remain outside it; this is evidence of a discrete 77-node integration bias at both points, not a continuous-domain claim. 305-node z12 interaction is now closed at index 7; index 3 remains at 153 nodes. Artifact/script `cr_positive_sobol_index3_frequency_refinement_z10_305_20261117.json` / `verify_sobol_index3_frequency_refinement_z10_305_20261117.py`.
- Existing opt-in Fast adaptive grid target root screen at Sobol 3/7 against directly converged 305-node z10 Reference anchors: target `3e-4` gives root errors `+3.349e-4/+3.483e-4` (fail); target `2e-4` gives `+1.517e-4/+1.623e-4` (root-only pass). Direct Reference validation of target `1e-4` on the exact Fast adaptive grids: 175/167 nodes, both updated-root residuals 0, Fast/Reference same-root same-grid max spectrum errors `7.202e-4/8.967e-4 dex`, no warnings and all raw nodes positive finite. Index 7 is close to the `1e-3 dex` limit. z12 adaptive-grid validation at Sobol indices 3 and 7 is now complete; matched runtime/memory A/B remains open, so no default or production change. Artifacts `sobol_fast_adaptive_target_root_screen_20261117.json`, `sobol_index3_adaptive_target1e4_reference_z10_20261117.json`, `sobol_index7_adaptive_target1e4_reference_z10_20261117.json` and scripts.
- Direct z12 Reference validation of the existing Fast adaptive target `1e-4` at Sobol index 7 completed on the exact Fast root and 167-point grid. The initial map residual is `+8.2411e-6`; a second direct map at the updated root returns the same DN_gw (residual 0), with zero pointwise spectrum change. Fast/Reference same-root same-grid spectrum max/P95 error is `8.9672e-4/3.3358e-4 dex`; 167/167 raw signal nodes are positive finite, no warnings, quadrature/interpolation estimates are `2.277e-66/1.823e-55`. Each Reference map costs `3017.56/3020.80 s`. The updated root differs from the 305-node z10 frequency anchor by `+1.1322e-4`; it is not a z12 frequency-convergence estimate. This is one discrete point, close to the spectrum gate, not domain certification. Artifact `sobol_index7_adaptive_target1e4_reference_z12_20261117.json`; exact executed runner snapshot SHA256 `8d8a0929…cb7e8a`, corrected source runner SHA256 `dc27d4d1…2cad823`.
- Independent Sobol index 3 adaptive target `1e-4` direct z12 Reference closure completed on 175 exact Fast nodes. Initial root residual `+8.8804e-6`, updated-root residual 0, same-grid spectrum max/P95 `7.2020e-4/5.0215e-4 dex`; 175 positive finite raw nodes, no warnings. The Fast root differs from the directly refined 305-node z10 anchor by `+1.0311e-4` (not a z12 frequency-converged comparison). Map runtimes `3103.18/3107.49 s`; quadrature/interpolation estimates `3.475e-63/2.763e-52`. Together with index 7, both discrete z12 checks pass root-anchor and spectrum gates, but index 7 has little spectrum margin. Artifact `sobol_index3_adaptive_target1e4_reference_z12_20261117.json`; runner/source hashes recorded.
- Sobol index 7 nested z12 Reference refinement from 153 to 305 direct nodes completed. At the fixed original Fast root, the map shifts `−4.0830e-5` from 153 and `−3.9198e-4` cumulatively from 77; all 153 common spectra match exactly and all 305 raw nodes are positive finite, no warnings. The 152 added midpoint direct spectra vs linear-signal PCHIP from 153 have P95/max `0.032914/0.247211 dex`. A second direct 305-node map closes at residual 0, zero pointwise change; runtimes `5184.10/5182.43 s`. The closed root is `1.654870115968357e-47`; adaptive target `1e-4` z12 Fast root differs by `+9.6773e-5`, within `2e-4`. Evidence is one discrete point and 153→305 only; 609/937 z12 and broader-domain evidence remain open. Artifact `cr_positive_sobol_index7_frequency_refinement_z12_305_20261117.json`; direct runner/source hash recorded.
- Independent Sobol index 3 z12 nested Reference grid closure from 77→153→305 direct nodes completed. 77→153 map shift `−2.3428e-4`; 153→305 shift `−2.6420e-5`; cumulative 77→305 shift `−2.6069e-4`. The 305-grid updated-root map closes at residual 0, all 305 raw nodes positive finite and no warnings. Direct solve runtimes `5299.46/5570.48 s`; PCHIP added-midpoint discrepancy P95/max `0.007661/0.240738 dex`. Closed Reference root `2.5255159374590494e-44`; adaptive target `1e-4` Fast root error `+9.8950e-5`. The 153-node grid’s root misses the 2e-4 screen, while the adaptive candidate passes against the 305 z12 root for this discrete holdout. 609/937 and broader-domain verification remain open. Artifact `cr_positive_sobol_index3_frequency_refinement_z12_305_20261117.json`; direct runner/source hashes recorded.
- Re-ran the underflow 235/469/937 nested-grid screen on the current Fast/Reference sources; all nodes were direct solves. Fast roots move `−2.2160e-6` (235→469), then `−8.7088e-8` (469→937). Direct Reference root-closure maps at 235/469/937 pass after one updated-root map with residuals `3.8e-15/−2.1e-15/0`; Fast vs closed Reference roots differ `8.35e-7/2.32e-7/4.28e-7`. Spectrum max errors `4.804e-4/4.804e-4/4.830e-4 dex`, all below gate. Reference map changes shrink from `−1.613e-6` to `−2.826e-7`; warnings persist and normalized zero counts are 13/25/49. Current artifacts `underflow_grid_235_469_937_current_20261117.json` and `underflow_reference_root_closure_235_469_937_current_20261117.json`.
- Current-source 937 Reference PCHIP diagnostic: 49 normalized zero nodes, 10 positive subnormals, four divide-overflow warnings. Log arithmetic on the same positive linear PCHIP polynomial gives relative integral difference `−1.776e-15`; zero-adjacent panels contribute `10^-322.714`. P95 finite midpoint discrepancy is `3.77e-14 dex`; weak-tail worst case is `0.30090 dex`, and 48 of 936 SciPy midpoint values become zero. Reference production code remains untouched. Artifact `reference_pchip_937_current_20261117.json`.


- 2026-11-17 Round 54: ran the same-grid one-thread amplitude A/B. The first diagnostic failed because its check mixed constructor and Fast evolution mesh sizes; inspected Fast's iteration/data flow, corrected the check, then reran successfully. Standard-vs-forced-log max spectrum `1.6315e-9 dex`, map relative `−1.5654e-14`; product-underflow-vs-standard ratio errors `1.6315e-9 dex` / `9.3259e-15`; exact zero remained `−inf`/zero-map. Confirmed all four returned background hashes match. Re-ran the existing 937 phase A/B artifact and recorded the one-point limits. New artifacts and runners are pending full validation, commit, push, and CI.
- 2026-11-17 Round 54 validation: full pytest `411 passed, 9 deselected` (two existing DeprecationWarnings); Ruff initially flagged one unused import and two diagnostics-only arrays left after fixing the mesh-comparison assertion, removed them and Ruff passed. Mypy passed; sdist/wheel built; installed-wheel smoke passed; Cobaya marker `1 passed, 419 deselected`. Script rerun regenerated the artifact with its final runner hash. Commit/push and remote CI are next.
- Round 54 delivery: commit `b8c8a26db4ffd1a472806771e146f25a3b408c1b` has Author/Committer `2966684515@qq.com`, was pushed to `fast/fast_v0.2`, and remote HEAD matches. GitHub Actions run `44046389606` completed successfully (all nine jobs). The separate final ledger update is being committed now.
- 2026-11-17 Round 55 started from clean current HEAD `c67b2c9`; remote `fast_v0.2` matches. Live Fast source hashes were re-read. Plan: 3 tilts × 4 amplitude-representation strata × 3 reheating temperatures × 3 `kappa10` values (`cr=0`), with exact zero as a separate control; first gather derived amplitude and current default frequency-grid feasibility only.
- Round 55 current-source preflight: 162 `cr=0` combinations across `n_t=17.30/17.38/17.45`, `T_re=10/2000/10000`, `kappa10=.001/.01/1`, six amplitudes (normal-edge, adjacent subnormal, interior subnormal, rounded-zero product, minimum-subnormal r, exact zero). All 162 background/grid preflights succeeded; exact zero uses 73–76 point `goal`, nonzero high-tilt cases choose `construct` at 100–236 points. Corrected the preflight to mirror the production fallback condition after source inspection. In guarded Fast execution at the center background, normal-edge, subnormal-edge and subnormal-interior all hit the existing `shared_Neff_guard`; r=1e-316 and r=minsub both converge (roots `2.74122595e-2` and `1.35434559e-9`), exact zero maps to zero. This motivates adding `n_t=0,1.6` controls for valid normal-amplitude propagation; a guard rejection is not a stability pass.
- Round 55 direct Reference: fixed `cr=0,n_t=1.6,T_re=2000,kappa10=.01`; solved unit-transfer Reference z8/z10 on exact 75/235 Fast grids and scaled nonzero amplitudes in log space. Same-grid Fast/Reference max spectral error `3.42e-4–3.90e-4 dex`; z10 map difference around `−6e-6`. 235-grid minimum-subnormal z10 map re-solve at updated root closes exactly (quadrature estimate `1.1e-14`, interpolation estimate `1.05e-3`, no warnings). Cross-grid direct maps shift `−4.97e-3` (75→235); 75→235 PCHIP new-node error max/P95 `0.542/0.0362 dex`, so 235→469 nested direct refinement is required. Runner's initial cross-grid shape mismatch was diagnosed and corrected before successful completion.
- Round 55 nested Reference complete: 235 shared `Ogw-Oj` nodes match exactly in the 469 direct z10/z12 result; 234 new modes were directly solved. 235→469 map shifts are `−3.23924e-4` (z10) and `−3.20302e-4` (z12), above the `2e-4` gate; 469 interpolation estimate is `3.20404e-4`. Added-node linear-PCHIP max/P95 discrepancy `0.43560/0.00876 dex` (z10) and `0.43560/0.00872 dex` (z12). z10→z12 map shift `−6.987e-7`; direct z12 updated-root closure at minsub is 0, quadrature estimate `1.11e-14`, no warnings. Same-grid 469 Fast/Reference max spectrum error `3.701e-4 dex`, map relative error about `−9.86e-7`. Frequency convergence remains open; next nested direct grid is 469→937.
- 2026-11-17 provenance repair: detected stale preflight runner hash from a post-generation unused-import cleanup. Re-ran 270 preflight rows (all identical except runner hash), reran 12 Fast samples (numerical data unchanged except runtime and lineage), and refreshed the two Reference artifact input-hash fields without touching any direct ODE arrays. All four current runner hashes and source hashes now validate.
- 2026-11-17 Round 55 final local gates: full pytest `411 passed, 9 deselected` (two existing DeprecationWarnings), Ruff all repo, Mypy, compileall, sdist/wheel build, installed-wheel smoke, Cobaya marker `1 passed`. No production source changed. Round 55 artifact lineage was revalidated after regeneration; commit/push and remote CI remain.
- Round 55 delivery verification: commit `1ed916a8716a9d5af192d3b1579bbd0411502833` has the required Author/Committer email, is pushed and verified at `fast_v0.2`; GitHub Actions run `44096271690` completed successfully across all nine jobs, including canonical coverage and schema checks.
- Round 55 CI concluded at current commit `2e9edd833f8467498219b250131fd6aa21baafdb`: GitHub Actions `44097623095` success; remote `fast_v0.2` matches and worktree is clean. 469→937 is next because the 235→469 Reference root shift is still `3.20e-4`–`3.24e-4`.
- Round 56 nested direct Reference complete: 469→937 fixed-root map changes `+8.7905e-5` (z10), `+8.7228e-5` (z12); cumulative 235→937 remains `−2.3310e-4`. 937 z10/z12 direct solves have zero warnings, positive finite raw signals, and exact equality on all 469 shared nodes. Newly solved 468 nodes vs 469-linear-signal PCHIP max/P95 `0.1933/0.00532 dex`. Same-grid 937 Fast/Reference max/P95 `3.743e-4/6.147e-6 dex`; z10→z12 shift `−1.376e-6`; 937 updated-root z12 closure residual 0 (quad `1.11e-14`, interpolation `8.722e-5`). At 2,001 sampled N values the actual Reference H2/sigma were bitwise identical between DN_eff=0 and tested roots; finite probe only. Runner/source hashes pass. Next is an independent T_re/kappa holdout.
- Round 56 local validation: full pytest `411 passed, 9 deselected` (two existing DeprecationWarnings); Ruff all repo, Mypy, compileall, sdist/wheel build, installed-wheel smoke, and Cobaya marker `1 passed` pass. The 937-node artifact source/runner hashes, shared-node exactness, warning-free status, and direct root closure checks pass. Commit/push and GitHub CI are next.
- Round 56 delivery: commit `546a024a9dad209383389a140d69af90fe386d74` was pushed and remote `fast_v0.2` SHA verified. GitHub Actions run `44215527410` completed successfully across all nine jobs, including canonical, static, Cobaya, package, and Python compatibility.
- Round 56 final ledger delivery: commit `b545f8503ea0a4a68584385ccf278f181cde706c` pushed and remote SHA verified; GitHub Actions run `44216116648` succeeded on all nine jobs. Current remote/worktree basis is now verified; Round 57 holdout begins next.
- Round 57 selection from current preflight: low-temperature/low-kappa point (`T_re=10,kappa10=.001,n_t=1.6`) has two amplitude holdouts (`r=1e-316`, minsub `r=5e-324`), both rounded-zero stored `A_t`, finite log amplitude, valid `N_inf=42.3905`, and same direct Fast 221-node construct grid hash. Propagation and Reference remain pending.

## 2026-11-19 — final-output underflow classification

- Confirmed the proposed regression existed: nonzero tensor model returned success with one zero linear spectrum bin, `log10OmegaGW=-inf`, and a divide-by-zero warning. The failing test reproduced the issue before the source change.
- Updated only final output classification in `fast_sgwb.py`: reject nonfinite and nonpositive final spectra explicitly, distinguish a zeroed linear output as `nonrepresentable_output/final_spectrum_underflow`, and restore input `DN_eff`/non-convergence state on rejection.
- Regression validation: focused behavior/zero-limit/log-amplitude trio `3 passed`; affected test files `30 passed, 5 deselected`; full default suite `423 passed, 9 deselected, 2 existing DeprecationWarnings`.
- Static/package validation: Ruff passed; Mypy reports no issues; compileall and `git diff --check` passed; sdist and wheel built; distribution boundary and installed-wheel smoke passed.
- This is a fail-closed output correctness fix, not a new Fast success or expansion result. Commit `55ec264cb7a18d9bb75dcf39909a448eb866f18a` was pushed to `fast_v0.2`; remote SHA matched. GitHub Actions `44474144299` completed successfully across all 9 jobs (canonical, static, package, Cobaya, Python 3.9–3.13). Worktree was clean after delivery verification.

## 2026-11-19 — cutoff/reheating overlap mesh experiment

- Audited current frequency-grid flow: formal `goal` uses a fixed 64-node sparse seed plus 17 feature nodes; the opt-in adaptive grid starts from a separate uniform support and refines once. Existing `goal` errors are concentrated where cutoff support overlaps the reheating feature.
- Prototyped a generic recursive midpoint selector using the existing adaptive-grid indicator `|PCHIP''| h²/8 × interval-integral-share`, with target `2e-5`, then solved only once on the final mesh. It does not inspect named parameter values.
- Fast-only 3x3 screen (`T_re={5,10,20}`, `kappa10={1e-8,1e-6,1e-4}`, `N_inf=740`, cutoff one dex below `f_re`) produced candidate meshes of 72–83 points versus 64 baseline points; all nine roots shifted. The two opposing training corners and the independent interior holdout have direct deep-tail Reference validation.
- One direct z8 Reference corner screen passed root (`1.318e-4` relative) and pointwise spectrum (`7.063e-4 dex`) gates. That exploratory run did not save raw Reference arrays, so the result is being rerun once with incremental artifact persistence. Stopped the remaining eight deep Reference closures to avoid redundant cost; will validate only selected cases/holdout and reuse existing deep-tail artifacts.
- Candidate remains provisional: no production change, test or expansion claim yet. Remaining gates are the reproducible raw direct-reference artifact, independent region holdout, a second valid 2D slice, stable old-region P50/P95, and actual end-to-end adaptive-path cost.

## 2026-11-19 — Fast-only mesh refinement and independent region checks

- Archived the recursive mesh screen on two distinct cutoff/reheating 2D slices: `T_re × kappa10` (3x3 with independent interior holdout) and `T_re × n_t` (3x3 with its own independent holdout), each under mathematically valid cr0 cutoff conditions.
- Same-grid z8 Reference validation is complete for two opposing corners plus the first slice holdout, and for the second slice holdout. All selected cases pass: root relative error max `1.678e-4`, pointwise spectrum max `7.343e-4 dex`, fixed-point residual max `1.612e-4`. Raw arrays and SHA256 lineage are in the corresponding `validate_*.json` files. A previously archived deep-tail anchor has the same T/kappa/n_t/A_s but different DN_re, so it is excluded from candidate validation.
- Rejected `DN_re × n_t` as an independent slice because changing DN_re did not change the solved outputs when the blue cutoff remained active; changed the independent dimensions to `T_re × n_t` instead.
- The un-gated local curvature indicator would refine many old red-tilt trusted points, so global adaptation is unsafe. Next: write a red regression for a condition based on an active positive-tilt cutoff overlapping the existing reheating reserve, implement minimal production selection/refinement, then repeat full source-level end-to-end tests and fixed-resource old-region timing.

## 2026-11-19 — candidate rejected on performance and dispatch fidelity

- Production integration was prototyped in the goal-grid path and covered by red/green regressions. Starting refinement on the first propagation was rejected: the high-root corner moved to `4.8913218094` versus direct z8 Reference `4.8903057849` (relative error `2.08e-4`, over gate). Delaying the mesh probe to outer iteration 1 restored the two root anchors, but did not meet cost requirements.
- With the stricter `2e-5` estimator target, paired warm Fast A/B across two training regions and two independent holdouts measured aggregate P50 ratio `1.0375` and P95 ratio `1.2375`; the high-root corner alone had P50 ratio `1.4428`. Cold-process solver ratios were about `1.03–1.06`, but stable warm P95 remains a required gate. This candidate is rejected.
- Relaxing the target to `1e-4` was also rejected by red root regressions: the `DN_re=740` anchor shifted by `3.33e-4` relative and the high-root corner by `3.34e-4`, both above `2e-4`.
- The z8 same-grid Reference artifacts validate the offline root-informed mesh cases only. They do not validate the later production dispatch's iteration-1 mesh, so are not used as production accuracy proof. `verify_goal_mesh_auto_dispatch.json` was generated from an earlier implementation and is stale.
- Restored all candidate production/test edits. Current code is unchanged from base SHA `872a70d`; only the research ledger and experiment artifacts remain uncommitted. No Fast capability or performance claim is made from this experiment.
- Follow-up exact-reuse prototype: after a red regression, added a selected-mode argument to the Numba kernel and reused shared frequency rows only when DN_eff and all sampled background/primitive arrays were bitwise equal. Incremental versus full propagation returned bitwise-identical output arrays. Repeated the 16-pair warm benchmark: aggregate P50/P95 ratios are 1.3363/1.3390; four cold solver ratios are 1.03–1.08. Dense all-frequency history assembly and expanded output-grid cost dominate enough that exact mode reuse does not meet performance gates. Restored the production source and tests again; base remains `872a70d`. Next low-cost prototype: separate quadrature-only terminal modes from the unchanged public output-frequency grid.

- A terminal-only quadrature prototype kept public nfreq at 64 and passed a high-root Reference root check, but its 16-pair warm P50/P95 ratios still rose roughly 17–22%; it was rejected and removed.
- A final Fast-only fixed-cardinality reallocation probe kept 64 public nodes by replacing the lowest-share non-feature nodes with estimator-selected midpoints. All four selected roots met the Reference root screen (max relative error 1.653e-4), but every baseline input already returned successfully. This is accuracy correction, explicitly out of scope as a domain-expansion result; same-grid Reference spectrum validation is also outstanding. No candidate source or tests remain modified.

## 2026-11-19 — S2 overflow identity candidate rejected

- TDD reproduced two existing candidate defects: `Omega_bh2=1e5, nt=0, T_re=2000` returned root `0.012241370684629289` against the already archived same-grid Reference root `0.012402231615199637` (1.297% relative error), and the integrated terminal `Phi` differed from the H2 identity by up to `384.1746`.
- Derived and checked the smooth-branch identity `sigma = 1 - (d ln H2/dN)/3`; sampled finite differences at `N_re+1,+5,+10,+20` and near the terminal transition agreed with `sigma_vec` (largest sampled difference `2.41e-4` around `N_re+10`, where the derivative is steep). The identity-based primitive removed the artificial terminal `Phi≈432` spike and yielded terminal `Phi≈47.8537` for the `Omega_bh2=1e5` probe.
- The identity candidate still failed the acceptance gate: its Fast root was `0.012034753827431605`, 2.962% from the archived Reference root. A separate `Omega_bh2=1e6` exploratory success root moved from the old forced-relative estimate `0.1074165` to `0.1024309`; that estimate has no independent Reference validation and is not treated as accuracy evidence.
- Therefore neither the relative-scale-only candidate nor the H2-identity integration candidate is retained in production. The screen baseline generated by monkeypatching Numba globals remains invalid and is excluded from A/B evidence. No domain expansion, precision, or speed claim is made from this round.
- Root candidate source/test hashes before rollback: `exact_background.py=F9CB8920B3C5165072159D27383C80CA337A45404EE50D53F066FD15031994E0`, `fast_sgwb.py=6168BDE9EFA57BC66E381B515B80F3EE9B86BA67137C59848472A0FAF554505E`, `test_long_reheating_domain.py=1C50CC92CF0944ADEBF8C8A0DC5657A036D7D17B71733CF1EAEDD05E8BCC51ED`.
- Next research target: isolate the terminal background feature against direct Reference primitives on the exact Fast grid, then determine whether the remaining root mismatch comes from Fast transfer/terminal handling or a Reference/Fast background mismatch. Do not restore automatic S2 scaling until this is resolved.

- Follow-up isolated the mismatch: at 229 common points spanning the full post-reheating interval and a logarithmic terminal cluster, `H2_vec` versus Reference `background_at` differed by at most `2.29e-16` relative and `sigma_vec` versus Reference `sigma` by at most `2.22e-16` absolute. The remaining 2.96% root error is therefore not a Fast/Reference background-formula discrepancy; it lies in the propagation/integration response to the narrow terminal feature or another downstream numerical stage. After restoring production source/tests to `872a70d`, `pytest -q tests/test_long_reheating_domain.py` passed (`11 passed`, 10.90 s).

- Propagation resolution probe used the archived same-grid Reference root without rerunning Reference. With an exploratory H2-identity primitive, uniform `h=0.005,0.0025,0.00125,0.000625,0.0003125` gave root relative errors `2.963%,2.962%,2.081%,1.402%,0.852%`; uniform refinement is far too costly to reach `2e-4` and still misses the gate at 152,873 N nodes.
- A validation-only local terminal refinement based on midpoint curvature of `log(H2)` plus the exact H2-integral primitive lowered root error to `0.273%` at curvature threshold `0.01` with 15,783 N nodes and about `1.62 s` after the first JIT run. Threshold `0.003` gave `0.279%`, 15,822 nodes, about `1.81 s`; this nonmonotonic shift is within the evidence that this candidate is not converged. It remains over the root gate and is slower than Fast, so rejected. No spectrum gate was claimed.
- The one harness error was calling the public `model.SGWB_iter` with an unsupported `config` keyword; the probe was rerun through `SGWB_iter_fast(..., config=...)`. No production source changed.
- Kernel inspection shows regular non-kink steps reconstruct `z_end = 2*z_mid - z_node` and then reset to `Phi_grid` at the next node. Thus midpoint primitive curvature alone may not bound accumulated state transfer error through a narrow sigma feature. Next candidate must estimate the transfer-step defect or use a mathematically equivalent locally scaled transfer, and must validate at fixed root before outer closure.

## Round 59 — 2026-11-19 continued

- [x] Derive and implement stable H2-identity background primitives with exact reheating anchoring and representability checks; preserve invalid H2 classification.
- [x] Rework the nonadiabatic WKB handoff scan into a fused Numba kernel after the first implementation failed the trusted-runtime gate.
- [x] Re-run same-batch base-vs-candidate Fast screens on two cr0 rectangles and a separate cr1 consistency rectangle; replay eight saved cr0 Reference points on the final fused kernel exactly.
- [x] Validate cr0 corners/independent holdouts with same-grid z8 Reference. All selected roots and spectra meet gates; the cr0 sampled regions are eligible for a finite-sample expansion claim.
- [x] Screen cr1 with model-derived tilt/cutoff constraints and select Reference checks. The high corner fails the root gate at z8 and z10 and the spectrum gate at z10. Keep it outside accepted precision coverage.
- [x] Measure trusted P50/P95 on four full Fast cases using 64 alternating warm pairs: candidate/base aggregate ratios `0.8850/0.8603`; exact outputs, with each case below the 5% regression gate.
- [x] Generate paired status/precision figure and raw replay/benchmark artifacts.
- [x] Finish current-source full pytest (`428 passed, 9 deselected`), Ruff, Mypy, compileall, sdist/wheel build, and installed-wheel smoke.
- [x] Refresh final source hashes, stage only this round’s relevant code/tests/evidence, verify Git identity, and commit with the required bilingual format.
- [x] Push commit `5895d29ddbb685f5ed3a04574bf87565cd91e249` to `fast_v0.2`, verify the remote SHA, and confirm GitHub Actions run `44541293060` passed all nine jobs.

### Current delivery boundary

The verified expansion claim is limited to the two sampled cr0 slices; no continuous-domain certification is asserted. The tested high-density cr1 region returns Fast results but fails strict accuracy at its selected high corner and remains unresolved. Do not represent it as an accepted capability until a generic propagation correction passes the same precision/performance gates.


## 2026-11-19 continuation resumed

- Confirmed active checkout remains clean in tracked files at `b236c5632fa06cbd113f11627efadee92947c243`; only prior untracked research leftovers are present. No other worktree was modified.
- Re-read the selected plan/findings/progress, exact `scaled_step` and `solve_kernel` code, CR1 z10 artifact, and prior tensor-trajectory diagnostic. No Reference run and no source change performed.
- Next: recreate CR1 exact-grid primitives and quantify per-panel `Phi[k+1]-2*Phi_mid[k]+Phi[k]` with mode phase exposure.

- Corrected the panel analysis: the prior script accidentally used the exact endpoint, and the first aggregate trigger counter mis-unpacked `kink_index`; neither result informed the decision. Reproduced exact production transfer rule at the full CR1 root and confirmed state error grows near terminal active panels.
- Diagnostic quadratic quarter-panel candidate improves two representative high-error trajectories by 3x–4x at `|delta_z|>1e-5`; no full solve or production edit yet. Next add one focused red test for the curved panel helper, then implement the minimal local dispatch.


## 2026-11-19 — CR1 lower-density region validation

- Tested a generic terminal panel-curvature transfer prototype with red/green coverage. It improved two selected mode-state errors but the full `Omega_bh2=Omega_ch2=1e6` Fast root remained `4.43e-4` from the existing z10 Reference map, slightly worse than `4.35e-4`; removed production code and test. No solver change was retained.
- Audited `Omega_bh2`/`Omega_ch2` source dependencies and focused the CR1 claim on `[1e5,3e5]^2`. Reused the original same-input paired screen: four rectangle points old Fast `0/4`, current Fast `4/4`; separate exact holdout old `0/1`, current `1/1`. All old failures are `nonfinite_background_primitive`.
- Direct same-grid z8 Reference validation (76 nodes) now covers lower corner, cross corner, upper corner, and independent `(2e5,2.5e5)` holdout. Maximum root proxy `4.59633e-5`; maximum spectrum difference `9.07808e-4 dex`; all pass. The reference source hash is unchanged. New raw files and hash manifest are in `docs/research/unified_domain_20261119/`.
- Updated `domain_expansion.md` and `h2_identity_domain_coverage.png` to identify this finite-sample CR1 subregion; the previously screened high-density corner stays excluded. Full pytest: `428 passed, 9 deselected, 2 warnings`; Ruff, Mypy, compileall, diff check, sdist and wheel passed.
- Pending: stage only this continuation's report, figure, scripts, JSON data, manifest, and ledgers; verify author/committer email; commit/push `fast_v0.2`; verify CI and remote SHA.

- Continuation: added four paired CR1 crossing points to enlarge the accepted finite-sample rectangle to `[1e5,5e5]^2`; old/current Fast `0/9 → 9/9` on the 3×3 lattice, plus independent holdout `0/1 → 1/1`. Direct z8 Reference at the newly added 800k total density passes root `3.44e-5`, spectrum `8.70e-4 dex`. Updated report/plot and will refresh manifest before the follow-up commit.

- CI for prior evidence commit `0a0ef4742dadbabd6b3c6c8db57a479e70bf9739` completed successfully: all 9 jobs passed. The 5e5 extension is a follow-up documentation/data commit; its source and Reference code are unchanged.

- CI run `44551497120` for commit `5dcc29326ba862c7f1a550e97f2b9442dc6b3aa3` passed all 9 jobs. Continued with an asymmetric CR1 rectangle and 800k same-sum holdout check; its paired screen is `0/8 → 8/8` plus holdout `0/1 → 1/1`, and Reference reuse passes after fixing the runner import path. New report/manifest/figure changes are staged for the next focused commit.

- CI run `44552905102` for `98486a9781a56c0b4abc4f8766f6bc31ea779fcc` passed all 9 jobs. Continued boundary search: 1.1M and 1.2M paired CR1 corners have direct z8 Reference and both pass; the additional 700k cross sample passes too. Expanded square lattice `{1e5,5e5,6e5}` has old/current `0/9→9/9`; holdout `(5.5e5,5.5e5)` old/current `0/1→1/1`, using verified same-sum Reference reuse. Latest report/plot/manifest and scripts/data are pending a focused commit.

- CI `44554892785` for commit `b71a55f8f7a27f376a30e8dc351dc65c8c9ce9c0` passed all 9 jobs. Direct z8 Reference and a second closure map at the 1.4M boundary point confirmed root error `2.0355e-4` (slightly above gate), so this Fast success is rejected from the accepted domain; spectrum max is `8.696e-4 dex`. Added this boundary evidence and updated the plot/manifest.

- Additional 1.3M Reference closure at `(7e5,6e5)` passes root `1.6731e-4`, spectrum `8.588e-4 dex`; Fast paired baseline fails/candidate succeeds. The adjacent 1.4M square corner remains excluded after closure at `2.0355e-4`. Updated the claimed rectangle and evidence manifest; 600k/700k geometry and 1.3M/1.4M boundary now appear in the plot.

- Commit `57f38761625075100f185f4511e80474f8d6688a` is pushed to `fast_v0.2`; remote SHA matches. GitHub Actions run `44558068590` passed all 9 jobs. Current accepted finite-sample CR1 rectangle is `[1e5,7e5] × [1e5,6e5]`; adjacent total 1.4M fails the root gate by a small margin and is excluded. Production source and Reference source were unchanged in these evidence commits. Next research step is a generic transfer/integration error decomposition around the 1.3M/1.4M boundary, with neighboring samples and no single-point threshold tuning.

- CI `44558068590` for `57f3876` passed all 9 jobs. Additional generic local-curvature propagation candidate was red/green tested and run at the 1.4M boundary. It improved root to `6.22e-5` but worsened max spectrum error to `1.431e-3 dex`; reverted production/test changes. PCHIP decomposition attributes the old root-map gap to spectral samples, concentrated near `log10(f)=-15.62…-14.04`, not quadrature. Exact candidate patch and raw diagnostics are recorded; production source remains unchanged.


## 2026-11-20 — CR0 density × negative-tilt sampled expansion

- On the CR0 `Omega_bh2∈[1e5,1e6] × n_t∈[-0.5,0]` rectangle (`Omega_ch2=0.1201075`, `r=.01`, `T_re=2000`, `kappa10=.01`), the paired 3×3 lattice changed old Fast `0/9` to current Fast `9/9`; two independent interior holdouts changed `0/2→2/2`. All old failures were `nonfinite_background_primitive`; all candidate rows passed the physical guard and had finite positive sampled H2 and finite N_inf.
- Eight independently generated scrambled Sobol points (`seed=20261120`) changed `0/8→8/8` with the same baseline failure. A positive-tilt screen was excluded because `shared_Neff_guard` rejected those points.
- Five direct same-grid Reference z8 checks (two lower-tilt density edges, an interior point, and both independent holdouts) all pass: max closed-root relative error `6.5708e-5`, max spectrum error `8.5869e-4 dex`, max Reference closure residual `5.1316e-8`. Reference source SHA remains unchanged.
- Current Fast warm timings on four representatives over 64 repetitions: P50 `6.26–7.39 ms`, P95 `6.97–8.31 ms`, with bitwise-identical roots, grids, and complete spectra over repetitions. Cold worker startup `1.30 s`; first solve `0.222 s` for a complete 75-frequency output.
- The matching current-source trusted A/B aggregate P50/P95 ratios are `1.0220/1.0441` versus `48bbb96`, but several individual trusted cases exceed 1.05 P95. This tail behavior remains disclosed; no production source changed in this evidence tranche.
- Artifacts and SHA256 manifest: `docs/research/unified_domain_20261119/cr0_density_negative_tilt_manifest.json`; figure `cr0_density_negative_tilt_domain.png`. The sampled box is coverage evidence only, not a continuous-domain certification.


## 2026-11-20 — CR0 blue-tilt amplitude strip

- At `Omega_bh2=1e6`, `Omega_ch2=0.1201075`, `cr=0`, fixed `T_re=2000`, `kappa10=.01`, paired baseline/current Fast is `0/11→11/11` on the `log10(r)∈[-12,-8]`, `n_t∈[.05,.1]` 3×3 lattice plus two interior samples. Every baseline failure is `nonfinite_background_primitive`; every current sample passes the physical guard.
- Eight scrambled log-r Sobol points (`seed=20261122`) give `0/8→8/8`. Direct z8 Reference checks at both amplitude edges for `n_t=.05/.1` and one Sobol holdout have max closed-root error `1.6119e-4`, max spectrum error `8.0253e-4 dex`, and max Reference closure residual `2.3212e-8`.
- A wider positive-tilt candidate was rejected on accuracy: direct root errors at `n_t=.2`, `.3`, and `.5` examples exceed `2e-4` while spectra pass. No CR0 blue-tilt range beyond `[.05,.1]` is claimed.
- Current-source warm Fast on four representatives (64 full solves each): P50 `4.32–4.47 ms`, P95 `4.95–5.13 ms`, outputs exact across repeats; first cold solver `0.252 s`, worker startup `1.45 s`. Trusted A/B limits remain as documented in the prior section.
- No production source changed. Data, Reference arrays, timings, plot, and source hashes are in `cr0_blue_tilt_low_strip_manifest.json`.

## 2026-11-20 — Round 60 production WKB3

- Added a production `_wkb_transfer.py` kernel and routed `_phase_segment` through a local WKB3 attempt only when the established phase estimate is at least 64 substeps; all rejected cases use the original phase-substep map. Reheating-kink splitting is unchanged.
- The first embedded defect used one actual state and differed for basis directions; a red test exposed state dependence. Replaced it with the Frobenius norm of the WKB2/WKB3 2x2 transfer-matrix difference, making the local acceptance metric state-independent.
- The complete CR0 domain A/B is on exact same source baseline `63aee9d`, Windows affinity `[0,1]`, Numba 2/workqueue, BLAS family 1, goal grid, PCHIP, 64 alternating warm pairs per case. Default/high-kappa trusted CR1 have P50/P95 candidate/base `0.9774/0.9833` and `0.9869/0.9395`; three high-tilt CR0 cases have `0.4170/0.4263`, `0.3699/0.3788`, `0.3692/0.3800`. Every paired root is exactly equal; candidate spectrum changes are bounded by `2.433e-9 dex`. No status or frequency-count changes.
- Isolated high-tilt CR0 cases are both old/new successes. This round is a solver cost optimization on the previously accepted region; do not count it as additional domain expansion. Historical z10 Reference errors plus the measured candidate delta remain under root/spectrum limits; Reference was not rerun.
- First call including JIT on fresh source paths: baseline `7.376 s`, candidate `10.173 s` (+38%). This one-time cold compile cost is material and recorded separately from stable warm runs.
- Local production break-even probe: at 64 subdivisions, phase transfer `1.772 us`, WKB3 matrix test `0.335 us`; at 32 subdivisions ratio is 2.85. The 64-substep work gate is cost-derived and input-independent.
- TDD: initial production import test failed because `_phase_segment_wkb3` was absent; implemented and focused suite passed. A second test for state-independent defect failed on the original state-local estimate and passed after matrix-norm change.
- Final full pytest on the final production source: `436 passed, 9 deselected, 2 warnings`. Ruff, `python -m mypy` (CI command), compileall, sdist/wheel build, and isolated installed-wheel high-tilt solve smoke passed. Direct `mypy stiffgwpy_fast` was mistakenly run outside the project-configured file scope and showed 286 existing package-wide errors; rerunning the actual CI command passed. Final checks and installed-wheel smoke were rerun on the shipped source and passed.
- A first benchmark attempt used prior commit `48bbb96`, so its high-density results mixed unrelated source changes; it was discarded. A second same-source attempt showed bitwise identity even in the high-phase case because Numba reused compiled globals from a shared source path; it was discarded. The accepted A/B uses clean isolated source paths, fresh Numba cache paths, same HEAD baseline, and the actual production candidate.
- Raw timing data/hash and figure: `docs/research/unified_domain_20261119/benchmark_wkb3_fullsolve_screen.json`, `wkb3_fullsolve_ab.png`; local timing probe: `benchmark_wkb3_break_even.json`; derivation and limits: `adaptive_wkb3_summary.md`.
- Remaining: final full test/build/smoke, update source manifest after final code, then focused commit/push and verify both fast refs and CI.

- Commit `8b1b53ff11ec98ea740a6bcbda1c774d0102fc67` was created with author/committer `jiuluo-if <2966684515@qq.com>` and pushed fast-forward to both `fast_v0.2` and `fast/fast_v0.2`; both remote refs were verified at that SHA.
- GitHub Actions: run `44635505675` (`fast_v0.2`) and `44635525376` (`fast/fast_v0.2`) each completed `success`, 9/9 jobs. Local final suite is `436 passed, 9 deselected`; Ruff, CI Mypy, compileall, sdist/wheel, and isolated installed-wheel high-tilt solve smoke pass.
- No new Fast success region was introduced by WKB3. The verified production change reduces warm end-to-end time in the accepted high-tilt CR0 cases while preserving roots and staying within the spectrum limit; the one-time cold JIT overhead remains +38% and is an explicit limitation.

## 2026-11-20 — Round 61 CR0 high-density continuation

- Reused the paired-screen harness with pinned old Fast baseline `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45` and current source `702d54347e3c489d1632142dc3151cbd54dd6657`. On the CR0 box `Omega_bh2∈[1e6,2e6]`, independent `n_t∈[-.5,0]`, old Fast `0/11→` current `11/11` on the 3×3 lattice plus two interior points; all old failures were `nonfinite_background_primitive`. Eight scrambled Sobol points (seed 20261123) were `0/8→8/8`.
- Candidate H2 was finite and positive at all 19 points (`min H2≈0.4532139932`), `N_inf` finite, and physical guard accepted every point.
- Four direct, unchanged z8 Reference checks cover the two `2e6` tilt corners, the `1e6,n_t=0` corner, and an independent Sobol holdout. Max closed-root error `7.2607e-5`, max spectrum error `8.0846e-4 dex`, max closure residual `7.7315e-8`; all pass.
- The same-process current-versus-pre-WKB (`63aee9d`) A/B used 64 alternating warm pairs per case, affinity `[0,1]`, Numba 2/workqueue, BLAS 1, goal/PCHIP. In the new region candidate P50/P95 are `26.7–27.4 / 27.8–28.6 ms`, with root differences <=`1.24e-15`, exact frequency hashes, max spectrum delta `5.45e-9 dex`, and candidate/base ratios `0.30–0.34`. The 750k existing-region control runs `25.6/26.4 ms`; thus these sampled new points remain in the same warm-time range. Fresh isolated-module first-call/JIT timing was `9.60 s` baseline vs `19.29 s` candidate and is a material one-time startup cost.
- A stage profile built by wrapping `profile_fast_breakdown.run_once` returned `nonpositive_h2` for otherwise accepted points because its mutable accuracy-mode/profile instrumentation did not reproduce the shared screen contract; this result is excluded. The direct persistent-worker A/B above uses `resolve_config`, matching `physical_guard=True`, and validates every output.
- Raw artifacts, source hashes and figure are `screen_cr0_high_density_extension{,_sobol}.json`, `validate_cr0_high_density_extension_z8.json`, `benchmark_cr0_wkb3_negative_tilt_ab.json`, `cr0_high_density_extension_domain.png`, and `cr0_high_density_extension_manifest.json`; summary at `cr0_high_density_extension_summary.md`.
- This is finite-sample expansion evidence, not a proof over the continuum. Reference checks cover four points only. Reference production code and thresholds were not modified.
- Next: commit/push this evidence tranche, verify both fast refs and Actions, then select another old-Fast failure slice rather than stopping at this region.

## 2026-11-20 — Round 61 CR0 density continuation

- Reused the paired-screen harness with pinned old Fast baseline `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45` and current source `702d54347e3c489d1632142dc3151cbd54dd6657`. On the CR0 box `Omega_bh2∈[1e6,2e6]`, independent `n_t∈[-.5,0]`, old Fast `0/11` and current Fast `11/11` on the 3×3 lattice plus two interior points; all old failures were `nonfinite_background_primitive`. Eight scrambled Sobol points (seed 20261123) were old/current `0/8→8/8`.
- Candidate H² was finite and positive at all 19 points (minimum sampled H² `0.4532139932`), `N_inf` finite, and physical guard accepted every candidate.
- Four direct z8 Reference checks cover `(2e6,0)`, `(2e6,-.5)`, `(1e6,0)`, and an independent Sobol holdout. Maximum closed-root error `7.2607e-5`, maximum spectrum error `8.0846e-4 dex`, maximum closure residual `7.7315e-8`; all pass unchanged gates and Reference code remains unchanged.
- Same-process 64-pair A/B compares WKB3 Fast at `702d543` against stable-background pre-WKB `63aee9d`. New region p50/p95 candidate/base ratios are `.299/.305` at 1.5M,-.5; `.328/.335` at 2M,-.5; `.329/.336` at 2M,0; `.307/.314` at Sobol. New-region candidate warm p50/p95 are `26.7–27.4/27.8–28.6 ms`, roots equal to roundoff, frequency hashes exact, max spectrum delta `5.45e-9 dex`. Existing-region control p50/p95 is `25.62/26.40 ms`. Isolated first calls include JIT and are slower (baseline 9.60s, candidate 19.29s); this remains disclosed.
- A first stage profile using `profile_fast_breakdown.run_once` with `apply_accuracy_mode` produced `nonpositive_h2` even though the same parameter inputs passed the established screen runner. That instrumentation run was excluded from all conclusions; paired persistent-process timing and direct screen/Reference runs use `resolve_config` and the established call contract.
- Raw paired screens, z8 Reference arrays, 64-pair timings, source hashes, manifest, summary and plot are in `docs/research/unified_domain_20261119/cr0_high_density_extension_*` and adjacent `benchmark_cr0_*high_density*` JSON/Python artifacts.
- This is finite-sample evidence supporting the adjacent density extension, not a continuous-domain proof. Current Fast source was not modified in this evidence tranche.
- Next: commit/push only the new evidence/ledger changes and verify the two Fast refs and CI; then choose another old-Fast failure slice.

- Round 62 diagnostic note corrected by the JIT trace in Round 63: the apparent terminal underflow came from a WKB basis determinant rounding to zero inside the parallel mode loop; the earlier long-state dynamic-range hypothesis was not the cause. The initial output-only retry was removed before identifying the actual propagation failure.

## 2026-11-20 — Round 63 WKB basis singularity guard and CR0 extension

- Added a TDD regression using the captured panel `z0=8.326819262766648`, `z1=13.006402132438978`, `h=5.853210263495612e-4`; it reproduced `ZeroDivisionError` before the guard and exact fallback after it.
- The 3M `n_t=0` JIT solve now completes all 75 modes with `DN_gw[-1]=0.3121571825`; direct z8 Reference at that Fast root differs by `3.42e-6` in root and `7.58e-4 dex` in spectrum. This replaces the partial failed solve's invalid `0.002371` value.
- Pinned baseline `872a70d`: the `[2e6,3e6] × n_t[-.5,0]` CR0 3×3 grid plus two interior points is `0/11→11/11`; independent Sobol (seed `20261124`) is `0/8→8/8`, all old failures `nonfinite_background_primitive`. H2 is positive finite throughout the 19 samples.
- Same-Fast-root z8 Reference at 3M `n_t=-.5` and the upper-density Sobol holdout has maximum closed-root relative difference `1.7505e-5`, max spectrum error `7.7100e-4 dex`, max map closure residual `9.1044e-7`. The separately reused 3M zero-tilt check is `3.4198e-6 / 7.5822e-4 dex`.
- Exact old-helper/current-helper 32-pair A/B on two trusted points and one new 2.5M point has P50 ratios `0.995–1.010` and P95 ratios `0.994–1.037`; roots, frequencies and spectra are bitwise identical. New-region screen warm P50/P95 excluding first JIT is `90.9/94.6 ms`.
- Local verification complete: WKB focused tests `6 passed`; full pytest `438 passed, 9 deselected, 2 warnings`; Ruff, Mypy, compileall, sdist/wheel build and installed-wheel guard smoke pass. Source/artifact manifest and figure are generated. The domain evidence is finite sampling, not a continuum certificate.
- Delivery verified: commit `bb4217ba1404679f4f1ecae75ab0e9032e8c6853` is present on both `fast_v0.2` and `fast/fast_v0.2`. GitHub Actions runs `44661121653` and `44661128642` both completed successfully; canonical, static, package, Cobaya, and Python 3.9–3.13 compatibility jobs all passed.

## 2026-11-20 — CR1 square-domain evidence at b175f0b

- 对 pinned baseline `872a70dc5ab3bb2ff9d3c0d78f05850341b42a45` 和候选 `b175f0b3ffceccd2d046808921cb2dceffda1611` 进行同输入 CR1 配对屏幕。`cr=1` 的实际派生倾斜均为 `-0.00125`；`[1e5,7e5]^2` 的 3×3 网格加两个内部点为 `0/11→11/11`，seed `20261120` 的独立 scrambled Sobol 为 `0/8→8/8`。旧版 19 点全为 `nonfinite_background_primitive`；候选 H² 均正有限，76 个输出频点完整且无警告。
- 当前源码与已有 CR1 Reference 工件哈希不同，且复放旧输入时 Fast 根/谱未逐位复现；没有复用旧 Reference 误差。固定当前 Fast 根和当前原始频率数组的新 z8 Reference 三点（高角点、中角点、独立低区留出）及第三次映射显示：最大 Fast 根/Reference 相对差 `1.0923e-4`，最大同根频谱差 `8.9352e-4 dex`，最大第二至第三映射闭合差 `2.9228e-7`，通过原门槛。Reference 源码哈希 `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536` 未变。
- CR1 两点以旧 WKB helper / 当前 helper 交替完整求解各 32 对，P50/P95 最大候选/旧核比率 `1.0221/1.0092`；根、频率、谱逐位相同。隔离新 Numba cache 的第一点冷启动基线/候选为 `21.44/20.08 s`；其余报告使用预热分布。该 A/B 检查的是 WKB 守卫开销，不能解释背景 identity 对旧背景失败的改善。
- 一次 Reference 复跑 Fast 根产生了轻微变化的动态频率网格；闭合验证已改为固定首轮 Reference 记录的同一频率数组，再逐次迭代 Reference map。图表经视觉检查，Ruff 与 `git diff --check` 通过，哈希 manifest 覆核通过。
- 报告、原始配对数据、三个 Reference 样本、第三次映射闭合、A/B、复现脚本及图在 `docs/research/unified_domain_20261119/cr1_b175_square_*`、`screen_cr1_b175_square*`、`validate_cr1_b175_square_z8*`、`close_reference_cr1_b175_square_z8*` 和 `benchmark_cr1_b175_wkb_guard_ab*`。
- 这是 CR1 方域的有限样本证据，不是连续认证；本轮未改生产源代码。下一步沿 CR1 失败边界做邻接矩形与独立留出，并继续核验未覆盖的有效切片。

## 2026-11-20 — CR1 邻接方块与精度边界

- 当前 `51c9030` 源码哈希下，将 CR1 配对采样扩展到邻接 `[7e5,8e5]^2`：3×3 网格 + 2 内部点 `0/11→11/11`，seed `20261122` 的 8 点独立 Sobol `0/8→8/8`。旧版 19/19 为 `nonfinite_background_primitive`，且 H² 正有限；候选 19/19 成功，76 频率点、无警告。CR1 实际倾斜 `-0.00125`，`N_inf=68.81383977` 有限。
- 同源码 z8 Reference 在对角和 `1.5M`、`1.6M` 与邻接 Sobol 高和留出三点逐次映射至第三次：最大 Fast 根相对误差 `1.8543e-4`，最大同 Fast 根谱差 `9.1571e-4 dex`，最大第二至第三 map 闭合差 `1.7022e-7`，均过门槛。最高和 Sobol 的根差 `1.7738e-4`、谱差 `9.1032e-4 dex`。
- 继续 screen `[7e5,1.3e6]^2` 得到 `0/19→19/19` 求解成功，但精度拒绝此整块：`(1.3M,1.3M)` 根/谱误差 `6.7009e-4 / 1.1025e-3 dex`；上方独立 Sobol `6.0345e-4 / 1.0878e-3 dex`；`(1M,1M)` 根误差 `4.4352e-4`。报告只接纳相邻 `[7e5,8e5]^2` 样本范围。
- CR1 上方区域 32 对 pilot 的高角点 P95 比 `1.052`；保留原始 JSON 后扩展到 100 对。Side-wise 候选/旧核最大 P50/P95 比为 `1.0041/1.0239`，根、频率和谱逐位相同。成对 P95 仍观测到尾部尖峰（Sobol `1.0702`），摘要中保留此计时波动风险；冷首调高角点基线/候选 `19.68/18.64 s`。
- 生产源码与 Reference 源码未改。图表、复现脚本、paired screen、当前源码 Reference、map closure、32/100 对 benchmark、双 manifest 与摘要归档在 `docs/research/unified_domain_20261119/cr1_b175_adjacent_*`、`screen_cr1_b175_{adjacent,upper}_*` 及相邻 validate/benchmark 文件中。有限采样，不构成连续域证明。
- 当前 paired CR1 邻接扩域精度与 side-wise 运行门槛均通过；上边界更宽方块因精度失败而拒绝。下一步继续寻找另一片旧 Fast 有效失败区域，并检查零极限与不可表示输出分类。

## 2026-11-20 — CR0 正倾斜 × 小幅度二维扩域

- 当前 `c3c6ab3` 同源码 paired CR0 screen：`log10(r)∈[-12,-8] × n_t∈[.05,.1]`，固定 `Omega_bh2=1e6`、`Omega_ch2=.1201075`。3×3 网格 + 2 内部点为旧版/当前版 `0/11→11/11`；独立 scrambled Sobol（seed `20261123`）`0/8→8/8`。旧版 19 点均 `nonfinite_background_primitive`，H²正有限；候选全成功，`N_inf` 有限、75 节点完整、无警告、物理 guard 接受。
- 现源码 z8 Reference 五点（四角 + Sobol）最大首映射根误差 `1.61197e-4`、同 Fast 根最大谱差 `8.02530e-4 dex`；高倾斜两端及 Sobol 的第三次 map 根误差最大 `1.61218e-4`，最大第二至第三次闭合 `3.89572e-8`。Reference SHA `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536`，源码未改。
- 100 对旧/当前 WKB helper 全求解 A/B：side-wise P50/P95 最大比 `1.0267/1.0145`，根/频率/谱逐位一致；成对时间比 P95 尖峰为 `1.156/1.182`，保留时间抖动风险。候选屏幕 P50 `77.3 ms`。
- 分类审计的 12 项定向测试通过；新增一条 underflow 输出断言，确认 14 个线性 `Ogw_today=0` 的 bins 对应有限 log 谱 `[-367.439,-325.375]`，小于最小 subnormal 的 log 下界 `-323.306`。精确零路径另测为全 `-inf`。新断言及 exact-zero/final-underflow 两项复跑 `3 passed`；Ruff 通过。
- 区域屏幕、原始 Reference arrays、map closure、100 对 benchmark、图和 manifest 已存于 `cr0_b175_amplitude_tilt_*`；manifest 覆核 current Fast/Reference/test 源码哈希。生产 Fast 与 Reference 方程未改。
- 这是 finite sampling，不证明连续矩形。下一步检查正向最终溢出与不可表示输出分类，再选另一片旧版 Fast 失败区域；保持当前精度/性能区域不退化。

## 2026-11-20 — final tensor-power output overflow classification

- Added the red regression for CR0 `n_t=1e308`, positive finite `r/A_s`, finite `N_inf`, and positive finite H². Before the fix Fast returned generic `nonfinite` after nonfinite DN map formation, with expected but noisy overflow warnings.
- Fast now detects a nonfinite `candidate_log_P_t` before propagation completion and returns `nonrepresentable_output` with detail `tensor_power_log_out_of_range`; `freq_adaptive` and the Fast power expression locally suppress expected IEEE overflow warnings while the caller classifies the result. Input `DN_eff` is restored and `SGWB_converge` remains false.
- The pressure-domain classifier now maps this reason to `unrepresentable_final_output` and retains `failure_detail`. The red test for reason mapping plus evaluator detail both failed before this script fix.
- Related verification after the fix: `29 passed`; Ruff, project `python -m mypy`, compileall and diff check passed. `Reference` unchanged. This repairs failure telemetry/typing; it is not counted as a domain expansion.
- Next: rerun the CR0 strip's affected-source replay after this Fast source change; then continue the distinct valid failure-domain search and keep positive final-overflow separate from intermediate overflow and `H²≤0`.

## 2026-11-20 — CR0 区域在输出分类修复后的源码复放

- 当前源码相对 CR0 条带证据的代码变更仅涉及张量 log-power 非有限分类与 pressure detail 输出。为了保留哈希边界，新增 after-guard replay artifact，不把 `c3c6ab3` Reference 直接贴到新源码哈希上。
- 初次 replay 脚本从主 worktree 导入了另一份 `stiffgwpy_fast`，产生错误 `4.056 dex` 谱差；复核发现脚本未把目标 repo 根加入 `sys.path`。修正为先插入当前 checkout 根后重跑：19 个 lattice/interior/Sobol 参数当前全部成功、根逐点相同；五个 Reference 参数的频率网格和 Fast log 谱均逐位相同，max delta `0.0`。旧错误 replay 输出已覆盖，不作为证据。
- 正确 replay 记录当前 worktree commit 和 `exact_background.py/fast_sgwb.py/freq_adaptive.py/_wkb_transfer.py/reference.py` 哈希；Reference 源码仍与前序 z8 工件相同。脚本、JSON、摘要与 CR0 manifest 的 post-classification replay source hash 已归档。

## 2026-11-20 — CR1 `T_re×kappa10` paired 区域证据（Round 68）

- Pinned 旧版与当前版同输入，`Omega_bh2=Omega_ch2=7e5`、`T_re∈[2000,2500]`、`kappa10∈[.01,.015]`、`cr=1`。实际 tilt `-.00125`，有限 `N_inf`、正有限 H²；lattice/interior `0/11→11/11`，独立 Sobol `0/8→8/8`，19 个旧失败全为 `nonfinite_background_primitive`。
- 当前源 z8 Reference 四点 max root proxy `1.54486e-4`、max same-root spectrum `8.07271e-4 dex`；两点 third-map max root delta `1.54487e-4`、second→third map closure `1.05017e-8`。
- 100 alternating full-solve pairs, two points: maximum side-wise P50/P95 ratios `0.9997/0.9887`; root/frequency/spectrum outputs match exactly. Paired P95 ratios `1.0717/1.0685` indicate measurable tail timing noise.
- Current-source replay of shared base corner is bitwise exact; saved Reference arrays are reusable. Reference source remains byte-identical.
- 5 source hashes and 20 artifact hashes verified; report, scripts, raw JSON, manifest and figure are in `docs/research/unified_domain_20261119/cr1_b175_reheat_kappa_*` and matching screen/validate/close/replay/benchmark files.
- Evidence covers finite samples only; no continuum claim. Production Fast code is unchanged in this tranche.
- 本轮交付门禁：`python -m pytest -q` 为 `441 passed, 9 deselected, 2 warnings`；仓库 CI 标准 Ruff、根目录无参数 `python -m mypy`、compileall 与 `python -m build` sdist/wheel 均通过。手动向 mypy 指定两个文件触发的额外告警不属于仓库标准配置，改用 CI 命令后 `Success: no issues found in 2 source files`。

## 2026-11-20 — CR1 低温高 κ 邻接区与 WKB 摘要复用优化

- 对 `cr=1`、实际派生 `n_t=-0.00125`、`Omega_bh2=Omega_ch2=7e5` 的 `T_re∈[1750,2000] × kappa10∈[.0125,.015]` 做同输入配对筛查：格点/内部点 `0/11→11/11`，独立 Sobol `1/8→8/8`，总计旧版 18 个 `nonfinite_background_primitive` 失败变为候选 19/19 成功；旧版另 1 点原本成功。
- 为去掉 WKB 尾部屏障对完整背景数组的重复扫描，在现有 `fast_phi_s2` 生成循环中累计 `Phi` 单调性、最大正向步长与最大绝热导数，并仅在原始背景路径将等价摘要传给 WKB kernel。H² 身份替代路径传无效标记并保留原数组扫描。同步用 `Psi` 与 `S2/S2inv` 可表示区间交集减少指数范围比较；未改变方程、守卫或输出。
- 测试先红后绿发现并修复了摘要漏计末节点左侧单边 sigma 的问题。最终 19 点相对 `fda343d` 根/频率摘要/频谱逐位复放；max root 和谱差均 0。Reference raw arrays 与 Reference 源码哈希 `4e2eb80709da944923206a840c0864d454faad9b51dc3fb9ebfda90abe460536` 经复核未变，按逐位输入输出复用已有 z8 与第三 map 核验。
- 精度最大根误差代理 `1.7512e-4`、谱误差 `9.0787e-4 dex`、第三 map 最大根差 `1.7513e-4`、second-to-third map 闭合 `6.1895e-8`。
- TBB、2 Numba threads、CPU `[0,1]`、BLAS1、goal grid、`physical_guard=True` 的隔离持久进程 A/B：相对本轮起点 fda，500 对可信 CR0 P50/P95 `0.9963/1.0110`，500 对扩域旧成功点 `0.9984/0.9929`；新增低角点与 Sobol 100 对分别 `0.9982/0.9980`、`0.9984/1.0079`。输出全部逐位相同。
- 最终源码全量 pytest `445 passed, 9 deselected, 2 warnings`；Ruff、Mypy、compileall、sdist/wheel build 通过。有限样本证据不构成连续域认证。复现脚本、逐位 replay、完整运行原始数据、Reference arrays、图与哈希见 `docs/research/unified_domain_20261119/cr1_b175_reheat_kappa_upper_candidate_manifest.json`。
