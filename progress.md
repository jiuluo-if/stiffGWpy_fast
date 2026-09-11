# Progress Log: stiffgwpy_fast 单一 fast 生产求解器

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
