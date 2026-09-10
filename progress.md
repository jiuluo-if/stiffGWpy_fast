# Progress Log: stiffgwpy_fast 单一 fast 生产求解器

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

### Errors

| Error | Resolution |
|-------|------------|
| 计划文件首次补丁未匹配 init-session 生成的简化模板 | 重新读取实际文件后以完整替换方式更新 |
| 直接查询 `numba.threading_layer()` 时线程层尚未初始化并抛出 `ValueError` | 将线程层初始化纳入基准 warmup，再记录实际 layer；不把该探测异常误判为 solver failure |
| PowerShell 内联 Python 装饰器/多行函数探测命令发生 `SyntaxError` | 不再用复杂 `python -c`；改用已有脚本或通过 `apply_patch` 创建可复现诊断脚本 |
| 固定 profiler 首次用 `from scripts import profile_fast_breakdown` 导入失败 | `scripts` 无包初始化文件；改为把脚本目录加入 `sys.path` 后直接导入 |
| 全仓 `ruff check stiffgwpy_fast scripts tests` 暴露 200+ 个既有 lint 错误 | 未对无关旧文件做大范围格式化；改为对本轮新增/修改 profiling scripts 做 scoped ruff，结果通过，并保留全仓门禁未通过事实 |
| quadrature 比较脚本先导入 SciPy/NumPy 后才设置 `NUMBA_THREADING_LAYER`，实际 layer 为 `tbb` | 该次输出作废；已将环境设置移到数值库导入前，重新运行时强制验证 `workqueue` |
