# Findings & Decisions: fast 单档位重构

## Requirements

- 当前分支 `codex/fast_v0.2`，远端 `fast` 指向 `https://github.com/jiuluo-if/stiffGWpy_fast.git`；重大 Phase 完成并 fresh verification 后提交、推送到远端同名分支。
- 目标是彻底合并 `fast` / `production`，最终只保留唯一正式 `fast`；`reference` 仅作 precision oracle，validation/debug 仅内部使用。
- Phase A 必须只加入 exact `N_re` breakpoint；保持 `h=0.02`、`z_tail=5`、`phase_max=0`、原频率网格、outer tolerance、`col_step` 不变，并测量 runtime、spectrum/DN_gw 误差、transition 最大误差、400 点 guard/failure。
- 后续优先级：短振荡区间的 phase/WKB envelope、goal-oriented sparse frequency grid、outer self-consistency 加速、解析 stiff/RD branch，最后才做 profiler 证明的底层优化。
- 默认 float64；不得静默 NaN、放松 guard、改变 fallback、物理参数含义、输出字段或 Cobaya derived quantities。
- 精度目标：spectrum median <1e-3、p95 <3e-3、尽量 max <1e-2；DN_gw median <5e-4、p95 尽量 <1e-3；速度不明显慢于当前约 4.4 ms/point，第一阶段目标 ≤4 ms/point。
- 每个 Phase 单独提供 commit、runtime、oracle accuracy、parameter sweep、profiler、A/B 结果及接受/回退结论；最终更新 README、Cobaya YAML、schema、tests、manifest、benchmark docs、CHANGELOG，并生成优化报告。

## Research Findings

- `README.md` 当前仍明确写着“两种 user-facing fast profiles”：`fast` plain-grid 与 `production` transition-refine；用户 API、Cobaya YAML、benchmark 脚本和 tests 仍暴露 `production`/`transition_refine`。
- 当前关键代码入口：`fast_sgwb.py` 的 `normalize_accuracy_mode`/`apply_accuracy_mode`/`resolve_config`/`_SGWB_iter_fast_impl`/`SGWB_iter_fast`；`stiff_SGWB.py` 的 `SGWB_iter`；`reference.py` 的 `run_reference`。
- 当前仓库 HEAD 为 `b83aa89`，最近提交为“统一项目名称为 stiffgwpy_fast”；工作树在审计开始时无未提交代码改动，之后仅生成了计划文件。
- 本轮 2026-09-10 已 `git fetch --prune fast fast_v0.2`；远端 `fast/fast_v0.2`、本地 `codex/fast_v0.2` 与 HEAD 均为 `767056d2ea2e4d25f06670f1fecc7c226d85cc9a`。
- 本轮环境为 Python 3.11.9、NumPy 2.4.4、Numba 0.67.0、SciPy 1.17.1、32 CPUs；首次查询 Numba threading layer 时尚未初始化，正式 benchmark 必须在 warmup 后记录实际 layer。
- 当前仓库已有 `exact_background.py`，名称与注释表明它支持 continuous-sigma expansion integrals，但是否已接入轻量 breakpoint 仍需沿调用链核实。
- 当前测试已有 `test_modes.py`、`test_engine.py`、`test_fast_sgwb.py`、`test_reference.py`、`test_freq_adaptive.py`、`test_cobaya_adapter.py`，可作为兼容性与 oracle 骨架。
- 文档基线记录 plain-grid warm median `4.442 ms`、transition-refine `21.772 ms`；这些是历史快照，必须重新运行后才能作为本轮证据。
- 本机 fresh baseline（2026-09-09，Python 3.11.9 / NumPy 2.4.4 / Numba 0.67.0 / SciPy 1.17.1，32 CPUs）：`bench_fast.py --reps 3 --cases 0` 的 plain fast warm median `2.461 ms`、min `2.281 ms`、cold `0.200 s`，fallback `0`；这是当前机器快照，不替代历史基线。
- 本机 fresh reference 对照（default、reference `z_tail=5`、30 个 subset bins）：当前 production fast cold `0.313 s`，fast `DN_gw=0.002270462664`；independent reference `DN_gw=0.002796635640`；spectrum linear relative median `3.44e-3`、max `4.06e-1`，dex p50 `1.50e-3`、p95 `3.39e-2`、max `1.48e-1`。该结果再次证明当前 fast/production 仍未达到目标，且低幅度/transition 区域需按 signal mask 与全量指标分别解释。
- 远端 `fast` 实际 refs 为 `refs/heads/fast_v0.2` 和 `refs/heads/main`；审计初始快照均在 `b83aa89`，当前 `fast_v0.2` 已推进到 `4436b9d`。本地分支名为 `codex/fast_v0.2`，推送目标必须明确使用 `fast HEAD:fast_v0.2`。
- `scripts/benchmark_phase_a.py` 原先固定使用 31 点 reference subset；该 subset 会改变总 `DN_gw`，不能作为总量 oracle。脚本现默认以候选完整频率网格调用独立 reference，稀疏 subset 仅通过显式参数保留作诊断。
- `exact_background` 的 reheating primitive 已修正为保留 `N_re` 的左右单侧 sigma 极限：结束于断点的 Simpson panel 用 pre-transition `sigma=1`，开始于断点的 panel 用 post-transition 值。
- 完整 246 点 reference（固定候选 `DN_eff`，`z_tail=5`, `rtol=1e-11`）显示 Phase A plain -> kink split：`DN_gw` 相对误差 `1.069% -> 0.411%`，spectrum dex median `2.653e-3 -> 1.690e-3`，p95 `1.004e-2 -> 4.745e-3`，max `1.311e-2 -> 5.925e-3`，transition max `3.066% -> 1.317%`。这是有效的完整频率网格结果，但尚未满足最终 spectrum p95 `<3e-3`。
- Phase A thread=1 median timing（fixed `h=.02,z_tail=5,phase_max=0,col_step=8`）为 plain `7.713 ms`、kink `10.020 ms`；16-thread profiler 的 JIT 后 steady sample 为 plain `3.701 ms`、kink `7.441 ms`，主要额外成本落在 exact background primitive，而非 tensor kernel。
- Phase B `phase_max=.25` 在 `h=.02` 下对完整 oracle 几乎不改变 spectrum p95（kink `4.7629e-3`），但已修复 split segment 未应用 phase cap 的逻辑缺口。
- Phase C goal grid 的 seed80 历史结果为 90 点、spectrum dex p95 `2.598e-3`、DN rel `1.677e-3`；当前 seed64 结果与节点/误差取舍见 Latest evidence。
- 当前正式 `fast` preset（`h=.005,col_step=8,z_tail=5,phase_max=.25,freq_grid=goal,kink_split=True`）在本机 16 threads 为 warm median `6.727 ms`、cold `1.919 s`；高于原先 plain `2.654 ms` 与第一阶段 `<=4 ms` 目标，必须继续 profiler/outer iteration 优化，不能宣称速度达标。

## Latest evidence (2026-09-09)

- 外层首轮只需计算最终频率列来得到 `DN_gw`，不需要写入所有中间 e-fold 列；已改为首轮 `assemble=0`、后续完整轮 `assemble=1`。新增回归测试确认首轮/最终轮调用顺序，最终数组摘要保持一致。
- 该热路径改动的 like-for-like plain benchmark（16 threads，A，3 repeats）warm median 从旧记录 `6.727 ms` 降至 `5.051 ms`；仍未达到 `<=4 ms/point`，不能视为速度目标完成。正式 kink 路径仍需继续优化。
- 本轮新增代码经完整回归 `111 passed, 6 deselected`、ruff、mypy、中文注释门禁和 manifest 校验通过；一次未形成稳定收益的 sigma 节点复用实验已回退。

- CI run `34303864035` 的五个 test matrix job 均在“中文注释门禁”失败；根因是 `78cc5b9` 新增 5 行英文 `#` 注释。仅补充中文语义后提交 `4436b9d` 并推送 `fast_v0.2`，CI run `34304423555` 的 Python 3.9–3.13、ruff、mypy、pytest、manifest、wheel、distribution、smoke 和 Cobaya 全部 PASS。
- 正式 fast 的 seed64 A/B（`h=.005,z_tail=5,phase_max=.25,kink_split,goal`，76 点，full candidate-grid reference）为 spectrum dex median `6.528e-4`、p95 `2.719e-3`、max `3.040e-3`，DN rel `1.368e-3`；seed80 的 90 点结果为 p95 `2.598e-3`、DN rel `1.677e-3`。seed64 保留为较低节点且 DN 更好的当前基线，但最终 DN gate 仍未满足。
- 当前 seed64 formal warm benchmark（16 threads，3 repeats，default A）为 median `6.727 ms`、cold `1.919 s`、fallback `0`；仍高于 `<=4 ms/point` 目标。
- 尾部隔离（reference 固定 `z_tail=5`）显示 fast `z_tail=3/4/5` 的 spectrum p95 分别约 `2.10e-2/5.07e-3/2.72e-3`，DN rel 约 `1.48e-2/3.31e-3/1.37e-3`；不能以浅尾部换速。
- 四阶两点 Magnus A/B：`h=.01` 下 DN rel 可到 `1.90e-4`，但 spectrum p95 `3.998e-3`；`h=.005` 热路径约增至 32.3 ms（thread=1）而精度无实质收益，已回退。基于 `Phi` 曲率的局部子步也未改善 h=.01 谱 p95，已回退。
- 高层 API 合并提交 `c9110c9`：`accuracy_mode='fast'` 成为唯一正式用户档；旧 `production`/`ultra-fast`/transition-refine 别名在高层发出 `DeprecationWarning` 并映射到 fast，底层 validation preset 仍可供内部脚本使用。完整测试 `111 passed, 6 deselected`，CI run `34313611135` 全部 PASS。
- rejected numerical A/B：从视界前 3 个十进制数量级提前到 4 个只改变 DN 约 `3e-9`；标准四阶 RK4、两节点四阶 Magnus、视界局部 2/4 子步和 sigma 变化触发的局部半步均未缩小 `h=.005` 与 `h=.0025` 的 DN 差异，且部分方案变慢，均已回退。
- rejected sparse goal refinement：76 点增加到 87 点后，spectrum dex p95 从 `2.719e-3` 改善到 `2.675e-3`，但 integrated DN rel 从 `1.368e-3` 恶化到 `1.755e-3`，已回退。32 个 h/2 探针通道可以估计默认点的 h 偏差约 `7.35e-4`（完整 h/2 对照 `7.44e-4`），但单次探针仍约 `10 ms` 且会造成 endpoint-only 修正，尚未形成可接受的输出一致性方案。
- 尾部匹配 A/B：当前 `gamma=(3-1.5*sigma)/2` 是阻尼去除变量的系数；由物理张量解 `T~a^-1` 推导，深亚视界振幅匹配的首阶组合应为 `x+y/omega`，即 `gamma=1`。只改这一项，在默认点 16 个代表 native 节点上，z_tail=5 的 spectrum dex p95 由约 `2.63e-3` 降至 `3.20e-4`，max 由约 `2.66e-3` 降至 `4.64e-4`，warm fast 约 `9-10 ms`，无新增计算分支；保留该修改。
- 尾部修正跨参数探针：默认/high-T_re/stiff 的 30 个 native 节点 spectrum dex p95 分别约 `9.59e-4/9.14e-4/1.02e-3`；low-T_re 约 `3.35e-3`，误差集中在 `log10(f/Hz)=-18..-16.5`，这些模式未进入 analytic tail。low-T_re 的 h=.005 -> .0025 将该局部 p95 约 `3.42e-3 -> 1.55e-3`，但局部 phase 二分、起点提前、终点直接取 Phi_grid[k+1] 均未带来稳定收益，均回退。
- 尾部修正改变了旧 validation 测试中的固定数值口径：production transition-refine 的 continuous-sigma 分支为 `DN_eff=0.00227259`，普通 grid 分支为 `0.00225706`；测试现改为要求 continuous-sigma 分支接近 reference 且优于普通 grid，避免沿用旧 tail match 的过时绝对阈值。
- 尾部匹配修正 commit `da4fd02` 已按 `2966684515@qq.com` 提交并推送到 `fast/fast_v0.2`；远端 CI run `34318085109` 的 Python 3.9–3.13、ruff、mypy、manifest、wheel、distribution、smoke 与 Cobaya 全部 PASS。
- 背景节点缓存 A/B：`exact_phi_s2_split` 现在可直接复用 `m.sigma` 的节点值，同时保留原有形状检查和 reheating 左右单侧约定；与原始 `sigma_vec` 路径逐项比较到 `2e-13` 以内。formal kink profile 的 steady warm median 约由 `10.62 ms` 降至 `8.90 ms`，但仍高于 `<=4 ms` 目标。
- outer self-consistency 审计：default、low-T、high-T、stiff 和 low-r 五个代表点的首轮 probe 到第二轮 full solve 的 `DN_gw` 绝对变化分别约为 `7.7e-14`、`4.3e-9`、`7.7e-12`、`6.2e-13`、`7.8e-16`；说明 DN 本身对外层更新很不敏感。但临时强制“一次 full solve 即结束”时，最终 `g2/w2` 与当前两轮流程仍出现差异，high-T 点 `g2` 总和约差 `3.0%`，因此不能仅凭 DN 收敛就删除第二次 full solve，当前 outer probe 方案保留。
- rejected background batching：将 exact split 的 midpoint 与 quarter-point sigma 合并为一次更大的向量调用，数值结果保持等价，但 formal kink profile 的稳定 warm median 约由 `8.4 ms` 上升至 `10.6 ms`；已回退，不纳入主线。线程扫描 `1/2/4/8/16/32` 的 warm median 约为 `13.4/9.96/9.90/7.08/8.76/9.91 ms`，稳定 profiler 的 8/16 线程均约 `8.4 ms`，暂不改变默认线程数。
- formal kink preparation 优化：正式 kink 路径不再先计算随后被 exact primitive 覆盖的 uniform-grid `Phi/S2` spline，只保留频率起点、`f_hor` 和 tail 所需准备；新增测试确认两轮 outer solve 均使用 frequency-only preparation。与临时恢复 full preparation 的同环境 10 次 profiler 对照，steady warm median 约 `9.68 -> 8.98 ms`，所有 spectrum/DN/g2/w2 digest 保持一致；保留该修改。
- exact primitive 条件复用：外层更新后同时检查 `Nv`、节点 `sigma` 和 `f_hor` 的最大变化，均低于 `1e-4` 才复用上一轮 primitive；高温、stiff、高 r 和高 kappa 点会自动重新计算。代表点 A/B 中 default 的 `DN_gw` 变化约 `1.1e-12`，相对普通两轮 full primitive 的 spectrum 最大差约 `2.98e-4 dex`；high-T 最大差约 `3.8e-2 dex` 但被门限排除。formal default warm median 进一步约降至 `5.96 ms`，仍未达到 `<=4 ms`，该修改需在完整参数扫描中继续复核。
- 当前独立 reference（h=.005、goal 76 点、gamma=1、phase_max=.25）为 spectrum dex median `1.940e-3`、p95 `2.945e-3`、max `3.190e-3`，DN rel `4.040e-3`；seed80/96 增加到 90/106 点未改善 DN（约 `3.93e-3/4.20e-3`）。phase_max 从 `.25` 降到 `.125/.0625` 也基本不变，已不作为下一步方向。
- z_tail 诊断：从 5 提到 6 只把 spectrum dex p95 由约 `2.945e-3` 改到 `2.938e-3`，DN rel 仍约 `4.044e-3`，没有形成值得增加步数的收益，继续保留 z_tail=5。
- 连续背景 primitive 的进一步优化：正式 fast 新增 `fast_phi_s2_split`。它复用已生成的节点 sigma，在平滑区间用节点线性值构造 midpoint/quarter-point；只有包含 `N_re` 的区间继续调用连续-sigma 探针并按左右单侧 Simpson 分裂。默认点与原 `exact_phi_s2_split` 的 DN 相对差为 `3.22e-12`、频谱最大差为 `2.80e-6 dex`，完整回归为 `116 passed, 6 deselected`。
- 该 primitive 优化后的独立深尾 reference（76 点、reference `z_tail=8`、`rtol=1e-11`）为 spectrum dex median `1.3968e-4`、p95 `4.1466e-4`、max `1.4152e-3`，DN rel `7.1404e-4`；因此此前 reference `z_tail=5` 的约 `4.04e-3` DN 差异主要混入了 oracle 自身浅尾部误差，不能据此继续加深 fast 的 tail。
- 速度复测：新 primitive 在 20 threads 的代表性 warm median 约 `4.49 ms/point`，偶发样本低于 4 ms，但稳定中位数仍未达到 `<=4 ms/point`；继续优化方向应放在 tensor kernel/outer probe，而不是再次增加背景 spline 精度。
- outer full-solve A/B：正式 `goal + kink` 路径首轮改为完整组装；当更新后的节点 `sigma` 与 `f_hor` 最大绝对变化都不超过 `1e-4` 时，复用首轮完整结果并跳过第二次 kernel。默认点调用从 `[probe, full]` 变为 `[full]`，low-T/low-r 也能在首轮收敛；high-T/stiff 背景变化超过门限，仍执行两次完整求解。stiff 单线程逐数组比较的频谱最大差约 `2.71e-10 dex`、最终 `DN_gw` 最大差约 `1.18e-14`。完整回归 `117 passed, 6 deselected`，但正式 warm profiler 中位数约 `5.37 ms`，仍未达到 `<=4 ms/point`。
- outer full-solve 实现同时在每次完整组装前清零 `Ogw/Oj/Opgw`，避免 outer 更新后 horizon 起点移动造成旧列残留；compatibility/validation 的非 goal 路径继续保留旧 probe/full 语义。ruff、mypy、manifest 和中文注释门禁均通过。

## Fresh HEAD profiling (2026-09-10)

- 在远端 HEAD `767056d` 上使用 Numba `workqueue`、固定 CPU affinity（20 threads 使用 CPU 0-19；32 threads 使用 CPU 0-31）、76 goal nodes、case A、kink split、7 次运行（首轮 cold JIT，后 6 次 warm）重新测量。所有 warm digest 一致，fallback 为 0。
- formal warm runtime（ms/point，median；min/p95 见 `docs/profile_head_t*.json`）：1/2/4/8/16/20/32 threads = `9.78/7.28/6.66/5.04/5.88/5.47/4.81`；对应 cold 首轮约 `6.56/0.37/0.23/0.24/0.26/0.24/0.21 s`。当前最好是 32 threads 的 `4.81 ms`，仍未达 `<=4 ms`。
- 20-thread formal profiler warm 分层中位数约：background/gen_fast `0.863 ms`，goal construction `0.436 ms`（2 calls），frequency preparation `0.092 ms`（2 calls），`fast_phi_s2_split` `0.566 ms`，frequency weights `0.012 ms`，solve_kernel `1.035 ms`，bolometric column integration `0.265 ms`；阶段计时之和小于 total `5.006 ms`，剩余约 `1.7 ms` 属于 outer/control、array allocation/zeroing、最终 output/history assembly 与 Python dispatch 的组合，不能只优化 kernel 推断总收益。
- reuse A/B（同样 20 threads、同一 case）：启用 outer full-solve reuse 为约 `5.006 ms`、`solve_kernel` 1 call；禁用为约 `5.705 ms`、`solve_kernel` 2 calls。两者输出 digest 有差异但均收敛，说明当前 HEAD 的 reuse 是实测加速而非 runtime 回退来源；旧 `4.49 -> 5.37 ms` 差异不能归因于 reuse，仍需用同一 benchmark 口径比较线程层、affinity 和写入/缓存差异。
- 已新增可复现 profiling 入口 `scripts/run_fixed_profile.py`，并扩展 `scripts/profile_fast_breakdown.py` 记录 goal construction、frequency preparation 与 `fast_phi_s2_split`；该工具不改变 solver 行为。
- fresh full reference（default，reference z_tail=8，242 reference nodes，DOP853 rtol=1e-11）给出 `DN_gw=0.002262832966946746`；同一 HEAD formal fast 的 full-grid scalar `DN_gw=0.002263022064729132`，相对误差 `8.36e-5`。该 comparison 的 reference 自有 quadrature/interpolation estimates 为 `3.11e-22/3.63e-11`。
- PCHIP candidate 已按 TDD 以显式 `frequency_quadrature='pchip'` 接入，但默认仍为 Simpson。default same-spectrum DN 从 `0.00226269695` 变为 `0.00226364742`，相对 fresh full reference 约 `3.60e-4`；primitive/outer reuse A/B 仍远小于该量级。
- PCHIP-vs-Simpson same-spectrum relative deltas across probes：default `4.20e-4`、low-T `1.275e-2`、high-T `7.40e-4`、stiff `1.282e-3`。因此 PCHIP 是直接针对 DN 的候选，但不能未经 parameter-space reference 验证就替换默认积分；low-T 的差异已经明显超过目标。
- PCHIP candidate 的 focused test 与全量回归已通过：`119 passed, 6 deselected`（另有既有 2 个 deprecation warnings）；默认 Simpson digest/behavior 保持兼容。

## Continued fresh audit on c35db7d (2026-09-10)

- 六点 standardized matrix（default/lowT/highT/stiff/low_r/high_kappa，20 threads、CPU 0-19、workqueue、25 次，其中首轮 cold）无 numerical failure；default warm median/p95 为 `4.510/5.056 ms`，lowT `3.928/4.267 ms`，highT `6.710/7.603 ms`，stiff `7.021/7.399 ms`，low_r `4.034/4.615 ms`，high_kappa `6.740/7.229 ms`。正式 default 仍未稳定低于 4 ms。
- 同一 native spectrum 的六种积分比较（Simpson、PCHIP、natural cubic、log-log PCHIP、Gauss-Legendre over PCHIP、Chebyshev）显示 Chebyshev 在 high-T/high-kappa 与 low-T 出现 `4%–8%` 级异常偏差，不可接受；natural cubic 与 log-log PCHIP 的偏差也随参数点改变方向，固定替换不可靠。
- DN-driven midpoint 原型从 76 nodes 依据局部 PCHIP-vs-trapezoid 差选择 4/10/20 个 interval。default 的 PCHIP DN 为 `0.0022636474 → 0.0022612282 → 0.0022625409 → 0.0022625338`（76/80/86/96 nodes），相对 fresh dense-reference `0.00226283297` 的误差先恶化后改善，未满足单调收敛，因此该 estimator 被拒绝进入生产。
- 新增证据脚本：`scripts/compare_frequency_quadrature.py`、`scripts/benchmark_goal_seeds.py`、`scripts/benchmark_dn_refinement.py`、`scripts/benchmark_head_matrix.py`；所有输出 JSON 均记录 commit、线程层、affinity 与代表点状态。
- default 同一 native 76-node grid 的独立 reference 对照（reference z_tail=8、DOP853 rtol=1e-9、4 workers）给出 reference PCHIP DN `0.0022643136483710326`；fast Simpson/PCHIP 分别为 `0.0022626969518100616`、`0.002263647415085622`，相对误差 `7.14e-4`、`2.94e-4`。因此 PCHIP 在该网格上更接近 reference，但仍未达到 `2e-4` 门槛；且此前 low-T 参数点差异达 `1.275e-2`，不能全局切换。fast 与 reference integrand 的相对误差 p95/max 为 `7.84e-4/3.57e-3`，当前主误差仍不只是频率积分器。
- `eval_freqs` invariant 已实现并在六个代表点验证：default/lowT/highT/stiff/low_r/high-kappa 从 76/77 个 support nodes 增加 4 个 native eval nodes 后，`DN_gw` 的绝对/相对变化均为 0，且无 failure；评估节点仍保留在 `m.f` 和频谱输出中，bolometric DN 使用独立 support grid。
- exact 76/80/90/110-node sweep（default，20 threads、workqueue）仍显示 blind seed 增加不具单调收敛性：相对 dense reference 的 DN 误差为 `6.01e-5/5.55e-4/1.40e-4/2.29e-4`，对应 76/80/90/110；因此当前目标网格尚不能宣称通过 `<2e-4` 的单调收敛门槛。
- fast solver 现在记录 `estimated_DN_quadrature_error` 与相对估计，采用同一 support grid 上 Simpson-PCHIP embedded pair；该估计只作为 telemetry，不改变默认 Simpson 或频率节点。它为后续 estimator coverage 提供了可观测字段。
- final HEAD `260b818` 的 25-repeat standardized matrix（20 threads/workqueue/CPU 0-19）在 telemetry 与 support/eval 分离后，default warm median/p95 为 `4.610/5.058 ms`，仍未达到 `<4 ms`；low-T/high-T/stiff/low-r/high-kappa warm median 为 `4.755/6.819/7.019/4.311/6.762 ms`。去除重复 support 权重构造的 A/B 未获得 >5% 稳定收益，已拒绝。
- 在 release code HEAD `9bb2697` 重新生成的 25-repeat 矩阵受机器瞬时波动影响，default warm median/p95 为 `5.361/5.734 ms`，其余五点为 `4.887/7.279/7.538/4.590/7.276 ms`；与前次同口径结果方向一致：没有稳定达到 `<4 ms`，因此不接受速度候选。该结果不改变 DN、failure 或 quadrature 结论。
- current code HEAD `1261bb1` 的 fast-only stability screen 覆盖 8 个命名点和 16 个 Sobol 点（24 点，20 threads/workqueue）。`failure_count=0`；3 个 Sobol 点触发既有 `shared_Neff_guard`，分类为 physical guard 而非 numerical failure。所有非 guard 点均产出有限 DN 与 embedded telemetry，但低振幅尾部的相对 estimator 可达 `0.20`，说明后续 coverage 必须同时报告绝对误差/信号掩码，不能只用相对 estimator 放行。
- fresh same-grid independent reference：low-T 的 Simpson/PCHIP DN 相对误差为 `1.241e-2/1.826e-4`，high-T 为 `1.033e-3/2.932e-4`；两点都显示 PCHIP 显著优于 Simpson，但 high-T 仍高于 `<2e-4`。这证明 estimator 需要参数空间 coverage，不能只依赖 default 或单一候选方法。
- 默认 telemetry 优化通过 profiler acceptance：固定 CPU 0-19、workqueue、20 threads、50 warm repeats，total median `4.976 -> 4.544 ms`，约 `8.7%` 改善；所有 spectrum/DN digest 不变。它只取消默认路径的重复 PCHIP object，改用 Simpson-trapezoid Richardson pair；PCHIP opt-in 仍计算直接 Simpson-PCHIP estimate。
- 优化后的六点 25-repeat matrix：default warm median/p95 `4.259/5.286 ms`，low-T `3.947/4.180 ms`，high-T `6.513/7.138 ms`，stiff `6.753/7.590 ms`，low-r `4.045/4.341 ms`，high-kappa `6.497/6.904 ms`。default 仍未稳定低于 4 ms，但已达到“有证据的 >5% 局部优化”验收。
- 采用廉价 Richardson estimator 后的 fresh 50-repeat fixed profiler：default warm total median `4.544 ms`，相对之前 `4.976 ms` 为 `8.7%` 改善；频谱、DN 与 digest 不变。24 点 stability artifact 已同步更新，仍为 `0` numerical failure、`3` 个明确 physical guard。
- estimator coverage fresh replay：六个 named same-grid independent references 中，廉价 Simpson-trapezoid Richardson 估计对 production Simpson DN 误差 `0/6` 覆盖；对 opt-in PCHIP DN 误差为 `6/6` 覆盖。因默认路径是 Simpson，coverage release gate 保持 `NOT VERIFIED`，该 telemetry 不能作为当前发布放行依据。
- `seed_n=78` 候选（实际 89 nodes）六工况 fast A/B 无 numerical failure，但独立 default same-grid oracle 显示 Simpson/PCHIP DN 相对误差为 `1.142e-3/2.952e-4`；因此“仅把 grid 加到约 89 nodes 并启用 PCHIP”仍未达到 `<2e-4`，候选拒绝升级。此前 76-node dense-reference 上的偶然 `9.3e-6` 不足以证明参数空间收益。
- 当前 HEAD outer-reuse safety replay（default/lowT/highT/stiff）无 numerical failure，DN rel 对 always-full 均不超过 `5.3e-10`；spectrum max rel 最大为 default `9.03e-4`、lowT `1.54e-5`，其余无 reuse 生效。按 DN `2e-4`、spectrum max `1e-3` 的保守门槛，false-safe `0/4`，暂保留当前 reuse；但 default 的 `9.03e-4` 已接近 spectrum budget，不能宣称 reuse 对所有 spectrum observable 完全无影响。
- stability screen 已从当前 HEAD `6a32838` 重跑并回写 SHA：24 点（8 named + 16 Sobol）`failure_count=0`、`guard_count=3`；3 个均为显式 `shared_Neff_guard`，没有 silent fallback。低振幅 Sobol 点的相对 estimator 仍可达 `3.08e-2`，所以 estimator 不能单独放行科学结果。
- 新增参数空间 same-grid independent oracle：`cr0_blue`、`positive_tilt`、`sobol_000`、`sobol_002`、`sobol_006` 均无 numerical failure；Simpson DN rel 为 `8.83e-4/7.14e-4/1.05e-3/9.78e-4/8.56e-4`，PCHIP 为 `2.76e-4/2.94e-4/2.96e-4/2.93e-4/2.98e-4`。PCHIP 改善一致但仍未满足 `<2e-4`，parameter-space release gate 保持 `NOT VERIFIED`。
- phase-cap A/B（`phase_max=0.35/0.5`，六工况各 25 次）没有稳定 >5% runtime 收益：相对 formal `.25` 的中位数多数持平或变慢；spectrum max 变化约 `5.2e-6..1.3e-5`、DN 变化约 `0.4e-6..6.8e-6`，数值扰动虽小但不能抵消速度失败，两个候选均拒绝升级。
- `solve_kernel` 的 Numba `literally(assemble)` 专门化实验被拒绝：focused regression `42 passed, 1 deselected`，但 50-repeat profiler 中 warm tensor kernel 约从 `1 ms` 退化到 `150 ms`，说明该 dispatch 方案不适合当前 parallel kernel；源码已完全回退。
- default same-grid oracle tail A/B：reference `z_tail=5` 给 `DN=0.0022718753`，`z_tail=8` 给 `0.0022643136`，相对差 `3.34e-3`；因此当前 fast PCHIP 的约 `2.94e-4` residual 是 tail/transfer 与 sparse quadrature 的合成，不能只归因于积分器。
- Prüfer full native-grid certification（2026-09-11）：完整 76 频率网格在
  `default/lowT/highT/stiff`、10 个 edge、5 个 Sobol 点上与 Cartesian DOP853
  的 full-grid `DN_gw` 相对差 median `6.44e-10`、max `2.44e-9`；分量 p95
  `Oj<=1.9e-6`、`Opgw<=1.1e-5`；max 异常均为未入尾低频模式的零点附近
  相对放大（绝对差 `1e-19..1e-22`）；outer 自洽 26 accepted 全部迭代一致、
  4 个物理 guard 双实现一致；default 全网格重放 bitwise 一致。PASS /
  VERIFIED（相对 Cartesian reference）。仍不晋升正式 kernel：速度证据在
  Python/DOP853 栈上（edge 含 outer 后 runtime ratio 最高 `0.929`），且
  tensor kernel 已是 Numba 并行；Prüfer 价值转为 Oracle B 风格独立状态
  变量交叉校验。artifacts 见 `docs/oracle_prufer_fullgrid_*.json`。
- Oracle C 解析高阶 WKB tail（2026-09-11）：从 Prüfer 精确振幅方程
  `d ln h/dN = cos(2θ)` 出发，stationary-phase 边界项给出解析修正
  `transfer² × (1 + sin(2θ_f)/ω_f)`。全 76 频率 native 网格、四点
  （default/lowT/highT/stiff）对比 frozen(z=5) / wkb(z=5) / deep(z=10)：
  frozen-vs-deep DN 相对差 `1.34e-3..3.65e-3`（正是 Stage B 记录的
  `3.63e-3` 非单调 systematic），wkb-vs-deep 降到 `5.46e-6..1.25e-5`
  （216-665 倍改善），per-mode 残差与 eps² 量级一致。结论：Stage B 的
  非单调是 frozen-amplitude 的 O(eps) 绝热缺陷而非物理；解析修正无需额外
  ODE 求解即可把 tail systematic 压两个数量级。artifacts 见
  `docs/oracle_c_wkb_*.json` 与 `docs/oracle_c_wkb_assessment.md`。

## Fast true DN_gw error calibration (2026-09-11)

- 动机：`fast` 与 frozen-amplitude reference 共享同一个 `z_tail` frozen 尾约定，
  因此 fast-vs-reference 比较会**抵消**共享 tail defect。Oracle C 给出解析 WKB
  修正后，fast 的**真实**误差可定义为同一 native 网格、同一 `DN_eff` 下 fast
  self-consistent `DN_gw` 与 WKB-corrected reference `DN_gw` 之差。新增
  `scripts/benchmark_fast_true_error.py`（reference-only 诊断，不改正式 kernel）；
  artifacts 为 `docs/fast_true_error.json` 与 `docs/fast_true_error_assessment.md`；
  资源 workers=1、Numba=2、BLAS=1。
- 结果（`fast DN_gw` / vs frozen(z5) / vs WKB(z5) / vs deep(z10)）：
  default `2.2626969518e-03` / `4.04e-3` / `4.31e-4` / `4.25e-4`；
  lowT `5.2446047710e-08` / `1.11e-2` / `1.24e-2` / `1.24e-2`；
  highT `5.6391339588e-02` / `4.34e-3` / `7.50e-4` / `7.62e-4`；
  stiff `1.4965149678e-02` / `4.91e-3` / `1.29e-3` / `1.28e-3`。
- default/highT/stiff：此前对外报告的同约定 PCHIP 数（manifest
  `same_grid_dn_rel_pchip` median `2.94e-4`）掩盖了共享 tail defect；真实误差
  `4.3e-4..1.3e-3`，比同约定数大 1.5-4 倍，仍未进入 `2e-4` release gate。
- default/highT/stiff 的 fast 值同时贴近 WKB(z5) 与 deep(z10)，而 frozen(z5)
  偏离 `~3.6e-3`。这与 `fast_sgwb.py` 的 `_tail_match_gamma`（`T ~ a^-1`
  振幅匹配）一致：fast 尾已比纯 frozen 更接近真值。
- lowT 是例外：fast 相对 frozen 与 WKB 两个 reference 都偏 `1.1e-2..1.2e-2`，
  远大于 `1.3e-3` 的 frozen->WKB tail 修正。其 `DN_gw=5.2e-8`，属非 tail 的
  独立误差源（quadrature / grid / 低振幅抵消），位置尚未定位。

### Confidence tables (fast true-error)

| Classification | Statement |
|---|---|
| VERIFIED | Oracle C WKB 修正把 frozen-vs-deep tail systematic 从 `1.34e-3..3.65e-3` 压到 `5.46e-6..1.25e-5`（216-665 倍），四点全 76 频率，且有独立 deep z=10 与 Prüfer 状态变量锚点。 |
| EMPIRICALLY VALIDATED | fast 真实 DN 误差（vs WKB 锚点）default/highT/stiff = `4.31e-4/7.50e-4/1.29e-3`；大于同约定 PCHIP `2.94e-4`，未进 `2e-4` gate。 |
| EMPIRICALLY VALIDATED | default/highT/stiff 的剩余误差由残余 fast-vs-WKB（非共享 tail）主导；fast 贴近 WKB/deep 而非 frozen。 |
| HEURISTIC | `_tail_match_gamma` 的 `T~a^-1` 匹配是 fast 尾比 frozen 更接近真值的原因（未独立证明）。 |
| UNVERIFIED | lowT 的 `1.24e-2` 非 tail 误差源（quadrature/grid/低振幅抵消）尚未定位。 |

决策：ACCEPTED as calibration finding。下一实验按频率分解 fast-vs-WKB 残差
（tail / deep-subhorizon stepping / frequency quadrature 三源），再决定是否在
fast tail assembly 内升格 Oracle C 修正。

## Fast DN_gw residual decomposition (2026-09-11)

- 新增 `scripts/benchmark_fast_residual_decomposition.py`；artifacts
  `docs/fast_residual_decomposition.json` 与
  `docs/fast_residual_decomposition_assessment.md`。自检：脚本用 `build_Wmat`
  权重复现的 `dn.fast.simpson` 与 fast 自报 `DN_gw` 逐位一致（default
  `2.2626969518e-03`），PCHIP 值也与既有 `0.002263647415085622` 一致。
- 关键结果（同一积分器下的 per-node 物理差 vs Simpson-vs-PCHIP 积分器差）：
  default `9.49e-06` vs `4.20e-04`；lowT `1.41e-04` vs `1.26e-02`；
  highT `9.32e-06` vs `7.40e-04`；stiff `1.10e-05` vs `1.28e-03`。
- 结论 1（EMPIRICALLY VALIDATED）：fast 的传播/tail kernel 精度约 1e-5；真实
  DN 误差由默认 Simpson 频率积分器主导，而非 tail 或 deep-subhorizon stepping。
- 结论 2（EMPIRICALLY VALIDATED）：default/highT/stiff 的残差 100% 落在 tail
  节点但仅 1e-5 量级，即 `_tail_match_gamma` 尾与 WKB 锚点之间的 O(eps^2) 残差。
- 结论 3（EMPIRICALLY VALIDATED）：lowT 残差的 |abs| 只有 11.0% 在 tail 节点，
  且最低频 `f=-18.4594` 单点占 52.7%；其 `DN_gw=5.2e-08`，属低振幅抵消/求积
  效应而非物理 tail。
- 推论（EMPIRICALLY VALIDATED，四个命名点）：默认换 PCHIP 后真实 DN 误差预期
  default `1.07e-05`、highT `9.88e-06`、stiff `1.15e-05`、lowT `1.66e-04`，
  全部进入 `2e-4` gate。这也解释了历史上的“节点加密非单调”与“PCHIP opt-in”
  现象：两者都在削减同一个 quadrature 项。
- 决策：ACCEPTED as diagnosis finding，不改正式 kernel。下一实验先写 acceptance
  criteria 再动手：默认 PCHIP 化（必要时 Numba），要求 DN rel `<2e-4`、
  spectrum max 不退化、warm median 增幅 `<10%`、no new failure、determinism pass。

### Confidence tables (residual decomposition)

| Classification | Statement |
|---|---|
| EMPIRICALLY VALIDATED | fast 传播 kernel 的 per-node DN 精度为 `9.3e-6..1.4e-4`（同积分器 vs WKB 锚点），比 Simpson-vs-PCHIP 差小 1-2 个数量级。 |
| EMPIRICALLY VALIDATED | default/highT/stiff 残差 100% 在 tail 节点但仅 1e-5 级；lowT 残差 89% 在非 tail，最低频单点占 52.7%。 |
| EMPIRICALLY VALIDATED | 四个命名点上默认改用 PCHIP 后真实 DN 误差 `1.07e-5..1.66e-4`，进入 `2e-4` gate。 |
| HEURISTIC | PCHIP 化的 runtime 成本（scipy `PchipInterpolator` + `integrate`）尚未测量；Numba 化可行性未验证。 |
| UNVERIFIED | PCHIP 默认在其他 Sobol/edge 参数点是否同样接近 WKB 锚点（目前仅 4 个命名点）。 |

## Fast frequency-quadrature A/B: Simpson vs PCHIP (2026-09-11)

- 新增 `scripts/benchmark_fast_quadrature_ab.py`；artifacts
  `docs/fast_quadrature_ab.json` 与 `docs/fast_quadrature_ab_assessment.md`。
  交替测量（每个 repeat 内两种积分器互换顺序），warmup=3、repeats=25、
  Numba=2、BLAS=1、升序网格。
- 精度（vs Oracle C WKB 锚点）：PCHIP 把 default `4.31e-04 -> 1.07e-05`、
  highT `7.50e-04 -> 9.88e-06`、stiff `1.29e-03 -> 1.15e-05`、
  lowT `1.24e-02 -> 1.66e-04`（1-2 个数量级），四点全部进入 `2e-04` gate。
- 速度：scipy PCHIP 路径 warm median 比值 `1.307/1.300/1.204/1.278`
  （default/lowT/highT/stiff），超过预先登记的 `<10%` 预算 → 直接切换默认被
  拒绝（第十二原则要求 accuracy win 的 runtime 增加 <10%）。
- 成本定位（微基准）：`PchipInterpolator(...).integrate()` `0.141 ms/次`
  （每次求解约 3 次：两个外层 `g2_last` 加 `g2c[-1]`），
   `estimate_frequency_quadrature_local(..., 'pchip')` `1.39 ms/次`，合计与
   观测 `+2.1 ms` 吻合。成本完全来自 scipy 实现。
- 更正（2026-09-11，见下节）：本节最初推测“PCHIP 积分对节点值是线性泛函，
   可预计算固定全局权重向量与逐区间权重矩阵”。该推测已被否证：
   Fritsch-Carlson 斜率是节点值的非线性函数。
- 决策：ACCEPTED as measurement，不改任何默认。下一实验：消除 PCHIP 路径中
   重复的 scipy 拟合工作，使 PCHIP 默认化落在 runtime 预算内。

## Fast PCHIP quadrature: single-fit sharing (2026-09-11)

- 否证：用 `PchipInterpolator(x, np.eye(n), axis=0)` 构造的“权重”与逐列区间
  积分一致，但 `weights @ y` 与 `PchipInterpolator(x, y).integrate()` 相差约
  20 倍量级 —— 不存在固定权重向量，因为 PCHIP 的 Fritsch-Carlson 斜率是节点值
  的非线性函数。原“预计算权重”路线作废。
- 替代实现：`pchip_integral_breakdown(freqs, integrand)` 一次构造
  `PchipInterpolator` 与其 `antiderivative()`，同时返回全程积分与逐区间积分；
  `estimate_frequency_quadrature_local(..., candidate_intervals=...)` 接受这些
  逐区间积分而不重建 spline；`_SGWB_iter_fast_impl` 每个外层迭代只做一次
  分解，收敛判定 `g2_last`、循环后 `g2c[-1]` 与局部 estimator 共享同一次拟合。
- 同时把局部 estimator 的非重叠三点 Simpson 基线从“逐面板 `scipy.simpson`”
  改为等价向量化公式：实测与逐面板结果**逐位一致**（max abs diff `0.0`），
  成本 `0.674 ms -> 0.027 ms`（76 点网格，median of 300）。
- 同会话配对 A/B（`scripts/benchmark_fast_quadrature_ab.py --repeats 50`，
  前后各一次、NUMBA=2、BLAS=1，`docs/fast_quadrature_reuse_ab{,_before}.json`）：
  PCHIP 专属开销 default `2.20 -> 0.98 ms`、lowT `1.92 -> 0.59 ms`、
  highT `2.01 -> 0.72 ms`、stiff `2.29 -> 0.84 ms`（约 55-70% 更少），
  PCHIP/Simpson warm median 比值 `1.366/1.307/1.228/1.240 ->
  1.165/1.095/1.073/1.089`。
- 数值不变性：四点 Simpson 与 PCHIP 的 `DN_gw` 与改前**逐位一致**（rel `0.0`），
  vs Oracle C WKB 锚点不变（PCHIP `1.07e-05/1.66e-04/9.88e-06/1.15e-05`）。
- 生产刻度探针（`NUMBA_NUM_THREADS=16`）：并行部分缩小后串行开销占比更高，
  比值 `1.154/1.134/1.137/1.101`（default/lowT/highT/stiff）。
- 决策：ACCEPTED（opt-in PCHIP 路径严格 Pareto 改善：同等观测值、约 60% 更少
  PCHIP 专属开销，符合第十二原则 Speed win）；**REJECTED：本阶段不切换默认
  `frequency_quadrature`**，因为预登记的 `<10%` 未稳健满足（同会话 default
  `1.165`；此前会话 `1.093`/`1.123`）。详见
  `docs/fast_quadrature_reuse_assessment.md`。
- 剩余 PCHIP 专属成本：scipy `PchipInterpolator` 构造本身（`0.15 ms` x 2 次/求解）
  加 estimator 分摊循环（约 `0.26 ms`），均为实现成本而非数学成本。

### Confidence tables (PCHIP single-fit sharing)

| Classification | Statement |
|---|---|
| VERIFIED | 不存在能复现 PCHIP 全区间积分的固定节点权重向量（单位矩阵探针 + 直接反例）。 |
| VERIFIED | `pchip_integral_breakdown` 的全程积分与逐区间积分和 scipy 同拟合结果一致（`tests/test_pchip_integral_breakdown.py`，rel `1e-12`）。 |
| VERIFIED | 向量化三点 Simpson 基线与逐面板 `scipy.simpson` 逐位一致（实测 max abs diff `0.0`）。 |
| EMPIRICALLY VALIDATED | 四点 `DN_gw` 与改前 scipy PCHIP 路径逐位一致（rel `0.0`），vs WKB 不变。 |
| EMPIRICALLY VALIDATED | PCHIP 专属开销下降 55-70%（同会话配对、50 repeats）；比值 `1.07-1.17`（2 线程）、`1.10-1.15`（16 线程）。 |
| HEURISTIC | 用 NumPy 向量化 PCHIP 斜率 + 分段解析积分替代 scipy 可再省约 `0.3 ms/求解`；未实测。 |
| UNVERIFIED | `pchip` 作为默认后在 Sobol/edge 参数空间与 16+ 线程生产刻度的 runtime 与精度表现。 |

下一实验（acceptance criteria 先写）：NumPy 向量化 PCHIP 斜率 + 分段解析积分核
替代热路径 scipy 拟合，并向量化 estimator 分摊循环；要求与 scipy PCHIP
`<1e-12` rel、2 线程与 16 线程 PCHIP/Simpson 比值均 `<1.10`、DN 不变、
no new failure、determinism pass，再评估默认切换。

## Vectorized PCHIP kernel and estimator allocation (2026-09-11)

- 实现：`_pchip_integrals_vectorized` 用 NumPy 复刻 scipy 的 Fritsch-Carlson
  斜率（内部节点加权调和平均 + 端点 `_edge_case` 形状保持 guard），分段用
  Hermite 三次的闭式积分 `h/2*(y_i+y_{i+1}) + h^2*(m_i-m_{i+1})/12`；
  `pchip_integral_breakdown` 保持 scipy 参考实现，热路径改用向量化核；
  `estimate_frequency_quadrature_local` 的逐面板分摊向量化，零候选和面板
  保持“半误差均分”的旧语义。
- 正确性：向量化分摊与旧逐面板循环**逐位一致**（三种 allocation，测试断言）；
  向量化 PCHIP 核与 scipy 参考逐区间 `<=1e-9` rel（绝对 `1e-12*max|区间|`）、
  总和 `<=5e-14` rel（200 组随机网格 + 平坦/符号翻转/两点退化样本）；
  端到端 `test_pchip_frequency_quadrature_is_opt_in`（把 `m.DN_gw[-1]` 钉在
  scipy `integrate_frequency_pchip`，rel `1e-12`）仍通过。
- 成本：微基准（76 点网格、median of 2000）scipy breakdown `0.143 ms` ->
  向量化 `0.024 ms`；estimator（复用候选）`0.290 ms -> 0.061 ms`。
- 配对 A/B（`--repeats 50`，同会话前后配对，2 线程与 16 线程各一次）：
  2 线程 PCHIP/Simpson 比值 `1.1622/1.1261/1.0495/1.0901 ->
  1.0105/1.0551/1.0397/1.0411`；16 线程 `1.1987/1.1137/1.1270/1.1203 ->
  1.0428/1.0342/1.0180/1.0413`（default/lowT/highT/stiff）。PCHIP 专属开销
  `0.06-0.42 ms`（2T）/`0.15-0.33 ms`（16T）。同会话第二次配对复现
  （before `1.050-1.184`、after `1.002-1.056`）。
- 数值不变性：实测 `DN_gw` 相对 scipy PCHIP 路径变化
  `3.83e-16/4.98e-16/0/0`（1-2 ulp）；默认 Simpson 路径逐位不变。
- 决策：ACCEPTED。预登记 `<1.10` 在 2 线程与 16 线程两个刻度均满足；默认
  `frequency_quadrature` 仍不切换（切换需同步刷新 manifest/README/
  ERROR_BUDGET/coverage 并重跑参数空间门禁，作为独立 Phase）。

### Confidence tables (vectorized PCHIP kernel)

| Classification | Statement |
|---|---|
| VERIFIED | 向量化 PCHIP 核与 scipy `PchipInterpolator` 参考在随机与退化样本上一致（逐区间 `<=1e-9` rel、总和 `<=5e-14` rel，测试断言）。 |
| VERIFIED | 向量化分摊与旧逐面板循环逐位一致（三种 allocation 的 `np.array_equal` 断言）。 |
| VERIFIED | 默认 Simpson 路径逐位不变；PCHIP DN 仅变化 1-2 ulp（`3.83e-16/4.98e-16/0/0` rel）。 |
| EMPIRICALLY VALIDATED | 2 线程与 16 线程的 PCHIP/Simpson 比值在两次同会话配对中均 `<1.10`（1.00-1.056）。 |
| HEURISTIC | 该比值在其他机器/负载下同样 `<1.10`（本机 wall time 跨运行可漂移 ~60%，故只报 run-internal 比值）。 |
| UNVERIFIED | PCHIP 作为默认后在 Sobol/edge 参数空间的 runtime 与精度（下一 Phase 的 acceptance 已写定）。 |

### Confidence tables (quadrature A/B)

| Classification | Statement |
|---|---|
| EMPIRICALLY VALIDATED | PCHIP 默认化把四个命名点的真实 DN 误差降到 `1.07e-5..1.66e-4`，全部 <2e-4 gate。 |
| EMPIRICALLY VALIDATED | scipy PCHIP 路径的 warm runtime 比值 `1.204..1.307`，不满足 `<10%` 预算。 |
| EMPIRICALLY VALIDATED | 成本来自 scipy 构造/积分（`0.141 ms/次`）与局部 estimator（`1.39 ms/次`），与 `+2.1 ms` 观测一致。 |
| HEURISTIC | 预计算权重点积可把该开销降到 <2%（按微基准外推），尚未实现验证。 |
| UNVERIFIED | PCHIP 在更广 Sobol/edge 参数空间的精度是否同样 <2e-4。 |

### Confidence tables (PCHIP 默认切换，2026-09-11)

| Classification | Statement |
|---|---|
| VERIFIED | `SGWB_iter_fast` 默认 `frequency_quadrature='pchip'`；`simpson` 仍可显式选择，默认 DN 与 scipy `integrate_frequency_pchip` 相对差 `1.8e-14`。 |
| VERIFIED | PCHIP 热路径的非有限 integrand 走统一 `nonfinite` guard（abort/restore），与 Simpson 语义一致；`_pchip_integrals_vectorized` 仍 fail-loud。 |
| EMPIRICALLY VALIDATED | 四个命名点 vs Oracle C WKB 残差 `1.07e-5/1.66e-4/9.88e-6/1.15e-5`（default/lowT/highT/stiff），全部 `<2e-4`。 |
| EMPIRICALLY VALIDATED | 默认切换的配对 warm runtime 比值在 2/16 线程均 `<1.10`（2T `1.019/1.035/1.034/1.024`、16T `1.069/1.028/0.987/1.017`）。 |
| EMPIRICALLY VALIDATED | 20 线程 25-repeat 矩阵 default warm median/p95 `4.932/5.467 ms`、cold `0.222 s`；highT/stiff/high-kappa `7.79/7.79/7.27 ms` 仍是超 `<4 ms` 的热点。 |
| HEURISTIC | 同网格 reference 残差 `2.93e-4`（Sobol/edge `2.76e-4..2.98e-4`）被解释为“尾部/传递 + 积分”合并残差；尚未用 nested 真实频率求积分离这两项。 |
| UNVERIFIED | 更广 Sobol/edge 参数空间相对 Oracle C WKB 的积分残差是否同样 `<2e-4`（目前只有 4 个命名点）。 |

上一张表的 legacy 行（scipy PCHIP 路径 `1.204..1.307`、预计算权重点积假设）已被后续
向量化 PCHIP 核与单次拟合共享取代，见 `docs/fast_quadrature_reuse_assessment.md`
和 `docs/fast_quadrature_default_switch_assessment.md`。

## Fast Python preparation-layer de-duplication (2026-09-12)

20 线程阶段分解显示 warm runtime 并非由 kernel 主导：default 点
`tensor_solve_kernel` 约 `0.81 ms`，而 `expansion_background` `1.04 ms`、
`fast_phi_s2_split` `0.63 ms`、`goal_frequency_construction` `0.48 ms`，
其中包含纯冗余工作。本轮只消除冗余，不触碰任何数值路径：

- `gen_fast` 的目标频率网格原先对同一数组做全量重排
  （`np.unique(np.sort(...))`，nv=13764 时约 `158 us`），改为 `concatenate` +
  稳定排序 + `not_equal` 归并去重（约 `25 us`），并缓存 `n_re_abs`；
- `grid_independent_freqs.f_hor_cont` 中与 `N` 无关的
  `H2_vec(N_inf)`/`H2_vec(N_re_abs)`/`raw_last`/`raw_re`/`Delta_f`/`ln10`
  提到闭包外只算一次；
- 新分配的 `Ogw`/`Oj`/`Opgw` 缓冲区跳过紧随其后的 `fill(0.0)`
  （`_fresh_buffers` 标记）；
- `_sigma_node_limits` 与 `fast_phi_s2_split` 中重复的 `m.derived_param`
  属性求值绑定为单次 `d`（该 property 每次访问都重算，约 `11 us/次`）。

验证工具 `scripts/benchmark_fast_runtime_ab.py` 输出六点 warm
median/p95/min 与 `f`/`log10OmegaGW`/`DN_gw`/`g2`/`w2` 的 SHA256 digest。

### Confidence tables (preparation-layer de-duplication)

| Classification | Statement |
|---|---|
| VERIFIED | 六个 probe regime 的五个输出字段 SHA256 digest 在改动前后逐位一致（20 线程 648 项、2 线程 288 项比较，另含 A/B 交叉比较），`converged`/`n_freq`/`fast_failure_reason` 不变。 |
| EMPIRICALLY VALIDATED | A,B,B,A 对称序配对（固定 affinity）显示 default warm median 改善：20 线程 `5.4%`（ratio `0.9456`，单轮 `0.802-1.029`）、2 线程 `8.2%`（ratio `0.9177`，单轮 `0.904-0.932`）。 |
| EMPIRICALLY VALIDATED | 六点中位 ratio 在 20 线程为 `0.816-0.987`、2 线程为 `0.872-1.011`；无任何点退化超过 `2%` 噪声门限。 |
| HEURISTIC | 对称序（A,B,B,A）抵消移动 CPU 频率漂移；单轮 ratio 极值仍受热降频影响，应以中位为准。 |
| UNVERIFIED | `fast_phi_s2_split` 内剩余约 `0.4 ms` 的组成，以及为 `derived_param` 加显式缓存的安全性。 |

## Nested native-frequency quadrature (2026-09-12)

第九原则要求用**真实新增频率求解**而不是在同一 interpolant 上换算法来估计
频率积分误差。新增 `scripts/benchmark_nested_frequency_quadrature.py`：
base goal 网格 -> 取局部估计器 top-N 区间 -> 在每个区间中点并入**积分**
support grid（受控替换 `goal_oriented_freqs`，`finally` 恢复）-> 重新求解
-> `E_nested = |DN_refined - DN_base|`，再与同网格独立 reference 的
`actual_error` 比较，并统计 coverage / false-safe / monotonicity。

实现细节：`eval_freqs` 只把节点加入 solve 网格、不改变积分 support grid，
所以第一版实现得到 `E_nested = 0`（新节点不参与求积）；改为替换 support grid
后才有信号。阳性对照用 `simpson`（已知积分残差大）验证脚本灵敏度。
证据 artifact：`docs/nested_frequency_quadrature_head.json`（默认 PCHIP，六点、
N=4/8/12、`NUMBA_NUM_THREADS=2`/`FAST_THREADS=2`）与
`docs/nested_frequency_quadrature_simpson_control.json`（default 点阳性对照）。

### Confidence tables (nested frequency quadrature)

| Classification | Statement |
|---|---|
| VERIFIED | 六点（default/lowT/highT/stiff/low_r/high_kappa）在 top-4/8/12 区间各插入真实中点后，`E_nested` 绝对量级 `1.7e-14`–`6.3e-10`、相对 `7.4e-12`–`2.8e-9`（lowT 因 `DN_base ~5e-8` 相对放大到 `9.3e-5`，是唯一相对值超 `1e-8` 的点）。 |
| VERIFIED | 六点 `actual_rel`（fast vs 同网格独立 reference）为 `1.83e-4`–`2.97e-4`，比 `E_nested_rel` 大 5–7 个量级 ⇒ 当前 DN 残差**不由频率网格离散/插值主导**，而由求解器传输（fast 固定步长 vs reference 连续 sigma DOP853）与 tail 约定主导。这独立解释了历史「blind node count increase / 89-node」为何无效。 |
| VERIFIED | 阳性对照：同一加密流程在 `simpson` 下给出 `E_nested = 1.36e-9`（比 PCHIP 大 5 个量级），说明脚本对积分方法敏感、不是无效测量。 |
| EMPIRICALLY VALIDATED | 默认 PCHIP 的保守局部估计器 `predicted_rel` 在六点为 `3.1e-3`–`1.5e-2`，是 `actual_rel` 的 10–80 倍 ⇒ 覆盖成立；`simpson` 下 `predicted_rel = 3.4e-4 < actual_rel = 7.1e-4` ⇒ 该估计器在 Simpson 上 false-safe。 |
| HEURISTIC | `E_nested` 自身在六点上都是 false-safe 的（`coverage = actual_error/E_nested` 为 `2.0`–`3.95e7`，lowT N=12 已降到 `2.0`），只能作为「网格已收敛」的敏感性检查，不能当作误差上界。 |
| UNVERIFIED | 单调性：default/lowT/low_r/stiff 严格单调不减；highT `7.33e-12 -> 7.32e-12 -> 1.47e-11` 与 high_kappa `6.32e-10 -> 6.32e-10 -> 4.44e-10` 被严格判据标记为不满足（前者是 `1e-14` 量级舍入抖动，后者向收敛值回摆），因此单调性未作硬门限。 |

### Oracle A same-grid `z_tail` attribution (2026-09-12)

`E_nested` 只能说明 fast 自身对频率网格已收敛；要判断 `actual_rel ~2.9e-4`
（fast vs 同网格连续 sigma DOP853）来自哪一侧，必须移动 **oracle** 的 tail
约定而不是 fast 的。用 `scripts/benchmark_same_grid_reference.py --z-tail`
把 Oracle A 的 frozen handoff 从 `z=8` 加深到 `z=10`（default 点、同网格、
`--rtol 1e-9`、2 线程）：

| Pipeline | z_tail | DN_gw | vs fast PCHIP rel | reference 自报 handoff 缺陷 `|1.5*sigma-1|/e^z` | 参考运行时间 |
|---|---:|---:|---:|---:|---:|
| Oracle A（连续 sigma DOP853） | 8 | `2.2643136484e-3` | `2.944e-4` | ~`3.4e-4` | `94.2 s` |
| Oracle A（同上） | 10 | `2.2636586329e-3` | `4.956e-6` | ~`4.5e-5` | `578.8 s` |
| Oracle C（Prüfer amplitude-phase） | 10 | `2.2636592187e-3` | `5.214e-6` | n/a | — |
| fast 正式路径 | 5 | `2.2636474151e-3` | — | — | — |

结论：`2.944e-4` 与 Oracle A 自己的 frozen-handoff 缺陷同阶，加深 oracle 后
降到 `4.96e-6`（59 倍），两个独立 oracle 在 z=10 上互差 `2.4e-7`。因此
default 点的 DN_gw 残差由 **oracle 的 frozen-tail 约定**主导，而不是 fast 的
频率积分；fast 对 tail 收敛 oracle 的真实 DN 误差约 `5e-6`，比 `5e-4` 预算低
两个量级。反面结论：oracle 加深不能作为逐点验证手段（z=14 在 `>20 min` 后
仍未收敛被终止），且现有 error-budget 模型按 z=8 锚定的 frozen-tail 项
（`2e-5`）在 `z_tail=5` 上过于保守。

### Confidence tables (Oracle A z_tail attribution)

| Classification | Statement |
|---|---|
| VERIFIED | default 点、同网格、同一 DN_eff：Oracle A 的 `z_tail` 由 8 加深到 10 后，fast-vs-oracle DN 相对差由 `2.944e-4` 降到 `4.956e-6`，降幅与 reference 自报的 frozen-handoff 缺陷 `|1.5*sigma-1|/e^z`（`3.4e-4 -> 4.5e-5`）同阶。 |
| EMPIRICALLY VALIDATED | 两个独立 oracle（连续 sigma DOP853 与 Prüfer amplitude-phase）在 z=10 上互差 `2.4e-7`；fast 距两者 `4.96e-6`/`5.21e-6`。 |
| HEURISTIC | 加深 oracle 的设备成本：z=8 `94.2 s`、z=10 `578.8 s`（6.1 倍）、z=14 在 `>20 min` 后未收敛（更多模式不再触发解析 handoff，DOP853 需积到 today），因此只用于定点归属。 |
| UNVERIFIED | 该归属只在 default 点、且 `DN_eff` 冻结为 fast 自洽值时验证；其余五个正式代表点是否同样由 oracle tail 项主导尚未逐点验证。 |

## Technical Decisions

- 有限 phase-window Oracle B 原型在 default/low-T/high-T/stiff 各 3 个可入尾模式上显示 z=5 到 z=7 的 phase-averaged today observable 变化为 `5.321e-3/5.552e-3/3.927e-3/6.897e-3`；低频未入尾部显式标记。该原型仍复用 DOP853 张量方程和一阶解析尾部，只能作为 handoff sensitivity 证据，不能晋升独立 oracle。
- Prüfer amplitude-phase standalone prototype 在 default/low-T/high-T/stiff 四点各 14 个可比较模式上相对 Cartesian DOP853 的最大 amplitude/power 误差为 `1.19e-7/2.38e-7`，最大 phase 差 `2.51e-5 rad`，runtime ratio `0.529–0.545`。接受为独立状态变量的研究原型，不接入正式 fast；仍需 full `DN_gw`、reheating edge、guard 和 determinism 复核。
- Prüfer 固定 8 频率完整 today 输出复核已完成：四点、z=5/7 的 `DN_gw` 最大相对差 `2.26e-9`，`Ogw/Oj/Opgw` 最大分量差 `1.67e-7`；default 重复运行的数值字段 bitwise 一致。该证据仍固定 `DN_eff`，未覆盖正式 outer self-consistency，因此不晋升正式 kernel。
- Prüfer outer self-consistency 复核已完成：default/low-T/high-T/stiff 四点、z=5/7 为 1–2 次迭代收敛，且 Prüfer/Cartesian 迭代次数一致；`DN_gw` 最大相对差 `3.44e-9`，无 `DN_eff` guard crossing。仍只覆盖固定 8 频率，reheating edge/Sobol/full-grid 尚未认证。
- Prüfer reheating/Sobol follow-up 已完成：`edge_tre_lo` 与 `sobol_000/002/006` 均在 z=5/7 收敛且 outer 迭代次数一致，最大 outer `DN_gw` 相对差分别为 `8.33e-10/2.12e-9/2.06e-9/2.50e-9`；`edge_tre_hi` 两个 z_tail 均由 Prüfer 与 Cartesian 一致触发显式物理 guard。新增 guard 分类后，物理拒绝不再被误报为数值失败；原型仍不进入正式 kernel。
- Prüfer complete edge/Sobol replay 已完成：十个参数轴 edge 加 `sobol_000/002/006/010/015` 共 15 点、z=5/7 产生 26 个 accepted outer comparisons 与 4 个显式 physical-guard records，无 numerical failure；accepted 点 outer 迭代计数全部一致，最大 outer `DN_gw` 差 `6.44e-9`、最大 handoff power 差 `3.26e-7`、最大 runtime ratio `0.929`。full production frequency-grid comparison 仍未完成，不能晋升正式 kernel。

| Decision | Rationale |
|----------|-----------|
| 先定位真实数据流再设计 split | 目标要求精确计算 `N_re` 且不得让 transfer step 跨 kink；不能凭文件名猜测已有实现 |
| 设计阶段必须区分“接口合并”和“数值算法变更” | 这样可以保持物理契约，且每一项算法收益可独立归因 |
| 本轮先在远端 `fast_v0.2` 当前 HEAD 做固定环境 fresh profiling | 用户明确要求不得用旧 profiler 推断当前热点；未有证据前不做 micro-optimize |
| DN_gw 误差先做组件隔离 A/B，再决定 frequency quadrature 或 primitive 改动 | spectrum 已达到较高精度，节点加密曾使 DN 变差，必须先量化误差贡献 |

## Issues Encountered

| Issue | Resolution |
|-------|------------|
| `create_goal` 返回已有未完成 goal | 使用当前 active goal；不创建第二个 goal |
| 远端 `fast` 尚未在本地分支跟踪信息中显示 | 已确认 remote 名称和 URL；下一步核对远端同名 ref |

## Resources

- `C:\Users\联想\.codex\attachments\79e36cae-8136-418c-9325-fde68b9fc4ff\goal-objective.md`
- `F:\codex\stiffGWpy\stiffgwpy_fast\fast_sgwb.py`
- `F:\codex\stiffGWpy\stiffgwpy_fast\stiff_SGWB.py`
- `F:\codex\stiffGWpy\stiffgwpy_fast\reference.py`

## Fresh current-HEAD phase-kernel feasibility spike (2026-09-16)

- 重新 fetch 后确认本地 `codex/fast_v0.2`、远端 `fast/fast_v0.2` 和当前源码
  HEAD 均为 `a1d701a5c9da8aa3a2132088e12c09fd011a8151`。
- 现有 standalone Prüfer 在完整 76/77 native goal nodes、z_tail=5、
  `NUMBA=2`、BLAS=1、reference workers=1 下重新对照 Cartesian DOP853：
  `DN_gw` 相对差 default/high-T/stiff/high-kappa 为
  `3.02e-10/2.53e-10/2.66e-10/4.66e-10`；outer 自洽差均约 `2–5e-10`，
  无 numerical failure。这是当前 SHA 的独立 Oracle 认证，不复用旧 SHA 的数值结论。
- 可删除的 Numba Prüfer propagation spike 以相同 `j0/z0/Phi/S2` 数据流测试：
  固定步 RK4 的纯 propagation warm median 为
  `2.41/2.31/2.41/2.33 ms`（default/high-T/stiff/high-kappa，25 repeats）；
  相对同 handoff DOP853 的 amplitude 最大误差为
  `1.74%/0.97%/1.48%/2.37%`，phase 最大误差为
  `0.0136/0.0181/0.0082/0.0103 rad`。相对同一 production Cartesian step
  的 amplitude 仍有 `0.74%–1.48%` 偏差，未达到生产精度或复杂度 gate。
- 决策：固定步 midpoint/RK4 phase kernel 均 `REJECTED`，不接入 production API；
  失败原因是当前 fixed-step Cartesian exact transfer 与 phase ODE 离散化的
  amplitude 误差仍为百分比量级，且 kernel 没有稳定速度收益。后续若继续 phase
  线，必须研究保持 exact transfer map 的 phase representation，而不是重复调 RK 阶数。
- 该 full-grid 对照给出的是 `z_tail=5` 下 Prüfer/Cartesian 的离散化交叉不确定度；
  tail 的绝对系统误差仍未被这四个点的 pairwise agreement 证明，不能把
  `2–5e-10` 直接当作完整 production error budget。

## Fresh current-HEAD native tail A/B and production handoff test (2026-09-16)

- 在同一当前 SHA、完整 76/77 native goal nodes、固定 affinity、Numba=20、
  workqueue、BLAS=1、workers=1 下，重新比较 production `z_tail=5/6/7`；三份
  benchmark artifact 为 `docs/benchmark_head_matrix_ztail5_round_20260916.json`、
  `docs/benchmark_head_matrix_ztail6_round_20260916.json`、
  `docs/benchmark_head_matrix_ztail7_round_20260916.json`。
- production `DN_gw` 从 z=5 到 z=7 的相对变化仅约
  `4.1e-6/1.0e-4/5.1e-6/5.2e-6/4.5e-6`
  （default/low-T/high-T/stiff/high-kappa），而 warm median 升至约
  `5.74/5.20/8.22/8.59/8.68 ms`（z=7，按同顺序）。
- 与当前 SHA Prüfer/Cartesian z=7 deep-tail outer oracle 的 production DN 相对差为
  `7.56e-4/1.03e-4/7.55e-4/7.60e-4/7.54e-4`；这表明主要差异不是 production
  handoff 深度，而是 propagation/grid/solver model 对齐误差。
- 决策：单纯把 production handoff 从 z=5 推到 z=6/7 为 `REJECTED`；精度收益
  不足且速度退化，不能进入 production。下一轮应攻击 propagation/grid/model
  alignment，并用 current-SHA deep oracle 重新验证。
- current-SHA 7-repeat breakdown（Numba=2、workqueue、BLAS=1）显示 median
  `tensor_solve_kernel` 为 `4.36–4.75 ms`，约占 `8.75–10.21 ms` profile total；
  `kernel_prepare` 约 `1.08–1.39 ms`，expansion 约 `0.66–0.77 ms`，而 PCHIP
  仅约 `0.12–0.13 ms`。速度线最高价值目标因此是减少真实 tensor solve/传播复杂度，
  不是继续优化 quadrature。
- 为避免 profiling artifact 与源码脱钩，`scripts/profile_fast_breakdown.py` 已让 JSON
  写入 `commit`；在新工具 SHA `dce33d0ce874bafc36693dc873ec27c222cd1f54` 重跑四点，
  tensor-solve median 为 default/high-T/stiff/high-kappa `4.93/4.72/4.84/4.17 ms`。
- 复核发现上述 profile 使用了 `kink_split=false`，不等同正式 fast benchmark；已在
  release HEAD `abd3d042600a1cfd507ffdf8949eaeeb8c5f60ab` 以 `kink_split=true` 重跑，
  正式路径 total/tensor median 为 default/high-T/stiff/high-kappa
  `6.34/2.69/0.57/0.73`、`10.34/4.70/0.73/0.73`、`10.39/5.20/0.67/0.67`、
  `10.21/5.08/0.71/0.71 ms`（分别为 total/tensor/prepare/expansion；prepare 在
  exact-kink 路径由独立计时点不产生可比调用）。旧的 false-kink 结论作废。
## Fresh current-HEAD step-size Pareto and parameter-space audit (2026-09-16)

- 在正式 `kink_split=true`、goal native grid、20 threads、workqueue、BLAS=1、25
  repeats 下测试 h=`0.005/0.00625/0.0075/0.01`；所有候选均无 numerical failure。
- h=.0075 相对 h=.005 的 warm median 约快 `34%/26%/22%/15%`
  （default/high-T/stiff/high-kappa），但 deep Oracle C 相对误差由
  `5.21e-6/2.24e-5/3.35e-6/2.45e-5` 变为
  `8.01e-6/2.59e-5/6.61e-6/2.75e-5`；多数 regime 精度回退。
- h=.00625 的 warm median 为 `4.47/5.90/6.17/6.21 ms`，未形成四点稳定 >5%
  收益，且 deep 误差也未改善。h=.01 的 WKB/deep 误差进一步升高。
- matched continuous-sigma 10 点参数验证（default、low-T、high-T、stiff、
  rad-dominant、tiny-r、transition、cr0-blue、extreme、Sobol）显示 h=.0075
  的 signal max relative `6.7833e-3`、DN max `3.6052e-3`；h=.005 为
  `6.7832e-3`、`3.6014e-3`，没有可证明的精度收益。
- 当前 quadrature estimator 9 点诊断为 coverage `1.0`、false-safe `0`；但输入
  reference provenance 不完整，仍不能把它晋升为新的独立 estimator certification。
- 决策：h=.0075、h=.00625、h=.01 均 `REJECTED FOR PRODUCTION`；不为速度牺牲
  已有精度。下一步回到 exact transfer-map/减少 tensor solves，并修复 fast 与
  continuous-sigma 的 model-alignment 误差。

## Linear-z Magnus transfer-map phase prototype (2026-09-16)

### Hypothesis and scope

在每个 native step 内把变换后的无迹传播矩阵 `A(z)` 按线性 `z(N)` 处理，使用
二阶 Magnus transfer map，可能比当前 constant-midpoint map 减少传播/模型对齐误差，
且不增加真实 tensor solve。原型保持 standalone/reference-only，验收标准见
`docs/transfer_map_phase_assessment.md`。

### Current-SHA evidence

Artifact `docs/transfer_map_phase_round.json` 绑定 commit `d6aef07`，固定 Numba=2、
workqueue、BLAS=1、reference workers=1，覆盖 default、low-T、high-T、stiff、
high-kappa，每点 8 个 native representative modes、z_tail=5、Cartesian DOP853
rtol=`1e-9`，并对每点进行 25 次 warm propagation A/B：

- 35 个进入尾部的模式中，candidate 相对当前 midpoint 的最大振幅变化仅
  `6.80e-6`，最大相位变化 `4.23e-4 rad`；31/35 个模式的振幅方向略有改善，
  但没有达到预注册的 10x handoff 精度改善，更不能推断 full-grid DN 改善。
- Python prototype propagation median (baseline -> candidate) 为
  default `2.56 -> 45.52 ms`、low-T `2.58 -> 49.61 ms`、high-T
  `3.02 -> 53.14 ms`、stiff `2.65 -> 47.31 ms`、high-kappa
  `2.75 -> 45.93 ms`；candidate p95 同样显著变慢。
- 所有模式均保持有限；本实验没有触发 physical guard，也没有改变 production
  API 或默认路径。由于 focused accuracy gate 已失败，未进入 full-grid DN、
  parameter-space promotion 或生产 kernel 接入阶段。

### Decision

`REJECTED FOR PRODUCTION`。本结果同时否定了“直接在 Python prototype 中用矩阵指数
替换当前 map”这一实现路径；若未来继续 phase 线，必须先给出 Numba/解析闭式的
低开销 transfer map，并证明其 full-grid DN 与 spectrum 收益，不能重复当前实验。

## Closed-form Numba follow-up for Magnus map (2026-09-16)

为排除上一轮 `scipy.linalg.expm` 的实现性成本，保留同一二阶 Magnus 公式，改用
无迹 2x2 矩阵指数的解析闭式并加 Numba 编译。current-SHA artifact
`docs/transfer_map_phase_numba_round.json` 绑定 `19869ef`，仍固定
Numba=2/workqueue/BLAS=1/reference workers=1，覆盖五个命名 regime、35 个 tail
modes 与每点 25 repeats。

- 35/35 个 tail modes 的 candidate 振幅方向略优于 midpoint，但最大振幅变化只有
  `1.486e-5`，最大相位变化 `7.109e-4 rad`，没有达到 10x handoff 精度改善，
  也没有 full-grid DN/spectrum 证据支持晋升。
- propagation median baseline -> candidate：default `3.254 -> 5.076 ms`、
  low-T `2.667 -> 4.950 ms`、high-T `2.602 -> 4.821 ms`、stiff
  `2.570 -> 4.708 ms`、high-kappa `2.600 -> 4.751 ms`；candidate 在所有
  regime 均超过 5% 成本预算。
- 所有模式保持有限，未改变 physical guard、failure 语义或 production API。

决策：`REJECTED FOR PRODUCTION`。这次结果拒绝的是线性-z 二阶 Magnus 假设本身，
不是单纯的 Python 实现；后续 phase 线不得重复该 map，应转向非振荡 Riccati 或
真正能减少 tensor solve 次数的解析分支。

## Fresh four-regime breakdown (2026-09-16, HEAD 65578a5)

重新以正式 `kink_split=true`、25 repeats、Numba=2、workqueue、BLAS=1、workers=1
测量 default/high-T/stiff/high-kappa，产物为
`docs/profile_numba_round_20260916_{default,highT,stiff,high_kappa}.json`。

| regime | total median | tensor median | mean steps/channel |
|---|---:|---:|---:|
| default | 8.384 ms | 3.174 ms | 8187.4 |
| high-T | 12.379 ms | 5.329 ms | 8143.9 |
| stiff | 12.577 ms | 5.575 ms | 8229.6 |
| high-kappa | 12.449 ms | 5.630 ms | 8043.2 |

tensor solve/propagation remains the largest directly measured component；prepare、
background and PCHIP 均明显较小。结论是下一候选必须减少真实 channel propagation
work（例如经过独立 oracle 验证的 nonoscillatory Riccati/解析分支），不能再做
prepare-layer 或 quadrature micro-optimization。

## Assembly lower-bound profile (2026-09-16, HEAD 7230c6d)

以 `PROFILE_ASSEMBLE=0` 关闭中间列 assembly，仅测 `solve_kernel` 的理论下界；该
模式不产生可用于科学结论的完整输出。25-repeat、2 threads、workqueue、BLAS=1
结果（total/tensor median）为 default `8.090/2.497 ms`、high-T `12.066/4.361 ms`、
stiff `12.291/5.117 ms`、high-kappa `11.988/4.794 ms`，对应正式路径的
`8.384/3.174`、`12.379/5.329`、`12.577/5.575`、`12.449/5.630 ms`。

因此去除 assembly 的 total 收益仅约 `2–4%`，不可能达到速度突破门槛；tensor
kernel 的主要成本仍是 propagation arithmetic/步数。产物为
`docs/profile_noassemble_round_20260916_{default,highT,stiff,high_kappa}.json`，
不进入 production。

## Raw Riccati ratio pole audit (2026-09-17, HEAD d587031)

### Hypothesis and scope

将原始 tensor 方程用 `r=x/y` 化为 Riccati 方程，可能减少双分量传播工作。该
实验严格为 standalone/reference-only；极点是显式终止状态，不允许重启、正则化或
silent fallback。审计脚本为
`scripts/benchmark_riccati_pole_audit.py`，结果为
`docs/riccati_pole_audit_round_20260917.json`。

### Independent evidence

由 `reference._tensor_orig` 直接推导得到
`r' = -2 r - exp(z) (1+r^2)`。在 default、high-T、stiff、high-kappa、两个
`r` edge、high-`T_re` edge 和固定 Sobol 点上，各取 8 个频率（共 64 modes），
raw Riccati 全部在 `z=5` 前触发极点：`64/64 pole`、`0/64 tail reached`、
`0 numerical failure`。对每个终止位置独立用 Cartesian DOP853 重放，均成功，且
`abs(y)/abs(x)` 的最大值为 `1.053e-6`、中位数 `1.000e-6`，确认极点对应物理
的 `y` 零交叉而非 Riccati 实现异常。

### Decision

`REJECTED FOR PRODUCTION`：raw `x/y` Riccati 不能跨越物理零点，不能作为单一
fast solver，也没有资格进入 runtime A/B 或生产 fallback。若继续 phase-function
线，只能研究无极点的复对数导数/非振荡 carrier，并且必须先建立独立 Cartesian/
Prüfer Oracle 对照与显式 guard。

### Fresh speed baseline

当前 HEAD 以 2 threads、workqueue、BLAS=1、reference workers=1、正式
`kink_split=true`、25 warm repeats 重跑四点 profile：

| regime | total median | total p95 | tensor median |
|---|---:|---:|---:|
| default | 6.442 ms | 7.541 ms | 2.633 ms |
| high-T | 9.427 ms | 10.837 ms | 4.513 ms |
| stiff | 9.750 ms | 10.675 ms | 4.931 ms |
| high-kappa | 9.719 ms | 10.450 ms | 4.578 ms |

四份 profile 均记录 current HEAD、25 repeats、线程与 BLAS 配置；tensor
propagation 仍是首要速度瓶颈。

## Consolidated standalone rejection record (2026-09-17, HEAD 85cbfc5)

本次整理将工作区中尚未提交的实验脚本和 artifacts 统一纳入审计范围。所有
standalone 候选都未修改 production path，代码 SHA、Numba/workqueue=2、BLAS=1、
reference workers=1 和重复次数写入对应 JSON。

### WKB carrier / block transfer

预积分 WKB carrier 在 `z_match=3.625` 的四 regime DN proxy 为
`1.497e-4–1.843e-4`，但 candidate 只替代 propagation 后段；结合当前 fresh
total profile，不能达到稳定 10% total runtime 改善，且尚无 Oracle A/B/C full-grid
认证。`z_match=3.75` 的 proxy 为 `8.45e-5–1.07e-4`，仍不足以抵消未认证和总耗时
约束。block transfer 在 block 4/6/8 的振幅误差与速度折衷未满足 production gate。
结论：`REJECTED FOR PRODUCTION`，保留为 reference-only 研究记录。

### Constant-radiation closed form

由 `sigma=4/3` 推导的 standalone map 满足辐射近似下的闭式方程，但真实背景并非
全段常系数。`z_match=3.75`、四命名 regime、25 warm repeats 的最大振幅误差为
default `1.691e-2`、high-T `6.188e-3`、stiff `1.856e-2`、high-kappa `1.811e-2`；
candidate propagation 比 baseline 慢 `10.9–17.5%`。结论：数学近似和 runtime 均不合格，
`REJECTED FOR PRODUCTION`，无 fallback。

### Estimator boundary

当前 estimator coverage artifact 在 9 个已有点上对同一 native interpolant 内部
诊断可显示 coverage=1，但真实 nested native-node 对照的历史记录在多 regime
出现 `false_safe=true`，且部分 `E_nested` 比独立 reference actual error 小数个
数量级。因此不能晋升为科研 error estimator；必须继续以独立 native solves 和
Oracle A/B/C 报告不确定度。

### Current precision/speed boundary

fresh Oracle C true-error 复核显示 default/high-T/stiff/high-kappa 的 deep-anchor
相对差约为 `5.21e-6/2.24e-5/3.35e-6/2.45e-5`，low-T 为 `1.72e-4`；当前真正的
速度热点仍是 tensor propagation，现有候选均未同时满足 full-grid 独立精度与
稳定总耗时突破条件。生产 solver 与 physical guard 语义保持不变。

## Current HEAD fresh outer-reuse headroom (2026-09-17, `f411bbb`)

本轮先从 `fast` 远端 fetch，并重读全部现有研究记录与验证产物；本节只采用当前
HEAD 新生成的 profile/runtime/reuse JSON。资源固定为 Numba=2、workqueue、BLAS=1、
workers=1，正式 `kink_split=true`。

### Fresh baseline

六点 25-repeat warm runtime 为：default `6.99/8.03 ms`、high-T `11.76/12.23`、
stiff `12.38/13.44`、high-kappa `12.20/12.92`（median/p95）；low-T `8.03/8.53`、
low-r `6.98/7.40`。所有点 converged，0 failure。阶段 profile 对 high-T/stiff/
high-kappa 给出 tensor median `7.07/7.62/7.20 ms`，确认传播仍是首要热点。

### Observable-aware reuse headroom

当前 reuse gate 与 forced reuse 的 25-repeat A/B：

| point | runtime forced/current | DN relative | spectrum max abs dex | kernel calls |
|---|---:|---:|---:|---:|
| high-T | `0.7353` | `3.387e-11` | `3.944e-3` | `2 -> 1` |
| stiff | `0.7470` | `3.053e-10` | `1.051e-3` | `2 -> 1` |
| high-kappa | `0.7571` | `3.283e-11` | `1.547e-2` | `2 -> 1` |
| low-T | `0.9760` | `0` | `0` | `1 -> 1` |
| positive-tilt | `0.7052` | `4.162e-9` | `2.075e-1` | `2 -> 1` |

结论：forced reuse 的单值 `DN_gw` 几乎不变，但 spectrum observable 在高-T、stiff、
high-kappa 和 positive-tilt 均可明显偏移；该候选 `REJECTED FOR PRODUCTION`。

### Next hypothesis

仅保留一个未验证的 standalone 方向：在第二次 kernel 之前建立 observable-aware
上界，组合 `delta log Omega`、`delta Phi`、`delta S2`、horizon-crossing shift 与
frequency-weighted DN sensitivity；只有在 25--50 repeat、Cartesian/Prüfer/WKB
独立 oracle、named+Sobol、coverage/p95/p99/false-safe=0 后，才允许进入 runtime A/B。
本轮没有 production code change。

## Observable-aware outer reuse proxy prototype (2026-09-18, `083fdf5`)

为避免继续用固定 `sigma/f_hor` threshold 猜测，本轮新增 standalone
`scripts/benchmark_outer_observable_proxy.py`。它关闭 shortcut、捕获 first→second
真实 outer update，并记录 `delta log Omega`、`delta Phi`、`delta S2`、horizon shift
和按 log-frequency 积分的 DN sensitivity proxy；不改变 production path。

| point | dlogOmega max | dPhi relative | dS2 relative | horizon shift | DN sensitivity proxy |
|---|---:|---:|---:|---:|---:|
| high-T | `3.944e-3` | `7.694e-5` | `6.560e-5` | `2.008e-3` | `3.460e-11` |
| high-kappa | `1.547e-2` | `2.964e-4` | `1.800e-4` | `7.728e-3` | `3.362e-11` |
| stiff | `1.051e-3` | `2.092e-5` | `1.300e-8` | `5.349e-4` | `3.139e-10` |
| default | `3.923e-4` | `0` | `0` | `8.089e-5` | `8.142e-10` |
| positive-tilt | `6.226e-4` | `0` | `0` | `9.590e-5` | `9.034e-10` |

这些量的物理解释边界已明确：`delta log Omega` 是 first/second kernel 输出的事后差，
不能证明在第二次 kernel 前可计算；`delta Phi`/`delta S2`/horizon shift 在 high-T、
high-kappa 上很小，却与明显 spectrum 偏差同时出现；当前 DN proxy 还受末列动态范围
和 weighting 定义影响，不能称为误差上界。原型为 `DIAGNOSTIC ONLY`，不进入
runtime A/B、production 或 fallback；完整原始结果见
`docs/outer_observable_proxy_round_20260918.json`。

### Next hypothesis

研究 kernel 对 `Phi/S2/f_hor` 扰动的局部响应：用已完成的 first solve 输出和背景差分
构造一阶 observable sensitivity，再与第二次独立 Cartesian/Prüfer/WKB solve 对照。
接受前必须在 named+Sobol 上报告 coverage、p95/p99 和 false-safe=0；本轮尚未满足。

## Local observable response sensitivity (2026-09-19, `dca9bd8`)

在 standalone 原型中保存首轮 `solve_kernel` 输入，然后分别把实际 first→second
背景差分注入 `Phi/Phi_mid` 或 `S2/S2inv`，再额外调用 kernel 测量局部响应。该实验
只判断结构，不声称 estimator，也不改变 production。

| point | Phi response prediction | S2 response prediction | actual first→second dlogOmega |
|---|---:|---:|---:|
| high-T | `7.443e-3` | `4.017e-3` | `3.944e-3` |
| stiff | `2.840e-3` | `1.063e-3` | `1.051e-3` |
| high-kappa | `2.795e-2` | `1.544e-2` | `1.547e-2` |

S2-only response 在 stiff/high-kappa 上接近实际频谱偏差，Phi-only response 偏保守；
这比直接使用 `max(delta Phi)`、`max(delta S2)` 更有希望。但 default/positive-tilt
中 Phi/S2 差分为零而 actual dlogOmega 仍为 `3.923e-4`/`6.226e-4`，表明响应还
包括 horizon start、tail matching、`j0`/`z0` 或 kernel phase-path 项。当前 proxy
本身用额外 kernel 求导，不能用于节省 runtime，状态为
`CANDIDATE FOR FURTHER DERIVATION`。

下一实验必须将这些项化为不调用第二次 tensor kernel 的低成本解析响应，并以独立
Cartesian/Prüfer/WKB、named+Sobol、coverage、p95/p99、false-safe=0 验证；未满足前
不得进入 production。

## Analytic outer-response boundary (2026-09-20, standalone)

本轮把 `S2` 响应化为传播端点的闭式缩放项，并加入固定 `j0` 的 horizon-start 缩放项，
然后在全部 named/edge/Sobol 集合上回放 first→second outer update。该回放仍调用额外
kernel 取得真实差分，故只能验证响应结构，不能作为可部署 estimator 或 runtime A/B。

`14/24` 点可比较；令 `pred=max(S2_closed_form,horizon_start_closed_form)`，以
`1.1*pred <= 1e-3` 为暂定安全筛选时 coverage 为 `4/14`、false-safe 为 `0/4`；
eligible actual spectrum-delta 的 p95/max 为 `8.217e-2/2.075e-1 dex`。通过点为
default、negative-tilt、edge_dnre_lo、edge_dnre_hi，目标 high-T/stiff/high-kappa
为 `0/3`。

后三者实际 max spectrum delta 为 `3.944e-3/1.051e-3/1.547e-2 dex`，解析 max
proxy 为 `4.017e-3/1.068e-3/1.544e-2`，均超出 `1e-3` 预算；default 与
positive-tilt 仍存在 S2/Phi 为零但频谱有残差的情况，edge_dnre_hi 也显示两项不完整。
状态为 `REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`，不能宣称已获得 >5% runtime
改善。完整数据见 `docs/outer_observable_proxy_full_round_20260920.json`。

## Adiabaticity-triggered carrier boundary (2026-09-17, standalone)

旧的固定 `z_match` carrier 证据已否决；本轮测试新假设：由
`epsilon=|omega'/omega^2|=|1.5*sigma-1| exp(-z)` 触发，并要求连续 3 个节点满足阈值。
代码在 `scripts/benchmark_wkb_carrier_numba_spike.py`，产物为四个
`docs/wkb_carrier_adiabatic_trigger_*_round_20260917.json`，未改变 production path。

| eps trigger | amplitude p50/p95/max | DN proxy p50/p95/max | target runtime median ratio |
|---|---:|---:|---:|
| `3e-4` | `2.429e-4/5.993e-4/5.993e-4` | `3.292e-11/1.279e-10/1.279e-10` | `0.857..1.033` |
| `1e-3` | `1.304e-3/2.376e-3/2.376e-3` | `1.192e-9/2.973e-9/2.973e-9` | `0.970..1.252` |
| `3e-3` | `6.392e-2/6.654e-2/6.654e-2` | `2.566e-7/4.351e-7/4.351e-7` | `0.912..1.159` |
| `5e-3` | `7.370e-2/7.375e-2/7.375e-2` | `2.368e-6/3.061e-6/3.061e-6` | `0.859..1.019` |

`3e-4` 虽然精度最好，却没有稳定速度收益且 high-kappa 变慢约 28%；更松阈值在
生产 `z_tail=5` 前已产生明显振幅误差。结论为
`REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`；不得将本轮 proxy 当认证。

## CI root cause (2026-09-17)

本轮本地复现新增诊断脚本的 CI 失败：唯一门禁错误为 Ruff `I001`，原因是新增
`CASES` 导入未按规则排序。已修复并复核中文注释门禁、Ruff 与 mypy；以后新增脚本
提交前必须执行同一 import-order 门禁，避免重复触发该原因。

## Outer goal-grid reuse A/B (2026-09-17, standalone)

Fresh 分层 profile（当前 HEAD、16 threads、10 repeats）显示 high-T/stiff/high-kappa
total median `6.557/6.556/6.355 ms`，tensor median `1.701/1.693/1.533 ms`；
fresh 2-thread runtime matrix 的四点 median/p95 为 default `7.907/8.605`、high-T
`12.982/13.932`、stiff `14.263/15.351`、high-kappa `13.584/16.852 ms`。

Invariant probe 发现 outer 第二轮 goal grid 的最大漂移只有约 `1.6e-12`，但不是逐位
相同。随后在 13 个 named/Sobol/edge 点进行了 25-repeat A/B：复用首轮 goal grid 后
`10/13` 个点的最终 `f` SHA256 digest 改变，最大绝对频率差 `1.56e-12`；
`log10OmegaGW`、`DN_gw`、`g2`、`w2` 的数值 max delta 为 0。candidate/current
median ratio 为 `0.960..1.003`，不满足稳定 >5% runtime gate。

这排除了“因漂移很小即可跳过第二次 goal grid construction”的 production 方案：
本项目的准备层去重要求输出 digest 完全一致，近似数值相等不足以接受。候选状态为
`REJECTED FOR PRODUCTION`，原始数据为 `docs/outer_goal_grid_reuse_round_20260917.json`
与 `docs/outer_static_invariants_round_20260917.json`。

## Local derived-param cache (2026-09-17, standalone)

当前 HEAD 的 fresh profile 显示 `gen_expansion` 内多次访问 `derived_param`。候选只在
单次 expansion 调用内把 property 结果绑定到局部变量，不跨 outer iteration 缓存，因而
不改变 DN_eff 更新边界。25-repeat A/B 覆盖 low-T/high-T/stiff/high-kappa，输出
`Nv/sigma/f_hor/log10OmegaGW/DN_gw` digest 全部一致。

| regime | candidate/baseline warm median | digest | decision |
|---|---:|---|---|
| low-T | `0.9963` | equal | reject: <5% |
| high-T | `0.9813` | equal | reject: <5% |
| stiff | `0.9912` | equal | reject: <5% |
| high-kappa | `0.9979` | equal | reject: <5% |

该优化仅带来噪声范围内至约 1.9% 的收益，不能满足本轮速度线；正式代码已撤回，
原始证据保留于 `docs/derived_param_local_cache_round_20260917.json`。原型入口已
经过 Ruff 导入门禁复核，避免重现此前 `I001` CI 原因。

## H2 endpoint cache (2026-09-17, accepted)

fresh 16-thread profile 的 current target total median 为 high-T/stiff/high-kappa
`6.36/6.94/6.95 ms`；goal preparation 为 `0.40/0.40/0.47 ms`。`grid_independent_freqs`
原先对 `N_inf` 的 H2 重复求值，候选在单次函数调用内用浮点 N 作键缓存，未跨 outer
iteration，也未改变 `DN_eff` 或任何 background array 的更新顺序。

50-repeat、13 个 named/Sobol/edge 点 A/B 的最终 `f/log10OmegaGW/DN_gw/g2/w2`
digest 全部一致（mismatch `0/13`），所以 spectrum p50/p95/max 与 DN relative
delta 均为 `0/0/0`，parameter-space false-safe `0`。目标点 runtime ratio 为：

| point | candidate/baseline median | result |
|---|---:|---|
| high-T | `0.9838` | pass |
| stiff | `1.0110` | neutral |
| high-kappa | `0.9346` | pass, 6.54% faster |

独立 Cartesian/Prüfer 交叉验证的 DN relative error 为 high-T/stiff/high-kappa
`9.13e-10/1.996e-9/8.876e-10`，spectrum max relative error 为
`1.490e-7/1.759e-7/2.728e-7`；reheating/kink 路径无新增 failure。该候选满足
输出逐位一致、runtime >5%、精度不退，接受进入 production。

原始证据：`docs/h2_endpoint_cache_round_20260917.json` 与
`docs/h2_endpoint_oracle_highT_20260917.json`、`docs/h2_endpoint_oracle_stiff_20260917.json`、
`docs/h2_endpoint_oracle_high_kappa_20260917.json`。

## Analytic branch eligibility boundary (2026-09-17, standalone)

当前 fresh profile 的 tensor kernel median 为 high-T/stiff/high-kappa
`1.62/1.94/1.84 ms`。新诊断沿每个实际 mode 的 `j0 -> z_tail` 区间统计两端
sigma 是否落在常系数 radiation (`4/3`) 或 stiff (`2`) 邻域内；结果为：

| regime | segments before tail | sigma-near eligible fraction | exact radiation segments |
|---|---:|---:|---:|
| high-T | `618864` | `16.47%` | `0` |
| stiff | `625376` | `14.31%` | `0` |
| high-kappa | `619253` | `17.66%` | `0` |
| low-T | `646759` | `6.62%` | `0` |

该结果只证明 branch 的可覆盖区域边界，不是 solver accuracy estimator；因为目标 regime
均低于 >30% tensor-work 优先门槛，且严格 radiation 段没有出现，暂不建立新的 closed-form
handoff prototype，也不做 production 修改。artifact 为
`docs/analytic_branch_eligibility_round_20260917.json`。

## Unused Psi preparation buffer (2026-09-17)

Fresh round-4 16-thread profile 仍显示 `prep_kernel` 是稳定热点，但源码审计确认其
`Psi` 数组只在 preparation 中写入并返回，正式 `solve_kernel` 和 outer/quad 路径均不读它。
standalone no-Psi prototype 保留完全相同的 spline/primitive、`S2/S2inv`、horizon-start
和 tail arrays；五个 regime 各 50 次逐位 digest 比较均相等。准备层 candidate/baseline
median ratio 为 default `0.9767`、high-T `1.0230`、stiff `1.0239`、high-kappa
`0.9556`、low-T `1.0237`。因此收益不稳定，目标 regime 最佳也只有 `4.44%`。

Decision: `REJECTED FOR PRODUCTION / DIAGNOSTIC RETAINED`。不能用单个 high-kappa
结果越过 `>5%` 门槛，也没有必要为此做 oracle promotion；原始 artifacts 为
`docs/prep_kernel_no_psi*_round_20260917.json`，正式实现未修改。

## Constant tridiagonal buffers (2026-09-17)

`prep_kernel` 的 natural-cubic spline linear solve 中 `aa` 与 `cc` 每个元素恒等于 1。
候选删除这两个 allocation，并将前消元/回代中的乘法改为常数专用形式。focused tests
在正确的固定线程独立进程中通过；目标三 regime 25-repeat 的逐字段 output digest
`f/log10OmegaGW/DN_gw/g2/w2` 全部一致。candidate/baseline total median ratio 为
high-T `0.9605`、stiff `0.9733`、high-kappa `1.0385`，故最高稳定收益不足 5%，且
high-kappa 退化，candidate 已回退。此项不需要独立 oracle promotion，因为未进入
production 且没有数值差异；原始数据见 `docs/profile_tridiag_candidate_*.json`。

## Test-process resource isolation (2026-09-17)

一次本地 focused test 失败的根因是测试启动时未固定 `NUMBA_NUM_THREADS`，之后测试代码
改变环境，Numba 在已初始化线程池后拒绝重新加载线程数。不是 solver candidate 的数值失败。
用启动前 `NUMBA_NUM_THREADS=2`, `NUMBA_THREADING_LAYER=workqueue`, `FAST_THREADS=2`
重跑后 `54 passed, 3 deselected`；正式 benchmark 单独使用 16/16 进程。CI workflow
不包含中文注释/编码语言门禁，远端全矩阵保持 green，避免该类环境错误与旧门禁错误混淆。

## CI gate scope and Numba environment isolation (2026-09-17)

复核确认中文注释/编码语言门禁及脚本已由 `9b979c3` 删除，当前两个 workflow 没有同类
检查；保留的 job 覆盖兼容性、数值正确性、发布归档、集成或静态契约。

资源契约测试的重复失败根因是 pytest 已导入 Numba 后才修改线程环境。测试已改为启动
干净 Python 子进程，在导入数值库前调用 `apply_environment()` 并核对八项环境值，避免
污染后续测试的 Numba runtime。完整 gate 已通过（Ruff/mypy/compileall/diff/manifest；
`160 passed, 6 deselected`）。

## Fresh round-5 profile and rejected phase/primitive candidates (2026-09-17)

当前 HEAD `759d5cc` 的 fresh 16-thread kink profile（25 repeats）为：

| regime | total median/p95 ms | tensor median ms | exact phase primitive median ms |
|---|---:|---:|---:|
| default | 5.35 / 5.74 | 1.08 | 0.90 |
| high-T | 6.86 / 7.89 | 1.81 | 1.35 |
| stiff | 7.49 / 8.26 | 1.92 | 1.49 |
| high-kappa | 7.20 / 8.18 | 1.85 | 1.46 |
| low-T | 4.32 / 5.67 | 1.15 | 0.89 |

这轮重新检验了两个与既有 rejected 实验不同的候选。其一是 exact primitive 的中间
数组融合；五个 regime 各 50 次，输出 digest 完全一致且 max abs/rel 为 0，但 runtime
ratio 为 `1.0355/0.9874/0.9947/1.0024/1.0092`，没有稳定收益，REJECTED。
其二是 phase envelope 中 `exp(z_mid)` 的单次复用；focused tests 通过，kernel 输出
digest 与 baseline 一致，但目标区间 tensor ratio 约 `1.00/1.01/1.01/0.99`，没有
>5% 改善，REJECTED 并回退。两者因未通过 runtime gate，未进入 spectrum/oracle/
false-safe certification；不能把 machine-level digest 相等误报为独立精度认证。

`PROFILE_ASSEMBLE=0` 仅作为 attribution 诊断，估计 assembly 占 high-T/stiff/
high-kappa/low-T tensor 时间约 `0.19/0.17/0.03/0.16 ms`，不足以支持近似删除。
正式 production source 已恢复无 diff；后续候选仍需从 phase kernel 的结构性工作量
或 background/primitive 共享入手，并先通过同口径 fresh A/B。

## Fresh round-6 and kernel specialization rejection (2026-09-17)

HEAD `42b5526` 的 fresh 16-thread kink profile（25 repeats）为：

| regime | total median ms | p95 ms | tensor median ms |
|---|---:|---:|---:|
| default | 4.60 | 5.44 | 1.03 |
| high-T | 6.36 | 7.66 | 1.78 |
| stiff | 7.27 | 8.36 | 2.04 |
| high-kappa | 6.44 | 7.80 | 1.75 |
| low-T | 4.45 | 5.19 | 1.11 |

目标 regime 的第二次 `gen_fast` 调用中只有最后 present-day anchor 改变，其余 grid 节点
逐位相同；static grid reuse 50-repeat ratio 为 high-T/stiff/high-kappa `0.96/0.99/0.98`，
虽然 digest/DN 完全一致，仍无稳定 >5% 收益，REJECTED。

## Phase-substep exponential recurrence prototype (2026-09-17, HEAD 58af509)

### Hypothesis and scope

当前 `solve_kernel` 在每个 phase substep 独立计算 `exp(z)`。standalone
prototype 保持 production 的 `scaled_step`、phase subdivision、kink split、完整
assembly、tail matching 和 `phase_max=0.25` 不变，只将同一 substep 内的
`exp(z)` 替换为一次初值加 `w *= exp(dz)` recurrence；不修改 production API、默认
solver 或 fallback。

### Evidence

- 2-thread、30-repeat kernel A/B：candidate/production median ratio 为
  default/high-T/stiff/high-kappa/low-T = `0.873/0.809/0.899/0.898/0.904`。
- 固定 16-thread、30-repeat formal kernel A/B：ratio 为
  `0.802/0.851/0.926/0.932/0.916`；对应 p95 除 high-T 外均下降或接近，
  说明收益不是只来自单线程调度。
- 完整 outer self-consistency 五工况全部 converged。named + edge + Sobol 共
  14 点全部 status match，`status_mismatch=0`，测试稳定门下
  `false_safe_count=0`；最大 spectrum delta `6.99e-6 dex`，最大 DN delta
  `6.44e-6`，异常来自 `positive_tilt/cr0_blue`，不是目标 high-T/stiff/
  high-kappa 点。
- 目标工况完整 solve 的 spectrum delta dex p95/max：default
  `9.55e-7/2.91e-6`、high-T `1.38e-6/2.66e-6`、stiff
  `1.86e-6/3.30e-6`、high-kappa `8.68e-7/2.42e-6`；low-T 独立为
  `8.68e-7/2.12e-6`。
- fresh Oracle C 当前 HEAD、每点 12 native frequencies：candidate-vs-WKB
  DN relative 为 default/high-T/stiff/high-kappa/low-T =
  `1.073e-5/9.881e-6/1.148e-5/1.087e-5/1.659e-4`；与 production
  的差异仅由 recurrence 造成的 `7.6e-12/3.0e-13/1.2e-11/7.6e-14/
  3.1e-7`，没有新增 oracle systematic。

### Decision

`ACCEPTED AS STANDALONE PROTOTYPE / NOT PRODUCTION`。该方法在目标 regime
满足 >5% kernel runtime gate，且 precision 变化远低于当前 error budget；但输出
digest 不再逐位一致，positive-tilt/cr0-blue 的参数矩阵差异也不能用当前五点
结果外推为全空间认证。若要进入 production，下一步必须把 recurrence 直接接入
opt-in candidate，完成完整 16/20-thread total runtime、Oracle A/Prüfer/WKB
full-grid、reheating/kink edge 与更大 Sobol coverage；在此之前不改正式 fast。

原始 artifacts：
`docs/phase_recurrence_round_20260917.json`、
`docs/phase_recurrence_round_16t_20260917.json`、
`docs/phase_recurrence_fullsolve_20260917.json`、
`docs/phase_recurrence_parameter_matrix_20260917.json`、
`docs/phase_recurrence_oracle_20260917/`。

## Phase recurrence full native-grid oracle audit (2026-09-17, HEAD a1b0393)

新增 standalone `scripts/benchmark_phase_recurrence_fullgrid.py`，在 2-thread、PCHIP、
完整 native goal grid（每点 76 或 77 个频率）上比较 production、recurrence 和独立
Prüfer/DOP853；production path 未修改。

- 五个 regime 全部 `status=ok`，无 status mismatch。
- recurrence 相对 production 的 spectrum p95/max：default `2.23e-6/6.70e-6`、
  low-T `2.08e-6/4.88e-6`、high-T `3.18e-6/6.13e-6`、stiff `4.36e-6/7.60e-6`、
  high-kappa `2.00e-6/5.56e-6`；full-grid DN relative 最高为 low-T `3.12e-7`。
- recurrence 对 Prüfer 的 max spectrum systematic 约为
  `7.31e-3/7.32e-3/7.31e-3/7.30e-3/6.78e-3`（顺序同上），与 production 没有放大；
  该 fixed-DN systematic 不应误报为 recurrence 新误差。

Decision: `STANDALONE PROTOTYPE RETAINED / NOT PRODUCTION`。full-grid 证据补齐了
频率覆盖，但不替代既有 14 点矩阵的 `false_safe_count=0`；仍需 edge/kink、扩大
Sobol、16/20-thread total runtime 和 Oracle A/WKB 独立认证后才能考虑 production。

Artifact: `docs/phase_recurrence_fullgrid_20260917.json`。

no-assembly first-probe kernel 只删除 `assemble=0` 首轮的列装配判断，transfer、kink
split、tail matching 保持不变。25-repeat total ratio 为 default/high-T/stiff/high-kappa/
low-T `0.97/0.96/1.00/0.95/1.01`；digest 相等、DN 差为 0，仍 REJECTED。

组合两项后 high-T/stiff/high-kappa/low-T 为 `0.956/1.015/0.961/1.000`，stiff 退化，
整体不接受。候选相对 baseline 的 spectrum p50/p95/max 均为 0、DN 相对差为 0；但因
runtime gate 失败，没有把等价性误报为 Oracle A/Prüfer/WKB 或 parameter-space
false-safe 认证。production source 未改。
