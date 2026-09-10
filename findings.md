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

## Technical Decisions

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
