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

## Phase recurrence edge/Sobol outer certification and total-runtime gate (2026-09-17, HEAD 3c631e5)

本轮将 recurrence standalone audit 扩展到 Prüfer 参数表的全部 24 点，并把
`edge_tre_hi`、`edge_nt_blue` 的 `shared_Neff_guard` 显式作为 physical guard 记录，
不再把它们误报为数值失败。

- full native-grid Prüfer/DOP853：`22/24` accepted、`2/24` physical guard、无
  numerical failure；accepted 点的 fixed-DN recurrence-vs-production 最大 spectrum
  relative 为 `1.609e-5`（positive-tilt），最大 DN relative 为 `6.44e-6`
  （cr0-blue）。
- full outer self-consistency：base/candidate 均 `22/24` accepted，两个 guard 点
  同步拒绝；`status_mismatch=0`、`false_safe_count=0`；最大 spectrum delta
  `6.988e-6 dex`，最大 DN relative `6.441e-6`。low-T outer delta 保持
  p95/max `8.68e-7/2.12e-6 dex`，DN relative `3.12e-7`，未由高频结果外推。
- formal total runtime 30-repeat、PCHIP、完整 outer solve：16-thread 重跑的
  high-T/stiff/high-kappa median ratio 为 `1.000/0.967/0.986`；20-thread 为
  `1.010/0.939/1.005`。跨线程与重复轮次均不能证明目标 regime 稳定 `>5%`
  total-runtime 改善；p95 也没有形成一致优势。
- accepted 目标点的 total A/B 精度仍为 high-T spectrum p95/max
  `1.38e-6/2.66e-6 dex`、stiff `1.86e-6/3.30e-6 dex`、high-kappa
  `8.68e-7/2.42e-6 dex`；DN relative 分别 `3.0e-13/1.18e-11/7.62e-14`。

Decision: `REJECTED FOR PRODUCTION / RETAINED AS STANDALONE PROTOTYPE`。recurrence
的精度与 guard 语义通过扩大参数空间验证，但 total-runtime gate 未通过；不修改
正式 fast，不声称已达到 default `<4 ms` 或目标 regime `4–5 ms`。后续转向新的
profiler-driven tensor/background 结构候选，不能重复本 recurrence 的同一 runtime
路径。

Artifacts：
`docs/phase_recurrence_fullgrid_all24_20260917.json`、
`docs/phase_recurrence_outer_matrix_20260917.json`、
`docs/phase_recurrence_total_runtime_16t_20260917.json`、
`docs/phase_recurrence_total_runtime_20t_20260917.json`。

no-assembly first-probe kernel 只删除 `assemble=0` 首轮的列装配判断，transfer、kink
split、tail matching 保持不变。25-repeat total ratio 为 default/high-T/stiff/high-kappa/
low-T `0.97/0.96/1.00/0.95/1.01`；digest 相等、DN 差为 0，仍 REJECTED。

组合两项后 high-T/stiff/high-kappa/low-T 为 `0.956/1.015/0.961/1.000`，stiff 退化，
整体不接受。候选相对 baseline 的 spectrum p50/p95/max 均为 0、DN 相对差为 0；但因
runtime gate 失败，没有把等价性误报为 Oracle A/Prüfer/WKB 或 parameter-space
false-safe 认证。production source 未改。

## Inline scalar transfer-map boundary (2026-09-17, standalone)

### Hypothesis and scope

fresh round-8 profile 的目标热点仍为 `tensor_solve_kernel`。本候选不改变数值公式，
只把 production `scaled_step` 的标量 2x2 transfer-map 计算内联到 phase loop，
保留 `exp/cos/sin/sqrt`、phase subdivision、reheating/kink split、tail matching
和 assembly；production fast source 未修改。

### Evidence

- 2-thread full-outer 25-repeat ratio 为 default/high-T/stiff/high-kappa/low-T
  `0.9656/0.9351/0.9488/0.9137/0.9669`；正式目标 ratio 为 16-thread
  `0.9856/0.9974/0.9127`、20-thread `1.0690/0.9876/0.9986`，只有 16-thread
  high-kappa 单点超过 5%，跨 formal 线程不稳定。
- 2-thread outer spectrum delta p50/p95/max（dex）：default
  `2.14e-11/9.55e-7/2.91e-6`，high-T `8.28e-12/1.38e-6/2.66e-6`，
  stiff `1.89e-11/1.86e-6/3.30e-6`，high-kappa `5.92e-13/8.68e-7/2.42e-6`，
  low-T `3.33e-10/8.68e-7/2.12e-6`。
- 对应 DN relative 为 default/high-T/stiff/high-kappa/low-T
  `7.51e-12/3.01e-13/1.18e-11/7.70e-14/3.12e-7`。candidate digest 不逐位
  相同，但变化属于舍入级；这不等于独立 oracle 或 parameter-space 认证。
- formal runtime gate 失败，因此没有进入 Oracle A/Prüfer/WKB、named+Sobol、
  coverage/p95/p99/false-safe gate，避免把局部 2-thread 收益误报为 production 证据。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。不修改 production fast，
不更新 error budget，不声称 default `<4 ms` 或目标 regime `4--5 ms` 已达成。

Artifacts（均绑定 HEAD `0671e1fc549bf45428496935c364956b17b64405`）：
`docs/kernel_inline_spike_round8_20260917.json`、
`docs/kernel_inline_outer_spike_round8_20260917.json`、
`docs/kernel_inline_outer_spike_16t_round8_20260917.json`、
`docs/kernel_inline_outer_spike_20t_round8_20260917.json`。

## Phase exponential hoist boundary (2026-09-17, standalone)

### Hypothesis and scope

fresh round9 tensor profile 显示 high-T/stiff/high-kappa 的 propagation 仍是主热点。
本候选只在 `n_sub=1` 分支复用 phase subdivision 已计算的 `exp(z_mid)`，不改变
transfer map、phase subdivision、kink split、tail matching、assembly 或 production
API；正式 fast source 未修改。

### Evidence

- 2-thread kernel 30-repeat ratio 为 default/high-T/stiff/high-kappa/low-T
  `0.9582/0.9238/0.9285/0.9911/0.9236`。
- 2-thread full-outer 25-repeat ratio 为
  `0.9786/0.9597/0.9610/0.9865/0.9375`；目标 high-T/stiff/high-kappa 均未达到
  稳定 >5% total-runtime 改善，故不进入 formal 16/20-thread gate。
- full-outer spectrum relative p95/max：high-T `3.18e-6/6.13e-6`、stiff
  `4.27e-6/7.60e-6`、high-kappa `2.00e-6/5.56e-6`；DN relative
  `3.01e-13/1.18e-11/7.70e-14`。digest 不逐位相同，不能替代独立 oracle。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。由于最小 full-outer
runtime gate 已失败，不继续 Oracle A/Prüfer/WKB、named+Sobol 或 false-safe 认证。

Artifacts（绑定 spike HEAD `06b7df1f6b1cdab406d7d7bc441dc29716a212a4`）：
`docs/phase_exp_hoist_spike_round9_20260917.json`、
`docs/phase_exp_hoist_outer_round9_20260917.json`。

## Round9 assembly attribution (2026-09-17, current HEAD)

在当前 HEAD `df6880a` 重新运行 `PROFILE_ASSEMBLE=0`，目标三点 tensor stage
median 为 high-T/stiff/high-kappa `3.971/4.284/4.312 ms`，而正式 assemble=1
分别为 `4.409/4.895/4.663 ms`。这只能给出“删除大部分中间列 assembly”的上限，
不是可部署算法；该诊断会丢失中间列，不能用于精度、oracle 或输出契约认证。

完整 total 没有形成稳定的对应 >5% 收益，且 low-T 仍需独立处理。因此 assembly-only
优化排除，不重复 no-assembly kernel 或其组合；下一候选必须减少真实 propagation
工作量，并先满足 hybrid/phase standalone 的精度与 runtime gate。

## Phase fastmath boundary (2026-09-17, standalone)

### Hypothesis and scope

在 phase-exp hoist twin 的 `scaled_step`、phase segment 和 parallel kernel 上启用
Numba `fastmath=True`，期望只改善 trig/标量 arithmetic 调度。该假设不改变 phase
subdivision、handoff、tail、assembly、输出 API 或 production fast source；由于
fastmath 放宽 IEEE 浮点语义，digest 不逐位相同，必须把误差与速度分别看待。

### Evidence

- 2-thread kernel 30-repeat ratio 为 default/high-T/stiff/high-kappa/low-T
  `0.8692/0.8468/0.9018/0.8695/0.9434`，但 full-outer ratio 只有
  `0.9527/0.9620/0.9061/0.9418/0.9957`。
- formal 16-thread ratio 为目标三点 high-T/stiff/high-kappa
  `0.9992/0.9840/0.9680`，20-thread 为 `1.0050/0.9649/0.9771`；收益不稳定且
  未达到 >5% total-runtime gate。default 在 20T 略慢 `1.0038`，low-T 在 16T 慢
  `1.0343`，不能把 2-thread kernel 收益外推。
- 2-thread full-outer 目标点 spectrum delta p50/p95/max（dex）分别为 high-T
  `8.28e-12/1.38e-6/2.66e-6`、stiff `1.89e-11/1.86e-6/3.30e-6`、
  high-kappa `5.92e-13/8.68e-7/2.42e-6`；DN relative 分别为
  `3.01e-13/1.18e-11/7.80e-14`。low-T spectrum max `2.12e-6`、DN relative
  `3.12e-7`。这些是 production 对照差异，不是独立 oracle 或 false-safe 证明。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。runtime gate 失败，
因此不继续 Oracle A/Prüfer/WKB、named+Sobol、coverage/p95/p99 或 false-safe 认证；
不修改 production、不更新 error budget。fastmath 只能作为局部编译实验，不能替代
真正减少 propagation 工作量的 hybrid solver。

Artifacts（均绑定 HEAD `164bd84f2984fce72755c3ad10c5f7698b33663f`）：
`docs/phase_exp_fastmath_spike_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_16t_round10_20260917.json`、
`docs/phase_exp_fastmath_outer_20t_round10_20260917.json`。

## Analytic branch eligibility boundary (2026-09-17, fresh HEAD)

在当前 HEAD `ba95ecc013dbd7c5423e6efdd3f614547a26cf15`、2 threads、BLAS=1、
workers=1 下，先重新统计实际 kernel segment 的背景可解析性，未修改 production。
严格 `sigma=4/3` 的 radiation segment 为 0；严格 `sigma=2` 的 stiff segment 在
high-T/stiff/high-kappa/low-T 分别为 `40202/28140/48294/0`。放宽到相邻节点
`|sigma-sigma_target|<1e-4` 后，eligible fraction 为
high-T/stiff/high-kappa/low-T `16.47%/14.31%/17.66%/6.62%`。

已有固定 `z_match=3.75` radiation closed-form spike 的 30-repeat kernel A/B
也在同一 HEAD 重跑：candidate/exact runtime ratio 为 default/high-T/stiff/high-kappa
`1.149/1.156/1.104/1.137`，且 amplitude max relative 为
`1.69%/0.62%/1.86%/1.81%`。因此固定 handoff 既没有速度收益，也不能作为
误差受控的 analytic branch；后续只测试按 sigma 连续区间局部跳跃，不能把该结果
外推为 branch 可行性。

Artifacts（绑定 HEAD `ba95ecc013dbd7c5423e6efdd3f614547a26cf15`）：
`docs/analytic_branch_eligibility_round11_20260917.json`、
`docs/radiation_closed_form_round11_20260917.json`。

## Radiation local-branch jump boundary (2026-09-17, standalone)

### Hypothesis and scope

如果仅在实际 `sigma≈4/3` 的连续 segment run 上合并调用 radiation exact map，
可以减少真实 propagation 调用；非候选段仍逐段使用 production `_phase_segment`。
prototype 只改变 standalone kernel，未触碰 production。

### Evidence

- tolerance `1e-8/1e-6/1e-4` 的 2-thread 30-repeat 结果显示，前两者几乎没有
  branch hit；`1e-4` 时 branch segment 总数仅约 `7723–8914`，每 mode median
  run length 为 0，说明 eligibility 不是长连续区间。
- `1e-4` candidate/baseline runtime ratio 为 default/high-T/stiff/high-kappa/low-T
  `1.022/1.070/1.020/1.045/1.195`，没有目标 regime 的稳定 >5% 改善。
- `1e-4` amplitude relative p95/max 为 high-T
  `1.51e-4/3.58e-3`、stiff `3.35e-4/2.09e-3`、high-kappa
  `1.51e-4/3.02e-3`；low-T `1.51e-4/5.67e-3`。这不是可接受的 production
  error budget，也没有独立 oracle 支持。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。局部 branch 的连续性与
速度门均失败；不继续 threshold tuning、Oracle 或 parameter-space certification，
转向基于 `|omega'/omega^2|` 的 adiabaticity-trigger prototype。

Artifact（绑定 HEAD `a7773d2e01dcdfd9c387429348b09cdea782edd2`）：
`docs/radiation_local_branch_round11_20260917.json`。

## Adiabaticity-trigger handoff boundary (2026-09-17, standalone)

### Hypothesis and scope

用 `epsilon=|omega'/omega^2|=|1.5*sigma-1| exp(-z)` 替代固定 `z_tail=5`，
在 epsilon 低于阈值时提前进入现有 analytic tail。prototype 保持 production
phase propagation、tail formula、assembly 和 API 不变，仅改变 standalone handoff。

### Evidence

- 2-thread 30-repeat kernel probe 的 threshold `0.02/0.01/0.005/0.001` 在目标
  high-T/stiff/high-kappa 上均比 baseline 慢；以 `0.001` 为例 runtime ratio
  分别为 `1.489/1.425/1.514`，low-T 为 `1.483`。
- 同一 `0.001` 阈值的 spectrum relative p50/p95/max：high-T
  `1.91e-11/2.46e-2/4.63e-1`、stiff `4.35e-11/1.98e-3/4.00e-1`、
  high-kappa `1.36e-12/2.22e-2/4.63e-1`、low-T
  `7.67e-10/1.97e-2/3.45e-1`；DN relative 分别
  `2.76e-8/1.24e-7/6.95e-9/8.45e-3`。
- 误差最大的 mode 在 `sigma≈2/3` 时因 `1.5*sigma-1≈0` 使一阶 epsilon
  虚假变小，提前 handoff；此时尚未满足 tail 的整体 phase/background 条件。
  这证明单独使用 `|omega'/omega^2|` 不是充分的 handoff criterion。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。accuracy gate 已失败，
不进入 16/20-thread runtime、Oracle A/Prüfer/WKB 或 false-safe certification；
不修改 production。后续若再研究 hybrid，必须加入至少一个独立的 phase/curvature
或 mode-local amplitude criterion，不能继续只调 epsilon threshold。

Artifact（绑定 HEAD `84c28e645a4d3e0bfd9329fce83355cb992092cb`）：
`docs/adiabaticity_handoff_round11_20260917.json`。

## Fresh round11 stage profile (2026-09-17, current HEAD)

在 `c5566a825af22342d134824a386a7eb7d32958fd`、Numba=2、BLAS=1、workers=1、
kink split 下，high-T/stiff/high-kappa 各独立 10-repeat fresh profile 的 warm
median（ms）为：total `9.696/10.052/10.410`，tensor kernel
`4.514/5.090/4.660`，`fast_phi_s2_split` `1.378/1.420/1.496`，background
construction `0.670/0.691/0.783`，column integration `0.231/0.181/0.344`。

tensor propagation 仍是最大单项，`fast_phi_s2_split` 是次要但可量化的准备层热点；
本轮不根据 p95 作结论，因为每个独立进程的首次编译/初始化 outlier 使 p95 达到
约 `132–141 ms`。fresh profile artifacts 已绑定该 HEAD；后续候选必须重新用这些
阶段比例设计，不能回到已拒绝的 phase micro-optimizations。

Artifacts：`docs/profile_fast_breakdown_round11_highT_20260917.json`、
`docs/profile_fast_breakdown_round11_stiff_20260917.json`、
`docs/profile_fast_breakdown_round11_high_kappa_20260917.json`。

## Second-order adiabaticity trigger boundary (2026-09-17, standalone)

### Hypothesis and scope

在一阶 `epsilon_1=|q|/omega` 之外加入
`epsilon_2=|q^2+q'|/omega^2`，其中 `q=1.5*sigma-1`，要求连续三个 native
node 同时满足两个界限后才启动 WKB carrier。该项用于抑制 `q≈0` 的假绝热判断；
只修改 standalone `scripts/benchmark_wkb_carrier_numba_spike.py`，不改变 production。

### Evidence

- `eps=3e-4`、2-thread focused 5-repeat：目标 high-T/stiff/high-kappa 的
  amplitude p50/p95/max 分别为 high-T `1.91e-11/9.6e-7/3.83e-6`、
  stiff `4.35e-11/1.15e-5/4.58e-5`、high-kappa
  `1.36e-12/7.3e-6/3.04e-5`；DN proxy 分别为
  `2.54e-13/5.71e-12/2.00e-12`。但 runtime ratio 为
  `1.071/1.026/0.998`，没有稳定的目标区收益。
- `eps=1e-3` 的目标 runtime ratio 为 `0.885/1.039/1.002`，amplitude max
  `3.52e-5/1.36e-4/1.21e-4`；`eps=3e-3` 的 ratio 为 `1.068/0.735/1.016`，
  high-T/high-kappa amplitude max `1.79e-3/1.06e-3`，已接近或超过既定 budget。
- 以上 DN 是同一 standalone carrier 的 proxy 对照，不是 Oracle A/Prüfer/WKB
  independent systematic；由于候选未通过稳定 runtime gate，未继续 full-grid
  oracle、coverage、false-safe 或 p95/p99 certification，避免把局部精度改善误报为
  production evidence。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。二阶项能显著收紧 handoff
误差，但没有产生目标 regime 的稳定速度收益；下一步回到真正减少 propagation
work 的 nonoscillatory phase/Riccati 结构，不再继续调该 trigger threshold。

Artifacts（均绑定 HEAD `db4685fa7cf6aa485b285cde3f1aa7fec2fd740c`）：
`docs/wkb_carrier_second_order_3e-4_round11_20260917.json`、
`docs/wkb_carrier_second_order_1e-3_round11_20260917.json`、
`docs/wkb_carrier_second_order_3e-3_round11_20260917.json`。

## Kink post-transition sigma probe cache boundary (2026-09-17, standalone)

### Hypothesis and scope

formal `kink_split` 的 `gen_fast` 不把 `N_re` 插入网格；因此每次
`fast_phi_s2_split` 都会在同一 outer solver 调用中重复计算相同的三 probe
`sigma_vec([N_re]` 左/点/右)。候选只在 standalone wrapper 中按 model、DN 和
probe 坐标缓存第二次结果，production source 未修改。

### Evidence

- 固定 Numba=2、BLAS=1、workers=1，五点各 50 repeats；`f`、`log10OmegaGW`、
  `DN_gw`、`g2`、`w2` digest 全部逐位一致，收敛状态与 failure reason 一致。
- candidate/baseline warm median ratio：default `0.996`、high-T `1.004`、
  stiff `1.010`、high-kappa `1.008`、low-T `0.916`。目标 high-T/stiff/high-kappa
  没有稳定收益，且 stiff/high-kappa 退化；不满足总耗时稳定改善 `>5%`。
- candidate 每次 solver 调用只保留一次 probe 计算，证明缓存命中；low-T 的局部收益
  不能外推为目标工况收益。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。该准备层去重在
bitwise 结果上安全，但收益不足且目标慢工况退化；不继续放大到 full-grid/oracle
认证，下一步继续围绕真实 tensor propagation 的结构性 standalone prototype。

Artifact（绑定当前 HEAD `b8e71bfffd98504f14084458c724493c03633340`）：
`docs/kink_sigma_cache_round12_20260917.json`。

## Observable-aware outer reuse boundary (2026-09-17, fresh HEAD)

### Hypothesis and scope

用 first-to-second outer update 的 `delta log Omega`、`delta Phi`、`delta S2`、
horizon shift 和频率加权 DN sensitivity 预测是否可以安全跳过第二次 tensor
kernel。实验只做 forced-reuse diagnostic，不改变正式 reuse 判据或 production source。

### Evidence

- 当前 HEAD `8d8102d`、Numba=2、BLAS=1、workers=1 下，forced reuse 的总耗时比现行
  判据为 high-T/stiff/high-kappa `0.748/0.766/0.729`；但 spectrum 最大偏差为
  `3.944e-3/1.051e-3/1.547e-2 dex`，不能作为安全复用。
- 对应 DN 相对差仅为 `3.39e-11/3.05e-10/3.28e-11`，说明 DN sensitivity proxy
  单独不足以约束 spectrum；`delta log Omega` 分别为
  `3.944e-3/1.051e-3/1.547e-2`。
- `delta Phi` 为 `7.69e-5/2.09e-5/2.96e-4`，`delta S2` 为
  `6.56e-5/1.30e-8/1.80e-4`，horizon shift 为
  `2.01e-3/5.35e-4/7.73e-3`；这些量可解释误差来源，但尚未形成无需第二次
  kernel 的充分安全判据。
- default/lowT 当前路径已经只调用一次 kernel；强制复用没有额外收益，不能用其
  结果外推目标慢工况。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS OBSERVABLE-AWARE DIAGNOSTIC`。任何能覆盖
 目标慢工况的宽松 guard 都会放大 spectrum 误差；DN-only 或单一 `Phi/S2` 阈值均不够。
不继续调阈值，下一步回到有新数学依据的 nonoscillatory phase/Riccati standalone
propagation prototype。

Artifacts（均绑定当前 HEAD `8d8102dbfcd9345f43d91b99fa6f656f9c24623f`）：
`docs/outer_observable_proxy_round12_20260917.json`、
`docs/outer_reuse_headroom_round12_20260917.json`。

## WKB-corrected early-tail handoff boundary (2026-09-17, standalone)

### Hypothesis and scope

Oracle C 给出了 frozen tail 的一阶边界修正
`transfer^2 *= 1 + sin(2 theta_handoff) / exp(z_handoff)`。候选将 tensor
Cartesian propagation 的 `z_tail` 从 5 提前到 4，并在 standalone kernel twin
中应用该修正；production solver 未修改。

### Focused evidence

- 当前 HEAD `b0cdf54`、default、Numba=2、BLAS=1、workers=1；本轮 focused
  accuracy/runtime 使用 5 warm repeats，candidate/kernel median ratio
  为 `0.708`，约减少 29% propagation kernel 时间。
- 与同一 native grid 的独立 Prüfer z=10 oracle 比较：现行 z=5 baseline 的
  DN 相对差为 `5.21e-6`，candidate z=4 的 DN 相对差为 `7.40e-3`；candidate
  spectrum 相对误差 p50/p95/max 为 `1.292e-2/1.833e-2/2.407e-2`。
- candidate 相对 production z=5 的 spectrum p50/p95/max 为
  `1.291e-2/1.832e-2/2.391e-2`；这不是可接受的舍入级变化，且已在 focused
  accuracy gate 失败，因此不进入 named/Sobol/full-grid 认证。
- 首次 prototype smoke 曾暴露缺失 endpoint assembly；修正后重跑，当前正式
  artifact 使用完整中间/末端 assembly，错误不再来自装配缺口。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。Oracle C 的一阶修正
不足以把 handoff 提前到 z=4；不继续调 z_tail 或同一 correction，转向真正的
nonoscillatory phase/carrier prototype，并保持 production `z_tail=5`。

Artifact（绑定当前 HEAD `b0cdf549cb0d9b2026a3f51fa3bfa67b81f5e695`）：
`docs/wkb_tail_handoff_round13_20260917.json`。

## WKB carrier preintegration handoff boundary (2026-09-17, fresh HEAD)

### Hypothesis and scope

在 standalone tensor propagation 中，保持 production 的 `z_tail=5`，仅把
`z_match=4` 到 `z_tail=5` 的 carrier 相位积分改为预积分 Simpson 累积，并保留
Cartesian exact propagation 到 handoff；这是与已拒绝的 early-tail truncation 不同的
工作量削减候选，production source 未修改。

### Evidence

- 当前 HEAD `cf183c4`、Numba=2、BLAS=1、workers=1，四工况 25-repeat kernel
  candidate/exact median ratio 为 default/high-T/stiff/high-kappa
  `0.910/0.943/1.051/0.928`。stiff 反而退化 `5.1%`，没有目标工况一致的稳定
  `>5%` 收益。
- candidate 相对 exact twin 的最大 amplitude 相对差为
  default/high-T/stiff/high-kappa `1.352e-4/1.282e-4/1.069e-4/2.071e-4`，tail
  integral 相对差为 `4.963e-5/5.141e-5/5.601e-5/5.911e-5`；四工况
  `numerical_failure_count=0`，但这只是 standalone twin 的局部等价性证据，尚未
  通过独立 Oracle A/Prüfer/WKB。
- smoke 与 formal run 均绑定 `cf183c4`；formal artifact 为
`docs/wkb_carrier_round14_20260917.json`。

## Primitive sigma-node no-copy boundary (2026-09-17, fresh HEAD)

### Hypothesis and scope

formal `kink_split` 网格通常不含精确 `N_re` 节点；此时
`_sigma_node_limits` 返回的 left/right node 数组逐位相同，却仍复制两份。standalone
candidate 只在无 exact kink node 时让两个视图共享输入 `sigma_nodes`，遇到 exact
kink 时保留原有 one-sided copy 和 branch convention；production source 未修改。

### Evidence

- TDD contract 覆盖 formal grid 与显式 exact-kink grid，primitive 四组数组及
  kink metadata 均逐位一致；测试结果 `2 passed`。
- 当前 HEAD `b23190a`、Numba=2、BLAS=1、workers=1，五工况 outer 25-repeat
  A/B 的 candidate/baseline median ratio 为
  default/high-T/stiff/high-kappa/low-T `0.997/0.982/0.987/0.997/0.989`。
- 五工况 `f/log10OmegaGW/DN_gw/g2/w2` digest 全部一致，spectrum max absolute
  dex diff 为 `0`，DN relative 为 `0`，两边均 converged 且 failure reason 为
  `null`；但总耗时改善最高仅 high-T `1.81%`，不满足 `>5%`。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。零复制语义安全但
收益过小，未进入 production；不再继续围绕同一 primitive allocation 做微调。

Artifacts：
`docs/fast_phi_nocopy_round15_default_20260917.json`、
`docs/fast_phi_nocopy_round15_highT_20260917.json`、
`docs/fast_phi_nocopy_round15_stiff_20260917.json`、
`docs/fast_phi_nocopy_round15_high_kappa_20260917.json`、
`docs/fast_phi_nocopy_round15_lowT_20260917.json`。

## FD interpolator lookup cache boundary (2026-09-17, fresh HEAD)

### Hypothesis and scope

`exact_background._fd_from_ref()` 在每次 `H2_vec/sigma_vec` 调用中重复解析同一对
FD interpolator；standalone candidate 在模块加载后缓存 callable identity，并只替换
lookup，不改变插值函数或浮点计算顺序，production source 未修改。

### Evidence

- TDD identity contract 通过：cached lookup 返回与 production lookup 相同的两个
  callable 对象。
- 当前 HEAD `2dd7547`、Numba=2、BLAS=1、workers=1；初始 25-repeat 结果曾出现
  default `0.914`，但扩大到 50 repeats 后 candidate/baseline median ratio 为
  default/high-T/stiff/high-kappa/low-T `1.011/1.004/0.999/1.016/1.012`。
- 50 repeats 中五个输出字段的 baseline/candidate digest unique count 均为 `1`；
  最终 digest 相等、spectrum max delta `0`、DN relative `0`、两边均 converged 且
  failure reason 为 `null`。runtime 结果显示 lookup cache 的 25-repeat 局部收益
  是测量噪声，不能稳定复现。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。不把一次较短重复数的
局部改善误认为速度突破，也不继续微调同一 Python import/cache 路径。

Artifacts：`docs/fd_lookup_cache_round16_*.json`、
`docs/fd_lookup_cache_round16_50_*.json`。

## Cross-channel exact transfer reuse boundary (2026-09-17, fresh HEAD)

### Hypothesis and scope

若多个频率通道共享 horizon start `j0` 或 tail endpoint，可计算一次线性 transfer
path 后复用，减少 tensor propagation work。只做 read-only overlap diagnostic，不
改变 kernel；并显式检查 `(j0, tail_index)` pair，而不是把单独的 `j0` 重合误当作
可复用证据。

### Evidence

- 当前 HEAD `780d70e`、Numba=2、BLAS=1、workers=1；round17 fresh profile 的
  tensor kernel median default/high-T/stiff/high-kappa 为 `2.34/4.29/4.49/4.19 ms`。
- 五工况 `j0` reuse fraction 仅 `9.1%–10.5%`，tail-index reuse fraction 约 `6.5%`；
  但每个 `(j0, tail_index)` pair 均唯一，pair reuse fraction 为 `0`，最大 pair
  multiplicity 为 `1`。
- formal `phase_max=.25` 下，即使 `j0` 相同，不同 `z0` 也可能产生不同的 phase
  subdivision；因此不能用单一 shared transfer map 保持逐位结果。

### Decision

`REJECTED AS EXACT REUSE CANDIDATE / DIAGNOSTIC RETAINED`。不建立简单的跨通道
transfer cache；下一步转向 outer allocation/assembly attribution，或提出能处理
mode-local phase subdivision 的独立数学 prototype。

Artifact：`docs/channel_overlap_round17_20260917.json`。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS STANDALONE SPIKE`。该 carrier 预积分在
default/high-T/high-kappa 有局部收益，但 stiff 退化且未达到跨目标 regime 的稳定
速度门槛；不进入独立 oracle 认证，也不修改 production。下一轮回到 fresh profile
后选择新的、具有独立数学依据的 propagation 工作量候选，避免重复调同一 handoff。
## Outer assembly attribution (2026-09-17, HEAD 9173f0d)

### Fresh diagnostic

本轮先以当前 HEAD 运行 `scripts/benchmark_outer_attribution.py`，固定
Numba=2、BLAS=1、workers=1、workqueue，并对 default、high-T、stiff、
high-kappa、low-T 各 warm 后重复 5 次。脚本只包裹现有 production path 的
stages，不修改计算结果或控制流；产物为
`docs/outer_attribution_round18_20260917.json`。

### Evidence

| regime | total median | primitive | solve per call | full solves/run |
|---|---:|---:|---:|---:|
| default | 6.182 ms | 0.847 ms | 2.450 ms | 1 |
| high-T | 9.094 ms | 0.588 ms | 2.201 ms | 2 |
| stiff | 9.921 ms | 0.690 ms | 2.452 ms | 2 |
| high-kappa | 9.214 ms | 0.636 ms | 2.115 ms | 2 |
| low-T | 6.323 ms | 0.787 ms | 3.123 ms | 1 |

所有观测到的 `solve_kernel` 调用均为 `assemble=1`；慢工况确实因 outer
收敛需要两次完整 assembly，但单纯移除 assembly 或强制 first-probe 的方向
已有独立 profile，收益仅约 2--4%，不足以成为新的速度候选。当前证据因此把
主要优化目标继续收敛到 tensor propagation arithmetic/steps，而不是 tail、
PCHIP 或 outer assembly bookkeeping。

### Decision

`DIAGNOSTIC CONFIRMED / NO PRODUCTION CHANGE`：不放宽 outer reuse，不新增
assembly shortcut；保留可复查的阶段归因脚本和 artifact，下一候选必须减少
真实 channel propagation 工作量，并先以 standalone 数学残差及
Cartesian/Prüfer/WKB 对照证明安全。

## Tail-factor cache spike (2026-09-17, HEAD 3921c8c)

### Hypothesis and scope

tail assembly 在每个 coarse slot 重复计算 `ev_minus[kk2] * fp_minus[kk2]`。
standalone candidate 预计算该 grid-only factor，并只用于 `xf` 的 tail
amplitude；`Th` 保留 production 的原始计算顺序。脚本为
`scripts/benchmark_tail_factor_cache.py`，TDD contract 为
`tests/test_tail_factor_cache_spike.py`，production source 未修改。

### Evidence

当前 HEAD、Numba=2、BLAS=1、workers=1、workqueue，25-repeat kernel A/B：

| regime | candidate/base | spectrum p50/p95/max | DN relative | digest |
|---|---:|---:|---:|---|
| default | 0.9633 | `1.13e-10/2.23e-6/6.70e-6` | `7.51e-12` | mismatch |
| high-T | 0.9763 | `1.91e-11/3.18e-6/6.13e-6` | `2.99e-13` | mismatch |
| stiff | 0.8911 | `1.79e-10/4.36e-6/7.60e-6` | `1.18e-11` | mismatch |
| high-kappa | 0.9610 | `2.48e-12/2.00e-6/5.56e-6` | `7.79e-14` | mismatch |
| low-T | 0.8892 | `3.65e-9/2.08e-6/4.88e-6` | `3.12e-7` | mismatch |

所有 A/B kernel 调用均完成且无 numerical failure；但 factor 重排改变了浮点
运算顺序，无法满足 bitwise digest gate。该结果只证明 tail assembly 有局部
算术冗余，不证明可安全改变科学结果；kernel-only speed 也不能直接外推为
full outer runtime。

### Decision

`REJECTED FOR PRODUCTION / RETAINED AS FAILED SPIKE`：不进入 Oracle A/Prüfer/WKB
或 outer gate，不修改 production；后续若再研究 tail，只能寻找保持原始运算
顺序的缓存边界，不能把近似数值一致误报为严格等价。

Artifact：`docs/tail_factor_cache_round19_20260917.json`。

## Phase loop-state reuse (2026-09-17, candidate source HEAD 6c35c1e)

### Hypothesis and scope

在 `solve_kernel` 的每个 mode/interval 中，`zz = z0 + Phi_grid[k] - Phi0`
已在进入 interval 前计算；原实现仍在 while 条件和 `z_node` 再次计算同一表达式。
候选只复用这个 loop state，不改变 `z_mid_step`、kink 分支 endpoint、transfer
map、assembly 或 tail expression。standalone twin 为
`scripts/benchmark_phase_z_reuse.py`，TDD contract 为
`tests/test_phase_z_reuse_spike.py`。

### Evidence

2-thread、BLAS=1、workers=1、workqueue，25-repeat kernel A/B：五工况 digest、
spectrum p50/p95/max、DN 全部逐位相同，kernel median ratio 为
default/high-T/stiff/high-kappa/low-T `0.880/0.849/0.845/0.895/0.861`。

Full outer 50-repeat（2 threads）结果：

| regime | candidate/base | digest | failure |
|---|---:|---|---|
| default | 0.9260 | equal | none |
| high-T | 0.9970 | equal | none |
| stiff | 0.9755 | equal | none |
| high-kappa | 0.9138 | equal | none |
| low-T | 0.9560 | equal | none |

固定正式资源的 25-repeat full outer 结果同样逐位一致且无 failure：16 threads
ratio 为 `0.984/1.016/0.977/0.985/0.940`，20 threads ratio 为
`0.988/0.978/0.996/0.952/1.007`（default/high-T/stiff/high-kappa/low-T）。
正式线程下收益较小，说明 Python/background preparation 已成为主要限制；未将
kernel 局部收益夸大为全路径 10% breakthrough。

### Decision

`ACCEPTED / STRICT-EQUIVALENCE PRODUCTION PATCH`：这是保持原始浮点运算顺序的
loop-state 去重，满足 bitwise、determinism、spectrum/DN、failure/guard gate，
接入 production；不触碰 phase threshold、outer reuse 或结构性数学路径。

Artifacts：
`docs/phase_z_reuse_round20_20260917.json`、
`docs/phase_z_reuse_outer_round20_20260917.json`、
`docs/phase_z_reuse_outer_round20_50_20260917.json`、
`docs/phase_z_reuse_outer_round20_t16_20260917.json`、
`docs/phase_z_reuse_outer_round20_t20_20260917.json`。

### Post-commit verification

Production patch commit `94144c14c0a79869b2043845b83f13c18ab7e4f1` 的 fresh
full pytest 为 `165 passed, 6 deselected, 2 warnings`。新 stage artifacts
`docs/profile_fast_breakdown_round20_postcommit_*.json` 均直接绑定该 SHA；
post-commit warm total median 为 default/high-T/stiff/high-kappa/low-T
`6.04/9.83/10.15/9.98/6.42 ms`，五点均 converged，failure reason 为 null。
## Round 21 — Cross-layer H2 endpoint reuse rejected (2026-09-17)

Fresh postcommit attribution at `3e5bf7ee82d50a5901d5c66aad886ef909819d80` checked whether `_correct_kink_background` should receive the endpoint values already evaluated by the frequency-grid preparation layer. The current code recomputes `H2_vec([n_re, Nv[-1]], ...)` once per invocation (`stiffgwpy_fast/fast_sgwb.py:785-809`).

Under the fixed study resources (Numba=2, BLAS=1, workers=1, 25 repeats), the complete `correct_kink_background` stage is only `0.10–0.18 ms` warm median against `6.04–10.15 ms` total. Even removing the whole stage would have a theoretical total-runtime ceiling below 2%, so it cannot satisfy the >5% acceptance gate. No production patch was made; changing the `gen_kernel` return contract would add risk without sufficient payoff.

Decision: **REJECTED — insufficient speed headroom**. Artifact: `docs/h2_endpoint_cross_layer_round21_20260917.json`.

## Round 22 — Observable-aware outer reuse screen (2026-09-17)

Fresh current-HEAD profiles at `7742f770c1ebd948612351a66ab47d430a4fd56c` kept the same regime split: high-T/stiff/high-kappa require two kernel calls, with warm kernel medians `4.48/4.74/4.68 ms` and total medians `10.42/10.50/10.28 ms` in the fixed 2-thread study environment. The fresh artifacts are `docs/profile_fast_breakdown_round22_{default,highT,stiff,high_kappa,lowT}_20260917.json`.

The existing outer reuse headroom A/B was rerun for five named cases plus `cr0_blue`, `positive_tilt`, `negative_tilt`, and four Sobol points. Forced reuse saves about 21–28% when it removes the second kernel, but the spectrum changes are `3.94e-3 dex` for high-T, `1.55e-2 dex` for high-kappa, `0.208 dex` for positive tilt, `8.36e-2/4.77e-2 dex` for two Sobol points, and `4.43e-3 dex` for cr0-blue. The small-update rows stiff and sobol_015 are `1.05e-3` and `9.78e-4 dex`, respectively, with tiny DN relative differences.

The S2/horizon observable proxies can screen this sample with a provisional `1.5e-3` threshold, but the proxies are explicitly diagnostic rather than error bounds and provide no independent DN safety certificate. This would be a non-bitwise algorithmic change requiring Oracle A/Prüfer/WKB and broader Sobol validation. **No production outer-reuse threshold was changed; candidate rejected for production pending a certified predictor.** Artifacts: `docs/outer_reuse_headroom_round22_20260917.json`, `docs/outer_reuse_headroom_round22_extended_20260917.json`, `docs/outer_observable_proxy_round22_20260917.json`, `docs/outer_observable_proxy_round22_extended_20260917.json`, and `docs/outer_observable_gate_round22_20260917.json`.

## Round 24 — Numba Phi/S2 primitive standalone spike (2026-09-17)

Fresh current HEAD `72231d76c47732b35d484c42ce07a0c9a6fa532e` was used for a standalone Numba rewrite of the existing `fast_phi_s2_split` array arithmetic. The prototype preserved the kink split and formulas, and its primitive outputs were bitwise identical in all five named cases; primitive warm ratios at 2 threads over 25 repeats were `0.826/0.794/0.866/0.858/0.866` for default/high-T/stiff/high-kappa/lowT.

The full solver did not meet the stable multi-regime gate. At 2 threads and 50 repeats, ratios were `0.950/0.968/0.976/1.007/1.003`; at 20 formal threads and 25 repeats they were `1.042/0.875/0.957/0.986/0.967`. Extended named/tilt/Sobol 25-repeat output digests all remained equal, but `sobol_006` showed a `1.6%` slowdown and the default 20-thread point regressed `4.2%`. **REJECTED_FOR_PRODUCTION**: the primitive win does not transfer stably to end-to-end runtime. Prototype and contract test remain standalone only.

Artifact: `docs/phi_s2_numba_spike_round24_summary_20260917.json` plus per-case `docs/phi_s2_numba_spike_round24_*.json`.

## Round 25 — Phi/S2 workspace reuse (2026-09-17)

The previous arithmetic-only Numba spike was not promoted because its end-to-end
benefit was unstable. This round isolated allocation/copy overhead instead: a
Numba fill helper writes the existing formulas into model-local buffers, and the
next primitive call reuses those buffers after the previous primitive has been
consumed by the outer iteration.

The standalone contract first failed on missing implementation, then passed with
array identity reuse. The production contract first failed on distinct output
identities, then passed after the minimal integration. The formal kink-grid
regressions and the existing fast split accuracy test also pass after correcting
the midpoint buffer length to `n-1`.

Fresh standalone A/B evidence at candidate source HEAD `b22042329f8da8f6aabcc986c89fa49231b43cca`:

| resources | default | high-T | stiff | high-kappa | low-T |
|---|---:|---:|---:|---:|---:|
| 2 threads, 50 repeats | 0.9480 | 0.9267 | 0.9247 | 0.9163 | 0.9171 |
| 16 threads, 25 repeats | 0.8768 | 0.9071 | 0.8957 | 0.9155 | 0.8856 |
| 20 threads, 25 repeats | 0.9181 | 0.9186 | 0.8892 | 0.9249 | 0.8943 |

All listed comparisons are digest-equal, converged, failure-free, with zero DN
relative difference. Eight additional tilt/Sobol cases at 2 threads and 25
repeats were also digest-equal, with ratios `0.9182–0.9770`. The complete
per-case artifacts are the `docs/phi_s2_workspace_spike_round25_*.json` files.

Decision: **ACCEPTED FOR PRODUCTION** pending a fresh post-commit full test and
profile gate. This is strict-equivalence allocation reuse, not an arithmetic or
scientific algorithm change; no Oracle promotion is required.

### Post-commit verification

Production commit `bd2cff602339b386efe9cfdcb3c19056d113e6b8` was pushed to
`fast_v0.2`. Fresh fixed-resource profiles at this SHA used 25 repeats and
Numba=2/BLAS=1/workers=1. Warm total medians were default/high-T/stiff/
high-kappa/low-T `5.64/7.68/8.55/7.87/5.52 ms`; every case converged with no
failure, and repeated spectrum/f/DN/g2/w2 digests were stable. The default
target remains `<4 ms`, so this round does not claim that final target has been
reached. Artifacts: `docs/profile_fast_breakdown_round25_postcommit_*.json`.

## Round 26 — Outer snapshot reference reuse (2026-09-17)

Fresh profile attribution showed that high-T, stiff, and high-kappa are slower
primarily because they execute two outer iterations. The candidate targeted only
allocation/copy overhead: because `gen_fast` replaces `m.Nv`, `m.sigma`, and
`m.f_hor` on the next iteration, the previous outer snapshots can theoretically
be retained by reference rather than copied. No reuse threshold, kernel formula,
or numerical expression was changed.

The TDD identity contract passed. A/B output digests for `f`, spectrum, `DN_gw`,
`g2`, and `w2` were equal in every named case; status and failure were equal and
DN relative difference was zero.

| resources | default | high-T | stiff | high-kappa | low-T |
|---|---:|---:|---:|---:|---:|
| 2 threads, 50 repeats | 0.9905 | 0.9697 | 0.9870 | 0.9713 | 0.9543 |
| 16 threads, 25 repeats | 0.9890 | 0.9312 | 0.9735 | 0.9566 | 1.0056 |
| 20 threads, 25 repeats | 0.9543 | 0.9351 | 0.9877 | 0.9826 | 0.9863 |

**REJECTED FOR PRODUCTION**: the candidate is strict-equivalence safe but fails
the stable >5% gate and has a small formal low-T regression. The production source
was restored; the A/B artifacts remain under
`docs/outer_snapshot_reference_round26_*.json` to prevent repeating this test.
## 2026-09-23 fresh baseline setup

- Repository state: local `codex/fast_v0.2` and remote `fast/fast_v0.2` both resolve to `465196c82b1c938af7bfce125933a4dfda22f4a2`; configured Git email is `2966684515@qq.com`.
- The workspace contains many pre-existing untracked 2026-09-17 experiment artifacts. They are preserved and are not implicitly staged.
- Resource environment for the fresh profile: Windows, affinity `[0, 1]`, Numba 2 threads, `workqueue`, BLAS/OMP/NUMEXPR budgets all 1, Python 3.11.9, NumPy 2.4.4, SciPy 1.17.1, Numba 0.67.0.
- First attempted profile artifacts `docs/profile_fast_breakdown_round28_fresh_20260923_{default,lowT,highT,stiff,high_kappa}.json` used `kink_split=False` because the wrapper invocation omitted `--kink-split`. They are non-canonical diagnostics only and must not be used for fast-profile performance claims; the corrected run must pass `--kink-split`.
- Corrected canonical profile artifacts are `docs/profile_fast_breakdown_round28_canonical_20260923_{default,lowT,highT,stiff,high_kappa}.json`, all bound to the baseline SHA and `kink_split=true`. At 2 threads, affinity `[0,1]`, 25 warm repeats, total medians are `6.361/7.239/9.961/10.288/12.397 ms` for default/lowT/highT/stiff/high-kappa; tensor-kernel medians are `3.259/4.160/5.884/6.529/6.349 ms`. Default/lowT use one kernel call; highT/stiff/high-kappa use two. The canonical tensor kernel remains the dominant stage.
- The first invalid profile exposed a configuration hazard: the profiler's default is not the canonical kink-split fast path. Future profile invocations must explicitly pass `--kink-split` and record that flag in the artifact.

## P0 strict repeated-exp elimination — 2026-09-23

- Candidate: `scripts/benchmark_phase_exp_hoist_strict_spike.py`; production source was unchanged. The candidate computes production `z_mid` and `w_mid = exp(z_mid)` once, reuses `w_mid` for `n_sub` and the `n_sub == 1` transfer step, and retains production `FS.scaled_step` for `n_sub > 1`; no `fastmath` is enabled.
- TDD contract: `tests/test_phase_exp_hoist_strict_spike.py` first failed at collection because the candidate module did not exist, then passed after the minimal twin was added.
- Kernel artifacts: `docs/phase_exp_hoist_strict_round28_kernel_20260923.json`, `docs/phase_exp_hoist_strict_round28_t16_20260923.json`, `docs/phase_exp_hoist_strict_round28_t20_20260923.json`. Named default/lowT/highT/stiff/high-kappa kernel outputs had equal `Ogw/Oj/Opgw/handoff_eps` digests and `DN_gw_relative=0` at 2, 16 and 20 threads; no first divergence was found.
- 2-thread, affinity `[0,1]`, 30-repeat kernel median ratios were `1.0000/0.9951/0.9777/0.9777/1.0045`; 25-repeat full-outer ratios were `0.9843/1.0140/0.9914/1.0028/0.9461` for default/lowT/highT/stiff/high-kappa.
- Formal 16-thread 20-repeat full-outer ratios were `1.0061/1.0033/1.0679/1.0036/1.0016`; formal 20-thread ratios were `1.0274/0.9894/1.0530/0.9821/1.0004`. All formal full-outer spectrum/DN digests, convergence and failure reasons matched, but runtime was not stable and several points regressed.
- Decision: `REJECTED_FOR_PRODUCTION` on end-to-end speed stability. Do not repeat this candidate or relax the bitwise gate; move to the independent P1 assembly-state hypothesis.

## P1 counted assembly-state — 2026-09-23

- Candidate: `scripts/benchmark_counted_assembly_spike.py`; TDD contract first failed at collection because the module did not exist, then passed after implementation. The candidate uses `next_output_k`/`output_slot` only for main-loop scheduling and records a trace; a separate production modulo/division reference trace verified the assembly slot sequence.
- Kernel artifacts: `docs/count_assembly_round28_kernel_20260923.json`, `docs/count_assembly_round28_t16_20260923.json`, `docs/count_assembly_round28_t20_20260923.json`. Named five-regime kernels were bitwise equal, `assembly_nodes_equal=true`, no first divergence, and `DN_gw_relative=0` at 2/16/20 threads.
- 2-thread 30-repeat kernel ratios were `0.9226/0.9466/0.9297/0.9242/0.9322`; 25-repeat full-outer ratios were `0.9792/0.9759/0.9400/0.9437/0.9171`.
- Formal 16-thread 20-repeat full-outer ratios were `1.0528/0.9563/0.9669/0.9838/0.9683`; formal 20-thread 20-repeat ratios were `1.0093/0.9554/1.0187/0.9896/0.9270`. The required canonical 20-thread 50-repeat rerun was `0.9977/0.9915/1.0549/1.0058/1.0064`; all output, convergence and failure gates still matched.
- Decision: `REJECTED_FOR_PRODUCTION` because the formal full-outer result is not stable across regimes; do not modify `solve_kernel` for this candidate. The large isolated kernel win does not establish a production win.

## P2 LLVM/ASM audit — 2026-09-23

- Initial inspection hit Numba's cache limitation (`Inspection disabled for cached code`); after root-cause confirmation, the audit calls `Dispatcher.recompile()` before reading LLVM/ASM. The corrected artifacts are `docs/fast_kernel_llvm_audit_round28_t2_20260923.json` and `docs/fast_kernel_llvm_audit_round28_t20_20260923.json`.
- The canonical `solve_kernel` signature is specialized with `h_arr=None` and `Sv` as a contiguous float64 array. LLVM still contains 9 integer division and 9 integer remainder operations; ASM is large (`115158` bytes) with 1474 `rsp` references under the current counting heuristic.
- `_phase_segment` LLVM contains two `exp` call sites in the inlined path (one for substep selection and one inside the one-step transfer); `scaled_step` has one LLVM exp call and separate sin/cos calls. The corrected audit reports no LLVM fast-math flags for these production functions.
- P2 confirms that generic constant/optional branches remain visible, so a guarded canonical-fast specialization is a justified next standalone hypothesis. It must preserve production modulo/division scheduling to remain independent of P1.

## P3 canonical-fast specialization — 2026-09-23

- Candidate: `scripts/benchmark_canonical_fast_specialized_spike.py`; it fixes only canonical fast constants and known optional facts while retaining production modulo/division scheduling. TDD collection first failed because the module was absent, then the default bitwise test passed.
- Artifact `docs/canonical_fast_specialized_round28_t2_20260923.json`: five named kernel digests matched at 2 threads; 30-repeat kernel ratios were `0.9392/0.9386/0.9450/0.9172/0.9518`; 25-repeat full-outer ratios were `0.9839/0.9577/0.9732/0.9481/0.9954`.
- Artifact `docs/canonical_fast_specialized_round28_t20_50_20260923.json`: five named kernel digests matched at 20 threads; canonical 50-repeat full-outer ratios were `0.9644/0.9934/0.9794/0.9741/1.0129`. Spectrum/DN/convergence/failure gates matched in every row, but high-kappa regressed by about 1.3%.
- Decision: `REJECTED_FOR_PRODUCTION`; no guarded dispatch or production source change is justified by a candidate with a repeatable formal-regime regression. The three candidates are kept as diagnostic evidence only.

## Verification note — 2026-09-23

- The new contracts `tests/test_phase_exp_hoist_strict_spike.py`, `tests/test_counted_assembly_spike.py`, and `tests/test_canonical_fast_specialized_spike.py` pass independently.
- The scoped existing run `tests/test_fast_sgwb.py tests/test_resource_budget.py` produced 3 failures / 45 passes. Two failures are stale monkeypatch tests that define `solve_kernel` doubles with 19–23 positional parameters while the unchanged baseline call site passes 26 parameters. The third `test_eval_freqs_are_native_grid_nodes` aborts in the unchanged transition-refine path with the existing `shared_Neff_guard`.
- No production source file was changed in this round, so these failures are recorded as baseline compatibility/guard issues and are not silently repaired within this performance-only task.

## Autonomous continuation bootstrap — 2026-09-23

- Fresh fetch completed before new work. Local `codex/fast_v0.2`, remote-tracking `fast/fast_v0.2`, and `HEAD` are all `1750772a1c70cb3f99edc054fd7841220bdf3601`; no tracked source changes are pending.
- Production path is `stiffgwpy_fast.fast_sgwb.SGWB_iter_fast` -> `gen_fast`/frequency preparation -> `solve_kernel` -> frequency quadrature and guards. Canonical profiling must call `kink_split=True`, use goal grid plus PCHIP, fixed affinity/resources, alternating baseline/candidate measurements, and report median/p95.
- Existing negative knowledge remains binding: repeated-exp hoist, counted assembly state, and canonical specialization all fail the formal full-outer stability gate; outer reuse and previous arithmetic/allocation micro-optimizations also lack stable end-to-end headroom.
- The user-authorized continuation loop is now active. Next selection must begin with Amdahl headroom from a fresh canonical profile, then use a mathematically distinct standalone prototype rather than mechanically reopening Round 28 candidates.
- Round 28 canonical medians remain the usable Amdahl baseline: total/tensor ms are default `6.36/3.26`, lowT `7.24/4.16`, highT `9.96/5.88`, stiff `10.29/6.53`, high-kappa `12.40/6.35`; tensor shares are `51.2%/57.5%/59.1%/63.5%/51.2%`. A candidate that cannot reduce this propagation work is screened out before implementation.
- Fresh Round 29 profile at `1750772a1c70cb3f99edc054fd7841220bdf3601` used five named cases, `kink_split=true`, 25 repeats, affinity `[0,1]`, Numba 2/workqueue and BLAS-family budgets of 1. Total/tensor medians are default `5.984/3.144 ms`, lowT `6.202/3.670 ms`, highT `9.503/5.587 ms`, stiff `9.360/6.133 ms`, high-kappa `9.183/5.694 ms`; tensor shares are `52.5%/59.2%/58.8%/65.5%/62.0%`, with tensor p95 `3.713/4.407/6.218/6.987/6.093 ms` in the same case order. Hard cases still execute two propagation calls; average propagation steps/channel are about `8.19k/8.51k/8.14k/8.23k/8.04k`.

## Literature scan — 2026-09-23

- Bremer's phase-function work shows that a slowly varying nonoscillatory phase can make cost largely independent of oscillation magnitude, but the construction uses a nonlinear Kummer/Riccati equation and needs explicit treatment of turning points; this is mathematically distinct from merely reusing `exp(z)` or changing assembly state. Sources: SIAM adaptive spectral phase method and arXiv turning-point phase method.
- Lorenz/Jahnke/Lubich's adiabatic midpoint and Magnus integrators support larger-than-period step sizes for time-varying high-frequency second-order systems, but their error analysis assumes a controlled adiabatic transformation. A safe prototype must measure the local residual/commutator and keep a hard fallback to the exact production transfer, not silently widen the step.
- SIMD literature supports SoA/AoSoA and padding for ensembles of independent ODE trajectories, while Numba documents `prange` as parallel rather than a guarantee of cross-mode SIMD. This makes a batched multi-frequency kernel plausible, but only after measuring LLVM vectorization and guarding lane divergence; it is higher implementation cost than a one-mode phase prototype.
- Initial hypothesis ordering: (H1) residual-controlled adiabatic/phase transfer with exact fallback, highest potential and information value; (H2) batched SoA/AoSoA propagation, potentially high throughput but likely lane-divergence risk; (H3) certified second-outer predictor/reuse, highest speed headroom but non-strict and already failed proxy-only gates, so requires an independent observable certificate before any implementation.
- Historical audit narrows H1: raw `x/y` Riccati hit poles in all tested modes, fixed-step Prüfer RK4 had 0.97–2.37% amplitude error, WKB/adiabatic handoff candidates either regressed accuracy or runtime, and linear-z Magnus was already rejected as a method family. Therefore H1 is retained only as a future global phase-function study, not the next mechanical spike.
- Next concrete experiment selected: H2 grouped SoA/AoSoA propagation prototype. It is mathematically strict (same Cartesian transfer per lane), changes only execution layout/order, and can be rejected cheaply if `j0`/tail divergence causes excess work or LLVM fails to vectorize. Before coding, measure the current `j0` distribution and bucket widths to set an honest work-overhead bound.
- TDD contract for `tests/test_grouped_soa_spike.py` first failed on the absent module, then passed after adding the standalone implementation. The five-case 2-thread/30-repeat kernel A/B at bucket width 32 was bitwise equal for `Ogw/Oj/Opgw/handoff_eps`, but candidate/base medians were default `1.3143`, lowT `1.3500`, highT `1.3217`, stiff `1.3293`, high-kappa `1.3418`; p95 ratios were `1.2903/1.4086/1.3471/1.3558/1.3848`. Measured lockstep work overhead was only `5.15–5.57%`, so the extra `31–35%` runtime is implementation/layout overhead rather than arithmetic work; this is a strong early negative signal, but a full-outer confirmation is still required by the contract.
- Full-outer 2-thread/25-repeat confirmation preserved spectrum, `DN_gw`, `g2`, `w2`, failure and convergence bitwise in all five cases, but candidate/base median ratios were default `1.2028`, lowT `1.2432`, highT `1.1864`, stiff `1.2323`, high-kappa `1.2080`; p95 ratios were `1.0974/1.0920/1.1171/1.3418/1.2220`. LLVM/ASM inspection found no `<N x double>` vector tokens in the grouped kernel despite loop-vectorization blocks (the visible vector loops are integer/layout machinery); LLVM had 16 exp, 7 sin and 7 cos references, with 258 assembly calls. Decision: `REJECTED_FOR_PRODUCTION`; no formal 16/20-thread gate or production change.
- A read-only first-iteration capture exposed a materially different outer hypothesis: the first full propagation's `DN_gw` is about `0.434x` of the final self-consistent value in all five named cases (`0.00098309 -> 0.00226365` default, `0.0245086 -> 0.0564331` highT, `0.00650761 -> 0.0149843` stiff, `0.0980231 -> 0.225706` high-kappa, `2.30674e-8 -> 5.31146e-8` lowT). Sparse frequency subsets do not fix this: 16-node first-pass estimates remain near the first-pass value, with `42–59%` relative error to the final value. The stable factor suggests an initial-map fixed-point predictor may have headroom, but it must be tested on edge/Sobol and with an actual one-pass full solve; do not assume the factor is universal.
- Fixed-gain one-pass predictor screen (`gain=2.30`, coarse count 32) failed the first named correctness gate before any performance promotion: default candidate `DN_gw=9.9067e-4` versus converged `2.26365e-3` (`56.2%` relative) and spectrum max difference `5.006e-3 dex`, although status and deterministic replay matched and the frequency grid was identical. Decision: reject this fixed-gain hypothesis; no gain tuning on the same sample.

## Round 30 correction — predictor harness audit and final rejection (2026-09-23)

- The first predictor result had a standalone harness defect, not a production
  numerical finding: `_evaluate_once` used the PCHIP integral in `d log10(f)`
  directly, while the production path multiplies the native integral by
  `ln(10)` to obtain `d ln(f)`. A TDD regression test reproduced the mismatch;
  the helper now uses the exact production vectorized PCHIP measure.
- After correction, the fixed-point first iterate was set to the only
  assumption justified by the data, `x1=x0+g(x0)` (`gain=1.0`). The corrected
  five-case screen remained deterministic and status-equivalent, but failed
  the spectrum gate in every named case: max differences were
  `5.006e-3/5.905e-3/4.905e-3/3.342e-3/5.044e-3 dex` for
  default/low-T/high-T/stiff/high-kappa, all above `1e-3 dex`.
- Root cause of the persistent mismatch is production semantics, not a cache:
  the fast outer loop may accept after one propagation and set the final
  `DN_eff` without regenerating the already-returned spectrum/background at
  that final value. The predictor's one full solve at `x1` therefore does not
  reproduce the production output at `x0`; forcing it to solve at `x0` removes
  the proposed work saving. The current output contract makes this predictor
  unsuitable for promotion, so no gain tuning or second-order secant variant
  will be reopened on this path.
- Artifact: `docs/coarse_fixed_point_round30_20260923.json`. The prototype and
  contract test remain standalone; production source is unchanged.
- Decision: **REJECTED_FOR_PRODUCTION**. Next concrete experiment is a
  residual-controlled local phase/adiabatic transfer prototype with exact
  production fallback, not another outer predictor or grouped-layout retune.

## Round 31 — residual-controlled two-step transfer (2026-09-23)

- Fresh current-HEAD canonical profile at `fc26e3408a507cbe8b9b4c56df5199773739895c`
  used `kink_split=true`, goal/PCHIP, fixed workqueue/BLAS budget and 25
  repeats. Total medians were `3.866/3.478/5.816/5.529/5.774 ms` for
  default/low-T/high-T/stiff/high-kappa; tensor propagation was the largest
  stage at `29.2%–33.7%`, and high-T/stiff/high-kappa still used two calls.
- The standalone candidate applies a two-native-step midpoint exponential only
  when a local variation/phase defect proxy passes; otherwise it reproduces the
  production `_phase_segment` path, including the exact kink split and tail
  handoff. TDD first exposed the expected missing-module red state, then a
  resource-initialization issue and two fallback-state mismatches; all were
  fixed before measurement. At threshold `0`, the candidate is bitwise equal
  with zero accepted blocks, confirming the fallback twin.
- At threshold `1e-3`, named kernel medians were `0.619–0.683x` of baseline
  and accepted roughly `75k–97k` blocks. Independent continuous reference
  screens on five named cases plus eight edge/Sobol cases showed the candidate
  changed production `DN_gw` by only `3.14e-5–3.60e-5` relative and did not
  materially worsen the existing reference residual. `edge_tre_hi` was
  explicitly recorded as the pre-existing `shared_Neff_guard` before the
  candidate ran; it is not a regression.
- Full outer, warm and alternating, 25 repeats at two threads had median
  ratios `0.898/0.876/0.867/0.838/0.876` for default/low-T/high-T/stiff/
  high-kappa; output/failure/convergence semantics matched, with max spectrum
  differences about `1.88e-5 dex` and DN relative differences
  `3.14e-5–3.60e-5`. Formal 16-thread ratios were
  `0.972/0.975/0.921/0.990/0.991`; formal 20-thread ratios were
  `1.000/0.969/0.981/0.982/1.021`. The high-kappa formal regression and lack
  of stable multi-resource >5% end-to-end gain reject production promotion,
  despite the low-thread kernel win.
- Artifacts: `docs/residual_block_transfer_round31_kernel_20260923.json`,
  `docs/residual_block_transfer_round31_oracle_20260923.json`,
  `docs/residual_block_transfer_round31_oracle_extended_20260923.json`,
  `docs/residual_block_transfer_round31_outer_20260923.json`,
  `docs/residual_block_transfer_round31_outer_t16_20260923.json`, and
  `docs/residual_block_transfer_round31_outer_t20_20260923.json`.
- Decision: **REJECTED_FOR_PRODUCTION**. The next concrete experiment is a
  fourth-order two-exponential commutator-free Magnus block over four native
  intervals, with an independently computed local defect gate and exact
  production fallback. This is a different integrator order/composition, not
  a threshold or bucket-width retune.

## Round 32 — fourth-order commutator-free Magnus screen (2026-09-23)

- The prototype implements the published Gaussian-node two-exponential CF4
  composition for the existing 2x2 Cartesian transfer matrix, with a local
  commutator defect gate and exact production fallback. Threshold zero passed
  full intermediate/final output bitwise equality, not only the final column.
- At threshold `1e-3`, the candidate's independent-oracle residual remained
  close to baseline on all five named cases: candidate DN-relative-to-oracle
  was `0.001514–0.003639` versus baseline `0.001514–0.003636`, while the
  candidate spectrum max-dex residual was `0.002934–0.003163`. This was a
  scientific screen only; no production source was changed.
- The full-assembly kernel screen did not provide a stable multi-regime win:
  at threshold `1e-3`, 2-thread candidate/base ratios were
  `0.913/1.069/0.812/0.824/0.790` for default/high-T/stiff/high-kappa/low-T;
  high-T regressed and the candidate's full-output path is not a formal
  end-to-end improvement across regimes. Lower threshold `1e-4` also left
  high-T and high-kappa above baseline. Since the core kernel prerequisite
  already fails the stable named-regime gate, no misleading full-outer timing
  expansion was run.
- Artifacts: `docs/cf4_magnus_round32_kernel_20260923.json` and
  `docs/cf4_magnus_round32_oracle_20260923.json`.
- Decision: **REJECTED_FOR_PRODUCTION**. Do not retune CF4 coefficients,
  threshold or phase cap on this sample. Next candidate: a phase-function
  Chebyshev representation feasibility spike with residual certification,
  explicitly avoiding the previously rejected raw Riccati-pole and
  first-order WKB-handoff implementations.

## Round 33 — phase-function feasibility bootstrap (2026-09-23)

- Fresh canonical profile after fetch at `00e2bfd5b27261363015a5f293e3a3ca118800ac`
  used the same goal/PCHIP and fixed resources. Total/tensor medians were
  `3.961/0.923` default, `3.389/1.005` low-T, `5.543/1.571` high-T,
  `5.702/1.876` stiff, and `5.745/1.783` high-kappa ms; tensor shares were
  `23.3%–32.9%` and hard cases still used two calls.
- A read-only Kummer-equation feasibility screen used the real prepared
  background, smooth post-horizon interval `z=2..5`, and WKB-initialized
  positive phase amplitude. The solve stayed finite/positive for the sampled
  modes in all five named regimes (`4/4` usable modes per regime; max rho
  about `0.949`). This is only a mathematical feasibility signal, not a
  correctness or speed result; the next standalone must fit the phase/amplitude
  representation and certify its Kummer residual against Cartesian/Prüfer.
- Next concrete experiment remains the residual-certified phase-function
  Chebyshev representation spike. No production source has been changed.

### Round 33 phase-function result

- The standalone Kummer representation was corrected to use the positive
  amplitude constructed from two independent scalar solutions, avoiding the
  unstable WKB-initialized nonlinear Kummer IVP. The positivity gate then
  passed `3/3` usable windows in every named regime.
- The useful-window compression gate failed: degree-12 Chebyshev fits over
  real `z=2..4` windows had residual maxima about `12.5–57.6`, relative rho
  fit errors about `4.0–4.3`, and phase fit errors about `0.93–1.10`. A very
  short `z=2..2.2` diagnostic window reached residual `4.1e-5`, but contains
  too little propagation work to justify a speed path.
- This rejects the specific fixed-window Chebyshev representation before any
  runtime or production test. It does not reject boundary-conditioned Kummer
  phase functions, which remain a materially different next prototype.
- Artifact: `docs/phase_function_chebyshev_round33_20260923.json`.

## Round 34 hypothesis — boundary-conditioned Kummer phase (2026-09-23)

The Round 33 amplitude was formed from an arbitrary fundamental basis, so its
Kummer solution inherited large oscillatory modulation even though the equation
itself was positive. The next standalone hypothesis is materially different:
construct a complex scalar solution by imposing an outgoing-WKB amplitude and
derivative at the smooth right boundary, integrate that solution backward, and
use its modulus as the Kummer amplitude. A valid result must show a long smooth
interval, positive amplitude, small Kummer residual, and agreement with the
Cartesian/Prüfer transfer reference. A short-window fit or a lower polynomial
degree is not sufficient evidence and must be rejected.

This follows the phase-function literature's boundary-conditioned special
solution idea, but the prototype remains diagnostic and does not alter
production code.

## Round 34 result — boundary-conditioned Kummer phase (2026-09-23)

Fresh canonical profile at `e40a2d0f5373843bfe3fc66a0173741be2346882` used
fixed 2-thread/workqueue/BLAS resources, 25 repeats, goal/PCHIP and
`kink_split=true`. Total medians were `5.528/5.557/10.170/9.533/10.319 ms`
for default/low-T/high-T/stiff/high-kappa; tensor medians were
`2.336/2.754/4.767/5.010/5.003 ms`, and hard cases used two propagations.
Artifacts are `docs/profile_fast_breakdown_round34_20260923_*.json`.

The right-boundary-WKB-conditioned complex solution stayed finite and positive
for the usable modes, but degree-12 compression was not a phase function in
the required numerical sense: Kummer residuals ranged from about `1.0` to
`41.9`, Cartesian transfer relative errors from `0.25` to `1.60`, and the
usable windows were only `Delta z≈2` (some `Delta N≈1`). The specific
boundary-conditioned WKB implementation is **REJECTED_FOR_PRODUCTION** and
was stopped before runtime; no degree/window/tolerance retuning is allowed.
Artifact: `docs/phase_function_boundary_round34_20260923.json`.

The next hypothesis is deliberately non-structural: a small-angle polynomial
transfer twin that preserves the exact production phase subdivision and step
sequence while replacing only `sin`/`cos` on the certified phase-cap branch.
It targets the LLVM-confirmed transcendental cost without reopening rejected
phase-function, recurrence, or grouped-layout experiments.

## Round 35 result — small-angle polynomial transfer (2026-09-23)

Fresh current-HEAD profiles at `27b9fc5db56b069945d3decfd17d8c0470cf79c0`
used fixed 2-thread resources, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `5.637/2.330` default, `5.479/2.887` low-T,
`9.053/4.250` high-T, `8.478/4.633` stiff and `9.721/4.462` high-kappa ms;
hard cases still executed two propagation calls. Artifacts are
`docs/profile_fast_breakdown_round35_20260923_*.json`.

The candidate preserved exact phase subdivision, assembly nodes, tail handoff,
failure and convergence semantics, and replaced only the `sin`/`cos` formulas
on the `|omega*h|<=0.25` branch with degree-10 Taylor forms. Local transfer
absolute error was `5.55e-17`. Named and extended independent-reference
screens passed; candidate spectrum errors to baseline were `2.1e-6–3.3e-6`
dex and the largest extended DN relative difference was `6.44e-6`.

The isolated 2-thread kernel improved, but full end-to-end did not satisfy the
production gate. 2-thread ratios default/high-T/stiff/high-kappa/low-T were
`0.976/0.966/1.000/0.923/0.969`; 16-thread ratios were
`1.001/1.003/0.986/0.896/0.965`; 20-thread ratios were
`1.013/1.045/0.962/0.969/0.946`. The implementation is
**REJECTED_FOR_PRODUCTION**; no formal promotion or trigonometric retuning.
Artifacts: `docs/small_angle_transfer_round35_{kernel,oracle,
oracle_extended,outer,outer_t16,outer_t20}_20260923.json`.

Next selected hypothesis is a separate range-reduced polynomial for the
per-transfer `exp(z)` value, while retaining the production exponential for
the phase-subdivision count. This isolates the other LLVM-confirmed
transcendental cost and will be rejected early if range-reduction overhead
outweighs the saved libm call.

## Round 36 result — range-reduced exponential transfer (2026-09-23)

Fresh profiles at `f94785e0a99d46331216d73a1177f6d91798965c` used the fixed
2-thread/workqueue/BLAS contract, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `6.209/2.750` default, `5.543/3.125` low-T,
`9.157/4.486` high-T, `9.082/4.908` stiff and `9.244/4.649` high-kappa ms.
Artifacts: `docs/profile_fast_breakdown_round36_20260923_*.json`.

The standalone candidate kept the production `exp` for phase-subdivision
counts and used a range-reduced degree-14 Taylor evaluation only for transfer
values. Its local errors were good (`5.13e-16` relative for exp and
`1.11e-16` absolute for the transfer), but the extra integer/range-reduction
and Horner work dominated: default kernel median/p95 were `4.576/4.873 ms`
versus baseline `2.304/2.529 ms` (`1.986x`). It failed the kernel prerequisite
and is **REJECTED_FOR_PRODUCTION** without oracle or full-outer expansion.
Artifact: `docs/range_exp_transfer_round36_default_20260923.json`.

The next hypothesis is a tangent/sensitivity-assisted outer correction. It is
not the rejected fixed-point predictor: it must produce a first-order response
for the spectrum/background and preserve the final output semantics, with an
early cost comparison against the second real propagation.

## Round 37 result — full-state tangent outer screen (2026-09-23)

Fresh profiles at `c2ff97a10b5ae31b2e623980feabd0ad6eaf784b` used the fixed
2-thread/workqueue/BLAS contract, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `5.604/2.462` default, `5.479/2.913` low-T,
`8.237/4.355` high-T, `8.246/4.515` stiff and `8.321/4.341` high-kappa ms.
Artifacts: `docs/profile_fast_breakdown_round37_20260923_*.json`.

The standalone screen derived the exact constant-z transfer Jacobian and
carried a tangent `(x,y)` state alongside the production state. The one-step
Jacobian matched a central finite difference to `2.62e-11`, including the
analytic `omega->0` limit. However, the long propagation state became
non-finite for default modes even without assembly/tail, and the median cost
was `110.8 ms` versus `2.69 ms` for the full production kernel (`41.1x`).
This naive full-state response is **REJECTED_FOR_PRODUCTION** on both numerical
stability and cost; no outer correctness gate was attempted.
Artifact: `docs/tangent_outer_round37_default_20260923.json`.

The next distinct response hypothesis is an observable/adjoint screen: avoid
retaining a full tangent vector per frequency and estimate only scalar
integrated-output sensitivity, with finite-difference verification before any
runtime promotion.

## Round 38 result — observable adjoint outer screen (2026-09-23)

Fresh profiles at `85cec7aa40e19840a762ce5f151dfb7c604ae736` used the fixed
2-thread/workqueue/BLAS contract, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `5.419/2.429` default, `5.247/2.902` low-T,
`8.280/4.321` high-T, `8.194/4.508` stiff and `8.329/4.343` high-kappa ms.
Artifacts: `docs/profile_fast_breakdown_round38_20260923_*.json`.

The standalone adjoint used the exact transpose of the constant-z transfer;
the one-step contract error was `0`. On the long default propagation, however,
the no-assembly/tail adjoint became non-finite and had median/p95
`65.1/71.1 ms` versus full baseline `2.68/3.05 ms` (`24.3x`). More
fundamentally, a scalar adjoint cannot reconstruct the API's full spectrum.
This route is **REJECTED_FOR_PRODUCTION** on API, numerical-stability and
Amdahl grounds; no outer benchmark was run.
Artifact: `docs/adjoint_outer_round38_default_20260923.json`.

The next selected hypothesis is a combined transfer-map lookup/interpolation
table. It is distinct from the rejected Taylor and range-reduced polynomial
forms: the table approximates the already-composed transfer coefficients and
will be rejected early if memory lookup/interpolation costs exceed libm.

## Round 39 result — combined transfer-map table (2026-09-23)

Fresh profiles at `0d0d77c73478cd4fa0f3dfaa1773ff215616647a` used the fixed
2-thread/workqueue/BLAS contract, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `5.337/2.489` default, `5.127/2.757` low-T,
`7.823/4.222` high-T, `8.608/4.723` stiff and `8.086/4.486` high-kappa ms.
Artifacts: `docs/profile_fast_breakdown_round39_20260923_*.json`.

The table feasibility screen stored the combined transfer coefficients for
`z∈[-1,5.1]` at 4097 nodes and 32 phase-substep rows, using about 2.1 MiB.
Cubic interpolation still reached only `9.73e-10` maximum absolute transfer
error, above the `1e-11` local prerequisite. It is
**REJECTED_FOR_PRODUCTION** before kernel timing; no table-density or
interpolation-order retuning is allowed.
Artifact: `docs/transfer_table_round39_20260923.json`.

The next selected mathematical direction is a Levin-collocation
nonoscillatory phase feasibility prototype. It is distinct from the rejected
arbitrary-basis Kummer fit and right-boundary-WKB solve, and must first pass
residual/transfer gates before any performance claim.

## Round 40 result — Levin interaction-picture envelope (2026-09-23)

Fresh profiles at `1e88d3161b32fa780bee10546eb4b1e42d7f84eb` used the fixed
2-thread/workqueue/BLAS contract, 25 repeats, goal/PCHIP and `kink_split=true`.
Total/tensor medians were `5.274/2.277` default, `5.192/2.760` low-T,
`7.952/4.262` high-T, `8.362/4.698` stiff and `8.106/4.330` high-kappa ms.
Artifacts: `docs/profile_fast_breakdown_round40_20260923_*.json`.

The standalone prototype factored the fast WKB carrier from the transformed
scalar equation and integrated the interaction-picture envelope, then fitted
its complex fundamental matrix with degree-12 Chebyshev coefficients. Across
the named regimes, envelope fit relative errors were approximately
`1.36e10–1.97e10`, and independent scalar-transfer relative errors were
`0.0058–0.0562`. The carrier transformation did not produce a slowly
compressible representation, so this implementation is
**REJECTED_FOR_PRODUCTION** before timing.
Artifact: `docs/levin_envelope_round40_20260923.json`.

Next selected direction is a fixed-width four-mode lockstep prototype with
explicit unrolled lanes. It is only a candidate if fresh LLVM confirms real
vector lanes or a materially different scheduling effect; the prior indirect
grouped SoA/AoSoA result remains negative knowledge.
