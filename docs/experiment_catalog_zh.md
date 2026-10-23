# 实验目录

[English catalogue](experiment_catalog.md) | 中文目录

状态：研究主题和运行产物导航

最近一次完整盘点：2026-10-22；2026-10-23 定向补核了本地 Prüfer 原始日志的来源信息。

代码版本：每个实验文件均绑定自己的 commit；validation manifest 日期为 `2026-09-17`，较新的 profile 和报告需按各自运行提交号读取。

本目录按科学问题归纳 `docs/` 中的实验材料。Markdown 报告说明假设、实验方案和决策；JSON/JSONL 文件保留测量值。某个候选方案在少量点上更快或数值相近，并不表示它已进入生产路径。

## 如何阅读证据

1. 先看[`benchmarks.md`](benchmarks.md)了解最新固定资源全流程耗时，并看[`accuracy.md`](accuracy.md)了解当前有限范围的精度证据。`fast_v02_audit_report.md` 和 `validation/validation_manifest.json` 是带日期的基线，阅读时需同时留意其覆盖缺口。
2. 再读下表中的主题报告，确认假设、对照、否定结果，以及接受或拒绝方法的原因。
3. 如需核验逐点数值、运行环境、代码提交、重复次数和重放信息，再查看对应的 JSON/JSONL 数据。
4. `SMOKE`、`SPIKE` 和 `NOT VERIFIED` 表示探索性或未完成结果。失败、显式物理保护和被拒绝的候选方案也是研究证据；除非报告如此分类，不要将其记作求解器数值故障。

## 研究主题

| 主题 | 研究问题与归纳 | 报告和数据 |
|---|---|---|
| 物理模型与求解器契约 | 定义宇宙学模型、输出量、求解档位、精度限制和 Cobaya 适配器。此组文档说明接口和物理约定，不代表整个参数空间都已通过认证。 | [`physics.md`](physics.md)、[`numerical_method.md`](numerical_method.md)、[`accuracy.md`](accuracy.md)、[`cobaya.md`](cobaya.md)、[`reproducibility.md`](reproducibility.md) |
| 参数空间验证与独立参照 | 将成功对比、物理保护和数值失败分开统计。已有 400 点 LHS 屏查记录 254 个成功点、146 个显式保护点和 0 个数值失败；240 点 Sobol 运行记录 212 个成功点和 28 个保护点。这些是带日期的抽样，不是当前 fast 档位的全参数空间认证。 | [`parameter_validation.md`](parameter_validation.md)、[`parameter_validation/parameter_validation_report.md`](parameter_validation/parameter_validation_report.md)、[`paramsweep_plain/validation_summary.md`](paramsweep_plain/validation_summary.md)、[`paramsweep_z8/validation_summary.md`](paramsweep_z8/validation_summary.md)、[`paramsweep_z8b/validation_summary.md`](paramsweep_z8b/validation_summary.md)、[`parameter_validation/`](parameter_validation/)、[`paramsweep_plain/`](paramsweep_plain/)、[`paramsweep_z8/`](paramsweep_z8/)、[`paramsweep_z8b/`](paramsweep_z8b/)、[`paramsweep_ref/`](paramsweep_ref/)、[`validation/`](validation/) |
| `DN_gw` 误差预算与频率求积 | 较早的 A/B 实验发现，SciPy PCHIP 精度更高，但开始时未达到预注册的耗时门槛。共享和向量化减少了 PCHIP 开销；之后以 WKB 为参照的误差分解指出频率求积是主要残差来源，后续的默认方案评估记录了切换到 PCHIP 的决策。2026-09-12 的四区域审计是历史基线；较新的两处 `DN_re` 边界精度扫描和转折细化结果见 `accuracy.md` 与 `oracle_c_wkb_assessment.md`。 | [`fast_true_error_assessment.md`](fast_true_error_assessment.md)、[`fast_residual_decomposition_assessment.md`](fast_residual_decomposition_assessment.md)、[`fast_quadrature_ab_assessment.md`](fast_quadrature_ab_assessment.md)、[`fast_quadrature_reuse_assessment.md`](fast_quadrature_reuse_assessment.md)、[`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md)、[`fast_v02_audit_report.md`](fast_v02_audit_report.md)、`quadrature_*`、`fast_quadrature_*` 和 `dn_*` JSON 数据 |
| 独立参照与尾部交接 | Oracle C 的高阶 WKB 修正在五个抽样全网格点上与深层参照吻合。正式 fast 求解器在其中三个点通过 `2e-4` 精度门限，但在两个 `DN_re` 边界点误差为 `2.5e-3..2.7e-3`。使用匹配 `DN_eff` 和频率网格重跑 Oracle 后，`h=0.0003125` 和转折细化候选均通过这两个点的门限。半径 `0.125` 的扫描保持 `0.0003125` 局部步长，在 `h=0.0025`、`0.005` 和 `0.02` 下都通过；最快的通过项是 `h=0.02/factor=64`，误差为 `1.16e-4` 和 `8.25e-5`，中位耗时 `7.63/9.49 ms`，仍比正式 fast 慢约 `1.34/2.09` 倍。继续将 `h` 加粗至 `0.04/factor=64` 后两点均未通过门限。候选都未推广；这些仍是抽样证据，不等于全参数空间认证。Oracle B 的相位平均和有限窗口原型尚未推广，因为它们还不能提供独立真值参照。Prüfer 保留为独立参考原型。 | [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md)、[`accuracy.md`](accuracy.md)、[`oracle_c_wkb_corner_20261022/`](oracle_c_wkb_corner_20261022/)、[`oracle_b_phase_averaged_assessment.md`](oracle_b_phase_averaged_assessment.md)、[`oracle_b_phase_window_assessment.md`](oracle_b_phase_window_assessment.md)、[`oracle_prufer_assessment.md`](oracle_prufer_assessment.md)、`oracle_*`、`tail_*` 和 `wkb_*` JSON 数据 |
| 数值方法候选方案 | 线性 z Magnus transfer-map 候选方案因运行时间高于 midpoint 基线且未达到实验门槛，被拒绝用于生产。其他 transfer、phase、Riccati、Bessel、Airy、稀疏频率和小角度记录均为具体候选方案的证据；阅读时须同时核对其状态和代码提交。 | [`transfer_map_phase_assessment.md`](transfer_map_phase_assessment.md)、`*transfer*`、`*phase*`、`*riccati*`、`*bessel*`、`*airy*`、`*frequency*` 和 `*small_angle*` JSON 数据 |
| 性能与热点实验 | 基准文档区分预热耗时、冷启动/JIT 开销、线程缩放和精度门槛。修正后的 2026-10-22 formal `fast` profile 使用求解器提交 `6707572e73d576a0ccf294172579cedc26e5e394`、预设容差 `tol=1e-6`，并对 default、high-T、stiff、high-kappa 各测 25 条 warm 样本。2026-10-21 的 `cdbdb387479e11cd87a38d76a0e58348ca13e70f` profile 漏传底层 `tol`，实际使用 `1e-4`，因此作为较低容差的历史数据保留。两组数据都单独记录首条进程调用，详见 `benchmarks.md`。同日 PTA 似然缓存实验测量了 IPTA 和 NANOGrav 适配器的单次直接调用耗时，并在固定输入和每个适配器 10 组变化参数上得到完全相同的输出；这些数据不代表完整 MCMC 加速。2026-10-21 的张量采样索引优化以提交 `68b17ac` 为基线，在六个场景交替运行 80 轮；原始记录和复跑脚本见 `benchmarks.md`。另有 TBB 双调用吞吐对照，见同页；该结果不代表单次求解加速。2026-10-22 的参考频谱求积实验在合成谱和 formal-fast 默认谱输入上，将 `integrate_spectrum` 的函数调用数减少约 25%，并修正偶数长度粗网格遗漏末端节点的问题；计时范围仅为积分函数，不代表求解器整体加速。截至 2026-10-22，受 Git 跟踪的 `profile_fast_breakdown_*.json` 共 228 个文件，记录不同轮次和参数区域的分阶段分析，不代表 228 个独立生产改进。`derived_param_*`、`phi_s2_*`、`fd_lookup_*`、`grouped_soa_*` 和 `outer_snapshot_*` 是缓存、工作区和内核候选方案的多轮实验；接受性能结论前要比较摘要/状态一致性和全部参数区域。 | [`benchmarks.md`](benchmarks.md)、[`benchmark_pta_likelihood_cache_20261021.json`](benchmark_pta_likelihood_cache_20261021.json)、[`benchmark_pta_likelihood_cache.py`](../scripts/benchmark_pta_likelihood_cache.py)、[`benchmark_reference_quad_breakpoints_20261022.json`](benchmark_reference_quad_breakpoints_20261022.json)、[`benchmark_reference_quad_breakpoints.py`](../scripts/benchmark_reference_quad_breakpoints.py)、[`performance_comparison_20260903.md`](performance_comparison_20260903.md)、[`baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md)、`benchmark_*`、`profile_*`、`derived_param_*`、`phi_s2_*`、`fast_phi_*`、`fd_lookup_*`、`grouped_soa_*`、`outer_snapshot_*` 和 `count_assembly_*` 数据 |
| 后验与外部模型对比 | 2026-09-24 的报告在同一 LVK 数据下比较三种 MCMC 方法。当前 fast 的采样单步约比 SageNet+ 快 2.9–4.9 倍；这衡量采样步骤耗时，不是独立求解器测速。低再加热温度场景对三种方法的 LVK 频段覆盖均为 0%，不能据此比较数据拟合表现。报告正式链文件包含每种方法/情景 4 条链、每条保留 10,000 个样本；`.gitignore` 下另有一份本机 3 条链 × 2,000 样本的较早独立运行，尚未记录完整来源，暂保留待核，不作为正式报告证据。 | [`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md)、[`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md)、[`mcmc_sagenet_compare/`](mcmc_sagenet_compare/)、[`mcmc_posterior/`](mcmc_posterior/) |
| 工程、测试与历史计划 | 这组是维护记录，不是物理结果报告。工程审计仍被根 README、变更日志和复现指南引用。2026-09-10 的 `DN_gw`/速度计划已成为历史计划；独立参照的同网格误差和稳定低于 4 ms 的速度门槛，在当前证据中仍未通过。2026-09-23 的张量传播计划也已成为历史记录：精确指数复用候选未被采用，后续 `next_assembly_k` 实现及测量结果见 2026-10-21 基准。 | [`engineering_audit.md`](engineering_audit.md)、[`test_coverage_matrix.md`](test_coverage_matrix.md)、[`test_duplication_audit.md`](test_duplication_audit.md)、[`superpowers/plans/2026-09-10-fast-v02-dn-gw.md`](superpowers/plans/2026-09-10-fast-v02-dn-gw.md)、[`superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`](superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md)、[`accuracy.md`](accuracy.md)、[`benchmark_kernel_assembly_stride_20261021.json`](benchmark_kernel_assembly_stride_20261021.json) |

## 2026-10-22 转折细化路径的 `H2_vec` 优化

转折细化路径过去会对每个网格节点单独调用 SciPy 密度插值器。批量求值后，所有求解输出数组与旧实现逐位一致；在两个 `DN_re` 边界点，预热后的中位耗时缩短了 `3.86–4.09` 倍。这只是在验证用转折细化路径上的实现速度改进，不改变正式精确 kink 档位，也不增加精度覆盖。详见[基准汇总](benchmarks.md#transition-refined-h2_vec-batching-2026-10-22)、[原始计时记录](benchmark_h2_vectorization_20261022.json)和[复跑脚本](../scripts/benchmark_h2_vectorization.py)。

## 2026-10-22 局部细化匹配扫描

三个半径 `0.125` 的配置都将转折区局部步长保持在 `0.0003125`，同时逐步加粗区外背景网格。`h=0.02/factor=64` 在抽查的两个 `DN_re` 边界上通过 `2e-4` 门限，中位耗时为 `7.63/9.49 ms`，但仍比正式 fast 慢约 `1.34/2.09` 倍。将 `h` 提高到 `0.04` 后耗时约降低三分之一，却在两条边界上都未通过精度门限。没有候选被推广。实验方案、匹配参照差异和原始记录见[Oracle C 报告](oracle_c_wkb_assessment.md)及链接的 `oracle_c_wkb_corner_20261022/transition_h*` 目录。

后续将 `h=0.04` 的半径增至 `0.25`：倍数 64 时两处都未通过；倍数 128 虽恢复 `0.0003125` 的局部步长，但低边界误差仍为 `4.51e-4`，仅高边界通过。倍数 128 的重复中位耗时为 `47.15/7.75 ms`，两点差异明显。两种配置均不推广，详见[Oracle C 报告](oracle_c_wkb_assessment.md)。

另将通过门槛的 `h=0.02/factor=64` 半径从 `0.125` 缩至 `0.0625`，两处误差降至 `1.14e-4/8.30e-5`，分次运行中位耗时为 `6.49/8.74 ms`。它仍慢于正式 fast，且尚未与其他半径做固定资源配对比较，因此只作为验证候选。

半径继续减至 `0.03125` 后，两处仍通过（`1.12e-4/8.29e-5`），中位耗时为 `6.41/8.09 ms`。这点耗时差来自不同的未锁定 CPU 运行，尚不能证明有稳定提速，仍不推广。

另将半径 `0.015625` 与 `0.0078125` 固定在 CPU `[0,1]` 上交错比较，每个边界测 30 组。更小半径的配对中位耗时比分别为 `0.989/0.979`，但 95% 区间都覆盖 `1.0`，因此没有证据证明它更快。原始样本见[配对计时记录](benchmark_transition_radius_paired_20261022.json)。

## 原始运行数据家族

| 位置或文件名前缀 | 内容 | 处理方式 |
|---|---|---|
| `paramsweep_*`、`parameter_validation/`、`validation/` | 同网格、参数空间、收敛性和参照独立性记录 | 保留：被报告、manifest 或复现说明引用 |
| `oracle_a_*`、`oracle_b_*`、`oracle_prufer_*`、`oracle_tail_*`、`tail_*`、`wkb_*` | 同网格 Oracle A/B/C、Prüfer 参考对比、相位/尾部敏感性、WKB 修正和诊断子集 | 作为有限范围的独立参照与误差归属证据保留；Prüfer 仍仅作参考，任何产物都不代表全参数空间认证 |
| `phase_recurrence_oracle_20260917/` | 默认、高 `kappa`、高温、低温和 stiff 五个情景的 Oracle C WKB 记录，每个文件仅含 12 个抽样频率，生成源码为提交 `58af509` | 作为稀疏相位递推敏感性子集保留；它不是 native 全网格对比，也不构成全参数空间认证 |
| `validation_h005_round_20260916/`、`validation_h0075_round_20260916/`、`validation_h0075_round_20260916b/` | 源码提交 `1bf840d` 下三组带日期的 10 点步长扫描；首轮 `h=0.0075` 的 10 点均为 `fast_exception`，完成的 `h=0.005` 与重跑 `h=0.0075` 均未通过信号、转折区和 `DN_gw` 验收门槛 | 作为被拒绝的历史证据保留；完成扫描的最大信号相对误差约 `6.78e-3`，最大 `DN_gw` 相对误差约 `3.60e-3`，不能代表当前精度 |
| `profile_fast_breakdown_*.json`、`benchmark_*.json` | 各轮性能、线程、稳定性和阶段耗时测量 | 作为带日期的运行数据保留；比较前核对对应报告和代码提交 |
| 其他 `profile_*.json`、`*_round*.json`、`*_head*.json`、`*_candidate*.json`、`*_spike*.json` | 各轮 profile 和独立实现候选；包含被拒绝或仅作参考的方案 | 原始记录随对应报告保留；状态以报告为准，不能按文件名或日期推断 |
| `ab_baseline_njit_*.json`、`ab_candidate_njit_*.json`、`production_ab_*.txt` | 2026-09-03 Numba 并行内核修改的基线/候选计时及 production 输出指纹；线程缩放仅保留在报告汇总中，原始 `scaling_njit_t*.json` 已由清理提交 `2619065` 移除 | 剩余原始记录随[`performance_comparison_20260903.md`](performance_comparison_20260903.md)保留；它们和报告汇总均不是当前性能基线 |
| `fast_true_error_*`、`fast_residual_*`、`error_budget_probe_*`、`fast_quadrature_*`、`quadrature_*` | 带日期的 Oracle 误差归属、误差估计校准和求积对比 | 作为有限范围精度证据保留；当前解释以 `accuracy.md` 和关联评估报告为准 |
| `derived_param_*.json`、`fd_lookup_*.json`、`grouped_soa_*.json`、`outer_snapshot_*.json` | 缓存、工作区、查表及外层求解候选方案 | 属于探索性结果；比较状态/摘要一致性和参数区域覆盖 |
| `gen_fast_*`、`h2_*`、`outer_*`、`kernel_*`、`prep_*`、`probe_*`、`canonical_fast_*`、`fast_pyoverhead_*` | 网格、背景、外层求解、内核和 Python 开销的分轮候选实验 | 保留已接受和被拒绝的对照；引用计时前核对评估报告及资源元数据 |
| `phi_s2_*`、`fast_phi_*`、`phase_*` | 相位和 `phi_s2` 实现候选方案，包括 smoke 运行 | 属于探索性结果；不能根据 smoke 结果推断已进入生产路径 |
| `*transfer*`、`*wkb*`、`*riccati*`、`*bessel*`、`*airy*`、`*small_angle*`、`*frequency*` | 数值方法原型、参照研究和后续实验轮次 | 按每条记录的状态保留；被拒绝和仅作参考的结果也是研究过程的一部分 |
| `mcmc_posterior/`、`mcmc_sagenet_compare/`、`mcmc/` | 后验验证、方法对比报告、图和本地采样链 | 保留当前报告与输入；大型/本地链数据仍遵循 `.gitignore`。正式对比链为 `mcmc_sagenet_compare/posterior_chains_20260924.npz`（每种方法/情景 4 × 10,000）；`mcmc/chains/sagenet_compare_20260924.npz` 是不同的 3 × 2,000 运行，不是重复副本。先留在本机，待补齐运行设置和来源后再决定是否纳入正式结论。 |

## 归档与清理记录

- [`paramsweep_ref/validation_summary.json`](paramsweep_ref/validation_summary.json) 是基于提交 `f87e969`、参数 `z_tail=7`、`rtol=1e-9`、`freq_res=2` 生成的 10 点汇总；信号和积分 `DN_gw` 的验收门槛均未通过。代码和基准脚本都不读取该文件，逐点输入仍保存在同目录的 JSONL 中。该文件作为这次未通过/未完成扫描的简明记录保留，不代表当前认证结果。
- [`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md) 是 2026-09-11 的早期双线程快照。求解器和基准脚本均不读取它；README 和目录中的链接仅用于查阅历史记录。该文件保留为当次运行的证据，不代表当前基线。
- [`archive/reference/head_20260910/benchmark_default.jsonl`](archive/reference/head_20260910/benchmark_default.jsonl) 是基于提交 `767056d` 的独立单点 fast/reference 记录；代码和基准脚本都不读取。它保留为带日期的单点证据，不替代后续双线程快照，也不代表当前基线。
- 三份日期为 2026-09-11 的 Prüfer 原始 `.log` 记录已从 `docs/` 根目录移至 `archive/oracle_prufer/raw_logs/`；每份日志开头的 JSON 记录都标明来源提交为 `91598e7`。它们被仓库通用的 `*.log` 规则忽略，且没有活动代码或基准脚本读取；后续 JSON 数据和评估报告是共享导航入口。由于这些记录来自较早提交且与后续 JSON 并非逐字节相同，仍在本机保留，不删除。新克隆不会包含这些本地日志，见[`archive/README.md`](archive/README.md)。
- 未删除 2026-09-02/03 的参数验证或后验数据：当前验证/报告链仍引用这些数据。也没有仅因为轮次较新或日期较早而删除当前实验输出。
- 2026-09-25 复查发现被 `.gitignore` 排除的 `mcmc/chains/sagenet_compare_20260924.npz`（3 条链 × 每条 2,000 个样本）与报告归档的 4 × 10,000 链不同；报告和脚本均未引用前者。由于它是近期独立运行且设置来源尚未查清，暂留本机并排除在正式对比结论之外。同次查漏还确认旧的 `mcmc_smoke2/chains/` 空目录无文件、历史或引用，已于 2026-09-25 清理。

## 版本控制内的 Markdown 文档清单

以下清单覆盖当前 `docs/` 中的研究说明和配套文档；各主题的关联方式见上表。

- 核心说明：`physics.md`、`numerical_method.md`、`accuracy.md`、`parameter_validation.md`、`cobaya.md`、`benchmarks.md`、`reproducibility.md`。
- 精度与求积：`fast_v02_audit_report.md`、`fast_true_error_assessment.md`、`fast_residual_decomposition_assessment.md`、`fast_quadrature_ab_assessment.md`、`fast_quadrature_reuse_assessment.md`、`fast_quadrature_default_switch_assessment.md`、`performance_comparison_20260903.md`。
- 参照方法：`oracle_b_phase_averaged_assessment.md`、`oracle_b_phase_window_assessment.md`、`oracle_c_wkb_assessment.md`、`oracle_prufer_assessment.md`、`transfer_map_phase_assessment.md`。
- 参数/后验：`parameter_validation/parameter_validation_report.md`、`paramsweep_plain/validation_summary.md`、`paramsweep_z8/validation_summary.md`、`paramsweep_z8b/validation_summary.md`、`mcmc_posterior/posterior_validation.md`、`mcmc_sagenet_compare/report.md`。
- 工程和计划：`engineering_audit.md`、`test_coverage_matrix.md`、`test_duplication_audit.md`、`superpowers/plans/2026-09-10-fast-v02-dn-gw.md`、`superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`。
- 导航和历史：`README.md`、`README_zh.md`、`experiment_catalog.md`、`experiment_catalog_zh.md`、`archive/README.md`、`archive/README_zh.md`、`archive/baselines/baseline_54d65e3.md`。
