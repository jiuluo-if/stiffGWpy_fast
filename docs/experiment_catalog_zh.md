# 实验目录

[English catalogue](experiment_catalog.md) | 中文目录

状态：截至 2026-09-24 的研究主题和运行产物导航

代码版本：每个实验文件均绑定自己的 commit；validation manifest 日期为 `2026-09-17`，较新的 profile 和报告需按各自运行提交号读取。

本目录按科学问题归纳 `docs/` 中的实验材料。Markdown 报告说明假设、实验方案和决策；JSON/JSONL 文件保留测量值。某个候选方案在少量点上更快或数值相近，并不表示它已进入生产路径。

## 如何阅读证据

1. 先看[`benchmarks.md`](benchmarks.md)了解最新固定资源全流程耗时，再把[`fast_v02_audit_report.md`](fast_v02_audit_report.md)和[`validation/validation_manifest.json`](validation/validation_manifest.json)作为带日期的精度证据及其明确缺口阅读。
2. 再读下表中的主题报告，确认假设、对照、否定结果，以及接受或拒绝方法的原因。
3. 如需核验逐点数值、运行环境、代码提交、重复次数和重放信息，再查看对应的 JSON/JSONL 数据。
4. `SMOKE`、`SPIKE` 和 `NOT VERIFIED` 表示探索性或未完成结果。失败、显式物理保护和被拒绝的候选方案也是研究证据；除非报告如此分类，不要将其记作求解器数值故障。

## 研究主题

| 主题 | 研究问题与归纳 | 报告和数据 |
|---|---|---|
| 物理模型与求解器契约 | 定义宇宙学模型、输出量、求解档位、精度限制和 Cobaya 适配器。此组文档说明接口和物理约定，不代表整个参数空间都已通过认证。 | [`physics.md`](physics.md)、[`numerical_method.md`](numerical_method.md)、[`accuracy.md`](accuracy.md)、[`cobaya.md`](cobaya.md)、[`reproducibility.md`](reproducibility.md) |
| 参数空间验证与独立参照 | 将成功对比、物理保护和数值失败分开统计。已有 400 点 LHS 屏查记录 254 个成功点、146 个显式保护点和 0 个数值失败；240 点 Sobol 运行记录 212 个成功点和 28 个保护点。这些是带日期的抽样，不是当前 fast 档位的全参数空间认证。 | [`parameter_validation.md`](parameter_validation.md)、[`parameter_validation/parameter_validation_report.md`](parameter_validation/parameter_validation_report.md)、[`paramsweep_plain/validation_summary.md`](paramsweep_plain/validation_summary.md)、[`paramsweep_z8/validation_summary.md`](paramsweep_z8/validation_summary.md)、[`paramsweep_z8b/validation_summary.md`](paramsweep_z8b/validation_summary.md)、[`parameter_validation/`](parameter_validation/)、[`paramsweep_plain/`](paramsweep_plain/)、[`paramsweep_z8/`](paramsweep_z8/)、[`paramsweep_z8b/`](paramsweep_z8b/)、[`paramsweep_ref/`](paramsweep_ref/)、[`validation/`](validation/) |
| `DN_gw` 误差预算与频率求积 | 较早的 A/B 实验发现，SciPy PCHIP 精度更高，但开始时未达到预注册的耗时门槛。共享和向量化减少了 PCHIP 开销；之后以 WKB 为参照的误差分解指出频率求积是主要残差来源，后续的默认方案评估记录了切换到 PCHIP 的决策。最新的详细精度审计日期为 2026-09-12，提供的是有限范围的历史证据。 | [`fast_true_error_assessment.md`](fast_true_error_assessment.md)、[`fast_residual_decomposition_assessment.md`](fast_residual_decomposition_assessment.md)、[`fast_quadrature_ab_assessment.md`](fast_quadrature_ab_assessment.md)、[`fast_quadrature_reuse_assessment.md`](fast_quadrature_reuse_assessment.md)、[`fast_quadrature_default_switch_assessment.md`](fast_quadrature_default_switch_assessment.md)、[`fast_v02_audit_report.md`](fast_v02_audit_report.md)、`quadrature_*`、`fast_quadrature_*` 和 `dn_*` JSON 数据 |
| 独立参照与尾部交接 | Oracle C 的高阶 WKB 修正已在报告列出的对比范围内验证，并解释了冻结尾部造成的缺陷。Oracle B 的相位平均和有限窗口原型尚未推广，因为它们还不能提供独立真值参照。Prüfer 作为独立参考原型通过了已记录的筛查，但在更广泛认证完成前仍只用于参考。 | [`oracle_c_wkb_assessment.md`](oracle_c_wkb_assessment.md)、[`oracle_b_phase_averaged_assessment.md`](oracle_b_phase_averaged_assessment.md)、[`oracle_b_phase_window_assessment.md`](oracle_b_phase_window_assessment.md)、[`oracle_prufer_assessment.md`](oracle_prufer_assessment.md)、`oracle_*`、`tail_*` 和 `wkb_*` JSON 数据 |
| 数值方法候选方案 | 线性 z Magnus transfer-map 候选方案因运行时间高于 midpoint 基线且未达到实验门槛，被拒绝用于生产。其他 transfer、phase、Riccati、Bessel、Airy、稀疏频率和小角度记录均为具体候选方案的证据；阅读时须同时核对其状态和代码提交。 | [`transfer_map_phase_assessment.md`](transfer_map_phase_assessment.md)、`*transfer*`、`*phase*`、`*riccati*`、`*bessel*`、`*airy*`、`*frequency*` 和 `*small_angle*` JSON 数据 |
| 性能与热点实验 | 基准文档区分预热耗时、冷启动/JIT 开销、线程缩放和精度门槛。`profile_fast_breakdown_*.json` 是最大的原始数据家族（本次 2026-09-24 清点到 220 条），表示不同轮次和参数区域的分阶段分析，不代表 220 个独立生产改进。`derived_param_*`、`phi_s2_*`、`fd_lookup_*`、`grouped_soa_*` 和 `outer_snapshot_*` 是缓存、工作区和内核候选方案的多轮实验；接受性能结论前要比较摘要/状态一致性和全部参数区域。 | [`benchmarks.md`](benchmarks.md)、[`performance_comparison_20260903.md`](performance_comparison_20260903.md)、[`baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md)、`benchmark_*`、`profile_*`、`derived_param_*`、`phi_s2_*`、`fast_phi_*`、`fd_lookup_*`、`grouped_soa_*`、`outer_snapshot_*` 和 `count_assembly_*` 数据 |
| 后验与外部模型对比 | 2026-09-24 的报告在同一 LVK 数据下比较三种 MCMC 方法。当前 fast 的采样单步约比 SageNet+ 快 2.9–4.9 倍；这衡量采样步骤耗时，不是独立求解器测速。低再加热温度场景对三种方法的 LVK 频段覆盖均为 0%，不能据此比较数据拟合表现。报告正式链文件包含每种方法/情景 4 条链、每条保留 10,000 个样本；`.gitignore` 下另有一份本机 3 条链 × 2,000 样本的较早独立运行，尚未记录完整来源，暂保留待核，不作为正式报告证据。 | [`mcmc_sagenet_compare/report.md`](mcmc_sagenet_compare/report.md)、[`mcmc_posterior/posterior_validation.md`](mcmc_posterior/posterior_validation.md)、[`mcmc_sagenet_compare/`](mcmc_sagenet_compare/)、[`mcmc_posterior/`](mcmc_posterior/) |
| 工程和测试覆盖 | 这组是维护审计，不是物理结果报告。工程审计仍被根 README、变更日志和复现指南引用，因此保留在当前文档树中。 | [`engineering_audit.md`](engineering_audit.md)、[`test_coverage_matrix.md`](test_coverage_matrix.md)、[`test_duplication_audit.md`](test_duplication_audit.md) |

## 原始运行数据家族

| 位置或文件名前缀 | 内容 | 处理方式 |
|---|---|---|
| `paramsweep_*`、`parameter_validation/`、`validation/` | 同网格、参数空间、收敛性和参照独立性记录 | 保留：被报告、manifest 或复现说明引用 |
| `oracle_prufer_fullgrid*.json`、`oracle_prufer_*_sobol*.json` | Prüfer 完整本征频率网格对比及边界/Sobol 样本 | 保留：支撑“仅作参考”的认证边界 |
| `profile_fast_breakdown_*.json`、`benchmark_*.json` | 各轮性能、线程、稳定性和阶段耗时测量 | 作为带日期的运行数据保留；比较前核对对应报告和代码提交 |
| `derived_param_*.json`、`fd_lookup_*.json`、`grouped_soa_*.json`、`outer_snapshot_*.json` | 缓存、工作区、查表及外层求解候选方案 | 属于探索性结果；比较状态/摘要一致性和参数区域覆盖 |
| `phi_s2_*`、`fast_phi_*`、`phase_*` | 相位和 `phi_s2` 实现候选方案，包括 smoke 运行 | 属于探索性结果；不能根据 smoke 结果推断已进入生产路径 |
| `*transfer*`、`*wkb*`、`*riccati*`、`*bessel*`、`*airy*`、`*small_angle*`、`*frequency*` | 数值方法原型、参照研究和后续实验轮次 | 按每条记录的状态保留；被拒绝和仅作参考的结果也是研究过程的一部分 |
| `mcmc_posterior/`、`mcmc_sagenet_compare/`、`mcmc/` | 后验验证、方法对比报告、图和本地采样链 | 保留当前报告与输入；大型/本地链数据仍遵循 `.gitignore`。正式对比链为 `mcmc_sagenet_compare/posterior_chains_20260924.npz`（每种方法/情景 4 × 10,000）；`mcmc/chains/sagenet_compare_20260924.npz` 是不同的 3 × 2,000 运行，不是重复副本。先留在本机，待补齐运行设置和来源后再决定是否纳入正式结论。 |

## 归档与清理记录

- [`archive/baselines/baseline_54d65e3.md`](archive/baselines/baseline_54d65e3.md) 是早期双线程快照，日期为 2026-09-11，仓库内没有活动路径引用。它作为历史证据保留，不再作为当前基线展示。
- 三份日期为 2026-09-11 的 Prüfer 原始 `.log` 记录已从 `docs/` 根目录移至 `archive/oracle_prufer/raw_logs/`。它们被仓库通用的 `*.log` 规则忽略，且没有活动路径引用；后续 JSON 数据和评估报告是共享导航入口。由于这些记录来自较早提交且与后续 JSON 并非逐字节相同，仍在本机保留，不删除。新克隆不会包含这些本地日志，见[`archive/README.md`](archive/README.md)。
- 未删除 2026-09-02/03 的参数验证或后验数据：当前验证/报告链仍引用这些数据。也没有仅因为轮次较新或日期较早而删除当前实验输出。
- 2026-09-25 复查发现被 `.gitignore` 排除的 `mcmc/chains/sagenet_compare_20260924.npz`（3 条链 × 每条 2,000 个样本）与报告归档的 4 × 10,000 链不同；报告和脚本均未引用前者。由于它是近期独立运行且设置来源尚未查清，暂留本机并排除在正式对比结论之外。另发现空的 `mcmc_smoke2/chains/` 目录，无文件、历史或引用；本机自动审核拒绝了目录清理操作，因此仍留在本机，不纳入版本控制。

## 版本控制内的 Markdown 文档清单

以下清单覆盖当前 `docs/` 中的研究说明和配套文档；各主题的关联方式见上表。

- 核心说明：`physics.md`、`numerical_method.md`、`accuracy.md`、`parameter_validation.md`、`cobaya.md`、`benchmarks.md`、`reproducibility.md`。
- 精度与求积：`fast_v02_audit_report.md`、`fast_true_error_assessment.md`、`fast_residual_decomposition_assessment.md`、`fast_quadrature_ab_assessment.md`、`fast_quadrature_reuse_assessment.md`、`fast_quadrature_default_switch_assessment.md`、`performance_comparison_20260903.md`。
- 参照方法：`oracle_b_phase_averaged_assessment.md`、`oracle_b_phase_window_assessment.md`、`oracle_c_wkb_assessment.md`、`oracle_prufer_assessment.md`、`transfer_map_phase_assessment.md`。
- 参数/后验：`parameter_validation/parameter_validation_report.md`、`paramsweep_plain/validation_summary.md`、`paramsweep_z8/validation_summary.md`、`paramsweep_z8b/validation_summary.md`、`mcmc_posterior/posterior_validation.md`、`mcmc_sagenet_compare/report.md`。
- 工程和计划：`engineering_audit.md`、`test_coverage_matrix.md`、`test_duplication_audit.md`、`superpowers/plans/2026-09-10-fast-v02-dn-gw.md`、`superpowers/plans/2026-09-23-fast-v02-tensor-propagation.md`。
- 导航和历史：`README.md`、`README_zh.md`、`experiment_catalog.md`、`experiment_catalog_zh.md`、`archive/README.md`、`archive/README_zh.md`、`archive/baselines/baseline_54d65e3.md`。
