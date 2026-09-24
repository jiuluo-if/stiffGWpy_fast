# stiffgwpy_fast — 文档索引

[English](README.md) | 中文说明

状态：当前文档导航

日期：2026-09-24

代码版本：见 `docs/validation/validation_manifest.json` 中的 `commit`

## 研究与实验导航

先看[实验目录](experiment_catalog_zh.md)，按主题查找研究线索、当前证据状态、报告和原始实验数据。目录会区分已接受的发现、未通过的候选方法、探索性运行，以及尚未验证的发布结论。

原始 JSON/JSONL 文件是单次运行记录。实验目录按文件名和文件夹归纳相近轮次，避免把 smoke 测试或候选方案的性能分析误当作正式生产结果。请保留记录中的数据路径，因为实验输出和验证清单会引用这些路径。

根目录 README 是面向使用者的主指南。本目录的文档详细说明两种 fast 使用档位。所有精度数值均从 `docs/validation/validation_manifest.json` 读取；该清单由已提交的验证产物（`docs/paramsweep_*`、`docs/mcmc_posterior/`）通过 `scripts/build_two_mode_manifest.py` 只读重放生成，不会重新运行物理计算。

后验模拟使用的独立参考真值保存在 `mcmc_posterior/oracle_truth.json`。它是必需的测试数据，不是临时基准结果。

| 文档 | 内容 | 状态 |
|---|---|---|
| `physics.md` | 背景演化、张量模和 `Delta N_eff` 闭合关系 | 当前 |
| `numerical_method.md` | fast 双档位方案、转变区细化、尾部处理和自适应网格 | 当前 |
| `accuracy.md` | 分层精度边界和有日期的验证证据 | 以 manifest 和实验目录为准 |
| `parameter_validation.md` | 参数定义和带日期的扫描结果 | 以 manifest 和实验目录为准 |
| `cobaya.md` | 适配器选项、模式映射和 `eval_freqs` | 当前 |
| `benchmarks.md` | 当前冷启动/预热运行时间及方法对比入口 | 当前 |
| `performance_comparison_20260903.md` | 优化前后对比、耗时分解、A/B 和门槛 | 当前参考 |
| `reproducibility.md` | 驱动脚本、验收门槛和环境元数据 | 当前 |
| `experiment_catalog_zh.md` | 研究主题、证据状态和原始实验数据家族 | 当前 |
| `archive/` | 已被替代的快照和仅保存在本机的诊断日志 | 历史 |

2026-09-11 的 `baseline_54d65e3.md` 双线程快照已保存在[归档目录](archive/baselines/baseline_54d65e3.md)。它记录早期基线，不代表当前发布结论。
