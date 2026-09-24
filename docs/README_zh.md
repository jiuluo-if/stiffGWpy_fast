# stiffgwpy_fast — 文档索引

[English](README.md) | 中文说明

状态：当前文档导航

日期：2026-09-24

代码版本：每个证据文件都记录自己的 commit；manifest 日期为 `2026-09-17`，不是当前 HEAD 的认证。

## 研究与实验导航

先看[实验目录](experiment_catalog_zh.md)，按主题查找研究线索、当前证据状态、报告和原始实验数据。目录会区分已接受的发现、未通过的候选方法、探索性运行，以及尚未验证的发布结论。

原始 JSON/JSONL 文件是单次运行记录。实验目录按文件名和文件夹归纳相近轮次，避免把 smoke 测试或候选方案的性能分析误当作正式生产结果。请保留记录中的数据路径，因为实验输出和验证清单会引用这些路径。

根目录 README 是面向使用者的主指南。当前只有一个用户可见的 `fast` 档位；历史验证中的 `plain-grid`、`production` 档位只适用于各自记录的实验，不是当前用户档位。精度数字必须结合实验 commit 和覆盖范围阅读：这里的 manifest 是较早已提交验证运行的只读汇总，不代表当前 HEAD 已通过全参数空间认证。

后验模拟使用的独立参考真值保存在 `mcmc_posterior/oracle_truth.json`。它是必需的测试数据，不是临时基准结果。

| 文档 | 内容 | 状态 |
|---|---|---|
| `physics.md` | 背景演化、张量模和 `Delta N_eff` 闭合关系 | 当前 |
| `numerical_method.md` | 唯一 fast 档位、精确 kink 拆分、尾部处理和 goal 网格 | 当前 |
| `accuracy.md` | 分层精度边界和有日期的验证证据 | 以 manifest 和实验目录为准 |
| `parameter_validation.md` | 参数定义和带日期的扫描结果 | 以 manifest 和实验目录为准 |
| `cobaya.md` | 适配器选项、模式映射和 `eval_freqs` | 当前 |
| `benchmarks.md` | 当前冷启动/预热运行时间及方法对比入口 | 当前 |
| `performance_comparison_20260903.md` | 优化前后对比、耗时分解、A/B 和门槛 | 带日期的历史快照 |
| `reproducibility.md` | 驱动脚本、验收门槛和环境元数据 | 当前 |
| `experiment_catalog_zh.md` | 研究主题、证据状态和原始实验数据家族 | 当前 |
| `archive/` | 已被替代的快照和仅保存在本机的诊断日志 | 历史 |
| `mcmc_sagenet_compare/` | 2026-09-24 LVK MCMC 速度、链诊断和后验能谱比较 | 带日期的实验 |
| `mcmc_posterior/` | 与连续 sigma reference 对比的重要性采样后验验证 | 带日期的验证 |

2026-09-11 的 `baseline_54d65e3.md` 双线程快照已保存在[归档目录](archive/baselines/baseline_54d65e3.md)。它记录早期基线，不代表当前发布结论。
