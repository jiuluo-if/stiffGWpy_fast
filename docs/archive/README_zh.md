# 历史文档与本地实验数据

[English](README.md) | 中文说明

此目录保留不再作为当前报告或基准输入的历史材料。README 和实验目录仍会链接部分文件以便查阅历史记录；这些链接不表示它们是当前证据。

- `baselines/baseline_54d65e3.md` 是 2026-09-11 的双线程快照。2026-09-12 的四区域精度记录见 `../fast_v02_audit_report.md`；较新的两处边界 Oracle C 证据见 `../accuracy.md` 和 `../oracle_c_wkb_assessment.md`。`../validation/validation_manifest.json` 日期为 2026-09-17。上述文件均不代表当前 HEAD 已通过全参数空间认证。
- `oracle_prufer/raw_logs/` 是仅保存在本机的归档目录，存放从 `docs/` 根目录移出的三份 2026-09-11 诊断记录。仓库的 `*.log` 忽略规则会阻止这些原始日志进入 Git；后续 JSON 记录和评估报告仍是共享、可导航的证据。新克隆不会包含这些被忽略的日志。
- [`reference/head_20260910/benchmark_default.jsonl`](reference/head_20260910/benchmark_default.jsonl)
  记录了基于源码提交 `767056d` 的单个默认参数点 fast/reference 对比。
  代码和基准脚本都不读取它。该记录不同于后续
  `baselines/baseline_54d65e3.md` 双线程报告，因此作为历史单点证据保留，不能视为当前性能基线。

归档材料只代表历史证据。除非核对其记录的代码提交和实验方案，否则不要将其视为当前生产结论。
