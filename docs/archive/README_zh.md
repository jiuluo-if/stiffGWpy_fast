# 历史文档与本地实验数据

[English](README.md) | 中文说明

此目录保留已经不再承担当前报告或主要导航作用、但仍有溯源价值的材料。

- `baselines/baseline_54d65e3.md` 是 2026-09-11 的双线程快照。当前验证状态请看 `../fast_v02_audit_report.md` 和 `../validation/validation_manifest.json`。
- `oracle_prufer/raw_logs/` 是仅保存在本机的归档目录，存放从 `docs/` 根目录移出的三份 2026-09-11 诊断记录。仓库的 `*.log` 忽略规则会阻止这些原始日志进入 Git；后续 JSON 记录和评估报告仍是共享、可导航的证据。新克隆不会包含这些被忽略的日志。
- `reference/head_20260910/` 保存较早的基准参考快照。

归档材料只代表历史证据。除非核对其记录的代码提交和实验方案，否则不要将其视为当前生产结论。
