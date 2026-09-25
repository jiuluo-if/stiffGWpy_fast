# 随包提供的 Cobaya 似然

[English](README.md) | 中文说明

本程序包包含 LVK 和 PTA 似然适配器，以及对应 YAML 配置引用的数据资源。上级[Cobaya 适配器说明](../README_zh.md)介绍 theory 接口；本页说明似然类与配置、数据目录的对应关系。

| 目录 | 适配器/配置 | 随包输入 |
|---|---|---|
| `LIGO_SGWB/` | `LVK_SGWB_CC` / `LVK_SGWB_CC.yaml` | `C_O1_O2_O3.dat` |
| `PTA/` | `IPTA` / `IPTA.yaml` | `EPTAdr2/` 下的 EPTA DR2 mock 样本和频率 bin |
| `PTA/` | `NANOGrav` / `NANOGrav.yaml` | `NANOGrav15yr/` 下的密度网格和频率网格 |

YAML 文件定义资源路径、频率点数、别名和似然参数。各目录 README 列出输入文件映射。本索引不改动数据内容或来源信息。
