# 脉冲星计时阵列（PTA）似然

[English](README.md) | 中文说明

本目录包含 IPTA 和 NANOGrav 似然适配器。类选项和默认资源路径分别定义在 `IPTA.yaml` 与 `NANOGrav.yaml`。

| 适配器 | 默认输入 | 计算方式 |
|---|---|---|
| `IPTA` | `EPTAdr2/EPTA_dr2new_mock.dat`、`EPTAdr2/freqs_dr2new.txt` | 使用配置的频率 bin、理论谱和 SMBHB nuisance 参数，逐频点计算样本密度似然。 |
| `NANOGrav` | `NANOGrav15yr/density_mock.npy`、`NANOGrav15yr/log10rhogrid.npy`、`NANOGrav15yr/freqs.npy` | 将理论谱与 SMBHB 分量合并后，插值随包提供的逐频点似然网格。 |

`Nfreqs` 指定适配器使用的频点数。这些文件作为包内资源读取；应以 YAML 路径为准，不要假设当前工作目录是仓库根目录。另见[似然包索引](../README_zh.md)和[`docs/cobaya.md`](../../../../docs/cobaya.md)。
