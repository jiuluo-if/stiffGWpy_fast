# Reference 下溢积分的 PCHIP 警告对照

日期：2026-11-15

## 样本和方法

使用 `cr=0,n_t=17.38`、冻结 `DN_eff=0.027412228089632874` 的 235 个独立 Reference 频点结果，覆盖 `r={1e-316,1e-317,1e-318,1e-320,5e-324}`。先按 Reference 原路径将对数积分核局部归一化，再调用原始 `integrate_spectrum`，记录完整 warning；随后从同一批有限 log 节点构造 `LogPchip`，用稳定对数算术积分其正值线性 PCHIP。Reference 生产代码未修改，原始调用仍是独立数值基准。

## 结果

- 五个振幅的归一化线性数组都含 13 个精确零和 2 个非零次正规节点；每次原始积分重现 4 条 `overflow encountered in divide`。Stable log-node PCHIP 不产生 warning。
- 稳定积分映射相对原始 PCHIP 映射变化范围 `−7.11e-15` 至 `+7.11e-15`。原始重算结果与保存的 Reference 积分在四个样本完全一致，最小正 `r` 样本相对差 `2.22e-16`。quadrature estimate 约 `1.11e-14`，half-density interpolation estimate 约 `1.5993e-5`。
- 原警告由构造线性 PCHIP 的倒数割线斜率算术触发；零节点仍在被积函数归一化时由 `exp(log_y−shift)` 下溢产生。稳定实现从有限 log 值还原正插值数据，避免除以极小割线，但也会在普通线性数组只能表示为零的节点恢复正值。因此此结果仅说明这些样本的积分几乎不变，不能单独证明新增节点谱偏差或所有参数下的插值误差。

## 判定

稳定对数算术能在这个 235 点下溢样本族消除 warning，并保持 PCHIP 映射到舍入误差；这不是简单屏蔽 warning。仍不修改 Reference 基准实现或默认切换插值器：需要在更多背景/网格上用直接新频率节点验证正值重建偏差，并与已记录的 log-PCHIP 大谱偏差区分。相关同曲线算术验证见[稳定线性 PCHIP 网格诊断](underflow_frequency_and_tail_screen_20261115.md#可复查文件与判定)及其 `235→469/469→937` 新节点结果。

脚本：[underflow_reference_pchip_warning_diagnostic_20261115.py](../../scripts/underflow_reference_pchip_warning_diagnostic_20261115.py)；逐点结果和源码哈希：[JSON](results/underflow_reference_pchip_warning_diagnostic_20261115.json)。此次同时修复并重算了振幅阶梯 JSON 的产物结构；之前变量覆盖导致完整 Reference 原始谱未写入文件，已由当前脚本重现并核验 6 个振幅行和来源 SHA-256。
