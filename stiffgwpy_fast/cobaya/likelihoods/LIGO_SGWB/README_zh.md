# LVK 交叉相关似然

[English](README.md) | 中文说明

`LVK_SGWB_CC.py` 提供 Cobaya 似然 `LVK_SGWB_CC`。默认选项见 `LVK_SGWB_CC.yaml`；其中 `CC_file` 指向 `C_O1_O2_O3.dat`。

适配器读取三列数据：频率、实测交叉相关量和不确定度。它将频率转换为 `log10(Hz)`，在理论谱支持范围内插值，并计算高斯卡方 log-likelihood。似然别名定义在 YAML 中。

更多入口见[似然包索引](../README_zh.md)、[Cobaya theory 说明](../../../../docs/cobaya.md)和[根目录 README](../../../../README_zh.md)。本页只说明代码和文件布局，不单独作出数据来源或模型拟合结论。
