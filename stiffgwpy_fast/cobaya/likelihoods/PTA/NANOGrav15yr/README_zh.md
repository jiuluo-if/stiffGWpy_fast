# NANOGrav 15 年 mock 数据资源

[English](README.md) | 中文说明

`../NANOGrav.yaml` 使用本目录中的三个 NumPy 资源：

- `density_mock.npy`：逐频点似然密度样本。
- `log10rhogrid.npy`：似然密度对应的 `log10(rho)` 网格。
- `freqs.npy`：单位为 Hz 的 PTA 频率 bin。

适配器按配置路径加载数组，并用 `Nfreqs` 选择参与计算的频点。请保持文件名及其与 YAML 的对应关系。本说明不新增数据来源结论。
