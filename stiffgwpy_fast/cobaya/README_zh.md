# Cobaya 适配器程序包

[English](README.md) | 中文说明

本程序包提供 Cobaya theory 适配器及随包安装的似然模块。选项和调用契约详见[`docs/cobaya.md`](../../docs/cobaya.md)。

## Theory 入口

在 Cobaya YAML 中使用完整类名：

```yaml
theory:
  stiffgwpy_fast.cobaya.stiffGW.stiffGW:
    engine: fast
    fallback: True
    accuracy_mode: fast
```

串行适配器只需安装 `cobaya` extra。只有运行环境确实使用 `mpi4py` 时才安装单独的 `mpi` extra。

## 数值契约

- 适配器按“`accuracy_mode` → 预设默认值 → 用户显式覆盖”的顺序解析配置。YAML 中的零值哨兵不会屏蔽所选预设。
- `fast` 是唯一正式的用户档位。只有启用回退后，数值故障才会尝试 LSODA；确定性的物理保护不会重试。
- `eval_freqs` 会把似然频率作为 fast 求解器的原生节点；这些节点不进入玻尔积分支撑网格。
- `engine_stats` 提供本次运行的引擎、故障、保护、回退和升级遥测。

本目录还包含 LVK、PTA 似然适配器及其数据资源。适配器与数据文件的对应关系见[似然与资源索引](likelihoods/README_zh.md)。发布 wheel 前请按复现和打包说明检查资源是否完整。
