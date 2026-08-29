# 三档精度模式审计（2026-08-30）

审计对象：`stiffgwpy.fast_sgwb` 的命名精度模式（`ACCURACY_MODES`：`reference` / `production` /
`ultra-fast`）vs `stiffgwpy.stiff_SGWB`（LSODA 参考）。
工具：`scripts/validate_modes.py`（三模式 × 三点 default/p1/p2 的 LSODA-vs-fast 验证 + 线程缩放，
支持断点续跑）；原始记录 `docs/modes/validate_modes.jsonl`、`docs/modes/thread_scaling.json`。

## 0. 结论先行

1. **production / ultra-fast 在三个测试点上的终值精度满足预期。**
   ΔN_eff 相对差 ≤ 5.3e-5（production ≤ 5.3e-5，ultra-fast ≤ 1.2e-5），DN_gw[-1] 与 κ_r 同量级；
   频谱 dex_max 为 3e-4–1.7e-1（见第 3 点口径说明），中位 dex 误差 3e-6–1.5e-4。
2. **warm 单点求解 4–10 ms，端到端相对 LSODA 加速 ~4e3–2e4×**（含冷启动的首次 JIT 除外，
   见 5.2 节）；多线程收益在 8 线程后饱和（8→16 线程不再变快）。
3. **p2 点的 dex_max / lin_max 大值是单通道极值，不能当典型误差**：p2
   （cr=0, r=1e-2, n_t=0.2, DN_re=5, T_re=5e3, kappa10=1e-3）的频谱在近零/边缘通道出现
   ~0.17 dex 的通道级偏差，而中位 dex 误差仅 1e-4–1.5e-4、ΔN_eff 相对差 2e-6–1e-5；
   阶段二扫描（`docs/audit_phase2.md`）与阶段三唯一完整点（`docs/audit_phase3.md`）观测到
   同样的现象，属于近零区的相对误差放大，与 fast 引擎无关（LSODA 自身在该区同样敏感）。
4. **reference 模式未在本轮跑 LSODA 对照**（LSODA h=0.00125 单点 >30 分钟，被终止）。
   其引擎差引用阶段二同网格数据：h=0.00125 时 fast-vs-LSODA 引擎差降至 2.6e-7（ΔN_eff），
   但该模式每点 ~0.1–0.4 s（8 线程），比 production 慢一个量级，**不建议用于 MCMC**。

## 1. 三档模式定义

| 模式 | h | col_step | z_tail | freq_res | 外层 tol | 目标用途 |
|---|---|---|---|---|---|---|
| `reference` | 0.00125 | 1 | 10 | 2.0 | 1e-8 | 最接近 LSODA 参考（最慢） |
| `production` | 0.01 | 4 | 7 | 1.0 | 1e-7 | 科研运行默认 |
| `ultra-fast` | 0.01 | 8 | 5 | 1.0 | 1e-6 | 快速探索性扫描 |

取值依据见 `docs/audit_phase2.md`：h=0.00125 时引擎差 2.6e-7；z_tail=7 将解析尾部误差降至
~2e-5（z_tail=10 降至 2.4e-7）；col_step 对 ΔN_eff 无影响（<1e-9）；freq_res=2.0 减半低频
陡尾欠采样误差。

## 2. 测试点

- `default`：r=1e-2, cr=1, T_re=2e3, kappa10=1e-2（12-case 网格 case A）；
- `p1`：r=5e-3, cr=1, T_re=1e2, kappa10=0.1（阶段二随机验证起点 P1）；
- `p2`：r=1e-2, cr=0, n_t=0.2, DN_re=5, T_re=5e3, kappa10=1e-3（stiff 主导 / 无一致性关系路径）。

## 3. 本文验证结果（LSODA rtol=1e-8, atol=[1e-12,1e-22,1e-22], 外层 tol=1e-7）

| 模式 | 点 | ΔN_eff 相对差 | dex_max | dex 中位 | lin_max | lin 中位 | t_LSODA (s) | t_fast (s) | 加速 |
|---|---|---|---|---|---|---|---|---|---|
| production | default | 5.3e-6 | 3.9e-4 | 9.1e-5 | 9.0e-4 | 2.1e-4 | 127.6 | 0.0063 | ~2.0e4× |
| production | p1 | 5.3e-5 | 2.9e-4 | 8.5e-5 | 6.7e-4 | 2.0e-4 | 132.7 | 0.0101 | ~1.3e4× |
| production | p2 | 2.3e-6 | 1.7e-1* | 1.5e-4 | 3.2e-1* | 3.4e-4 | 61.7 | 0.0063 | ~9.8e3× |
| ultra-fast | default | 1.1e-5 | 5.7e-3 | 3.8e-6 | 1.3e-2 | 8.8e-6 | 28.6 | 0.224** | ~128× |
| ultra-fast | p1 | 1.2e-5 | 3.0e-3 | 3.0e-6 | 6.8e-3 | 7.0e-6 | 27.8 | 0.0050 | ~5.5e3× |
| ultra-fast | p2 | 1.0e-5 | 1.7e-1* | 8.3e-6 | 3.2e-1* | 1.9e-5 | 16.8 | 0.0041 | ~4.1e3× |
| reference | default/p1/p2 | fast_only（引用阶段二 h=0.00125 引擎差 2.6e-7） | | | | | | 0.14–0.40 | — |

*：单通道极值，位于频谱近零/边缘区（与阶段二/阶段三现象一致），典型误差看中位列。
**：ultra-fast default 首次含 JIT 编译（冷启动 0.224 s）；warm 各点 4–10 ms，见 5.2 节。

## 4. 线程缩放（production / default，5 次取中位）

| 线程 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| t_fast (ms) | 16.35 | 10.16 | 6.79 | 4.80 | 5.25 |

8 线程后不再增益（16 线程略回升）。默认 `threads=8` 是该机型的合理上限；
`reference` 模式同理，不必盲目加大线程数。

## 5. 复现

```bash
python scripts/validate_modes.py --out docs/modes/validate_modes.jsonl   # 三模式 × 三点
python scripts/validate_modes.py --include-reference-lsoda               # 可选：reference 的昂贵 LSODA 对照
python scripts/validate_modes.py --thread-scaling --out docs/modes/thread_scaling.json
```

注意：LSODA 参考路径按 `SGWB_POOL_SIZE` 开进程池（默认 4；MPI 下默认 1）。本机内存脆弱，
大批量运行时建议 `SGWB_POOL_SIZE=2`。

## 6. NOT CERTIFIED 清单

- **Cobaya MCMC 后验对比未完成**（LSODA chain vs fast chain、ΔlogL、后验偏移）；adapter
  （engine/fallback/threads/h/col_step/z_tail/freq_res/accuracy_mode）已实现并通过单测，
  但需在装有 cobaya 的环境做真实采样对比。
- **1000 点参数空间扫描未完成**（`docs/audit_phase3.md`：本机提交内存墙）；本轮仅补了
  三模式 × 三点的 LSODA 对照。
- 因此三档模式是**经验推荐配置**，不是全面认证；跨参数空间的最坏情况误差仍以
  阶段三全量扫描为准（未完成）。

## 7. 使用建议

- 科研运行 / MCMC：`production`（若链太长可退 `ultra-fast`，误差代价约 +6e-6 相对）；
- `reference` 仅用于局部交叉验证，不做 MCMC（单点 ~0.1–0.4 s，比 production 慢 ~50–100×，
  且 LSODA 对照更贵）；
- 在 Cobaya yaml 中通过 `accuracy_mode: production` 一键选择；显式 `h/col_step/z_tail/threads`
  仍可覆盖预设。
