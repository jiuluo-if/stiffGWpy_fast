# 第三阶段审计：大规模参数空间验证（2026-08-30）

审计对象：`stiffgwpy.fast_sgwb`（fast）vs `stiffgwpy.stiff_SGWB`（LSODA 参考）。
目标：建立 ≥1000 个确定性参数点的 LSODA-vs-fast 交叉验证，输出 worst cases / 误差分布 / 参数误差图。
工具：`scripts/param_sweep.py`、`scripts/plot_param_sweep.py`（均已通过冒烟验证与 ruff）。

**结论：第三阶段 PARTIALLY CERTIFIED（91% 完成；用户要求停止长跑，非内存墙）。**
扫描已中止于 **918/1030 点（91%）**：sobol 400 + lhs 350 + edge 168 全部落盘
（`docs/paramsweep/sweep_phase3.jsonl`），extreme 80 点未启动。中止原因是**用户要求停止
长跑**，不是内存墙（本轮以 `SGWB_POOL_SIZE=2` 约束池大小，全程无提交内存墙复现）。
已完成 918 点给出误差分布、worst cases、参数误差图与失败区画像；严格 ≥1000 点认证
（含 extreme 角点）仍未完成。项目整体维持 **PARTIALLY CERTIFIED**（阶段一物理等价性、
阶段二 default 点收敛、阶段三 91% 扫描证据）。

## 1. 设计与协议

- 采样（确定性，seed=20260830，共 1030 点）：Sobol 400 + LHS 350 + 边界过采样 LHS 200
  （unit 坐标向 0/1 推 `|u-0.5|^0.6`）+ 手选极端/角点 80（单参数上下界 22、联合极端 22、
  cr=1 的 (r,T_re) 角格 24、cr=0 的 (n_t,DN_re) 角格 12）。
- 参数覆盖：Omega_bh2 [0.018,0.026]、Omega_ch2 [0.09,0.15]、H0 [60,76]、DN_eff [0,2]、
  A_s [1e-9,4e-9]（log）、r [1e-4,1e-1]（log）、n_t [-0.5,0.5]、cr {0,1}、T_re [10,1e6]（log）、
  DN_re [0,30]、kappa10 [1e-4,1]（log）。对数参数 log 采样，线性参数 linear 采样。
- 每点：同一 σ 网格上 LSODA（rtol=1e-8, atol=[1e-12,1e-22,1e-22], outer tol=1e-7）vs fast
  （h=0.01, col_step=4, z_tail=5, freq_res=1.0, outer tol=1e-7）；记录 status、耗时、ΔN_eff_final、
  DN_gw[-1]、κ_r、频谱 dex 误差（max/p50/p95/p99）、信号区线性 Ω 相对误差（max/median）、
  DN_gw 曲线相对误差（max/median）。
- 绘图脚本产出：`summary.json`、`worst_cases.json`（worst 20 × 2 维度）、`error_distribution.png`、
  `parameter_error_map/`（11 张，含 Spearman 相关）、`points.json`、`sweep_phase3.jsonl`（918 条）。

## 2. 执行与中止：918/1030（91%，用户要求停止）

- 执行：`--workers 2` + `SGWB_POOL_SIZE=2`，checkpoint 续跑；日志停在 `[895/1010]`，
  jsonl 实落盘 918 条 = sobol 400 + lhs 350 + edge 168（edge 余 32 点、extreme 80 点未跑）。
- 中止原因：**用户要求停止长跑**（LSODA 单点中位 28 s，全量 1030 点约 8–10 小时级），
  非内存墙；本轮未复现早期提交内存墙。

## 3. 已完成 918 点的结果

### 3.1 完成率与失败率

- ok=697 / fast_failed=221，失败率 **24.1%**（sobol 85、lhs 83、edge 53）。
- 失败区特征（failed vs ok 参数中位数）：r 0.0083 vs 0.0022、kappa10 0.056 vs 0.0054、
  T_re 1.6e5 vs 762、n_t 0.23 vs -0.08、cr=0（ok 组中位 cr=1）。
  → 失败集中在 **cr=0 且高 T_re、高 r/kappa10、较高 n_t** 的区域（fast 在这些点数值发散），
  与物理配置相关而非随机。

### 3.2 ok 点（n=697）误差分布

| 指标 | median | p95 | max |
|---|---|---|---|
| DN_gw_last_rel | 9.9e-6 | 1.9e-5 | **2.0e-4** |
| dex_max | 3.0e-3 | 7.7e-2 | 0.237 |
| lin_max | 6.8e-3 | 0.162 | 0.421 |

- DN_gw 终值相对差全样本 max 2.0e-4，在 `validate_fast.py` 的 1e-3 门禁内；
  但频谱 dex / 线性 Ω 尾部误差显著（p95 达 0.08/0.16），集中出现在失败区边界的 ok 点。

### 3.3 最差 ok 点（worst_cases.json）

- 按 DN_gw_last_rel：`edge-0109`（r=0.0187, n_t=0.5, cr=0, T_re=132, kappa10=1e-4）：
  DN_gw_last_rel=2.0e-4、dex_max=0.070、lin_max=0.150。
- 按 dex_max：`lhs-0300`（r=0.0082, n_t=0.14, cr=0, T_re=6e3, kappa10=0.0096）：
  dex_max=0.237、lin_max=0.421、DN_gw_last_rel=1.1e-5。

### 3.4 参数敏感性（Spearman，log10 DN_gw_last_rel vs 参数）

- 主导：n_t +0.38、cr +0.26、T_re −0.19、kappa10 −0.15；其余 |ρ|<0.05（r≈0、
  DN_eff≈0.01、Omega_bh2≈0.01 等）。即误差主要随谱指数/重加热/耗散参数变化，与
  失败区画像一致。

## 4. 未完成清单与恢复路径

- 未完成：edge 余 32 点 + extreme 80 点（共 112 点，9%）；严格 ≥1000 点认证、extreme
  角点失败率与最坏误差（阶段四输入）。
- 恢复：`python scripts/param_sweep.py --out docs/paramsweep --workers 1 --retry-failed --no-warmup`
  （checkpoint 续跑，`--retry-failed` 重跑失败点），完成后 `python scripts/plot_param_sweep.py`
  刷新汇总与图。
- 低内存主机建议：将 `stiff_SGWB.run_SGWB` 的硬编码 `mp.Pool(4)` 改为可配置池大小
  （env var，默认 4；纯进程管理改动，不改变任何数值结果），或延迟化 `global_param` 的
  astropy 顶层 import（每子进程私有提交可降至 ~0.8 GB）。
