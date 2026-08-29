# Cobaya MCMC 对比审计（短链，N=20）

**日期：** 2026-08-30
**运行环境：** Windows 本机，cobaya 3.6.2，Python 3.11（miniconda3），`SGWB_POOL_SIZE=2`，
`fast_threads=8`，accuracy mode = `production`。
**run-info：** `stiffgwpy/cobaya/mcmc_compare.yaml`（NANOGrav + LVK_SGWB_CC，仓库内数据）。
**命令：**
`python scripts/mcmc_compare.py --samples 20 --n-eval 10 --fast-mode production --fast-threads 8 --seed 20260830`

## 1. 真实 MCMC 对比结果

两条链同一 seed、同一 run-info，仅 `engine` 不同（lsoda / fast），各 20 个接受样本。

| 指标 | 值 |
|---|---|
| LSODA 链墙钟 | 526.6 s（约 26 s/样本，含低接受率） |
| fast 链墙钟 | 0.90 s（0.047 s/样本） |
| **MCMC 墙钟加速** | **585×** |
| 逐点评估（10 点）LSODA | 167.8 s（16.8 s/点） |
| 逐点评估（10 点）fast | 0.158 s（0.016 s/点） |
| **逐点评估加速** | **1061×** |
| LSODA 链 failure rate | 0/20 |
| fast 链 failure rate | 0/20 |
| 逐点对照 n_both_ok | 10/10（两引擎均有限） |

## 2. 逐点 ΔlogL（同一批点，fast − lsoda，logL = −minuslogpost 方向）

| 统计量 | 值 |
|---|---|
| n | 10 |
| median | +0.029 |
| mean | +0.028 |
| max / max_abs | +0.090 / 0.090 |
| p95 绝对值 | 0.087 |
| frac \|ΔlogL\| > 0.5 | 0 |
| frac \|ΔlogL\| > 1.0 | 0 |

逐点一致性很强：两个引擎在完全相同参数点上的 log-likelihood 差异最大约 **0.09**
（远小于 1），这是“引擎级差异”最可靠的度量。

## 3. 后验偏移（仅供参考，不做结论）

`(mean_fast − mean_lsoda)/std_lsoda`（15 个参数，含派生量）：

| 参数 | shift | 参数 | shift |
|---|---|---|---|
| log10r | −5.27 | A_BBH | −1.54 |
| n_t | +6.42 | gamma_BBH | −0.41 |
| log10T_re | −15.41 | h | 0.00（常数） |
| DN_re | −0.28 | r | −3.53 |
| log10kappa10 | +1.80 | T_re | −14.44 |
| kappa10 | +2.36 | Delta_Neff_GW | −2.48 |
| Delta_Neff_total | −2.48 | log10hc_prim_fyr | −6.24 |
| f_end | +16.80 | | |

**该表不得用于判定 fast 造成 posterior bias。** 理由（如实说明）：

- N=20 远未收敛（7 维参数空间）；链内 std 反映的是局部游走步长而非后验宽度，
  “shift/std”在未收敛链上没有统计意义。
- 同 seed 下两条链在首次被拒后即分叉，各自停留在不同邻域；均值差主要来自
  采样噪声与未混合，而不是引擎偏差。
- 能证明引擎一致性的指标是 §2 的逐点 ΔlogL（max 0.09）以及 failure rate=0，
  而不是未收敛链的均值差。

## 4. 推荐参数（本机）

- **MCMC 主运行：`accuracy_mode: production`，`fast_threads: 8`**（8 线程后无增益，
  见 `docs/audit_modes.md` 线程标定）。预计 2000 接受样本约 1–2 分钟（fast 链），
  逐点 LSODA 对照按 `--n-eval 50` 约 14 分钟。
- 内存受限时设 `SGWB_POOL_SIZE=2`；MPI 环境该变量默认 1。

## 5. 结论与 NOT CERTIFIED

- **已验证：** 端到端 Cobaya MCMC 链路（`run` → `sampler.products()` → 逐点
  `model.logpost` 对照）在本机跑通；production 模式下 MCMC 墙钟加速 ~585×，
  逐点评估加速 ~1061×，两引擎 failure rate 均为 0，逐点 max|ΔlogL| = 0.09。
- **NOT CERTIFIED：**
  1. “fast 不造成 posterior bias”尚无收敛链证据——需 `--samples >= 2000` 的
     收敛链（LSODA 链约 2–15 小时，本机）才能给出可信的 posterior mean/std/CI/MAP
     偏移结论；
  2. 1000 点全参数空间扫描仍未完成（`docs/audit_phase3.md`），worst-case 误差
     与失败区域仍以阶段三全量扫描为准；
  3. 三档模式（ultra-fast/production/reference）仍是经验推荐，不是全面认证。
