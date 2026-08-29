# stiffGWpy_fast 审计与推荐总结（2026-08-30）

本页只保留对“能否在 posterior 不显著偏移的前提下把真实 Cobaya MCMC 跑快”最重要的数据。
详细证据见 `audit_phase1.md`（物理等价性）、`audit_phase2.md`（收敛）、`audit_phase3.md`
（918 点参数扫描）、`audit_modes.md`（三档模式）、`audit_mcmc.md`（短链 MCMC 对比）。

## 1. 实际 MCMC speedup（production 模式，短链 N=20，同 seed 同 run-info）

| 指标 | 值 |
|---|---|
| MCMC 墙钟加速 | **585×**（LSODA 526.6 s vs fast 0.90 s） |
| 逐点评估加速 | **1061×**（16.8 s/点 vs 0.016 s/点） |

单点求解（warm，8 线程）：production 4–10 ms；ultra-fast 4–5 ms。LSODA 单点 17–28 s。

## 2. 主要 observable 最大误差（vs LSODA 参考）

| 量 | 全参数空间（918 点 ok=697） | MCMC 主区（DN_eff≤2，n=624） |
|---|---|---|
| ΔN_eff 相对差 | ≤ 5.3e-5（三测试点） | 同左 |
| DN_gw 终值相对差 | **max 2.0e-4**（p95 1.9e-5） | max 2.0e-4 |
| 频谱 dex 误差 | p95 0.077 / max 0.237 | **p95 0.019 / max 0.128** |
| 信号区线性 Ω 相对误差 | p95 0.162 / max 0.421 | p95 0.044 / max 0.255 |

频谱大误差（0.1–0.24 dex）只出现在紧邻 `DN_eff>5` 护栏边界的 ok 点（最差 8 点全部
`DN_eff_lsoda ≥ 2.95`）；终值误差（ΔN_eff / DN_gw / κ_r）在任何子区都 ≤ 2e-4。

## 3. 最大 ΔlogL（同一批点，fast − lsoda，production）

| 统计量 | 值 |
|---|---|
| max / max_abs | **+0.090 / 0.090** |
| median / mean | +0.029 / +0.028 |
| frac \|ΔlogL\|>0.5 / >1.0 | 0 / 0 |

## 4. posterior shift

**尚无收敛链证据，不作 bias 结论。** N=20 短链的均值差无统计意义（见 `audit_mcmc.md` §3）；
需要 `--samples ≥ 2000` 的收敛链（LSODA 链本机约 2–15 小时）才能给出可信的
posterior mean/std/CI/MAP 与 shift。

## 5. failure rate

- 物理区（LSODA `DN_eff ≤ 5`，n=697）：**fast 失败率 0.00%**（697/697 ok）。
- 全参数空间 24.1% 的 `fast_failed` 是**双引擎共享** `DN_eff>5` 护栏拦截
  （生产 `engine='lsoda'` 同样返回 None，已实测），扫描参考循环无护栏才返回发散非物理值；
  这些点被任一引擎拒绝（loglike → −∞），不会污染后验。

## 6. 推荐参数（本机标定）

| 模式 | h | col_step | z_tail | freq_res | 外层 tol | 用途 |
|---|---|---|---|---|---|---|
| `ultra-fast` | 0.01 | 8 | 5 | 1.0 | 1e-6 | 快速探索 / 早期扫描 |
| `production`（推荐 MCMC） | 0.01 | 4 | 7 | 1.0 | 1e-7 | 科研运行默认 |
| `reference` | 0.00125 | 1 | 10 | 2.0 | 1e-8 | 最接近 LSODA；慢，不建议 MCMC |

线程：8 线程后无增益（1→8 线程 16.4→4.8 ms；16 线程 5.3 ms）。MCMC 用法：
`accuracy_mode: production` + `fast_threads: 8` + `fallback: true`（护栏触发自动回退 LSODA）。
内存受限设 `SGWB_POOL_SIZE=2`。

内存优化（2026-08-30，已实施）：`global_param` 移除 `astropy.cosmology` 顶层导入
（`TCMB` 硬编码为 Planck18 精确值 2.7255 K）；A/B 实测每个子进程 `import global_param`
RSS 增量 −24.7 MB（95.8→71.1 MB），数值逐位不变、测试全绿。

## 7. 验证等级（诚实声明）

- **未声称**“完全等价 / 无损 / 全参数空间验证通过”。
- 参数空间扫描完成 **918/1030（91%）**：edge 余 32 点 + extreme 80 点未跑；extreme 角点
  失败率与最坏误差未知（阶段四输入）。
- MCMC 对比为 **N=20 短链**：逐点 ΔlogL 与 failure rate 可信，posterior shift 不可信。
- 工程化（2026-08-30）：`pip wheel` 构建 + 全新目录安装 + 导入/快速求解冒烟通过
  （wheel 33 文件 4.0 MB，无已删除的陈旧数据）；CI 门禁（pytest + ruff + wheel build）
  已加入 `.github/workflows/ci.yml`。

### 续跑命令（补齐后刷新本页）

```bash
python scripts/param_sweep.py --out docs/paramsweep --workers 1 --retry-failed --no-warmup
python scripts/plot_param_sweep.py
python scripts/mcmc_compare.py --samples 2000 --fast-mode production --n-eval 50
```
