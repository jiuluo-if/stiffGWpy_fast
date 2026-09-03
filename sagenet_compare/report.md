# stiffGWpy(收口版) vs SageNet+ — SGWB 谱模拟精度 / 误差 / 执行时间对比

Status: current
Date: 2026-09-03
Code version: stiffgwpy `git_commit = c79ee7c`；SageNet 源码来自 `github.com/ML4GW/SageNet`（`best_gw_model_Transformer.pth` 权重）。

> 诚实声明：本节是**独立可控复现**的实测。所有 stiffgwpy（fast production / plain-grid）与
> SageNet（Transformer 神经网络）的谱，都在同一个参数点上算，并与**独立连续-sigma reference**
> 真值比较。数值均读取自 `sagenet_compare/data/report.json` / `quick.json`（本仓库内，可由
> `sagenet_compare/run_compare.py` 重放）。SageNet 侧是"该权重 + 该 venv"的直接推理结果，不是论文声称值。

---

## 1. 要回答的问题

对一个宇观学参数点，四种方法输出 SGWB 能量密度谱 `log10 Omega_GW(f)`：

| 方法 | 类型 | 说明 |
|---|---|---|
| `reference` | stiffgwpy 独立连续-σ DOP853 | **真值锚点**（signal 区频率子集，rtol=1e-6，z_tail=8） |
| `production` | stiffgwpy fast（transition-refine） | 生产档，谱内收敛 ~7e-4 |
| `plain_grid` | stiffgwpy fast（plain-grid） | 速度优先，已知系统偏置 |
| `sagenet` | SageNet+ Transformer NN | 深度学习模拟器 |

**精确率**这里不是分类准确率（无分类任务），而是**回归忠实度**：谱点 `|Δlog10 Ω_GW|`（dex）、
Omega 空间相对误差 `|ΔΩ|/Ω`（rel），以及 SageNet 作者口径的 `calculate_area_difference`
（对数空间面积相对差）。

## 2. 计算环境

* stiffGWpy 在独立 venv `F:\codex\sagenet_env` 跑：`numpy 1.26.4`，`scipy 1.14.1`，
  `numba 0.60.0`，`astropy 6.1.0`，`torch 2.3.1+cpu`，`scikit-learn 1.6.1`。
* torch 在本机 base 环境因 `c10.dll` DLL 初始化失败（WinError 1114，新版 torch 的
  c10 DllMain 不兼容）无法直接用；venv 内通过 `KMP_DUPLICATE_LIB_OK` +
  `torch\lib`+`C:\miniconda3\Library\bin` 加入 DLL 搜索路径解决。SageNet 权重由 venv 加载。
* 全频 `reference` 每点约 `6+ min`（低频亚视界模式需积分到今天）；本报告真值用
  **signal 区频率子集**（`log10 f ∈ [-6, 1]`，均为已重入、可快速触发尾部事件的模式），
  单点 23–55 s，覆盖驱动可观测量的物理区。

## 3. 参数点（均在 SageNet+ 合法盒内）

| 标签 | r | n_t | kappa10 | T_re [GeV] | DN_re | Ω_bh2 | Ω_ch2 | H0 | A_s | 备注 |
|---|---|---|---|---|---|---|---|---|---|---|
| sage_center | 3.96e-5 | 1.011 | 110.4 | 0.175 | 39.37 | 0.0224 | 0.1201 | 67.32 | 2.10e-9 | SageNet README 示例点 |
| neutral | 1e-3 | 0.0 | 1.0 | 1e3 | 20 | 0.0224 | 0.1201 | 67.32 | 2.10e-9 | 中等刚度 |
| red_tilt | 1e-1 | -0.5 | 5.0 | 5e2 | 10 | 0.0224 | 0.1201 | 67.32 | 2.10e-9 | 大 r、红移、弱 stiff |
| stiff | 5e-3 | 0.0 | 50.0 | 1e2 | 30 | 0.0224 | 0.1201 | 67.32 | 2.10e-9 | 高 stiff |
| high_Tre | 1e-2 | 0.2 | 0.05 | 5e4 | 5 | 0.0224 | 0.1201 | 67.32 | 2.10e-9 | 高温重加热；fast 触发共享 N_eff 守卫 |

## 4. 谱误差（相对真值，signal 区 `log10 f ∈ [-6,1]`）

| 点 | 方法 | dex_max | rel_max (Ω) | area_diff（SageNet 口径） |
|---|---|---|---|---|
| sage_center | production | **3.3e-4** | **7.5e-4** | 6.7e-4 |
| sage_center | plain-grid | 2.17e-2 | 5.1e-2 | 3.96e-2 |
| sage_center | SageNet | 2.26e-1 | 4.1e-1 | 1.87e-1 |
| neutral | production | **3.4e-4** | **7.8e-4** | 4.6e-4 |
| neutral | plain-grid | 7.5e-3 | 1.7e-2 | 7.1e-3 |
| neutral | SageNet | 1.30e-1 | 3.5e-1 | 1.20e-1 |
| red_tilt | production | **2.6e-4** | **6.0e-4** | 3.5e-4 |
| red_tilt | plain-grid | 7.6e-3 | 1.7e-2 | 7.1e-3 |
| red_tilt | SageNet | 2.74e-1 | 8.8e-1 | 2.52e-1 |
| stiff | production | **3.3e-4** | **7.5e-4** | 2.4e-4 |
| stiff | plain-grid | 7.5e-3 | 1.7e-2 | 7.1e-3 |
| stiff | SageNet | 1.22e-1 | 3.2e-1 | 1.36e-1 |
| high_Tre | SageNet（无 fast，fast 被物理守卫拒绝） | 8.0e-1 | 8.4e-1 | 9.2e-1 |

## 5. 执行时间（本机 warm，单点）

| 点 | production [s] | plain-grid [s] | SageNet [s] | reference(signal) [s] |
|---|---|---|---|---|
| sage_center | 3.12 | 0.33 | 0.033 | 22.9 |
| neutral | 2.00 | 0.30 | 0.021 | 35.0 |
| red_tilt | 1.60 | 0.12 | 0.020 | 34.7 |
| stiff | 2.27 | 0.34 | 0.020 | 34.9 |
| high_Tre | (守卫拒绝) | (守卫拒绝) | 0.020 | 55.0 |

速度对比（signal 子集真值）：

* **SageNet vs production** ≈ **60–150×**（0.02 s vs 1.6–3.1 s）。
* **SageNet vs plain-grid** ≈ **6–17×**（0.02 s vs 0.12–0.34 s）。
* **SageNet vs reference** ≈ **1000–1500×**（0.02 s vs 23–55 s）；对**全频 reference (~360 s)** 达
  **~1.5e4×**，符合 SageNet 声称的 ~10000× 量级。
* **production vs reference** ≈ **12–22×**（signal 子集）；对全频 reference ≈ **150×**。

## 6. 关键结论与本轮答案

1. **stiffgwpy 两个 fast 档的"精密 vs 高速"边界是真实的**：
   `production` 谱误差稳健在 **6–8e-4（rel）/ ~3e-4（dex）**，与上一轮 matched-z8 的一致；
   `plain-grid` 落到 **1.7e-2（rel）/ 7.5e-3（dex）**（高 stiff 点 5.1e-2）——快 ~10× 但精度差 ~20×，
   因此 plain-grid 仅适合作探索/筛查外层，不适合科学结论。这与 stiffGWpy 本仓库的认证一致。
2. **SageNet（Transformer）在本机实测的 signal 区谱误差显著大于 stiffgwpy 生产档**：
   `dex_max 0.12–0.80`、`rel_max 30–88%`、`area_diff 0.12–0.92`。
   主峰附近两者接近（~0.1 dex），**误差集中于低频 / 过渡段**（stiffgwpy 峰 −7.5 vs SageNet −7.4，
   但低频端 stiffgwpy −13.7 vs SageNet −13.1，差 ~0.6 dex）。
3. **"纯速度"的代价**：SageNet 以 ≈0.02 s/点 换取 ~1.5e4× 速度，但在本机本权重下，
   **其精度在信号区远未达到 stiffgwpy production 的 1e-3 级**，也低于 SageNet 论文声称的 ~1%。
4. **为何 SageNet 误差偏大（诚实 caveat）**：
   * SageNet 权重/训练数据**不是本仓库连续-σ reference 生成的**；若训练基于旧版 stiffGWpy
     （LSODA 路径），两者在低频/过渡段可能有系统性差异。
   * 这些参数点可能落在 SageNet 训练分布边缘或低频分辨率不足（256 点固定网格）。
   * 我无法验证 5400 样本在 log10 f ∈ [-6,1] 的覆盖与模型版本，因此**这是"该权重+该点"的直接
     实测**，不应外推为"SageNet 全局弱"。峰的绝对定位（~0.1 dex）证明其整体缩放正确。
5. **参考真值的使用范围**：本报告真值是 signal 区频率子集（`log10 f ∈ [-6,1]`），
   `DN_gw` 为子集积分；全频 `reference` 因低频模式积分成本过高（~360 s/点）未使用。

## 7. 复现

```bash
# 1) 建独立 venv 并装依赖（numpy<2 以兼容 torch 2.3.1；astropy 6 兼容 numpy 1.26）
F:\codex\sagenet_env\Scripts\python.exe -m pip install "numpy==1.26.4" "scipy==1.14.1" \
    "astropy==6.1.0" "numba==0.60.0" pyyaml scikit-learn==1.6.1
F:\codex\sagenet_env\Scripts\python.exe -m pip install "torch==2.3.1+cpu" \
    --index-url https://download.pytorch.org/whl/cpu

# 2) 在仓库根，把本地 stiffgwpy 接到 SageNet 子模块目录，并克隆 SageNet 源码
#    (sagenetgw/stiffGWpy -> F:\codex\stiffGWpy\stiffgwpy；SAGE_SRC 指向克隆目录)

# 3) 运行
F:\codex\sagenet_env\Scripts\python.exe sagenet_compare\run_compare.py --phase quick
F:\codex\sagenet_env\Scripts\python.exe sagenet_compare\run_compare.py --phase reference --points neutral stiff red_tilt sage_center high_Tre
F:\codex\sagenet_env\Scripts\python.exe sagenet_compare\run_compare.py --phase report
```

数据：`sagenet_compare/data/quick.json`（各方法谱 + 时间）、`sagenet_compare/data/report.json`
（误差汇总）。脚本：`sagenet_compare/run_compare.py`。

## 8. 主要限制

* `reference` 仅在 signal 区（`log10 f ∈ [-6,1]`）作为真值；低频超视界尾部（静态 Ω_GW 定义不清、
  无 ΔN_eff 权重）被排除。
* `high_Tre` 点 stiffgwpy fast 两种档都被共享 `ΔN_eff > 5` 物理守卫拒绝（不是求解器失败），
  故该点只有 SageNet-vs-reference，未给 fast 误差。
* SageNet 侧在本机 venv 用单一 `Transformer` 模型；未测 LSTM/RNN/CosmicNet2。sklearn scaler
  版本已与 checkpoint 对齐（1.6.1）。
* 这是"该权重"的独立复现，**未**复现 SageNet 论文训练流程；论文声称的 ~1% 需要其原始训练数据/流程。
