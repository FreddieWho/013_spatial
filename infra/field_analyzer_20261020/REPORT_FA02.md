# FA-02 交付：合成 kill-test 生成器与配对诊断

**日期：** 2026-10-20
**节点：** FA-02（计算验证；服务 AF-01/02/03，本身不是生物学证据）
**口径：** 06_SYNTHETIC_KILL_TESTS.md + `docs/FA_CONTRACT_ADDENDUM_20261020.md`

## 1. 两类生成器

- **G1 解析/样条**：高斯峰、环、角向加厚环（可精确匹配径向均值）、双峰、Bezier 管道、Y 形分支、冗余环。几何真值显式。
- **G2 异机制**：网格图上的**各向异性扩散**（conductance_x ≠ conductance_y，带衰减），真值来自传播而不是几何模板。测试验证了它是各向异性的（x 方向二阶矩约为 y 的 1.3 倍），且当 conductance 相等时各向同性（相对差 <5%）。默认场景未使用 G2，作为机制错配检验保留。

两者都不从"待测模型"的公式抽样，避免生成器同构造成虚假优势。

## 2. 场景与配对诊断（实测）

64×64 规则格，pitch = 1（`synthetic_pitch`）。每场 12 个 marker，另 128 个背景基因。计数为负二项（dispersion 8）、对数正态深度（sd 0.25），S0 另加共享技术深度梯度。

| 场景 | 构造 | 配对是否真的匹配（实测） |
|---|---|---|
| S0 | 4 个独立平滑场 + 共享深度梯度 | 场间相关 <0.6，无跨场结构 |
| S1 | 圆环 vs 角向加厚环；单峰 vs 双峰 | 径向剖面 max\|Δ\| = **2.2e-3 / 2.0e-10**；形变幅度 >0.1；双峰质心差 **1.1e-4**、质量差 <2e-14、半高以上峰数 = 2 |
| S2 | 同点相关匹配的伙伴 vs 空间偏移伙伴 | 同点相关 **0.822 vs 0.822**（二分求解混合权重 w=0.855），实际峰位移 **7.0 pitch**，标 `separable=True` |
| S3 | 连续桥 / 断桥（中段 8 pitch 置零）/ 错端 / Y 分支 / 冗余环 + 三个端点区 | 断桥中段最大强度 = 0；错端与左端点峰值位置不同 |
| S4 | 同一真桥 + 0/4/16/32 个无关平滑场 | 真桥数组逐位不变；场数 8→12→24→40 |
| S5 | 可见计数逐位相同，隐藏圆盘（半径 8 pitch）内桥存在/不存在 | 可见计数 max\|Δ\| = **0.0**（两个种子）；隐藏区真值明显不同 |

S2 的设计值得单独记一笔：**"同点相关相同"与"空间峰位移"天然冲突**。混合权重越大峰越偏移、同点相关越低。本实现用二分法把同点相关精确匹配到目标值，并**如实报告**因此只能达到 7 pitch 的位移（请求 8）。若换成更小位移，配对会不可分——这就是后面评价 B 路线时必须先看配对诊断的原因。

## 3. 数据组织（物理分离）

每个样本一个目录：

```
synth/samples/<scenario>_<base_id>/
  observed_inputs/counts_visible.npz   模型可读：计数（隐藏行为 NaN）、gene、xy、visible_mask
  latent_truth/fields.npz              评分用真值（字段图）
  evaluation_only/answers.npz          mask、结构参数、配对答案
  meta.json                            场景、种子、角色、字段数、点/基因数
```

模型侧读取器只拿到第一项；第三项含结构参数，训练器不得读取。

种子角色互斥且完整：train 0–7、development 8–11、evaluation 12–27。本次贯通只用 **0 和 8**（train 10 个样本、development 6 个），**evaluation 种子 12–27 未被消耗**，留给后续确认性比较。同一 `base_id` 的配对变体同角色。

## 4. 产物规模

```
synth/scenario_manifest.tsv          16 行（base_id / scenario / pair_id / seed / role / n_fields / 路径）
synth/pair_matching_diagnostics.tsv  38 行（含 S1/S2/S5 配对残差）
synth/generator_parameters.json      生成参数与说明
synth/summary.json                   计数与关键诊断
synth/samples/**                     16 个样本，共 11 MB
```

生成耗时 7.6 s，峰值内存约 200 MB，纯 CPU。

## 5. 边界与未做

- 合成成功不等于生物真实；本轮**没有任何路线跑过这些数据**（那是 FA-03）。
- S1 的径向匹配残差 2.2e-3 不是零——报告里给的是实测值，不声称精确匹配。
- S3 的"断开"用 8 pitch 宽的置零带实现，真实病理的断点更窄；窄于一个 pitch 的断桥应视为分辨率边界，不能要求算法恢复。
- S5 只检验"应当承认看不出"，不检验预测是否正确（也不应该要求猜对）。
- S4 只准备材料，**尚未验证"场数增加时不增加路径预算"**——那是 C 路线的评价规则，须在 FA-03C 里作为独立断言实现。
- G2 未进入默认 manifest；用于后续"换机制是否仍成立"的复核。

## 6. 对 AF 假设的影响

无直接支持或否定。FA-02 只提供可判定性与反例材料：它使"A 优于 B / C 有独立价值"这类陈述第一次具备可证伪的形式。按 roadmap 规则，生成器本身记为基础设施。

## 7. 复现

```bash
/opt/anaconda3/bin/python -m field_analyzer.cli synth --seeds 0,8
/opt/anaconda3/bin/python -m pytest tests/test_fa_synthetic.py -q   # 10 passed
```
