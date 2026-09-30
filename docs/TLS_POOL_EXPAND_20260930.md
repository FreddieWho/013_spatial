# TLS-signature pool expansion：跨队列/平台/组织（2026-09-30）

只数池子规模与诊断，不拟合梯度、不反推、不算 p。合同见
`infra/tls_pool_expand_20260930/run_contract_v3.json`，脚本
`scripts/tls_pool_expand.py`，产物 `pool_diagnostics.tsv` /
`pool_components.tsv` / `pool_summary.json`。

## 口径（一句话）

同一签名（B 轴 49 + 浆细胞 7 + CXCL13/CCL19/CCL21/LTB，共 57 符号）、
同一 AUCell 算子、同一 top-10% 分位阈值、同一六边连通组件规则，
按平台换几何（Visium 系 100 μm；VisiumHD/Xenium 先聚到 55 μm 伪格）。
重复符号跨行加总（D-148 路线）。每 section 另跑 10 次值置换 null。

## 池子规模（55/55 section 全过）

| 队列 | section 数 | 组件总数 | 备注 |
|---|---|---|---|
| USZ（肾/肺，Visium v1/v2） | 8 | 430 | 标签验证见下 |
| GSE175540（ccRCC，Visium v1） | 18 | 1008 | 标签验证见下 |
| ST-CRC（CRC，Visium） | 14 | 364 | 无 TLS 标签；IC aggregate 类 FAIL-CLOSED |
| HTAN-CRC（淋巴滤泡验证 8 片） | 8 | 452 | 标签验证见下 |
| Cervilla CRC 对（Visium v1 + CytAssist v2） | 2 | 76 | 无标签 |
| 10x 公开 CRC（Parent Visium + CytAssist 11mm） | 2 | 250 | 无标签 |
| VisiumHD 结肠 FF（16 μm→55 μm 伪格） | 1 | 125 | 无标签 |
| Xenium 肾/肺（55 μm 伪格，panel 受限） | 2 | 153 | 无标签，见下 |
| 合计 | 55 | 2859 | |

## 签名覆盖

Visium 系全部 50–56/57。有缺口的是 Xenium（20/57：缺 CXCL13/CCL21/LTB/
JCHAIN/SDC1/CD74/HLA-DRA/COL1A1 等）——Xenium 是独立的 panel 受限层，
只报告、不与 Visium 拼数。

## 标签验证（有标签的队列）

- USZ：TLS vs NOR，AUC 0.63–0.94；LC4 0.38、LC5 0.18（已知标签异常，D-157）。
- GSE175540：TLS vs NO_TLS，AUC 0.56–0.96（17/18 ≥ 0.70）。
- HTAN 8 片：lymphoid_follicle vs 其余，AUC 0.66–0.999（7003_AS_5 0.999）。

## 置换 null 对照（每 section 10 次）

- Visium 系：obs/null 中位约 1.3–3.5（HD 12.5）。高分点成团显著多于随机。
- Cervilla CytAssist v2：0.95，基本等于 null（探针/化学敏感性备注）。
- Xenium：0.72/0.60——观测到的组件数**少于** null 中位，即大团块而非散点；
  panel 受限、无标签验证，不解读为 TLS 证据。

## 不能说的话

- 这些是 B-签名热点，不是已验证 TLS（D-157 同样结论）。
- 数量不能跨平台直接比（Xenium 缺基因、CytAssist 探针不同、HD 是伪格）。
- ST-CRC 的 IC aggregate 类维持 FAIL-CLOSED，不计入任何 TLS 口径。
- LC4/LC5（USZ）与 GSE175540 frozen_c_5（AUC 0.56）等弱验证行保留原样，
  不删不修。
