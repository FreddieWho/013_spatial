# T3 triage log（2026-09-18，主 agent 记录决策理由；零入围也有证据链）

## 输入
- laneB 191 全记录（`laneB_residual_20260917.json`）
- Tier-2 13 复现组 markers（R16P003/004/005/007/008/009/010/013/014/015/025/056/058）
- 联合臂 A/B（216 候选）、G（356）、H（83）：npat>=2 者 top10 markers
- 六轴定义（T/B/Mye/Epi/Stromal/Plasma；ILC 退役）

## 去重（candidate_union.tsv：556 基因）
- 191 中 140 个在 A/B 出现、113 在 H、96 在 G、48 在 Tier-2 十三组——各方法捞同一批，
  投票数不作独立证据数（审计 A17 已警告口径差异：H 78/83 vs G 32/356 是口径差异非优劣）。
- 191 中 42 个是轴定义基因（leave-out 在程序级处理）；33 个 laneB 独有（无 joint/tier2 支持）。

## 共表达（coexpression_R.npz + modules.tsv）
- 患者内 median → 30 患者等权平均；平均链接 1-|r| 切 0.5。
- 结果：1 个巨模块 M213（159 基因，off-diag |r|中位 0.59）＋ Ig 模块 M93（11）＋ 基质模块 M132（9）
  ＋ 平滑肌 M139（3）＋ MHC2 M123（2）＋ 应激 M197（2）＋ MT M187（2）＋ 368 单例。
- M213 回归六轴后残余协同中位 0.18（原 0.59）：大部分是组成搭便车，但残余非零——
  按 04-B 分层：identity/state-compatible 保留，不作新场解读，不判死。

## 程序冻结（6 个，program_definitions.json）
| 程序 | 基因 | 解释标签 | 可分性 | 去留理由 |
|---|---|---|---|---|
| P-epi | 16 上皮分化 | identity/state | SEPARABLE（Epi轴loo 0.996） | 留：跨30病人QC残余0.34，杯状/肠上皮两态 |
| P-stromal | 9 基质 | identity/state | SEPARABLE（Stromal轴loo 0.971） | 留：QC残余0.22，非线性敏感（quad砍半）如实记 |
| P-plasma | 11 浆/Ig | identity | NOT_SEPARABLE（Plasma轴loo 0.815） | 留但降级：与Plasma轴不可分，T4只验分子再现，不争增量 |
| P-smmhc | 3 平滑肌 | identity/state | SEPARABLE | 留：小而干净（|r|0.52），DES/MYH11联合臂独立支持 |
| P-mhc2 | 2 CD74/HLA-DRA | state | SEPARABLE | 留：唯一非丰度候选（抗原呈递活性），QC残余0.14最弱，如实记 |
| P-stress | 2 FOS系 | state/technical | SEPARABLE（无轴重叠） | 留：FOS/FOSB跨25+病人残余相干；技术梯度未排除，两标签并列 |

- 未入围：M187（MT1E/MT2A，纯技术线粒体，按04-B技术伪影排除）；368 单例（无协同证据）；
  33 laneB独有基因中非持家者（如 AREG、EREG 缺席——实际无强单例，名单见 union 表 n_sources=1 行）。
- HIST1H 簇（组蛋白，M213内）：归入 P-epi 的增殖亚态，不单独成程序（增殖轴未建，方向待定）。

## 非线性敏感性（nonlinear_sensitivity.json）
- QC_lin→QC_quad：P-stromal 被砍最多（0.35→0.17），P-stress 几乎不动（0.28→0.25 / 0.61→0.60）——
  基质残余部分是线性模型欠拟合，应激残余不是。这是"模型缺陷≠发现"的活证据，已记入程序卡。

## 患者效应（program_patient_effects.tsv：180 行=6程序x30病人，全表无截断）
- QC残余>0.1 病人数：epi 30、plasma 29、stress 25、stromal 24、smmhc 20、mhc2 20。
- readout 基因划分：T4 待定（与 input 分离，机器可检查；冻结时填）。
