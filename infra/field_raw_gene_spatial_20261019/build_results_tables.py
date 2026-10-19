#!/usr/bin/env python3
"""Concise result tables for the descriptive spatial contrast (deliverable, pending review).

Reads only the finished run tables; performs no fitting. Applies the C-01 direction rule
(the signed radius deviation carries no direction) and surfaces the C-02 requirements
(estimable-unit counts, per-source detail, cross-organ marker).
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / 'infra/field_raw_gene_spatial_20261019'
KNOWN_SINGLE_CASE_VARIABILITY = 2.0e-4   # source_0/EMT point case, 1/2 vs 8 effective threads
PRIMARY = [('rec_local_pattern_pearson', 'higher_is_better'),
           ('rec_jaccard_q75', 'higher_is_better'),
           ('rec_jaccard_q80', 'higher_is_better (sensitivity)'),
           ('rec_jaccard_q90', 'higher_is_better (sensitivity)'),
           ('rec_centroid_displacement_pitch', 'lower_is_better'),
           ('rec_patch75_largest_radius_abs_delta', 'lower_is_better'),
           ('rec_sd_ratio', 'descriptive_only'),
           ('rec_rmse', 'lower_is_better'),
           ('rec_centered_rmse', 'lower_is_better')]
MODELS = ['point_ridge', 'point_ridge_smooth_half', 'local_ridge', 'background_ridge',
          'permuted_ridge']
BASELINES = ['point_ridge', 'point_ridge_smooth_half', 'background_ridge', 'permuted_ridge']
SIGNED_RADIUS = ['rec_patch75_largest_radius_pitch_delta',
                 'rec_patch90_largest_radius_pitch_delta']


def rows(path):
    with (W / path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def num(value):
    try:
        return float(value) if value not in (None, '', 'None', 'nan') else None
    except ValueError:
        return None


def fmt(value, digits=4):
    return 'NA' if value is None else f'{value:.{digits}f}'


def main():
    macro = rows('macro_summary.tsv')
    paired = rows('macro_paired.tsv')
    unit_metrics = rows('unit_metrics.tsv')
    unit_paired = rows('unit_paired.tsv')
    denom = rows('denominators.tsv')
    claims = rows('readout_claim_summary.tsv')
    splits = json.loads((W / 'splits.json').read_text())
    cross_organ = {s['name'] for s in splits
                   if s['arm'] == 'source' and s.get('same_organ_train_sections', 1) == 0}
    held_source = {s['name']: s.get('held_source', '') for s in splits}
    per_source = {}
    for r in unit_paired:
        if r['arm'] != 'source' or r['baseline'] != 'point_ridge':
            continue
        key = (r['target'], r['source_split_key'], r['split'])
        per_source.setdefault(key, []).append(r)
    source_rows = []
    for (target, source, split), group in sorted(per_source.items()):
        row = dict(target=target, source=source, split=split, n_units=len(group))
        for metric in ['delta_rec_local_pattern_pearson', 'delta_rec_jaccard_q75',
                       'delta_rec_centroid_displacement_pitch',
                       'delta_rec_patch75_largest_radius_abs_delta', 'delta_rec_rmse']:
            values = [num(g.get(metric)) for g in group]
            values = [v for v in values if v is not None]
            row[metric] = (sum(values) / len(values)) if values else None
            row[f'{metric}__n'] = len(values)
        row['cross_organ_extrapolation'] = split in cross_organ
        source_rows.append(row)
    arm_index = {(r['arm'], r['target'], r['model'], r['metric']): r for r in macro
                 if r['aggregate'] == 'arm'}
    pair_index = {(r['arm'], r['target'], r['baseline'], r['metric']): r for r in paired
                  if r['aggregate'] == 'arm'}
    out = []
    out.append('# 同原始基因输入空间对照：描述性结果表（REVIEWED：只允许描述性报告，不得声称真实邻域排列带来恢复；见 review/result_comparison_review.json）\n')
    out.append('来源：`infra/field_raw_gene_spatial_20261019/` 的 evaluate/summarize 产物；'
               '同一批 500 个预测包、5 个模型；内部 65 单位 / 来源 47 单位；无显著性、无 PASS 阈值。')
    out.append('已知数值边界：single-case（source_0/EMT point）1/2 线程 vs 8 线程的逐点差约 2.1e-4、'
               '每片 RMSE 差约 1.9e-4；该证据只覆盖这一个 case 的 point 模型，'
               '未覆盖 candidate/background/permuted，也不假定同批误差相消。'
               '下表中 |Δ| ≤ 2e-4 的配对差标 `~num`（处在已知数值波动内，结论未确定）。\n')

    out.append('## 1. 覆盖与分母（完整 / 部分 / 缺失单位组）\n')
    out.append('| 臂 | 模型 | 单位组 | 主读出可评 | 部分 | 缺失 |')
    out.append('|---|---|---:|---:|---:|---:|')
    for arm in ['internal', 'source']:
        for model in MODELS:
            sel = [r for r in denom if r['arm'] == arm and r['model'] == model
                   and r['metric'] == 'rec_local_pattern_pearson']
            out.append(f"| {arm} | {model} | {sel[0]['unit_rows']} | {sel[0]['complete']} | "
                       f"{sel[0]['partial']} | {sel[0]['missing']} |")
    out.append('')

    out.append('## 2. 主读出与次读出（臂级单位等权宏平均；内部/来源分开）\n')
    header = ('| 目标 | 臂 | 模型 | local-detail Pearson | q75 Jaccard | q80 / q90 Jaccard | '
              '质心位移(pitch) | 半径绝对偏差(pitch) | sd比 | RMSE | 去均值RMSE |')
    out.append(header)
    out.append('|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for target in sorted({r['target'] for r in macro}):
        for arm in ['internal', 'source']:
            for model in MODELS:
                def v(metric):
                    r = arm_index.get((arm, target, model, metric))
                    return fmt(num(r['mean']) if r else None)
                out.append(f"| {target} | {arm} | {model} | {v('rec_local_pattern_pearson')} | "
                           f"{v('rec_jaccard_q75')} | {v('rec_jaccard_q80')} / {v('rec_jaccard_q90')} | "
                           f"{v('rec_centroid_displacement_pitch')} | "
                           f"{v('rec_patch75_largest_radius_abs_delta')} | {v('rec_sd_ratio')} | "
                           f"{v('rec_rmse')} | {v('rec_centered_rmse')} |")
    out.append('')
    out.append('注：有符号半径差（`rec_patch75/90_largest_radius_pitch_delta`）只作描述，'
               '负值不代表范围改善（C-01）；范围结论一律以“每片绝对值再按完整单位聚合”的读数为准。\n')

    out.append('## 3. 候选减四基线（臂级单位等权配对差；C-01 方向标注）\n')
    out.append('| 目标 | 臂 | 基线 | 判据 | local-detail Δ | q75 Jaccard Δ | 质心位移 Δ | '
               '半径绝对偏差 Δ | RMSE Δ |')
    out.append('|---|---|---|---|---:|---:|---:|---:|---:|')
    for target in sorted({r['target'] for r in macro}):
        for arm in ['internal', 'source']:
            for base_name in BASELINES:
                cells = []
                flags = []
                for metric in ['rec_local_pattern_pearson', 'rec_jaccard_q75',
                               'rec_centroid_displacement_pitch',
                               'rec_patch75_largest_radius_abs_delta', 'rec_rmse']:
                    r = pair_index.get((arm, target, base_name, f'delta_{metric}'))
                    value = num(r['mean']) if r else None
                    if value is None:
                        cells.append('NA')
                    else:
                        flags.append(abs(value) <= KNOWN_SINGLE_CASE_VARIABILITY)
                        cells.append(f"{value:+.4f}" + (' ~num' if abs(value) <= KNOWN_SINGLE_CASE_VARIABILITY else ''))
                n_est = num(pair_index[(arm, target, base_name, 'delta_rec_rmse')]['n_estimable']) \
                    if (arm, target, base_name, 'delta_rec_rmse') in pair_index else None
                criterion = f"{'~num' if any(flags) else 'beyond-2e-4'}（n_est={int(n_est) if n_est else 'NA'}）"
                out.append(f"| {target} | {arm} | {base_name} | {criterion} | " +
                           ' | '.join(cells) + ' |')
    out.append('')
    out.append('`~num`：该差落在单一 case 已知数值波动（2e-4）之内，按审核要求标为未确定；'
               '`beyond-2e-4` 仅表示数值上超出该单 case 波动参考，不是显著性。\n')

    out.append('## 4. 来源分层（来源臂，按 held_source；乳腺/肝脏为器官外推）\n')
    out.append('派生自 `unit_paired.tsv`（候选减同点，单位等权；每个来源内单位平均），')
    out.append('四目标 × 6 来源的完整行见 `RESULTS_SOURCE_TABLE.tsv`。\n')
    out.append('| 目标 | 来源 | 单位组(可评) | local-detail Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ | 器官外推 |')
    out.append('|---|---|---:|---:|---:|---:|---:|---|')
    for r in source_rows:
        label = r['source'][:44] if r['source'] else '(empty)'
        out.append(f"| {r['target']} | {label} | {r['n_units']}"
                   f"({r['delta_rec_local_pattern_pearson__n']}) | "
                   f"{fmt(r['delta_rec_local_pattern_pearson'], 5)} | "
                   f"{fmt(r['delta_rec_jaccard_q75'], 5)} | "
                   f"{fmt(r['delta_rec_centroid_displacement_pitch'], 5)} | "
                   f"{fmt(r['delta_rec_patch75_largest_radius_abs_delta'], 5)} | "
                   f"{'是' if r['cross_organ_extrapolation'] else '否'} |")
    out.append('')
    with (W / 'RESULTS_SOURCE_TABLE.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(source_rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(source_rows)
    out.append('注：分来源必须在多目标、多基线下逐一核对（完整分来源汇总可再按需派生），'
               '不得把这个分层当成显著性或成功率（C-02）。\n')

    out.append('## 5. 描述性方向合取（不是阈值、不是显著性）\n')
    out.append('| 臂 | 目标 | 读出 | 方向 | 同时优于四个对照 | 最小可评单位数 |')
    out.append('|---|---|---|---|---|---:|')
    for c in claims:
        out.append(f"| {c['arm']} | {c['target']} | {c['readout']} | {c['direction']} | "
                   f"{c['favourable_vs_all_four']} | {c.get('min_n_estimable_units', 'NA')} |")
    out.append('')
    out.append('## 6. 限制与未确定项\n')
    out.append('- 未知患者：内部 72 片中 14 片（8 单位）为研究内确认身份，其余 58 片（57 单位）未知；'
               '切片数不冒称患者数。')
    out.append('- 来源身份不足的 19 片只保留内部结果，不进入来源臂；来源臂 47 单位覆盖 6 个来源。')
    out.append('- 乳腺（source_5）与肝脏（source_0）为器官外推，未与同器官来源混合解释。')
    out.append('- 冻结基线复现：11/44 case 超出旧逐点 1e-4 / 每片 RMSE 1e-5 界，'
               '按 D-275 保留为已披露的失败，不改写成通过、不放宽界限。')
    out.append('- 数值波动证据只覆盖 source_0/EMT 的 point 模型（单 case）；'
               '不推广到 candidate/background/permuted，也不假定同批误差相消。')
    out.append('- 本轮 NOT_RUN：人工 TLS、独立 marker、正式显著性、物理尺度、视野外推断、'
               'shuffle 参考（本轮关闭）。')
    (W / 'RESULTS_TABLES.md').write_text('\n'.join(out) + '\n')
    print('\n'.join(out[:40]))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
