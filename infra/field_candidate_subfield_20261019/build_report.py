#!/usr/bin/env python3
"""Build REPORT.md for the two-component candidate-field probe from the summary tables."""
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'infra/field_candidate_subfield_20261019'


def read_tsv(path):
    import csv
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def fmt(value, digits=3):
    return f'{float(value):.{digits}f}' if value not in ('', None) and np.isfinite(float(value)) else '–'


def main():
    macro = read_tsv(OUT / 'summary_macro.tsv')
    unit = read_tsv(OUT / 'summary_per_unit.tsv')
    section = read_tsv(OUT / 'summary_per_section.tsv')
    status = json.loads((OUT / 'summary_status.json').read_text())
    cells = sorted({(r['target'], r['split']) for r in macro} | {(k.split('|')[0], k.split('|')[1]) for k in status})
    lines = ['# 两分量候选场探针（D-279）：结果报告', '']
    lines += ['## 口径与结论上限', '',
              '- 四目标全报；候选为逐折固定 K=2 的两分量非负分解，A/B 按基因名与种子 0 拆分。',
              '- A、B 图共用“排除本目标整个面板”的外部文库分母；训练单位冻结可测 A、等量抽样与缩放；测试 B 仅事后评价。',
              '- 普通参照（A 均值、K=1、单基因、全 A Ridge）与候选同折同 A 同训练抽样；全 A Ridge 是更强参照。',
              '- **即使稳定、非冗余、受 B 支持，也只称候选场**：不是独立病理真值，不是独立生物场，不是 HF-01；',
              '  旧 2026-10-06 FC-03 分数与氧化磷酸化优势不计入本轮。', '']
    counts = Counter(v.get('status') for v in status.values())
    lines += ['## 运行状态', '', f'- 目标×划分单元：{len(cells)}（期望 {len(status)}）',
              f'- 状态分布：{dict(counts)}', '']
    not_testable = [(k, v) for k, v in status.items() if v.get('status') != 'OK']
    if not_testable:
        lines += ['| 单元 | 状态 | 可用 A | 说明 |', '|---|---|---|---|']
        for key, value in sorted(not_testable):
            lines.append(f"| {key} | {value.get('status')} | {value.get('usable_A')} | "
                         f"`{value.get('nmf_k2', {}).get('status', '')}` |")
        lines.append('')
    lines += ['## 各目标：候选与参照（完整单位等权，B 基因逐单位汇总）', '']
    for target in sorted({r['target'] for r in macro} | {c[0] for c in cells}):
        lines += [f'### {target}', '',
                  '| 臂 | 可评单位数（范围） | 基因×划分行 | Pearson 中位 | Jaccard 中位 | Pearson>0 比例 |',
                  '|---|---|---|---|---|---|']
        subset = [r for r in macro if r['target'] == target]
        for arm in ['A_mean', 'K1', 'single_gene', 'K2', 'full_A_ridge']:
            rows = [r for r in subset if r['arm'] == arm]
            pears = [float(r['pearson']) for r in rows if r['pearson']]
            jacs = [float(r['jaccard']) for r in rows if r['jaccard']]
            n_units = sorted({int(r['n_units']) for r in rows}) if rows else []
            pos = np.mean([p > 0 for p in pears]) if pears else float('nan')
            unit_span = (f'{min(n_units)}–{max(n_units)}' if n_units else '–')
            lines.append(f'| {arm} | {unit_span} | {len(rows)} | '
                         f'{fmt(np.median(pears)) if pears else "–"} | '
                         f'{fmt(np.median(jacs)) if jacs else "–"} | '
                         f'{fmt(pos, 2) if pears else "–"} |')
        lines.append('')
    lines += ['## K=2 相对参照的配对差（同目标同划分同 B 基因）', '',
              '| 目标 | 对比 | 配对数 | 中位差 | K2 更高比例 |', '|---|---|---|---|---|']
    index = {(r['target'], r['split'], r['gene'], r['arm']): r for r in macro}
    for target in sorted({r['target'] for r in macro}):
        for reference in ('A_mean', 'full_A_ridge'):
            deltas = []
            for (t, sp, gene, arm) in list(index):
                if t != target or arm != 'K2':
                    continue
                other = index.get((t, sp, gene, reference))
                here = index[(t, sp, gene, 'K2')]
                if not other or not other['pearson'] or not here['pearson']:
                    continue
                deltas.append(float(here['pearson']) - float(other['pearson']))
            if deltas:
                lines.append(f'| {target} | K2 − {reference} | {len(deltas)} | '
                             f'{fmt(np.median(deltas))} | {fmt(np.mean([d > 0 for d in deltas]), 2)} |')
    lines += ['', '## 内部与来源留出（完整单位等权，Pearson 中位）', '',
              '| 目标 | 臂 | 内部 5 折 | 来源 6 屏 |', '|---|---|---|---|']
    for target in sorted({r['target'] for r in macro}):
        for arm in ['A_mean', 'K1', 'K2', 'full_A_ridge', 'single_gene']:
            row = []
            for prefix in ('internal', 'source'):
                vals = [float(r['pearson']) for r in macro if r['target'] == target and r['arm'] == arm
                        and r['split'].startswith(prefix) and r['pearson']]
                row.append(fmt(np.median(vals)) if vals else '–')
            lines.append(f'| {target} | {arm} | {row[0]} | {row[1]} |')
    lines += ['', '## 分量稳定性与技术相关', '']
    stab = defaultdict(list)
    sim = defaultdict(list)
    for key, value in status.items():
        target = key.split('|')[0]
        for item in value.get('stability', []) or []:
            if item.get('status') == 'OK':
                stab[target].append(float(item['cos_sum']) / 2.0)
        if value.get('similarity', {}).get('k2_component_cosine') is not None:
            sim[target].append(float(value['similarity']['k2_component_cosine']))
    lines += ['| 目标 | 折内重抽样配对余弦（中位，按分量均） | 两分量载荷余弦（中位） |', '|---|---|---|']
    for target in sorted(set(stab) | set(sim)):
        lines.append(f'| {target} | {fmt(np.median(stab[target])) if stab[target] else "–"} | '
                     f'{fmt(np.median(sim[target])) if sim[target] else "–"} |')
    tech = []
    for path in sorted((OUT / 'results/main').glob('*/*/technical.tsv')):
        for row in read_tsv(path):
            row['target'] = path.parent.parent.name
            tech.append(row)
    lines += ['', '| 目标 | 分量 | 深度相关（中位） | 检出率相关（中位） |', '|---|---|---|---|']
    for target in sorted({r['target'] for r in tech}):
        for comp in ('0', '1'):
            rows = [r for r in tech if r['target'] == target and r['component'] == comp]
            depths = [abs(float(r['depth_corr'])) for r in tech if r['target'] == target
                      and r['component'] == comp and r['depth_corr'] not in ('', 'nan')]
            dets = [abs(float(r['detection_corr'])) for r in tech if r['target'] == target
                    and r['component'] == comp and r['detection_corr'] not in ('', 'nan')]
            lines.append(f'| {target} | {comp} | {fmt(np.median(depths)) if depths else "–"} | '
                         f'{fmt(np.median(dets)) if dets else "–"} |')
    lines += ['', '## 失败与不可评', '',
              f'- 逐片行 {len(section)}，逐单位行 {len(unit)}，完整单位汇总行 {len(macro)}。',
              f'- 不可评（单位级）：{sum(1 for r in unit if r["status"] != "OK")}。',
              '- 稀疏/并列导致的逐片不可评见 `summary_per_section.tsv` 的 `status` 列。', '']
    lines += ['## 边界', '',
              '- 片内同向与 B 支持不等于生物协同、独立病理真值或 HF-01；片数不是患者数。',
              '- 本探针不检验场间空间恢复；隐藏场恢复需在候选冻结后另做，仍保留同点与整片背景对照。',
              '- 人工 TLS 标签、旧人工评价、HookNet 全轮封存。', '']
    (OUT / 'REPORT.md').write_text('\n'.join(lines) + '\n')
    print(f'wrote {OUT / "REPORT.md"} ({len(lines)} lines)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
