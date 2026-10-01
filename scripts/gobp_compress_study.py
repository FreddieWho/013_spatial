#!/usr/bin/env python3
"""GOBP set compression study: size floor + ontology hierarchy + gene overlap.

Inputs (all local, no download):
  --gmt   /home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt  (MSigDB v2025.1)
  --obo   /home/huyudi/012_conference/aaai2027/data/raw/ontology/go-basic.obo (releases/2026-06-15)
  --defs  infra/r16/recovery_20260918/programs_go/definitions.json (6987 cohort-covered programs)

Reports compression curves only; it selects no frozen panel by itself.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/gobp_compress_20261001"
DEFAULT_GMT = Path("/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt")
DEFAULT_OBO = Path("/home/huyudi/012_conference/aaai2027/data/raw/ontology/go-basic.obo")
DEFAULT_DEFS = ROOT / "infra/r16/recovery_20260918/programs_go/definitions.json"
ROOT_TERM = "GO:0008150"


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def norm(name: str) -> str:
    s = name.strip().upper()
    s = re.sub(r"[^A-Z0-9]+", "_", s)
    return re.sub(r"_+", "_", s).strip("_")


def read_gmt(path):
    sets = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = sorted({g for g in parts[2:] if g})
    return sets


def parse_obo(path, namespace="biological_process"):
    """Return (terms, parents) for non-obsolete terms of the requested namespace."""
    terms, parents, cur = {}, defaultdict(set), None

    def flush():
        nonlocal cur
        if cur and cur.get("id"):
            terms[cur["id"]] = cur
        cur = None

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("["):
                flush()
                if line == "[Term]":
                    cur = {"id": None, "name": None, "namespace": None, "obsolete": False,
                           "alt_ids": [], "replaced_by": [], "consider": []}
                continue
            if cur is None or not line:
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                cur["id"] = val.strip()
            elif key == "name":
                cur["name"] = val.strip()
            elif key == "namespace":
                cur["namespace"] = val.strip()
            elif key == "is_obsolete" and val.strip() == "true":
                cur["obsolete"] = True
            elif key == "alt_id":
                cur["alt_ids"].append(val.strip())
            elif key == "replaced_by":
                cur["replaced_by"].append(val.strip())
            elif key == "consider":
                cur["consider"].append(val.strip())
            elif line.startswith("is_a:"):
                parents[cur["id"]].add(line.split()[1])
            elif line.startswith("relationship: part_of"):
                parents[cur["id"]].add(line.split()[2])
    if cur and cur["id"]:
        terms[cur["id"]] = cur
    keep = {t for t, d in terms.items()
            if d["namespace"] == namespace and not d["obsolete"] and d["id"].startswith("GO:")}
    parents = {t: {p for p in parents.get(t, set()) if p in keep} for t in keep}
    return {t: terms[t] for t in keep}, parents


def depths_and_ancestors(parents, root=ROOT_TERM):
    """depth = shortest path to root; ancestors = all transitive parents."""
    depth, anc = {}, {}

    def walk(t, seen):
        if t in depth:
            return depth[t], anc[t]
        if t in seen:
            return 10 ** 6, set()
        seen = seen | {t}
        ps = parents.get(t, set())
        if not ps:
            depth[t], anc[t] = (0 if t == root else 10 ** 6), set()
            return depth[t], anc[t]
        best, acc = 10 ** 6, set()
        for p in ps:
            dp, ap = walk(p, seen)
            best = min(best, dp + 1)
            acc |= ap | {p}
        depth[t], anc[t] = best, acc
        return best, anc[t]

    for t in parents:
        walk(t, frozenset())
    return depth, anc


def compress2(names, genes, anc_by_prog, goid_by_prog, depth_by_prog, size_floor, jaccard,
              order="desc", depth_window=None, size_cap=None):
    """Generalised greedy reduction.

    order='desc' keeps broad terms and removes their descendants (keep-general).
    order='asc'  keeps the most specific term in each lineage (keep-specific).
    depth_window=(lo, hi) restricts candidates to an ontology depth band.
    """
    cands = [n for n in names
             if (depth_window is None
                 or (depth_by_prog.get(n) is not None
                     and depth_window[0] <= depth_by_prog[n] <= depth_window[1]))
             and (size_cap is None or len(genes[n]) <= size_cap)]
    sign = -1 if order == "desc" else 1
    order_list = sorted(cands, key=lambda n: (sign * len(genes[n]), n))
    fsets = {n: frozenset(genes[n]) for n in names}
    kept, kept_fsets, reasons = [], [], {}
    kept_ids, kept_anc_union, inv = set(), set(), defaultdict(list)
    for n in order_list:
        g = fsets[n]
        if len(g) < size_floor:
            continue
        tid = goid_by_prog.get(n)
        a = anc_by_prog.get(n, set())
        if order == "desc":
            if a & kept_ids:
                reasons[n] = "descendant"
                continue
        else:
            if tid in kept_anc_union:
                reasons[n] = "ancestor"
                continue
        drop = None
        for j in _cand_idx(g, inv):
            h = kept_fsets[j]
            if len(g) * jaccard > len(h) or len(h) * jaccard > len(g):
                continue
            if len(g & h) >= jaccard * len(g | h):
                drop = f"jaccard>={jaccard}"
                break
        if drop:
            reasons[n] = drop
            continue
        kept.append(n)
        kept_fsets.append(g)
        if tid:
            kept_ids.add(tid)
            kept_anc_union |= a
        for gene in g:
            inv[gene].append(len(kept) - 1)
    return kept, reasons


def _cand_idx(g, inv):
    out = set()
    for gene in g:
        out.update(inv[gene])
    return out


def compress(names, genes, anc_by_prog, goid_by_prog, depth, size_floor, jaccard,
             keep_general=True):
    """Greedy reduction: larger sets first. O(n * (|anc| + overlaps))."""
    order = sorted(names, key=lambda n: (-len(genes[n]), n))
    fsets = {n: frozenset(genes[n]) for n in names}
    kept, kept_fsets, reasons = [], [], {}
    kept_ids, kept_anc_union, inv = set(), set(), defaultdict(list)
    for n in order:
        g = fsets[n]
        if len(g) < size_floor:
            reasons[n] = f"size<{size_floor}"
            continue
        tid = goid_by_prog.get(n)
        a = anc_by_prog.get(n, set())
        if keep_general:
            if a & kept_ids:
                reasons[n] = "descendant"
                continue
        else:
            if tid in kept_anc_union:
                reasons[n] = "ancestor"
                continue
        cand = set()
        for gene in g:
            cand.update(inv[gene])
        drop = None
        for j in cand:
            h = kept_fsets[j]
            if len(g) * jaccard > len(h) or len(h) * jaccard > len(g):
                continue
            if len(g & h) >= jaccard * len(g | h):
                drop = f"jaccard>={jaccard}"
                break
        if drop:
            reasons[n] = drop
            continue
        kept.append(n)
        kept_fsets.append(g)
        if tid:
            kept_ids.add(tid)
            kept_anc_union |= a
        for gene in g:
            inv[gene].append(len(kept) - 1)
    return kept, reasons


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gmt", type=Path, default=DEFAULT_GMT)
    ap.add_argument("--obo", type=Path, default=DEFAULT_OBO)
    ap.add_argument("--defs", type=Path, default=DEFAULT_DEFS)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    gmt = {n: g for n, g in read_gmt(args.gmt).items() if n.startswith("GOBP_")}
    print(f"gobp sets={len(gmt)}", flush=True)
    defs = json.loads(args.defs.read_text())
    covered = {pid: int(d["covered_size"]) for pid, d in defs.items()}
    terms, parents = parse_obo(args.obo)
    print(f"obo bp terms={len(terms)}", flush=True)
    depth, anc = depths_and_ancestors(parents)
    n_fin = sum(1 for v in depth.values() if v < 10 ** 6)
    print(f"depths finite={n_fin}/{len(depth)}", flush=True)

    by_norm = defaultdict(list)
    for tid, d in terms.items():
        by_norm[norm(d["name"])].append(tid)
    matched, unmatched, ambiguous = {}, [], 0
    for name in gmt:
        cand = by_norm.get(norm(name.replace("GOBP_", "", 1)))
        if not cand:
            unmatched.append(name)
            continue
        if len(cand) > 1:
            ambiguous += 1
        matched[name] = sorted(cand)[0]

    prov = {"gmt": str(args.gmt), "gmt_sha256": sha256(args.gmt),
            "obo": str(args.obo), "obo_sha256": sha256(args.obo),
            "obo_data_version": next((l.split(": ", 1)[1].strip()
                                      for l in open(args.obo, encoding="utf-8")
                                      if l.startswith("data-version")), ""),
            "defs": str(args.defs), "defs_sha256": sha256(args.defs)}
    n_bp = len(terms)
    res = {"provenance": prov, "gmt_sets": len(gmt), "obo_bp_terms": n_bp,
           "matched": len(matched), "unmatched": len(unmatched), "ambiguous_names": ambiguous,
           "unmatched_examples": unmatched[:20]}
    (OUT / "match_report.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))

    # depth distribution of matched sets
    dvals = [depth[t] for t in matched.values() if depth.get(t, 10 ** 6) < 10 ** 6]
    if not dvals:
        raise SystemExit(f"depth computation produced no finite values (matched={len(matched)}, bp_terms={len(terms)})")
    hist = defaultdict(int)
    for v in dvals:
        hist[v] += 1
    with (OUT / "level_distribution.tsv").open("w") as f:
        f.write("depth\tn\n")
        for k in sorted(hist):
            f.write(f"{k}\t{hist[k]}\n")

    # size-floor curve (raw MSigDB size)
    rows = []
    universe = set()
    for g in gmt.values():
        universe |= set(g)
    for floor in (3, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100):
        keep = [n for n in gmt if len(gmt[n]) >= floor]
        u = set()
        for n in keep:
            u |= set(gmt[n])
        cov_ok = sum(1 for n in keep if covered.get(n, 0) >= 10)
        rows.append({"rule": f"size>={floor}", "n_kept": len(keep),
                     "gene_union": len(u), "gene_union_frac": round(len(u) / len(universe), 4),
                     "kept_with_covered10": cov_ok})
    with (OUT / "size_floor_curve.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # combined rules
    anc_by_prog = {p: anc.get(t, set()) for p, t in matched.items()}
    combos = []
    for floor in (5, 10, 15, 20, 30):
        for jac in (0.5, 0.7, 0.9):
            keep, reasons = compress(list(gmt), gmt, anc_by_prog, matched, depth, floor, jac)
            n_size = sum(1 for r in reasons.values() if r.startswith("size<"))
            n_hier = sum(1 for r in reasons.values() if r.startswith("descendant"))
            n_jac = sum(1 for r in reasons.values() if r.startswith("jaccard"))
            combos.append({"size_floor": floor, "jaccard": jac, "n_kept": len(keep),
                           "removed_total": len(gmt) - len(keep),
                           "removed_size": n_size, "removed_hierarchy": n_hier,
                           "removed_jaccard": n_jac,
                           "kept_covered10": sum(1 for n in keep if covered.get(n, 0) >= 10),
                           "median_raw_size": int(sorted(len(gmt[n]) for n in keep)[len(keep) // 2])})
    with (OUT / "compression_curve.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(combos[0]), delimiter="\t")
        w.writeheader()
        w.writerows(combos)

    # detail for one default combination (reported as the study default, not frozen)
    # rule matrix: keep-general vs keep-specific vs depth-window
    depth_by_prog = {p: depth.get(t) for p, t in matched.items()}
    rules = []
    settings = [("general", "desc", None, None), ("specific", "asc", None, None),
                ("window_3_6", "desc", (3, 6), None),
                ("bal_spec_d5_c500", "asc", (5, 20), 500),
                ("bal_gen_d5_c500", "desc", (5, 20), 500),
                ("bal_spec_d5_c300", "asc", (5, 20), 300)]
    for rname, order, win, cap in settings:
        for floor in (10, 15, 20, 30):
            keep, _ = compress2(list(gmt), gmt, anc_by_prog, matched, depth_by_prog,
                                floor, 0.7, order=order, depth_window=win, size_cap=cap)
            sizes = sorted(len(gmt[n]) for n in keep)
            in_defs = [n for n in keep if n in defs]
            c10 = sum(1 for n in in_defs if covered.get(n, 0) >= 10)
            dep = sorted(depth_by_prog[n] for n in keep if depth_by_prog.get(n) is not None)
            rules.append({"rule": rname, "size_floor": floor, "size_cap": cap or "", "jaccard": 0.7,
                          "n_kept": len(keep), "median_size": sizes[len(sizes) // 2],
                          "p90_size": sizes[int(0.9 * (len(sizes) - 1))],
                          "median_depth": dep[len(dep) // 2] if dep else "",
                          "in_defs": len(in_defs), "covered10": c10,
                          "kept_frac_of_7583": round(len(keep) / len(gmt), 3),
                          "kept_frac_of_6987": round(len(in_defs) / max(1, len(defs)), 3)})
    with (OUT / "rule_matrix.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rules[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rules)
    print("\nrule matrix (jaccard 0.7):")
    print(f"{'rule':>18} {'floor':>5} {'cap':>4} {'kept':>6} {'med':>5} {'p90':>5} {'dep':>4} {'in_defs':>7} {'cov10':>6} {'%7583':>6} {'%6987':>6}")
    for r in rules:
        print(f"{r['rule']:>18} {r['size_floor']:5} {str(r['size_cap']):>4} {r['n_kept']:6} {r['median_size']:5} "
              f"{r['p90_size']:5} {str(r['median_depth']):>4} {r['in_defs']:7} {r['covered10']:6} "
              f"{r['kept_frac_of_7583']:6} {r['kept_frac_of_6987']:6}")
    keep_win, reasons_win = compress2(list(gmt), gmt, anc_by_prog, matched, depth_by_prog,
                                      15, 0.7, order="desc", depth_window=(3, 6))
    keep_bal, reasons_bal = compress2(list(gmt), gmt, anc_by_prog, matched, depth_by_prog,
                                      20, 0.7, order="asc", depth_window=(5, 20), size_cap=500)
    keep_spec, reasons_spec = compress2(list(gmt), gmt, anc_by_prog, matched, depth_by_prog,
                                        20, 0.7, order="asc", depth_window=None, size_cap=None)
    with (OUT / "survivors_specific_floor20_jac70.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "raw_size", "covered_size", "in_definitions", "depth", "rule"])
        for n in sorted(keep_spec):
            w.writerow([n, len(gmt[n]), covered.get(n, ""), int(n in defs),
                        depth_by_prog.get(n, ""), "keep_specific_asc_floor20_jac70"])
    sp_sizes = sorted(len(gmt[n]) for n in keep_spec)
    print(f"specific panel: n={len(keep_spec)} min={sp_sizes[0]} med={sp_sizes[len(sp_sizes)//2]} "
          f"p90={sp_sizes[int(0.9*(len(sp_sizes)-1))]} max={sp_sizes[-1]}")
    with (OUT / "survivors_balanced_spec_d5_c500_floor20.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "raw_size", "covered_size", "in_definitions", "depth", "n_ancestors"])
        for n in sorted(keep_bal, key=lambda x: -len(gmt[x])):
            w.writerow([n, len(gmt[n]), covered.get(n, ""), int(n in defs),
                        depth_by_prog.get(n, ""), len(anc_by_prog.get(n, set()))])
    with (OUT / "survivors_window3_6_floor15_jac70.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "raw_size", "covered_size", "in_definitions", "depth"])
        for n in sorted(keep_win, key=lambda x: -len(gmt[x])):
            w.writerow([n, len(gmt[n]), covered.get(n, ""), int(n in defs),
                        depth_by_prog.get(n, "")])

    keep15, reasons15 = compress(list(gmt), gmt, anc_by_prog, matched, depth, 15, 0.7)
    goids_kept = {matched[o] for o in keep15 if o in matched}
    with (OUT / "survivors_size15_jac70.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "raw_size", "covered_size", "depth", "n_ancestors", "ancestor_kept"])
        for n in sorted(keep15, key=lambda x: -len(gmt[x])):
            tid = matched.get(n)
            a = anc.get(tid, set()) if tid else set()
            tek = [x for x in a if x in goids_kept]
            w.writerow([n, len(gmt[n]), covered.get(n, ""),
                        depth.get(tid, "") if tid else "", len(a), len(tek)])
    with (OUT / "drop_reasons_size15_jac70.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "raw_size", "reason"])
        for n, r in sorted(reasons15.items(), key=lambda kv: -len(gmt[kv[0]])):
            w.writerow([n, len(gmt[n]), r])

    print(json.dumps({"gmt_sets": len(gmt), "matched": len(matched),
                      "unmatched": len(unmatched), "ambiguous": ambiguous,
                      "bp_terms": n_bp}, ensure_ascii=False))
    print("\ncombined rules (rows = size floor x jaccard):")
    print(f"{'floor':>5} {'jac':>4} {'kept':>6} {'removed':>8} {'by_size':>8} {'by_hier':>8} {'by_jac':>7}")
    for r in combos:
        print(f"{r['size_floor']:5} {r['jaccard']:4} {r['n_kept']:6} {r['removed_total']:8} "
              f"{r['removed_size']:8} {r['removed_hierarchy']:8} {r['removed_jaccard']:7}")


if __name__ == "__main__":
    main()
