#!/usr/bin/env python3
"""Cross-cohort / cross-platform / cross-tissue TLS-signature pool expansion.

Frozen operator: infra/tls_pool_expand_20260930/run_contract_v3.json.
Pool sizes and diagnostics ONLY. No gradients, no inversion, no p-values.

Usage: --smoke (loader sanity on first sections) | --run (all contract sections)
Env: LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import tarfile
import zipfile
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_pool_expand_20260930"
CONTRACT_PATH = OUT / "run_contract_v3.json"


def decode(values):
    return [v.decode() if isinstance(v, bytes) else v for v in values]


def contract():
    return json.loads(CONTRACT_PATH.read_text())


def signature_symbols():
    proxy = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
    from r16 import axes as A
    genes = set(proxy["classes"]["B"]["voted_genes"]) | set(A.PLASMA_GENES)
    genes |= {"CXCL13", "CCL19", "CCL21", "LTB"}
    return sorted(genes)


def aucell(X, cols, max_frac=0.05, chunk=600):
    """Rank 1 = highest expression (D-157 implementation, copied verbatim)."""
    from scipy.stats import rankdata
    n, k = X.shape[0], len(cols)
    out = np.full(n, np.nan)
    if k == 0:
        return out
    max_rank = max(1, int(np.floor(max_frac * X.shape[1])))
    for a in range(0, n, chunk):
        b = min(n, a + chunk)
        block = X[a:b].toarray() if sparse.issparse(X) else X[a:b]
        r = rankdata(-block, axis=1, method="average")
        pos = np.clip(max_rank - r[:, cols] + 1.0, 0.0, None)
        out[a:b] = pos.sum(axis=1) / (k * max_rank)
    return out


def hex_xy(array_row, array_col, pitch_um):
    array_row = np.asarray(array_row, dtype=float)
    array_col = np.asarray(array_col, dtype=float)
    return np.column_stack([array_col + 0.5 * (array_row % 2),
                            array_row * (3 ** 0.5) / 2]) * pitch_um


def components(xy, mask, pitch_um, min_size=2):
    if len(np.asarray(mask)) != len(xy):
        raise ValueError(f"components(): mask length {len(np.asarray(mask))} != coordinate rows {len(xy)}")
    idx = np.where(np.asarray(mask, dtype=bool))[0]
    if len(idx) == 0:
        return []
    parent = np.arange(len(idx))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(xy[idx]).query_pairs(1.3 * pitch_um):
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a
    groups = {}
    for i in range(len(idx)):
        groups.setdefault(find(i), []).append(int(idx[i]))
    return [g for g in groups.values() if len(g) >= min_size]


def read_visium_v1_positions(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8-sig") as handle:
            lines = [ln.rstrip("\n") for ln in handle if ln.strip()]
    except Exception as exc:
        return None, f"positions unreadable: {type(exc).__name__}"
    if not lines:
        return None, "positions file empty"
    header = [h.strip() for h in lines[0].split(",")]
    named = {"barcode", "in_tissue", "array_row", "array_col",
             "pxl_row_in_fullres", "pxl_col_in_fullres"}
    if set(header) >= named:
        idx = {h: i for i, h in enumerate(header)}
        out = []
        for ln in lines[1:]:
            parts = ln.split(",")
            try:
                out.append((parts[idx["barcode"]].strip(),
                            float(parts[idx["array_row"]]),
                            float(parts[idx["array_col"]])))
            except (ValueError, IndexError):
                continue
        if not out:
            return None, "positions header present but no parseable rows"
        return out, None
    if len(lines[0].split(",")) == 6:
        out = []
        for ln in lines:
            parts = ln.split(",")
            try:
                out.append((parts[0].strip(), float(parts[2]), float(parts[3])))
            except (ValueError, IndexError):
                continue
        if not out:
            return None, "positions 6-col layout but no parseable rows"
        return out, "positions header absent; used fixed 6-col v1 layout"
    return None, "positions layout unknown (need header or 6-col v1 layout)"


def parse_positions_text(text):
    if isinstance(text, (bytes, bytearray)):
        text = bytes(text).decode("utf-8", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None, "positions text empty"
    header = [h.strip() for h in lines[0].split(",")]
    named = {"barcode", "in_tissue", "array_row", "array_col",
             "pxl_row_in_fullres", "pxl_col_in_fullres"}
    if set(header) >= named:
        idx = {h: i for i, h in enumerate(header)}
        bc, rows, cols = [], [], []
        for ln in lines[1:]:
            parts = ln.split(",")
            try:
                bc.append(parts[idx["barcode"]].strip())
                rows.append(float(parts[idx["array_row"]]))
                cols.append(float(parts[idx["array_col"]]))
            except (ValueError, IndexError):
                continue
        if not bc:
            return None, "positions header present but no parseable rows"
        return bc, (rows, cols)
    if len(lines[0].split(",")) == 6:
        bc, rows, cols = [], [], []
        for ln in lines:
            parts = ln.split(",")
            try:
                bc.append(parts[0].strip())
                rows.append(float(parts[2]))
                cols.append(float(parts[3]))
            except (ValueError, IndexError):
                continue
        if not bc:
            return None, "positions 6-col layout but no parseable rows"
        return bc, (rows, cols)
    return None, "positions layout unknown (need header or 6-col v1 layout)"


def stream_tar_member(tar_path, predicate, byte_budget=8 * 1024 ** 3):
    try:
        handle = tarfile.open(tar_path, "r")
    except Exception:
        return None
    try:
        for member in handle:
            if member.isfile() and predicate(member.name):
                if member.size > byte_budget:
                    return None
                extracted = handle.extractfile(member)
                if extracted is None:
                    return None
                return extracted.read()
    except Exception:
        return None
    finally:
        try:
            handle.close()
        except Exception:
            pass
    return None


def read_10x_h5_matrix(path):
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(v) for v in matrix["shape"][()])
        values = matrix["data"][()]
        indices = matrix["indices"][()]
        indptr = matrix["indptr"][()]
        feats = matrix["features"]
        ids = decode(feats["id"][()] if "id" in feats else feats["name"][()])
        names = decode(feats["name"][()] if "name" in feats else feats["id"][()])
        ftypes = decode(feats["feature_type"][()]) if "feature_type" in feats else ["?"] * len(names)
        genomes = decode(feats["genome"][()]) if "genome" in feats else ["?"] * len(names)
        barcodes = decode(matrix["barcodes"][()])
    csc = sparse.csc_matrix((values.astype(float), indices, indptr), shape=shape)
    return csc.T.tocsr(), names, ids, ftypes, genomes, barcodes


def gene_index_from_names(names, ftypes, genomes):
    index = {}
    for i, (name, ftype, genome) in enumerate(zip(names, ftypes, genomes)):
        if ftype == "Gene Expression" and genome == "GRCh38":
            index.setdefault(name, []).append(i)
    return index


def gene_index_from_ensembl(names, ids):
    index = {}
    for i, (name, ens) in enumerate(zip(names, ids)):
        index.setdefault(name, []).append(i)
        index.setdefault(ens, []).append(i)
    return index


def load_cervilla_visium():
    p = ROOT / "data/other_sources/cervilla_2026/visium_Colorectal/filtered_feature_bc_matrix.h5"
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(p)
    pos, reason = read_visium_v1_positions(
        ROOT / "data/other_sources/cervilla_2026/visium_Colorectal/spatial/tissue_positions_list.csv")
    if pos is None:
        return None, reason
    keep = {b for b, _, _ in pos}
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut = {b: (r, c) for b, r, c in pos}
    rows = [lut[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 100)
    return ("cervilla-Visium-v1", X[order], names, xy, 100,
            {"vocab": "none", "y": None, "positive_values": None},
            {g: v for g, v in gene_index_from_names(names, ftypes, genomes).items()})


def load_cervilla_cytassist():
    p = ROOT / "data/other_sources/cervilla_2026/cytassist_Colorectal/filtered_feature_bc_matrix.h5"
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(p)
    pos, reason = read_visium_v1_positions(
        ROOT / "data/other_sources/cervilla_2026/cytassist_Colorectal/spatial/tissue_positions.csv")
    if pos is None:
        return None, reason
    keep = {b for b, _, _ in pos}
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut = {b: (r, c) for b, r, c in pos}
    rows = [lut[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 100)
    return ("cervilla-CytAssist-v2", X[order], names, xy, 100,
            {"vocab": "none", "y": None, "positive_values": None},
            {g: v for g, v in gene_index_from_names(names, ftypes, genomes).items()})


def load_parent_crc():
    p = ROOT / ("data/other_sources/10x_genomics/Parent_Visium_Human_ColorectalCancer/"
                "Parent_Visium_Human_ColorectalCancer_filtered_feature_bc_matrix.h5")
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(p)
    pos = stream_tar_member(
        ROOT / ("data/other_sources/10x_genomics/Parent_Visium_Human_ColorectalCancer/"
                "Parent_Visium_Human_ColorectalCancer_spatial.tar.gz"),
        lambda n: n.endswith("tissue_positions_list.csv"))
    if pos is None:
        return None, "parent spatial tar lacks tissue_positions_list.csv"
    bc, arr = parse_positions_text(pos)
    if bc is None:
        return None, arr
    keep = set(bc)
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut = {b: (r, c) for b, r, c in zip(bc, arr[0], arr[1])}
    rows = [lut[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 100)
    return ("parent-CRC-Visium", X[order], names, xy, 100,
            {"vocab": "none", "y": None, "positive_values": None},
            {g: v for g, v in gene_index_from_ensembl(names, ids).items()})


def load_cytassist11_crc():
    p = ROOT / ("data/other_sources/10x_genomics/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/"
                "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_filtered_feature_bc_matrix.h5")
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(p)
    pos = stream_tar_member(
        ROOT / ("data/other_sources/10x_genomics/CytAssist_11mm_FFPE_Human_Colorectal_Cancer/"
                "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_spatial.tar.gz"),
        lambda n: n.endswith("tissue_positions.csv"))
    if pos is None:
        return None, "cytassist11 spatial tar lacks tissue_positions.csv"
    bc, arr = parse_positions_text(pos)
    if bc is None:
        return None, arr
    keep = set(bc)
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut = {b: (r, c) for b, r, c in zip(bc, arr[0], arr[1])}
    rows = [lut[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 100)
    return ("cytassist-CRC-11mm", X[order], names, xy, 100,
            {"vocab": "none", "y": None, "positive_values": None},
            {g: v for g, v in gene_index_from_ensembl(names, ids).items()})


def load_usz(alias):
    path = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed" / f"{alias}.h5ad"
    with h5py.File(path, "r") as handle:
        X = handle["X"][:]
        genes = decode(handle["var"]["gene_name"][:])
        obs = handle["obs"]
        cats = decode(obs["ground_truth"]["categories"][:])
        codes = np.asarray(obs["ground_truth"]["codes"][:])
        lab = np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
        xa = np.asarray(obs["x_array"][:], float)
        ya = np.asarray(obs["y_array"][:], float)
    xy = hex_xy(xa, ya, 100)
    index = {}
    for i, g in enumerate(genes):
        index.setdefault(g, []).append(i)
    return (f"usz-{alias}", sparse.csr_matrix(X), genes, xy, 100,
            {"vocab": "ground_truth TLS vs NOR", "y": lab,
             "positive_values": ["TLS"]}, index)


def load_gse175540(gsm):
    base = ROOT / "data/GEO/GSE175540/raw"
    matrix = base / f"{gsm}_filtered_feature_bc_matrix.h5"
    posp = base / f"{gsm}_tissue_positions_list.csv.gz"
    annp = base / f"{gsm}_TLS_annotation.csv.gz"
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(matrix)
    pos, reason = read_visium_v1_positions(posp)
    if pos is None:
        return None, reason
    keep = {b for b, _, _ in pos}
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut = {b: (r, c) for b, r, c in pos}
    rows = [lut[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 100)
    lab = np.full(len(order), "", dtype=object)
    try:
        import csv as _csv
        with gzip.open(annp, "rt", encoding="utf-8-sig") as handle:
            amap = {}
            for row in _csv.DictReader(handle):
                amap[row.get("Barcode", "").strip()] = row.get("TLS_2_cat", "").strip()
    except FileNotFoundError:
        return None, f"annotation file missing: {annp.name}"
    ordered = [barcodes[i] for i in order]
    for j, b in enumerate(ordered):
        if b in amap:
            lab[j] = amap[b]
    return (f"gse175540-{gsm}", X[order], names, xy, 100,
            {"vocab": "TLS_2_cat: TLS vs NO_TLS", "y": lab,
             "positive_values": ["TLS"]},
            {g: v for g, v in gene_index_from_names(names, ftypes, genomes).items()})


def load_stcrc(sample):
    base = ROOT / "data/derived/r04/st_crc_cms" / sample
    matrix = base / "filtered_feature_bc_matrix.h5"
    if not matrix.exists():
        return None, f"derived matrix missing: {sample} (no zip streaming for STCRC)"
    X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(matrix)
    posp = ROOT / "data/other_sources/zenodo_st_crc_cms/tissue_positions" / f"{sample}_tissue_positions_list.csv"
    pos, reason = read_visium_v1_positions(posp)
    if pos is None:
        return None, reason
    ann = ROOT / "data/other_sources/zenodo_st_crc_cms/Pathology_SpotAnnotations" / f"Pathologist_Annotations_{sample}.csv"
    if not ann.exists():
        return None, f"pathology annotation missing: {ann.name}"
    import csv as _csv
    with open(ann, encoding="utf-8-sig") as handle:
        rows = list(_csv.DictReader(handle))
    labcol = [c for c in rows[0].keys() if c != "Barcode"][0]
    lut = {r["Barcode"].strip(): r[labcol].strip() for r in rows}
    keep = {b for b, _, _ in pos}
    order = [i for i, b in enumerate(barcodes) if b in keep]
    lut2 = {b: (r, c) for b, r, c in pos}
    rows2 = [lut2[b] for b in (barcodes[i] for i in order)]
    xy = hex_xy([r for r, _ in rows2], [c for _, c in rows2], 100)
    lab = np.array([lut.get(barcodes[i], "") for i in order], dtype=object)
    return (f"stcrc-{sample}", X[order], names, xy, 100,
            {"vocab": f"pathologist {labcol}; IC aggregate_* FAIL-CLOSED",
             "y": lab, "positive_values": []},
            {g: v for g, v in gene_index_from_names(names, ftypes, genomes).items()})


def load_htan(capture):
    path = ROOT / "data/other_sources/htan/HTAN_Vanderbilt_CRC_Visium_OSF_hftq2" / f"{capture}_filtered_trimmed.h5ad"
    with h5py.File(path, "r") as handle:
        barcodes = decode(handle["obs"]["_index"][:])
        genes = decode(handle["var"]["_index"][:])
        matrix = handle["X"]
        if isinstance(matrix, h5py.Group):
            X = sparse.csr_matrix((matrix["data"][()], matrix["indices"][()],
                                   matrix["indptr"][()]),
                                  shape=(len(barcodes), len(genes)))
        else:
            X = sparse.csr_matrix(np.asarray(matrix[()], dtype=float))
        xa = np.asarray(handle["obs"]["array_row"][:], float)
        ya = np.asarray(handle["obs"]["array_col"][:], float)
    xy = hex_xy(xa, ya, 100)
    index = {}
    for i, g in enumerate(genes):
        index.setdefault(g, []).append(i)
    return (f"htan-{capture}", X, genes, xy, 100,
            {"vocab": "Heiser pathology CSV joined on barcode (fail-closed)",
             "y": np.full(len(barcodes), "", dtype=object),
             "positive_values": ["lymphoid_follicle"],
             "__barcodes__": barcodes}, index)


def attach_htan_labels(loader_out, key_short):
    section_id, X, genes, xy, pitch, labels, index = loader_out
    barcodes = labels.pop("__barcodes__")
    import csv as _csv
    path = ROOT / "data/other_sources/htan/spatial_CRC_atlas_repo/resources/ST" / f"{key_short}_pathology_annotation.csv"
    seen, dup = {}, []
    with open(path, encoding="utf-8-sig") as handle:
        for row in _csv.DictReader(handle):
            b = row["Barcode"].strip()
            v = row["pathology_annotation"].strip()
            if b in seen:
                dup.append(b)
            seen[b] = v
    y = np.array([seen.get(b, "") for b in barcodes], dtype=object)
    labels = dict(labels)
    labels["y"] = y
    labels["vocab"] += f"; duplicates_in_csv={len(dup)}"
    return section_id, X, genes, xy, pitch, labels, index


def load_visiumhd_6p5():
    p = ROOT / ("data/other_sources/10x_genomics/Visium_HD_6p5mm_Human_Colon_Cancer/"
                "Visium_HD_6p5mm_Human_Colon_Cancer_binned_outputs.tar.gz")
    member = "binned_outputs/square_016um/filtered_feature_bc_matrix.h5"
    raw = stream_tar_member(p, lambda n: n == member, byte_budget=8 * 1024 ** 3)
    if raw is None:
        return None, "HD 16um h5 member unreadable or over budget"
    mf = Path("/tmp/tls_pool_hd16.h5")
    mf.write_bytes(raw)
    try:
        X, names, ids, ftypes, genomes, barcodes = read_10x_h5_matrix(mf)
        pos = stream_tar_member(
            p, lambda n: n == "binned_outputs/square_016um/spatial/tissue_positions.parquet",
            byte_budget=8 * 1024 ** 3)
        if pos is None:
            return None, "HD 16um positions member unreadable or over budget"
        import pyarrow.parquet as pq
        table = pq.read_table(io.BytesIO(pos))
        cols = {n: table.column(n).to_pylist() for n in table.column_names}
        lut = {b: (r, c) for b, r, c in zip(cols["barcode"], cols["array_row"], cols["array_col"])}
        keep = [b for b in barcodes if b in lut]
        order = [i for i, b in enumerate(barcodes) if b in keep]
        rows = [lut[b] for b in (barcodes[i] for i in order)]
        xy = hex_xy([r for r, _ in rows], [c for _, c in rows], 16)
        return ("hd-colon-16um", X[order], names, xy, 55,
                {"vocab": "none", "y": None, "positive_values": None},
                {g: v for g, v in gene_index_from_ensembl(names, ids).items()})
    finally:
        try:
            mf.unlink(missing_ok=True)
        except Exception:
            pass


def load_xenium(alias):
    zp = ROOT / f"data/other_sources/10x_genomics/{alias}/{alias}_outs.zip"
    try:
        archive = zipfile.ZipFile(zp)
    except Exception as exc:
        return None, f"xenium zip unreadable: {type(exc).__name__}"
    try:
        panel = json.loads(archive.read("gene_panel.json"))
        targets = panel["payload"]["targets"]
        genes = [t["type"]["data"]["name"] for t in targets]
    except Exception as exc:
        return None, f"xenium panel unreadable: {type(exc).__name__}"
    try:
        raw = archive.read("cell_feature_matrix.h5")
    except KeyError:
        return None, "xenium cell_feature_matrix.h5 missing"
    mf = Path(f"/tmp/tls_pool_xenium_{alias.split('_')[2].lower()}.h5")
    mf.write_bytes(raw)
    try:
        with h5py.File(mf, "r") as handle:
            m = handle["matrix"]
            shape = tuple(int(v) for v in m["shape"][()])
            csc = sparse.csc_matrix((m["data"][()], m["indices"][()], m["indptr"][()]), shape=shape)
            X = csc.T.tocsr()
            barcodes = decode(m["barcodes"][()])
        import pyarrow.parquet as pq
        cells = pq.read_table(io.BytesIO(archive.read("cells.parquet")),
                              columns=["cell_id", "x_centroid", "y_centroid"]).to_pandas()
        lut = {str(r.cell_id): (float(r.x_centroid), float(r.y_centroid)) for r in cells.itertuples()}
        order = [i for i, b in enumerate(barcodes) if b in lut]
        if not order:
            return None, "xenium cell_id join produced zero rows"
        xy = np.array([lut[barcodes[i]] for i in order])
        gx = np.asarray(np.floor(np.asarray(xy[:, 0]) / 55.0), dtype=int)
        gy = np.asarray(np.floor(np.asarray(xy[:, 1]) / 55.0), dtype=int)
        grid = {}
        for k, (xx, yy) in enumerate(zip(gx, gy)):
            grid.setdefault((int(xx), int(yy)), []).append(k)
        keep_cells, bx, by = [], [], []
        for (ggx, ggy), members in grid.items():
            keep_cells.append(members)
            bx.append((ggx + 0.5) * 55)
            by.append((ggy + 0.5) * 55)
        agg = sparse.lil_matrix((len(keep_cells), X.shape[1]))
        for bi, members in enumerate(keep_cells):
            agg[bi] = X[[order[m] for m in members]].sum(axis=0)
        Xb = agg.tocsr()
        # snap each 55um bin to the nearest hex lattice node so components() sees hex adjacency
        hx = np.round(np.asarray(bx) / 55.0 - 0.5).astype(int)
        hy = np.round(np.asarray(by) / (55.0 * (3 ** 0.5) / 2)).astype(int)
        xyb = hex_xy(hy, hx, 55)
        index = {}
        for i, g in enumerate(genes):
            index.setdefault(g, []).append(i)
        return (f"xenium-{alias.split('_')[2].lower()}", Xb, genes, xyb, 55,
                {"vocab": "none", "y": None, "positive_values": None}, index)
    finally:
        try:
            mf.unlink(missing_ok=True)
        except Exception:
            pass


def gene_index_from_ensembl(names, ids):
    index = {}
    for i, (name, ens) in enumerate(zip(names, ids)):
        index.setdefault(name, []).append(i)
        index.setdefault(ens, []).append(i)
    return index


def run_section(cohort, section_id, X, genes, xy, pitch_um, labels, gene_index,
                diag_rows, comp_rows, bag):
    sig = signature_symbols()
    cols, missing = [], []
    for g in sig:
        if g in gene_index:
            cols.extend(gene_index[g] if isinstance(gene_index[g], list) else [gene_index[g]])
        else:
            missing.append(g)
    coverage = (len(sig) - len(missing), len(sig))
    diag = {"cohort": cohort, "section_id": section_id,
            "n_spots": int(X.shape[0]), "n_genes": int(X.shape[1]),
            "sig_coverage": f"{coverage[0]}/{coverage[1]}",
            "sig_missing": ";".join(sorted(missing)),
            "pitch_um": pitch_um, "status": "ok", "reason": "",
            "mask_frac": "", "n_components": 0, "n_comp_spots": 0,
            "auc_vs_labels": "", "null_median": "", "obs_over_null": "",
            "label_vocab": labels.get("vocab", "") if labels else "",
            "extent_um": ""}
    try:
        if X.shape[0] != len(xy):
            diag.update(status="skipped",
                        reason=f"matrix/coordinate row mismatch {X.shape[0]} vs {len(xy)}")
            diag_rows.append(diag)
            return
        if X.shape[0] < 10 or not cols:
            diag.update(status="skipped",
                        reason="too few spots or zero signature genes present")
            diag_rows.append(diag)
            return
        scores = aucell(X, cols, 0.05)
        finite = np.isfinite(scores)
        if finite.sum() < 10:
            diag.update(status="skipped", reason="fewer than 10 finite scores")
            diag_rows.append(diag)
            return
        thr = float(np.quantile(scores[finite], 0.9))
        mask = finite & (scores >= thr)
        extent = (float(np.ptp(xy[:, 0])) if len(xy) else 0.0,
                  float(np.ptp(xy[:, 1])) if len(xy) else 0.0)
        diag["extent_um"] = f"{extent[0]:.0f}x{extent[1]:.0f}"
        comps = components(xy, mask, pitch_um)
        diag["mask_frac"] = round(float(mask.mean()), 4)
        diag["n_components"] = len(comps)
        diag["n_comp_spots"] = int(sum(len(c) for c in comps))
        rng = np.random.default_rng(abs(hash(section_id)) % (2 ** 32))
        null_counts = []
        for _ in range(10):
            pmask = finite & (scores[rng.permutation(len(scores))] >= thr)
            null_counts.append(len(components(xy, pmask, pitch_um)))
        nm = float(np.median(null_counts))
        diag["null_median"] = nm
        diag["obs_over_null"] = round(len(comps) / nm, 3) if nm > 0 else "inf_if_null_zero"
        lab = labels.get("y") if labels else None
        pos = labels.get("positive_values") if labels else None
        if lab is not None and pos:
            from scipy.stats import rankdata as _rd
            y = np.asarray([1 if v in pos else 0 for v in lab])
            if y.sum() > 0 and (1 - y).sum() > 0 and finite.sum() > 0:
                r = _rd(-scores[finite])
                order = np.argsort(r, kind="stable")
                yy = y[finite][order]
                tp = np.cumsum(yy == 1)
                fp = np.cumsum(yy == 0)
                tpr = tp / max(yy.sum(), 1)
                with np.errstate(invalid="ignore"):
                    xs = fp / max((1 - yy).sum(), 1)
                    auc = float(np.sum((tpr[1:] + tpr[:-1]) / 2 * np.diff(xs)))
                diag["auc_vs_labels"] = round(auc, 3)
        bag["expected_components"] = bag.get("expected_components", 0) + len(comps)
        for ci, c in enumerate(comps):
            cx, cy = float(xy[c, 0].mean()), float(xy[c, 1].mean())
            comp_rows.append({
                "cohort": cohort, "section_id": section_id,
                "component_id": f"{section_id}::c{ci}",
                "cx_um": round(cx, 1), "cy_um": round(cy, 1),
                "n_spots": len(c), "mean_score": round(float(scores[c].mean()), 4),
            })
        diag_rows.append(diag)
    except Exception as exc:
        diag.update(status="skipped", reason=f"{type(exc).__name__}: {str(exc)[:160]}")
        diag_rows.append(diag)


def smoke_loaders():
    for name, fn in [("usz-LC3", lambda: load_usz("LC3")),
                     ("gse175540-ffpe", lambda: load_gse175540("GSM5924030_ffpe_c_2")),
                     ("stcrc-Rep1", lambda: load_stcrc("SN048_A121573_Rep1")),
                     ("htan-6723_KL_1", lambda: attach_htan_labels(load_htan("6723_KL_1"), "6723_1")),
                     ("cervilla-visium", load_cervilla_visium),
                     ("hd-6p5-16um", load_visiumhd_6p5),
                     ("xenium-kidney", lambda: load_xenium("Xenium_V1_hKidney_cancer_section"))]:
        try:
            out = fn()
            if out is None or (isinstance(out, tuple) and out[0] is None):
                print(f"{name}: SKIP", out[1] if isinstance(out, tuple) else "")
            else:
                sid, X, genes, xy, pitch, labels, index = out
                print(f"{name}: ok spots={X.shape[0]} genes={X.shape[1]} pitch={pitch} "
                      f"extent={np.ptp(xy[:, 0]):.0f}x{np.ptp(xy[:, 1]):.0f}")
        except Exception as exc:
            print(f"{name}: ERROR {type(exc).__name__}: {str(exc)[:200]}")


def run_pool():
    run_c = contract()
    diag_rows, comp_rows, bag = [], [], {}
    section_ids = []

    def emit(out_tuple, cohort):
        if out_tuple is None or (isinstance(out_tuple, tuple) and out_tuple[0] is None):
            section_ids.append((cohort, "<unresolved>", "skipped: loader unresolved"))
            return
        sid, X, genes, xy, pitch, labels, index = out_tuple
        section_ids.append((cohort, sid, ""))
        run_section(cohort, sid, X, genes, xy, pitch, labels, index, diag_rows, comp_rows, bag)

    for alias in ["KC1", "KC2", "KC3", "LC1", "LC2", "LC3", "LC4", "LC5"]:
        emit(load_usz(alias), "USZ")
    qual = json.loads((ROOT / "infra/gobp_external_qualification_20260925/run_contract.json").read_text())
    gsm_done = set()
    for src in qual["selected_sources"]:
        gsm = src["gsm"]
        if gsm in gsm_done:
            continue
        gsm_done.add(gsm)
        stem = Path(src["matrix_path"]).stem.replace("_filtered_feature_bc_matrix", "")
        emit(load_gse175540(stem), "GSE175540")
    assert len(gsm_done) == 18, f"expected 18 GSMs, got {len(gsm_done)}"
    manifest = json.loads((ROOT / "infra/r04/role_manifests/internal_validation_manifest.json").read_text())
    assert len(manifest["rows"]) == 14
    for row in manifest["rows"]:
        emit(load_stcrc(row["section_id"].split("::")[-1]), "STCRC")
    for cap, key in [("6723_KL_4", "6723_4"), ("8578_AS_4", "8578_4"),
                     ("7319_AS_2", "7319_2"), ("8899_AS_5", "8899_5"),
                     ("8899_AS_6", "8899_6"), ("8270_AS_2", "8270_2"),
                     ("7003_AS_5", "7003_5"), ("8899_AS_8", "8899_8")]:
        emit(attach_htan_labels(load_htan(cap), key), "HTAN")
    emit(load_cervilla_visium(), "Cervilla")
    emit(load_cervilla_cytassist(), "Cervilla")
    emit(load_parent_crc(), "Parent")
    emit(load_cytassist11_crc(), "Parent")
    emit(load_visiumhd_6p5(), "VisiumHD")
    emit(load_xenium("Xenium_V1_hKidney_cancer_section"), "Xenium")
    emit(load_xenium("Xenium_V1_hLung_cancer_section"), "Xenium")

    with (OUT / "pool_diagnostics.tsv").open("w", newline="") as handle:
        cols = ["cohort", "section_id", "n_spots", "n_genes", "sig_coverage",
                "sig_missing", "pitch_um", "status", "reason", "mask_frac",
                "n_components", "n_comp_spots", "auc_vs_labels", "null_median",
                "obs_over_null", "label_vocab", "extent_um"]
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(diag_rows)
    with (OUT / "pool_components.tsv").open("w", newline="") as handle:
        cols = ["cohort", "section_id", "component_id", "cx_um", "cy_um",
                "n_spots", "mean_score"]
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(comp_rows)
    import datetime
    summary = {
        "sections_attempted": len(section_ids),
        "sections_ok": sum(1 for d in diag_rows if d["status"] == "ok"),
        "sections_skipped": sum(1 for d in diag_rows if d["status"] != "ok"),
        "expected_components": bag.get("expected_components", 0),
        "finished_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    (OUT / "pool_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def expected_section_ids():
    run_c = contract()
    ids = []
    for s in run_c["sections"]:
        sid = s["id"]
        if sid == "gse175540:18-R02-frozen":
            ids.append("gse175540-<18 GSMs from run_contract.json>")
        elif sid == "stcrc:14-manifest-sections":
            ids.append("stcrc-<14 samples from internal_validation_manifest.json>")
        elif sid == "htan:8-lymphoid-verifier-captures":
            ids.append("htan-<8 captures from lymphoid table>")
        elif sid == "usz:8-rerun-under-v3":
            ids.append("usz-<KC1..LC5>")
        else:
            ids.append(sid)
    return ids


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="Load the first section of the first four cohorts and print loader summaries only.")
    ap.add_argument("--run", action="store_true", help="Run all contract sections.")
    args = ap.parse_args()
    run_c = contract()
    assert run_c["status"] == "FROZEN_BEFORE_SCORING"
    if args.smoke:
        smoke_loaders()
    elif args.run:
        run_pool()
    else:
        raise SystemExit("pass --smoke or --run")


if __name__ == "__main__":
    main()
