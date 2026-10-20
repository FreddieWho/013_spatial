#!/usr/bin/env python3
"""A fixed, CPU-only molecular recovery experiment. No human TLS inputs.

Run through run_field_discovery_20261018.sh to use the existing clean CPU environment.
The input guard covers repository Python opens; HDF5 reads additionally pass an
explicit file/field allowlist. This is an auditable launcher, not an OS sandbox.
"""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import os
import pickle
import shutil
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import h5py
import numpy as np
import scipy
import sklearn
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/field_discovery_20261018"
HEST = ROOT / "data/other_sources/hest1k"
ORGANS = ("Bowel", "Breast", "Lung", "Liver")
PLATFORMS = ("Visium", "Visium HD")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def tsv_write(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        path.write_text(""); return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def finite(value):
    return float(value) if np.isfinite(value) else None


def reserve(config):
    free = shutil.disk_usage(WORK).free
    if free < config["min_free_bytes"] + config["estimated_max_new_bytes"]:
        raise RuntimeError(f"Storage reserve violated: free={free}")


class InputGuard:
    def __init__(self, root, allowed, output):
        self.root = Path(root).resolve()
        self.allowed = {Path(p).resolve() for p in allowed}
        self.output = Path(output).resolve()
        self.events = Counter()
        self.records = {}

    def check(self, path):
        path = Path(path).resolve()
        if path not in self.allowed:
            raise PermissionError(f"Input allowlist rejected {path}")
        return path

    def record(self, path, fields, purpose):
        path = self.check(path)
        key = str(path.relative_to(self.root))
        if key not in self.records:
            self.records[key] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
        self.records[key].update(fields=fields, purpose=purpose)
        return path

    def install(self):
        def audit(event, args):
            if event != "open" or not isinstance(args[0], (str, bytes, os.PathLike)):
                return
            path = Path(os.fsdecode(args[0])).resolve()
            if not path.is_relative_to(self.root):
                return  # Python libraries and OS runtime resources, not project data.
            if path.is_relative_to(self.output):
                return
            mode = args[1]
            flags = args[2] if len(args) > 2 else 0
            writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or (
                isinstance(flags, int) and bool(flags & (os.O_WRONLY | os.O_RDWR)))
            if writing:
                raise PermissionError(f"Write outside experiment output: {path}")
            self.check(path)
            self.events[str(path.relative_to(self.root))] += 1
        sys.addaudithook(audit)

    def receipt(self):
        return {"inputs": self.records, "python_read_events": dict(self.events),
                "human_tls_files_opened": [],
                "scope": "repository Python-open guard plus explicit HDF5 path/field whitelist; not an OS sandbox"}


def load_config():
    config = json.loads((WORK / "contract.json").read_text())
    split = ROOT / config["split_file"]
    if sha256(split) != config["split_sha256"]:
        raise ValueError("Frozen split hash changed")
    rows = list(csv.DictReader(split.open(), delimiter="\t"))
    for row in rows:
        row["fold"] = int(row["fold"])
    assert_grouped(rows)
    allowed = {WORK / "contract.json", split, ROOT / config["panel_file"], Path(__file__)}
    for row in rows:
        sid = row["section_id"]
        allowed.update({HEST / f"st/{sid}.h5ad", HEST / f"metadata/{sid}.json"})
    guard = InputGuard(ROOT, allowed, WORK)
    guard.install()
    for path in [WORK / "contract.json", split, ROOT / config["panel_file"], Path(__file__)]:
        guard.record(path, ["whole file"], "frozen experiment definition")
    panel = json.loads((ROOT / config["panel_file"]).read_text())
    return config, rows, panel, guard


def assert_grouped(rows):
    groups = {}
    for row in rows:
        unit, fold = row["unit_id"], int(row["fold"])
        if unit in groups and groups[unit] != fold:
            raise ValueError(f"Unit crosses folds: {unit}")
        groups[unit] = fold


def strings(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def read_matrix(node, nobs, nvar):
    if isinstance(node, h5py.Dataset):
        x = sparse.csr_matrix(node[:])
    else:
        encoding = node.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        shape = tuple(node.attrs.get("shape", (nobs, nvar)))
        parts = (node["data"][:], node["indices"][:], node["indptr"][:])
        if encoding == "csr_matrix":
            x = sparse.csr_matrix(parts, shape=shape)
        elif encoding == "csc_matrix":
            x = sparse.csc_matrix(parts, shape=shape).tocsr()
        else:
            nind = len(parts[2]) - 1
            if nobs == nvar:
                raise ValueError("Ambiguous square sparse matrix without encoding-type")
            if nind == nobs:
                x = sparse.csr_matrix(parts, shape=(nobs, nvar))
            elif nind == nvar:
                x = sparse.csc_matrix(parts, shape=(nobs, nvar)).tocsr()
            else:
                raise ValueError("Sparse shape/encoding mismatch")
    if x.shape != (nobs, nvar):
        if x.shape == (nvar, nobs):
            x = x.T.tocsr()
        else:
            raise ValueError(f"Expression shape mismatch: {x.shape}, {(nobs, nvar)}")
    if not np.all(np.isfinite(x.data)) or np.any(x.data < 0):
        raise ValueError("Expression must be finite nonnegative source counts")
    return x


def normalized_subset(x, columns, library):
    sub = x[:, columns].astype(np.float64).tocsr()
    scale = 10000.0 / np.maximum(library, 1.0)
    sub.data = np.log1p(sub.data * np.repeat(scale, np.diff(sub.indptr)))
    return sub


def target_readout(x, genes, members):
    mask = np.isin(genes, list(members))
    if not np.any(mask):
        raise ValueError("No target genes")
    library = np.asarray(x.sum(axis=1)).ravel()
    return np.asarray(normalized_subset(x, mask, library).mean(axis=1)).ravel().astype(np.float32)


def visible_features(x, genes, panel, target, names):
    hidden = set(panel[target])
    definitions = [set(panel[name]) - hidden for name in names]
    union = set().union(*definitions)
    if hidden & union:
        raise AssertionError("Target genes entered visible input")
    mask = np.isin(genes, list(union))
    visible_genes = genes[mask]
    counts = x[:, mask].tocsr()
    library = np.asarray(counts.sum(axis=1)).ravel()
    z = normalized_subset(x, mask, library)
    columns = []
    for members in definitions:
        select = np.isin(visible_genes, list(members))
        if not np.any(select):
            raise ValueError("No genes remain in input module")
        columns.append(np.asarray(z[:, select].mean(axis=1)).ravel())
    columns += [np.log1p(library), np.asarray((counts > 0).sum(axis=1)).ravel() / max(len(visible_genes), 1)]
    return np.column_stack(columns).astype(np.float32)


def neighbor_graphs(xy):
    tree = cKDTree(xy)
    distances = tree.query(xy, k=2)[0][:, 1]
    if np.any(distances <= 0) or not np.all(np.isfinite(distances)):
        raise ValueError("Duplicate/nonfinite spatial coordinates")
    pitch = float(np.median(distances))
    pairs = tree.query_pairs(3.1 * pitch, output_type="ndarray")
    dd = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1) / pitch
    graphs = []
    for select in [dd <= 1.1, dd > 1.1]:
        edges = pairs[select]
        graph = sparse.csr_matrix((np.ones(2 * len(edges), dtype=np.float32),
                                  (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])),
                                 shape=(len(xy), len(xy)))
        degrees = np.asarray(graph.sum(axis=1)).ravel()
        graph = sparse.diags(1 / np.maximum(degrees, 1)) @ graph
        graphs.append(graph.tocsr())
    return *graphs, pitch


def neighbor_mean(graph, values):
    mean = graph @ values
    missing = np.diff(graph.indptr) == 0
    mean[missing] = values[missing]
    return mean


def input_names(config, target, view):
    return config["small_inputs"] if view == "small5" else [n for n in config["dictionary"] if n != target]


def metadata_value(value):
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return None
    if isinstance(value, str) and value.strip().lower() in ("", "nan", "none", "unknown"):
        return None
    return value


def prepare():
    config, rows, panel, guard = load_config()
    reserve(config)
    if (WORK / "input_manifest.json").exists():
        raise RuntimeError("Prepared inputs exist; do not overwrite a run")
    (WORK / "cache").mkdir(exist_ok=True)
    manifest = []; definitions = {}
    # Qualification reads presence/coordinates, never target strength or human labels.
    for row in rows:
        sid = row["section_id"]; row = dict(row)
        mp = guard.record(HEST / f"metadata/{sid}.json", config["metadata_fields"], "non-target sample context")
        raw_meta = json.loads(mp.read_text())
        meta = {k: metadata_value(raw_meta.get(k)) for k in config["metadata_fields"]}
        row.update(metadata=meta, patient_known=meta["patient"] is not None,
                   source_key=meta["study_link"] or "title:" + str(meta["dataset_title"]))
        path = guard.record(HEST / f"st/{sid}.h5ad", config["allowed_h5_fields"], "molecular measurements and geometry")
        with h5py.File(path) as f:
            genes = strings(f["var/_index"][:]); found = set(genes)
            keep = np.flatnonzero(f["obs/in_tissue"][:] == 1)
            xy = np.column_stack((f["obs/pxl_col_in_fullres"][:][keep], f["obs/pxl_row_in_fullres"][:][keep]))
        failures = []; coverage = {}
        for target in config["targets"]:
            td = {}; definitions[target] = {}
            for view in config["input_sizes"]:
                names = input_names(config, target, view)
                visible = {name: sorted(set(panel[name]) - set(panel[target])) for name in names}
                definitions[target][view] = visible
                td[view] = {}
                for name, members in {target: panel[target], **visible}.items():
                    n = len(found & set(members)); cov = n / len(members)
                    td[view][name] = {"n": n, "coverage": cov}
                    q = config["qualification"]
                    if cov < q["min_gene_coverage"] or n < q["min_genes_per_module"]:
                        failures.append(f"{target}/{view}/{name}: {n}/{len(members)} genes")
            coverage[target] = td
        if len(keep) < config["qualification"]["min_tissue_spots"]:
            failures.append("TOO_FEW_TISSUE_SPOTS")
        if not np.all(np.isfinite(xy)):
            failures.append("NONFINITE_COORDINATES")
        if meta["organ"] != row["organ"]:
            failures.append("ORGAN_METADATA_CONFLICT")
        if meta["st_technology"] not in PLATFORMS:
            failures.append("UNSUPPORTED_PLATFORM")
        row.update(status="EXCLUDED" if failures else "QUALIFIED", reasons=sorted(set(failures)),
                   n_spots=len(keep), n_genes=len(genes), coverage=coverage,
                   expression_sha256=guard.records[str(path.relative_to(ROOT))]["sha256"])
        manifest.append(row)
    assert_grouped(manifest)
    json_write(WORK / "input_manifest.json", manifest)
    json_write(WORK / "visible_definitions.json", definitions)
    json_write(WORK / "environment.json", {"python": sys.version, "executable": sys.executable,
               "numpy": np.__version__, "scipy": scipy.__version__, "h5py": h5py.__version__,
               "sklearn": sklearn.__version__, "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
               "contract_sha256": sha256(WORK / "contract.json"), "code_sha256": sha256(__file__)})
    # All qualification choices are persisted before reading expression values.
    for row in manifest:
        if row["status"] != "QUALIFIED":
            continue
        reserve(config)
        sid = row["section_id"]; path = guard.check(HEST / f"st/{sid}.h5ad")
        with h5py.File(path) as f:
            genes = strings(f["var/_index"][:]); barcodes = strings(f["obs/_index"][:])
            keep = np.flatnonzero(f["obs/in_tissue"][:] == 1)
            xy = np.column_stack((f["obs/pxl_col_in_fullres"][:][keep], f["obs/pxl_row_in_fullres"][:][keep]))
            x = read_matrix(f["X"], len(barcodes), len(genes))[keep]
        if len(set(barcodes[keep])) != len(keep):
            raise ValueError(f"Duplicate source barcodes in {sid}")
        near, outer, pitch = neighbor_graphs(xy)
        arrays = {"barcodes": barcodes[keep], "source_row": keep, "xy_source_pixels": xy,
                  "y": np.column_stack([target_readout(x, genes, panel[t]) for t in config["targets"]])}
        for ti, target in enumerate(config["targets"]):
            for view in config["input_sizes"]:
                arrays[f"x{ti}_{view}"] = visible_features(x, genes, panel, target, input_names(config, target, view))
        np.savez_compressed(WORK / f"cache/{sid}.npz", **arrays)
        sparse.save_npz(WORK / f"cache/{sid}_near.npz", near)
        sparse.save_npz(WORK / f"cache/{sid}_outer.npz", outer)
        row.update(empirical_pitch_source_pixels=pitch,
                   estimated_pitch_um=finite(pitch * (row["metadata"]["pixel_size_um_estimated"] or np.nan)),
                   near_isolated_spots=int(np.sum(np.diff(near.indptr) == 0)))
        print(json.dumps({"prepared": sid, "n_spots": len(keep)}, ensure_ascii=False), flush=True)
    json_write(WORK / "input_manifest.json", manifest)
    json_write(WORK / "input_receipt_prepare.json", guard.receipt())
    (WORK / "STATUS").write_text("INPUTS_PREPARED_NO_MODELS_FIT\n")
    print("QUALIFIED", Counter(r["organ"] for r in manifest if r["status"] == "QUALIFIED"), flush=True)


def sample_indices(n, sid, seed, limit):
    hashed = hashlib.sha256(f"{seed}|{sid}".encode()).digest()
    rng = np.random.default_rng(int.from_bytes(hashed[:8], "little"))
    return np.sort(rng.choice(n, min(n, limit), replace=False))


def load_cache(row):
    sid = row["section_id"]
    with np.load(WORK / f"cache/{sid}.npz", allow_pickle=False) as z:
        cache = {key: z[key] for key in z.files}
    cache["near"] = sparse.load_npz(WORK / f"cache/{sid}_near.npz")
    cache["outer"] = sparse.load_npz(WORK / f"cache/{sid}_outer.npz")
    return cache


def feature_views(row, cache, ti, view):
    x = cache[f"x{ti}_{view}"]
    meta = np.array([float(row["organ"] == o) for o in ORGANS] +
                    [float(row["metadata"]["st_technology"] == p) for p in PLATFORMS], dtype=np.float32)
    broadcast = lambda a: np.broadcast_to(a, (len(x), len(a)))
    point = np.c_[x, broadcast(meta)]
    bag = np.c_[point, broadcast(x.mean(axis=0)), broadcast(x.std(axis=0))]
    local = np.c_[neighbor_mean(cache["near"], x), neighbor_mean(cache["outer"], x),
                  np.log1p(np.diff(cache["near"].indptr)), np.log1p(np.diff(cache["outer"].indptr))]
    return {"point": point.astype(np.float32), "bag": bag.astype(np.float32),
            "space": np.c_[bag, local].astype(np.float32)}


def training_data(rows, caches, features, ti, seed, config):
    sections_per_unit = Counter(r["unit_id"] for r in rows)
    indices = {r["section_id"]: sample_indices(r["n_spots"], r["section_id"], seed,
               config["training"]["max_spots_per_section"]) for r in rows}
    y = []; weights = []; xs = defaultdict(list); section_means = {}
    for row in rows:
        sid = row["section_id"]; idx = indices[sid]
        yy = caches[sid]["y"][idx, ti]; y.append(yy)
        weights.append(np.full(len(idx), 1 / (len(idx) * sections_per_unit[row["unit_id"]])))
        section_means[sid] = float(yy.mean())
        for level, array in features[sid].items():
            xs[level].append(array[idx])
    weight = np.concatenate(weights); weight /= weight.mean()
    return {k: np.concatenate(v) for k, v in xs.items()}, np.concatenate(y), weight, indices, section_means


def priors(rows, section_means):
    nu = Counter(r["unit_id"] for r in rows)
    accum = defaultdict(list)
    for row in rows:
        value = section_means[row["section_id"]]; w = 1 / nu[row["unit_id"]]
        for key in ["global", "organ:" + row["organ"], "source:" + row["source_key"]]:
            accum[key].append((value, w))
    return {k: float(np.average([v for v, w in pairs], weights=[w for v, w in pairs])) for k, pairs in accum.items()}


def correlation(a, b):
    a = a - a.mean(); b = b - b.mean()
    den = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return finite(np.dot(a, b) / den) if den > 1e-12 else None


def metrics(y, pred):
    return {"mse": float(np.mean((y - pred) ** 2)),
            "centered_mse": float(np.mean(((y - y.mean()) - (pred - pred.mean())) ** 2)),
            "pearson": correlation(y, pred)}


def stability(a, b):
    k = max(1, math.ceil(len(a) * 0.1))
    top = lambda x: set(np.lexsort((np.arange(len(x)), -x))[:k])
    aa, bb = top(a), top(b)
    return {"repeat_pearson": correlation(a, b), "repeat_top_decile_jaccard": len(aa & bb) / len(aa | bb)}


def run_fold(fold):
    config, _, _, guard = load_config(); reserve(config)
    rows = [r for r in json.loads((WORK / "input_manifest.json").read_text()) if r["status"] == "QUALIFIED"]
    assert_grouped(rows)
    train = [r for r in rows if r["fold"] != fold]; test = [r for r in rows if r["fold"] == fold]
    if {r["unit_id"] for r in train} & {r["unit_id"] for r in test}:
        raise ValueError("Train/test unit overlap")
    out = WORK / f"folds/fold{fold}"; out.mkdir(parents=True, exist_ok=True)
    if (out / "COMPLETE").exists():
        return {"fold": fold, "status": "already complete"}
    caches = {r["section_id"]: load_cache(r) for r in rows}
    records = []; start = time.time()
    json_write(out / "group_receipt.json", {"fold": fold, "train": [r["section_id"] for r in train],
               "test": [r["section_id"] for r in test], "test_units": sorted({r["unit_id"] for r in test}),
               "contract_sha256": sha256(WORK / "contract.json"), "code_sha256": sha256(__file__)})
    with threadpool_limits(limits=config["training"]["threads_per_fold"]):
        for ti, target in enumerate(config["targets"]):
            for view in config["input_sizes"]:
                features = {r["section_id"]: feature_views(r, caches[r["section_id"]], ti, view) for r in rows}
                x, y, weight, indices, section_means = training_data(train, caches, features, ti, config["seed"], config)
                mean = priors(train, section_means); models = {}
                for level in ["point", "bag", "space"]:
                    ridge = make_pipeline(StandardScaler(), Ridge(alpha=config["training"]["ridge_alpha"]))
                    ridge.fit(x[level], y, standardscaler__sample_weight=weight, ridge__sample_weight=weight)
                    models[level + "_ridge"] = ridge
                    hgb = HistGradientBoostingRegressor(**config["training"]["hgb"])
                    hgb.fit(x[level], y, sample_weight=weight)
                    models[level + "_hgb"] = hgb
                if view == "full15":
                    rx, ry, rw, ri, _ = training_data(train, caches, features, ti, config["seed"] + 1, config)
                    repeat = HistGradientBoostingRegressor(**config["training"]["hgb"])
                    repeat.fit(rx["space"], ry, sample_weight=rw)
                    models["space_hgb_repeat"] = repeat
                    np.savez_compressed(out / f"subsample_repeat_t{ti}.npz", **ri)
                np.savez_compressed(out / f"subsample_t{ti}_{view}.npz", **indices)
                with (out / f"models_t{ti}_{view}.pkl").open("wb") as f:
                    pickle.dump({"models": models, "priors": mean, "target": target, "view": view,
                                 "feature_dimensions": {k: v.shape[1] for k, v in x.items()},
                                 "train_sections": [r["section_id"] for r in train]}, f, protocol=5)
                for row in test:
                    sid = row["section_id"]; cache = caches[sid]; xx = features[sid]
                    pred = {name: model.predict(xx[name.split("_")[0]]).astype(np.float32) for name, model in models.items()}
                    organ_mean = mean.get("organ:" + row["organ"], mean["global"])
                    for name, value in [("train_global_mean", mean["global"]), ("train_organ_mean", organ_mean),
                                        ("train_source_mean", mean.get("source:" + row["source_key"], organ_mean))]:
                        pred[name] = np.full(row["n_spots"], value, dtype=np.float32)
                    smooth = neighbor_mean(cache["near"], pred["bag_hgb"])
                    pred["bag_hgb_smooth_half"] = (pred["bag_hgb"] + smooth) * 0.5
                    pred["bag_hgb_smooth_all"] = smooth
                    shuffled = xx["space"].copy(); nbag = xx["bag"].shape[1]
                    perm = sample_indices(len(shuffled), f"layout|{sid}|{ti}|{view}", config["seed"], len(shuffled))
                    # sample_indices sorts indices; generate a genuine fixed permutation here.
                    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(f"layout|{config['seed']}|{sid}".encode()).digest()[:8], "little"))
                    perm = rng.permutation(perm)
                    shuffled[:, nbag:] = shuffled[perm, nbag:]
                    pred["space_hgb_shuffled_context"] = models["space_hgb"].predict(shuffled).astype(np.float32)
                    yy = cache["y"][:, ti]
                    stab = stability(pred["space_hgb"], pred["space_hgb_repeat"]) if view == "full15" else {}
                    np.savez_compressed(out / f"{sid}_t{ti}_{view}.npz", y=yy, barcodes=cache["barcodes"],
                                        source_row=cache["source_row"], xy_source_pixels=cache["xy_source_pixels"], **pred)
                    base = {"section_id": sid, "unit_id": row["unit_id"], "fold": fold, "organ": row["organ"],
                            "platform": row["metadata"]["st_technology"], "source_key": row["source_key"],
                            "disease_state": row["metadata"]["disease_state"], "patient_known": row["patient_known"],
                            "target": target, "view": view, "n_spots": len(yy), **stab}
                    records.extend({**base, "model": name, **metrics(yy, values)} for name, values in pred.items())
                print(json.dumps({"fold": fold, "target": target, "view": view,
                                  "elapsed_seconds": round(time.time() - start, 1)}), flush=True)
                del x, features, models
    tsv_write(out / "section_metrics.tsv", records)
    json_write(out / "input_receipt.json", guard.receipt())
    (out / "COMPLETE").write_text("ALL_TARGETS_AND_INPUT_SIZES_COMPLETE\n")
    return {"fold": fold, "seconds": round(time.time() - start, 1), "test_sections": len(test)}


def numeric_rows(path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        for key in ["mse", "centered_mse", "pearson", "repeat_pearson", "repeat_top_decile_jaccard"]:
            row[key] = float(row[key]) if row.get(key) not in (None, "") else None
        row["fold"] = int(row["fold"])
    return rows


def summarize():
    config = json.loads((WORK / "contract.json").read_text())
    rows = []
    for fold in config["training"]["folds"]:
        out = WORK / f"folds/fold{fold}"
        if not (out / "COMPLETE").exists():
            raise RuntimeError(f"Fold {fold} incomplete")
        receipt = json.loads((out / "group_receipt.json").read_text())
        if receipt["code_sha256"] != sha256(__file__) or receipt["contract_sha256"] != sha256(WORK / "contract.json"):
            raise RuntimeError("Code/contract changed across folds")
        rows += numeric_rows(out / "section_metrics.tsv")
    tsv_write(WORK / "section_metrics.tsv", rows)
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["target"], row["view"], row["unit_id"], row["model"])].append(row)
    units = defaultdict(dict)
    for (target, view, unit, model), values in grouped.items():
        row = values[0]
        val = {"rmse": math.sqrt(np.mean([v["mse"] for v in values])),
               "centered_rmse": math.sqrt(np.mean([v["centered_mse"] for v in values]))}
        units[(target, view, unit)][model] = val
        units[(target, view, unit)]["meta"] = {k: row[k] for k in ["fold", "organ", "source_key", "patient_known", "platform"]}
        if model == "space_hgb":
            for key in ["repeat_pearson", "repeat_top_decile_jaccard"]:
                vv = [v[key] for v in values if v[key] is not None]
                units[(target, view, unit)][key] = float(np.mean(vv)) if vv else None
    comparisons = []
    for (target, view, unit), models in units.items():
        ordinary = {name: models[name] for name in config["baselines"]}
        best = min(ordinary, key=lambda name: ordinary[name]["rmse"])
        base = ordinary[best]["rmse"]
        cbase = min(v["centered_rmse"] for v in ordinary.values())
        space = models["space_hgb"]
        comparisons.append({"target": target, "view": view, "unit_id": unit, **models["meta"],
            "n_sections": sum(r["target"] == target and r["view"] == view and r["unit_id"] == unit and r["model"] == "space_hgb" for r in rows),
            "best_ordinary_model": best, "ordinary_envelope_rmse": base, "space_hgb_rmse": space["rmse"],
            "relative_rmse_gain": (base - space["rmse"]) / max(base, 1e-12),
            "relative_centered_rmse_gain": (cbase - space["centered_rmse"]) / max(cbase, 1e-12),
            "gain_vs_bag_hgb": 1 - space["rmse"] / max(models["bag_hgb"]["rmse"], 1e-12),
            "gain_vs_point_hgb": 1 - space["rmse"] / max(models["point_hgb"]["rmse"], 1e-12),
            "shuffled_context_relative_rmse_penalty": models["space_hgb_shuffled_context"]["rmse"] / max(space["rmse"], 1e-12) - 1,
            "repeat_pearson": models["repeat_pearson"], "repeat_top_decile_jaccard": models["repeat_top_decile_jaccard"]})
    tsv_write(WORK / "unit_comparisons.tsv", comparisons)
    summaries = []; strata = []
    rule = config["screening_rule"]
    for target in config["targets"]:
        for view in config["input_sizes"]:
            rr = [r for r in comparisons if r["target"] == target and r["view"] == view]
            gain = float(np.mean([r["relative_rmse_gain"] for r in rr])); positive = float(np.mean([r["relative_rmse_gain"] > 0 for r in rr]))
            center = float(np.mean([r["relative_centered_rmse_gain"] for r in rr]))
            signs = [float(np.mean([r["relative_rmse_gain"] for r in rr if r["fold"] == f])) for f in config["training"]["folds"]]
            repeats = [r["repeat_pearson"] for r in rr if r["repeat_pearson"] is not None]
            repeat = float(np.median(repeats)) if repeats else None
            passed = (gain >= rule["min_macro_unit_relative_rmse_gain"] and positive >= rule["min_positive_unit_fraction"]
                      and sum(g > 0 for g in signs) >= rule["min_positive_folds"] and center > 0
                      and repeat is not None and repeat >= rule["full15_min_median_repeat_prediction_correlation"])
            summaries.append({"target": target, "view": view, "n_units": len(rr),
                "macro_relative_rmse_gain": gain, "positive_unit_fraction": positive,
                "macro_relative_centered_rmse_gain": center, "positive_folds": sum(g > 0 for g in signs),
                "fold_gains": signs, "median_repeat_pearson": repeat,
                "median_top_decile_jaccard": float(np.median([r["repeat_top_decile_jaccard"] for r in rr])) if repeats else None,
                "status": ("INTERNAL_CANDIDATE" if passed else "NOT_RETAINED_BY_FIXED_RULE") if view == "full15" else "INPUT_LADDER_ONLY_NO_REPEAT",
                "interpretation": "fixed molecular proxy; internal grouped exploration; not independent TLS truth or patient/source validation"})
            for key in ["organ", "source_key", "patient_known", "platform", "fold"]:
                for val in sorted({str(r[key]) for r in rr}):
                    subset = [r for r in rr if str(r[key]) == val]
                    strata.append({"target": target, "view": view, "stratum": key, "value": val, "n_units": len(subset),
                                   "macro_relative_rmse_gain": float(np.mean([r["relative_rmse_gain"] for r in subset])),
                                   "positive_unit_fraction": float(np.mean([r["relative_rmse_gain"] > 0 for r in subset]))})
    json_write(WORK / "summary.json", {"results": summaries, "not_run": config["not_run_by_design"],
                                      "primary_comparison": config["screening_rule"]["comparison"]})
    tsv_write(WORK / "relations.tsv", [{**r, "fold_gains": json.dumps(r["fold_gains"])} for r in summaries])
    tsv_write(WORK / "stratified_results.tsv", strata)
    manifest = json.loads((WORK / "input_manifest.json").read_text()); qualified = [r for r in manifest if r["status"] == "QUALIFIED"]
    report = ["# 无人工TLS标签的分子场恢复：固定首轮结果", "",
              f"合格 {len(qualified)}/{len(manifest)} 片；{len({r['unit_id'] for r in qualified})} 个冻结划分单位；全部留出点参与评价。",
              "结果来自既有来源的内部五折，不构成独立患者/来源验证。B-axis不是TLS结构真值。", "",
              "| 目标 | 可见模块 | 对最强普通基线的单位平均RMSE改善 | 正增益单位 | 正向折 | 重复相关中位数 | 判定 |",
              "|---|---:|---:|---:|---:|---:|---|"]
    for r in summaries:
        report.append(f"| {r['target']} | {r['view']} | {100*r['macro_relative_rmse_gain']:.2f}% | {100*r['positive_unit_fraction']:.1f}% | {r['positive_folds']}/5 | {r['median_repeat_pearson']} | {r['status']} |")
    report += ["", "普通基线包络按每单位最小误差取值，是偏向普通模型的保守参照；不是可部署的选模策略。",
               "候选门槛预先冻结。去均值误差另比较空间图形，重复只衡量训练子采样稳定性，不是校准不确定性。",
               "全部预测、分子定义、输入指纹、模型、分组、失败和分层结果保留。", "", "## 未完成的工作", ""]
    report += [f"- NOT_RUN: {v}" for v in config["not_run_by_design"]]
    report += ["", "## 资格排除", ""]
    for r in manifest:
        if r["status"] != "QUALIFIED":
            report.append(f"- {r['section_id']}: " + "; ".join(r["reasons"]))
    (WORK / "REPORT.md").write_text("\n".join(report) + "\n")
    (WORK / "STATUS").write_text("MOLECULAR_DISCOVERY_FROZEN_HUMAN_TLS_NOT_READ\n")
    frozen = {}
    for path in sorted(WORK.rglob("*")):
        if path.is_file() and path.name != "freeze_manifest.json" and path.suffix != ".log":
            frozen[str(path.relative_to(WORK))] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    json_write(WORK / "freeze_manifest.json", {"code_sha256": sha256(__file__), "files": frozen,
               "human_tls_labels_read": False, "final_tls_validation": "NOT_RUN"})
    print(json.dumps(summaries, ensure_ascii=False, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prepare", "run", "fold", "summarize"])
    parser.add_argument("--fold", type=int)
    args = parser.parse_args()
    if os.environ.get("CUDA_VISIBLE_DEVICES") not in (None, "", "-1"):
        raise RuntimeError("GPU disabled")
    if args.command == "prepare":
        prepare()
    elif args.command == "fold":
        print(run_fold(args.fold), flush=True)
    elif args.command == "summarize":
        summarize()
    else:
        config = json.loads((WORK / "contract.json").read_text()); reserve(config)
        if not (WORK / "input_receipt_prepare.json").exists():
            raise RuntimeError("Prepare must finish before training")
        (WORK / "STATUS").write_text("MODELS_RUNNING_NO_HUMAN_TLS_LABELS\n")
        with ProcessPoolExecutor(max_workers=config["training"]["fold_workers"],
                                 mp_context=multiprocessing.get_context("spawn")) as pool:
            futures = [pool.submit(run_fold, fold) for fold in config["training"]["folds"]]
            for future in as_completed(futures):
                print(json.dumps(future.result()), flush=True)
        summarize()


if __name__ == "__main__":
    main()
