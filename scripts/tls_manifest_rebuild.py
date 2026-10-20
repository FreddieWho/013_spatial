"""Rebuild four-cancer manifest from scores npz files on disk."""
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path("/home/huyudi/013_spatial")
OUT = ROOT / "infra/tls_fourcancer_20261006"
rows = list(csv.DictReader(open(ROOT / "data/other_sources/hest1k/selected_samples.tsv"), delimiter="\t"))
organ_of = {r["id"]: r["organ"] for r in rows}

manifest = []
for sp in sorted((OUT / "scores").glob("hest-*.npz")):
    sid = sp.stem
    raw = sid.removeprefix("hest-")
    z = np.load(sp)
    valid = z["valid"] if "valid" in z.files else None
    manifest.append({
        "section_id": sid,
        "organ": organ_of.get(raw, "?"),
        "status": "OK",
        "n_spots": int(z["scores"].shape[0]),
        "n_valid_sets": None if valid is None else int(np.asarray(valid).sum()),
    })
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
ok = sum(1 for m in manifest if m["status"] == "OK")
(OUT / "STATUS").write_text(f"SCORED {ok}/{len(manifest)}\n")
print(f"rebuilt {len(manifest)} ok={ok}", flush=True)
