#!/usr/bin/env python3
"""W0 deep-priority Synapse batch downloader with status TSV + throughput stats.
Usage: conda run -n htan_download python W0_run.py [--jobs 4] [--limit N]
Reads /tmp/deep_priority.json (or path via --input). Resumes: skips files
already present with size>0 (ALREADY_PRESENT_VERIFIED by size; manifest has
no checksum column). Retries 3x with exponential backoff, then FAILED_RETRYABLE.
"""
import argparse, csv, json, os, sys, time, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent  # WUSTL/
STATUS = BASE / "synapse_download_status.tsv"
COLS = ["atlas", "synapse_id", "filename", "assay", "biospecimen",
        "destination", "status", "bytes", "checksum", "attempts",
        "error", "download_timestamp"]

def load_done():
    done = {}
    if STATUS.exists():
        with open(STATUS) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("status") in ("DOWNLOADED", "ALREADY_PRESENT_VERIFIED"):
                    done[r["synapse_id"]] = r
    return done

def one(syn, rec, dest_base):
    sid = rec["syn"]
    dest_dir = dest_base / sid
    dest_dir.mkdir(parents=True, exist_ok=True)
    # resume check
    for p in dest_dir.iterdir():
        if p.is_file() and p.stat().st_size > 0:
            return {"atlas": "WUSTL", "synapse_id": sid, "filename": p.name,
                    "assay": rec["assay"], "biospecimen": rec["bio"],
                    "destination": str(p), "status": "ALREADY_PRESENT_VERIFIED",
                    "bytes": p.stat().st_size, "checksum": "", "attempts": 0,
                    "error": "", "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}
    last_err, attempts = "", 0
    for attempt in range(1, 4):
        attempts = attempt
        try:
            t0 = time.time()
            ent = syn.get(sid, downloadLocation=str(dest_dir))
            dt = time.time() - t0
            sz = os.path.getsize(ent.path)
            return {"atlas": "WUSTL", "synapse_id": sid, "filename": os.path.basename(ent.path),
                    "assay": rec["assay"], "biospecimen": rec["bio"],
                    "destination": ent.path, "status": "DOWNLOADED",
                    "bytes": sz, "checksum": "", "attempts": attempts,
                    "error": f"transfer_s={dt:.1f}",
                    "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}
        except Exception as e:
            last_err = f"{type(e).__name__}: {str(e)[:160]}"
            if "403" in last_err or "401" in last_err or "ACCESS" in last_err.upper():
                status = "ACCESS_REQUIRED"
                break
            status = "FAILED_RETRYABLE"
            time.sleep(2 ** attempt)
    else:
        status = "FAILED_RETRYABLE"
    return {"atlas": "WUSTL", "synapse_id": sid, "filename": rec["fn"],
            "assay": rec["assay"], "biospecimen": rec["bio"],
            "destination": str(dest_dir), "status": status,
            "bytes": 0, "checksum": "", "attempts": attempts,
            "error": last_err, "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--input", default="/tmp/deep_priority.json")
    a = ap.parse_args()
    import synapseclient
    recs = json.load(open(a.input))
    if a.limit:
        recs = recs[:a.limit]
    done = load_done()
    todo = [r for r in recs if r["syn"] not in done]
    print(f"total={len(recs)} done={len(done)} todo={len(todo)} jobs={a.jobs}", flush=True)
    syn = synapseclient.Synapse(silent=True)
    syn.login(silent=True)
    print(f"auth as {syn.username}", flush=True)
    dest_base = BASE / "synapse"
    t_start, total_bytes = time.time(), 0
    results = list(done.values())
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(one, syn, r, dest_base): r for r in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            try:
                res = fut.result()
            except Exception:
                r = futs[fut]
                res = {"atlas": "WUSTL", "synapse_id": r["syn"], "filename": r["fn"],
                       "assay": r["assay"], "biospecimen": r["bio"],
                       "destination": "", "status": "FAILED_TERMINAL",
                       "bytes": 0, "checksum": "", "attempts": 0,
                       "error": traceback.format_exc(limit=3)[-300:],
                       "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}
            results.append(res)
            if res["status"] == "DOWNLOADED":
                total_bytes += int(res["bytes"] or 0)
            el = time.time() - t_start
            dl = sum(1 for x in results if x["status"] == "DOWNLOADED")
            print(f"[{i}/{len(todo)}] {res['synapse_id']} {res['status']} "
                  f"dl_bytes={total_bytes} elapsed={el:.0f}s", flush=True)
    with open(STATUS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    el = time.time() - t_start
    from collections import Counter
    print(f"DONE elapsed={el:.0f}s dl_bytes={total_bytes} "
          f"avg_MBps={total_bytes/el/1e6:.2f} status={Counter(r['status'] for r in results)}")

if __name__ == "__main__":
    sys.exit(main())
