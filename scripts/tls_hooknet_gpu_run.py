#!/usr/bin/env python3
"""Run HookNet over the staged input, with a wall-clock budget and resume.

Design constraints, because the GPU window is rented by the hour:
  * stop before the budget expires and leave a clean resume point
  * write each slide's output atomically, so a kill never leaves a partial file
  * decide slide order up front and never restart a finished slide
  * run unchanged on CPU (CUDA_VISIBLE_DEVICES=-1) so the same code path is
    tested before the window opens

Output per slide: a float16 probability mosaic for all classes on a declared
output grid, plus a JSON sidecar with the affine transform back to source pixels.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import time
from pathlib import Path

import numpy as np

PATCH = 1244
CONTEXT_FACTOR = 4
# HookNet's own writer convention (hooknet/inference/writing.py +
# inference/apply.py): TILE_SIZE = 1024 output pixels at 0.5 um/px, placed at
# (info.x - 512, info.y - 512), and the training code centre-crops the label to
# the model's output shape. Both together mean the CENTRE 1024 of the 1030-wide
# prediction corresponds to the centre of the input patch. Reproducing that
# exactly avoids guessing any offset.
TILE = 1024
TILE_OFFSET = TILE // 2          # 512
PRED_HW = 1030


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Budget:
    def __init__(self, seconds: float):
        self.t0 = time.time()
        self.limit = seconds
        self.stop_requested = False

    def elapsed(self) -> float:
        return time.time() - self.t0

    def remaining(self) -> float:
        return self.limit - self.elapsed()

    def expired(self) -> bool:
        return self.stop_requested or self.remaining() <= 0


def load_planes(entry: dict, root: Path) -> tuple[np.ndarray, np.ndarray]:
    """Target plane at 0.5 um/px and context plane at 2.0 um/px, decoded once.

    The earlier `bad_alloc` was from the CPU-fallback activation footprint, not
    from holding one extra downsampled plane (16x smaller). Per-position PIL
    resampling was slow enough to starve the GPU, so the context plane is built
    once here.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(root / entry["image"]).convert("RGB")
    target = np.array(im, dtype=np.uint8)
    cw, ch = max(1, im.size[0] // CONTEXT_FACTOR), max(1, im.size[1] // CONTEXT_FACTOR)
    context = np.asarray(im.resize((cw, ch), Image.BILINEAR), dtype=np.uint8)
    im.close()
    del im
    return target, context


def crop_target(arr: np.ndarray, cx: int, cy: int, size: int = PATCH) -> np.ndarray:
    h, w = arr.shape[:2]
    x0, y0 = cx - size // 2, cy - size // 2
    cv = np.full((size, size, 3), 255, dtype=np.uint8)
    a, b = max(0, x0), max(0, y0)
    c, d = min(w, x0 + size), min(h, y0 + size)
    if c > a and d > b:
        cv[b - y0:d - y0, a - x0:c - x0] = arr[b:d, a:c]
    return cv


def crop_from_plane(plane: np.ndarray, cx: int, cy: int, size: int = PATCH) -> np.ndarray:
    """Plain slice from a precomputed plane; the caller must pick the right
    plane for the branch (target at 0.5 um/px, context at 2.0 um/px)."""
    return crop_target(plane, cx, cy, size)


def crop(arr: np.ndarray, cx: int, cy: int, size: int) -> np.ndarray:
    """Crop centred on (cx, cy), padding with white outside the image."""
    h, w = arr.shape[:2]
    x0, y0 = cx - size // 2, cy - size // 2
    canvas = np.full((size, size, 3), 255, dtype=np.uint8)
    sx0, sy0 = max(0, x0), max(0, y0)
    sx1, sy1 = min(w, x0 + size), min(h, y0 + size)
    if sx1 > sx0 and sy1 > sy0:
        canvas[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = arr[sy0:sy1, sx0:sx1]
    return canvas


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="staged pack directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--budget-min", type=float, default=50.0,
                    help="wall-clock budget for the inference loop")
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--output-px-um", type=float, default=4.0,
                    help="resolution of the saved mosaic")
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--max-slides", type=int, default=0)
    ap.add_argument("--no-wait", action="store_true",
                    help="do not re-scan for slides that are still uploading")
    ap.add_argument("--wait-timeout-min", type=float, default=0.0)
    ap.add_argument("--cpu", action="store_true",
                    help="force CPU; default uses the GPU when one is visible")
    ap.add_argument("--weights", default=None,
                    help="path to weights.h5; if unset, search a few layouts")
    args = ap.parse_args()

    root = Path(args.input)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "maps").mkdir(exist_ok=True)
    man = json.loads((root / "manifest.json").read_text())

    if args.cpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    import tensorflow as tf
    from tensorflow.keras.layers import LeakyReLU
    from hooknet.models.tensorflowmodel import HookNet

    gpus = tf.config.list_physical_devices("GPU")
    log(f"tensorflow {tf.__version__}  GPU={gpus if gpus else 'none (CPU)'}")
    for g in gpus:
        try:
            tf.config.experimental.set_memory_growth(g, True)
        except Exception:
            pass
    model = HookNet(input_shape=[[PATCH, PATCH, 3], [PATCH, PATCH, 3]],
                    hook_indices=[2, 4], n_classes=3, batch_norm=True,
                    n_filters=32, activation=LeakyReLU(alpha=0.1))
    if args.weights:
        wpath = Path(args.weights)
    else:
        here = Path(__file__).resolve()
        candidates = [
            here.parents[1] / "infra/tls_weaklabel_probe_20261006/weights/weights.h5",
            here.parent / "weights.h5",
            here.parent.parent / "work/weights.h5",
            Path("/root/work/weights.h5"),
        ]
        wpath = next((c for c in candidates if c.exists()), candidates[0])
    if not wpath.exists():
        raise FileNotFoundError(f"weights not found; tried {candidates}")
    model.load_weights(str(wpath))
    log(f"weights loaded from {wpath}")

    order = [s["section_id"] for s in man["slides"]]
    if args.sections:
        order = [s for s in order if s in args.sections]
    if args.max_slides:
        order = order[:args.max_slides]

    budget = Budget(args.budget_min * 60)

    def on_signal(signum, frame):
        budget.stop_requested = True
        log(f"signal {signum}: finishing the current slide then stopping")
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)

    entry_of = {s["section_id"]: s for s in man["slides"]}
    done, skipped, not_yet = [], [], []
    summed = [False]
    t_infer = 0.0

    # One pass over the manifest; if some slides have not finished uploading
    # yet, re-scan so the GPU keeps working while the transfer streams in.
    done: list[str] = []
    skipped: list[str] = []
    not_yet: list[str] = []
    summed = [False]
    t_infer = 0.0
    passes = 0

    while True:
        passes += 1
        not_yet = []
        for sid in order:
            if sid in done:
                continue
            if budget.expired():
                skipped.append(sid)
                continue
            fin = out / "maps" / f"{sid}.npz"
            if fin.exists():
                done.append(sid)
                continue
            e = entry_of[sid]
            if not ((root / e["image"]).exists()
                    and (root / f"{sid}.positions.json").exists()):
                not_yet.append(sid)
                continue

            t_load = time.time()
            try:
                target, context = load_planes(e, root)
                pos = json.loads((root / f"{sid}.positions.json").read_text())
            except (MemoryError, OSError) as exc:
                log(f"{sid}: LOAD FAILED ({type(exc).__name__}); skipping, resumable")
                skipped.append(sid)
                continue
            load_s = time.time() - t_load

            step = args.output_px_um / e["um_per_px"]
            ow = int(np.ceil(e["width"] / step))
            oh = int(np.ceil(e["height"] / step))
            acc = np.zeros((oh, ow, 3), dtype=np.float32)
            cnt = np.zeros((oh, ow, 1), dtype=np.float32)

            t0 = time.time()
            n = 0
            for b0 in range(0, len(pos), args.batch):
                if budget.stop_requested and b0 > 0:
                    break
                chunk = pos[b0:b0 + args.batch]
                tb = np.stack([crop_from_plane(target, int(x), int(y)) for x, y in chunk])
                cb = np.stack([crop_from_plane(
                    context, int(x) // CONTEXT_FACTOR, int(y) // CONTEXT_FACTOR)
                    for x, y in chunk])
                # The model has its own Rescaling(scale=1/255) layer, so it wants
                # RAW 0-255 input. Dividing here too feeds it near-black and
                # yields a constant prediction; verified empirically.
                pred = model.predict([tb.astype(np.float32), cb.astype(np.float32)],
                                     verbose=0, batch_size=args.batch)
                pred = pred.reshape(len(chunk), PRED_HW, PRED_HW, -1)
                if not summed[0]:
                    s_ = pred.sum(axis=-1)
                    log(f"output shape {pred.shape}; class sum min {s_.min():.4f} "
                        f"max {s_.max():.4f} "
                        f"({'softmax-like' if abs(s_.mean()-1) < 0.05 else 'NOT softmax'})")
                    summed[0] = True
                n += len(chunk)
                crop0 = (PRED_HW - TILE) // 2
                block = int(round(step))
                for k, (x, y) in enumerate(chunk):
                    t = pred[k, crop0:crop0 + TILE, crop0:crop0 + TILE]
                    if block > 1:
                        t = t[:TILE // block * block, :TILE // block * block]
                        t = t.reshape(TILE // block, block,
                                      TILE // block, block, -1).mean((1, 3))
                    gx0 = int(round((x - TILE_OFFSET) / step))
                    gy0 = int(round((y - TILE_OFFSET) / step))
                    th_, tw_ = t.shape[0], t.shape[1]
                    a0, b0i = max(0, gx0), max(0, gy0)
                    a1, b1i = min(ow, gx0 + tw_), min(oh, gy0 + th_)
                    if a1 <= a0 or b1i <= b0i:
                        continue
                    acc[b0i:b1i, a0:a1] += t[b0i - gy0:b1i - gy0, a0 - gx0:a1 - gx0]
                    cnt[b0i:b1i, a0:a1, 0] += 1
            dt = time.time() - t0
            t_infer += dt
            del target, context
            import gc
            gc.collect()
            mean_map = np.divide(acc, np.maximum(cnt, 1), dtype=np.float32)
            log(f"{sid}: {len(pos)} positions in {dt:.0f}s "
                f"({dt/max(1,len(pos)):.3f} s/pos, load {load_s:.0f}s) "
                f"| elapsed {budget.elapsed()/60:.1f} min, left {budget.remaining()/60:.1f} min")

            try:
                tmp = out / "maps" / f".{sid}.tmp.npz"
                np.savez_compressed(tmp,
                                    prob=np.round(mean_map * 255).astype(np.uint8),
                                    coverage=np.minimum(cnt[:, :, 0], 255).astype(np.uint8))
                tmp.replace(fin)
            except (MemoryError, OSError) as exc:
                log(f"{sid}: SAVE FAILED ({type(exc).__name__}); not counted as done")
                skipped.append(sid)
                del target, context
                gc.collect()
                continue
            (out / f"{sid}.json").write_text(json.dumps({
                "section_id": sid,
                "input_image": e["image"],
                "input_um_per_px": e["um_per_px"],
                "output_um_per_px": args.output_px_um,
                "grid_step_px": step,
                "output_shape": [oh, ow],
                "n_positions": len(pos),
                "positions_completed": n,
                "batch": args.batch,
                "seconds": round(dt, 1),
                "seconds_per_position": round(dt / max(1, len(pos)), 4),
                "complete": bool(n == len(pos) and not budget.stop_requested),
                "class_order": "channel 2 = TLS per hooknettls inv_label_map {2: tls, 3: gc}",
                "placement": ("HookNet writer convention: centre 1024 px of the "
                              "1030 prediction at 0.5 um/px, placed at info.x-512"),
                "tile_px": TILE, "tile_offset_px": TILE_OFFSET,
                "mosaic_block_px": int(round(step)),
                "prob_dtype": "uint8", "prob_scale": 255,
            }, indent=1) + "\n")
            done.append(sid)

        if not not_yet or args.no_wait or budget.remaining() <= 90:
            break
        if args.wait_timeout_min and budget.elapsed() > args.wait_timeout_min * 60:
            log("wait timeout reached")
            break
        log(f"pass {passes}: {len(done)} done, waiting on {len(not_yet)} "
            f"still uploading ({not_yet[:4]})")
        time.sleep(20)

    summary = {
        "input": str(root), "out": str(out),
        "budget_min": args.budget_min,
        "elapsed_min": round(budget.elapsed() / 60, 2),
        "slides_completed": done,
        "slides_skipped_out_of_time": skipped,
        "slides_not_yet_uploaded": not_yet,
        "positions": sum(entry_of[s]["n_positions"] for s in done),
        "inference_seconds": round(t_infer, 1),
        "stopped_by_signal": budget.stop_requested,
        "resume_hint": "rerun the same command; finished slides are skipped",
    }
    (out / "run_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    log(f"completed {len(done)}, skipped {len(skipped)}")
    log(f"resume with the same command; done slides are skipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())