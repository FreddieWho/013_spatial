#!/usr/bin/env python3
"""Whole-slide polygon annotator for the four-cancer TLS pilot.

RNA-free by construction: this program never opens an expression file, a score
file or any registration overlay. It reads only the whole-slide TIFF and the
per-section pixel scale, and writes two GeoJSON files in full-resolution pixel
coordinates.

    DISPLAY=:0 python3 scripts/tls_annotate_gui.py
    DISPLAY=:0 python3 scripts/tls_annotate_gui.py --sections NCBI642 ZEN44

Keys:
  1..4   class: TLS / lymphoid_aggregate / non_TLS / unresolved
  5      class: reviewed (coverage)
  Enter  close polygon      Esc  cancel polygon      u  undo last
  Space+drag or right-drag  pan        wheel  zoom
  n / p  next / previous section       Ctrl+S  save        h  help
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from PIL import Image
from PyQt5.QtCore import Qt, QPoint, QPointF, QRect, QRectF
from PyQt5.QtGui import QColor, QFont, QImage, QPainter, QPen, QPolygonF
from PyQt5.QtWidgets import (QApplication, QFileDialog, QMainWindow, QMessageBox,
                             QWidget)

Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
PACK = ROOT / "infra/tls_four_cancer_labels_20261008/annotator_pack"
ANN = ROOT / "infra/tls_four_cancer_labels_20261008/annotation/incoming"

CLASSES = {
    "1": ("TLS", QColor(220, 30, 30)),
    "2": ("lymphoid_aggregate", QColor(255, 140, 0)),
    "3": ("non_TLS", QColor(140, 140, 140)),
    "4": ("unresolved", QColor(60, 120, 220)),
    "5": ("reviewed", QColor(20, 160, 60)),
}
INSTANCE_CLASSES = {"TLS", "lymphoid_aggregate", "non_TLS", "unresolved"}
MAX_PLANE_MP = 1400.0
AUTOSAVE_SECONDS = 30


def pil_to_qimage(im: Image.Image) -> QImage:
    """PIL images cannot be handed to QPainter; QImage does not own the buffer,
    so keep an explicit copy alive."""
    im = im.convert("RGB")
    data = im.tobytes("raw", "RGB")
    qimg = QImage(data, im.size[0], im.size[1], im.size[0] * 3, QImage.Format_RGB888)
    return qimg.copy()


def load_sections() -> list[str]:
    f = PACK / "pilot_sections.tsv"
    with f.open() as fh:
        return [r["section_id"] for r in csv.DictReader(fh, delimiter="\t")]


class Slide:
    """Pyramid TIFF opened for random level access."""

    def __init__(self, sid: str, max_plane_mp: float = MAX_PLANE_MP):
        self.sid = sid
        self.max_plane_mp = max_plane_mp
        meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
        self.w = int(meta["fullres_px_width"])
        self.h = int(meta["fullres_px_height"])
        self.px_um = float(meta["pixel_size_um_estimated"])
        self.path = HEST / "wsis" / f"{sid}.tif"
        self.im = Image.open(self.path)
        self.n = getattr(self.im, "n_frames", 1)
        self.level_sizes = []
        for i in range(self.n):
            self.im.seek(i)
            self.level_sizes.append(self.im.size)
        # Finest level we will actually decode, bounded by memory.
        self.base_level = 0
        for i in range(self.n):
            lw, lh = self.level_sizes[i]
            if lw * lh / 1e6 <= self.max_plane_mp:
                self.base_level = i
                break
            self.base_level = i
        self._cache: tuple[int, Image.Image] | None = None

    def level_scale(self, level: int) -> float:
        return self.level_sizes[level][0] / self.w

    def um_per_px(self, level: int) -> float:
        return self.px_um / self.level_scale(level)

    def covering_level(self, viewport_w: int) -> int:
        """Coarsest level that still spans the viewport width."""
        for i in range(self.n - 1, -1, -1):
            if self.level_sizes[i][0] >= viewport_w:
                return i
        return 0

    def pick_level(self, viewport_w: int, zoom: float) -> int:
        """Covering level, then descend while a finer level is affordable and
        still not upscaled.

        Searching all levels matters when zoomed out: level 0 is hundreds of
        megapixels and handing that to the painter is unusably slow.
        """
        lvl = self.covering_level(viewport_w)
        while lvl > 0 and lvl > self.base_level:
            if self.level_scale(lvl - 1) <= zoom:
                lvl -= 1
            else:
                break
        return max(lvl, 0)

    def frame(self, level: int) -> Image.Image:
        if self._cache is None or self._cache[0] != level:
            self.im.seek(level)
            self._cache = (level, self.im.convert("RGB"))
        return self._cache[1]

    def crop_level(self, level: int, box, target_px: int | None = None):
        """Crop in level pixels, pre-shrinking with PIL when the target is small.

        The painter cannot cheaply downscale a hundred-megapixel source, so the
        shrink happens in PIL first and the painter gets an image near 1:1.
        """
        img = self.frame(level).crop(box)
        if target_px:
            tw, th = img.size
            if tw > 3 * target_px:
                f = target_px / tw
                img = img.resize((max(1, int(tw * f)), max(1, int(th * f))),
                                 Image.BOX)
        return img


class Store:
    """Per-section GeoJSON store in full-resolution pixel coordinates."""

    def __init__(self, sid: str):
        self.sid = sid
        self.dir = ANN / sid
        self.dir.mkdir(parents=True, exist_ok=True)
        self.instances: list[dict] = []
        self.reviewed: list[dict] = []
        self.counter = 0
        self._load()

    def _path(self, kind: str) -> Path:
        return self.dir / ("instances.geojson" if kind == "instances"
                           else "reviewed_regions.geojson")

    def _load(self) -> None:
        for kind, target in (("instances", "instances"), ("reviewed", "reviewed")):
            p = self._path(kind)
            if not p.exists():
                continue
            try:
                feats = json.loads(p.read_text()).get("features", [])
            except (json.JSONDecodeError, OSError):
                continue
            for f in feats:
                geom = f.get("geometry") or {}
                if geom.get("type") != "Polygon":
                    continue
                ring = geom["coordinates"][0]
                props = f.get("properties") or {}
                entry = {
                    "pts": [(float(x), float(y)) for x, y in ring],
                    "structure_class": props.get("structure_class"),
                    "confidence": props.get("confidence"),
                    "maturity": props.get("maturity"),
                    "note": props.get("note", ""),
                }
                (self.instances if target == "instances" else self.reviewed).append(entry)
                if target == "instances":
                    self.counter = max(self.counter, self._num(props.get("instance_id")))

    @staticmethod
    def _num(instance_id) -> int:
        try:
            return int(str(instance_id).split("_")[-1])
        except (TypeError, ValueError):
            return 0

    def add(self, kind: str, pts: list[tuple[float, float]], cls: str) -> None:
        if kind == "reviewed":
            self.reviewed.append({
                "pts": pts, "structure_class": "reviewed_region",
                "confidence": "high", "maturity": "unknown", "note": "",
            })
            return
        self.counter += 1
        self.instances.append({
            "pts": pts, "structure_class": cls, "confidence": "medium",
            "maturity": "unknown", "note": "",
        })

    def undo(self) -> bool:
        if self.instances:
            self.instances.pop()
            return True
        if self.reviewed:
            self.reviewed.pop()
            return True
        return False

    def save(self) -> None:
        for kind, feats in (("instances", self.instances),
                            ("reviewed_regions", self.reviewed)):
            out = []
            for i, e in enumerate(feats, 1):
                iid = f'{self.sid}_{kind[:3].upper()}_{i:04d}'
                out.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[[round(x, 2), round(y, 2)] for x, y in e["pts"]]],
                    },
                    "properties": {
                        "section_id": self.sid,
                        "instance_id": iid,
                        "structure_class": e["structure_class"],
                        "evidence_level": "independent_reference",
                        "reviewer_credential": "expert_human_review",
                        "double_read": "false",
                        "confidence": e["confidence"],
                        "maturity": e["maturity"],
                        "reviewer": "project_owner",
                        "annotation_protocol": "blind_first_pass_single_reader",
                        "note": e["note"],
                        "annotation_version": "v1",
                    },
                })
            doc = {"type": "FeatureCollection",
                   "name": f"{self.sid}_{kind}",
                   "features": out}
            self._path("instances" if kind == "instances"
                       else "reviewed_regions").write_text(
                json.dumps(doc, indent=1) + "\n")


class Canvas(QWidget):
    def __init__(self, win: "Win"):
        super().__init__()
        self.win = win
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)
        self.cx = 0.0            # view centre, full-resolution pixels
        self.cy = 0.0
        self.zoom = 0.02          # screen px per fullres px
        self.panning = False
        self.last = None
        self.draft: list[tuple[float, float]] = []
        self.fitted = False

    def fit(self):
        s = self.win.slide
        if s is None or self.width() < 10:
            return
        self.cx = s.w / 2
        self.cy = s.h / 2
        self.zoom = min(self.width() / s.w, self.height() / s.h) * 1.05
        self.fitted = True

    def showEvent(self, _):
        # A widget that has not been shown reports its default 640x480 size, so
        # fitting here is what makes the initial zoom correct.
        if not self.fitted:
            self.fit()
        super().showEvent(_)

    # ---- coordinate helpers ----
    def to_screen(self, x, y):
        s = self.win.slide
        sx = (x - self.cx) * self.zoom + self.width() / 2
        sy = (y - self.cy) * self.zoom + self.height() / 2
        return sx, sy

    def to_full(self, sx, sy):
        x = (sx - self.width() / 2) / self.zoom + self.cx
        y = (sy - self.height() / 2) / self.zoom + self.cy
        return x, y

    # ---- painting ----
    def paintEvent(self, _):
        s = self.win.slide
        if s is None:
            return
        p = QPainter(self)
        p.fillRect(self.rect(), QColor(24, 24, 28))
        lvl = s.pick_level(self.width(), self.zoom)
        lscale = s.level_scale(lvl)
        z = self.zoom / lscale          # level px -> screen px
        # visible box in level pixels
        lx0, ly0 = self.to_full(0, 0)
        lx1, ly1 = self.to_full(self.width(), self.height())
        box = (max(0, int(lx0 * lscale)), max(0, int(ly0 * lscale)),
               min(s.level_sizes[lvl][0], int(lx1 * lscale) + 2),
               min(s.level_sizes[lvl][1], int(ly1 * lscale) + 2))
        if box[2] > box[0] and box[3] > box[1]:
            cw, ch = box[2] - box[0], box[3] - box[1]
            target_px = max(64, int(cw * z))
            try:
                img = s.crop_level(lvl, box, target_px=target_px)
                qimg = pil_to_qimage(img)
                # Draw into an explicit rect: the crop may have been pre-shrunk by
                # PIL, so a 1:1 draw would put it at the wrong size and position.
                sx = box[0] * z + (self.width() / 2 - self.cx * self.zoom)
                sy = box[1] * z + (self.height() / 2 - self.cy * self.zoom)
                p.drawImage(QRectF(sx, sy, cw * z, ch * z), qimg)
            except Exception as exc:  # decode failure must not kill the session
                p.setPen(QColor(255, 120, 120))
                p.drawText(20, 40, f"decode error: {exc}")

        def poly(pts, color, dashed=False, width=2):
            if len(pts) < 2:
                return
            pf = QPolygonF([QPointF(*self.to_screen(x, y)) for x, y in pts])
            pen = QPen(color, width)
            if dashed:
                pen.setStyle(Qt.DashLine)
            p.setPen(pen)
            p.setBrush(Qt.NoBrush)
            p.drawPolygon(pf)

        st = self.win.store
        for e in st.reviewed:
            poly(e["pts"], QColor(20, 160, 60, 110), dashed=True, width=2)
        for e in st.instances:
            cls = e["structure_class"]
            col = dict((v[0], v[1]) for v in CLASSES.values()).get(cls, QColor(255, 0, 255))
            poly(e["pts"], col, width=2)
        if self.draft:
            col = dict((v[0], v[1]) for v in CLASSES.values()).get(
                self.win.current_class, QColor(255, 0, 255))
            poly(self.draft, col, width=2)
            if len(self.draft) >= 2:
                sx, sy = self.to_screen(*self.draft[-1])
                p.setPen(QPen(QColor(120, 120, 120), 1, Qt.DotLine))
                p.drawLine(self.to_screen(*self.draft[0]), QPoint(int(sx), int(sy)))

        # overlays
        p.setFont(QFont("monospace", 9))
        p.setPen(QColor(255, 255, 255))
        # Fixed ladder, ascending: take the smallest bar that is still readable.
        # A halving loop can never terminate when zoomed far out.
        bar_um = next((v for v in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
                       if v / s.px_um * self.zoom >= 40), 2000)
        bar_px = min(bar_um / s.px_um * self.zoom, self.width() - 30)
        p.fillRect(QRect(10, self.height() - 34,
                         int(bar_px) + 6, 26), QColor(0, 0, 0, 160))
        p.drawLine(10, self.height() - 12, int(10 + bar_px), self.height() - 12)
        p.drawText(12, self.height() - 38, f"scale bar {bar_um:.0f} um")
        p.fillRect(QRect(0, 0, self.width(), 20), QColor(0, 0, 0, 170))
        p.drawText(6, 14, self.win.status_text())

    # ---- interaction ----
    def wheelEvent(self, e):
        f = 1.25 if e.angleDelta().y() > 0 else 0.8
        mx, my = e.position().x(), e.position().y()
        fx, fy = self.to_full(mx, my)
        self.zoom = max(1e-5, min(4.0, self.zoom * f))
        self.cx += fx - (mx - self.width() / 2) / self.zoom
        self.cy += fy - (my - self.height() / 2) / self.zoom
        self.update()

    def mousePressEvent(self, e):
        self.setFocus()
        if e.button() == Qt.RightButton or (e.button() == Qt.LeftButton
                                            and self.win.space_down):
            self.panning = True
        elif e.button() == Qt.LeftButton:
            x, y = self.to_full(e.position().x(), e.position().y())
            self.draft.append((x, y))
        self.last = e.position()

    def mouseMoveEvent(self, e):
        if self.panning and self.last is not None:
            dx = e.position().x() - self.last.x()
            dy = e.position().y() - self.last.y()
            self.cx -= dx / self.zoom
            self.cy -= dy / self.zoom
            self.last = e.position()
            self.update()
        elif self.draft:
            self.update()

    def mouseReleaseEvent(self, _):
        self.panning = False
        self.last = None

    def resizeEvent(self, _):
        self.update()


class Win(QMainWindow):
    def __init__(self, sections: list[str], start: str | None,
                 max_plane_mp: float = MAX_PLANE_MP):
        super().__init__()
        self.sections = sections
        self.max_plane_mp = max_plane_mp
        self.idx = sections.index(start) if start in sections else 0
        self.space_down = False
        self.current_class = "TLS"
        self.last_save = 0.0
        self.slide = None
        self.store = None
        self.canvas = Canvas(self)
        self.setCentralWidget(self.canvas)
        self.resize(1400, 950)
        self.setWindowTitle("TLS annotation - RNA free")
        self.open(self.sections[self.idx])

    # ---- section handling ----
    def open(self, sid: str):
        self.store = Store(sid)
        self.slide = Slide(sid, self.max_plane_mp)
        self.canvas.cx = self.slide.w / 2
        self.canvas.cy = self.slide.h / 2
        self.canvas.fitted = False
        self.canvas.fit()
        self.setWindowTitle(f"{self.slide.sid}  [{self.idx + 1}/{len(self.sections)}]"
                            f"  RNA free")
        self.canvas.update()

    def next_section(self, delta: int):
        self.store.save()
        self.idx = (self.idx + delta) % len(self.sections)
        self.open(self.sections[self.idx])

    # ---- status ----
    def status_text(self):
        s, st = self.slide, self.store
        if s is None:
            return ""
        inst = len(st.instances)
        rev = len(st.reviewed)
        base_um = s.um_per_px(s.base_level)
        cur = base_um / self.canvas.zoom if self.canvas.zoom > 0 else base_um
        return (f"{s.sid} | {s.w}x{s.h}px | native {s.px_um:.3f} um/px | "
                f"finest loaded {base_um:.3f} um/px | cursor {cur:.3f} um/px | "
                f"class {self.current_class} | draft {len(self.canvas.draft)} | "
                f"instances {inst} | reviewed {rev}")

    def autosave(self):
        if time.time() - self.last_save > AUTOSAVE_SECONDS:
            self.store.save()
            self.last_save = time.time()

    def save(self):
        self.store.save()
        self.last_save = time.time()
        QMessageBox.information(self, "saved", f"{self.slide.sid} saved")

    def help(self):
        QMessageBox.information(
            self, "keys",
            "1 TLS | 2 lymphoid_aggregate | 3 non_TLS | 4 unresolved | 5 reviewed\n"
            "Left click: add polygon point | Enter: close | Esc: cancel\n"
            "Right drag or Space+drag: pan | wheel: zoom | u: undo | Ctrl+S: save\n"
            "n / p: next / previous section | h: this help\n\n"
            "Do not open score or heatmap files. Annotate tissue structure, "
            "not Visium spots.")

    def keyPressEvent(self, e):
        # Map key codes as well as text, so the tool can be driven by tests that
        # synthesise key codes without text.
        code_map = {
            Qt.Key_1: "1", Qt.Key_2: "2", Qt.Key_3: "3",
            Qt.Key_4: "4", Qt.Key_5: "5",
            Qt.Key_U: "u", Qt.Key_N: "n", Qt.Key_P: "p", Qt.Key_H: "h",
            Qt.Key_Escape: "escape",
        }
        k = e.text().lower() or code_map.get(e.key(), "")
        if k in CLASSES:
            self.current_class = CLASSES[k][0]
            self.canvas.draft.clear()
            self.canvas.update()
            self.autosave()
        elif e.key() in (Qt.Key_Return, Qt.Key_Enter):
            if len(self.canvas.draft) >= 3:
                kind = ("reviewed" if self.current_class == "reviewed"
                        else "instances")
                self.store.add(kind, list(self.canvas.draft), self.current_class)
            self.canvas.draft.clear()
            self.canvas.update()
            self.autosave()
        elif k == "escape":
            self.canvas.draft.clear()
            self.canvas.update()
        elif k == "u":
            self.store.undo()
            self.canvas.update()
            self.autosave()
        elif k == "n":
            self.next_section(1)
        elif k == "p":
            self.next_section(-1)
        elif k == "h":
            self.help()
        elif e.key() == Qt.Key_Space:
            self.space_down = True
        elif e.modifiers() & Qt.ControlModifier and k == "s":
            self.save()
        else:
            super().keyPressEvent(e)

    def keyReleaseEvent(self, e):
        if e.key() == Qt.Key_Space:
            self.space_down = False
        super().keyReleaseEvent(e)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--start", default=None)
    ap.add_argument("--max-plane-mp", type=float, default=MAX_PLANE_MP,
                    help="memory ceiling per decoded pyramid level; lower it if "
                         "the viewer stutters, raise it for finer detail")
    args = ap.parse_args()
    sections = args.sections or load_sections()
    app = QApplication(sys.argv)
    win = Win(sections, args.start, args.max_plane_mp)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()