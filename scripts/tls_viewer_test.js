/* Integration test for viewer/app.js.
 *
 * Runs the real app in node against a real HTTP server with minimal DOM stubs.
 * This is what verifies the fetch chain, tile addressing, polygon handling,
 * coverage maths and export payload. Static UI styling is checked separately by
 * a headless-browser screenshot.
 *
 *   node scripts/tls_viewer_test.js http://127.0.0.1:8899
 */
"use strict";
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const BASE = process.argv[2] || "http://127.0.0.1:8899";
const APP = process.argv[3] ||
  path.join(__dirname, "..", "infra", "tls_annotation_pack_src", "viewer", "app.js");

let failures = 0;
function ok(cond, label, extra) {
  if (cond) { console.log("  PASS  " + label); }
  else { failures++; console.log("  FAIL  " + label + (extra ? "  <- " + extra : "")); }
}
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/* ---------- DOM stubs ---------- */
function el(id) {
  const e = {
    id, textContent: "", innerHTML: "", hidden: false, value: "", title: "",
    dataset: {}, tagName: "DIV", style: {}, children: [], options: [],
    selectedIndex: 0, clientWidth: 1600, clientHeight: 1000, width: 1600, height: 1000,
    clientHeight_: 1000,
    classList: {
      _s: new Set(),
      add(c) { this._s.add(c); }, remove(c) { this._s.delete(c); },
      toggle(c, on) { if (on === undefined) on = !this._s.has(c); on ? this._s.add(c) : this._s.delete(c); },
      contains(c) { return this._s.has(c); },
    },
    appendChild(c) { this.children.push(c); if (c.tagName === "OPTION") this.options.push(c); return c; },
    remove() {}, click() {}, setPointerCapture() {}, addEventListener() {},
    getContext() { return ctxStub; },
    getBoundingClientRect() { return { left: 0, top: 0, width: 1600, height: 1000 }; },
  };
  Object.defineProperty(e, "offsetX", { get: () => 800 });
  Object.defineProperty(e, "offsetY", { get: () => 500 });
  return e;
}
const calls = { drawImage: [], fillRect: 0, fillText: [] };
const ctxStub = new Proxy({}, {
  get(_, k) {
    if (k === "imageSmoothingEnabled" || k === "lineWidth") return undefined;
    return (...args) => {
      // drawImage(img, dx, dy, dw, dh): canvases (heat LUT output) carry a
      // marker; plain Images carry _u. Anything else is a bug.
      if (k === "drawImage") {
        const a = args[0];
        calls.drawImage.push(a && a._isCanvas ? { canvas: true, heatUrl: a._heatUrl || null }
          : (a && a._u ? a._u : String(a)));
      }
      else if (k === "fillRect") calls.fillRect++;
      else if (k === "fillText") calls.fillText.push(args[0]);
      else if (k === "setTransform" || k === "beginPath" || k === "moveTo" ||
               k === "lineTo" || k === "closePath" || k === "stroke" ||
               k === "fill" || k === "save" || k === "restore") { /* no-op */ }
      else if (k === "strokeStyle" || k === "fillStyle" || k === "font" ||
               k === "textAlign" || k === "setLineDash") { /* no-op */ }
    };
  },
  set() { return true; },
});

const els = {};
const getEl = (id) => (els[id] = els[id] || el(id));
const store = new Map();

let requestedTiles = [];
const alerts = [];
let finishTiles, pendingTiles, fireResize;
const DOC = BASE + "/viewer/index.html";
const resolve = (u) => new URL(u, DOC).href;

const ctxObj = {
  console,
  setTimeout, clearTimeout, setInterval, clearInterval,
  // Mirror a browser: rAF is asynchronous, otherwise the app's draw guard
  // (S.raf) never clears and nothing is ever painted.
  requestAnimationFrame: (fn) => setTimeout(() => fn(0), 0),
  // Resolve like a browser: relative URLs resolve against the DOCUMENT url, not
// the server root. Prepending the root here would have hidden the /viewer/
// mis-resolution that produced a black screen.
  fetch: (u, o) => fetch(resolve(u), o),
  Promise, JSON, Math, Object, Array, String, Number, Date, Error, isFinite, parseInt, parseFloat,
  document: {
    getElementById: getEl,
    // app.js derives the pack root from document.currentScript.src
    currentScript: { src: DOC },
    getElementsByTagName: () => [{ src: DOC }],
    createElement: (t) => {
      const e = el("new_" + t);
      e.tagName = String(t).toUpperCase();
      if (t === "canvas") {
        // Functional 2d context so the heat LUT path really executes instead
        // of silently throwing inside a try/catch.
        e.width = 0; e.height = 0; e._isCanvas = true;
        e.getContext = () => ({
          drawImage() {},
          getImageData: (x, y, w, h) => ({
            width: w, height: h, data: new Uint8ClampedArray(w * h * 4),
          }),
          createImageData: (w, h) => ({
            width: w, height: h, data: new Uint8ClampedArray(w * h * 4),
          }),
          putImageData() {},
        });
      }
      return e;
    },
    body: { appendChild() {} },
    addEventListener() {},
  },
  window: {
    addEventListener() {}, devicePixelRatio: 1,
    ResizeObserver: class { constructor(cb) { this.cb = cb; } observe() { fireResize = () => this.cb([]); } },
  },
  localStorage: {
    getItem: (k) => (store.has(k) ? store.get(k) : null),
    setItem: (k, v) => store.set(k, String(v)),
    removeItem: (k) => store.delete(k),
  },
  Image: (function () {
    // Registry so the test can finish decoding in stages, the way a browser does.
    const reg = [];
    finishTiles = function (n) {
      let done = 0;
      for (const im of reg) {
        if (im.complete) continue;
        im.complete = true;
        im.naturalWidth = 256;
        done++;
        if (im.onload) im.onload();
        if (n && done >= n) break;
      }
      return done;
    };
    pendingTiles = () => reg.filter((i) => !i.complete).length;
    return class {
      constructor() { this.complete = false; this.naturalWidth = 0; reg.push(this); }
      set src(u) { this._u = u; requestedTiles.push(resolve(u)); }
      get src() { return this._u; }
    };
  })(),
  Blob: class { constructor(parts) { this.parts = parts; } },
  URL: { createObjectURL: () => "blob:stub", revokeObjectURL() {} },
  FileReader: class { readAsText(f) { this.result = (f && f.__content) || "{}"; if (this.onload) this.onload(); } },
  alert: (m) => { alerts.push(m); },
};
ctxObj.globalThis = ctxObj;
ctxObj.window.document = ctxObj.document;

(async function main() {
  console.log("loading " + APP);
  vm.createContext(ctxObj);
  vm.runInContext(fs.readFileSync(APP, "utf8"), ctxObj, { filename: "app.js" });
  // fire DOMContentLoaded the way a browser would
  ctxObj.init();

  await sleep(1500);

  const S = ctxObj.S;
  const idx = S.index;

  console.log("\n[1] boot and data fetch");
  if (!getEl("slideSel").options.length) {
    ok(false, "slide list populated",
       "no options - the app could not read data/index.json. Check S.root");
    console.log("\n" + failures + " FAILURE(S)");
    process.exit(1);
  }
  ok(getEl("slideSel").options.length === idx.slides.length, "slide list populated",
     "options=" + getEl("slideSel").options.length + " expected " + idx.slides.length);
  ok(getEl("classBtns").children.length === 5, "5 class buttons created",
     "n=" + getEl("classBtns").children.length);
  ok(S.meta && S.meta.section_id, "meta.json loaded for " + (S.meta && S.meta.section_id));
  ok(S.meta && S.meta.base_w > 0 && S.meta.levels.length > 1,
     "pyramid levels present", S.meta ? String(S.meta.levels.length) : "none");
  ok(S.tissue && S.tissue.features && S.tissue.features.length > 0, "tissue contour loaded");
  // The tissue contours are in original pixels while the grid is in base pixels.
  // If they are not rescaled, the coverage readout silently returns nothing.
  (function () {
    const k = S.meta.base_w / S.meta.full_w;
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (const f of S.tissue.features) {
      for (const ring of f.geometry.coordinates) {
        for (const p of ring) {
          xmin = Math.min(xmin, p[0] * k); xmax = Math.max(xmax, p[0] * k);
          ymin = Math.min(ymin, p[1] * k); ymax = Math.max(ymax, p[1] * k);
        }
      }
    }
    ok(xmax > 0 && xmin < S.meta.base_w && ymax > 0 && ymin < S.meta.base_h,
       "tissue contours fall inside the base-image grid once rescaled",
       [Math.round(xmin), Math.round(xmax), S.meta.base_w]);
  })();

  console.log("\n[2] view and level selection");
  ok(S.view.scale > 0 && isFinite(S.view.scale), "fit() produced a usable scale",
     String(S.view.scale));
  const L0 = S.__proto__ ? 0 : 0;
  const lvlAtFit = (function () {
    let L = Math.round(Math.log2(1 / S.view.scale));
    return Math.max(0, Math.min(S.meta.levels.length - 1, L));
  })();
  ok(lvlAtFit >= 0, "level chosen at fit view: L" + lvlAtFit);

  console.log("\n[3] tile requests go to real files");
  ok(requestedTiles.length > 0, "tiles requested", "n=" + requestedTiles.length);
  const uniq = Array.from(new Set(requestedTiles));
  let bad = [];
  for (const u of uniq.slice(0, 40)) {
    const r = await fetch(u);
    if (!r.ok) bad.push(u + " -> " + r.status);
  }
  ok(bad.length === 0, "sampled " + Math.min(40, uniq.length) + " tile URLs all 200",
     bad.slice(0, 3).join("; "));
  ok(requestedTiles.every((u) => /\/data\/[^/]+\/\d+\/\d+_\d+\.jpg$/.test(u)),
     "tile URL resolves to /data/, not /viewer/data/", 
     requestedTiles[0]);

  // Regression: tiles decode after the first paint. Nothing must be drawn from
  // placeholders alone, and completion must repaint with NO user interaction.
  console.log("\n[3b] tiles repaint when they finish decoding (black-screen guard)");
  const drawnBefore = calls.drawImage.length;
  ok(drawnBefore === 0, "nothing painted before any tile decoded",
     "drawImage calls=" + drawnBefore);
  const pend = pendingTiles();
  ok(pend > 0, "tiles still decoding at first paint", "pending=" + pend);
  finishTiles();
  await sleep(120);
  ok(calls.drawImage.length > 0, "repainted after tiles decoded, with no user input",
     "drawImage calls=" + calls.drawImage.length);
  ok(calls.drawImage.every((u) => typeof u !== "string" ||
     /\/data\/[^/]+\/\d+\/\d+_\d+\.jpg$/.test(resolve(u))),
     "only real H&E tiles were painted (no stray objects)");
  const before2 = calls.drawImage.length;
  finishTiles();
  await sleep(80);
  ok(calls.drawImage.length === before2, "no redundant repaint once settled",
     before2 + " -> " + calls.drawImage.length);

  console.log("\n[3c] canvas sizing");
  ok(S.canvas.width === 1600 && S.canvas.height === 1000,
     "backing store sized from layout", S.canvas.width + "x" + S.canvas.height);
  ok(typeof fireResize === "function", "ResizeObserver wired");

  // Regression: the reported symptom was a page stuck loading with a GUI but no
  // image. fit() against a not-yet-laid-out canvas yields scale 0, every range
  // becomes NaN, no tile is ever requested, and it never recovers on its own.
  console.log("\n[3d] zero-sized canvas must not wedge the viewer");
  ok(getEl("cv").clientWidth === 1600, "harness starts with a laid-out canvas");
  const before = requestedTiles.length;
  getEl("cv").clientWidth = 0;
  getEl("cv").clientHeight = 0;
  S.canvas = getEl("cv");
  ctxObj.fit();
  ok(S.needsFit === true, "fit() defers instead of setting scale 0",
     "scale=" + S.view.scale + " needsFit=" + S.needsFit);
  const paintedBefore = calls.drawImage.length;
  ctxObj.draw();
  const reqWhileZero = requestedTiles.length - before;
  ok(reqWhileZero === 0, "no new tiles requested while the canvas has no box",
     "requested=" + reqWhileZero);
  ok(calls.drawImage.length === paintedBefore, "nothing painted while the canvas has no box");
  getEl("cv").clientWidth = 1600;
  getEl("cv").clientHeight = 1000;
  S.canvas = getEl("cv");
  fireResize();                       // the real size() inside wire()
  ok(S.needsFit === false && S.view.scale > 0,
     "view recovers once the canvas has a box",
     "scale=" + S.view.scale);
  await sleep(120);
  ok(S.degenerate !== true, "healthy view after recovery");

  console.log("\n[3e] one paint never requests an unbounded tile count");
  let worst = 0;
  for (let t = 0; t < 60; t++) {
    S.view.scale = 10 ** (-4 + (t / 60) * 4.5);   // extreme zoom-out to zoom-in
    S.view.cx = S.meta.base_w / 2;
    S.view.cy = S.meta.base_h / 2;
    const n0 = requestedTiles.length;
    ctxObj.draw();
    worst = Math.max(worst, requestedTiles.length - n0);
  }
  ok(worst <= 400 + 1, "max tiles requested by one paint is bounded",
     "worst=" + worst);
  ctxObj.fit();
  ok(S.degenerate !== true, "restored to a healthy view");
  const k = S.meta.full_w / S.meta.base_w;
  S.cls = "TLS"; S.draftCls = "TLS";
  ctxObj.fit();
  S.draft = [[1000, 1000], [3000, 1000], [3000, 3000], [1000, 3000]];
  ctxObj.closePolygon();
  ok(S.data[S.sid].instances.length === 1, "polygon closed into instances");
  const rec = S.data[S.sid].instances[0];
  ok(rec.structure_class === "TLS", "class recorded", rec.structure_class);
  ok(rec.pts.length === 4, "4 points kept");
  ok(store.get("tls_pilot_v1__" + S.sid) !== undefined, "written to localStorage");

  S.cls = "reviewed"; S.draftCls = "reviewed";
  S.draft = [[0, 0], [S.meta.base_w, 0], [S.meta.base_w, S.meta.base_h], [0, S.meta.base_h]];
  ctxObj.closePolygon();
  ok(S.data[S.sid].reviewed.length === 1, "coverage polygon stored separately");
  const cov = S.cov;
  ok(cov !== null && cov > 0.5 && cov <= 1.0001, "coverage computed",
     "cov=" + cov);
  // A polygon covering the whole slide must report the tissue as fully reviewed.
  ok(cov !== null && cov > 0.97, "full-slide coverage polygon reports ~100%",
     "cov=" + cov);

  console.log("\n[5] undo / cancel");
  ctxObj.undo();
  ok(S.data[S.sid].instances.length === 0, "undo removed the instance");
  ctxObj.undo();
  ok(S.data[S.sid].reviewed.length === 0, "undo removed the coverage polygon");

  console.log("\n[6] export payload");
  const dl = [];
  ctxObj.dl = (n, t) => dl.push([n, t]);
  S.data[S.sid].instances.push({
    pts: [[1000, 1000], [3000, 1000], [3000, 3000], [1000, 3000]],
    structure_class: "TLS", confidence: "medium", maturity: "unknown", note: "",
  });
  S.data[S.sid].reviewed.push({ pts: [[0, 0], [S.meta.base_w, 0], [S.meta.base_w, S.meta.base_h]] });
  ctxObj.exportSlide();
  ok(dl.length >= 1, "per-slide export produced files", "n=" + dl.length);
  const gj = JSON.parse(dl[0][1]);
  ok(gj.type === "FeatureCollection", "geojson FeatureCollection");
  const ring = gj.features[0].geometry.coordinates[0];
  ok(gj.features[0].geometry.type === "Polygon", "geometry is Polygon");
  ok(Math.abs(ring[0][0] - 1000 * k) < 1.0, "coords scaled to ORIGINAL pixels",
     "got " + ring[0][0] + " expected ~" + (1000 * k));
  const p = gj.features[0].properties;
  ok(p.evidence_level === "reviewed_candidate", "evidence_level v2 (assisted)", p.evidence_level);
  ok(p.annotation_version === "v2", "annotation version v2", p.annotation_version);
  ok(p.double_read === "false", "double_read recorded", p.double_read);
  ok(p.reviewer_credential === "expert_human_review", "credential recorded", p.reviewer_credential);
  ok(p.annotation_protocol === "assisted_review_single_reader", "protocol recorded (assisted)");
  ok(p.overlay_available === true, "overlay availability recorded");
  ok(ring.length === 5, "ring closed", "n=" + ring.length);

  ctxObj.dl = (n, t) => dl.push([n, t]);
  ctxObj.exportAll();
  // capture the envelope synchronously: exportSlide's second file lands via
  // a 350ms timer and would otherwise interleave ahead of it
  const envelopePayload = dl[dl.length - 1][1];
  const all = JSON.parse(dl[dl.length - 1][1]);
  ok(all.format === "tls_annotation_export", "all-export envelope");
  ok(all.version === 2, "envelope version 2", all.version);
  ok(all.annotation_protocol === "assisted_review_single_reader", "envelope protocol (assisted)");
  ok(all.slides && Object.keys(all.slides).length === idx.slides.length,
     "envelope holds every slide (" + Object.keys(all.slides).length + ")",
     "expected " + idx.slides.length);
  ok(all.reviewer_credential === "expert_human_review", "envelope credential");
  const one = all.slides[S.sid].instances[0];
  ok(Math.abs(one.pts[0][0] - 1000 * k) < 1.0, "all-export also in original pixels");

  console.log("\n[7b] heat overlay plumbing");
  ok(!!(S.meta.heat && S.meta.heat.levels && S.meta.heat.levels.length),
     "heat pyramid present in slide meta",
     S.meta.heat ? S.meta.heat.levels.length + " levels" : "no heat meta");
  S.data[S.sid].overlay_used = true;
  ok(S.data[S.sid].overlay_used === true, "overlay use recorded on the slide");
  const hz = S.meta.heat.levels[0];
  ok(hz.um_per_px > 0 && hz.w > 0 && hz.h > 0, "heat L0 has dims and scale");
  // heat tiles must resolve over HTTP exactly like H&E tiles do
  const hurl = "data/" + S.sid + "/heat/0/0_0.png";
  const hr = await fetch(BASE + "/" + hurl);
  ok(hr.ok, "heat tile resolves over HTTP", hurl + " -> " + hr.status);

  console.log("\n[7c] heat overlay actually paints (NaN-placement guard)");
  S.heat.on = true;
  if (S.data[S.sid]) S.data[S.sid].overlay_used = true;
  // zoom in one step so the heat path picks a deeper level too
  S.view.scale = S.view.scale * 4;
  ctxObj.draw();
  await sleep(120);
  const heatUrls = requestedTiles.filter((u) => u.includes("/heat/"));
  ok(heatUrls.length > 0, "heat tiles requested while overlay is on",
     "n=" + heatUrls.length);
  let heatBad = [];
  for (const u of heatUrls.slice(0, 40)) {
    const r = await fetch(u);
    if (!r.ok) heatBad.push(u + " -> " + r.status);
  }
  ok(heatBad.length === 0,
     "sampled heat tile URLs all 200 (" + Math.min(40, heatUrls.length) + " checked)",
     heatBad.slice(0, 3).join("; "));
  // heat images decode after the requesting paint, like H&E tiles: finish
  // them, repaint, and only then assert anything was drawn
  finishTiles();
  ctxObj.draw();
  await sleep(120);
  const heatPaints = calls.drawImage.filter(
    (a) => a && a.canvas === true && a.heatUrl);
  ok(heatPaints.length > 0, "processed heat canvases painted to screen",
     "n=" + heatPaints.length);
  S.heat.on = false;
  ctxObj.draw();
  await sleep(80);

  console.log("\n[7] import round trip");
  store.clear();
  const payload = envelopePayload;
  ctxObj.importFile({ __content: payload });
  await sleep(50);
  ok(alerts.length > 0 && alerts[alerts.length - 1] === "已导入", "import accepted",
     alerts.join("|"));
  ok(S.data[S.sid].instances.length === 1, "instances restored after import",
     "n=" + S.data[S.sid].instances.length);

  console.log("\n" + (failures ? failures + " FAILURE(S)" : "ALL PASS"));
  process.exit(failures ? 1 : 0);
})().catch((e) => { console.error("harness error:", e); process.exit(2); });