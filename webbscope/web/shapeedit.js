"use strict";
// Selecting and editing the boxes on a slide preview (shared by the single-run and the batch import pages).
// A shape is {x0,y0,x1,y1,name,keep} (rectangle) or {pts:[[x,y],...],name,keep} (free-hand), in data coordinates (micrometres).
// A selection is {kind:"shape", i} (index into the shapes) or {kind:"comp", id} (a detected tissue box, id of the piece).
const SE = (() => {
  const norm = (r) => { if (!r.pts) { const a = Math.min(r.x0, r.x1), b = Math.max(r.x0, r.x1), c = Math.min(r.y0, r.y1), d = Math.max(r.y0, r.y1); r.x0 = a; r.x1 = b; r.y0 = c; r.y1 = d; } };
  const handles = (r) => (r.pts ? r.pts : [[r.x0, r.y0], [r.x1, r.y0], [r.x1, r.y1], [r.x0, r.y1]]);
  const inPoly = (p, x, y) => { let c = false; for (let i = 0, j = p.length - 1; i < p.length; j = i++) if ((p[i][1] > y) !== (p[j][1] > y) && x < ((p[j][0] - p[i][0]) * (y - p[i][1])) / (p[j][1] - p[i][1]) + p[i][0]) c = !c; return c; };
  const segDist = (x, y, a, b) => { const dx = b[0] - a[0], dy = b[1] - a[1], L = dx * dx + dy * dy, t = L ? Math.max(0, Math.min(1, ((x - a[0]) * dx + (y - a[1]) * dy) / L)) : 0; return Math.hypot(x - a[0] - t * dx, y - a[1] - t * dy); };
  const hitShape = (r, x, y, tol) => {
    if (!r.pts) return x >= Math.min(r.x0, r.x1) - tol && x <= Math.max(r.x0, r.x1) + tol && y >= Math.min(r.y0, r.y1) - tol && y <= Math.max(r.y0, r.y1) + tol;
    if (inPoly(r.pts, x, y)) return true;
    for (let i = 0; i < r.pts.length; i++) if (segDist(x, y, r.pts[i], r.pts[(i + 1) % r.pts.length]) <= tol) return true;
    return false;
  };
  const hitHandle = (r, x, y, tol) => { const h = handles(r); let best = -1, bd = tol; h.forEach((q, i) => { const d = Math.hypot(q[0] - x, q[1] - y); if (d <= bd) { bd = d; best = i; } }); return best; };
  const move = (r, dx, dy) => { if (r.pts) { r.pts = r.pts.map((q) => [q[0] + dx, q[1] + dy]); if (r.holes) r.holes = r.holes.map((h) => h.map((q) => [q[0] + dx, q[1] + dy])); } else { r.x0 += dx; r.x1 += dx; r.y0 += dy; r.y1 += dy; } };
  // the outline of a free-hand shape as one path: its outer ring and its holes (a tissue box loaded from a GeoJSON can have holes); stroke it afterwards
  function ringsPath(ctx, r, toPx) {
    ctx.beginPath();
    for (const ring of [r.pts, ...(r.holes || [])]) { ring.forEach((q, k) => { const [x, y] = toPx(q[0], q[1]); k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }); ctx.closePath(); }
  }
  const dragHandle = (r, h, x, y) => { if (r.pts) r.pts[h] = [x, y]; else { if (h === 0 || h === 3) r.x0 = x; else r.x1 = x; if (h === 0 || h === 1) r.y0 = y; else r.y1 = y; } };
  // the shape under the point (the last drawn wins), else the smallest detected box around it
  function pick(shapes, boxes, x, y, tol) {
    for (let i = shapes.length - 1; i >= 0; i--) if (hitShape(shapes[i], x, y, tol)) return { kind: "shape", i };
    let best = null, ba = Infinity;
    for (const b of boxes) { const [x0, y0, x1, y1] = b.bbox; if (x >= x0 - tol && x <= x1 + tol && y >= y0 - tol && y <= y1 + tol) { const a = (x1 - x0) * (y1 - y0); if (a < ba) { ba = a; best = { kind: "comp", id: b.id }; } } }
    return best;
  }
  // a readable label: dark outline, then the text, kept inside the canvas
  function label(ctx, text, x, y, color) { ctx.save(); ctx.font = "600 12px system-ui"; ctx.lineJoin = "round"; ctx.lineWidth = 3.5; ctx.strokeStyle = "rgba(0,0,0,.85)"; ctx.strokeText(text, x, y); ctx.fillStyle = color; ctx.fillText(text, x, y); ctx.restore(); }
  // detected tissue boxes: dashed (red when switched off), labelled with the tissue's name in the tissue's colour; the selected one solid yellow
  // info(id) -> {name, keep, color?}
  // cap a shape at the slide b = [x0, y0, x1, y1] in place: a rectangle is clamped, a free-hand outline is clipped to the rectangle (Sutherland-Hodgman);
  // returns false when nothing of the shape is left on the slide
  function clip(r, b) {
    if (!b) return true;
    if (!r.pts) {
      norm(r); r.x0 = Math.max(r.x0, b[0]); r.y0 = Math.max(r.y0, b[1]); r.x1 = Math.min(r.x1, b[2]); r.y1 = Math.min(r.y1, b[3]);
      return r.x1 > r.x0 && r.y1 > r.y0;
    }
    const clipRing = (ring) => {
      let p = ring;
      const edges = [[0, b[0], 1], [0, b[2], -1], [1, b[1], 1], [1, b[3], -1]];   // [axis, limit, +1 keep >= limit / -1 keep <= limit]
      for (const [k, lim, sg] of edges) {
        const inn = (q) => sg * (q[k] - lim) >= 0, out = [];
        for (let i = 0; i < p.length; i++) {
          const a = p[(i + p.length - 1) % p.length], c = p[i];
          if (inn(c) !== inn(a)) { const t = (lim - a[k]) / (c[k] - a[k]), q = [a[0] + t * (c[0] - a[0]), a[1] + t * (c[1] - a[1])]; q[k] = lim; out.push(q); }
          if (inn(c)) out.push(c);
        }
        p = out; if (!p.length) break;
      }
      p = p.filter((q, i) => { const n = p[(i + 1) % p.length]; return q[0] !== n[0] || q[1] !== n[1]; });
      return p.length >= 3 ? p : null;
    };
    const p = clipRing(r.pts); if (!p) return false; r.pts = p;
    if (r.holes) { r.holes = r.holes.map(clipRing).filter(Boolean); if (!r.holes.length) delete r.holes; }
    return true;
  }
  // polys = [[ring, ...], ...] (rings of [x, y] in data units): the tissue's box after the overlap with its neighbours was removed
  // (what the export writes); strokes it and returns its top-left in pixels, or null when there is nothing to draw
  function strokePolys(ctx, polys, toPx) {
    if (!polys || !polys.length) return null;
    let tx = Infinity, ty = Infinity; ctx.beginPath();
    for (const poly of polys) for (const ring of poly) ring.forEach((q, k) => { const [x, y] = toPx(q[0], q[1]); if (y < ty || (y === ty && x < tx)) { tx = x; ty = y; } k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
    ctx.closePath(); ctx.stroke(); return [tx, ty];
  }
  // polysOf(id, name) -> the tissue's carved box (see strokePolys) or null to draw the plain bounding box; a tissue made of several pieces is outlined once
  function paintBoxes(ctx, boxes, toPx, info, sel, polysOf) {
    ctx.save(); ctx.font = "11px system-ui"; const done = new Set();
    for (const b of boxes) {
      const i = info(b.id); if (!i) continue; const isel = sel && sel.kind === "comp" && sel.id === b.id;
      const col = isel ? "#ffd54a" : i.keep ? (i.color || "#ffffff") : "#ff8a80";
      ctx.lineWidth = isel ? 2.5 : 1.2; ctx.setLineDash(isel ? [] : [4, 2]); ctx.strokeStyle = col;
      const pp = i.keep && polysOf ? polysOf(b.id, i.name) : null;
      if (pp && done.has(i.name) && !isel) continue;
      let at = pp ? strokePolys(ctx, pp, toPx) : null;
      if (pp) done.add(i.name);
      if (!at) { const [ax, ay] = toPx(b.bbox[0], b.bbox[1]), [bx, by] = toPx(b.bbox[2], b.bbox[3]); ctx.strokeRect(ax, ay, bx - ax, by - ay); at = [ax, ay]; }
      label(ctx, i.name + (i.keep ? "" : " (off)"), at[0] + 3, at[1] < 14 ? at[1] + 13 : at[1] - 4, col);
    }
    ctx.restore();
  }
  // the selected drawn shape: thick outline and a handle at every corner / vertex
  function paintSel(ctx, shapes, sel, toPx) {
    if (!sel || sel.kind !== "shape" || !shapes[sel.i]) return; const r = shapes[sel.i];
    ctx.save(); ctx.setLineDash([]); ctx.lineWidth = 2.5; ctx.strokeStyle = "#ffd54a";
    if (r.pts) { ringsPath(ctx, r, toPx); ctx.stroke(); }
    else { const [ax, ay] = toPx(Math.min(r.x0, r.x1), Math.min(r.y0, r.y1)), [bx, by] = toPx(Math.max(r.x0, r.x1), Math.max(r.y0, r.y1)); ctx.strokeRect(ax, ay, bx - ax, by - ay); }
    ctx.fillStyle = "#ffd54a"; ctx.lineWidth = 1; ctx.strokeStyle = "#000";
    for (const q of handles(r)) { const [x, y] = toPx(q[0], q[1]); ctx.fillRect(x - 3.5, y - 3.5, 7, 7); ctx.strokeRect(x - 3.5, y - 3.5, 7, 7); }
    ctx.restore();
  }
  // api: { tool() -> drawing tool or null, view() -> {x0,y0,s}, bounds() -> the slide [x0,y0,x1,y1] shapes are capped at, shapes(), boxes() -> [{id,bbox}], getSel(), setSel(sel), redraw(), commit(), active() }
  function attach(cv, api) {
    let drag = null;
    const data = (e) => { const b = cv.getBoundingClientRect(), v = api.view(); return [v.x0 + (e.clientX - b.left - 4) / v.s, v.y0 + (e.clientY - b.top - 4) / v.s]; };
    const tol = () => 7 / api.view().s;
    cv.addEventListener("mousedown", (e) => {
      if (api.tool()) return;
      const [x, y] = data(e), s = api.shapes(), sel = api.getSel();
      if (sel && sel.kind === "shape" && s[sel.i]) { const h = hitHandle(s[sel.i], x, y, tol()); if (h >= 0) { norm(s[sel.i]); drag = { mode: "handle", sh: s[sel.i], h, moved: false }; e.preventDefault(); return; } }
      const p = pick(s, api.boxes(), x, y, tol());
      api.setSel(p);
      if (p && p.kind === "shape") { norm(s[p.i]); drag = { mode: "move", sh: s[p.i], x, y, moved: false }; e.preventDefault(); }
    });
    cv.addEventListener("mousemove", (e) => {
      if (api.tool()) { cv.style.cursor = "crosshair"; return; }
      const [x, y] = data(e);
      if (drag) { if (drag.mode === "move") { move(drag.sh, x - drag.x, y - drag.y); drag.x = x; drag.y = y; } else dragHandle(drag.sh, drag.h, x, y); drag.moved = true; api.redraw(); return; }
      const s = api.shapes(), sel = api.getSel();
      cv.style.cursor = sel && sel.kind === "shape" && s[sel.i] && hitHandle(s[sel.i], x, y, tol()) >= 0 ? "grab" : pick(s, api.boxes(), x, y, tol()) ? "pointer" : "default";
    });
    window.addEventListener("mouseup", () => {
      if (!drag) return; const m = drag.moved, sh = drag.sh; norm(sh); drag = null; if (!m) return;
      if (api.bounds && !clip(sh, api.bounds())) { const s = api.shapes(), k = s.indexOf(sh); if (k >= 0) s.splice(k, 1); api.setSel(null); }   // capped at the slide
      api.commit();
    });
    window.addEventListener("keydown", (e) => {
      if ((e.key !== "Delete" && e.key !== "Backspace") || !api.active() || /INPUT|SELECT|TEXTAREA/.test((document.activeElement || {}).tagName || "")) return;
      const sel = api.getSel();
      if (sel && sel.kind === "shape") { e.preventDefault(); api.shapes().splice(sel.i, 1); api.setSel(null); api.commit(); }
      else if (sel && sel.kind === "comp" && api.offComp) { e.preventDefault(); api.offComp(sel.id); }
    });
  }
  return { attach, paintBoxes, paintSel, pick, label, strokePolys, clip, ringsPath };
})();
