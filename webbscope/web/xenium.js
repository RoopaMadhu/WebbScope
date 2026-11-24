"use strict";
const $ = (s) => document.querySelector(s);
const PAL = ["#4e9cf5", "#f5a142", "#5ec26a", "#e05d8b", "#a37ff0", "#3fc7c0", "#e6d44a", "#d9694a", "#8fb85c", "#c27dd1", "#7aa2c9", "#c9a27a"];
const X = { info: null, xy: null, spec: { gap: 200, min_cells: 500, min_counts: 1, comps: {}, rects: [] }, prev: null, draw: false, tok: 0, ann: null };
const post = (u, b) => fetch(u, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(b) });
const colorOf = (name) => { const i = X.prev.names.indexOf(name); return i < 0 ? "#555" : PAL[i % PAL.length]; };
const fmt = (n) => n.toLocaleString();

// ---- 1: scan
function files() { const o = {}; for (const k of ["matrix", "cells", "boundaries", "clusters", "umap", "pca", "experiment", "transcripts", "if_image"]) { const v = $("#f_" + k).value.trim(); if (v) o[k] = v; } return o; }
async function pollScan() {
  const s = await (await fetch("/api/xenium/status")).json();
  $("#log").textContent = s.log.join("\n"); $("#log").scrollTop = 1e9;
  if (s.state === "scanning") return setTimeout(pollScan, 700);
  $("#scan").disabled = false;
  if (s.state === "error") { $("#err").textContent = s.error; return; }
  if (s.state === "ready") loaded();
}
$("#scan").onclick = async () => {
  const path = $("#path").value.trim(), f = files();
  if (!path && !Object.keys(f).length) { $("#err").textContent = "Give the .tar.gz or .zip, the output folder, or the files."; return; }
  $("#err").textContent = ""; $("#scan").disabled = true; $("#log").textContent = "";
  for (const s of ["#s2", "#s3", "#s4"]) $(s).hidden = true;
  const r = await post("/api/xenium/scan", { path, files: f, transcripts: $("#wanttx").checked, ifimg: $("#wantif").checked, he_czi: $("#f_hecz").value.trim(), he_ann: $("#f_heann").value.trim() });
  if (!r.ok) { $("#err").textContent = await r.text(); $("#scan").disabled = false; return; }
  pollScan();
};
$("#path").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#scan").click(); });
// Browse buttons: the native window's file dialog when there is one, else the local server's (macOS) dialog; if neither works the path can still be typed.
async function pickPath(kind, types) {
  const api = window.pywebview && window.pywebview.api;
  try {
    if (api) return kind === "dir" ? await api.pick_folder() : await api.pick_file(types && types.length ? types : null);
    const r = await fetch(`/api/pick?kind=${kind}&types=${(types || []).join(",")}`);
    if (!r.ok) throw new Error(await r.text());
    return (await r.json()).path;
  } catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return null; }
}
const bind = (btn, target, kind, types) => { btn.addEventListener("click", async () => { const p = await pickPath(kind, types); if (p) { target.value = p; target.dispatchEvent(new Event("change")); } }); };
for (const b of document.querySelectorAll("button[data-pick-for]")) bind(b, $("#" + b.dataset.pickFor), b.dataset.kind, (b.dataset.types || "").split(",").filter(Boolean));
for (const inp of document.querySelectorAll("input[data-pick]")) {          // wrap the field and add a Browse… button after it
  const w = document.createElement("div"), b = document.createElement("button"); w.className = "pick"; b.type = "button"; b.textContent = "Browse…";
  inp.replaceWith(w); w.append(inp, b); bind(b, inp, inp.dataset.pick, (inp.dataset.types || "").split(",").filter(Boolean));
}

async function loaded() {
  X.info = await (await fetch("/api/xenium/info")).json();
  X.xy = new Float32Array(await (await fetch("/api/xenium/xy")).arrayBuffer());
  const i = X.info, p = i.panel;
  $("#runinfo").textContent = `${i.name} — ${fmt(i.n_cells)} cells, ${fmt(i.n_genes)} panel genes` + (p.panel_name ? ` (${p.panel_name}, ${p.panel_organism || ""} ${p.panel_tissue_type || ""})` : "")
    + `; ${fmt(i.zero_count_cells)} cells have no transcripts.` + (i.has_transcripts ? " Transcripts will be included." : " Transcripts are not included.");
  $("#lab10x").textContent = i.has_clusters ? `10x graph clusters (${i.n_clusters} clusters)` : "none found in this output (all cells get one label)";
  $("#pname").value = i.name; $("#outdir").value = i.default_output;
  $("#herow").hidden = !i.he; $("#hethrow").hidden = !i.he;
  if (i.he) { const c = i.he.czi, a = i.he.annotations; $("#hetxt").textContent = `Include the H&E (${c.file}: ${c.scenes.length} scene(s), ${c.px_um.toFixed(3)} µm/pixel` + (a ? `; ${a.n_regions} annotation polygons in ${a.names.filter((n) => n.n_regions).length} annotation names` : "") + "); building the tiles takes about a minute and about 0.5 GB"; }
  for (const s of ["#s2", "#s3", "#s4"]) $(s).hidden = false;
  X.sel = null; X.pending = null;
  X.spec = { gap: +$("#p_gap").value, min_cells: +$("#p_minc").value, min_counts: +$("#p_mint").value, detect: "cells", tx_min: 0, comps: {}, rects: [], cut_off: $("#cutoff").checked }; $("#p_detect").value = "cells"; $("#p_txmin").value = ""; txUi();
  sizeCanvas(); refresh();
}

// ---- 2: preview
function sizeCanvas() {
  const [x0, y0, x1, y1] = X.info.bounds, dx = x1 - x0, dy = y1 - y0, H = 640, W = Math.max(260, Math.min(620, Math.round(H * dx / dy)));
  const cv = $("#cv"); cv.width = W; cv.height = H; cv.style.width = W + "px"; cv.style.height = H + "px";
  X.view = { x0, y0, s: Math.min((W - 8) / dx, (H - 8) / dy) };
}
const toPx = (x, y) => [4 + (x - X.view.x0) * X.view.s, 4 + (y - X.view.y0) * X.view.s];
const toData = (px, py) => [X.view.x0 + (px - 4) / X.view.s, X.view.y0 + (py - 4) / X.view.s];
async function refresh() {
  const tok = ++X.tok;
  const r = await post("/api/xenium/preview", X.spec);
  if (!r.ok) { let m = ""; try { m = JSON.parse(await r.text()).detail; } catch (e) { /* plain text */ } $("#cnote").textContent = m || "preview failed"; if (X.spec.detect === "transcripts") { X.spec.detect = "cells"; $("#p_detect").value = "cells"; txUi(); $("#cnote").textContent += " (back to cell density)"; } return; }
  const d = await r.json(); if (tok !== X.tok) return;
  X.prev = d; draw(); table(); xSelBar();
  if (X.spec.detect === "transcripts") $("#txnote").textContent = `Transcript detection: a bin counts as tissue with at least ${d.tx_threshold ?? "?"} transcripts` + (X.spec.tx_min ? "." : " (chosen automatically).");
}
// transcript density: one tinted square per bin (brightness = log count, colour = the tissue of the bin, grey = no tissue), written into the pixels before the cells
function fillBins(u, W, H, grid, view, colorOf) {
  const bs = Math.max(1, Math.ceil(grid.bin * view.s)), q = grid.q, bg = [16, 18, 22];
  for (let i = 0; i < q.length; i += 4) {
    const col = q[i + 3] > 0 ? colorOf(q[i + 3]) : null, b = q[i + 2] / 255, a = col ? 0.15 + 0.6 * b : 0.05 + 0.3 * b, c = col || [138, 147, 166];
    const v = (255 << 24) | (((bg[2] + (c[2] - bg[2]) * a) | 0) << 16) | (((bg[1] + (c[1] - bg[1]) * a) | 0) << 8) | ((bg[0] + (c[0] - bg[0]) * a) | 0);
    const px0 = (4 + (grid.ox + q[i] * grid.bin - view.x0) * view.s) | 0, py0 = (4 + (grid.oy + q[i + 1] * grid.bin - view.y0) * view.s) | 0;
    for (let y = Math.max(0, py0); y < Math.min(H, py0 + bs); y++) for (let x = Math.max(0, px0); x < Math.min(W, px0 + bs); x++) u[y * W + x] = v;
  }
}
function hex(c) { return [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)]; }
function draw() {
  const cv = $("#cv"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height, img = ctx.createImageData(W, H), u = new Uint32Array(img.data.buffer);
  const pal = X.prev.names.map((_, i) => { const [r, g, b] = hex(PAL[i % PAL.length]); return (255 << 24) | (b << 16) | (g << 8) | r; }), grey = (255 << 24) | (85 << 16) | (85 << 8) | 85;
  const L = X.prev.labels, v = X.view;
  if (X.prev.grid) {                               // transcript density under the cells
    const rgb = {};
    for (const c of X.prev.comps) { const e = X.spec.comps[c.id] || {}, k = X.prev.names.indexOf(e.name ?? c.name); rgb[c.id] = (e.keep ?? c.keep) && k >= 0 ? hex(PAL[k % PAL.length]) : null; }
    fillBins(u, W, H, X.prev.grid, v, (id) => rgb[id]);
  }
  for (let i = 0; i < L.length; i++) {
    const l = L[i]; if (l === -2) continue;
    const px = (4 + (X.xy[2 * i] - v.x0) * v.s) | 0, py = (4 + (X.xy[2 * i + 1] - v.y0) * v.s) | 0;
    if (px >= 0 && px < W && py >= 0 && py < H) u[py * W + px] = l < 0 ? grey : pal[l % pal.length];
  }
  ctx.putImageData(img, 0, 0);
  ctx.lineWidth = 1.5; ctx.font = "11px system-ui"; ctx.setLineDash([5, 3]);
  for (const r of X.spec.rects) {
    ctx.strokeStyle = r.keep ? "#ffffff" : "#ff8a80"; ctx.fillStyle = ctx.strokeStyle;
    if (r.pts) {
      SE.ringsPath(ctx, r, toPx); ctx.stroke();
      const [x, y] = toPx(r.pts[0][0], r.pts[0][1]); SE.label(ctx, r.name + (r.keep ? "" : " (removed)"), x + 3, y < 14 ? y + 13 : y - 4, ctx.strokeStyle);
    } else {
      const [ax, ay] = toPx(Math.min(r.x0, r.x1), Math.min(r.y0, r.y1)), [bx, by] = toPx(Math.max(r.x0, r.x1), Math.max(r.y0, r.y1));
      ctx.strokeRect(ax, ay, bx - ax, by - ay); SE.label(ctx, r.name + (r.keep ? "" : " (removed)"), ax + 3, ay + 13, ctx.strokeStyle);
    }
  }
  SE.paintBoxes(ctx, xBoxes(), toPx, xInfo, X.sel, (id, name) => polysByName(X.prev, name)); SE.paintSel(ctx, X.spec.rects, X.sel, toPx);
  ctx.setLineDash([]);
  $("#cnote").textContent = `${fmt(X.prev.n_kept)} cells in ${X.prev.names.length} tissue(s)` + (X.prev.n_below_min ? `; ${fmt(X.prev.n_below_min)} below the transcript minimum are hidden` : "") + "; grey = switched off";
}
function table() {
  const tb = $("#ttab tbody"); tb.innerHTML = "";
  const row = (cls, keep, name, n, onKeep, onName, onDel, selected, onSel) => {
    const tr = document.createElement("tr"); tr.className = cls + (selected ? " selrow" : ""); tr.addEventListener("click", (e) => { if (!e.target.closest("input, .x")) onSel(); });
    tr.innerHTML = `<td></td><td><input type="checkbox" ${keep ? "checked" : ""}></td><td><span class="sw" style="background:${colorOf(name)}"></span> <input type="text" value="${name.replace(/"/g, "&quot;")}"></td><td>${n}</td><td>${onDel ? '<span class="x" title="remove">✕</span>' : ""}</td>`;
    tr.querySelector("input[type=checkbox]").onchange = (e) => onKeep(e.target.checked);
    tr.querySelector("input[type=text]").addEventListener("change", (e) => onName(e.target.value));
    if (onDel) tr.querySelector(".x").onclick = () => { X.sel = null; onDel(); };
    tb.appendChild(tr);
  };
  for (const c of X.prev.comps) {
    const e = X.spec.comps[c.id] || {}, name = e.name ?? c.name, keep = e.keep ?? c.keep;
    row("", keep, name, fmt(c.n), (k) => { (X.spec.comps[c.id] ||= {}).keep = k; refresh(); }, (nm) => { (X.spec.comps[c.id] ||= {}).name = nm; refresh(); }, undefined, X.sel && X.sel.kind === "comp" && X.sel.id === c.id, () => xSetSel({ kind: "comp", id: c.id }));
  }
  X.spec.rects.forEach((r, i) => row("rect", r.keep, r.name, r.pts ? "(free-hand)" : "(rectangle)", (k) => { r.keep = k; refresh(); }, (nm) => { r.name = nm; refresh(); }, () => { X.spec.rects.splice(i, 1); refresh(); }, X.sel && X.sel.kind === "shape" && X.sel.i === i, () => xSetSel({ kind: "shape", i })));
  $("#tsum").textContent = "Project tissues (samples): " + (Object.entries(X.prev.counts).map(([n, c]) => `${n} (${fmt(c)})`).join(", ") || "none — switch at least one on");
}
// ---- tissue boxes from a GeoJSON an earlier import exported: they become the shapes that separate the tissues, the detected tissues start switched off.
// (Only for telling tissues apart on the slide; anatomical regions within a tissue are a different input.)
$("#tisloadgo").onclick = async () => {
  const path = $("#tisload").value.trim(), msg = $("#tismsg"); if (!path) { msg.textContent = "Give the path of a .geojson file."; return; }
  const r = await post("/api/xenium/tissues/load", { path, name_col: "" });
  if (!r.ok) { let m = "could not load the file"; try { m = JSON.parse(await r.text()).detail; } catch (e) { /* plain text */ } msg.textContent = m; return; }
  const d = await r.json(), sb = slideB(), shapes = d.shapes.filter((q) => SE.clip(q, sb));
  X.spec.rects = shapes; X.spec.comps = {}; X.spec.auto_off = true; X.sel = null; X.pending = null; xSelBar();
  msg.textContent = `Loaded ${d.n_tissues} tissue box(es)` + (d.n_removed ? ` and ${d.n_removed} switched-off shape(s)` : "") + ` from ${d.source.split("/").pop()} (names from “${d.name_col}”)` + (shapes.length < d.shapes.length ? `; ${d.shapes.length - shapes.length} shape(s) off the slide were dropped` : "") + ". The detected tissues are switched off; tick one to use it again.";
  refresh();
};
$("#cutoff").addEventListener("change", () => { X.spec.cut_off = $("#cutoff").checked; refresh(); });
// ---- export the tissue boxes now (as edited), without converting
async function tisExport(overwrite) {
  const msg = $("#tisoutmsg"), folder = $("#tisout").value.trim(); $("#tisover").hidden = true;
  if (!folder) { msg.textContent = "Choose a folder for the file."; return; }
  if (!X.prev || !X.prev.tissue_polys_raw) { msg.textContent = "No tissue boxes yet: read the run first."; return; }
  const name = $("#pname").value.trim() || "slide";     // from the outlines this page holds, as on the batch page: works even if the server has lost the run
  const r = await post("/api/xenium/batch/tissues/export", { runs: [{ name, shapes: X.spec.rects, tissues: X.prev.tissue_polys_raw }], folder, cut_off: $("#cutoff").checked, name_col: $("#regcol").value.trim() || "sample_id_unblinded", overwrite });
  let d = null; try { d = JSON.parse(await r.text()); } catch (e) { /* not JSON */ }
  if (r.status === 409) { msg.textContent = (d && d.detail) + ". Overwrite it?"; $("#tisover").hidden = false; return; }
  if (!r.ok) { msg.textContent = (d && d.detail) || "export failed"; return; }
  if (d.errors[name]) { msg.textContent = d.errors[name]; return; }
  const ov = d.overlaps[name] || [], em = d.empty[name] || [];
  msg.textContent = "Written: " + d.written.join(", ") + (ov.length ? `. Overlapping geometries found (written as they are): ${ov.map((q) => q.join(" & ")).join(", ")}` : "") + (em.length ? `. Left out because nothing of them remained after the cut: ${em.join(", ")}` : "");
}
$("#tisexport").onclick = () => tisExport(false); $("#tisover").onclick = () => tisExport(true);
function txUi() { document.querySelectorAll(".txonly").forEach((e) => (e.hidden = X.spec.detect !== "transcripts")); if (X.spec.detect !== "transcripts") $("#txnote").textContent = ""; }
for (const id of ["p_gap", "p_minc", "p_mint", "p_detect", "p_txmin"]) $("#" + id).addEventListener("change", () => {
  X.spec.gap = +$("#p_gap").value; X.spec.min_cells = +$("#p_minc").value; X.spec.min_counts = +$("#p_mint").value;
  X.spec.detect = $("#p_detect").value; X.spec.tx_min = Math.max(0, Math.round(+$("#p_txmin").value || 0)); X.spec.comps = {}; txUi(); refresh();
});
// ---- selecting a box on the slide: drag a drawn shape, drag its handles, Delete removes it; a detected tissue box can be made editable
const slideB = () => (X.prev && X.prev.slide_bounds) || (X.info && X.info.bounds) || null;   // the slide: shapes are capped at it
const polysByName = (prev, name) => { const t = prev && prev.tissue_polys && prev.tissue_polys.find((q) => q.name === name); return t ? t.polys : null; };
const xBoxes = () => (X.prev && X.prev.boxes ? Object.entries(X.prev.boxes).map(([id, b]) => ({ id: +id, bbox: b })) : []);
function xInfo(id) { const c = X.prev && X.prev.comps.find((q) => q.id === id); if (!c) return null; const e = X.spec.comps[id] || {}; const name = e.name ?? c.name, k = X.prev.names.indexOf(name); return { name, keep: e.keep ?? c.keep, color: k >= 0 ? PAL[k % PAL.length] : "#ffffff" }; }
function xSetSel(p) { X.sel = p; if (X.prev) { draw(); table(); } xSelBar(); }
function xSelBar() {
  const bar = $("#selbar"), sel = X.sel, btn = (t, f) => { const b = document.createElement("button"); b.textContent = t; b.onclick = f; return b; }; bar.replaceChildren();
  if (!sel) { bar.textContent = "Click a box on the slide to select it. A drawn shape can then be dragged, reshaped by its handles, or removed with Delete; a detected tissue box can be made editable."; return; }
  if (sel.kind === "shape") { const r = X.spec.rects[sel.i]; if (!r) return; bar.append(`Selected: ${r.name} (${r.pts ? "free-hand" : "rectangle"}): drag it to move, drag a handle to reshape.`, btn("Remove", () => { X.spec.rects.splice(sel.i, 1); X.sel = null; refresh(); })); return; }
  const i = xInfo(sel.id); if (!i) return;
  bar.append(`Selected detected tissue ${i.name}${i.keep ? "" : " (off)"}.`,
    btn("Make editable", () => {          // a rectangle copy of the box takes its place: the tissue is switched off, the shape keeps its name
      const b = X.prev.boxes[sel.id]; (X.spec.comps[sel.id] ||= {}).keep = false; X.spec.rects.push({ x0: b[0], y0: b[1], x1: b[2], y1: b[3], name: i.name, keep: true }); X.sel = { kind: "shape", i: X.spec.rects.length - 1 }; refresh();
    }),
    btn("Replace with rectangle", () => xReplace(sel.id, "rect")), btn("Replace with free-hand", () => xReplace(sel.id, "free")),
    btn(i.keep ? "Switch off" : "Switch on", () => { (X.spec.comps[sel.id] ||= {}).keep = !i.keep; refresh(); }));
}
// replace a detected tissue by a new shape that takes its name: the tissue is switched off, the next shape drawn is named like it (cancelling the tool puts the tissue back)
function xReplace(id, tool) {
  const i = xInfo(id), e = (X.spec.comps[id] ||= {}); X.pending = { id, name: i.name, prevKeep: e.keep };
  e.keep = false; X.sel = null; X.draw = null; setTool(tool); refresh();
}
function xCancelPending() { const p = X.pending; X.pending = null; if (!p) return; const e = (X.spec.comps[p.id] ||= {}); if (p.prevKeep === undefined) delete e.keep; else e.keep = p.prevKeep; refresh(); }
function xNewName() { if (X.pending) { const nm = X.pending.name; X.pending = null; return nm; } return "region_" + (X.spec.rects.length + 1); }
SE.attach($("#cv"), { tool: () => X.draw, view: () => X.view, bounds: () => slideB(), shapes: () => X.spec.rects, boxes: xBoxes, getSel: () => X.sel, setSel: xSetSel, redraw: draw, commit: () => refresh(), active: () => !$("#s2").hidden, offComp: (id) => { (X.spec.comps[id] ||= {}).keep = false; refresh(); } });
function setTool(t) {
  X.draw = X.draw === t ? null : t; if (X.draw) X.sel = null; else if (X.pending) xCancelPending();
  $("#draw").textContent = X.draw === "rect" ? "Rectangle: drag on the slide (click to cancel)" : "Rectangle";
  $("#drawfree").textContent = X.draw === "free" ? "Free-hand: trace on the slide (click to cancel)" : "Free-hand";
}
$("#draw").onclick = () => setTool("rect"); $("#drawfree").onclick = () => setTool("free");
function simplify(pts, tol) {                    // Douglas-Peucker on the screen path
  if (pts.length < 3) return pts;
  const [ax, ay] = pts[0], [bx, by] = pts[pts.length - 1], dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1;
  let m = -1, k = 0; for (let i = 1; i < pts.length - 1; i++) { const d = Math.abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / L; if (d > m) { m = d; k = i; } }
  if (m <= tol) return [pts[0], pts[pts.length - 1]];
  return simplify(pts.slice(0, k + 1), tol).slice(0, -1).concat(simplify(pts.slice(k), tol));
}
(function () {
  const cv = $("#cv"); let a = null, path = null;
  const pos = (e) => { const b = cv.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
  cv.addEventListener("mousedown", (e) => { if (X.draw === "rect") a = pos(e); else if (X.draw === "free") path = [pos(e)]; });
  cv.addEventListener("mousemove", (e) => {
    const ctx = cv.getContext("2d"), q = pos(e);
    if (a) { draw(); ctx.strokeStyle = "#ffd54a"; ctx.setLineDash([4, 3]); ctx.strokeRect(a[0], a[1], q[0] - a[0], q[1] - a[1]); ctx.setLineDash([]); }
    else if (path) {
      const l = path[path.length - 1]; if (Math.hypot(q[0] - l[0], q[1] - l[1]) < 2) return; path.push(q);
      draw(); ctx.strokeStyle = "#ffd54a"; ctx.beginPath(); path.forEach((p, k) => k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.stroke();
    }
  });
  window.addEventListener("mouseup", (e) => {
    const done = () => { setTool(X.draw); refresh(); };
    if (a) {
      const b = pos(e), p0 = toData(a[0], a[1]), p1 = toData(b[0], b[1]); a = null;
      if (Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) > 50) { { const r = { x0: p0[0], y0: p0[1], x1: p1[0], y1: p1[1], name: xNewName(), keep: true }; if (SE.clip(r, slideB())) { X.spec.rects.push(r); done(); } else draw(); } } else draw();
    } else if (path) {
      let far = 0, fd = -1; path.forEach((q, i) => { const d = Math.hypot(q[0] - path[0][0], q[1] - path[0][1]); if (d > fd) { fd = d; far = i; } });   // the path is a closed loop: split it at the point farthest from the start
      const pts = (fd < 3 ? [] : simplify(path.slice(0, far + 1), 1.2).slice(0, -1).concat(simplify(path.slice(far).concat([path[0]]), 1.2).slice(0, -1))).map((p) => toData(p[0], p[1])); path = null;
      { const cr = { pts }; if (pts.length >= 3 && SE.clip(cr, slideB())) pts.splice(0, pts.length, ...cr.pts); else pts.length = 0; }   // capped at the slide
      let ar = 0; pts.forEach((p, i) => { const q = pts[(i + 1) % pts.length]; ar += p[0] * q[1] - q[0] * p[1]; });
      if (pts.length >= 3 && Math.abs(ar) / 2 > 2500) { X.spec.rects.push({ pts, name: xNewName(), keep: true }); done(); } else draw();
    }
  });
})();

// ---- 3: labels
document.querySelectorAll("input[name=lab]").forEach((r) => r.addEventListener("change", () => { $("#annbox").hidden = document.querySelector("input[name=lab]:checked").value !== "ann"; }));
function fill(sel, cols, blank) { sel.innerHTML = (blank ? '<option value="">(none)</option>' : "") + cols.map((c) => `<option>${c}</option>`).join(""); }
$("#annload").onclick = async () => {
  $("#annmsg").textContent = "";
  const r = await post("/api/xenium/annotation", { path: $("#annpath").value.trim() });
  if (!r.ok) { X.ann = null; $("#anncols").hidden = true; $("#annmsg").textContent = JSON.parse(await r.text()).detail || "cannot read the file"; return; }
  const d = await r.json(); X.ann = d; $("#anncols").hidden = false;
  fill($("#a_id"), d.columns); fill($("#a_type"), d.columns, true); fill($("#a_state"), d.columns, true);
  $("#a_id").value = d.id_col; $("#a_type").value = d.type_col || ""; $("#a_state").value = d.state_col || "";
  const upd = () => { const k = $("#a_id").value, m = d.match[k] || 0; $("#annmsg").textContent = `${fmt(d.n_rows)} rows; ${(100 * m).toFixed(1)}% of them match this run's cell ids; ${fmt(d.cells_covered[k] || 0)} of the run's ${fmt(d.n_cells)} cells get a label (the rest become \"unlabelled\").` + ($("#a_type").value ? "" : " Pick the cell type column."); };
  $("#a_id").onchange = upd; $("#a_type").onchange = upd; upd();
};

// ---- 4: convert
async function pollConv() {
  const s = await (await fetch("/api/status")).json();
  $("#log2").textContent = s.log.join("\n"); $("#log2").scrollTop = 1e9;
  if (s.state === "importing") return setTimeout(pollConv, 600);
  $("#go").disabled = false;
  if (s.state === "error") { $("#err2").textContent = s.error; return; }
  if (s.state === "ready") location.href = "/";
}
$("#go").onclick = async () => {
  $("#err2").textContent = "";
  let ann = null;
  if (document.querySelector("input[name=lab]:checked").value === "ann") {
    if (!X.ann || !$("#a_type").value) { $("#err2").textContent = "Load the annotation file and pick its cell type column (or use the 10x clusters)."; return; }
    ann = { path: $("#annpath").value.trim(), id_col: $("#a_id").value, type_col: $("#a_type").value, state_col: $("#a_state").value || null };
  }
  if (!Object.keys(X.prev.counts).length) { $("#err2").textContent = "No tissue is switched on."; return; }
  $("#go").disabled = true; $("#log2").textContent = "";
  const r = await post("/api/xenium/convert", { spec: X.spec, output: $("#outdir").value.trim(), name: $("#pname").value.trim(), annotation: ann || {}, transcripts: !!(X.info.has_transcripts && $("#wanttx").checked), he: !!(X.info.he && $("#wanthe").checked), he_threads: Math.max(1, Math.round(+$("#hethreads").value || 4)),
    he_reg_folder: $("#f_hereg").value.trim(), he_reg_index: $("#f_heregidx").value.trim(), imc_file: $("#f_imc").value.trim(), export_regions: $("#wantreg").checked, regions_name_col: $("#regcol").value.trim() || "sample_id_unblinded", regions_dir: $("#regdir").value.trim() });
  if (!r.ok) { $("#err2").textContent = (JSON.parse(await r.text()).detail) || "failed"; $("#go").disabled = false; return; }
  pollConv();
};
(async () => { const s = await (await fetch("/api/xenium/status")).json(); if (s.state === "scanning") { $("#scan").disabled = true; pollScan(); } else if (s.state === "ready") { $("#log").textContent = s.log.join("\n"); loaded(); } })();
