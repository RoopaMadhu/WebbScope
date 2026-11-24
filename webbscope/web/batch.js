"use strict";
const $ = (s) => document.querySelector(s);
const post = (u, b) => fetch(u, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(b) });
let REG = {}, COMPS = {}, THUMB = {}, EDITPOLY = {}, EDITRAW = {}, AUTOOFF = {}, IMPORTED = {};   // EDITRAW: the editor's tissue outlines before the switched-off shapes are cut out (for the export)
// AUTOOFF / IMPORTED: runs whose tissue boxes were loaded from a GeoJSON -> true / [{name, polys}]
//   // run name -> drawn shapes; -> {piece id: {name, keep}} edits of the detected tissues; -> thumbnail data of the gallery
let RUNS = [], PARENT_MERGED = null, HE = null, ANN = { cols: null, values: null, match: null };   // HE = {czi: [{name, path}], annotations: [...], match: {run: {czi, annotations, why}}} once an H&E folder was matched
async function detail(r) { try { return JSON.parse(await r.text()).detail; } catch (e) { return "failed"; } }
$("#find").onclick = async () => {
  $("#err").textContent = "";
  const r = await post("/api/xenium/batch/scan", { path: $("#path").value.trim() });
  if (!r.ok) { $("#err").textContent = await detail(r); return; }
  const d = await r.json();
  RUNS = d.runs; PARENT_MERGED = d.merged || null;
  if (!RUNS.length) { $("#err").textContent = "No Xenium runs found directly inside this folder (looked for Xenium output sub-folders with cell_feature_matrix/, cells and cell_boundaries, .tar.gz or .zip bundles, and project folders already written by this app)." + (d.skipped.length ? " Skipped sub-folders: " + d.skipped.join(", ") : ""); return; }
  $("#runsum").textContent = `${RUNS.length} run(s)` + (d.skipped.length ? `; ${d.skipped.length} sub-folder(s) skipped because they are not Xenium outputs: ${d.skipped.join(", ")}` : "");
  REG = {}; COMPS = {}; THUMB = {}; EDITPOLY = {}; EDITRAW = {}; AUTOOFF = {}; IMPORTED = {}; HE = null; $("#tismsg").textContent = ""; $("#hemsg").textContent = ""; document.querySelectorAll(".hecol").forEach((e) => (e.hidden = true));
  $("#ttab tbody").innerHTML = RUNS.map((x, i) => `<tr><td><input type="checkbox" data-i="${i}" checked></td><td></td><td>${x.kind}</td><td>${x.kind === "project" ? '<span class="sub">n/a</span>' : `<button data-d="${i}">Draw…</button> <span class="sub" data-n="${i}"></span>`}</td><td class="hecol" hidden></td><td class="hecol" hidden></td><td class="anncol" hidden></td></tr>`).join(""); ANN.match = null; annShow();
  document.querySelectorAll("#ttab tbody tr").forEach((tr, i) => { tr.children[1].textContent = RUNS[i].name; });
  $("#outdir").value = d.default_output;
  $("#s2").hidden = false; $("#s3").hidden = false; checkMerged(); tisState();
};
// An earlier merged project (<folder>/_merged) is offered for opening as it is: nothing is converted or merged again.
async function checkMerged() {
  const found = PARENT_MERGED ? [PARENT_MERGED] : [], out = $("#outdir").value.trim(), box = $("#mergedbox"); box.replaceChildren();
  if (out) { try { const r = await fetch("/api/xenium/batch/merged?folder=" + encodeURIComponent(out)); const m = r.ok ? await r.json() : null; if (m && !found.some((f) => f.path === m.path)) found.push(m); } catch (e) { /* none */ } }
  const ticked = new Set([...document.querySelectorAll("#ttab input[type=checkbox]:checked")].map((c) => RUNS[+c.dataset.i].name));
  for (const m of found) {
    const same = m.runs.length === ticked.size && m.runs.every((r) => ticked.has(r)), d = document.createElement("div"), t = document.createElement("div"), b = document.createElement("button");
    d.style.cssText = "border:1px solid #29303b;border-radius:6px;padding:8px 10px;margin:6px 0";
    t.textContent = `An earlier merged project exists: ${m.path} — ${m.runs.length} run(s), ${m.n_tissues} tissue(s)` + (m.n_cells ? `, ${m.n_cells.toLocaleString()} cells` : "") + (m.transcripts ? ", with transcripts" : "") + `. It lists ${same ? "the same runs as the ticked ones" : "runs that differ from the ticked ones"}.`;
    b.textContent = "Open it (no merging)"; b.className = "primary"; b.style.marginTop = "6px";
    b.onclick = async () => {
      $("#err2").textContent = ""; b.disabled = true; $("#log2").textContent = "";
      const r = await post("/api/open", { path: m.path });
      if (!r.ok) { $("#err2").textContent = await r.text(); b.disabled = false; return; }
      pollConv();
    };
    d.append(t, b); box.append(d);
  }
}
$("#outdir").addEventListener("change", checkMerged);
$("#ttab").addEventListener("change", (e) => { if (e.target.type === "checkbox") { checkMerged(); tisState(); } });
// ---- tissue boxes from the <run>_regions.geojson files of an earlier import (only to tell tissues apart, not anatomical regions): each found file's shapes replace the run's shapes, its detected tissues start switched off
$("#tisfind").onclick = async () => {
  const msg = $("#tismsg"), folder = $("#tisfolder").value.trim(); msg.textContent = "";
  if (!folder) { msg.textContent = "Enter the folder with the <run>_regions.geojson files."; return; }
  const runs = RUNS.filter((x) => x.kind !== "project").map((x) => x.name);
  const r = await post("/api/xenium/batch/tissues", { folder, runs, name_col: "" });
  if (!r.ok) { msg.textContent = await detail(r); return; }
  const d = await r.json(), names = Object.keys(d.found);
  for (const n of names) { REG[n] = d.found[n].shapes; COMPS[n] = {}; AUTOOFF[n] = true; IMPORTED[n] = d.found[n].tissues; delete EDITPOLY[n]; delete EDITRAW[n]; }
  const errs = Object.entries(d.errors).map(([n, e]) => `${n}: ${e}`);
  msg.textContent = `Loaded tissue boxes for ${names.length} of ${runs.length} run(s)` + (d.missing.length ? `; no file for ${d.missing.length}: ${d.missing.join(", ")} (they keep their detected tissues)` : "") + (errs.length ? `; could not read ${errs.length}: ${errs.join("; ")}` : "") + ".";
  refreshNotes(); tisState(); if (!$("#gal").hidden) galRender();
};
// ---- export the tissue boxes of the ticked runs now (as edited here), from the geometry this page holds: no run is read
// A run's tissue outlines before the switched-off shapes are cut out: null = its kept shapes are the tissues (tissue boxes loaded from a GeoJSON, no
// detected tissue switched back on); else [{name, polys}] as the editor or the slide review built them; undefined = not built yet.
const hasCompEdits = (n) => Object.values(COMPS[n] || {}).some((e) => e.name || e.keep !== undefined);
function builtTissues(n) {
  if (AUTOOFF[n] && !Object.values(COMPS[n] || {}).some((e) => e.keep === true)) return null;
  if (EDITRAW[n] && EDITRAW[n].key === skey()) return EDITRAW[n].polys;
  const t = THUMB[n];
  if (!(REG[n] || []).length && !hasCompEdits(n) && t && t.d && t.key === skey()) return t.d.comps.filter((c) => c.id > 0 && c.keep && t.d.polys && t.d.polys[c.id]).map((c) => ({ name: c.name, polys: t.d.polys[c.id] }));
  return undefined;
}
// the export is greyed out until every ticked run has its tissue boxes (loaded from GeoJSON, or built in the slide review / the editor)
function tisState() {
  const runs = galRuns(), miss = runs.filter((x) => builtTissues(x.name) === undefined), ok = runs.length > 0 && !miss.length;
  $("#tisexport").disabled = !ok; $("#tisout").disabled = !ok; if (!ok) $("#tisover").hidden = true;
  $("#tisready").textContent = !runs.length ? "No run is ticked." : !miss.length ? "" : `Tissue boxes are not built yet for ${miss.length} of ${runs.length} ticked run(s): `
    + miss.slice(0, 6).map((x) => x.name).join(", ") + (miss.length > 6 ? ", …" : "") + ". Load them from GeoJSON, or build them in Review all slides… (every page; a run edited before the detection settings changed has to be opened again).";
}
async function tisExport(overwrite) {
  const msg = $("#tisoutmsg"), folder = $("#tisout").value.trim(); $("#tisover").hidden = true;
  if (!folder) { msg.textContent = "Choose a folder for the files."; return; }
  const runs = galRuns(); if (!runs.length) { msg.textContent = "No run is ticked."; return; }
  const body = runs.map((x) => ({ name: x.name, shapes: REG[x.name] || [], tissues: builtTissues(x.name) }));
  if (body.some((x) => x.tissues === undefined)) { tisState(); return; }
  msg.textContent = `writing ${runs.length} file(s)…`; $("#tisexport").disabled = true;
  const r = await post("/api/xenium/batch/tissues/export", { runs: body, folder, cut_off: $("#cutoff").checked, overwrite, name_col: $("#regcol").value.trim() || "sample_id_unblinded" });
  tisState();
  let d = null; try { d = JSON.parse(await r.text()); } catch (e) { /* not JSON */ }
  if (r.status === 409 && !overwrite) { msg.textContent = (d && d.detail) + ". Overwrite them?"; $("#tisover").hidden = false; return; }
  if (!r.ok) { msg.textContent = (d && d.detail) || "export failed"; return; }
  const errs = Object.entries(d.errors).map(([n, e]) => `${n}: ${e}`), ov = Object.entries(d.overlaps || {}), em = Object.entries(d.empty || {});
  msg.textContent = `Wrote ${d.written.length} file(s) to ${folder}` + (errs.length ? `; failed for ${errs.length}: ${errs.join("; ")}` : "") + "."
    + (ov.length ? ` Overlapping geometries found (written as they are) in ${ov.length} run(s): ` + ov.map(([n, ps]) => `${n}: ${ps.map((p) => p.join(" & ")).join(", ")}`).join("; ") + "." : "")
    + (em.length ? ` Left out because nothing of them remained after the cut: ` + em.map(([n, ts]) => `${n}: ${ts.join(", ")}`).join("; ") + "." : "");
}
$("#tisexport").onclick = () => tisExport(false); $("#tisover").onclick = () => tisExport(true);
$("#cutoff").addEventListener("change", () => { EDITPOLY = {}; if (!$("#gal").hidden) galRender(); });
// ---- H&E folder: suggested scan (and annotation file) per run, editable in the table
function heSelect(items, chosen, cls) {
  const sel = document.createElement("select"); sel.className = cls; sel.style.maxWidth = "260px";
  const o0 = document.createElement("option"); o0.value = ""; o0.textContent = "(none)"; sel.append(o0);
  for (const it of items) { const o = document.createElement("option"); o.value = it.path; o.textContent = it.name; o.title = it.path; sel.append(o); }
  sel.value = chosen || ""; return sel;
}
$("#hefind").onclick = async () => {
  const msg = $("#hemsg"); msg.textContent = "";
  const folder = $("#hepath").value.trim(); if (!folder) { msg.textContent = "Enter the folder with the H&E scans."; return; }
  const r = await post("/api/xenium/batch/he", { folder, runs: RUNS.map((x) => x.name) });
  if (!r.ok) { msg.textContent = await detail(r); return; }
  HE = await r.json();
  if (!HE.czi.length) { msg.textContent = "No H&E scans (.czi, .ome.tif, .ome.tiff) found in that folder (looked two levels down)."; HE = null; return; }
  document.querySelectorAll(".hecol").forEach((e) => (e.hidden = false));
  const rows = document.querySelectorAll("#ttab tbody tr"); let n = 0;
  RUNS.forEach((run, i) => {
    const m = HE.match[run.name] || {}, tr = rows[i]; if (m.czi) n++;
    const a = heSelect(HE.czi, m.czi, "hesel"), b = heSelect(HE.annotations, m.annotations, "hean"), why = document.createElement("div");
    why.className = "sub"; why.style.cssText = "font-size:11px;margin:2px 0 0"; why.textContent = m.why || "";
    a.addEventListener("change", () => { why.textContent = a.value ? "chosen by you" : ""; });
    tr.children[4].replaceChildren(a, why); tr.children[5].replaceChildren(b);
  });
  msg.textContent = `${HE.czi.length} scan(s) and ${HE.annotations.length} annotation file(s) found; ${n} of ${RUNS.length} run(s) matched by name. Runs left at (none) get no H&E.`;
};
// ---- cell labels from the user's annotation file (batch.ann_inspect / ann_match); applied when the runs are merged
const ID_PREFIX = "__id_prefix__";
const annEsc = (v) => String(v).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/"/g, "&quot;");
function annFill(sel, cols, blank, extra) { sel.innerHTML = (blank ? '<option value="">(none)</option>' : "") + (extra || []).map(([v, t]) => `<option value="${annEsc(v)}">${annEsc(t)}</option>`).join("") + cols.map((c) => `<option value="${annEsc(c)}">${annEsc(c)}</option>`).join(""); }
function annShow() {   // the "cell labels" column of the runs table
  const on = !!ANN.cols, byRun = on && !!$("#a_run").value;
  document.querySelectorAll(".anncol").forEach((e) => (e.hidden = !on));
  if (!on) return;
  document.querySelectorAll("#ttab tbody tr").forEach((tr, i) => {
    const td = tr.children[6], run = RUNS[i]; if (!td || !run) return;
    const note = (t) => td.replaceChildren(Object.assign(document.createElement("span"), { className: "sub", textContent: t }));
    if (!byRun) return note(RUNS.length > 1 ? "choose the run column" : "all rows (one run)");
    if (!ANN.match) return note("press Match to runs");
    const m = ANN.match[run.name] || {}, sel = document.createElement("select"); sel.className = "annsel"; sel.style.maxWidth = "220px";
    sel.innerHTML = '<option value="">(none: 10x clusters)</option>' + Object.entries(ANN.values).map(([v, n]) => `<option value="${annEsc(v)}">${annEsc(v)} (${n.toLocaleString()} cells)</option>`).join("");
    sel.value = m.value || ""; const why = document.createElement("div"); why.className = "sub"; why.style.cssText = "font-size:11px;margin:2px 0 0"; why.textContent = m.why || "";
    sel.addEventListener("change", () => { why.textContent = sel.value ? "chosen by you" : ""; }); td.replaceChildren(sel, why);
  });
}
$("#annload").onclick = async () => {
  const msg = $("#annmsg"); ANN.cols = null; ANN.match = null; $("#anncols").hidden = $("#annmatchrow").hidden = true;
  const path = $("#annpath").value.trim(); if (!path) { msg.textContent = "Enter the path of the annotation file."; annShow(); return; }
  msg.textContent = "reading…";
  const r = await post("/api/xenium/batch/ann", { path, runs: RUNS.map((x) => x.name) });
  if (!r.ok) { msg.textContent = await detail(r); annShow(); return; }
  const d = await r.json(); ANN.cols = d;
  annFill($("#a_id"), d.columns); annFill($("#a_run"), d.columns, true, d.id_prefix ? [[ID_PREFIX, "the prefix of the cell id"]] : []);
  annFill($("#a_type"), d.columns, true); annFill($("#a_state"), d.columns, true); annFill($("#a_donor"), d.columns, true);
  $("#a_id").value = d.id_col; $("#a_run").value = d.run_col || ""; $("#a_type").value = d.type_col || ""; $("#a_state").value = d.state_col || ""; $("#a_donor").value = d.donor_col || "";
  $("#anncols").hidden = false;
  const upd = () => { $("#annmatchrow").hidden = !$("#a_run").value; ANN.match = null; annShow(); }; $("#a_run").onchange = upd; $("#a_id").onchange = upd; upd();
  msg.textContent = `${d.n_rows.toLocaleString()} rows. Check the columns` + (d.run_col ? "." : RUNS.length > 1 ? "; with several runs choose the run column (Xenium cell ids repeat between runs)." : ".");
  if (d.run_col) await $("#annmatch").onclick();
};
$("#annmatch").onclick = async () => {
  const msg = $("#annmsg"); msg.textContent = "matching…";
  const r = await post("/api/xenium/batch/ann/match", { path: $("#annpath").value.trim(), runs: RUNS.map((x) => x.name), run_col: $("#a_run").value, id_col: $("#a_id").value });
  if (!r.ok) { msg.textContent = await detail(r); ANN.match = null; annShow(); return; }
  const d = await r.json(); ANN.values = d.values; ANN.match = d.match; annShow();
  const n = Object.values(d.match).filter((m) => m.value).length;
  msg.textContent = `${Object.keys(d.values).length} value(s) in the run column; ${n} of ${RUNS.length} run(s) matched. Check the cell labels column of the runs table; runs left at (none) keep their 10x clusters.`;
};
function annChoice() {   // {path, columns, runs: {run: value}} for the ticked runs, or {} when no file is loaded
  if (!ANN.cols || !$("#annpath").value.trim()) return {};
  const runs = {}, byRun = !!$("#a_run").value;
  document.querySelectorAll("#ttab tbody tr").forEach((tr, i) => { if (!tr.querySelector("input[type=checkbox]").checked) return; const s = tr.querySelector(".annsel"); if (byRun && s && s.value) runs[RUNS[i].name] = s.value; });
  return { path: $("#annpath").value.trim(), id_col: $("#a_id").value, run_col: $("#a_run").value, type_col: $("#a_type").value, state_col: $("#a_state").value, donor_col: $("#a_donor").value, runs };
}
function heChoice() {   // {run: {czi, annotations}} for the ticked runs that have a scan chosen
  const out = {}; if (!HE) return out;
  document.querySelectorAll("#ttab tbody tr").forEach((tr, i) => {
    const c = tr.querySelector(".hesel"); if (c && c.value && tr.querySelector("input[type=checkbox]").checked) out[RUNS[i].name] = { czi: c.value, annotations: tr.querySelector(".hean").value || "" };
  });
  return out;
}
$("#path").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#find").click(); });
async function pollConv() {
  const s = await (await fetch("/api/status")).json();
  $("#log2").textContent = s.log.join("\n"); $("#log2").scrollTop = 1e9;
  if (s.state === "importing") return setTimeout(pollConv, 1000);
  $("#go").disabled = false;
  if (s.state === "error") { $("#err2").textContent = s.error; return; }
  if (s.state === "ready") location.href = "/";
}
$("#go").onclick = async () => {
  $("#err2").textContent = "";
  const runs = [...document.querySelectorAll("#ttab input[type=checkbox]:checked")].map((c) => RUNS[+c.dataset.i]);
  if (!runs.length) { $("#err2").textContent = "No run is ticked."; return; }
  $("#go").disabled = true; $("#log2").textContent = "";
  const r = await post("/api/xenium/batch/run", { runs, output: $("#outdir").value.trim(), merge: $("#domerge").checked, transcripts: $("#wanttx").checked, tx_table: $("#wanttab").checked, workers: +$("#workers").value || 0, he_threads: Math.max(1, Math.round(+$("#hethreads").value || 4)), ifimg: $("#wantif").checked, he_reg_folder: $("#hereg_folder").value.trim(), he_reg_index: $("#hereg_index").value.trim(),
    ...settings(), he: heChoice(), annotation: annChoice(),
    export_regions: $("#wantreg").checked, regions_name_col: $("#regcol").value.trim() || "sample_id_unblinded", regions_dir: $("#regdir").value.trim(),
    regions: Object.fromEntries(runs.filter((x) => (REG[x.name] || []).length).map((x) => [x.name, REG[x.name]])),
    comps: Object.fromEntries(runs.filter((x) => Object.keys(COMPS[x.name] || {}).length).map((x) => [x.name, COMPS[x.name]])),
    auto_off: Object.fromEntries(runs.filter((x) => AUTOOFF[x.name]).map((x) => [x.name, true])), cut_off: $("#cutoff").checked });
  if (!r.ok) { $("#err2").textContent = await detail(r); $("#go").disabled = false; return; }
  pollConv();
};
(async () => { const s = await (await fetch("/api/status")).json(); if (s.state === "importing") { $("#s3").hidden = false; $("#go").disabled = true; pollConv(); } })();

// ---- slide gallery: thumbnails of the runs with their tissue bounding boxes, a click opens one run for editing (shapes, renaming or switching off tissues)
const PAL = ["#4e9cf5", "#f5a142", "#5ec26a", "#e05d8b", "#a37ff0", "#3fc7c0", "#e6d44a", "#d9694a", "#8fb85c", "#c27dd1", "#7aa2c9", "#c9a27a"];
const fmt = (n) => n.toLocaleString();
const settings = () => ({ gap: +$("#p_gap").value, min_cells: +$("#p_minc").value, min_counts: +$("#p_mint").value, detect: $("#p_detect").value, tx_min: Math.max(0, Math.round(+$("#p_txmin").value || 0)) });
const skey = () => JSON.stringify(settings());
const shapes = (n) => (REG[n] ||= []);
const compEdit = (n, id) => ((COMPS[n] ||= {})[id] ||= {});
const galRuns = () => [...document.querySelectorAll("#ttab input[type=checkbox][data-i]:checked")].map((c) => RUNS[+c.dataset.i]).filter((r) => r.kind !== "project");
function simplify(pts, tol) {                    // Douglas-Peucker on the screen path
  if (pts.length < 3) return pts;
  const [ax, ay] = pts[0], [bx, by] = pts[pts.length - 1], dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1;
  let m = -1, k = 0; for (let i = 1; i < pts.length - 1; i++) { const d = Math.abs(dy * (pts[i][0] - ax) - dx * (pts[i][1] - ay)) / L; if (d > m) { m = d; k = i; } }
  if (m <= tol) return [pts[0], pts[pts.length - 1]];
  return simplify(pts.slice(0, k + 1), tol).slice(0, -1).concat(simplify(pts.slice(k), tol));
}
function noteFor(i) {   // the small text next to a run's button
  const n = RUNS[i].name, sh = (REG[n] || []).length, cp = Object.values(COMPS[n] || {}).filter((e) => e.name || e.keep !== undefined).length;
  return [AUTOOFF[n] ? "tissue boxes from GeoJSON" : "", sh ? `${sh} shape(s)` : "", cp ? `${cp} tissue edit(s)` : ""].filter(Boolean).join(", ");
}
function refreshNotes() { document.querySelectorAll("[data-n]").forEach((sp) => { sp.textContent = noteFor(+sp.dataset.n); }); }
// is (x, y) inside one of the polygons [[outer ring, hole, ...], ...]?
function inPolys(polys, x, y) {
  const inRing = (p) => { let c = false; for (let i = 0, j = p.length - 1; i < p.length; j = i++) if ((p[i][1] > y) !== (p[j][1] > y) && x < ((p[j][0] - p[i][0]) * (y - p[i][1])) / (p[j][1] - p[i][1]) + p[i][0]) c = !c; return c; };
  return polys.some((pl) => inRing(pl[0]) && !pl.slice(1).some(inRing));
}
// for a run whose tissue boxes were loaded from a GeoJSON: the tissue (index into tl) of every thumbnail cell, -1 = none; kept with the thumbnail data
function tissueOfCells(d, tl) {
  if (d._tl === tl) return d._lab;
  const bb = tl.map((t) => { let a = Infinity, b = Infinity, c = -Infinity, e = -Infinity; for (const pl of t.polys) for (const q of pl[0]) { a = Math.min(a, q[0]); b = Math.min(b, q[1]); c = Math.max(c, q[0]); e = Math.max(e, q[1]); } return [a, b, c, e]; });
  const n = d.comp.length, lab = new Int16Array(n).fill(-1);
  for (let i = 0; i < n; i++) { const x = d.xy[2 * i], y = d.xy[2 * i + 1]; for (let k = 0; k < tl.length; k++) { const q = bb[k]; if (x >= q[0] && x <= q[2] && y >= q[1] && y <= q[3] && inPolys(tl[k].polys, x, y)) { lab[i] = k; break; } } }
  d._tl = tl; d._lab = lab; return lab;
}
// draw one slide: cells coloured by detected tissue (or, for a run whose tissue boxes were loaded from a GeoJSON, by loaded tissue box), the tissues' boxes with names, drawn shapes
function paintSlide(cv, d, name, tall) {
  const ctx = cv.getContext("2d"), W = cv.width, H = cv.height, [x0, y0, x1, y1] = d.bounds, sc = Math.min((W - 8) / (x1 - x0 || 1), (H - 8) / (y1 - y0 || 1)), ox = 4 + (W - 8 - (x1 - x0) * sc) / 2, oy = 4 + (H - 8 - (y1 - y0) * sc) / 2;
  const imp = !!AUTOOFF[name], P = (x, y) => [ox + (x - x0) * sc, oy + (y - y0) * sc], ed = COMPS[name] || {}, byId = {};
  const keepOf = (c) => (ed[c.id] && ed[c.id].keep !== undefined ? ed[c.id].keep : imp ? false : c.keep);
  const nameOf = (c) => (ed[c.id] && ed[c.id].name) || c.name, order = [];   // colours follow the tissues' order of names, as in the editor
  // the tissue outlines: after the overlap with neighbouring tissues was removed (as exported): the editor's when this run was edited, else the loaded ones / the thumbnail's
  const ep = EDITPOLY[name] && EDITPOLY[name].key === skey() ? EDITPOLY[name].polys : null, tl = imp ? (ep || IMPORTED[name] || []) : null;
  for (const c of d.comps) if (keepOf(c) && !order.includes(nameOf(c))) order.push(nameOf(c));
  if (tl) for (const t of tl) if (!order.includes(t.name)) order.push(t.name);
  const colOf = (nm) => PAL[order.indexOf(nm) % PAL.length];
  d.comps.forEach((c) => { byId[c.id] = { col: keepOf(c) ? colOf(nameOf(c)) : "#555", c }; });
  ctx.fillStyle = "#101216"; ctx.fillRect(0, 0, W, H);
  if (d.grid) {                                   // transcript detection: the density per bin, tinted by tissue, instead of the cells
    const g = d.grid, q = g.q, bs = g.bin * sc + 0.6;
    for (let i = 0; i < q.length; i += 4) {
      const b = q[i + 3] > 0 && keepOf(byId[q[i + 3]] ? byId[q[i + 3]].c : {}) ? byId[q[i + 3]] : null, br = q[i + 2] / 255;
      ctx.globalAlpha = b ? 0.2 + 0.8 * br : 0.06 + 0.3 * br; ctx.fillStyle = b ? b.col : "#8a93a6";
      const [px, py] = P(g.ox + q[i] * g.bin, g.oy + q[i + 1] * g.bin); ctx.fillRect(px, py, bs, bs);
    }
    ctx.globalAlpha = 1;
  } else {
    const lab = tl ? tissueOfCells(d, tl) : null;
    for (let i = 0; i < d.comp.length; i++) {
      const b = byId[d.comp[i]]; ctx.fillStyle = lab && lab[i] >= 0 ? colOf(tl[lab[i]].name) : b ? b.col : "#555";
      const [px, py] = P(d.xy[2 * i], d.xy[2 * i + 1]); ctx.fillRect(px, py, 1.3, 1.3);
    }
  }
  ctx.lineWidth = 1; ctx.font = (tall ? 11 : 9) + "px system-ui"; ctx.setLineDash([4, 2]);
  const done = new Set();
  if (tl) for (const t of tl) { ctx.strokeStyle = colOf(t.name); const at = SE.strokePolys(ctx, t.polys, P); done.add(t.name); if (at) SE.label(ctx, t.name, at[0] + 2, at[1] < 12 ? at[1] + 11 : at[1] - 3, colOf(t.name)); }
  const polysOf = (c) => (ep ? (ep.find((q) => q.name === nameOf(c)) || {}).polys : d.polys && d.polys[c.id]) || null;
  for (const c of d.comps) {
    if (!c.bbox) continue; const keep = keepOf(c); if (imp && !keep) continue;   // a loaded run: its switched-off detected pieces are not drawn
    const col = keep ? byId[c.id].col : "#ff8a80"; ctx.strokeStyle = col;
    const pp = keep ? polysOf(c) : null; if (keep && done.has(nameOf(c))) continue;
    let at = pp ? SE.strokePolys(ctx, pp, P) : null; if (pp) done.add(nameOf(c));
    if (!at) { const [ax, ay] = P(c.bbox[0], c.bbox[1]), [bx, by] = P(c.bbox[2], c.bbox[3]); ctx.strokeRect(ax, ay, bx - ax, by - ay); at = [ax, ay]; }
    SE.label(ctx, nameOf(c) + (keep ? "" : " (off)"), at[0] + 2, at[1] < 12 ? at[1] + 11 : at[1] - 3, col);
  }
  ctx.strokeStyle = "#ffd54a"; ctx.fillStyle = "#ffd54a";
  for (const r of shapes(name)) {
    if (r.pts) { SE.ringsPath(ctx, r, P); ctx.stroke(); }
    else { const [ax, ay] = P(Math.min(r.x0, r.x1), Math.min(r.y0, r.y1)), [bx, by] = P(Math.max(r.x0, r.x1), Math.max(r.y0, r.y1)); ctx.strokeRect(ax, ay, bx - ax, by - ay); }
  }
  ctx.setLineDash([]);
}
const G = { page: 0, per: 20, tok: 0, timer: null };
function galRender() {
  const runs = galRuns(), pages = Math.max(1, Math.ceil(runs.length / G.per)); G.page = Math.min(G.page, pages - 1);
  const slice = runs.slice(G.page * G.per, (G.page + 1) * G.per), grid = $("#ggrid"); grid.replaceChildren();
  $("#gpage").textContent = `page ${G.page + 1} of ${pages} (${runs.length} run(s))`; $("#gprev").disabled = G.page === 0; $("#gnext").disabled = G.page >= pages - 1;
  for (const run of slice) {
    const card = document.createElement("div"), cv = document.createElement("canvas"), cap = document.createElement("div"), note = document.createElement("div");
    card.className = "gcard"; card.dataset.run = run.name; cv.width = 240; cv.height = 180; cap.className = "gname"; cap.textContent = run.name; cap.title = run.name; note.className = "sub";
    card.append(cv, cap, note); card.onclick = () => openEditor(run); grid.append(card);
  }
  galLoad();
}
async function galLoad() {
  const tok = ++G.tok; clearTimeout(G.timer);
  const runs = galRuns().slice(G.page * G.per, (G.page + 1) * G.per), key = skey(), need = runs.filter((r) => !(THUMB[r.name] && THUMB[r.name].key === key && THUMB[r.name].d));
  if (need.length) {
    const r = await post("/api/xenium/batch/thumbs", { runs: need, ...settings() });
    if (tok !== G.tok) return;
    if (!r.ok) { $("#gmsg").textContent = await detail(r); return; }
    const got = await r.json();
    for (const [n, e] of Object.entries(got)) { if (e.state === "ready") THUMB[n] = { key, d: e.d }; else if (e.state === "error") THUMB[n] = { key, error: e.error }; }
  }
  if (tok !== G.tok) return;
  let pending = 0;
  for (const run of runs) {
    const card = document.querySelector(`.gcard[data-run="${CSS.escape(run.name)}"]`); if (!card) continue;
    const t = THUMB[run.name], cv = card.querySelector("canvas"), note = card.querySelector(".sub");
    if (t && t.d && t.key === key) { paintSlide(cv, t.d, run.name, false); note.textContent = `${fmt(t.d.n)} cells, ${AUTOOFF[run.name] ? (IMPORTED[run.name] || []).length : t.d.comps.filter((c) => c.keep).length} tissue(s)` + (noteFor(RUNS.indexOf(run)) ? ", " + noteFor(RUNS.indexOf(run)) : ""); }
    else if (t && t.error) note.textContent = t.error;
    else { pending++; note.textContent = "reading…"; }
  }
  $("#gmsg").textContent = pending ? `reading ${pending} run(s)…` : ""; tisState();
  if (pending) G.timer = setTimeout(galLoad, 1500);
}
$("#galopen").onclick = () => { $("#gal").hidden = false; G.page = 0; galRender(); };
$("#gclose").onclick = () => { $("#gal").hidden = true; clearTimeout(G.timer); G.tok++; refreshNotes(); };
$("#gprev").onclick = () => { G.page--; galRender(); }; $("#gnext").onclick = () => { G.page++; galRender(); };
$("#gper").onchange = () => { G.per = +$("#gper").value; G.page = 0; galRender(); };
for (const id of ["p_gap", "p_minc", "p_mint", "p_detect", "p_txmin"]) $("#" + id).addEventListener("change", () => { document.querySelectorAll(".txonly").forEach((e) => (e.hidden = $("#p_detect").value !== "transcripts")); COMPS = {}; THUMB = {}; refreshNotes(); tisState(); if (!$("#gal").hidden) galRender(); });   // piece ids belong to the detection settings

// ---- editor of one run (the single-run page's tools and preview)
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
const D = { run: null, info: null, xy: null, prev: null, draw: null, tok: 0, view: null, sel: null, pending: null };
const dToPx = (x, y) => [4 + (x - D.view.x0) * D.view.s, 4 + (y - D.view.y0) * D.view.s];
const dToData = (px, py) => [D.view.x0 + (px - 4) / D.view.s, D.view.y0 + (py - 4) / D.view.s];
async function openEditor(run, btn) {
  const old = btn && btn.textContent; if (btn) { btn.disabled = true; btn.textContent = "reading…"; }
  const r = await post("/api/xenium/batch/draw/open", { run, detect: $("#p_detect").value });
  if (btn) { btn.disabled = false; btn.textContent = old; }
  if (!r.ok) { (($("#gal").hidden) ? $("#err") : $("#gmsg")).textContent = await detail(r); return; }
  $("#err").textContent = ""; D.run = run; D.info = await r.json(); D.draw = null; D.prev = null; D.sel = null; D.pending = null;
  D.xy = new Float32Array(await (await fetch("/api/xenium/batch/draw/xy?run=" + encodeURIComponent(run.name))).arrayBuffer());
  const [x0, y0, x1, y1] = D.info.bounds, dx = x1 - x0, dy = y1 - y0, H = 560, W = Math.max(260, Math.min(560, Math.round(H * dx / dy))), cv = $("#cv");
  cv.width = W; cv.height = H; cv.style.width = W + "px"; cv.style.height = H + "px"; D.view = { x0, y0, s: Math.min((W - 8) / dx, (H - 8) / dy) };
  $("#dtitle").textContent = `${run.name} — ${fmt(D.info.n_cells)} cells`; $("#dlg").hidden = false; setTool(null); dRefresh();
}
$("#ttab").addEventListener("click", (e) => { const b = e.target.closest("button[data-d]"); if (b) openEditor(RUNS[+b.dataset.d], b); });
$("#ddone").onclick = () => {
  $("#dlg").hidden = true; refreshNotes();
  if (!$("#gal").hidden) galLoad();                // redraw this run's card (and any other still reading)
};
async function dRefresh() {
  const tok = ++D.tok, n = D.run.name; delete EDITRAW[n];   // until this preview arrives the run's outlines are not known
  const r = await post("/api/xenium/batch/draw/preview", { run: n, spec: { ...settings(), comps: COMPS[n] || {}, rects: shapes(n), auto_off: !!AUTOOFF[n], cut_off: $("#cutoff").checked } });
  if (!r.ok) { $("#derr").textContent = await detail(r); return; }
  const d = await r.json(); if (tok !== D.tok) return; $("#derr").textContent = ""; D.prev = d; EDITPOLY[n] = { key: skey(), polys: d.tissue_polys || [] }; EDITRAW[n] = { key: skey(), polys: d.tissue_polys_raw || [] }; tisState(); dDraw(); dTable();
}
function dDraw() {
  const cv = $("#cv"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height, img = ctx.createImageData(W, H), u = new Uint32Array(img.data.buffer), v = D.view, L = D.prev.labels;
  const hex = (c) => [parseInt(c.slice(1, 3), 16), parseInt(c.slice(3, 5), 16), parseInt(c.slice(5, 7), 16)];
  const pal = D.prev.names.map((_, i) => { const [r, g, b] = hex(PAL[i % PAL.length]); return (255 << 24) | (b << 16) | (g << 8) | r; }), grey = (255 << 24) | (85 << 16) | (85 << 8) | 85;
  if (D.prev.grid) {                               // transcript density under the cells
    const ed = COMPS[D.run.name] || {}, rgb = {};
    for (const c of D.prev.comps) { const e = ed[c.id] || {}, k = D.prev.names.indexOf(e.name ?? c.name); rgb[c.id] = (e.keep ?? c.keep) && k >= 0 ? hex(PAL[k % PAL.length]) : null; }
    fillBins(u, W, H, D.prev.grid, v, (id) => rgb[id]);
  }
  for (let i = 0; i < L.length; i++) {
    const l = L[i]; if (l === -2) continue;
    const px = (4 + (D.xy[2 * i] - v.x0) * v.s) | 0, py = (4 + (D.xy[2 * i + 1] - v.y0) * v.s) | 0;
    if (px >= 0 && px < W && py >= 0 && py < H) u[py * W + px] = l < 0 ? grey : pal[l % pal.length];
  }
  ctx.putImageData(img, 0, 0); ctx.lineWidth = 1.5; ctx.font = "11px system-ui"; ctx.setLineDash([5, 3]);
  for (const r of shapes(D.run.name)) {
    ctx.strokeStyle = r.keep ? "#ffffff" : "#ff8a80"; ctx.fillStyle = ctx.strokeStyle;
    if (r.pts) { SE.ringsPath(ctx, r, dToPx); ctx.stroke(); const [x, y] = dToPx(r.pts[0][0], r.pts[0][1]); SE.label(ctx, r.name + (r.keep ? "" : " (off)"), x + 3, y < 14 ? y + 13 : y - 4, ctx.strokeStyle); }
    else { const [ax, ay] = dToPx(Math.min(r.x0, r.x1), Math.min(r.y0, r.y1)), [bx, by] = dToPx(Math.max(r.x0, r.x1), Math.max(r.y0, r.y1)); ctx.strokeRect(ax, ay, bx - ax, by - ay); SE.label(ctx, r.name + (r.keep ? "" : " (off)"), ax + 3, ay + 13, ctx.strokeStyle); }
  }
  SE.paintBoxes(ctx, dBoxes(), dToPx, dInfo, D.sel, (id, name) => polysByName(D.prev, name)); SE.paintSel(ctx, shapes(D.run.name), D.sel, dToPx);
  ctx.setLineDash([]); $("#cnote").textContent = `${fmt(D.prev.n_kept)} cells in ${D.prev.names.length} tissue(s); grey = switched off`;
}
const swatch = (name) => { const k = D.prev.names.indexOf(name); return k >= 0 ? PAL[k % PAL.length] : "#555"; };
function dTable() {
  const tb = $("#dtab tbody"), n = D.run.name, s = shapes(n); tb.innerHTML = "";
  const row = (cls, keep, name, kind, onKeep, onName, onDel, selected, onSel) => {
    const tr = document.createElement("tr"); tr.className = cls + (selected ? " selrow" : ""); tr.addEventListener("click", (e) => { if (!e.target.closest("input, .x")) onSel(); }); tr.innerHTML = `<td><input type="checkbox" ${keep ? "checked" : ""}></td><td><span class="sw" style="background:${swatch(name)}"></span><input type="text"></td><td>${kind}</td><td>${onDel ? '<span class="x" title="remove">✕</span>' : ""}</td>`;
    tr.querySelector("input[type=text]").value = name;
    tr.querySelector("input[type=checkbox]").onchange = (e) => { onKeep(e.target.checked); dRefresh(); };
    tr.querySelector("input[type=text]").onchange = (e) => { onName(e.target.value); dRefresh(); };
    if (onDel) tr.querySelector(".x").onclick = () => { D.sel = null; onDel(); dRefresh(); }; tb.append(tr);
  };
  for (const c of D.prev.comps) { const e = (COMPS[n] || {})[c.id] || {}; row("", e.keep ?? c.keep, e.name ?? c.name, `detected, ${fmt(c.n)} cells`, (k) => { compEdit(n, c.id).keep = k; }, (v) => { compEdit(n, c.id).name = v; }, undefined, D.sel && D.sel.kind === "comp" && D.sel.id === c.id, () => dSetSel({ kind: "comp", id: c.id })); }
  s.forEach((r, i) => row("rect", r.keep, r.name, r.pts ? "free-hand" : "rectangle", (k) => { r.keep = k; }, (v) => { r.name = v; }, () => { s.splice(i, 1); }, D.sel && D.sel.kind === "shape" && D.sel.i === i, () => dSetSel({ kind: "shape", i })));
  $("#dsum").textContent = Object.entries(D.prev.counts).map(([k, c]) => `${k} (${fmt(c)})`).join(", ") || "no tissue switched on";
}
// ---- selecting a box on the slide: drag a drawn shape, drag its handles, Delete removes it; a detected tissue box can be made editable
const slideB = () => (D.prev && D.prev.slide_bounds) || (D.info && D.info.bounds) || null;   // the slide: shapes are capped at it
const polysByName = (prev, name) => { const t = prev && prev.tissue_polys && prev.tissue_polys.find((q) => q.name === name); return t ? t.polys : null; };
const dBoxes = () => (D.prev && D.prev.boxes ? Object.entries(D.prev.boxes).map(([id, b]) => ({ id: +id, bbox: b })) : []);
function dInfo(id) { const c = D.prev && D.prev.comps.find((q) => q.id === id); if (!c) return null; const e = (COMPS[D.run.name] || {})[id] || {}; const name = e.name ?? c.name, k = D.prev.names.indexOf(name); return { name, keep: e.keep ?? c.keep, color: k >= 0 ? PAL[k % PAL.length] : "#ffffff" }; }
function dSetSel(p) { D.sel = p; if (D.prev) { dDraw(); dTable(); } dSelBar(); }
function dSelBar() {
  const bar = $("#dsel"), n = D.run.name, sel = D.sel, btn = (t, f) => { const b = document.createElement("button"); b.textContent = t; b.onclick = f; return b; }; bar.replaceChildren();
  if (!sel) { bar.textContent = "Click a box on the slide to select it. A drawn shape can then be dragged, reshaped by its handles, or removed with Delete; a detected tissue box can be made editable."; return; }
  if (sel.kind === "shape") { const r = shapes(n)[sel.i]; if (!r) return; bar.append(`Selected: ${r.name} (${r.pts ? "free-hand" : "rectangle"}): drag it to move, drag a handle to reshape.`, btn("Remove", () => { shapes(n).splice(sel.i, 1); D.sel = null; dRefresh(); })); return; }
  const i = dInfo(sel.id); if (!i) return;
  bar.append(`Selected detected tissue ${i.name}${i.keep ? "" : " (off)"}.`,
    btn("Make editable", () => {          // a rectangle copy of the box takes its place: the tissue is switched off, the shape keeps its name
      const b = D.prev.boxes[sel.id]; compEdit(n, sel.id).keep = false; shapes(n).push({ x0: b[0], y0: b[1], x1: b[2], y1: b[3], name: i.name, keep: true }); D.sel = { kind: "shape", i: shapes(n).length - 1 }; dRefresh();
    }),
    btn("Replace with rectangle", () => dReplace(sel.id, "rect")), btn("Replace with free-hand", () => dReplace(sel.id, "free")),
    btn(i.keep ? "Switch off" : "Switch on", () => { compEdit(n, sel.id).keep = !i.keep; dRefresh(); }));
}
// replace a detected tissue by a new shape that takes its name: the tissue is switched off, the next shape drawn is named like it (cancelling the tool puts the tissue back)
function dReplace(id, tool) {
  const i = dInfo(id), e = compEdit(D.run.name, id); D.pending = { id, name: i.name, prevKeep: e.keep };
  e.keep = false; D.sel = null; D.draw = null; setTool(tool); dRefresh();
}
function dCancelPending() { const p = D.pending; D.pending = null; if (!p) return; const e = compEdit(D.run.name, p.id); if (p.prevKeep === undefined) delete e.keep; else e.keep = p.prevKeep; dRefresh(); }
function dNewName(n) { if (D.pending) { const nm = D.pending.name; D.pending = null; return nm; } return "region_" + (shapes(n).length + 1); }
SE.attach($("#cv"), { tool: () => D.draw, view: () => D.view, bounds: () => slideB(), shapes: () => shapes(D.run.name), boxes: dBoxes, getSel: () => D.sel, setSel: dSetSel, redraw: dDraw, commit: () => dRefresh(), active: () => !$("#dlg").hidden, offComp: (id) => { compEdit(D.run.name, id).keep = false; dRefresh(); } });
function setTool(t) {
  D.draw = D.draw === t ? null : t; if (D.draw) D.sel = null; else if (D.pending) dCancelPending();
  $("#draw").textContent = D.draw === "rect" ? "Rectangle: drag on the slide (click to cancel)" : "Rectangle";
  $("#drawfree").textContent = D.draw === "free" ? "Free-hand: trace on the slide (click to cancel)" : "Free-hand";
}
$("#draw").onclick = () => setTool("rect"); $("#drawfree").onclick = () => setTool("free");
(function () {
  const cv = $("#cv"); let a = null, path = null;
  const pos = (e) => { const b = cv.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
  cv.addEventListener("mousedown", (e) => { if (D.draw === "rect") a = pos(e); else if (D.draw === "free") path = [pos(e)]; });
  cv.addEventListener("mousemove", (e) => {
    const ctx = cv.getContext("2d"), q = pos(e);
    if (a) { dDraw(); ctx.strokeStyle = "#ffd54a"; ctx.setLineDash([4, 3]); ctx.strokeRect(a[0], a[1], q[0] - a[0], q[1] - a[1]); ctx.setLineDash([]); }
    else if (path) { const l = path[path.length - 1]; if (Math.hypot(q[0] - l[0], q[1] - l[1]) < 2) return; path.push(q); dDraw(); ctx.strokeStyle = "#ffd54a"; ctx.beginPath(); path.forEach((p, k) => k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1])); ctx.stroke(); }
  });
  window.addEventListener("mouseup", (e) => {
    if (!a && !path) return;
    const s = shapes(D.run.name), done = () => { setTool(D.draw); dRefresh(); };
    if (a) {
      const b = pos(e), p0 = dToData(a[0], a[1]), p1 = dToData(b[0], b[1]); a = null;
      if (Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) > 50) { { const r = { x0: p0[0], y0: p0[1], x1: p1[0], y1: p1[1], name: dNewName(D.run.name), keep: true }; if (SE.clip(r, slideB())) { s.push(r); done(); } else dDraw(); } } else dDraw();
    } else {
      let far = 0, fd = -1; path.forEach((q, i) => { const d = Math.hypot(q[0] - path[0][0], q[1] - path[0][1]); if (d > fd) { fd = d; far = i; } });   // the path is a closed loop: split it at the point farthest from the start
      const pts = (fd < 3 ? [] : simplify(path.slice(0, far + 1), 1.2).slice(0, -1).concat(simplify(path.slice(far).concat([path[0]]), 1.2).slice(0, -1))).map((p) => dToData(p[0], p[1])); path = null;
      { const cr = { pts }; if (pts.length >= 3 && SE.clip(cr, slideB())) pts.splice(0, pts.length, ...cr.pts); else pts.length = 0; }   // capped at the slide
      let ar = 0; pts.forEach((p, i) => { const q = pts[(i + 1) % pts.length]; ar += p[0] * q[1] - q[0] * p[1]; });
      if (pts.length >= 3 && Math.abs(ar) / 2 > 2500) { s.push({ pts, name: dNewName(D.run.name), keep: true }); done(); } else dDraw();
    }
  });
})();
