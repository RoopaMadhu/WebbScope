"use strict";
const $ = (s) => document.querySelector(s);
let opening = false;                       // true once the user pressed Open on this page: only then redirect to the viewer when ready
async function status() { const r = await fetch("/api/status"); return r.json(); }
async function poll() {
  const s = await status();
  $("#log").textContent = s.log.join("\n");
  if (s.state === "importing") return setTimeout(poll, 400);
  $("#open").disabled = false;
  if (s.state === "error") { $("#err").textContent = s.error; return; }
  if (s.state === "ready" && opening) location.href = "/";
}
let askedFor = null;                         // the path whose columns the user has been asked about
let INFO = null;                             // what /api/columns said about that path
let forceAsk = false;                        // "Change columns…": show the form even when the answers are known
const opt = (sel, value, text) => { const o = document.createElement("option"); o.value = value; o.textContent = text; sel.append(o); return o; };
const nuniqText = (info, c) => { const n = info.n_unique && info.n_unique[c]; return n != null ? `${n.toLocaleString()} different value${n === 1 ? "" : "s"}` : ""; };
const libNone = () => $("#colvalues").value === "lognorm" ? "(none)" : "(none: the sum of each cell's counts)";
// One select per role: the likely columns first, then every other column; optional roles start with a "none" entry.
function roleRow(info, cols) {
  const row = document.createElement("div"); row.className = "row crow";
  const lab = document.createElement("label"); lab.textContent = info.label; lab.htmlFor = "col_" + info.role; if (info.required) lab.className = "req";
  const sel = document.createElement("select"); sel.id = "col_" + info.role; sel.dataset.role = info.role;
  const n = document.createElement("span"); n.className = "sub";
  if (!info.required) opt(sel, "", info.role === "lib_size" ? libNone() : info.none || "(none)");
  for (const c of info.candidates) opt(sel, c, c);
  const rest = cols.filter((c) => !info.candidates.includes(c));
  if (rest.length) { const g = document.createElement("optgroup"); g.label = info.candidates.length ? "other columns" : "columns"; for (const c of rest) opt(g, c, c); sel.append(g); }
  sel.value = info.default != null ? info.default : "";
  const upd = () => { n.textContent = nuniqText(info, sel.value); };
  sel.onchange = upd; upd();
  row.append(lab, sel, n);
  return row;
}
function genesFor() {                       // gene-name choices of the gene table of the log-normalised matrix (else the counts'); raw/X has its own
  const key = $("#colln").value || $("#colcounts").value;
  const sel = $("#colgenes"), keep = sel.value;
  if (!key) { sel.innerHTML = ""; $("#colgenesn").textContent = ""; return; }
  const m = INFO.matrices.find((x) => x.key === key), g = INFO.genes.options[m.var];
  sel.innerHTML = ""; opt(sel, "", "the index of the gene table");
  for (const c of g.columns) opt(sel, c, c);
  sel.value = sel.dataset.touched && [...sel.options].some((o) => o.value === keep) ? keep : INFO.genes.default[m.var] || "";
  const ex = () => { const e = g.examples[sel.value] || []; $("#colgenesn").textContent = e.length ? "e.g. " + e.join(", ") : ""; };
  sel.onchange = () => { sel.dataset.touched = "1"; ex(); }; ex();
}
const hasCounts = () => INFO.matrices ? !!$("#colcounts").value : INFO.files ? !!INFO.files.counts : true;    // tables: told by the file names
const hasLn = () => INFO.matrices ? !!$("#colln").value : INFO.files ? !!INFO.files.lognorm : false;
function valuesNote() {
  const d = INFO.values.detected;
  $("#colvaln").textContent = d ? (d.integer ? "the stored values are whole numbers" : "the stored values are not whole numbers") + (d.negative ? ", some are negative (scaled data cannot be used)" : "") : "";
  const c = hasCounts(), l = hasLn();
  $("#scalerow").hidden = !INFO.full || !c || l;          // the scale is used only when the values are computed here
  if (INFO.matrices) $("#colmnote").textContent = c && l ? "Both: the log-normalised values are shown as they are; the counts give the raw-counts layer and the QC numbers."
    : c ? "Raw counts only: WebbScope log-normalises them." : l ? "Log-normalised values only: shown as they are; there is no raw-counts layer." : "Choose at least one of the two.";
  const lib = $("#col_lib_size"); if (lib && lib.options[0] && lib.options[0].value === "") lib.options[0].textContent = c ? "(none: the sum of each cell's counts)" : "(none)";
}
function showForm(info) {
  INFO = info;
  $("#colsintro").innerHTML = info.full
    ? `${info.h5ad ? `<b>${info.h5ad.replace(/[<>&]/g, "")}</b>${info.file_only ? " (only this file is read: cell outlines and other files in its folder are not; give the folder to use them)" : ""}: ` : ""}say which part of the ${info.h5ad ? "file" : "cell table"} is which. The choices were guessed from the names${info.remembered ? " or are the ones made the last time this project was opened" : ""}. Starred ones are needed. Then press Open.`
    : `Which columns of the cell table identify the <b>samples</b> and the <b>donors</b>? They were found by name (a column with "sample" or "donor" in it); change them if they are wrong, then press Open.`;
  const box = $("#colrows"); box.innerHTML = "";
  for (const r of info.roles) box.append(roleRow(r, info.columns));
  $("#h5rows").hidden = !info.matrices;
  if (info.matrices) {
    for (const [id, dflt] of [["#colcounts", info.counts_matrix], ["#colln", info.lognorm_matrix]]) {
      const sel = $(id); sel.innerHTML = ""; opt(sel, "", "(none)");
      for (const m of info.matrices) opt(sel, m.key, `${m.key} (${m.shape[1].toLocaleString()} genes, ${m.kind === "dense" ? "dense" : "sparse"}; ${m.integer ? "whole numbers" : "decimals"}${m.negative ? ", negative values" : ""})`);
      sel.value = dflt || "";
      sel.onchange = () => { genesFor(); valuesNote(); };
    }
    genesFor();
  }
  const f = info.files;
  $("#colfiles").hidden = !f;
  if (f) $("#colfiles").textContent = "Files: " + [f.counts && `raw counts = ${f.counts}`, f.lognorm && `log-normalised values = ${f.lognorm}`, f.cells && `cell metadata = ${f.cells}`].filter(Boolean).join(" · ");
  $("#valrow").hidden = !info.values.show;
  $("#colvalues").value = info.values.default;
  $("#colvalues").onchange = valuesNote;
  $("#colscale").value = info.norm_scale || "";
  const pj = info.project_json || {};
  $("#saverow").hidden = !info.full || info.file_only;     // a file alone: nothing is written next to it (the choices stay in the cache)
  $("#savepj").disabled = !pj.writable;
  $("#savepj").checked = pj.writable && (!pj.exists || pj.ours);
  $("#savelab").textContent = !pj.writable ? "project.json cannot be saved: the folder is not writable (the choices are kept in WebbScope's cache)"
    : pj.exists && !pj.ours ? "Save these choices as project.json in the folder (replaces the project.json that is there)"
    : "Save these choices as project.json in the folder, so they travel with the data (otherwise they are kept in WebbScope's cache)";
  valuesNote();
  $("#colsbox").hidden = false;
}
function formAnswers() {
  if ($("#colsbox").hidden || !INFO) return {};
  const val = (role) => { const e = $("#col_" + role); return e ? e.value : undefined; };
  const out = { sample_col: val("sample"), donor_col: val("donor") };
  if (INFO.values.show) out.values = $("#colvalues").value;
  if (!INFO.full) return out;
  out.roles = {};
  for (const r of INFO.roles) if (r.role !== "sample" && r.role !== "donor") out.roles[r.role] = val(r.role);
  if (INFO.matrices) { out.counts_matrix = $("#colcounts").value; out.lognorm_matrix = $("#colln").value; out.genes = $("#colgenes").value; }
  if (!$("#scalerow").hidden) out.norm_scale = $("#colscale").value.trim();
  out.save_project_json = $("#savepj").checked && !$("#savepj").disabled;
  return out;
}
async function openPath(p) {
  $("#err").textContent = "";
  if (askedFor !== p) {                      // first press for this path: ask which columns are which (projects this app wrote itself are opened directly)
    let info = null;
    try { const r = await fetch("/api/columns?path=" + encodeURIComponent(p)); if (r.ok) info = await r.json(); else { $("#err").textContent = (await r.json()).detail || ""; return; } } catch (e) { info = null; }
    askedFor = p;
    if (info && info.ask && (!info.known || forceAsk)) { forceAsk = false; showForm(info); return; }   // the user checks the columns and presses Open again
    if (info && info.known) $("#log").textContent = "";
    $("#colsbox").hidden = true; INFO = null;            // answers known (project.json or the last opening): open directly
  }
  const ans = formAnswers();
  if (INFO && INFO.matrices && !$("#colcounts").value && !$("#colln").value) { $("#err").textContent = "Choose the raw counts, the log-normalised values, or both."; return; }
  const unset = INFO && INFO.full ? INFO.roles.filter((r) => r.required && !(r.role === "sample" ? ans.sample_col : (ans.roles || {})[r.role])) : [];
  if (unset.length) { $("#err").textContent = "Choose a column for: " + unset.map((r) => r.label).join(", "); return; }
  if (ans.norm_scale && !(+ans.norm_scale > 0)) { $("#err").textContent = "The normalisation scale must be a number above 0, or empty."; return; }
  $("#open").disabled = true; $("#log").textContent = ""; opening = true;
  const r = await fetch("/api/open", { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify({ ...ans, path: p, roi_object: $("#roi").value.trim(), if_image: $("#ifimg").value.trim(), imc_file: $("#imcfile").value.trim(), he_reg_folder: $("#hereg").value.trim(), he_reg_index: $("#heregidx").value.trim() }) });
  if (!r.ok) { $("#err").textContent = await r.text(); $("#open").disabled = false; return; }
  poll();
}
$("#open").onclick = () => { const p = $("#path").value.trim(); if (p) openPath(p); };
$("#chcols").onclick = (e) => { e.preventDefault(); const p = $("#path").value.trim(); if (!p) { $("#err").textContent = "Type the project path first."; return; } askedFor = null; forceAsk = true; openPath(p); };
$("#path").addEventListener("input", () => { askedFor = null; INFO = null; $("#colsbox").hidden = true; });   // another path: ask again
$("#path").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });   // addEventListener: a handler that returns false would cancel typing
$("#browse").hidden = false; $("#browseh5").hidden = false; $("#roibrowse").hidden = false; $("#heregbrowse").hidden = false; $("#heregidxbrowse").hidden = false; $("#ifbrowse").hidden = false; $("#imcbrowse").hidden = false;     // native dialog in the app window, else the local server's macOS dialog; the path can always be typed
$("#browse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api ? await api.pick_folder() : (await (await fetch("/api/pick?kind=dir")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a folder dialog here; type the path instead."; return; }
  if (p) { $("#path").value = p; $("#open").click(); }
};
$("#browseh5").onclick = async () => {        // an .h5ad file itself (the folder dialog cannot select files)
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api && api.pick_file ? await api.pick_file(["h5ad"]) : (await (await fetch("/api/pick?kind=file")).json()).path; }   // macOS' own dialog: no type filter (.h5ad has no registered type, so a filter could grey it out)
  catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return; }
  if (p && !/\.h5ad$/i.test(p)) { $("#err").textContent = `Not an .h5ad file: ${p}`; return; }
  if (p) { $("#path").value = p; askedFor = null; INFO = null; $("#colsbox").hidden = true; $("#open").click(); }
};
$("#roibrowse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api && api.pick_file ? await api.pick_file(["rds", "qs2"]) : (await (await fetch("/api/pick?kind=file&types=rds,qs2")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return; }
  if (p) $("#roi").value = p;
};
$("#heregbrowse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api ? await api.pick_folder() : (await (await fetch("/api/pick?kind=dir")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a folder dialog here; type the path instead."; return; }
  if (p) $("#hereg").value = p;
};
$("#heregidxbrowse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api && api.pick_file ? await api.pick_file(["json"]) : (await (await fetch("/api/pick?kind=file&types=json")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return; }
  if (p) $("#heregidx").value = p;
};
$("#ifbrowse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api && api.pick_file ? await api.pick_file(["tif", "tiff"]) : (await (await fetch("/api/pick?kind=file&types=tif,tiff")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return; }
  if (p) $("#ifimg").value = p;
};
$("#imcbrowse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api && api.pick_file ? await api.pick_file(["mcd"]) : (await (await fetch("/api/pick?kind=file&types=mcd")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a file dialog here; type the path instead."; return; }
  if (p) $("#imcfile").value = p;
};
$("#ifimg").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });
$("#roi").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });
(async () => {
  const s = await status();
  if (s.open && s.project) { $("#current").hidden = false; $("#cur").textContent = s.project; }
  if (s.state === "importing") { opening = true; poll(); }
})();
