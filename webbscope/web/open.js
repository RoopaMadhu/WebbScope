"use strict";
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
let forceAsk = false;                        // "Change columns…": show the form even when the answers are known
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
  const ans = formAnswers(), bad = formProblem(ans);
  if (bad) { $("#err").textContent = bad; return; }
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
  if (s.state === "importing") { opening = true; poll(); return; }
  const want = new URLSearchParams(location.search).get("path");   // from the start page's Recent list: open it (asks only if its columns are not known)
  if (want) { $("#path").value = want; history.replaceState(null, "", "/open"); openPath(want); }
})();
