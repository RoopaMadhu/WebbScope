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
const fillCols = (sel, info, none) => {
  const rest = info.all.filter((c) => !info.candidates.includes(c));
  sel.innerHTML = (none ? `<option value="">(none: use the sample id)</option>` : "") + info.candidates.map((c) => `<option value="${c}">${c}</option>`).join("") + (rest.length ? `<optgroup label="other columns">${rest.map((c) => `<option value="${c}">${c}</option>`).join("")}</optgroup>` : "");
  sel.value = info.default || (none ? "" : info.candidates[0] || info.all[0]);
};
const nuniq = (info, sel, lab) => { const n = info.n_unique && info.n_unique[sel.value]; lab.textContent = n != null ? `${n.toLocaleString()} different value${n === 1 ? "" : "s"}` : ""; };
async function openPath(p) {
  $("#err").textContent = "";
  if (askedFor !== p) {                      // first press for this path: ask which columns are the sample id and the donor id (projects this app wrote itself are opened directly)
    let info = null;
    try { const r = await fetch("/api/columns?path=" + encodeURIComponent(p)); if (r.ok) info = await r.json(); else { $("#err").textContent = (await r.json()).detail || ""; return; } } catch (e) { info = null; }
    askedFor = p;
    if (info && info.ask) {
      const mk = (role) => ({ ...info[role], all: info.columns });
      const s = mk("sample"), d = mk("donor");
      fillCols($("#colsample"), s, false); fillCols($("#coldonor"), d, true);
      nuniq(s, $("#colsample"), $("#colsamplen")); nuniq(d, $("#coldonor"), $("#coldonorn"));
      $("#colsample").onchange = () => nuniq(s, $("#colsample"), $("#colsamplen")); $("#coldonor").onchange = () => nuniq(d, $("#coldonor"), $("#coldonorn"));
      $("#colsbox").hidden = false; return;   // the user checks the two columns and presses Open again
    }
    $("#colsbox").hidden = true;
  }
  $("#open").disabled = true; $("#log").textContent = ""; opening = true;
  const cols = $("#colsbox").hidden ? {} : { sample_col: $("#colsample").value, donor_col: $("#coldonor").value };
  const r = await fetch("/api/open", { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify({ ...cols, path: p, roi_object: $("#roi").value.trim(), if_image: $("#ifimg").value.trim(), he_reg_folder: $("#hereg").value.trim(), he_reg_index: $("#heregidx").value.trim() }) });
  if (!r.ok) { $("#err").textContent = await r.text(); $("#open").disabled = false; return; }
  poll();
}
$("#open").onclick = () => { const p = $("#path").value.trim(); if (p) openPath(p); };
$("#path").addEventListener("input", () => { askedFor = null; $("#colsbox").hidden = true; });   // another path: ask again
$("#path").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });   // addEventListener: a handler that returns false would cancel typing
$("#browse").hidden = false; $("#roibrowse").hidden = false; $("#heregbrowse").hidden = false; $("#heregidxbrowse").hidden = false; $("#ifbrowse").hidden = false;     // native dialog in the app window, else the local server's macOS dialog; the path can always be typed
$("#browse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api ? await api.pick_folder() : (await (await fetch("/api/pick?kind=dir")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a folder dialog here; type the path instead."; return; }
  if (p) { $("#path").value = p; $("#open").click(); }
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
$("#ifimg").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });
$("#roi").addEventListener("keydown", (e) => { if (e.key === "Enter") $("#open").click(); });
(async () => {
  const s = await status();
  if (s.open && s.project) { $("#current").hidden = false; $("#cur").textContent = s.project; }
  if (s.state === "importing") { opening = true; poll(); }
})();
