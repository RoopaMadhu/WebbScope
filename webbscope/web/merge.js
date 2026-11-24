"use strict";
// The merge page: find the datasets in a folder, tick them, answer the column form once for those whose columns are not known, merge and open.
let SCAN = null;                             // /api/merge/scan of the folder
let FORMFOR = null;                          // path of the dataset the column form was filled from
FORM = { verb: "Merge and open", save: false, intro: "These answers are used for every ticked dataset without saved choices, so their columns should have the same names. " };
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const ticked = () => SCAN ? SCAN.datasets.filter((d, i) => $(`#use${i}`).checked) : [];

async function find() {
  $("#err").textContent = ""; $("#err2").textContent = "";
  const p = $("#path").value.trim(); if (!p) return;
  const r = await fetch("/api/merge/scan?path=" + encodeURIComponent(p));
  if (!r.ok) { $("#err").textContent = (await r.json()).detail || "cannot read the folder"; return; }
  SCAN = await r.json();
  const tb = $("#ttab tbody"); tb.innerHTML = "";
  SCAN.datasets.forEach((d, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td><input type="checkbox" id="use${i}" checked></td><td>${esc(d.name)}</td><td>${esc(d.kind)}</td>`
      + `<td>${d.known ? "saved" : "<b>from the form below</b>"}</td><td class="sub">${esc(d.path)}</td>`;
    tb.append(tr);
    $(`#use${i}`).onchange = refresh;
  });
  $("#skipped").innerHTML = SCAN.skipped.length ? "Not datasets to merge here: " + SCAN.skipped.map((s) => `<b>${esc(s.name)}</b> (${esc(s.reason)})`).join("; ") : "";
  $("#earlier").innerHTML = SCAN.merged.length ? "Earlier merges in this folder: " + SCAN.merged.map((m) => `<a href="/open?path=${encodeURIComponent(m)}">${esc(m.split("/").pop())}</a>`).join(" · ") + " (click to open)" : "";
  $("#s2").hidden = false;
  if (!SCAN.datasets.length) { $("#err").textContent = "No datasets found in this folder."; $("#s3").hidden = $("#s4").hidden = true; return; }
  if (!$("#outdir").value) $("#outdir").value = SCAN.parent.replace(/\/$/, "") + "/_merged_datasets";
  FORMFOR = null;
  await refresh();
}

async function refresh() {                   // the column form is needed when a ticked dataset has no saved choices; it is filled from the first such dataset
  const need = ticked().filter((d) => !d.known);
  $("#s4").hidden = false;
  $("#s4h").textContent = need.length ? "4 · Merge" : "3 · Merge";
  if (!need.length) { $("#s3").hidden = true; INFO = null; FORMFOR = null; return; }
  $("#s3").hidden = false;
  if (FORMFOR === need[0].path) return;
  const r = await fetch("/api/columns?path=" + encodeURIComponent(need[0].path));
  if (!r.ok) { $("#err2").textContent = `${need[0].name}: ` + ((await r.json()).detail || "cannot read its columns"); return; }
  const info = await r.json();
  info.full = true;                           // every column is asked (also for a folder whose project.json names only some)
  if (!info.roles.some((x) => x.role === "cell_id")) { $("#err2").textContent = `${need[0].name}: open it alone first to choose its columns`; return; }
  FORMFOR = need[0].path;
  showForm(info);
  $("#colsintro").innerHTML = `Filled from <b>${esc(need[0].name)}</b>; ${need.length > 1 ? `used for all ${need.length} datasets without saved choices (${need.map((d) => esc(d.name)).join(", ")}), so their columns should have the same names. ` : ""}`
    + `The choices were guessed from the names. Starred ones are needed. A dataset that differs can be opened alone first (<a href="/open">One dataset</a>): its choices are then remembered and used here.`;
}

async function poll() {
  const s = await (await fetch("/api/status")).json();
  $("#log2").textContent = s.log.join("\n"); $("#log2").scrollTop = 1e9;
  if (s.state === "importing") return setTimeout(poll, 500);
  $("#go").disabled = false;
  if (s.state === "error") { $("#err2").textContent = s.error; return; }
  if (s.state === "ready") location.href = "/";
}

async function go() {
  $("#err2").textContent = "";
  const ds = ticked();
  if (ds.length < 2) { $("#err2").textContent = "Tick at least two datasets."; return; }
  let ans = {};
  if (!$("#s3").hidden) {
    if (!INFO) { $("#err2").textContent = "The column form could not be filled; see the message above."; return; }
    ans = formAnswers(); delete ans.save_project_json;
    const bad = formProblem(ans); if (bad) { $("#err2").textContent = bad; return; }
  }
  const out = $("#outdir").value.trim(); if (!out) { $("#err2").textContent = "Give the output folder."; return; }
  $("#go").disabled = true; $("#log2").textContent = "";
  const r = await fetch("/api/merge", { method: "POST", headers: { "content-type": "application/json" },
    body: JSON.stringify({ ...ans, datasets: ds.map((d) => ({ name: d.name, path: d.path })), out }) });
  if (!r.ok) { $("#err2").textContent = (await r.json()).detail || "could not start"; $("#go").disabled = false; return; }
  poll();
}

$("#find").onclick = find;
$("#path").addEventListener("keydown", (e) => { if (e.key === "Enter") find(); });
$("#go").onclick = go;
$("#browse").hidden = false;
$("#browse").onclick = async () => {
  let p = null;
  try { const api = window.pywebview && window.pywebview.api; p = api ? await api.pick_folder() : (await (await fetch("/api/pick?kind=dir")).json()).path; }
  catch (e) { $("#err").textContent = "Cannot open a folder dialog here; type the path instead."; return; }
  if (p) { $("#path").value = p; find(); }
};
(async () => { const s = await (await fetch("/api/status")).json(); if (s.state === "importing") { $("#s4").hidden = false; $("#go").disabled = true; poll(); } })();
