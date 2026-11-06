"use strict";
const $ = (s) => document.querySelector(s);
const post = (u, b) => fetch(u, { method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify(b) });
let RUNS = [], PARENT_MERGED = null, HE = null;   // HE = {czi: [{name, path}], annotations: [...], match: {run: {czi, annotations, why}}} once an H&E folder was matched
async function detail(r) { try { return JSON.parse(await r.text()).detail; } catch (e) { return "failed"; } }
$("#find").onclick = async () => {
  $("#err").textContent = "";
  const r = await post("/api/xenium/batch/scan", { path: $("#path").value.trim() });
  if (!r.ok) { $("#err").textContent = await detail(r); return; }
  const d = await r.json();
  RUNS = d.runs; PARENT_MERGED = d.merged || null;
  if (!RUNS.length) { $("#err").textContent = "No Xenium runs found directly inside this folder (looked for Xenium output sub-folders with cell_feature_matrix/, cells and cell_boundaries, .tar.gz or .zip bundles, and project folders already written by this app)." + (d.skipped.length ? " Skipped sub-folders: " + d.skipped.join(", ") : ""); return; }
  $("#runsum").textContent = `${RUNS.length} run(s)` + (d.skipped.length ? `; ${d.skipped.length} sub-folder(s) skipped because they are not Xenium outputs: ${d.skipped.join(", ")}` : "");
  HE = null; $("#hemsg").textContent = ""; document.querySelectorAll(".hecol").forEach((e) => (e.hidden = true));
  $("#ttab tbody").innerHTML = RUNS.map((x, i) => `<tr><td><input type="checkbox" data-i="${i}" checked></td><td></td><td>${x.kind}</td><td class="hecol" hidden></td><td class="hecol" hidden></td></tr>`).join("");
  document.querySelectorAll("#ttab tbody tr").forEach((tr, i) => { tr.children[1].textContent = RUNS[i].name; });
  $("#outdir").value = d.default_output;
  $("#s2").hidden = false; $("#s3").hidden = false; checkMerged();
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
$("#ttab").addEventListener("change", (e) => { if (e.target.type === "checkbox") checkMerged(); });
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
    tr.children[3].replaceChildren(a, why); tr.children[4].replaceChildren(b);
  });
  msg.textContent = `${HE.czi.length} scan(s) and ${HE.annotations.length} annotation file(s) found; ${n} of ${RUNS.length} run(s) matched by name. Runs left at (none) get no H&E.`;
};
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
  const r = await post("/api/xenium/batch/run", { runs, output: $("#outdir").value.trim(), merge: $("#domerge").checked, transcripts: $("#wanttx").checked, tx_table: $("#wanttab").checked, workers: +$("#workers").value || 0, ifimg: $("#wantif").checked, he_reg_folder: $("#hereg_folder").value.trim(), he_reg_index: $("#hereg_index").value.trim(),
    gap: +$("#p_gap").value, min_cells: +$("#p_minc").value, min_counts: +$("#p_mint").value, he: heChoice() });
  if (!r.ok) { $("#err2").textContent = await detail(r); $("#go").disabled = false; return; }
  pollConv();
};
(async () => { const s = await (await fetch("/api/status")).json(); if (s.state === "importing") { $("#s3").hidden = false; $("#go").disabled = true; pollConv(); } })();
