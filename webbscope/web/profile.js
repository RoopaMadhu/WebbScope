"use strict";
const $ = (s) => document.querySelector(s);
let LAST = null;
const f1 = (x) => (x == null ? "–" : (+x).toFixed(1)), f2 = (x) => (x == null ? "–" : (+x).toFixed(2));
const cls = (s) => (s >= 10 ? "vslow" : s >= 2 ? "slow" : "");
const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" })[c]);
const tab = (id, head, rows) => { $(id).innerHTML = `<tr>${head.map((h, i) => `<th class="${i ? "n" : ""}">${h}</th>`).join("")}</tr>` + rows.join(""); };
async function load() {
  const d = LAST = await (await fetch("/api/profile")).json();
  $("#head").textContent = d.project ? `${d.project} — opened in ${f1(d.open_seconds)} s; process peak memory ${f2(d.peak_gb)} GB` + (d.rss_gb != null ? `, now ${f2(d.rss_gb)} GB` : "") : "No project has been opened in this session yet.";
  tab("#steps", ["step", "seconds", "peak GB"], d.steps.map((s) => `<tr><td>${esc(s.name)}</td><td class="n ${cls(s.seconds)}">${f2(s.seconds)}</td><td class="n">${f2(s.peak)}</td></tr>`));
  tab("#events", ["message", "t (s)", "+dt (s)", "RSS GB", "peak GB"], d.events.map((e) => `<tr><td>${esc(e.msg)}</td><td class="n">${f1(e.t)}</td><td class="n ${cls(e.dt)}">${f2(e.dt)}</td><td class="n">${f2(e.rss)}</td><td class="n">${f2(e.peak)}</td></tr>`));
  tab("#reqs", ["endpoint", "calls", "total ms", "mean ms", "max ms", "last ms", "MB sent"], d.requests.map((r) => `<tr><td>${esc(r.endpoint)}</td><td class="n">${r.n}</td><td class="n">${r.total_ms}</td><td class="n">${r.mean_ms}</td><td class="n ${cls(r.max_ms / 1000)}">${r.max_ms}</td><td class="n">${r.last_ms}</td><td class="n">${r.mb}</td></tr>`));
}
$("#refresh").onclick = load;
$("#reset").onclick = async () => { await fetch("/api/profile/reset", { method: "POST" }); load(); };
$("#copy").onclick = () => {
  const d = LAST; if (!d) return;
  const t = [`project: ${d.project}`, `open: ${f1(d.open_seconds)} s, peak ${f2(d.peak_gb)} GB`, "", "steps:", ...d.steps.map((s) => `  ${f2(s.seconds)} s  ${s.name}`), "", "log:",
    ...d.events.map((e) => `  ${f1(e.t).padStart(7)} s  +${f2(e.dt).padStart(6)}  rss ${f2(e.rss)}  peak ${f2(e.peak)}  ${e.msg}`), "", "requests:", ...d.requests.map((r) => `  ${String(r.n).padStart(5)} x  total ${r.total_ms} ms  mean ${r.mean_ms}  max ${r.max_ms}  ${r.mb} MB  ${r.endpoint}`)].join("\n");
  navigator.clipboard.writeText(t).then(() => { $("#copy").textContent = "Copied"; setTimeout(() => ($("#copy").textContent = "Copy as text"), 1500); }, () => { $("#head").textContent = "Could not copy; select the tables instead."; });
};
let T = null;
$("#auto").onchange = () => { clearInterval(T); if ($("#auto").checked) T = setInterval(load, 2000); };
load();
