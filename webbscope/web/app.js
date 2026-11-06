/* Xenium spatial browser — single-state vanilla JS client (deck.gl for WebGL rendering). */
const $ = (id) => document.getElementById(id);
const MISSING = 65535;
const ROLE_ORDER = ["cell_type", "cell_state", "region_lvl1", "region_lvl2", "region_lvl3", "glom_quality"];
// Categorical palette (dark-surface steps of the validated reference palette, slots 1-8) + 3 extension colours for fields with up to 11
// levels. Checked with validate_palette.js against the viewer surface #0f1218: all 11 pass the lightness band, chroma floor, adjacent-pair CVD
// (worst 8.4 >= 8), adjacent normal-vision (worst 18.1 >= 15) and 3:1 contrast. Only the first 3 are all-pairs safe (see legends/tooltips for identity).
const PAL = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767", "#9033bd", "#9e4b58", "#5457e5"];
const STATUS = { good: "#0ca30c", critical: "#d03b3b" };  // fixed status palette: passed / failed qc (always shown with ✓/✕ + label)
const GRAY = [60, 66, 78];      // gene not detected
const DIM = [105, 110, 122];    // cell filtered out (drawn translucent so context stays visible)

// ---------------------------------------------------------------- state
const S = {
  meta: null, sample: null, n: 0,
  xy: null, ids: null, cats: {},          // per-sample arrays (cats: key -> Uint16Array)
  nums: {}, genes: {},                    // numeric metadata / gene vectors loaded for this sample
  filters: {},                            // key -> Set(level index); empty/undefined = no restriction
  exprFilters: [], exprV: [],             // gene / signature filters {name, genes, op, thr}; their per-cell scores for the current tissue
  hiddenCats: {},                         // colour key -> Set(level index) hidden via legend
  colorBy: "cell_type", geneA: "", geneB: "", sigGenes: [], sigVec: null, sigKey: "", sigGenes2: [], sigVec2: null, sigKey2: "", sigName2: "", smooth: 0, zscope: "none", within: "tissue", layer: "logcpx", clip: 99.5, numScale: { log: false, cmin: null, cmax: null, fmin: null, fmax: null },
  psize: 5, opacity: 0.9, selected: new Set(), cellSearch: "",
  vis: null, visIdx: null, drawIdx: null, // passes-filter mask; passing cells; all cells (filtered-out first, so highlights draw on top)
  ovLegend: false, ctxs: {}, boxMode: false, showPath: false, shape: "auto", polyZoom: 1, polyMax: 60000, mode: "pt", poly: null, colors: null, viewState: null, compField: "cell_type", paths: null,
};
let deckgl = null;

// ---------------------------------------------------------------- helpers
const loading = (on) => { $("loading").style.display = on ? "block" : "none"; };
async function jget(url) { const r = await fetch(url); if (!r.ok) throw new Error(`${url}: ${(await r.text()).slice(0, 200)}`); return r.json(); }
async function bget(url, T) { const r = await fetch(url); if (!r.ok) throw new Error(`${url}: ${(await r.text()).slice(0, 200)}`); return new T(await r.arrayBuffer()); }
const debounce = (f, ms = 250) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => f(...a), ms); }; };
const hex2rgb = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
const fmt = (n) => n.toLocaleString();
const catMeta = (k) => S.meta.cats[k];
const display = (k, i) => (catMeta(k).display ? catMeta(k).display[i] : catMeta(k).levels[i]);
// Stable colour per level (never depends on filters). Fields: slot = global frequency rank. Cell states: slot = position inside the parent cell
// type (natural order), so the states of one cell type always get distinct colours; states of different cell types may share a colour, which is why the
// legend / filter lists group states under their cell type. Fields with more levels than slots reuse slots (rare, e.g. region ids).
// Multi-membership fields (pathologist regions) store a bit mask per cell: a cell can lie in several nested regions, and a filter on such a
// field keeps the cells inside ANY of the ticked regions.
const isMulti = (k) => !!S.meta?.cats[k]?.multi;
const selMask = (set) => { let m = 0; for (const i of set) m |= 1 << i; return m; };
// Cell types use the Tableau 20 palette. Its ten saturated colours go to the most frequent cell types and the pale partners to the rarest, so the
// big cell types never get a dark/light pair of the same hue; Tableau's grey is skipped (grey already means "filtered out" / "not detected" here).
// The order of the first 11 is the one whose neighbouring pairs are easiest to tell apart (checked with the palette validator: worst adjacent
// colour-vision difference 12.0, normal-vision 16.9), e.g. orange and green are never neighbours.
const TAB20_CT = ["#1f77b4", "#ff7f0e", "#d62728", "#9467bd", "#2ca02c", "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#ffbb78", "#aec7e8",
  "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5"];
// Pathologist regions get their own neon palette (fixed per region name), kept away from the Tableau hues used for cell types, so an outline or
// region swatch is never mistaken for a cell type (the pen colours in the files repeat: Cortex and Perivascular are both red, Medulla and Subcapsular orange).
const REGION_COL = { "Cortex": "#00ffd5", "Medulla": "#ff2bd6", "IFTA": "#ffffff", "All gloms": "#b6ff00", "CMJ": "#ffe600", "Subcapsular Fibrosis": "#ff4d6d",
  "Perivascular Fibrosis": "#7dffb3", "Gloms w/endocapillary hypercellularity": "#c77dff", "Gloms w/periglomerular fibrosis": "#4dd2ff", "Sclerotic Gloms": "#ff9f1c", "Non-sclerotic Gloms": "#b6ff00",
  "Artifacts": "#a0a8b8", "Intraglom Crescents": "#f15bb5", "Poorly Preserved Gloms": "#d4ff8f", "(outside all outlines)": "#6b7482" };
const regionColor = (name, pen) => REGION_COL[String(name || "").trim()] || pen || "#ffffff";
const catColor = (k, i) => {
  if (k === "cell_type") return TAB20_CT[i % TAB20_CT.length];
  if (k === "path_region") return regionColor(catMeta(k).levels[i], catMeta(k).colors[i]);
  if (k === "glom_quality") { const n = catMeta(k).levels[i]; return n === "passed qc" ? STATUS.good : n === "failed qc" ? STATUS.critical : PAL[i % PAL.length]; }
  if (k === "cell_state" && S.stateInfo) return PAL[S.stateInfo.slot[i] % PAL.length];
  return PAL[i % PAL.length];
};
// Cell states grouped by parent cell type (cell-type filter order), natural sort inside a group (PT0, PT1, ... M2, M12).
function buildStateInfo() {
  S.stateInfo = null; const cs = S.meta.cats.cell_state, ct = S.meta.cats.cell_type; if (!cs?.parent || !ct) return;
  const nat = (a, b) => cs.levels[a].localeCompare(cs.levels[b], undefined, { numeric: true });
  const groups = [], slot = new Array(cs.levels.length).fill(0);
  for (let t = -1; t < ct.levels.length; t++) {
    const m = cs.levels.map((_, i) => i).filter((i) => cs.parent[i] === t).sort(nat); if (!m.length) continue;
    m.forEach((i, r) => { slot[i] = r; }); groups.push({ type: t, name: t < 0 ? "Unassigned" : ct.levels[t], levels: m });
  }
  S.stateInfo = { groups, slot, rank: Object.fromEntries(groups.flatMap((g, gi) => g.levels.map((i, r) => [i, [gi, r]]))) };
}
const ROLE_LABEL = { cell_type: "Cell type", cell_state: "Cell state", region_lvl1: "Region level 1",
  region_lvl2: "Region level 2", region_lvl3: "Region level 3", glom_quality: "Glom quality", path_region: "Pathologist region (outlines)" };
const labelOf = (k) => ROLE_LABEL[k] || (k.startsWith("meta:") ? k.slice(5) : k);

// viridis-like LUT (9 anchor stops, linear interpolation)
const VIR = [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],[31,158,137],[53,183,121],[109,205,89],[253,231,37]];
const LUT = Array.from({ length: 256 }, (_, i) => {
  const t = (i / 255) * (VIR.length - 1), a = Math.floor(t), b = Math.min(a + 1, VIR.length - 1), f = t - a;
  return VIR[a].map((v, j) => Math.round(v + (VIR[b][j] - v) * f));
});

const DIV = [[59, 76, 192], [141, 176, 254], [221, 221, 221], [244, 154, 123], [180, 4, 38]];
const DIVL = Array.from({ length: 256 }, (_, i) => {
  const t = (i / 255) * (DIV.length - 1), a = Math.floor(t), b = Math.min(a + 1, DIV.length - 1), f = t - a;
  return DIV[a].map((v, j) => Math.round(v + (DIV[b][j] - v) * f));
});

async function ensureCat(k) {
  if (!S.cats[k]) S.cats[k] = await bget(`/api/sample/${encodeURIComponent(S.sample)}/cat/${encodeURIComponent(k)}`, Uint16Array);
  return S.cats[k];
}
async function ensureNum(k) {
  if (!S.nums[k]) S.nums[k] = await bget(`/api/sample/${encodeURIComponent(S.sample)}/num/${encodeURIComponent(k)}`, Float32Array);
  return S.nums[k];
}
// Expression vectors come either straight from the matrix (/gene) or, when smoothing and/or z-scoring is on, from the
// server (/expr): mean over the genes -> Gaussian smoothing (sigma µm) -> z-score within the chosen scope.
// "Compute over": smoothing and the z-score reference can be limited to the glomeruli (each glom smoothed on its own) or to the region picked in the
// Regions window; cells outside that subset come back as NaN and are drawn grey.
const RG = { src: "path_region", level: -1, open: false };   // region chosen in the Regions window (also used by "Compute over: that region")
const withinSpec = () => (S.within === "glom" ? "glom" : S.within === "region" && RG.level >= 0 ? `cat:${RG.src}:${RG.level}` : "tissue");
const rgName = () => (RG.level >= 0 ? display(RG.src, RG.level) : "");
const transformed = () => S.smooth > 0 || S.zscope !== "none" || withinSpec() !== "tissue";
const gkey = (genes) => (transformed() ? `E|${genes.join(",")}|${S.layer}|${S.smooth}|${S.zscope}|${withinSpec()}` : genes[0] + "|" + S.layer);
const exprUrl = (genes) => `/api/sample/${encodeURIComponent(S.sample)}/expr?genes=${encodeURIComponent(genes.join(","))}&layer=${S.layer}&sigma=${S.smooth}&scope=${S.zscope}&within=${encodeURIComponent(withinSpec())}`;
const SCOPE_LABEL = { sample: "within each donor", subtype: "across samples of the same disease subtype", type: "across samples of the same disease type", all: "across all samples", nocontrol: "against the non-control samples (No LN + LN; Control excluded from the mean/SD)" };
const withinNote = () => { const w = withinSpec(); return w === "glom" ? "glomerulus cells only (each glom smoothed on its own; z-score mean/SD from glomerulus cells)" : w !== "tissue" ? `cells of the region “${rgName()}” only (smoothing and z-score mean/SD from that region)` : ""; };
const transformNote = () => [S.smooth > 0 ? `Gaussian-smoothed, σ = ${S.smooth} µm` : "", S.zscope !== "none" ? `z-score ${SCOPE_LABEL[S.zscope]}` : "", withinNote()].filter(Boolean).join("; ");
async function ensureGene(g) {
  const key = gkey([g]);
  if (!S.genes[key]) S.genes[key] = await bget(transformed() ? exprUrl([g]) : `/api/sample/${encodeURIComponent(S.sample)}/gene/${encodeURIComponent(g)}?layer=${S.layer}`, Float32Array);
  return S.genes[key];
}
const geneVec = (g) => S.genes[gkey([g])];

// ---------------------------------------------------------------- init
async function init() {
  loading(true);
  S.meta = await jget("/api/meta");
  if (S.meta.default_grey) { S.grey = true; $("greycells").checked = true; }   // merged batch projects: cell types are per tissue, so start uncoloured
  const m = S.meta; buildStateInfo(); rgnInit();
  for (const x of m.samples) if (!x.donor) x.donor = String(x.id).replace(/_fragment[-_]?\d+$/, "");   // older servers do not send the donor: derive it from the sample id
  $("dsinfo").textContent = `${fmt(m.n_cells)} Xenium cells (QC'ed) · ${fmt(m.n_genes)} genes · ${donTxt(nDonors(m.samples))} · matrix ${m.matrix.orientation}`;
  $("sample").innerHTML = m.samples.map((s) => `<option value="${s.id}">${dlabel(s.id)} — ${subtypeOf(s)} (${fmt(s.n)} cells)</option>`).join("");
  const colorOpts = [];
  for (const k of ROLE_ORDER) if (m.cats[k]) colorOpts.push([k, ROLE_LABEL[k]]);
  colorOpts.push(["__gene", "Gene expression (one gene)"], ["__dual", "Two genes (dual colour)"], ["__sig", "Gene signature (mean of genes)"], ["__sigdual", "Two gene signatures (dual colour)"]);
  for (const k of Object.keys(m.cats)) if (k.startsWith("meta:")) colorOpts.push([k, "meta: " + k.slice(5)]);
  for (const k of Object.keys(m.numeric)) colorOpts.push([k, "numeric: " + k.slice(5)]);
  $("colorby").innerHTML = colorOpts.map(([k, l]) => `<option value="${k}">${l}</option>`).join("");
  $("compfield").innerHTML = ROLE_ORDER.filter((k) => m.cats[k]).map((k) => `<option value="${k}">${ROLE_LABEL[k]}</option>`).join("");
  buildFilterUI();
  bindUI();
  await loadSample(m.samples[0].id);
}

async function loadSample(id) {
  loading(true);
  S.sample = id; S.cats = {}; S.nums = {}; S.genes = {}; S.exprV = []; S.selected = new Set(); S.pathLayer = null; S.pathInfo = null;
  const info = S.meta.samples.find((s) => s.id === id);
  S.n = info.n;
  const [xy, ids] = await Promise.all([bget(`/api/sample/${encodeURIComponent(id)}/xy`, Float32Array), jget(`/api/sample/${encodeURIComponent(id)}/ids`)]);
  S.xy = xy; S.ids = ids; S.poly = null;
  if (S.meta.has_polygons) {
    const q = `/api/sample/${encodeURIComponent(id)}/poly/`;
    const [off, pxy] = await Promise.all([bget(q + "offsets", Int32Array), bget(q + "xy", Float32Array)]);
    for (let k = 1; k < pxy.length; k += 2) pxy[k] = -pxy[k];  // same y flip as centroids
    S.poly = { off, xy: pxy };
  }
  await Promise.all([...ROLE_ORDER, "path_region"].filter((k) => S.meta.cats[k]).map(ensureCat));
  if (S.showPath) await loadPaths();
  syncPathUI(); S.cmpOverride = ""; renderCmp();
  S.viewState = null;
  await refresh(true);
  loading(false);
}

// ---------------------------------------------------------------- filters
function buildFilterUI() {
  const keys = ROLE_ORDER.filter((k) => S.meta.cats[k]), gi = keys.indexOf("glom_quality");
  keys.splice(gi < 0 ? keys.length : gi, 0, "path_region");   // after the region levels
  const pathBox = `<label class="chk"><input type="checkbox" id="showpath"> draw the outlines on the viewer</label><div id="pathstat" class="muted small"></div>
    <div class="muted small">${S.meta.cats.path_region ? "Cells inside any ticked region (regions can be nested, e.g. gloms lie inside the cortex)." : "Region membership has not been built (run scripts/build_paths.py)."}</div>`;
  $("filters").innerHTML = keys.map((k) => `<details class="fgroup" data-k="${k}" ${k === "cell_type" ? "open" : ""}><summary>${ROLE_LABEL[k]} <span class="muted" id="fcount-${k}"></span></summary>
    ${k === "path_region" ? pathBox : ""}<input class="fs" placeholder="search…" data-k="${k}"><div class="opts" id="fopts-${k}"></div></details>`).join("") +
    exprGroupHTML();
  keys.forEach(renderFilterOpts);
  document.querySelectorAll(".fs").forEach((el) => el.addEventListener("input", (e) => renderFilterOpts(e.target.dataset.k)));
  bindExprUI(); renderExprList();
}

// Gene / signature filters: keep the cells whose score (the gene's value, or the mean over several genes, with the Layer / Smoothing / Z-score /
// Compute-over settings currently chosen) is > thr ("positive") or <= thr ("negative"). 0 on the raw layers = detected. They combine (AND) with
// the other filters and are applied in the viewer, the windows and the server-side analyses (neighbours, glom heterogeneity, signature summary).
const exprPass = (e, x) => (e.op === "gt" ? x > e.thr : x <= e.thr);   // NaN (outside the "compute over" subset) fails both
const exprSpec = (e) => ({ name: e.name, genes: e.genes, op: e.op, thr: e.thr, layer: e.layer, sigma: e.sigma, scope: e.scope, within: e.within });   // each filter carries its own transform
const exprSpecs = () => S.exprFilters.map(exprSpec);
const exprDesc = (e) => { const t = [e.layer === "counts" ? "raw counts" : "", e.sigma > 0 ? `σ ${e.sigma} µm` : "", e.scope !== "none" ? "z-score" : "", e.within !== "tissue" ? "compute-over subset" : ""].filter(Boolean).join(", ");
  return `${e.name} ${e.op === "gt" ? "positive" : "negative"} (score ${e.op === "gt" ? ">" : "≤"} ${e.thr}${t ? "; " + t : ""})`; };
function exprGroupHTML() {
  return `<details class="fgroup" data-k="__expr"><summary>Gene / signature positive <span class="muted" id="fcount-__expr"></span></summary><div class="xf">
    <input id="xfgenes" list="genelist" placeholder="gene, or several genes (mean = signature)…" autocomplete="off">
    <select id="xfpreset" title="use a built-in signature instead"><option value="">— or a built-in signature —</option></select>
    <div class="xfrow"><select id="xfop"><option value="gt">positive: score &gt;</option><option value="le">negative: score ≤</option></select><input id="xfthr" type="number" step="any" value="0" title="threshold on the score; 0 on the raw layers means detected"><button id="xfadd" class="sm">Add</button></div>
    <div class="xfrow"><label class="small">Layer</label><select id="xflayer"><option value="logcpx">log-normalised</option><option value="counts">raw counts</option></select></div>
    <div class="xfrow"><label class="small">Smooth σ</label><input type="range" id="xfsigma" min="0" max="100" step="5" value="0"><span id="xfsigmalab" class="small">off</span></div>
    <div class="xfrow"><label class="small">Z-score</label><select id="xfscope"><option value="none">none (raw values)</option><option value="sample">within each donor</option><option value="subtype">across samples of the same disease subtype</option><option value="type">across samples of the same disease type</option><option value="all">across all samples</option><option value="nocontrol">across all non-control samples</option></select></div>
    <div class="xfrow"><label class="small">Compute over</label><select id="xfwithin"><option value="tissue">the whole tissue</option><option value="glom">glomeruli only</option><option value="region">the region chosen in the Regions window</option></select></div>
    <div id="xfmsg" class="muted small"></div><div id="xflist"></div>
    <div class="muted small">Score = the gene's value (mean over the genes for several) with the Layer / Smoothing / Z-score / Compute-over settings chosen here for this filter (independent of the viewer's colouring settings), so 0 on the raw layers means “detected”. Several filters and the other filters combine with AND. The number is the cells of this tissue passing that filter alone.</div></div></details>`;
}
function renderExprList() {
  const box = $("xflist"); if (!box) return;
  const sp = exprSpecs();
  box.innerHTML = sp.map((e, q) => { const v = S.exprV[q], n = v && v.length === S.n ? v.reduce((c, x) => c + (exprPass(e, x) ? 1 : 0), 0) : null;
    return `<div class="xfi"><span class="l" title="${exprDesc(e)} — genes: ${e.genes.join(", ")}">${exprDesc(e)}</span><span class="n">${n == null ? "" : fmt(n)}</span><button class="sm" data-q="${q}" title="remove this filter">✕</button></div>`; }).join("");
  box.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => { S.exprFilters.splice(+b.dataset.q, 1); refresh(); }));
  $("fcount-__expr").textContent = sp.length ? `(${sp.length} active)` : "";
}
async function addExprFilter() {
  const msg = $("xfmsg"), op = $("xfop").value, thr = parseFloat($("xfthr").value), preset = $("xfpreset").value; msg.textContent = "";
  if (!Number.isFinite(thr)) { msg.textContent = "The threshold must be a number."; return; }
  if (S.exprFilters.length >= 6) { msg.textContent = "At most 6 gene / signature filters."; return; }
  let genes, name;
  try {
    if (preset) { await loadSigList(); const sg = SIG.byId[preset]; genes = [...sg.genes]; name = sg.name; }
    else {
      const names = $("xfgenes").value.split(/[\s,;]+/).filter(Boolean);
      if (!names.length) { msg.textContent = "Enter a gene (or several), or choose a signature."; return; }
      const r = await jget("/api/genes/resolve?names=" + encodeURIComponent(names.join(",")));
      if (!r.found.length) { msg.textContent = `Not in the panel: ${r.missing.join(", ")}`; return; }
      genes = r.found; name = genes.length === 1 ? genes[0] : genes.length <= 3 ? genes.join("+") : `Signature (${genes.length} genes)`;
      if (r.missing.length) msg.textContent = `Not in the panel (left out): ${r.missing.join(", ")}`;
    }
  } catch (e) { msg.textContent = "Error: " + e.message; return; }
  const wv = $("xfwithin").value;
  if (wv === "region" && RG.level < 0) { msg.textContent = "Choose a region in the Regions window first."; return; }
  const tf = { layer: $("xflayer").value, sigma: +$("xfsigma").value, scope: $("xfscope").value, within: wv === "glom" ? "glom" : wv === "region" ? `cat:${RG.src}:${RG.level}` : "tissue" };
  const dup = S.exprFilters.some((e) => e.op === op && e.thr === thr && e.genes.join(",") === genes.join(",") && e.layer === tf.layer && e.sigma === tf.sigma && e.scope === tf.scope && e.within === tf.within);
  if (!dup) S.exprFilters.push({ name, genes, op, thr, ...tf });
  $("xfgenes").value = ""; $("xfpreset").value = ""; refresh();
}
function bindExprUI() {
  $("xfadd").addEventListener("click", addExprFilter);
  $("xfsigma").addEventListener("input", (e) => { $("xfsigmalab").textContent = +e.target.value ? e.target.value + " µm" : "off"; });
  $("xfgenes").addEventListener("keydown", (e) => { if (e.key === "Enter") addExprFilter(); });
  $("xfpreset").addEventListener("change", () => { if ($("xfpreset").value) $("xfgenes").value = ""; });
  $("xfgenes").addEventListener("input", debounce(async (e) => { const q = e.target.value; if (!q || /[\s,;]/.test(q)) return; const l = await jget("/api/genes/search?q=" + encodeURIComponent(q)); $("genelist").innerHTML = l.map((g) => `<option value="${g}">`).join(""); }));
}
const fltText = (spec) => [...Object.entries(spec.filters).map(([k, ids]) => `${labelOf(k)} = ${ids.map((i) => display(k, i)).join(", ")}`), ...(spec.expr || []).map(exprDesc)].join(" · ");
const UMX = {};
async function exprMaskAll() {   // server: 1 = the cell passes every gene / signature filter, for all cells in cache order (the UMAP shows all tissues)
  const body = JSON.stringify({ expr: exprSpecs() });
  if (UMX.key !== body) { const r = await fetch("/api/expr_mask", { method: "POST", headers: { "Content-Type": "application/json" }, body }); if (!r.ok) throw new Error((await r.text()).slice(0, 200)); UMX.arr = new Uint8Array(await r.arrayBuffer()); UMX.key = body; }
  return UMX.arr;
}

function renderFilterOpts(k) {
  const el = $("fopts-" + k); if (!el || !S.meta.cats[k]) return;
  const q = (document.querySelector(`.fs[data-k="${k}"]`)?.value || "").toLowerCase();
  const sel = S.filters[k] || new Set();
  const avail = Array.from({ length: catMeta(k).levels.length }, (_, i) => i);
  const counts = S.counts?.[k];
  const row = (i) => `<label class="opt" title="${display(k, i)}"><input type="checkbox" data-k="${k}" data-i="${i}" ${sel.has(i) ? "checked" : ""}>
      <span class="sw" style="background:${catColor(k, i)}"></span><span class="l">${display(k, i)}</span><span class="n">${counts ? fmt(counts[i] || 0) : ""}</span></label>`;
  const match = (i) => display(k, i).toLowerCase().includes(q);
  el.innerHTML = k === "cell_state" && S.stateInfo   // all states of the first cell type, then the next cell type, ...
    ? S.stateInfo.groups.map((g) => { const m = g.levels.filter(match); return m.length ? `<div class="ghead">${g.name}</div>` + m.map(row).join("") : ""; }).join("")
    : avail.filter(match).slice(0, 400).map(row).join("");
  $("fcount-" + k).textContent = sel.size ? `(${sel.size} selected)` : "";
  el.querySelectorAll("input").forEach((c) => c.addEventListener("change", () => {
    const set = S.filters[k] || (S.filters[k] = new Set());
    c.checked ? set.add(+c.dataset.i) : set.delete(+c.dataset.i);
    refresh();
  }));
}

// ---------------------------------------------------------------- visibility, colours
function computeVisible() {
  const n = S.n, vis = new Uint8Array(n).fill(1);
  for (const [k, set] of Object.entries(S.filters)) {
    if (!set.size || !S.cats[k]) continue;
    const c = S.cats[k], mu = isMulti(k), sm = mu ? selMask(set) : 0;
    for (let i = 0; i < n; i++) if (vis[i] && !(mu ? (c[i] & sm) !== 0 : set.has(c[i]))) vis[i] = 0;
  }
  S.exprFilters.forEach((e, q) => { const v = S.exprV[q]; if (!v || v.length !== n) return; for (let i = 0; i < n; i++) if (vis[i] && !exprPass(e, v[i])) vis[i] = 0; });
  const hid = S.hiddenCats[S.colorBy];
  if (hid?.size && S.cats[S.colorBy]) { const c = S.cats[S.colorBy]; for (let i = 0; i < n; i++) if (vis[i] && hid.has(c[i])) vis[i] = 0; }
  if (S.cellSearch) { const q = S.cellSearch.toLowerCase(); for (let i = 0; i < n; i++) if (vis[i] && !S.ids[i].toLowerCase().includes(q)) vis[i] = 0; }
  S.vis = vis;
  const idx = []; for (let i = 0; i < n; i++) if (vis[i]) idx.push(i);
  S.visIdx = Int32Array.from(idx);
  const draw = new Int32Array(n); let a = 0;
  for (let i = 0; i < n; i++) if (!vis[i]) draw[a++] = i;
  draw.set(S.visIdx, a); S.drawIdx = draw;
}

function quantile(arr, q) { if (!arr.length) return 0; const a = Float32Array.from(arr).sort(); return a[Math.min(a.length - 1, Math.floor(q * (a.length - 1)))]; }

// Colour for every cell in S.drawIdx. Cells that fail the filters are drawn translucent grey; stats/clipping use passing cells only.
function geneList() {  // [label, vector] pairs for the active gene colouring mode
  if (S.colorBy === "__sig") return S.sigVec ? [[sigLabel(), S.sigVec]] : [];
  if (S.colorBy === "__sigdual") return [S.sigVec ? [sigLabel(), S.sigVec] : null, S.sigVec2 ? [sigLabel2(), S.sigVec2] : null].filter(Boolean);
  if (S.colorBy !== "__gene" && S.colorBy !== "__dual") return [];
  return (S.colorBy === "__dual" ? [S.geneA, S.geneB] : [S.geneA]).filter((g) => g && geneVec(g)).map((g) => [g, geneVec(g)]);
}
// A "bundle" holds one sample's per-cell arrays: { n, pass (Uint8Array), codes | num | A [, B] }.
// makeScale() turns the current Colour-by settings plus the passing cells of one or more bundles into
// { info (legend description), color(bundle, i) -> [r,g,b] }. The main view uses one bundle; the Gallery and
// Glomeruli windows pool all the samples they show so colours are comparable across tiles.
// Values labelled along a continuous colour bar: 5 labels, t = position 0..1. A log-scaled numeric colouring is spaced in log space, so its labels are too.
const fmtTick = (v) => { const a = Math.abs(v); if (a === 0) return "0"; if (a >= 1000) return Math.round(v).toLocaleString("en-US"); if (a >= 100) return v.toFixed(0); if (a >= 10) return v.toFixed(1).replace(/\.0$/, ""); if (a >= 1) return v.toFixed(2).replace(/0$/, "").replace(/\.$/, ""); return v.toPrecision(2); };
function scaleTicks(info, n = 5) {
  let f;
  if (info.type === "num") { const L = (v) => (info.log ? Math.log1p(Math.max(v, 0)) : v), a = L(info.lo), b = L(info.hi); f = (t) => (info.log ? Math.expm1(a + (b - a) * t) : a + (b - a) * t); }
  else if (info.type === "gene") f = info.div ? (t) => -info.hiA + 2 * info.hiA * t : (t) => info.hiA * t;
  else return null;
  return Array.from({ length: n }, (_, i) => ({ t: i / (n - 1), v: f(i / (n - 1)) }));
}
const tickRowHTML = (ticks) => `<div class="tickrow">${ticks.map((k) => `<span style="left:${(k.t * 100).toFixed(1)}%">${fmtTick(k.v)}</span>`).join("")}</div>`;
const numLegend = (ci) => `${ci.log ? "log scale, log(1 + value) · " : ""}ends: ${ci.autoLo ? "p1" : "set by you"} and ${ci.autoHi ? "p99" : "set by you"}${ci.filt ? " · cells outside the keep-range are grey" : ""}; values beyond the ends take the end colours`;
const Z_LIM = 3;   // z-score colour range is fixed at ±3 SD so the choice of z-score scope stays visible (an auto-ranged scale would cancel it)
function makeScale(bundles) {
  const mode = S.colorBy;
  const pool = (get, keep) => { const o = []; for (const b of bundles) { const v = get(b); if (v) for (let i = 0; i < b.n; i++) if (b.pass[i] && v[i] === v[i] && (!keep || keep(v[i]))) o.push(v[i]); } return o; };   // NaN = outside the "compute over" subset
  if (mode === "__sigdual") {   // two signatures at once: A blue, B orange, both high -> near white (colour-blind-safe pair, additive mixing)
    const hasA = bundles.some((b) => b.A), hasB = bundles.some((b) => b.B), z = S.zscope !== "none";
    if (!hasA && !hasB) return { info: { type: "none" }, color: () => GRAY };
    const lim = (get) => { const v = pool(get), hi = z ? Z_LIM : (v.length ? quantile(v, S.clip / 100) : 1), lo = z ? 0 : (v.length ? quantile(v, 0.01) : 0); return { lo, hi: Math.max(hi, lo + 1e-6) }; };
    const lA = hasA ? lim((b) => b.A) : null, lB = hasB ? lim((b) => b.B) : null, cA = hex2rgb(SIG_A), cB = hex2rgb(SIG_B), t = (v, l) => Math.max(0, Math.min(1, (v - l.lo) / (l.hi - l.lo)));
    return { info: { type: "sigdual", z, labelA: sigLabel(), labelB: hasB ? sigLabel2() : null },
      color: (b, i) => { if ((b.A && b.A[i] !== b.A[i]) || (b.B && b.B[i] !== b.B[i])) return GRAY; const tA = lA && b.A ? t(b.A[i], lA) : 0, tB = lB && b.B ? t(b.B[i], lB) : 0, m = Math.max(tA, tB);
        return [0, 1, 2].map((q) => Math.min(255, Math.round(GRAY[q] * (1 - m) + tA * cA[q] + tB * cB[q]))); } };
  }
  if (mode === "__gene" || mode === "__dual" || mode === "__sig") {
    const has = (k) => bundles.some((b) => b[k]), hasA = has("A"), hasB = mode === "__dual" && has("B");
    if (!hasA) return { info: { type: "none" }, color: () => GRAY };
    const label = mode === "__sig" ? sigLabel() : S.geneA, sig = mode === "__sig";
    if (transformed() && !hasB) {  // smoothed and/or z-scored: continuous scale (sequential, or diverging around 0 for z-scores)
      const Z = S.zscope !== "none", arr = pool((b) => b.A);
      const hiq = quantile(arr, S.clip / 100), loq = quantile(arr, 1 - S.clip / 100);
      const mx = Z ? Z_LIM : Math.max(Math.abs(hiq), Math.abs(loq), 1e-6), lo = Z ? -mx : 0, hi = Z ? mx : Math.max(hiq, 1e-6), L = Z ? DIVL : LUT;
      return { info: { type: "gene", cont: true, div: Z, lo, hiA: hi, label, sig }, color: (b, i) => (b.A[i] !== b.A[i] ? GRAY : L[Math.round(Math.max(0, Math.min(1, (b.A[i] - lo) / (hi - lo))) * 255)]) };
    }
    const top = (get) => Math.max(quantile(pool(get, (v) => v > 0), S.clip / 100), 1e-6), hA = top((b) => b.A), hB = hasB ? top((b) => b.B) : 1;
    const info = { type: hasB ? "dual" : "gene", hiA: hA, hiB: hB, label, sig };
    if (!hasB) return { info, color: (b, i) => { const v = b.A[i]; return !(v > 0) ? GRAY : LUT[Math.round(Math.min(1, v / hA) * 255)]; } };
    return { info, color: (b, i) => {  // dual: A→red, B→green, both→yellow; zero (or negative z) = neutral gray
      const va = b.A[i], vb = b.B ? b.B[i] : 0; if (!(va > 0) && !(vb > 0)) return GRAY;
      let r = Math.round(255 * Math.min(1, (va > 0 ? va : 0) / hA)), g = Math.round(255 * Math.min(1, (vb > 0 ? vb : 0) / hB));
      if (r < 40 && g < 40) { r = Math.max(r, 60); g = Math.max(g, 60); }
      return [r, g, 0];
    } };
  }
  if (S.meta.cats[mode]) {
    const cl = catMeta(mode).levels.map((_, k) => hex2rgb(catColor(mode, k)));
    return { info: { type: "cat" }, color: (b, i) => { const c = b.codes[i]; return c === MISSING ? GRAY : cl[c]; } };
  }
  if (bundles.some((b) => b.num)) {   // numeric column: optional log scale, colour-bar limits (values beyond them take the end colours) and a keep-range (others grey)
    const ns = S.numScale, tf = ns.log ? (v) => Math.log1p(Math.max(v, 0)) : (v) => v, inF = (v) => (ns.fmin == null || v >= ns.fmin) && (ns.fmax == null || v <= ns.fmax);
    const arr = pool((b) => b.num, (v) => !Number.isNaN(v) && inF(v)), lo = ns.cmin ?? (arr.length ? quantile(arr, 0.01) : 0), hi = Math.max(ns.cmax ?? (arr.length ? quantile(arr, 0.99) : 1), lo + 1e-9);
    const tlo = tf(lo), sp = Math.max(tf(hi) - tlo, 1e-9);
    return { info: { type: "num", lo, hi, log: ns.log, autoLo: ns.cmin == null, autoHi: ns.cmax == null, filt: ns.fmin != null || ns.fmax != null }, color: (b, i) => { const x = b.num[i]; return Number.isNaN(x) || !inF(x) ? GRAY : LUT[Math.round(Math.max(0, Math.min(1, (tf(x) - tlo) / sp)) * 255)]; } };
  }
  return { info: { type: "none" }, color: () => GRAY };
}
function mainBundle() {
  const mode = S.colorBy, b = { n: S.n, pass: S.vis };
  if (mode === "__gene" || mode === "__dual") { if (S.geneA) b.A = geneVec(S.geneA); if (mode === "__dual" && S.geneB) b.B = geneVec(S.geneB); }
  else if (mode === "__sig") b.A = S.sigVec;
  else if (mode === "__sigdual") { b.A = S.sigVec; b.B = S.sigVec2; }
  else if (S.cats[mode]) b.codes = S.cats[mode];
  else if (S.nums[mode]) b.num = S.nums[mode];
  return b;
}
const GREY = [150, 155, 162], GREY_LINE = [128, 134, 143, 235];   // "Grey cells": plain fill and outline, so transcripts and other overlays stand out
function computeColors() {
  const di = S.drawIdx, m = di.length, out = new Uint8Array(m * 4), a = Math.round(S.opacity * 255), ad = Math.round(a * 0.45);
  const b = mainBundle(), sc = makeScale([b]); S.colorInfo = sc.info;
  for (let j = 0; j < m; j++) {
    const i = di[j], c = S.vis[i] ? (S.grey ? GREY : sc.color(b, i)) : DIM, al = S.vis[i] ? a : ad;
    out[j * 4] = c[0]; out[j * 4 + 1] = c[1]; out[j * 4 + 2] = c[2]; out[j * 4 + 3] = al;
  }
  return out;
}

// ---------------------------------------------------------------- deck.gl
function bounds() {
  const inf = S.meta.samples.find((s) => s.id === S.sample); return { x: inf.x, y: inf.y };
}
// deck.gl only applies initialViewState when it differs from the previous value, so a repeated reset (e.g. clicking Fit
// after panning/zooming) would be ignored. A changing counter makes every call register; layers are rebuilt for the new view.
function setView(vs) {
  S.viewState = vs; S.viewN = (S.viewN || 0) + 1;
  deckgl.setProps({ initialViewState: { ...vs, _n: S.viewN }, ...(S.colors ? { layers: layers(S.colors) } : {}) });
  updateScaleBar(); if (HEW.main.side && HEW.main.link && !scanMode()) heSideView(vs);
}
function fitView() {
  const b = bounds(), el = $("deckwrap"), w = el.clientWidth || 800, h = el.clientHeight || 600;
  const dx = b.x[1] - b.x[0], dy = b.y[1] - b.y[0];
  S.viewState = { target: [(b.x[0] + b.x[1]) / 2, -(b.y[0] + b.y[1]) / 2, 0], zoom: Math.log2(Math.min(w / dx, h / dy) * 0.92) };
  if (deckgl) setView(S.viewState);
  updateScaleBar();
  if (HEW.main.side) { if (scanMode()) scanFit(); else heSideView(S.viewState); }   // the image beside the cells is reset too, linked or not
}

// Indices j (into S.drawIdx) of cells whose centroid lies in the current viewport (+10% margin); null if too many for polygons.
function cellsInView() {
  const vs = S.viewState, el = $("deckwrap"), sc = Math.pow(2, vs.zoom);
  const hw = (el.clientWidth || 800) / 2 / sc * 1.1, hh = (el.clientHeight || 600) / 2 / sc * 1.1;
  const x0 = vs.target[0] - hw, x1 = vs.target[0] + hw, y0 = vs.target[1] - hh, y1 = vs.target[1] + hh;
  const out = [], di = S.drawIdx;
  for (let j = 0; j < di.length; j++) {
    const x = S.xy[di[j] * 2], y = -S.xy[di[j] * 2 + 1];
    if (x >= x0 && x <= x1 && y >= y0 && y <= y1) { out.push(j); if (out.length > S.polyMax) return null; }
  }
  return out;
}
function usePolys() {
  if (!S.poly || S.shape === "centroid" || !S.viewState) return null;
  if (S.shape === "auto" && S.viewState.zoom < S.polyZoom) return null;
  return cellsInView();  // null when too many cells in view to draw as polygons
}
const withPaths = (L) => { if (S.showPath && S.pathLayer) L.push(S.pathLayer()); return L; };  // pathologist outlines on top of the cells
// H&E under the cells in the main viewer: one BitmapLayer per visible pyramid tile (y is flipped: deck y = -um y).
function heLayers(vs = S.viewState, el = $("deckwrap"), redo = heMainRedo, list0 = null, w = "main") {   // BitmapLayers of the image tiles visible in view vs (list0: pyramids to draw, e.g. the slide scenes; default = what window w has ticked in its Images menu: H&E first, then IF)
  const list = list0 || heSel(S.sample, w); if (!list || !list.length || !vs) return [];
  const sc = 2 ** vs.zoom, wd = (el.clientWidth || 800) / sc, h = (el.clientHeight || 600) / sc, cx = vs.target[0], cy = -vs.target[1], out = [];
  const layer = (im, z, tx, ty, img) => { const sp = im.tsize * im.px * 2 ** z, X0 = im.x0 + tx * sp, Y0 = im.y0 + ty * sp, f = im.px * 2 ** z;
    const o = { id: `he-${im.key}-${z}-${tx}-${ty}`, image: img, bounds: [X0, -(Y0 + img.naturalHeight * f), X0 + img.naturalWidth * f, -Y0], pickable: false };
    if (im.key.startsWith("if:")) o.parameters = { blend: true, blendFunc: [1, 769] };   // IF adds its light to what is under it ("screen": ONE, ONE_MINUS_SRC_COLOR)
    return o; };
  for (const im of list) {
    const back = new Map(), front = [];
    const T = im.tsize, u = im.px, z = Math.max(0, Math.min(im.nz - 1, Math.floor(Math.log2(1 / (u * sc))))), span = T * u * 2 ** z;
    const t0 = Math.max(0, Math.floor((cx - wd / 2 - im.x0) / span)), t1 = Math.min(Math.ceil(im.w * u / span) - 1, Math.floor((cx + wd / 2 - im.x0) / span));
    const r0 = Math.max(0, Math.floor((cy - h / 2 - im.y0) / span)), r1 = Math.min(Math.ceil(im.h * u / span) - 1, Math.floor((cy + h / 2 - im.y0) / span));
    for (let ty = r0; ty <= r1; ty++) for (let tx = t0; tx <= t1; tx++) {
      const img = heTile(im.key, z, tx, ty, redo);
      if (img) { front.push(layer(im, z, tx, ty, img)); continue; }
      for (let zz = z + 1; zz < im.nz; zz++) {   // exact tile still loading: show the whole coarser parent underneath
        const d = zz - z, p = heTile(im.key, zz, tx >> d, ty >> d, null);
        if (p) { const c = layer(im, zz, tx >> d, ty >> d, p); back.set(c.id, c); break; }
      }
    }
    out.push(...back.values(), ...front);
  }
  return out.map((c) => new deck.BitmapLayer(c));
}
const heMainRedo = () => { clearTimeout(S.heT); S.heT = setTimeout(() => { if (HEW.main.on && S.colors) deckgl.setProps({ layers: layers(S.colors) }); }, 60); };
function layers(colors) { let L = layersCore(colors); if (!HEW.main.cells) L = L.filter((l) => !/^(cells|sel)-/.test(l.id)); const base = HEW.main.on ? [...heLayers(), ...L] : L, tl = txLayer(), pl = ptLayer(); ptUi(); return [...base, ...(tl ? [tl] : []), ...(pl ? [pl] : [])]; }
function layersCore(colors) {
  const di = S.drawIdx, m = di.length, L = [];
  const pl = usePolys();
  if (pl) {
    S.mode = "poly";
    const { off, xy } = S.poly;
    L.push(new deck.PolygonLayer({
      id: "cells-poly", data: pl, getPolygon: (j) => xy.subarray(off[di[j]] * 2, off[di[j] + 1] * 2),
      getFillColor: (j) => colors.subarray(j * 4, j * 4 + 4), filled: !HEW.main.outl,
      getLineColor: S.grey ? GREY_LINE : HEW.main.outl ? (j) => { const o = j * 4; return [colors[o], colors[o + 1], colors[o + 2], Math.max(200, colors[o + 3])]; } : [20, 20, 20, 140],
      stroked: true, lineWidthUnits: "pixels", getLineWidth: S.grey ? 0.9 : HEW.main.outl ? 1.4 : 0.6, positionFormat: "XY", pickable: true,
      updateTriggers: { getPolygon: [S.sample, S.drawIdx, pl.length, pl[0], pl[pl.length - 1]], getFillColor: [colors], getLineColor: [colors, HEW.main.outl, S.grey], getLineWidth: [S.grey, HEW.main.outl] },
    }));
    if (S.selected.size) {
      const sel = [...S.selected].filter((i) => S.vis[i]);
      L.push(new deck.PolygonLayer({ id: "sel-poly", data: sel, getPolygon: (i) => xy.subarray(off[i] * 2, off[i + 1] * 2), positionFormat: "XY",
        stroked: true, filled: false, getLineColor: [255, 255, 255, 255], lineWidthUnits: "pixels", getLineWidth: 2, pickable: false }));
    }
    return withPaths(L);
  }
  S.mode = "pt";
  const pos = new Float32Array(m * 2);
  for (let j = 0; j < m; j++) { pos[j * 2] = S.xy[di[j] * 2]; pos[j * 2 + 1] = -S.xy[di[j] * 2 + 1]; }
  L.push(new deck.ScatterplotLayer({
    id: "cells-pt", data: { length: m, attributes: { getPosition: { value: pos, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true, type: "unorm8" } } },
    getRadius: S.psize / 2, radiusUnits: "common", radiusMinPixels: 0.8, radiusMaxPixels: 60, pickable: true, antialiasing: true,
  }));
  if (S.selected.size) {
    const sel = [...S.selected].filter((i) => S.vis[i]);
    L.push(new deck.ScatterplotLayer({
      id: "sel-pt", data: sel, getPosition: (i) => [S.xy[i * 2], -S.xy[i * 2 + 1]], radiusUnits: "common", getRadius: (S.psize / 2) * 1.15,
      stroked: true, filled: false, getLineColor: [255, 255, 255, 255], lineWidthUnits: "pixels", getLineWidth: 1.5, radiusMinPixels: 3, pickable: false,
    }));
  }
  return withPaths(L);
}

// Transcript overlay (projects imported from Xenium output with transcripts): individual decoded transcripts of a few genes, only for the part in view.
const TXPAL = [[255, 70, 70], [60, 200, 255], [255, 210, 40], [120, 255, 120], [255, 120, 255], [255, 150, 40], [160, 140, 255], [255, 255, 255]];
const TX = { genes: [], warn: "", pos: null, col: null, n: 0, sample: null, key: "", tok: 0, t: null };
function txLayer() {
  if (!S.meta.has_transcripts || !$("txon").checked || !TX.n || TX.sample !== S.sample) return null;
  return new deck.ScatterplotLayer({ id: "tx", data: { length: TX.n, attributes: { getPosition: { value: TX.pos, size: 2 }, getFillColor: { value: TX.col, size: 4, normalized: true, type: "unorm8" } } },
    radiusUnits: "pixels", getRadius: +$("txsize").value, pickable: false });
}
const txRedraw = () => { if (deckgl && S.colors) deckgl.setProps({ layers: layers(S.colors) }); };
function txSoon() { if (!S.meta?.has_transcripts) return; clearTimeout(TX.t); TX.t = setTimeout(txFetch, 250); }
function txNote(msg) { const n = $("txnote"); n.hidden = !msg; n.innerHTML = msg || ""; }
async function txFetch() {
  if (!S.meta.has_transcripts || !deckgl) return;
  if (!$("txon").checked || !TX.genes.length) { TX.n = 0; $("txrow2").hidden = !$("txon").checked; if (!$("txon").checked) txNote(""); txRedraw(); return; }
  $("txrow2").hidden = false;
  const vp = deckgl.getViewports()[0]; if (!vp) return;
  const [bx0, by0, bx1, by1] = vp.getBounds(), pad = 0.1 * Math.max(bx1 - bx0, by1 - by0), r10 = (v) => Math.round(v / 10) * 10;
  const x0 = r10(bx0 - pad), x1 = r10(bx1 + pad), y0 = r10(-(by1 + pad)), y1 = r10(-(by0 - pad)), mx = $("txmax").value;          // deck y is the negated cell y
  const key = [S.sample, TX.genes.join(), x0, y0, x1, y1, mx].join("|"); if (key === TX.key) return; TX.key = key;
  const tok = ++TX.tok, genes = TX.genes.slice();
  const r = await fetch(`/api/sample/${encodeURIComponent(S.sample)}/tx?genes=${encodeURIComponent(genes.join(","))}&x0=${x0}&y0=${y0}&x1=${x1}&y1=${y1}&max=${mx}`);
  if (tok !== TX.tok) return;
  if (!r.ok) { TX.n = 0; txNote("Transcripts: " + (await r.text())); txRedraw(); return; }
  const cnt = JSON.parse(r.headers.get("X-Counts")), a = new Float32Array(await r.arrayBuffer()); if (tok !== TX.tok) return;
  const n = a.length / 3, pos = new Float32Array(2 * n), col = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) { pos[2 * i] = a[3 * i]; pos[2 * i + 1] = -a[3 * i + 1]; const c = TXPAL[a[3 * i + 2] % TXPAL.length]; col[4 * i] = c[0]; col[4 * i + 1] = c[1]; col[4 * i + 2] = c[2]; col[4 * i + 3] = 235; }
  TX.pos = pos; TX.col = col; TX.n = n; TX.sample = S.sample;
  const tot = cnt.in_view.reduce((u, v) => u + v, 0);
  txNote((TX.warn ? TX.warn + "<br>" : "") + cnt.genes.map((g, j) => `<span style="color:rgb(${TXPAL[j % TXPAL.length]})">■</span> ${g}: ${cnt.in_view[j].toLocaleString()} in view (${cnt.total[j].toLocaleString()} in this tissue)`).join("<br>")
    + (n < tot ? `<br>showing ${n.toLocaleString()} of ${tot.toLocaleString()} in view (thinned; zoom in or raise the limit)` : ""));
  txRedraw();
}
function greyInit() { $("greycells").addEventListener("change", (e) => { S.grey = e.target.checked; render(); for (const w of ["gal", "rgn", "glom", "gcmp", "nbh", "rnb"]) HE_RENDER[w](); }); }   // the windows draw their own tiles: redraw the open ones
function txInit() {
  $("txrow").hidden = !S.meta.has_transcripts; if (!S.meta.has_transcripts) return;
  $("txon").addEventListener("change", () => { TX.key = ""; txFetch(); });
  $("txsize").addEventListener("input", txRedraw); $("txmax").addEventListener("change", () => { TX.key = ""; txFetch(); });
  const set = async () => {
    const names = $("txgenes").value.split(/[\s,;]+/).filter(Boolean); if (!names.length) { TX.genes = []; TX.warn = ""; TX.n = 0; txNote(""); txRedraw(); return; }
    const d = await jget("/api/genes/resolve?names=" + encodeURIComponent(names.join(",")));
    const all = d.found.concat(d.excluded || []).slice(0, 8), changed = all.join() !== TX.genes.join();   // genes left out of the signature scores can still be shown as transcripts
    TX.genes = all; TX.key = ""; if (all.length && changed) $("txon").checked = true;                       // only when the list changed: a blur must not undo the user's untick
    TX.warn = (d.missing.length ? `Not on the panel: ${d.missing.join(", ")}` : "") + (d.found.length + (d.excluded || []).length > 8 ? "; only the first 8 genes are shown" : "");
    if (TX.warn) txNote(TX.warn);
    txFetch();
  };
  $("txgenes").addEventListener("change", set); $("txgenes").addEventListener("keydown", (e) => { if (e.key === "Enter") set(); });   // addEventListener: a handler returning false would cancel typing
}

// H&E slide window: the scanned H&E slide (he_slide.json: one pyramid per scanned section) with its annotations, in the slide's own frame (um, y down).
// Not registered to the Xenium coordinates, so it is independent of the cell viewer; tiles come through the same loader / layer function as the registered H&E.
const HS = { idx: null, ann: null, deck: null, vs: null, hide: new Set(), t: null, open: false };
const hesRedo = () => { clearTimeout(HS.t); HS.t = setTimeout(hesDraw, 50); };
const hesRgb = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
function hesFit(box) {
  const el = $("hesdeck"), w = el.clientWidth || 700, h = el.clientHeight || 450, [x0, y0, x1, y1] = box;
  HS.vs = { target: [(x0 + x1) / 2, -(y0 + y1) / 2, 0], zoom: Math.log2(Math.min(w / (x1 - x0), h / (y1 - y0)) * 0.94), _n: (HS.vs?._n || 0) + 1 };
  HS.deck.setProps({ initialViewState: HS.vs }); hesDraw();
}
function hesDraw() {
  if (!HS.deck || $("heswin").hidden || !HS.vs) return;
  const op = +$("hesop").value, L = heLayers(HS.vs, $("hesdeck"), hesRedo, HS.idx.scenes).map((l) => l.clone({ opacity: op }));
  if (HS.ann && $("hesann").checked) {
    const feats = HS.ann.features.filter((f) => !HS.hide.has(f.properties.annotation_name));
    L.push(new deck.GeoJsonLayer({ id: "hes-ann", data: { type: "FeatureCollection", features: feats }, stroked: true, filled: false, pickable: true, lineWidthUnits: "pixels",
      getLineWidth: (f) => (f.properties.negative_roa ? 1.2 : 2.4), getLineColor: (f) => [...hesRgb(f.properties.line_color_hex), f.properties.negative_roa ? 170 : 255] }));
  }
  HS.deck.setProps({ layers: L });
}
function hesPanel() {
  const P = $("hespanel"); P.replaceChildren();
  const row = (sw, text, n, onClick) => { const d = document.createElement("div"); d.className = "hr"; d.innerHTML = `<span class="sw" style="background:${sw}"></span><span>${text}</span><span class="n">${n}</span>`; d.addEventListener("click", onClick); return d; };
  const h = (t) => { const e = document.createElement("h4"); e.textContent = t; return e; };
  P.append(h("Scanned sections (click to zoom)"));
  for (const sc of HS.idx.scenes) P.append(row("#8aa4c8", sc.name, `${(sc.w * sc.px / 1000).toFixed(1)} × ${(sc.h * sc.px / 1000).toFixed(1)} mm`, () => hesFit([sc.x0, sc.y0, sc.x0 + sc.w * sc.px, sc.y0 + sc.h * sc.px])));
  const a = HS.idx.annotations;
  if (a) {
    P.append(h("Annotations (click to show / hide)"));
    for (const n of a.names) {
      const r = row(n.color, n.name, n.n_regions || "—", () => { HS.hide.has(n.name) ? HS.hide.delete(n.name) : HS.hide.add(n.name); r.style.opacity = HS.hide.has(n.name) ? 0.4 : 1; hesDraw(); });
      if (!n.n_regions) r.style.opacity = 0.4; P.append(r);
    }
  }
}
// A batch project can hold several scans (one per slide; several runs may share one): the window then has a scan selector, set to the scan of the tissue being viewed.
async function hesLoadScan(key) {
  const q = key ? "?scan=" + encodeURIComponent(key) : "";
  HS.idx = await jget("/api/heslide" + q); HS.ann = null; HS.hide = new Set(); HS.key = key || "";
  if (HS.idx.annotations) { try { HS.ann = await jget("/api/heslide/annotations" + q); for (const f of HS.ann.features) f.geometry.coordinates = f.geometry.coordinates.map((r) => r.map(([x, y]) => [x, -y])); } catch (e) { HS.ann = null; } }   // deck y is the negated slide y
  $("hestitle").textContent = `H&E slide — ${HS.idx.source}`; hesPanel();
  if (HS.deck) hesFit(HS.idx.bounds_um);
}
function hesScanOfSample() {   // the scan of the run of the current tissue in a multi-scan project ("" when that run was given none)
  const m = HS.multi, sm = S.meta.samples.find((x) => x.id === S.sample), run = sm ? subtypeOf(sm) : "";
  return m.runs[run] || "";
}
function hesNoScanNote(k) { if (!k) { const sm = S.meta.samples.find((x) => x.id === S.sample); $("hesnote").textContent = `No H&E scan was given for ${sm ? subtypeOf(sm) : "this run"}; showing ${HS.idx.source}.`; } }
async function hesToggle(on) {
  HS.open = on ?? !HS.open; $("heswin").hidden = !HS.open; $("hesbtn").classList.toggle("on", HS.open); if (!HS.open) return;
  if (!HS.idx) {
    const d = await jget("/api/heslide");
    if (d.multi) {
      HS.multi = d; const sel = $("hesscan"); sel.replaceChildren();
      for (const sc of d.scans) { const o = document.createElement("option"); o.value = sc.key; o.textContent = `${sc.source} · ${sc.runs.length} run${sc.runs.length === 1 ? "" : "s"}`; o.title = sc.runs.join("\n"); sel.append(o); }
      sel.hidden = false; const k0 = hesScanOfSample(); sel.value = k0 || d.scans[0].key; await hesLoadScan(sel.value); hesNoScanNote(k0);
    } else {
      HS.idx = d; HS.ann = null;
      if (HS.idx.annotations) { try { HS.ann = await jget("/api/heslide/annotations"); for (const f of HS.ann.features) f.geometry.coordinates = f.geometry.coordinates.map((r) => r.map(([x, y]) => [x, -y])); } catch (e) { HS.ann = null; } }
      $("hestitle").textContent = `H&E slide — ${HS.idx.source}`; hesPanel();
    }
    HS.deck = new deck.Deck({ parent: $("hesdeck"), views: new deck.OrthographicView({ flipY: false }), initialViewState: { target: [0, 0, 0], zoom: -4 }, controller: { dragPan: true, scrollZoom: true, doubleClickZoom: true }, layers: [],
      onViewStateChange: ({ viewState }) => { HS.vs = viewState; hesRedo(); },
      onHover: (i) => { const f = i.object?.properties; $("hesnote").textContent = f ? `${f.annotation_name}${f.negative_roa ? " (negative region: excluded from the annotation)" : ""} · ${f.region_id}` : ""; },
      getCursor: ({ isDragging }) => (isDragging ? "grabbing" : "grab") });
    hesFit(HS.idx.bounds_um);
  } else {
    if (HS.multi) { const k = hesScanOfSample(); if (k && k !== HS.key) { $("hesscan").value = k; await hesLoadScan(k); return; } hesNoScanNote(k); }   // reopened on another tissue: show its slide
    hesDraw();
  }
}
function hesInit() {
  $("hesbtn").hidden = !S.meta.has_he_slide; if (!S.meta.has_he_slide) return;
  initFloating($("heswin"), $("heshead"), "xb_hes");
  $("hesbtn").addEventListener("click", () => hesToggle()); $("hesx").addEventListener("click", () => hesToggle(false));
  $("hesscan").addEventListener("change", () => hesLoadScan($("hesscan").value));
  $("hesfit").addEventListener("click", () => HS.idx && hesFit(HS.idx.bounds_um));
  $("hesann").addEventListener("change", hesDraw); $("hesop").addEventListener("input", hesDraw);
  if (window.ResizeObserver) new ResizeObserver(() => { if (HS.open) hesRedo(); }).observe($("hesdeck"));
}

// H&E beside the cells in the main viewer: a second deck.gl view of the H&E only, optionally locked to the main view; hovering one marks the same spot in the other.
// ---- Pointers: arrows the user places on the tissue. They are stored in the tissue's own frame (the deck coordinates of the cell view), so the same arrow is drawn in the cell view and in the
// side pane whenever that pane shows a registered H&E or IF image (they share the cells' frame; it does not matter whether zoom is linked). Drawn as pixel-sized icons: the size does not change with zoom.
const PT = { on: false, list: {}, hist: [] };       // hist = the order arrows were placed in (Undo removes the last one, in any tissue)
const PT_SVG = '<svg xmlns="http://www.w3.org/2000/svg" width="48" height="48" viewBox="0 0 48 48"><g transform="rotate(-45 24 24)"><path d="M0 24 L16 10 L16 19 L46 19 L46 29 L16 29 L16 38 Z" fill="#e11d48" stroke="#ffffff" stroke-width="2" stroke-linejoin="round"/></g></svg>';
const PT_ICON = { url: "data:image/svg+xml;charset=utf-8," + encodeURIComponent(PT_SVG), width: 48, height: 48, anchorX: 7, anchorY: 41 };   // the arrow's tip is at (7, 41) of the icon
const ptList = () => (PT.list[S.sample] ||= []);
const ptAllowedFor = (w) => !!(HEW[w] && HEW[w].side && !(w === "main" && scanMode()) && ["he", "if"].some((k) => HEW[w].sides.has(k)));   // that window's side pane shows an image in the cells' frame
const ptAllowed = () => ptAllowedFor("main");
function ptPush(sample, x, y) { (PT.list[sample] ||= []).push({ x, y }); PT.hist.push(sample); }
function ptUndo() { const smp = PT.hist.pop(); if (smp && PT.list[smp]) PT.list[smp].pop(); }
function ptLayer() { const l = PT.list[S.sample]; return l && l.length ? new deck.IconLayer({ id: "pointers", data: l.slice(), getIcon: () => PT_ICON, getPosition: (d) => [d.x, d.y, 0], getSize: 42, sizeUnits: "pixels", pickable: false }) : null; }
function ptDraw() { if (S.colors && deckgl) deckgl.setProps({ layers: layers(S.colors) }); heSideLayers(); }
function ptClick(info) { if (!PT.on || !ptAllowed() || !info || !info.coordinate) return false; ptPush(S.sample, info.coordinate[0], info.coordinate[1]); ptUi(); ptDraw(); return true; }
// The other windows draw tissue tiles on 2-D canvases (cv._xf: canvas px = s * um + offset, y down), so the arrows are drawn from that transform after each tile is drawn, and a click on
// either canvas of a side-by-side pair places an arrow. The arrows are shared per tissue: one placed in any window shows in the others.
const PT_IMG = new Image(); PT_IMG.src = PT_ICON.url;
function ptCanvasDraw(cv, sample) {
  const l = PT.list[sample]; if (!l || !l.length || !cv._xf || !PT_IMG.complete) return;
  const ctx = cv.getContext("2d"), f = cv._xf, k = (34 / 48) * (cv.width / (cv.clientWidth || 700));          // about 34 screen px
  for (const p of l) ctx.drawImage(PT_IMG, f.s * p.x + f.ox - PT_ICON.anchorX * k, f.s * -p.y + f.oy - PT_ICON.anchorY * k, 48 * k, 48 * k);
}
function ptBindCanvas(el, w) {
  el.addEventListener("pointerdown", (e) => { el._pd = [e.clientX, e.clientY]; });
  el.addEventListener("click", (e) => {
    if (!PT.on || !el._xf || !el._ptSample || !ptAllowedFor(w)) return;
    if (el._pd && Math.hypot(e.clientX - el._pd[0], e.clientY - el._pd[1]) > 4) return;                          // a drag (pan), not a click
    e.stopPropagation(); e.preventDefault();
    const r = el.getBoundingClientRect(), px = (e.clientX - r.left) * el.width / r.width, py = (e.clientY - r.top) * el.height / r.height, f = el._xf;
    ptPush(el._ptSample, (px - f.ox) / f.s, -((py - f.oy) / f.s)); ptUi(); ptDraw(); HE_RENDER[w] && HE_RENDER[w]();
  }, true);
}
function ptWrap(fn, sampleOf) {   // draw the arrows on the canvases of a side-by-side pair after the tile is drawn
  return function (...a) { const r = fn.apply(this, a), cv = a[0]; try { if (cv && (cv._heOnly || (cv.parentNode && cv.parentNode.classList.contains("hepair")))) ptCanvasDraw(cv, sampleOf(a)); } catch (e) { /* no arrows */ } return r; };
}
function ptUi() {
  const row = $("ptrow"); if (!row) return; const ok = ptAllowed(); if (!ok && !HE_SIDE.some((k) => ptAllowedFor(k))) PT.on = false;   // the tool stays on while some window can use it
  if (row.hidden === ok) row.hidden = !ok;
  $("ptadd").classList.toggle("on", PT.on); document.querySelectorAll('.ptwbtn[data-act="tool"]').forEach((b) => b.classList.toggle("on", PT.on)); const n = (PT.list[S.sample] || []).length, t = n ? `${n} placed` : ""; if ($("ptn").textContent !== t) $("ptn").textContent = t;
}
function ptInit() {
  $("ptadd").addEventListener("click", () => { PT.on = !PT.on; ptUi(); });
  $("ptundo").addEventListener("click", () => { ptUndo(); ptUi(); ptDraw(); for (const k of HE_SIDE) HEW[k].side && HE_RENDER[k] && HE_RENDER[k](); });
  $("ptclr").addEventListener("click", () => { PT.list[S.sample] = []; PT.hist = PT.hist.filter((x) => x !== S.sample); ptUi(); ptDraw(); for (const k of HE_SIDE) HEW[k].side && HE_RENDER[k] && HE_RENDER[k](); });
  // the same tool in the header of every window that can show an image beside the cells
  document.addEventListener("click", (e) => {
    const b = e.target.closest(".ptwbtn"); if (!b) return; const w = b.dataset.w;
    if (b.dataset.act === "tool") PT.on = !PT.on;
    else if (b.dataset.act === "undo") ptUndo();
    else { PT.list = {}; PT.hist = []; }
    ptUi(); ptDraw(); for (const k of HE_SIDE) HEW[k].side && HE_RENDER[k] && HE_RENDER[k]();
  });
}
// ---- Slide scan (unregistered H&E) beside the cells: the side pane shows the whole scan with its own pan / zoom (nothing links it to the cells); the
// orientation buttons turn / flip it by transforming each tile's corners (the scan's own frame: micrometres from its top-left, y down).
const SC = { idx: null, multi: null, key: "", sample: null, orient: {}, loading: null };
const scanMode = () => HEW.main.sides.has("scan");
const scOri = () => (SC.orient[SC.key] ||= { rot: 0, flip: false });
function scXf(o, W, H, x, y) { if (o.flip) x = W - x; let w = W, h = H; for (let i = 0; i < o.rot; i++) { [x, y] = [h - y, x]; [w, h] = [h, w]; } return [x, y]; }   // scan frame -> displayed frame (mirror, then quarter turns clockwise)
function scIxf(o, W, H, x, y) {   // displayed frame -> scan frame
  const dims = [[W, H]]; for (let i = 0; i < o.rot; i++) { const [w, h] = dims[i]; dims.push([h, w]); }
  for (let i = o.rot - 1; i >= 0; i--) { const h = dims[i][1]; [x, y] = [y, h - x]; }
  return [o.flip ? W - x : x, y];
}
function scanLayers(vs, el, redo) {
  const idx = SC.idx; if (!idx || !vs) return [];
  const o = scOri(), W = idx.bounds_um[2], H = idx.bounds_um[3], k = 2 ** vs.zoom, hw = (el.clientWidth || 800) / k / 2, hh = (el.clientHeight || 600) / k / 2, cx = vs.target[0], cy = vs.target[1];
  const q = [[cx - hw, cy - hh], [cx + hw, cy - hh], [cx + hw, cy + hh], [cx - hw, cy + hh]].map(([x, y]) => scIxf(o, W, H, x, -y));   // the view's corners in the scan frame
  const vx0 = Math.min(...q.map((a) => a[0])), vx1 = Math.max(...q.map((a) => a[0])), vy0 = Math.min(...q.map((a) => a[1])), vy1 = Math.max(...q.map((a) => a[1]));
  const P = (x, y) => { const [a, b] = scXf(o, W, H, x, y); return [a, -b]; }, out = [];
  for (const im of idx.scenes) {
    const T = im.tsize, u = im.px, z = Math.max(0, Math.min(im.nz - 1, Math.floor(Math.log2(1 / (u * k))))), span = T * u * 2 ** z;
    const t0 = Math.max(0, Math.floor((vx0 - im.x0) / span)), t1 = Math.min(Math.ceil(im.w * u / span) - 1, Math.floor((vx1 - im.x0) / span));
    const r0 = Math.max(0, Math.floor((vy0 - im.y0) / span)), r1 = Math.min(Math.ceil(im.h * u / span) - 1, Math.floor((vy1 - im.y0) / span));
    const layer = (zz, tx, ty, img) => { const f = u * 2 ** zz, sp = T * f, X0 = im.x0 + tx * sp, Y0 = im.y0 + ty * sp, w = img.naturalWidth * f, h = img.naturalHeight * f;
      return { id: `sc-${im.key}-${zz}-${tx}-${ty}`, image: img, bounds: [P(X0, Y0 + h), P(X0, Y0), P(X0 + w, Y0), P(X0 + w, Y0 + h)], pickable: false }; };   // corners: bottom-left, top-left, top-right, bottom-right of the picture
    const back = new Map(), front = [];
    for (let ty = r0; ty <= r1; ty++) for (let tx = t0; tx <= t1; tx++) {
      const img = heTile(im.key, z, tx, ty, redo);
      if (img) { front.push(layer(z, tx, ty, img)); continue; }
      for (let zz = z + 1; zz < im.nz; zz++) { const d = zz - z, pimg = heTile(im.key, zz, tx >> d, ty >> d, null); if (pimg) { const c = layer(zz, tx >> d, ty >> d, pimg); back.set(c.id, c); break; } }
    }
    out.push(...back.values(), ...front);
  }
  return out.map((c) => new deck.BitmapLayer(c));
}
function scanFit() {
  if (!HEX.deck || !SC.idx) return;
  const o = scOri(), W = SC.idx.bounds_um[2], H = SC.idx.bounds_um[3], odd = o.rot % 2 === 1, w = odd ? H : W, h = odd ? W : H, el = $("hewrap");
  HEX.vs = { target: [w / 2, -h / 2, 0], zoom: Math.log2(Math.min((el.clientWidth || 600) / w, (el.clientHeight || 500) / h) * 0.94), _n: ++HEX.n };
  HEX.deck.setProps({ initialViewState: HEX.vs }); heSideLayers();
}
async function scanLoad() {   // the scan of the tissue being viewed (a batch project can hold several: the one of this tissue's run, else the first)
  const d = await jget("/api/heslide");
  if (d.multi) {
    SC.multi = d; const sm = S.meta.samples.find((x) => x.id === S.sample), run = sm ? subtypeOf(sm) : "", key = d.runs[run] || SC.key || d.scans[0].key, sel = $("scansel");
    sel.replaceChildren(); for (const sc of d.scans) { const op = document.createElement("option"); op.value = sc.key; op.textContent = `${sc.source} · ${sc.runs.length} run${sc.runs.length === 1 ? "" : "s"}`; sel.append(op); }
    sel.value = key; if (key !== SC.key || !SC.idx) { SC.key = key; SC.idx = await jget("/api/heslide?scan=" + encodeURIComponent(key)); SC.sample = S.sample; scanUi(); return scanFit(); }
  } else if (!SC.idx) { SC.idx = d; SC.key = ""; SC.sample = S.sample; scanUi(); return scanFit(); }
  SC.sample = S.sample; scanUi(); heSideLayers();
}
function scanUi() {
  const on = scanMode() && HEW.main.side; $("scanrow").hidden = !on; $("scansel").hidden = !(on && SC.multi);
  $("scannote").hidden = !on; $("scannote").textContent = SC.idx ? `${SC.idx.source}: not aligned to the cells; pan and zoom here are independent of the cell view.` : "loading the scan…";
}
function scanInit() {
  const turn = (f) => () => { if (!SC.idx) return; f(scOri()); scanFit(); };
  $("scrotl").addEventListener("click", turn((o) => (o.rot = (o.rot + 3) % 4))); $("scrotr").addEventListener("click", turn((o) => (o.rot = (o.rot + 1) % 4)));
  $("scflipx").addEventListener("click", turn((o) => (o.flip = !o.flip)));
  $("scflipy").addEventListener("click", turn((o) => { o.flip = !o.flip; o.rot = (o.rot + 2) % 4; }));   // top-bottom mirror = left-right mirror + half turn
  $("scfit").addEventListener("click", () => scanFit());
  $("scansel").addEventListener("change", async () => { SC.key = $("scansel").value; SC.idx = await jget("/api/heslide?scan=" + encodeURIComponent(SC.key)); scanFit(); });
}
const HEX = { deck: null, vs: null, n: 0 };
function heSideLayers() {
  if (!HEX.deck) return;
  if (HEW.main.side && scanMode()) {                                                    // the slide scan instead of the registered images
    if (!SC.idx || (SC.multi && SC.sample !== S.sample)) { if (!SC.loading) SC.loading = scanLoad().finally(() => (SC.loading = null)); return; }
    if (!HEX.vs) return scanFit();
    HEX.deck.setProps({ layers: scanLayers(HEX.vs, $("hewrap"), () => { clearTimeout(HEX.t); HEX.t = setTimeout(heSideLayers, 60); }) }); return;
  }
  const on = HEW.main.side && HEX.vs;
  HEX.deck.setProps({ layers: on ? [...heLayers(HEX.vs, $("hewrap"), () => { clearTimeout(HEX.t); HEX.t = setTimeout(heSideLayers, 60); }, heSel(S.sample, "main", true)), ...(ptLayer() ? [ptLayer()] : [])] : [] });
}
function heSideView(vs) { if (!HEX.deck || !vs) return; HEX.vs = vs; HEX.deck.setProps({ initialViewState: { ...vs, _n: ++HEX.n } }); heSideLayers(); }   // main -> H&E (linked)
function heMainView(vs) {   // H&E -> main (linked)
  S.viewState = vs; S.viewN = (S.viewN || 0) + 1; deckgl.setProps({ initialViewState: { ...vs, _n: S.viewN } }); updateScaleBar();
  if (S.colors) { clearTimeout(S.vt); S.vt = setTimeout(() => deckgl.setProps({ layers: layers(S.colors) }), 120); }
}
function heSideLayout() {
  const on = HEW.main.side; $("viewer").classList.toggle("he2", on);
  if (on && !HEX.deck) {
    HEX.deck = new deck.Deck({
      parent: $("hewrap"), views: new deck.OrthographicView({ flipY: false }), initialViewState: S.viewState || { target: [0, 0, 0], zoom: 0 }, controller: { dragPan: true, scrollZoom: true, doubleClickZoom: true }, layers: [],
      onViewStateChange: ({ viewState }) => { HEX.vs = viewState; if (HEW.main.link && !scanMode()) heMainView(viewState); heSideLayers(); },
      onHover: (info) => hePtr(info, "he"), onClick: (info) => ptClick(info), getCursor: ({ isDragging }) => (isDragging ? "grabbing" : PT.on && ptAllowed() ? "crosshair" : "default"),
    });
    $("hewrap").addEventListener("mouseleave", () => hePtr(null)); $("deckwrap").addEventListener("mouseleave", () => hePtr(null));
  }
  if (!on) { hePtr(null); heSideLayers(); return; }
  setTimeout(() => { window.dispatchEvent(new Event("resize")); if (scanMode()) heSideLayers(); else if (!HEX.vs || HEW.main.link) heSideView(S.viewState); else heSideLayers(); }, 30);   // the viewer just changed width
}
function hePtr(info, from) {   // mark the point under the cursor of one pane in the other
  const m = [$("heptr1"), $("heptr2")]; m.forEach((e) => (e.hidden = true));
  if (!info || !info.coordinate || !HEW.main.side || !HEX.deck || !deckgl || scanMode()) return;   // no hover marker for the unaligned scan
  const dst = from === "he" ? deckgl : HEX.deck, vp = dst.getViewports()[0], el = from === "he" ? $("deckwrap") : $("hewrap"); if (!vp) return;
  const [x, y] = vp.project([info.coordinate[0], info.coordinate[1], 0]), mk2 = from === "he" ? m[0] : m[1];
  if (x < 0 || y < 0 || x > el.clientWidth || y > el.clientHeight) return;
  mk2.style.left = x + "px"; mk2.style.top = y + "px"; mk2.hidden = false;
}
function render() {
  computeVisible();
  const colors = S.colors = computeColors();
  const props = { layers: layers(colors) };
  if (!deckgl) {
    deckgl = new deck.Deck({
      parent: $("deckwrap"), views: new deck.OrthographicView({ flipY: false }), initialViewState: S.viewState,
      controller: { dragPan: true, scrollZoom: true, doubleClickZoom: true }, ...props,
      onViewStateChange: ({ viewState }) => {
        S.viewState = viewState; updateScaleBar(); txSoon(); if (HEW.main.side && HEW.main.link && !scanMode()) heSideView(viewState);
        if (HEW.main.on && S.colors) { clearTimeout(S.heV); S.heV = setTimeout(() => deckgl.setProps({ layers: layers(S.colors) }), 90); }
        if (S.shape === "centroid" || !S.colors) return;
        // polygon mode depends on the viewport, so refresh (debounced) while in it or when entering it
        if (S.poly && (S.mode === "poly" || S.shape === "polygon" || viewState.zoom >= S.polyZoom)) {
          clearTimeout(S.vt); S.vt = setTimeout(() => deckgl.setProps({ layers: layers(S.colors) }), 120);
        }
      },
      onHover: onHover, onClick: onClick, getCursor: ({ isDragging }) => (S.boxMode || (PT.on && ptAllowed()) ? "crosshair" : isDragging ? "grabbing" : "default"),
    });
  } else deckgl.setProps(props);
  if (HEW.main.side) heSideLayers();
  txSoon();
}

function updateScaleBar() {
  if (!S.viewState) return;
  const ppu = Math.pow(2, S.viewState.zoom);  // pixels per µm
  const target = 110 / ppu, p = Math.pow(10, Math.floor(Math.log10(target)));
  const nice = [1, 2, 5, 10].map((k) => k * p).reduce((a, b) => (Math.abs(b - target) < Math.abs(a - target) ? b : a));
  $("sbline").style.width = nice * ppu + "px";
  $("sblabel").textContent = nice >= 1000 ? nice / 1000 + " mm" : nice + " µm";
}

// ---------------------------------------------------------------- interaction
function tipHTML(i) {
  const rows = [`<b>${S.ids[i]}</b>`];
  for (const k of ROLE_ORDER) if (S.cats[k]) { const c = S.cats[k][i]; if (k === "glom_quality" && c === MISSING) continue; rows.push(`${ROLE_LABEL[k]}: ${c === MISSING ? "<i>missing</i>" : display(k, c)}`); }
  if (S.cats.path_region) { const v = S.cats.path_region[i], lv = catMeta("path_region").levels, r = lv.filter((_, b) => (v >> b) & 1 && b < lv.length - 1); rows.push(`Pathologist regions: ${r.length ? r.join(", ") : "<i>outside all outlines</i>"}`); }
  if (S.cats.glom_id && S.cats.glom_id[i] !== MISSING) rows.push(`Glomerulus: ${catMeta("glom_id").levels[S.cats.glom_id[i]]}`);
  for (const [g, v] of geneList()) rows.push(`${g}: ${v[i] !== v[i] ? "n/a (outside the compute-over subset)" : v[i].toFixed(2)} (${S.layer === "counts" ? "counts" : "log-norm"}${S.colorBy === "__sig" || S.colorBy === "__sigdual" ? ", mean" : ""}${S.zscope !== "none" ? ", z-score" : ""}${S.smooth > 0 ? `, smoothed σ=${S.smooth}µm` : ""})`);
  return rows.join("<br>");
}
function cellAt(info) {  // original cell index for a pick on the "cells" layer (scatter: index; polygons: object = j)
  if (!(info.index >= 0) || !/^cells-/.test(info.layer?.id || "")) return -1;
  return S.drawIdx[S.mode === "poly" ? info.object : info.index];
}
function onHover(info) {
  hePtr(info, "main");
  const t = $("tip"), i = cellAt(info);
  if (i >= 0) { t.innerHTML = tipHTML(i); t.style.display = "block";
    t.style.left = Math.min(info.x + 14, $("viewer").clientWidth - 280) + "px"; t.style.top = Math.max(0, info.y + 14) + "px";
  } else t.style.display = "none";
}
async function onClick(info, ev) {
  if (ptClick(info)) return;                                          // the pointer tool is on: place an arrow here instead of selecting a cell
  const i = cellAt(info);
  if (S.boxMode || i < 0) return;
  if (ev?.srcEvent?.shiftKey) S.selected.has(i) ? S.selected.delete(i) : S.selected.add(i); else S.selected = new Set([i]);
  setHighlight(S.sample, S.selected); await showCell(i); updateAll();
}

async function showCell(i) {
  const el = $("cell"); el.innerHTML = "loading…";
  try {
    const r = await jget(`/api/cell/${encodeURIComponent(S.sample)}/${i}`);
    const md = r.metadata;
    const core = Object.entries(md).filter(([k, v]) => v !== null && v !== "").map(([k, v]) => `<tr><td class="muted">${k}</td><td>${typeof v === "number" ? +v.toFixed(3) : v}</td></tr>`).join("");
    const ex = r.expression.map((e) => `<tr class="g"><td>${e.gene}</td><td>${e.counts}</td><td>${e.logcpx.toFixed(2)}</td></tr>`).join("");
    el.innerHTML = `<div><b>${S.ids[i]}</b></div><div class="small muted">x=${S.xy[i * 2].toFixed(1)} µm, y=${S.xy[i * 2 + 1].toFixed(1)} µm</div>
      <details open><summary>Metadata</summary><table class="expr">${core}</table></details>
      <details open><summary>Expression (${r.n_genes_detected} genes detected)</summary><input placeholder="filter genes…" id="gfilt">
      <table class="expr"><tr><th>gene</th><th>counts</th><th>log-norm</th></tr>${ex}</table></details>`;
    $("gfilt").addEventListener("input", (e) => { const q = e.target.value.toLowerCase(); el.querySelectorAll("tr.g").forEach((tr) => { tr.style.display = tr.firstChild.textContent.toLowerCase().includes(q) ? "" : "none"; }); });
  } catch (e) { el.textContent = "Error: " + e.message; }
}

// box selection
function bindBox() {
  const v = $("viewer"), r = $("boxrect"); let s = null;
  v.addEventListener("mousedown", (e) => { if (!S.boxMode) return; e.stopPropagation(); const b = v.getBoundingClientRect(); s = [e.clientX - b.left, e.clientY - b.top]; r.style.display = "block"; }, true);
  window.addEventListener("mousemove", (e) => {
    if (!s) return; const b = v.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
    Object.assign(r.style, { left: Math.min(s[0], x) + "px", top: Math.min(s[1], y) + "px", width: Math.abs(x - s[0]) + "px", height: Math.abs(y - s[1]) + "px" });
  });
  window.addEventListener("mouseup", (e) => {
    if (!s) return; const b = v.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top;
    const vp = deckgl.getViewports()[0], p1 = vp.unproject([s[0], s[1]]), p2 = vp.unproject([x, y]);
    s = null; r.style.display = "none";
    const x0 = Math.min(p1[0], p2[0]), x1 = Math.max(p1[0], p2[0]), y0 = Math.min(p1[1], p2[1]), y1 = Math.max(p1[1], p2[1]);
    const add = e.shiftKey, sel = add ? new Set(S.selected) : new Set();
    for (const i of S.visIdx) { const px = S.xy[i * 2], py = -S.xy[i * 2 + 1]; if (px >= x0 && px <= x1 && py >= y0 && py <= y1) sel.add(i); }
    S.selected = sel; setHighlight(S.sample, sel); updateAll();
  });
  v.addEventListener("mousedown", (e) => { if (S.boxMode) e.preventDefault(); });
}

// ---------------------------------------------------------------- panels
function countsFor(k, subset) {
  const c = S.cats[k], n = catMeta(k).levels.length, out = new Array(n).fill(0); let miss = 0;
  const f = (i) => { const v = c[i]; v === MISSING ? miss++ : out[v]++; };
  subset ? subset.forEach(f) : S.visIdx.forEach(f);
  return { out, miss };
}
function updateStats() {
  S.counts = {}; for (const k of ROLE_ORDER) if (S.cats[k]) S.counts[k] = countsFor(k).out;
  if (S.cats.path_region) { const m = S.cats.path_region, nl = catMeta("path_region").levels.length, out = new Array(nl).fill(0);
    for (const i of S.visIdx) { const v = m[i]; for (let b = 0; b < nl; b++) if ((v >> b) & 1) out[b]++; } S.counts.path_region = out; }
  const nsel = [...S.selected].filter((i) => S.vis[i]).length;
  const rows = [["Cells in dataset", fmt(S.meta.n_cells)], ["Cells in this donor tissue", fmt(S.n)], ["Passing filters", fmt(S.visIdx.length)], ["Selected", fmt(nsel)]];
  const cur = S.meta.samples.find((x) => x.id === S.sample);
  if (cur && cur.tx_tissue_total != null) rows.splice(2, 0, ["Transcripts in this tissue", fmt(cur.tx_tissue_total)]);   // rows of the transcripts file that belong to this tissue (all qualities)
  if (cur && cur.tx_slide_total) rows.splice(cur && cur.tx_tissue_total != null ? 3 : 2, 0, ["Transcripts in this slide", fmt(cur.tx_slide_total)]);   // all rows of the slide's transcripts file (the run this tissue is from)
  for (const k of ROLE_ORDER) if (S.cats[k]) rows.push([`Distinct ${ROLE_LABEL[k].toLowerCase()}`, fmt(S.counts[k].filter((x) => x > 0).length)]);
  rows.push(["Genes available", fmt(S.meta.n_genes)]);
  for (const [g, v] of geneList()) {
    let e = 0, sum = 0, nv = 0; for (const i of S.visIdx) { if (v[i] > 0) e++; if (v[i] === v[i]) { sum += v[i]; nv++; } }
    const lay = S.layer === "counts" ? "raw counts" : `log1p(CP${fmt(S.meta.norm_scale)})`;
    if (transformed()) {
      const n = Math.max(1, nv), mu = sum / n; let ss = 0; for (const i of S.visIdx) if (v[i] === v[i]) ss += (v[i] - mu) ** 2;
      rows.push([`${g}: mean ± sd, ${withinSpec() !== "tissue" ? `${fmt(nv)} of ${fmt(S.visIdx.length)} filtered cells (those inside the subset)` : "filtered cells"}`, `${mu.toFixed(3)} ± ${Math.sqrt(ss / n).toFixed(3)} <span class="muted">${S.zscope !== "none" ? "z" : lay}</span>`]); continue;
    }
    rows.push([`${g}: expressing (>0)`, `${fmt(e)} (${(100 * e / Math.max(1, S.visIdx.length)).toFixed(1)}%)`]);
    rows.push([`${g}: mean, all filtered cells`, `${(sum / Math.max(1, S.visIdx.length)).toFixed(3)} <span class="muted">${lay}</span>`]);
  }
  $("stats").innerHTML = rows.map(([a, b]) => `<tr><td>${a}</td><td>${b}</td></tr>`).join("");
}

function scaleLegendHTML(ci) {  // legend for continuous scales (gene / two genes / signature / numeric)
  const lay = S.layer === "counts" ? "raw counts" : `log1p(counts/total×${fmt(S.meta.norm_scale)})`;
  const grad = (L) => `<div class="grad" style="background:linear-gradient(90deg,${L.filter((_, i) => i % 32 === 0).map((c) => `rgb(${c})`).join(",")})"></div>`;
  const note = transformed() ? `<div class="small muted">${transformNote()}</div>` : "";
  if (ci.type === "gene" && ci.cont) return `<div>${ci.label} — ${ci.sig ? "mean of " : ""}${lay}</div>${grad(ci.div ? DIVL : LUT)}${tickRowHTML(scaleTicks(ci))}
      <div class="small muted">${ci.div ? "z-score, 0 = average · " : ""}ends: ${ci.div ? "±" : "0 and "}p${S.clip}</div>${note}`;
  if (ci.type === "sigdual") return `<div class="small"><span style="color:${SIG_A}">■</span> ${ci.labelA} &nbsp; ${ci.labelB ? `<span style="color:${SIG_B}">■</span> ${ci.labelB} &nbsp; <span style="color:#ffe0ff">■</span> both high` : "(choose signature B)"}</div>
      <div class="small muted">Intensity = score ${ci.z ? "above the average" : "from the 1st percentile"} up to p${S.clip}, scaled separately for each signature; grey = below. ${lay}${transformed() ? ". " + transformNote() : ""}</div>`;
  if (ci.type === "gene") return `<div>${ci.label} — ${ci.sig ? "mean of " : ""}${lay}</div>${grad(LUT)}${tickRowHTML(scaleTicks(ci))}
      <div class="small muted">0 = gray: ${ci.sig ? "no signature gene" : "not"} detected · top = p${S.clip}</div>${note}`;
  if (ci.type === "dual") return `<div><span style="color:#d33">■</span> ${S.geneA} (to ${ci.hiA.toFixed(2)}) &nbsp; <span style="color:#3a3">■</span> ${S.geneB} (to ${ci.hiB.toFixed(2)}) &nbsp; <span style="color:#cc0">■</span> both</div>
      <div class="small muted">${lay}; gray = neither detected; each gene clipped at p${S.clip}${transformed() ? "; positive part shown. " + transformNote() : ""}</div>`;
  if (ci.type === "num") return `<div>${labelOf(S.colorBy)}</div>${grad(LUT)}${tickRowHTML(scaleTicks(ci))}<div class="small muted">${numLegend(ci)}</div>`;
  return `<div class="muted">${S.colorBy === "__sig" ? "Enter genes and press “Apply signature”." : "Search a gene in “Gene A”."}</div>`;
}
// Entries for a categorical legend; for cell states they are grouped under their cell type (same order as the filter list).
function legendEntries(mode, order) {
  if (mode !== "cell_state" || !S.stateInfo) return order.map(([i, c]) => ({ i, c }));
  const out = [], by = new Map(order.map(([i, c]) => [i, c]));
  for (const g of S.stateInfo.groups) { const m = g.levels.filter((i) => by.has(i)); if (m.length) { out.push({ head: g.name }); m.forEach((i) => out.push({ i, c: by.get(i) })); } }
  return out;
}
const legendRowHTML = (mode, e, hid) => e.head ? `<div class="ghead">${e.head}</div>`
  : `<div class="opt ${hid.has(e.i) ? "off" : ""}" data-i="${e.i}" title="click: grey out/show · double-click: highlight only this"><span class="sw" style="background:${catColor(mode, e.i)}"></span><span class="l">${display(mode, e.i)}</span><span class="n">${fmt(e.c)}</span></div>`;
function updateLegend() {
  const el = $("legend"), mode = S.colorBy, ci = S.colorInfo;
  $("legtitle").textContent = "— " + ({ __gene: "gene", __dual: "two genes", __sig: "signature", __sigdual: "two signatures" }[mode] || labelOf(mode));
  if (ci.type !== "cat") { el.innerHTML = scaleLegendHTML(ci); return; }
  const q = $("legsearch").value.toLowerCase(), counts = countsAll(mode), hid = S.hiddenCats[mode] || new Set();
  const order = counts.map((c, i) => [i, c]).filter(([, c]) => c > 0).sort((a, b) => b[1] - a[1]);
  el.innerHTML = legendEntries(mode, order.filter(([i]) => display(mode, i).toLowerCase().includes(q))).map((e) => legendRowHTML(mode, e, hid)).join("");
  el.querySelectorAll(".opt").forEach((o) => {
    o.addEventListener("click", () => { const h = S.hiddenCats[mode] || (S.hiddenCats[mode] = new Set()); const i = +o.dataset.i; h.has(i) ? h.delete(i) : h.add(i); refresh(); });
    o.addEventListener("dblclick", () => { S.hiddenCats[mode] = new Set(order.map(([i]) => i).filter((i) => i !== +o.dataset.i)); refresh(); });
  });
}
// Legend floating over the viewer (toggle in the toolbar): same content and click-to-grey behaviour as the side-panel legend.
function updateOvLegend() {
  const el = $("ovleg"); el.hidden = !S.ovLegend; if (!S.ovLegend) return;
  const mode = S.colorBy, ci = S.colorInfo; $("ovtitle").textContent = colorName();
  const body = $("ovbody"), pl = pathLegendHTML(), plBlock = pl ? `<div class="ghead">Pathologist regions</div>${pl}` : "";
  if (ci.type !== "cat") { body.innerHTML = scaleLegendHTML(ci) + plBlock; return; }
  const counts = countsAll(mode), hid = S.hiddenCats[mode] || new Set();
  const order = counts.map((c, i) => [i, c]).filter(([, c]) => c > 0).sort((a, b) => b[1] - a[1]);
  const f = S.n - S.visIdx.length;
  body.innerHTML = legendEntries(mode, order).map((e) => legendRowHTML(mode, e, hid)).join("")
    + (f ? `<div class="small muted" style="margin-top:6px">${fmt(f)} cells filtered out (grey)</div>` : "") + plBlock;
  body.querySelectorAll(".opt[data-i]").forEach((o) => {
    o.addEventListener("click", () => { const h = S.hiddenCats[mode] || (S.hiddenCats[mode] = new Set()); const i = +o.dataset.i; h.has(i) ? h.delete(i) : h.add(i); refresh(); });
    o.addEventListener("dblclick", () => { S.hiddenCats[mode] = new Set(order.map(([i]) => i).filter((i) => i !== +o.dataset.i)); refresh(); });
  });
}
// Floating legend: drag by the title bar, resize from the bottom-right corner (CSS resize), A−/A+ for text size; position/size persist.
// Generic floating panel: drag by its header, resize with the CSS corner handle; position and size persist.
function initFloating(el, head, KEY) {
  const view = $("viewer"), saved = (() => { try { return JSON.parse(localStorage.getItem(KEY) || "{}"); } catch (e) { return {}; } })();
  if (saved.w) el.style.width = saved.w + "px"; if (saved.h) el.style.height = saved.h + "px";
  const place = (l, t) => { el.style.right = "auto"; el.style.bottom = "auto"; el.style.left = l + "px"; el.style.top = t + "px"; };   // switch from the CSS default corner to explicit left/top
  if (saved.l != null) place(saved.l, saved.t);
  const save = () => { const r = el.getBoundingClientRect(), v = view.getBoundingClientRect(); localStorage.setItem(KEY, JSON.stringify({ w: Math.round(r.width), h: Math.round(r.height), l: Math.round(r.left - v.left), t: Math.round(r.top - v.top) })); };
  let drag = null;
  head.addEventListener("pointerdown", (e) => { if (e.target.closest("button, select, input, label")) return; const r = el.getBoundingClientRect(), v = view.getBoundingClientRect();
    drag = { dx: e.clientX - r.left, dy: e.clientY - r.top, v }; place(r.left - v.left, r.top - v.top); head.setPointerCapture(e.pointerId); });
  head.addEventListener("pointermove", (e) => { if (!drag) return;
    el.style.left = Math.max(0, Math.min(drag.v.width - 60, e.clientX - drag.dx - drag.v.left)) + "px"; el.style.top = Math.max(0, Math.min(drag.v.height - 36, e.clientY - drag.dy - drag.v.top)) + "px"; });
  head.addEventListener("pointerup", () => { if (drag) { drag = null; save(); } });
  const grip = el.querySelector(".fgrip"); let rs = null;   // own resize corner (the native one is covered by the canvas)
  if (grip) {
    grip.addEventListener("pointerdown", (e) => { const r = el.getBoundingClientRect(), v = view.getBoundingClientRect(); place(r.left - v.left, r.top - v.top);
      rs = { x: e.clientX, y: e.clientY, w: r.width, h: r.height, v, l: r.left - v.left, t: r.top - v.top }; grip.setPointerCapture(e.pointerId); e.preventDefault(); e.stopPropagation(); });
    grip.addEventListener("pointermove", (e) => { if (!rs) return;
      el.style.width = Math.max(260, Math.min(rs.v.width - rs.l, rs.w + e.clientX - rs.x)) + "px"; el.style.height = Math.max(220, Math.min(rs.v.height - rs.t, rs.h + e.clientY - rs.y)) + "px"; });
    grip.addEventListener("pointerup", () => { if (rs) { rs = null; save(); } });
  }
  if (window.ResizeObserver) { let t; new ResizeObserver(() => { if (!el.hidden) { clearTimeout(t); t = setTimeout(save, 300); } }).observe(el); }
}
function initOvLegend() {
  const el = $("ovleg"), view = $("viewer"), KEY = "xb_ovleg";
  const saved = (() => { try { return JSON.parse(localStorage.getItem(KEY) || "{}"); } catch (e) { return {}; } })();
  S.ovFont = saved.font || 16; el.style.fontSize = S.ovFont + "px";
  if (saved.w) el.style.width = saved.w + "px"; if (saved.h) el.style.height = saved.h + "px";
  if (saved.l != null) { el.style.left = saved.l + "px"; el.style.top = saved.t + "px"; el.style.right = "auto"; }
  const save = () => { const r = el.getBoundingClientRect(), v = view.getBoundingClientRect();
    localStorage.setItem(KEY, JSON.stringify({ font: S.ovFont, w: Math.round(r.width), h: el.style.height ? Math.round(r.height) : 0, l: Math.round(r.left - v.left), t: Math.round(r.top - v.top) })); };
  const font = (d) => { S.ovFont = Math.max(10, Math.min(40, S.ovFont + d)); el.style.fontSize = S.ovFont + "px"; save(); };
  $("ovsm").addEventListener("click", () => font(-2)); $("ovbg").addEventListener("click", () => font(2));
  $("ovx").addEventListener("click", () => { S.ovLegend = false; $("legbtn").classList.remove("on"); updateOvLegend(); });
  let drag = null;
  $("ovhead").addEventListener("pointerdown", (e) => {
    if (e.target.closest("button")) return;
    const r = el.getBoundingClientRect(), v = view.getBoundingClientRect();
    drag = { dx: e.clientX - r.left, dy: e.clientY - r.top, v }; el.style.right = "auto"; el.style.left = r.left - v.left + "px"; el.style.top = r.top - v.top + "px";
    $("ovhead").setPointerCapture(e.pointerId);
  });
  $("ovhead").addEventListener("pointermove", (e) => {
    if (!drag) return; const r = el.getBoundingClientRect();
    el.style.left = Math.max(0, Math.min(drag.v.width - 60, e.clientX - drag.dx - drag.v.left)) + "px";
    el.style.top = Math.max(0, Math.min(drag.v.height - 36, e.clientY - drag.dy - drag.v.top)) + "px";
  });
  $("ovhead").addEventListener("pointerup", () => { if (drag) { drag = null; save(); } });
  let rs = null; const grip = $("ovgrip");
  grip.addEventListener("pointerdown", (e) => { const r = el.getBoundingClientRect(), v = view.getBoundingClientRect();
    el.style.right = "auto"; el.style.left = r.left - v.left + "px"; el.style.top = r.top - v.top + "px";
    rs = { x: e.clientX, y: e.clientY, w: r.width, h: r.height, v, l: r.left - v.left, t: r.top - v.top }; grip.setPointerCapture(e.pointerId); e.preventDefault(); e.stopPropagation(); });
  grip.addEventListener("pointermove", (e) => { if (!rs) return;
    el.style.width = Math.max(200, Math.min(rs.v.width - rs.l, rs.w + e.clientX - rs.x)) + "px"; el.style.height = Math.max(70, Math.min(rs.v.height - rs.t, rs.h + e.clientY - rs.y)) + "px"; });
  grip.addEventListener("pointerup", () => { if (rs) { rs = null; save(); } });
  if (window.ResizeObserver) new ResizeObserver(() => { if (!el.hidden) { clearTimeout(S.ovT); S.ovT = setTimeout(save, 300); } }).observe(el);
}
function countsAll(k) { // counts over the filtered population ignoring legend-hiding (so hidden items stay listed)
  const c = S.cats[k], n = catMeta(k).levels.length, out = new Array(n).fill(0); const keep = new Uint8Array(S.n).fill(1);
  for (const [fk, set] of Object.entries(S.filters)) if (set.size && S.cats[fk]) { const mu = isMulti(fk), sm = mu ? selMask(set) : 0; for (let i = 0; i < S.n; i++) if (keep[i] && !(mu ? (S.cats[fk][i] & sm) !== 0 : set.has(S.cats[fk][i]))) keep[i] = 0; }
  S.exprFilters.forEach((e, q) => { const v = S.exprV[q]; if (v && v.length === S.n) for (let i = 0; i < S.n; i++) if (keep[i] && !exprPass(e, v[i])) keep[i] = 0; });
  for (let i = 0; i < S.n; i++) if (keep[i] && c[i] !== MISSING) out[c[i]]++;
  return out;
}

function updateComp() {
  const k = S.compField, sel = [...S.selected].filter((i) => S.vis[i]), subset = sel.length ? sel : null;
  const { out, miss } = countsFor(k, subset), tot = out.reduce((a, b) => a + b, 0) + miss;
  $("comphead").textContent = `Composition of ${subset ? "selected cells" : "filtered cells"} (n=${fmt(tot)})`;
  const order = out.map((c, i) => [i, c]).filter(([, c]) => c > 0).sort((a, b) => b[1] - a[1]).slice(0, 25), mx = order[0]?.[1] || 1;
  $("comp").innerHTML = order.map(([i, c]) => `<div class="bar" data-i="${i}" title="click to filter to this category"><span class="b" style="width:${Math.max(2, 120 * c / mx)}px;background:${catColor(k, i)}"></span>
    <span class="t">${display(k, i)}</span><span class="muted small">${fmt(c)} · ${(100 * c / tot).toFixed(1)}%</span></div>`).join("") +
    (miss ? `<div class="small muted">${fmt(miss)} cells with missing annotation</div>` : "") + `<div class="small muted">% of ${subset ? "selected" : "filtered"} cells (denominator ${fmt(tot)}). Descriptive only — no statistical test.</div>`;
  $("comp").querySelectorAll(".bar").forEach((b) => b.addEventListener("click", () => { S.filters[k] = new Set([+b.dataset.i]); refresh(); }));
}

function updateSelInfo() {
  const sel = [...S.selected].filter((i) => S.vis[i]);
  $("selinfo").innerHTML = sel.length ? `${fmt(sel.length)} cells selected (of ${fmt(S.visIdx.length)} passing filters)` : "No cells selected.";
}

function updateCrumbs() {  // active region filters (levels are independent; they combine with AND)
  const parts = [];
  for (const k of ["region_lvl1", "region_lvl2", "region_lvl3"]) {
    const s = S.filters[k]; if (s?.size) parts.push(`${ROLE_LABEL[k]}: ${[...s].map((i) => display(k, i)).join(", ")}`);
  }
  for (const e of exprSpecs()) parts.push(exprDesc(e));
  $("crumbs").textContent = parts.join(" · ");
}

function updateLegNote() {
  const f = S.n - S.visIdx.length;
  $("legnote").textContent = f ? `${fmt(f)} cells are filtered out and shown in grey; ${fmt(S.visIdx.length)} highlighted.` : "";
}
function updateAll() { clearTimeout(S.gzT); S.gzT = setTimeout(refreshGlomZoom, 300); clearTimeout(S.ghT); S.ghT = setTimeout(() => { if (!$("ghet").hidden) renderGH(); }, 400); clearTimeout(S.umT); S.umT = setTimeout(() => updateUmap(), 120); renderGallery(); renderSlide(); renderGlom(); clearTimeout(S.rgT); S.rgT = setTimeout(renderRgn, 450); clearTimeout(S.gcT); S.gcT = setTimeout(renderGlomCmp, 600); clearTimeout(S.rnT); S.rnT = setTimeout(renderRnb, 500); clearTimeout(S.nbT); S.nbT = setTimeout(renderNbr, 400); clearTimeout(S.nbhT); S.nbhT = setTimeout(renderNbh, 450); render(); updateLegNote(); updateStats(); updateLegend(); updateOvLegend(); updateComp(); updateSelInfo(); updateCrumbs(); renderExprList(); [...ROLE_ORDER, "path_region"].filter((k) => S.meta.cats[k]).forEach(renderFilterOpts); }

function syncTransformUI() {   // the main panel and the Glom heterogeneity window show the same gene / signature settings
  const set = (id, v) => { const e = $(id); if (e && e.value !== String(v)) e.value = v; }, sm = S.smooth ? S.smooth + " µm" : "off";
  set("smooth", S.smooth); $("smoothlab").textContent = sm; set("zscope", S.zscope); set("layer", S.layer); set("zregion", S.within);
  set("ghsm", S.smooth); $("ghsmlab").textContent = sm; set("ghzs", S.zscope); set("ghexl", S.layer); set("ghwi", S.within);
  const t = RG.level >= 0 ? `the region “${rgName()}” (Regions window)` : "the region chosen in the Regions window (none yet)";
  $("zregionopt").textContent = t; $("ghwi").options[2].textContent = t;
}
async function refresh(fit = false) {
  loading(true); syncTransformUI();
  try {
    if (S.colorBy === "__gene" || S.colorBy === "__dual") { if (S.geneA) await ensureGene(S.geneA); if (S.colorBy === "__dual" && S.geneB) await ensureGene(S.geneB); }
    else if (S.colorBy === "__sig") await loadSignature();
    else if (S.colorBy === "__sigdual") await Promise.all([loadSignature(""), loadSignature("2")]);
    else if (S.meta.cats[S.colorBy]) await ensureCat(S.colorBy);
    else if (S.meta.numeric[S.colorBy]) await ensureNum(S.colorBy);
    if (S.compField) await ensureCat(S.compField);
    S.exprV = await Promise.all(S.exprFilters.map((e) => galExpr(S.sample, e.genes, e)));
    if (fit || !S.viewState) fitView();
    updateAll();
  } catch (e) { alert("Error: " + e.message); }
  loading(false);
}

// Signature score = mean over signature genes of the selected layer (log-normalised by default).
async function loadSignature(sfx = "") {   // sfx "" = signature A, "2" = signature B
  const genes = S["sigGenes" + sfx], key = [S.sample, S.layer, S.smooth, S.zscope, withinSpec(), ...genes].join("|");
  if (!genes.length) { S["sigVec" + sfx] = null; S["sigKey" + sfx] = ""; return; }
  if (key === S["sigKey" + sfx] && S["sigVec" + sfx]) return;
  let vec;
  if (transformed()) vec = await bget(exprUrl(genes), Float32Array);  // server: mean -> smooth -> z-score
  else {
    const vs = await Promise.all(genes.map(ensureGene));
    vec = new Float32Array(S.n);
    for (const v of vs) for (let i = 0; i < S.n; i++) vec[i] += v[i];
    for (let i = 0; i < S.n; i++) vec[i] /= vs.length;
  }
  S["sigVec" + sfx] = vec; S["sigKey" + sfx] = key;
}

// ---------------------------------------------------------------- gallery: all tissues of a disease subtype
// Donors: a donor can have two tissue pieces (sample ids ending _fragment-1 / _fragment-2). Windows show the blinded donor id and count donors;
// a piece of a two-piece donor is labelled "<donor> · fragment N" so the two tiles can be told apart.
const donorOfId = (id) => S.meta.samples.find((s) => s.id === id)?.donor || String(id).replace(/_fragment[-_]?\d+$/, "");
const dlabel = (id) => { const d = donorOfId(id), f = /_fragment[-_]?(\d+)$/.exec(String(id)); return f && S.meta.samples.filter((s) => s.donor === d).length > 1 ? `${d} · fragment ${f[1]}` : d; };
const nDonors = (list) => new Set(list.map((s) => (typeof s === "string" ? donorOfId(s) : s.donor))).size;
const donTxt = (n) => `${n} donor${n === 1 ? "" : "s"}`;
const subtypeOf = (s) => (s.disease_subtype && s.disease_subtype !== "nan" ? s.disease_subtype : s.disease_type) || "unknown";
const GAL_BOX = 300;
const SUB_ORDER = ["Control", "No LN", "Class I/II", "Class III", "Class V"];  // display order; any other subtype follows
const sortedSubtypes = () => [...new Set(S.meta.samples.map(subtypeOf))].sort((a, b) => { const r = (t) => (SUB_ORDER.includes(t) ? SUB_ORDER.indexOf(t) : 99); return r(a) - r(b) || a.localeCompare(b); });
const galData = {};  // sample -> { xy: Promise, cats: {key: Promise} } (per-sample arrays are cached for the session)
function galGet(sample, key) {
  const c = galData[sample] || (galData[sample] = { cats: {} });
  const q = encodeURIComponent(sample);
  if (key === "xy") return c.xy || (c.xy = bget(`/api/sample/${q}/xy`, Float32Array));
  return c.cats[key] || (c.cats[key] = bget(`/api/sample/${q}/cat/${encodeURIComponent(key)}`, Uint16Array));
}
function galNum(sample, key) {
  const c = galData[sample] || (galData[sample] = { cats: {} }); c.nums = c.nums || {};
  return c.nums[key] || (c.nums[key] = bget(`/api/sample/${encodeURIComponent(sample)}/num/${encodeURIComponent(key)}`, Float32Array));
}
function galExpr(sample, genes, o = { layer: S.layer, sigma: S.smooth, scope: S.zscope, within: withinSpec() }) {  // mean over genes (server-side); default = smoothed / z-scored exactly as in the main view, filters pass their own settings
  const c = galData[sample] || (galData[sample] = { cats: {} }); c.ex = c.ex || {};
  const key = [...genes, o.layer, o.sigma, o.scope, o.within].join("|");
  return c.ex[key] || (c.ex[key] = bget(`/api/sample/${encodeURIComponent(sample)}/expr?genes=${encodeURIComponent(genes.join(","))}&layer=${o.layer}&sigma=${o.sigma}&scope=${o.scope}&within=${encodeURIComponent(o.within)}`, Float32Array));
}
// One sample's arrays for the current Colour-by + filters (see makeScale): any colour mode works in the windows.
async function winBundle(sid) {
  const mode = S.colorBy, fkeys = Object.keys(S.filters).filter((k) => S.filters[k].size && S.meta.cats[k]);
  const isCat = !!S.meta.cats[mode], hid = isCat ? (S.hiddenCats[mode] || new Set()) : new Set();
  const A = mode === "__sig" || mode === "__sigdual" ? S.sigGenes : (mode === "__gene" || mode === "__dual") && S.geneA ? [S.geneA] : [];
  const Bg = mode === "__dual" && S.geneB ? [S.geneB] : mode === "__sigdual" ? S.sigGenes2 : [];
  const [xy, fc, codes, num, a, bv, ev] = await Promise.all([galGet(sid, "xy"), Promise.all(fkeys.map((k) => galGet(sid, k))), isCat ? galGet(sid, mode) : null,
    !isCat && S.meta.numeric[mode] ? galNum(sid, mode) : null, A.length ? galExpr(sid, A) : null, Bg.length ? galExpr(sid, Bg) : null, Promise.all(S.exprFilters.map((e) => galExpr(sid, e.genes, e)))]);
  const n = xy.length / 2, pass = new Uint8Array(n).fill(1);
  fkeys.forEach((k, q) => { const set = S.filters[k], mu = isMulti(k), sm = mu ? selMask(set) : 0; for (let i = 0; i < n; i++) if (pass[i] && !(mu ? (fc[q][i] & sm) !== 0 : set.has(fc[q][i]))) pass[i] = 0; });
  S.exprFilters.forEach((e, q) => { for (let i = 0; i < n; i++) if (pass[i] && !exprPass(e, ev[q][i])) pass[i] = 0; });
  if (hid.size) for (let i = 0; i < n; i++) if (pass[i] && hid.has(codes[i])) pass[i] = 0;
  return { n, xy, pass, codes, num, A: a, B: bv, sample: sid };
}
// Scale shared by every tile of a window: categorical colours need no pooling; continuous ones are pooled over all shown samples.
async function winScale(bundlePromises, tok, cur) {
  if (S.meta.cats[S.colorBy]) return makeScale([]);
  const all = await Promise.all(bundlePromises);
  return tok() === cur ? makeScale(all) : null;
}
const colorName = () => $("colorby").selectedOptions[0]?.text || S.colorBy;
const winNote = () => { const fk = Object.keys(S.filters).filter((k) => S.filters[k].size).length + S.exprFilters.length; return fk || (S.hiddenCats[S.colorBy] || new Set()).size ? " · filters applied (others grey)" : ""; };
const rgbStr = (c) => `rgb(${c[0]},${c[1]},${c[2]})`;

function closeAllWindows() {
  for (const id of ["slide", "gallery", "glom", "glomcmp", "glomzoom", "nbr", "nbrtip", "nbh", "cmp", "sigwin", "vio", "ghet", "rgn", "rnb"]) if ($(id)) $(id).hidden = true;
}
function openGallery() {
  closeAllWindows();
  const subs = sortedSubtypes();
  $("galsub").innerHTML = `<option value="__all">All donors (${nDonors(S.meta.samples)})</option>` +
    subs.map((t) => `<option value="${t}">${t} (${nDonors(S.meta.samples.filter((s) => subtypeOf(s) === t))})</option>`).join("");
  $("galsub").value = S.galSub && (S.galSub === "__all" || subs.includes(S.galSub)) ? S.galSub : "__all";
  $("gallery").hidden = false; renderGallery();
}
function closeGallery() { $("gallery").hidden = true; if ($("cmp")) $("cmp").hidden = true; if ($("sigwin")) $("sigwin").hidden = true; if ($("vio")) $("vio").hidden = true; }   // opening any other window also closes the comparison window
let galToken = 0;
async function renderGallery() {
  if ($("gallery").hidden) return;
  const tok = ++galToken, sub = $("galsub").value; S.galSub = sub;
  const sel = S.galSel || (S.galSel = new Set()), only = $("galonly").checked, box0 = +$("galsz").value;
  const groups = sortedSubtypes().filter((t) => sub === "__all" || t === sub).map((t) => ({ t, samples: S.meta.samples.filter((s) => subtypeOf(s) === t && (!only || sel.has(s.id))) })).filter((g) => g.samples.length);
  const samples = groups.flatMap((g) => g.samples);
  const nT = samples.length, box = box0 > 0 ? box0 : GAL_BOX;
  $("galsz").disabled = only; $("galgrid").classList.toggle("cmpgrid", only); if (!only) $("galgrid").removeAttribute("style");
  $("galn").textContent = sel.size;
  $("galinfo").textContent = `${donTxt(nDonors(samples))} · colour: ${colorName()}${winNote()} · ${only ? "" : "scroll for more · "}${only ? "window split per tissue · wheel = zoom (cell outlines appear when zoomed in), drag = pan, double-click = reset, click a name to open the tissue" : "click a tissue to open it · tick tissues to compare them"}`;
  $("galzoomctl").hidden = !only;
  if (only && !samples.length) { $("galgrid").innerHTML = `<div class="muted">No ticked tissue${sub === "__all" ? "" : " in this group"}. Untick “only ticked” and tick the tissues to compare.</div>`; $("galleg").innerHTML = ""; return; }
  $("galgrid").innerHTML = groups.map((g) => (only ? "" : `<div class="galsec"><div class="sechead">${g.t} · ${donTxt(nDonors(g.samples))}</div><div class="secgrid">`) +
    g.samples.map((s) => `<div class="tile ${s.id === S.sample ? "cur" : ""}" data-id="${s.id}" title="Open ${dlabel(s.id)}"><input type="checkbox" class="gsel" title="tick to compare this tissue" ${sel.has(s.id) ? "checked" : ""}><span class="zb"></span><canvas width="${box}" height="60"></canvas>
    <div class="cap" style="max-width:${box}px"><b>${dlabel(s.id)}</b><br><span class="muted">${fmt(s.n)} cells · ${s.disease_type} · ${subtypeOf(s)}</span></div></div>`).join("") + (only ? "" : `</div></div>`)).join("");
  $("galgrid").querySelectorAll(".tile").forEach((t) => {
    t.classList.toggle("sel", sel.has(t.dataset.id));
    t.addEventListener("click", async (e) => { if (only && e.target.closest("canvas")) return; $("sample").value = t.dataset.id; closeGallery(); await loadSample(t.dataset.id); });   // compare view: the picture is for zooming; click the caption to open the tissue
    const cb = t.querySelector(".gsel"); cb.addEventListener("click", (e) => e.stopPropagation());
    cb.addEventListener("change", () => { cb.checked ? sel.add(t.dataset.id) : sel.delete(t.dataset.id); t.classList.toggle("sel", cb.checked); $("galn").textContent = sel.size; if (only) renderGallery(); });
  });
  const bp = samples.map((s) => winBundle(s.id)), heo = heOpt("gal"), wantPoly = only || HEW.gal.outl;
  const sc = await winScale(bp, () => galToken, tok).catch((e) => { console.error(e); return null; });
  if (!sc || tok !== galToken) return;
  $("galleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  S.ctxs.gallery = { bp, sc, ids: samples.map((s) => s.id) };
  let dims = null;
  if (only) {   // every chunk gets an equal share of the window below the header; the tissue is centred in it and fits it
    dims = cmpLayout("gallery", "galgrid", nT);
  }
  samples.forEach((s, k) => Promise.all([bp[k], wantPoly ? galPoly(s.id).catch(() => null) : null]).then(([b, poly]) => {
    const cv = $("galgrid").querySelector(`.tile[data-id="${CSS.escape(s.id)}"] canvas`);
    if (tok !== galToken || !cv) return;
    const t = cv.closest(".tile"), cache = {}, he = heo.on && HE.data?.[s.id] ? { a: heo.a, sample: s.id, redo: () => heRedo(t) } : null;
    const cv2 = only ? hePane(t, cv, "gal", s.id) : null, dm = cv2 ? hePairDims(dims, "gal") : dims, he2 = cv2 ? { a: 1, sample: s.id, redo: () => heRedo(t) } : null;
    const rec = t._dt = { he, last: null, last2: null, bb: { x0: s.x[0], y0: s.y[0], w: s.x[1] - s.x[0] || 1, h: s.y[1] - s.y[0] || 1 }, draw: (vw, vh) => { rec.last = vw; rec.last2 = vh; drawTile(cv, s, b, sc, box, vw, cache, dm, poly, he, "gal"); if (cv2) drawTile(cv2, s, b, sc, box, HEW.gal.link ? vw : vh, null, dm, null, he2, "gal"); } }, v = only ? galCfg.views()[galCfg.sync() ? "__all" : s.id] || null : null;
    rec.draw(v, only ? galCfg.views()[(galCfg.sync() ? "__all" : s.id) + "#h"] || null : null); if (only) { attachZoom(t, rec, galCfg); if (v && v.z > 1.01) t.querySelector(".zb").textContent = "×" + v.z.toFixed(1); }
  }).catch((e) => console.error(e)));
}
// ---------------------------------------------------------------- slide: every tissue of one slide (run) at its true position, each in a labelled box
// Tissues of a slide share one coordinate frame (the run's), so they are grouped by run (= disease subtype / run name in imported projects).
function tissueName(s, grp) {   // "<run>_tissue_2" -> "tissue_2": the slide name is dropped, the window is already one slide
  let n = String(s.id); if (n.startsWith(grp + "_")) return n.slice(grp.length + 1);
  const m = /(tissue_\d+|small_piece_\d+)$/.exec(n); return m ? m[1] : n;
}
let slideTok = 0, slideObs = null;
const slideGroups = () => [...new Set(S.meta.samples.map(subtypeOf))].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }));
function openSlide() {
  closeAllWindows();
  const groups = slideGroups();
  $("slidesel").innerHTML = `<option value="__all">All slides (${groups.length})</option>` + groups.map((g) => `<option value="${g}">${g} (${S.meta.samples.filter((s) => subtypeOf(s) === g).length} tissues)</option>`).join("");
  const cur = S.meta.samples.find((s) => s.id === S.sample); $("slidesel").value = S.slideSel && (S.slideSel === "__all" || groups.includes(S.slideSel)) ? S.slideSel : cur ? subtypeOf(cur) : groups[0];
  $("slide").hidden = false; renderSlide();
}
const closeSlide = () => { $("slide").hidden = true; };
function slideGeom(samples) {   // extent of a slide in µm (plus a margin)
  const x0 = Math.min(...samples.map((s) => s.x[0])), x1 = Math.max(...samples.map((s) => s.x[1])), y0 = Math.min(...samples.map((s) => s.y[0])), y1 = Math.max(...samples.map((s) => s.y[1]));
  const pad = 0.03 * Math.max(x1 - x0, y1 - y0, 1); return { x0, y0, pad, w: x1 - x0 + 2 * pad, h: y1 - y0 + 2 * pad };
}
// Draw one slide (all its tissues at their true positions, each in a labelled box) into cv at k CSS px per µm; dpr = backing-store scale. Returns the boxes for hit-testing.
const SLIDE_TOP = 24, SLIDE_BOT = 24;   // strips above / below the tissues: labels (and the scale bar) never sit on a tissue
function drawSlideInto(cv, g, samples, bs, sc, k, geom, dpr, ts = 1, opt = {}) {   // ts = text scale (exports draw at a high dpr, so their text is scaled down); opt = { imgs: draw the images ticked in Images, imgOnly: no cells (the pane beside), redo }
  const Wd = Math.round(geom.w * k), Hd = Math.round(geom.h * k) + SLIDE_TOP + SLIDE_BOT;
  cv.width = Math.round(Wd * dpr); cv.height = Math.round(Hd * dpr); cv.style.width = Wd + "px"; cv.style.height = Hd + "px";
  const ctx = cv.getContext("2d"); ctx.scale(dpr, dpr); ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, Wd, Hd);
  const X = (x) => (x - geom.x0 + geom.pad) * k, Y = (y) => (y - geom.y0 + geom.pad) * k + SLIDE_TOP, px = Math.max(1.2, 8 * k);
  if (opt.imgs) for (const s of samples) drawHE(ctx, { sample: s.id, redo: opt.redo, w: "slide", side: !!opt.imgOnly }, X, Y, Wd, Hd, k);   // registered H&E / IF of every tissue, at its true place
  const cellsOn = !opt.imgOnly && hw("slide").cells; ctx.globalAlpha = opt.imgs && hw("slide").on ? hw("slide").a : 1;
  if (cellsOn) for (const dim of [true, false]) bs.forEach((b) => {   // filtered-out cells first, coloured ones on top
    const groups = new Map();
    for (let i = 0; i < b.n; i++) {
      const col = b.pass[i] ? (S.grey ? rgbStr(GREY) : rgbStr(sc.color(b, i))) : "rgba(105,110,122,0.45)";
      if (col.startsWith("rgba") !== dim) continue; (groups.get(col) || groups.set(col, []).get(col)).push(i);
    }
    for (const [col, list] of groups) { ctx.fillStyle = col; for (const i of list) ctx.fillRect(X(b.xy[i * 2]) - px / 2, Y(b.xy[i * 2 + 1]) - px / 2, px, px); }
  });
  ctx.globalAlpha = 1;
  const boxes = samples.map((s) => ({ id: s.id, s, x: X(s.x[0]) - 5, y: Y(s.y[0]) - 5, w: (s.x[1] - s.x[0]) * k + 10, h: (s.y[1] - s.y[0]) * k + 10 })), hits = boxes.map((q) => ({ id: q.id, x: q.x, y: q.y, w: q.w, h: q.h }));
  const hit = (r, but) => boxes.some((q) => q !== but && r.x < q.x + q.w && r.x + r.w > q.x && r.y < q.y + q.h && r.y + r.h > q.y);   // does rectangle r cover another tissue's box?
  ctx.textBaseline = "alphabetic";
  for (const q of boxes) {
    const cur = q.id === S.sample; ctx.lineWidth = cur ? 2.5 : 1.5; ctx.strokeStyle = cur ? "#4c8dff" : "rgba(215,221,230,0.75)"; ctx.strokeRect(q.x, q.y, q.w, q.h);
    const nm = tissueName(q.s, g), cnt = ` · ${fmt(q.s.n)} cells`, fs = 13 * ts; ctx.font = `600 ${fs}px ${FONT}`; const wn = ctx.measureText(nm).width; ctx.font = `${12 * ts}px ${FONT}`;
    const lw = wn + ctx.measureText(cnt).width, lx = Math.max(2, Math.min(q.x, Wd - lw - 4));   // slid left when it would run off the picture
    let ly = q.y - 5; const above = { x: lx, y: ly - fs, w: lw, h: fs + 3 };   // label above the box; below it when that would cover another tissue
    if (hit(above, q) || ly - fs < 0) { const below = { x: lx, y: q.y + q.h + 2, w: lw, h: fs + 3 }; if (!hit(below, q) && below.y + below.h <= Hd) ly = q.y + q.h + fs + 1; }
    ctx.font = `600 ${fs}px ${FONT}`; ctx.fillStyle = cur ? "#8fb4ff" : "#e6ebf2"; ctx.fillText(nm, lx, ly); ctx.font = `${12 * ts}px ${FONT}`; ctx.fillStyle = "#8b95a5"; ctx.fillText(cnt, lx + wn, ly);
  }
  const target = 0.15 * Wd / k, p = 10 ** Math.floor(Math.log10(target)), nice = [1, 2, 5, 10].map((m) => m * p).reduce((a, c) => (Math.abs(c - target) < Math.abs(a - target) ? c : a));
  ctx.fillStyle = "#fff"; ctx.fillRect(10, Hd - 14 * ts, nice * k, 3 * ts); ctx.font = `${12 * ts}px ${FONT}`; ctx.fillText(nice >= 1000 ? nice / 1000 + " mm" : nice + " µm", 10, Hd - 20 * ts);
  return hits;
}
async function renderSlide() {
  if ($("slide").hidden) return;
  const tok = ++slideTok, g = $("slidesel").value, all = g === "__all"; S.slideSel = g;
  $("slidecv").hidden = all; $("slideall").hidden = !all; $("slidecolsl").hidden = !all; $("slidecv2").hidden = true;
  if (HEW.slide.on || HEW.slide.side) await heReady();
  if (all) return renderSlideAll(tok);
  const samples = S.meta.samples.filter((s) => subtypeOf(s) === g); if (!samples.length) return;
  const bp = samples.map((s) => winBundle(s.id)), sc = await winScale(bp, () => slideTok, tok).catch((e) => { console.error(e); return null; });
  if (!sc || tok !== slideTok) return;
  const bs = await Promise.all(bp); if (tok !== slideTok) return;
  const total = samples.reduce((a, s) => a + s.n, 0);
  $("slideleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  $("slideinfo").textContent = `${samples.length} tissue${samples.length === 1 ? "" : "s"} · ${fmt(total)} cells · colour: ${colorName()}${winNote()} · click a tissue to open it`;
  const geom = slideGeom(samples), win = $("slide"), cv = $("slidecv"), cv2 = $("slidecv2"), side = HEW.slide.side;
  const Wc = Math.max(300, win.clientWidth - 30 - (side ? 8 : 0)) / (side ? 2 : 1), availH = Math.max(300, win.clientHeight - (cv.getBoundingClientRect().top - win.getBoundingClientRect().top + win.scrollTop) - 16);
  const k = Math.min(Wc / geom.w, (availH - SLIDE_TOP - SLIDE_BOT) / geom.h), dpr = window.devicePixelRatio || 1, he = { imgs: HEW.slide.on, redo: slideHeRedo };
  cv._hits = drawSlideInto(cv, g, samples, bs, sc, k, geom, dpr, 1, he);
  cv2.hidden = !side;
  if (side) cv2._hits = drawSlideInto(cv2, g, samples, bs, sc, k, geom, dpr, 1, { imgs: true, imgOnly: true, redo: slideHeRedo });   // the same slide again: the ticked images, no cells; a click opens a tissue in either view
}
let slideHeT = 0; const slideHeRedo = () => { clearTimeout(slideHeT); slideHeT = setTimeout(renderSlide, 160); };   // an image tile arrived
// All slides: one scrolling page, common scale (same µm per pixel everywhere) and one colour scale over every tissue. Each slide is drawn when it scrolls into view.
async function renderSlideAll(tok) {
  const groups = slideGroups(), bySlide = groups.map((g) => ({ g, samples: S.meta.samples.filter((s) => subtypeOf(s) === g) })), cols = +$("slidecols").value, win = $("slide"), gap = 14;
  const colW = Math.max(250, Math.floor((win.clientWidth - 30 - gap * (cols - 1)) / cols)), geoms = bySlide.map((o) => slideGeom(o.samples)), k = colW / Math.max(...geoms.map((q) => q.w));
  const cache = new Map(), getB = (id) => cache.get(id) || (cache.set(id, winBundle(id)), cache.get(id)), N = S.meta.samples.length;
  let sc;
  if (S.meta.cats[S.colorBy]) sc = makeScale([]);
  else {   // a gene / signature / numeric colour: the shared scale needs every tissue's values first
    let done = 0; const bsAll = await Promise.all(S.meta.samples.map((s) => getB(s.id).then((b) => { if (tok === slideTok) $("slideinfo").textContent = `loading values for the shared colour scale: ${++done} of ${N} tissues…`; return b; })));
    if (tok !== slideTok) return; sc = makeScale(bsAll);
  }
  $("slideleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  $("slideinfo").textContent = `${groups.length} slides · ${N} tissues · ${fmt(S.meta.samples.reduce((a, s) => a + s.n, 0))} cells · colour: ${colorName()}${winNote()} · same scale and colours on every slide · click a tissue to open it` + (HEW.slide.side ? " · choose one slide to see the images beside the cells" : "");
  const host = $("slideall"); host.style.gridTemplateColumns = `repeat(${cols}, ${colW}px)`;
  host.innerHTML = bySlide.map((o, i) => `<div class="slsec"><div class="slhead">${o.g} · ${o.samples.length} tissue${o.samples.length === 1 ? "" : "s"} · ${fmt(o.samples.reduce((a, s) => a + s.n, 0))} cells</div><canvas data-i="${i}" style="width:${Math.round(geoms[i].w * k)}px;height:${Math.round(geoms[i].h * k) + SLIDE_TOP + SLIDE_BOT}px"></canvas></div>`).join("");
  S.slideAll = { bySlide, geoms, k, sc, getB, tok };
  if (slideObs) slideObs.disconnect();
  const drawOne = async (cv) => { const i = +cv.dataset.i, o = bySlide[i], bs = await Promise.all(o.samples.map((s) => getB(s.id))); if (tok !== slideTok) return; cv._hits = drawSlideInto(cv, o.g, o.samples, bs, sc, k, geoms[i], window.devicePixelRatio || 1, 1, { imgs: HEW.slide.on, redo: slideHeRedo }); };
  slideObs = new IntersectionObserver((ents) => { for (const e of ents) if (e.isIntersecting && !e.target._drawn) { e.target._drawn = true; drawOne(e.target).catch((er) => console.error(er)); } }, { root: win, rootMargin: "600px" });
  host.querySelectorAll("canvas").forEach((cv) => slideObs.observe(cv));
}
// Save PNG in the all-slides view: a zip with one picture per slide (title, tissues in boxes, legend), same scale and colours throughout.
async function saveSlideZip() {
  const A = S.slideAll; if (!A || A.tok !== slideTok) return;
  const files = [], dprE = 2800 / (Math.max(...A.geoms.map((q) => q.w)) * A.k), scratch = document.createElement("canvas").getContext("2d");
  for (let i = 0; i < A.bySlide.length; i++) {
    const { g, samples } = A.bySlide[i]; $("slideinfo").textContent = `writing the picture of slide ${i + 1} of ${A.bySlide.length}…`;
    const bs = await Promise.all(samples.map((s) => A.getB(s.id))), cv = document.createElement("canvas"); drawSlideInto(cv, g, samples, bs, A.sc, A.k, A.geoms[i], dprE, Math.min(1, 2.2 / dprE));
    let info = A.sc.info, items = null;
    if (S.grey) info = { type: "none" };
    else if (info.type === "cat") {
      const cnt = new Map(); let dimmed = 0;
      for (const b of bs) for (let q = 0; q < b.n; q++) { if (!b.pass[q]) { dimmed++; continue; } cnt.set(b.codes[q], (cnt.get(b.codes[q]) || 0) + 1); }
      items = [...cnt.entries()].sort((a, b) => b[1] - a[1]).map(([c, n]) => [c === MISSING ? "missing" : display(S.colorBy, c), c === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, c), n]);
      if (dimmed) items.push(["filtered out", "rgba(105,110,122,0.6)", dimmed]);
    }
    const W = cv.width, u = Math.max(1, W / 1100), X = 24 * u, lh = pngLegend(scratch, info, items, X, 0, W - 2 * X, u, 6), head = 60 * u, out = document.createElement("canvas");
    out.width = W; out.height = Math.round(head + cv.height + lh + 24 * u); const x = out.getContext("2d");
    x.fillStyle = "#0f1218"; x.fillRect(0, 0, out.width, out.height); x.drawImage(cv, 0, head);
    x.fillStyle = "#d7dde6"; x.font = `600 ${22 * u}px ${FONT}`; x.fillText(g, X, 38 * u); const tw = x.measureText(g).width;
    x.font = `${14 * u}px ${FONT}`; x.fillStyle = "#8b95a5"; x.fillText(`${samples.length} tissue${samples.length === 1 ? "" : "s"} · colour: ${colorName()}${S.grey ? " · cells grey" : ""}`, X + tw + 16 * u, 38 * u);
    pngLegend(x, info, items, X, head + cv.height + 24 * u, W - 2 * X, u, 6);
    const blob = await new Promise((r) => out.toBlob(r, "image/png"));
    files.push({ name: `slide_${String(i + 1).padStart(2, "0")}_${g.replace(/[^\w-]+/g, "_").slice(0, 60)}.png`, data: new Uint8Array(await blob.arrayBuffer()) });
  }
  const a = document.createElement("a"); a.href = URL.createObjectURL(zipBlob(files)); a.download = "slides.zip"; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 4000);
  $("slideinfo").textContent = `saved slides.zip (${files.length} pictures)`;
}
function slideBind() {
  const win = $("slide"), hit = (e) => { const cv = e.target; if (cv.tagName !== "CANVAS" || !cv._hits) return null; const r = cv.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top; return cv._hits.find((h) => x >= h.x && x <= h.x + h.w && y >= h.y && y <= h.y + h.h); };
  win.addEventListener("mousemove", (e) => { if (e.target.tagName === "CANVAS") e.target.style.cursor = hit(e) ? "pointer" : "default"; });
  win.addEventListener("click", async (e) => { const h = hit(e); if (!h) return; $("sample").value = h.id; closeSlide(); await loadSample(h.id); });
  $("slidesel").addEventListener("change", renderSlide); $("slidecols").addEventListener("change", renderSlide); $("slideclose").addEventListener("click", closeSlide); $("slidebtn").addEventListener("click", () => { closeGlom(); openSlide(); });
  let rt = 0; window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(renderSlide, 200); });
}
// ---- registered H&E underlay (tile pyramids from tools/he_pyramid.py, served by /api/he/<image>/<level>/<x>_<y>.jpg; placement in the same um frame as the cells)
const HE = { data: null, req: null, imgs: new Map(), clock: 0, crops: {}, cropImgs: new Map() };
function heReady() { return HE.req || (HE.req = Promise.all([jget("/api/he").catch(() => ({})), jget("/api/heglom").catch(() => ({})), jget("/api/if").catch(() => ({}))]).then(([d, g, f]) => {
  HE.crops = g; HE.cropMap = {}; IFS.he = d; IFS.idx = f; IFS.has = Object.keys(f).length > 0;
  ifApply(); return HE.data;
})); }
// ---- IF (Xenium morphology) image: a second source for the same underlay / beside views. HE.data holds the list of whichever source is selected, so every window
// that draws the H&E draws the IF image the same way (tiles come from /api/if, rendered on request in the chosen colour and brightness).
const IFS = { idx: {}, he: {}, cur: {}, has: false, sel: new Set([0]), cfg: [], g: 1, ver: 0 };   // sel = the channels shown (tick several to overlay them); cfg[i] = { color, br } of channel i
const IF_COLORS = ["blue", "green", "red", "yellow", "cyan", "magenta", "orange", "white"];
const ifCfg = (i) => (IFS.cfg[i] ||= { color: IF_COLORS[i % IF_COLORS.length], br: 1 });
function ifApply() {
  IFS.cur = {}; const o = {};
  for (const [smp, l] of Object.entries(IFS.he)) o[smp] = [...l];
  for (const [smp, l] of Object.entries(IFS.idx)) l.forEach((e, i) => { IFS.cur[e.key] = e; if (IFS.sel.has(i)) (o[smp] = o[smp] || []).push(e); });
  HE.data = o;   // per sample: its registered H&E images, then the IF channels ticked in the left panel (drawn on top of each other, added as light); each window draws what its Images menu ticks
}
function ifTile(key, z, x, y, redo) {
  const e = IFS.cur[key]; if (!e) return null; const [, img, ch] = key.split(":"), cf = ifCfg(+ch), lo = e.lo, hi = Math.max(lo + 1, Math.round(e.hi / cf.br));
  return heLoad(HE.imgs, `${key}/${z}/${x}_${y}/${IFS.ver}`, `/api/if/${encodeURIComponent(img)}/${ch}/${z}/${x}_${y}.jpg?c=${cf.color}&lo=${lo}&hi=${hi}&g=${IFS.g}`, redo, 400);
}
function ifRerender() {   // every window that shows the image redraws; tiles are requested again with the new settings
  IFS.ver++;
  for (const k of Object.keys(HEW)) if (HEW[k].on || HEW[k].side) { HE_T[k] && clearTimeout(HE_T[k]); HE_T[k] = setTimeout(() => HE_RENDER[k] && HE_RENDER[k](), 80); }
}
function ifUi() {   // the IF channel / colour / brightness controls (shared by every window) and the Images menus' availability
  $("ifbox").hidden = !IFS.has;
  const first = Object.values(IFS.idx)[0] || [], box = $("iflist"); box.replaceChildren();
  first.forEach((c, i) => {   // one row per channel: show it, its colour, its brightness
    const row = document.createElement("div"), lab = document.createElement("label"), cb = document.createElement("input"), nm = document.createElement("span"), sel = document.createElement("select"), br = document.createElement("input");
    row.className = "ifch"; cb.type = "checkbox"; cb.dataset.ifsel = i; cb.checked = IFS.sel.has(i); nm.textContent = c.name; lab.title = "show this channel; tick several to overlay them (their light is added)"; lab.append(cb, " ", nm);
    sel.dataset.ifcol = i; sel.title = "colour of this channel"; for (const k of IF_COLORS) { const o = document.createElement("option"); o.value = k; o.textContent = k; sel.append(o); } sel.value = ifCfg(i).color;
    br.type = "range"; br.dataset.ifbr = i; br.min = -100; br.max = 100; br.value = Math.round(50 * Math.log2(ifCfg(i).br)); br.title = "brightness of this channel";
    row.append(lab, sel, br); box.append(row);
  });
  for (const k of Object.keys(HEW)) heSync(k);
}
async function ifInit() {
  await heReady(); ifUi();
  $("iflist").addEventListener("change", (e) => {
    const t = e.target;
    if (t.dataset.ifsel !== undefined) { t.checked ? IFS.sel.add(+t.dataset.ifsel) : IFS.sel.delete(+t.dataset.ifsel); ifApply(); ifRerender(); }
    else if (t.dataset.ifcol !== undefined) { ifCfg(+t.dataset.ifcol).color = t.value; ifRerender(); }
  });
  $("iflist").addEventListener("input", (e) => { const t = e.target; if (t.dataset.ifbr !== undefined) { ifCfg(+t.dataset.ifbr).br = 2 ** (t.value / 50); ifRerender(); } });
  $("ifgamma").addEventListener("input", () => { IFS.g = +$("ifgamma").value; $("ifgammalab").textContent = IFS.g.toFixed(2); ifRerender(); });
}
// Display state, one per window that draws cells: H&E underlay, opacity of the cells over it, and cell outlines instead of filled cells.
// main = the main viewer; gal / rgn = Gallery / Regions; glom = Glomeruli; gcmp = Compare glomeruli; nbh = Neighborhoods; gh = Glom heterogeneity detail; gz = glom zoom popup; sig = Signature maps (no outlines there).
const HEW = {}; ["main", "gal", "rgn", "glom", "gcmp", "nbh", "gh", "gz", "sig", "rnb", "slide"].forEach((k) => (HEW[k] = { on: false, imgs: new Set(), sides: new Set(), a: 0.55, outl: false, side: false, link: true, cells: true }));   // cells = draw the Xenium cells (off: H&E and outlines only));
const hw = (w) => HEW[w || "main"];
const heOpt = (w) => ({ on: hw(w).on, a: hw(w).on ? hw(w).a : 1 });
const heSel = (sample, w, side) => { const all = HE.data && HE.data[sample], st = HEW[w || "main"]; return all && st ? all.filter((e) => (side ? st.sides : st.imgs).has(e.key.startsWith("if:") ? "if" : "he")) : []; };   // side: the images of the pane beside the cells (else those under them)   // the images window w shows for a sample: H&E first, then IF
const heFor = (sample, w) => (hw(w).on && HE.data && HE.data[sample] ? sample : null);
const HE_SIDE = ["gal", "rgn", "glom", "gcmp", "gh", "rnb"];   // windows that can show the H&E beside the cells ("main" has its own control in the sidebar)
const HE_RENDER = {   // redraw only the window whose setting changed
  main: () => { heSideLayout(); render(); }, slide: () => renderSlide(), gal: () => renderGallery(), rgn: () => renderRgn(), glom: () => renderGlom(), gcmp: () => renderGlomCmp(), nbh: () => renderNbh(),
  gh: () => ghDraw(), rnb: () => renderRnb(), gz: () => { if (S.gz) drawGlomZoom(); },
  sig: () => { if ($("sigwin") && !$("sigwin").hidden) renderSigPane(); },
};
// Report pages (QC report, Profile, validation JSON) open in a new tab; the native app window cannot do that, so it asks the app to open them in the default browser.
function openTab(path, ev) {
  const api = window.pywebview && window.pywebview.api;
  if (api && api.open_url) { if (ev) ev.preventDefault(); api.open_url(path); return false; }
  if (ev) return true;                      // a plain link: the browser follows it
  window.open(path, "_blank"); return false;
}
const IMG_NAMES = { he: "H&E", if: "IF", scan: "Slide scan (not aligned)" };
const heImgsHTML = (w) => `<button type="button" class="heimgsbtn" data-w="${w}" title="which images this window shows: under the cells (overlay) and / or in a pane beside them">Images ▾</button>`;
function heAvail() { return { he: Object.keys(IFS.he).length > 0, if: IFS.has, scan: !!(S.meta && S.meta.has_he_slide) }; }
function heMenu(btn) {   // one menu for every window: ticks = that window's images; fixed position so no window's overflow can clip it
  let m = $("heimgsmenu");
  if (!m) {
    m = document.createElement("div"); m.id = "heimgsmenu"; m.hidden = true;
    m.innerHTML = `<div class="imgrow imghead"><span></span><span>Under</span><span>Beside</span></div>` + Object.entries(IMG_NAMES).map(([k, n]) => `<div class="imgrow"><span>${n} <span class="muted small" data-none="${k}"></span></span><input type="checkbox" data-img="${k}" data-mode="under" title="draw it under the cells"><input type="checkbox" data-img="${k}" data-mode="side" title="show it in the pane beside the cells"></div>`).join("");
    document.body.appendChild(m);
    m.addEventListener("change", (e) => { const i = e.target.closest("input[data-img]"); if (i && m._w) heSet(m._w, "img", { src: i.dataset.img, on: i.checked, side: i.dataset.mode === "side" }); });
    document.addEventListener("click", (e) => { if (!m.hidden && !e.target.closest("#heimgsmenu") && !e.target.closest(".heimgsbtn")) m.hidden = true; });
  }
  if (!m.hidden && m._w === btn.dataset.w) { m.hidden = true; return; }
  m._w = btn.dataset.w; const av = heAvail(), st = hw(m._w);
  m.querySelectorAll("input[data-img]").forEach((i) => { const sc = i.dataset.img === "scan"; i.checked = (i.dataset.mode === "side" ? st.sides : st.imgs).has(i.dataset.img); i.disabled = !av[i.dataset.img] || (sc && (m._w !== "main" || i.dataset.mode === "under")); i.style.visibility = sc && i.dataset.mode === "under" ? "hidden" : ""; });
  m.querySelectorAll("[data-none]").forEach((e) => (e.textContent = !av[e.dataset.none] ? "(none in this project)" : (e.dataset.none === "scan" && m._w !== "main" ? "(main view only)" : "")));
  const r = btn.getBoundingClientRect(); m.style.left = Math.min(r.left, innerWidth - 190) + "px"; m.style.top = r.bottom + 2 + "px"; m.hidden = false;
}
function heCtlHTML(w) {   // the H&E controls placed in each window's header
  if (w === "slide") return `<span class="hectl" title="this window only: images of the tissues of the slide, under the cells and / or in a second view beside them"><label title="draw the Xenium cells (untick to see the images alone)"><input type="checkbox" data-he="cells" data-hew="slide" checked> Xenium</label> ${heImgsHTML("slide")} <label class="small heop" data-heop="slide" title="how solid the filled cells are over the images drawn under them (0 = see-through, 100 = solid); no effect with outlines ticked or no image under the cells">cell opacity <input type="range" data-he="a" data-hew="slide" min="0" max="100" value="55" style="width:70px;vertical-align:middle"></label></span>`;
  return `<span class="hectl" title="this window only: overlay the images ticked in Images (H&amp;E, IF) under the cells${w === "sig" ? "" : "; outlines draws cell boundaries unfilled"}"><label title="draw the Xenium cells (untick to see the images and outlines alone)"><input type="checkbox" data-he="cells" data-hew="${w}" checked> Xenium</label> ${heImgsHTML(w)} <label class="small heop" data-heop="${w}" title="how solid the filled cells are over the images drawn under them (0 = see-through, 100 = solid); no effect with outlines ticked or no image under the cells">cell opacity <input type="range" data-he="a" data-hew="${w}" min="0" max="100" value="55" style="width:70px;vertical-align:middle"></label>${w === "sig" ? "" : ` <label><input type="checkbox" data-he="outl" data-hew="${w}"> outlines</label>`}${HE_SIDE.includes(w) ? ` <span class="ptwin" data-ptw="${w}" hidden title="arrows you place on the tissue: the same arrow appears in the cells and in the image beside them"><button type="button" class="ptwbtn" data-act="tool" data-w="${w}">Pointer</button><button type="button" class="ptwbtn" data-act="undo" data-w="${w}">Undo</button><button type="button" class="ptwbtn" data-act="clear" data-w="${w}" title="remove the arrows of all tissues">Clear</button></span> <label class="small" data-hlink="${w}" hidden title="zoom and pan the images and the cells together"><input type="checkbox" data-he="link" data-hew="${w}" checked> link zoom</label>` : ""}</span>`;
}
function heSync(w) {   // show window w's settings in its controls
  const st = hw(w);
  document.querySelectorAll(`[data-hew="${w}"]`).forEach((el) => { const k = el.dataset.he; if (el.type === "checkbox") el.checked = !!st[k]; else el.value = Math.round(st[k] * 100); });
  document.querySelectorAll(`[data-hlink="${w}"]`).forEach((e) => (e.hidden = !st.side));
  document.querySelectorAll(`[data-heop="${w}"]`).forEach((e) => { const off = !st.on || (st.outl && w !== "slide"); e.classList.toggle("dim", off); e.querySelector("input").disabled = off; });   // the opacity of filled cells only matters with an image under them and no outlines
  document.querySelectorAll(`[data-ptw="${w}"]`).forEach((e) => (e.hidden = !ptAllowedFor(w)));
  const names = [...new Set([...st.imgs, ...st.sides])].sort().map((k) => (k === "scan" ? "scan" : IMG_NAMES[k]));      // "he" before "if"
  document.querySelectorAll(`.heimgsbtn[data-w="${w}"]`).forEach((b) => { b.textContent = names.length ? `Images: ${names.join(" + ")} ▾` : "Images ▾"; b.classList.toggle("on", names.length > 0); });
  if ($("heimgsmenu") && $("heimgsmenu")._w === w && !$("heimgsmenu").hidden) $("heimgsmenu").querySelectorAll("input[data-img]").forEach((i) => (i.checked = (i.dataset.mode === "side" ? st.sides : st.imgs).has(i.dataset.img)));
  if (w === "main" && $("scanrow")) { scanUi(); ptUi(); if (HEX.deck) ptDraw(); }
}
async function heSet(w, key, val) {
  if (key === "img" && val.on) await heReady();
  const st = hw(w);
  if (key === "img") {
    const set = val.side ? st.sides : st.imgs; val.on ? set.add(val.src) : set.delete(val.src);
    st.on = st.imgs.size > 0; st.side = st.sides.size > 0;   // on = something under the cells, side = something in the pane beside them
  }
  else st[key] = key === "a" ? val / 100 : val;
  heSync(w);
  if (w === "main" && key === "outl") for (const k of Object.keys(HEW)) if (k !== "main" && k !== "sig") { HEW[k].outl = val; heSync(k); HE_T[k] && clearTimeout(HE_T[k]); HE_T[k] = setTimeout(() => HE_RENDER[k] && HE_RENDER[k](), 0); }   // the main view's Outlines setting carries over to the windows; each window can still be changed on its own afterwards
  const none = (HEW.main.on || HEW.main.side) && !scanMode() && !(HE.data && HE.data[S.sample]); document.querySelectorAll(".henote").forEach((e) => (e.textContent = none ? (S.meta && S.meta.has_he_slide ? "the slide H&E is not aligned to the cells yet: use the toolbar “H&E slide”" : "none for this tissue") : ""));
  HE_T[w] && clearTimeout(HE_T[w]); HE_T[w] = setTimeout(() => HE_RENDER[w] && HE_RENDER[w](), key === "a" ? 150 : 0);
}
const HE_T = {};
function bindHE() {
  document.querySelectorAll("[data-hectl]").forEach((el) => { el.outerHTML = heCtlHTML(el.dataset.hectl); });
  document.querySelectorAll("[data-heimgs]").forEach((el) => { el.outerHTML = heImgsHTML(el.dataset.heimgs); });   // the main view's Images button in the left panel
  document.addEventListener("click", (e) => { const b = e.target.closest(".heimgsbtn"); if (b) heMenu(b); });
  $("ghcv")._hw = "gh";
  document.querySelectorAll("[data-he]").forEach((el) => el.addEventListener(el.type === "range" ? "input" : "change", () => heSet(el.dataset.hew, el.dataset.he, el.type === "checkbox" ? el.checked : +el.value)));
}
// H&E beside the cells in a tile: a second canvas next to cv that draws the H&E only (cv2._heOnly), with the same geometry; hovering one marks the same spot in the other.
const hePairDims = (dims, w) => (dims ? { W: Math.floor((dims.W - 6) / 2), H: dims.H } : dims);
function hePane(t, cv, w, sample) {   // idempotent: windows whose canvas persists (glom detail) call it on every draw
  cv._ptSample = sample; if (t._he2) t._he2._ptSample = sample;
  if (!cv._ptb) { cv._ptb = true; ptBindCanvas(cv, w); }
  if (!hw(w).side || !HE.data?.[sample]) { if (t._he2) { const pr = cv.parentNode; if (pr.classList.contains("hepair")) pr.replaceWith(cv); t._he2 = null; } return null; }
  if (t._he2 && t._he2.isConnected) return t._he2;
  const pair = document.createElement("div"); pair.className = "hepair"; cv.replaceWith(pair);
  const cv2 = document.createElement("canvas"), m1 = document.createElement("i"), m2 = document.createElement("i"); cv2._heOnly = true; cv2._hw = w; cv2.className = cv.className; cv2.title = "H&E (no cells)";
  m1.className = m2.className = "hept"; m1.hidden = m2.hidden = true; pair.append(cv, cv2, m1, m2); t._he2 = cv2; t._hw = w;
  const bind = (src, dst, m) => {
    src.addEventListener("pointermove", (e) => {
      if (!src._xf || !dst._xf) return; const r = src.getBoundingClientRect(), x = ((e.clientX - r.left) * src.width / r.width - src._xf.ox) / src._xf.s, y = ((e.clientY - r.top) * src.height / r.height - src._xf.oy) / src._xf.s;
      const qx = dst._xf.s * x + dst._xf.ox, qy = dst._xf.s * y + dst._xf.oy; if (qx < 0 || qy < 0 || qx > dst.width || qy > dst.height) { m.hidden = true; return; }
      m.style.left = dst.offsetLeft + qx * dst.clientWidth / dst.width + "px"; m.style.top = dst.offsetTop + qy * dst.clientHeight / dst.height + "px"; m.hidden = false;
    });
    src.addEventListener("pointerleave", () => (m.hidden = true));
  };
  bind(cv, cv2, m2); bind(cv2, cv, m1); cv2._ptSample = sample; ptBindCanvas(cv2, w); return cv2;
}
function heRedoCv(cv) { if (cv._hr) return; cv._hr = requestAnimationFrame(() => { cv._hr = 0; if (cv._redraw) cv._redraw(); }); }   // an H&E tile arrived for a canvas that can redraw itself
function heRedo(el) { if (el._hr) return; el._hr = requestAnimationFrame(() => { el._hr = 0; if (el._dt) el._dt.draw(el._dt.last, el._dt.last2); }); }   // an H&E tile arrived
function heLru(map, cap, drop) {   // evict the least recently used LOADED entries (pending ones are never dropped, so a tile cannot vanish before it is drawn)
  if (map.size <= cap) return;
  [...map.entries()].filter(([, e]) => e.ok || e.bad).sort((p, q) => p[1].t - q[1].t).slice(0, drop).forEach(([k]) => map.delete(k));
}
function heLoad(map, id, url, redo, cap) {   // -> loaded <img> or null (starts the fetch; redo() is called when it arrives)
  let e = map.get(id);
  if (!e) {
    heLru(map, cap, Math.ceil(cap * 0.4));
    e = { img: new Image(), ok: false, bad: false, w: new Set(), t: 0 }; map.set(id, e);
    e.img.onload = () => { e.ok = true; e.t = ++HE.clock; const ws = [...e.w]; e.w.clear(); for (const f of ws) f(); }; e.img.onerror = () => { e.bad = true; e.w.clear(); };
    e.img.src = url;
  }
  e.t = ++HE.clock;
  if (!e.ok && !e.bad && redo) e.w.add(redo);
  return e.ok ? e.img : null;
}
const heTile = (key, z, x, y, redo) => key.startsWith("if:") ? ifTile(key, z, x, y, redo) : heLoad(HE.imgs, `${key}/${z}/${x}_${y}`, `/api/he/${encodeURIComponent(key)}/${z}/${x}_${y}.jpg`, redo, 400);
// Per-glomerulus crop (tools/he_gloms.py): one image per glom instead of several pyramid tiles. Looked up by the glom's first cell index.
function heCrop(sample, cells) {
  const m = (HE.cropMap && HE.cropMap[sample]) || (HE.cropMap && (HE.cropMap[sample] = new Map((HE.crops[sample] || []).map((r) => [r.i, r])))); if (!m || !m.size) return null;
  let f = Infinity; for (const i of cells) if (i < f) f = i; return m.get(f) || null;
}
const heCropImg = (sample, r, redo, small) => { const f = small && r.fs ? r.fs : r.f; return heLoad(HE.cropImgs, `${sample}/${f}`, `/api/heglom/${encodeURIComponent(sample)}/${f}`, redo, 80); };   // small = half-resolution copy, enough for the ~330 px tiles
// he = { sample, redo }: X / Y map um to canvas pixels, kz = canvas px per um. Tiles not loaded yet are drawn from a coarser level when one is available.
function drawHE(ctx, he, X, Y, W, H, kz) {
  const list = heSel(he.sample, he.w, he.side); if (!list.length) return;
  ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "medium";
  const ux = (sx) => (sx - X(0)) / (X(1) - X(0)), uy = (sy) => (sy - Y(0)) / (Y(1) - Y(0));
  if (!list.some((im) => !im.key.startsWith("if:"))) for (const im of list) {   // IF alone: light added to black (screen over a light page would stay white)
    ctx.globalCompositeOperation = "source-over"; ctx.fillStyle = "#000"; ctx.fillRect(X(im.x0), Y(im.y0), im.w * im.px * (X(1) - X(0)), im.h * im.px * (Y(1) - Y(0))); break;
  }
  for (const im of list) {
    ctx.globalCompositeOperation = im.key.startsWith("if:") ? "screen" : "source-over";   // IF adds its light to what is under it
    const T = im.tsize, u = im.px, z = Math.max(0, Math.min(im.nz - 1, Math.floor(Math.log2(1 / (u * kz))))), span = T * u * 2 ** z;
    const t0 = Math.max(0, Math.floor((ux(0) - im.x0) / span)), t1 = Math.min(Math.ceil(im.w * u / span) - 1, Math.floor((ux(W) - im.x0) / span));
    const r0 = Math.max(0, Math.floor((uy(0) - im.y0) / span)), r1 = Math.min(Math.ceil(im.h * u / span) - 1, Math.floor((uy(H) - im.y0) / span));
    for (let ty = r0; ty <= r1; ty++) for (let tx = t0; tx <= t1; tx++) {
      let img = heTile(im.key, z, tx, ty, he.redo), lz = z, sx = 0, sy = 0, sw = 0, sh = 0;
      if (!img && he.miss) he.miss.n++;   // a figure waits until every tile it needs has arrived
      if (img) { sw = img.naturalWidth; sh = img.naturalHeight; }
      else for (let zz = z + 1; zz < im.nz; zz++) {   // coarser stand-in: the sub-rectangle of the parent tile that this tile covers
        const d = zz - z, ptx = tx >> d, pty = ty >> d, p = heTile(im.key, zz, ptx, pty, null); if (!p) continue;
        sx = (tx - (ptx << d)) * (T >> d); sy = (ty - (pty << d)) * (T >> d); sw = Math.min(T >> d, p.naturalWidth - sx); sh = Math.min(T >> d, p.naturalHeight - sy);
        if (sw > 0 && sh > 0) { img = p; lz = zz; break; }
      }
      if (!img) continue;
      const f = u * 2 ** lz;   // um per source pixel
      ctx.drawImage(img, sx, sy, sw, sh, X(im.x0 + tx * span), Y(im.y0 + ty * span), sw * f * kz + 0.6, sh * f * kz + 0.6);
    }
  }
  ctx.globalCompositeOperation = "source-over";
}
// Compare views draw real cell outlines once few enough cells are in view (zoomed in); otherwise squares at the cell centres.
const POLY_MAX = 9000;
function drawCells(ctx, order, b, X, Y, W, H, px, poly, a = 1, ol = false, q = 1, pmax = POLY_MAX) {   // q = render scale (line widths); px is already scaled   // a = opacity of filled cells (< 1 lets the H&E show through)
  const m = 40; let vis = 0;
  if (poly) for (const [, list] of order) for (const i of list) { const sx = X(b.xy[i * 2]), sy = Y(b.xy[i * 2 + 1]); if (sx > -m && sy > -m && sx < W + m && sy < H + m && ++vis > pmax) break; }
  if (poly && vis <= pmax) {
    ctx.lineJoin = "round";
    for (const [col, list] of order) {
      if (ol) { ctx.globalAlpha = 1; ctx.lineWidth = 1.6 * q; ctx.strokeStyle = col.startsWith("rgba") ? "rgba(150,155,168,0.55)" : col; }
      else { ctx.globalAlpha = a; ctx.lineWidth = (S.grey ? 0.9 : 0.8) * q; ctx.strokeStyle = S.grey ? "rgba(128,134,143,0.92)" : "rgba(10,10,10,0.75)"; ctx.fillStyle = col; }
      for (const i of list) {
        const sx = X(b.xy[i * 2]), sy = Y(b.xy[i * 2 + 1]); if (sx < -m || sy < -m || sx > W + m || sy > H + m) continue;
        ctx.beginPath(); for (let q = poly.off[i] * 2, f = true; q < poly.off[i + 1] * 2; q += 2, f = false) { const x = X(poly.xy[q]), y = Y(poly.xy[q + 1]); f ? ctx.moveTo(x, y) : ctx.lineTo(x, y); }
        ctx.closePath(); if (!ol) ctx.fill(); ctx.stroke();
      }
    }
    ctx.globalAlpha = 1; return true;
  }
  ctx.globalAlpha = a;
  for (const [col, list] of order) { ctx.fillStyle = col; for (const i of list) { const sx = X(b.xy[i * 2]), sy = Y(b.xy[i * 2 + 1]); if (sx < -px || sy < -px || sx > W || sy > H) continue; ctx.fillRect(sx, sy, px, px); } }
  ctx.globalAlpha = 1; return false;
}
// view = { z, rx, ry }: zoom factor and the centre as a fraction of the tissue's bounding box (null = whole tissue). cache keeps the colour groups between redraws.
function drawTile(cv, s, b, sc, box = GAL_BOX, view = null, cache = null, dims = null, poly = null, he = null, w8 = "gal", q = 1) {   // dims = { W, H }: fill a fixed rectangle instead of the tissue's own aspect ratio; q = render scale of the export (points and lines grow with it)
  const w = s.x[1] - s.x[0] || 1, h = s.y[1] - s.y[0] || 1, k = dims ? Math.min((dims.W - 8) / w, (dims.H - 8) / h) : Math.min((box - 8) / w, (box - 8) / h);
  cv.width = dims ? dims.W : Math.max(40, Math.ceil(w * k) + 8); cv.height = dims ? dims.H : Math.max(40, Math.ceil(h * k) + 8);
  const ctx = cv.getContext("2d"), n = b.n, z = view ? view.z : 1, W = cv.width, H = cv.height, kz = k * z;
  const cx = s.x[0] + (view ? view.rx : 0.5) * w, cy = s.y[0] + (view ? view.ry : 0.5) * h, px = (n > 30000 ? 1.4 : 2) * Math.min(z, 5) * q;
  ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, W, H);
  const X = (x) => W / 2 + (x - cx) * kz, Y = (y) => H / 2 + (y - cy) * kz;
  cv._xf = { s: kz, ox: W / 2 - cx * kz, oy: H / 2 - cy * kz };
  if (he) drawHE(ctx, { ...(he.w ? he : { ...he, w: w8 }), side: !!cv._heOnly }, X, Y, W, H, kz);
  if (cv._heOnly || !hw(w8).cells) return;   // the H&E beside the cells, or Xenium switched off: no cells
  let order = cache?.order;
  if (!order) {
    const groups = new Map();  // colour -> cells (fewer fillStyle switches)
    for (let i = 0; i < n; i++) {
      const col = b.pass[i] ? (S.grey ? rgbStr(GREY) : rgbStr(sc.color(b, i))) : "rgba(105,110,122,0.45)";
      (groups.get(col) || groups.set(col, []).get(col)).push(i);
    }
    order = [...groups.entries()].sort((p, q) => (p[0].startsWith("rgba") ? -1 : 0) - (q[0].startsWith("rgba") ? -1 : 0));  // grey first
    if (cache) cache.order = order;
  }
  drawCells(ctx, order, b, X, Y, W, H, px, poly, he ? he.a : 1, hw(w8).outl, q);   // cells sit on the H&E at the chosen opacity
}
// Compare views (Gallery and Regions, "only ticked"): the window is split into one chunk per ticked tissue; wheel = zoom at the cursor, drag = pan,
// double-click = reset; "link zoom" applies one view to every tile. rec = { bb: {x0, y0, w, h}, draw(view) }; cfg = { views(), sync(), tiles() }.
function cmpLayout(winId, gridId, n, capH = 62) {
  const win = $(winId), grid = $(gridId), gap = 10, cols = n <= 3 ? Math.max(1, n) : Math.ceil(Math.sqrt(n)), rows = Math.max(1, Math.ceil(n / cols));
  const top = grid.getBoundingClientRect().top - win.getBoundingClientRect().top + win.scrollTop, availH = Math.max(240, win.clientHeight - top - 14), cw = grid.clientWidth;
  grid.style.cssText = `display:grid;gap:${gap}px;grid-template-columns:repeat(${cols},minmax(0,1fr));grid-template-rows:repeat(${rows},${Math.floor((availH - gap * (rows - 1)) / rows)}px);height:${availH}px`;
  return { W: Math.max(120, Math.floor((cw - gap * (cols - 1)) / cols) - 16), H: Math.max(100, Math.floor((availH - gap * (rows - 1)) / rows) - capH) };
}
const ID_VIEW = { z: 1, rx: 0.5, ry: 0.5 };
function attachZoom(t, rec, cfg) {
  const clamp = (q) => Math.max(0, Math.min(1, q)), key = () => (cfg.sync() ? "__all" : t.dataset.id), hk = (pane) => (pane === 1 && !hw(t._hw).link ? key() + "#h" : key());   // pane 1 = the H&E beside the cells
  const redraw = (el) => { const r = el._dt; if (!r) return; const k = cfg.sync() ? "__all" : el.dataset.id, v = cfg.views()[k] || ID_VIEW; r.draw(v, cfg.views()[k + "#h"] || ID_VIEW); el.querySelector(".zb").textContent = v.z > 1.01 ? "×" + v.z.toFixed(v.z < 10 ? 1 : 0) : ""; };
  const redrawSoon = (el) => { if (el._raf) return; el._raf = requestAnimationFrame(() => { el._raf = 0; redraw(el); }); };
  const bind = (cv, pane) => {
    if (cv._zb) return; cv._zb = 1;   // handlers are bound once per canvas
    const view = () => cfg.views()[hk(pane)] || ID_VIEW, apply = (v) => { cfg.views()[hk(pane)] = v; if (cfg.sync()) cfg.tiles().forEach(redrawSoon); else redrawSoon(t); };
    const geom = () => { const bb = t._dt.bb; return { w: bb.w, h: bb.h, k: Math.min((cv.width - 8) / bb.w, (cv.height - 8) / bb.h) }; };
    const pos = (e) => { const r = cv.getBoundingClientRect(); return { X: (e.clientX - r.left) * cv.width / r.width, Y: (e.clientY - r.top) * cv.height / r.height }; };
    cv.style.cursor = "grab"; cv.style.touchAction = "none";
    cv.addEventListener("wheel", (e) => {
      e.preventDefault(); const v = view(), g = geom(), { X, Y } = pos(e), z2 = Math.max(1, Math.min(60, v.z * Math.exp(-e.deltaY * 0.0018)));
      const wx = (X - cv.width / 2) / (g.k * v.z) + v.rx * g.w, wy = (Y - cv.height / 2) / (g.k * v.z) + v.ry * g.h;   // tissue point under the cursor (relative to the box corner)
      apply({ z: z2, rx: z2 === 1 ? 0.5 : clamp((wx - (X - cv.width / 2) / (g.k * z2)) / g.w), ry: z2 === 1 ? 0.5 : clamp((wy - (Y - cv.height / 2) / (g.k * z2)) / g.h) });
    }, { passive: false });
    let drag = null;
    cv.addEventListener("pointerdown", (e) => { t._moved = false; drag = { x: e.clientX, y: e.clientY, v: view() }; cv.setPointerCapture(e.pointerId); cv.style.cursor = "grabbing"; });
    cv.addEventListener("pointermove", (e) => { if (!drag) return; if (Math.hypot(e.clientX - drag.x, e.clientY - drag.y) > 3) t._moved = true; const g = geom(), sc = cv.width / cv.getBoundingClientRect().width, v = drag.v;
      apply({ z: v.z, rx: clamp(v.rx - ((e.clientX - drag.x) * sc) / (g.k * v.z * g.w)), ry: clamp(v.ry - ((e.clientY - drag.y) * sc) / (g.k * v.z * g.h)) }); });
    cv.addEventListener("pointerup", () => { drag = null; cv.style.cursor = "grab"; });
    cv.addEventListener("dblclick", () => apply({ ...ID_VIEW }));
  };
  bind(t.querySelector("canvas"), 0); if (t._he2) bind(t._he2, 1);
}
function resetZoom(cfg) { for (const k of Object.keys(cfg.views())) delete cfg.views()[k]; cfg.tiles().forEach((el) => { if (el._dt) { el._dt.draw(null); el.querySelector(".zb").textContent = ""; } }); }
const galCfg = { views: () => S.galView || (S.galView = {}), sync: () => $("galsync").checked, tiles: () => $("galgrid").querySelectorAll(".tile") };
const rgnCfg = { views: () => S.rgnView || (S.rgnView = {}), sync: () => $("rgnsync").checked, tiles: () => $("rgngrid").querySelectorAll(".tile") };

// ---------------------------------------------------------------- glomerulus window
// One tile per glomerulus ("All gloms#N <tissue>") showing its real cell boundaries, coloured by the current colour-by / filters.
function galPaths(sample) {   // pathologist region polygons of a sample (null when none were matched)
  const c = galData[sample] || (galData[sample] = { cats: {} });
  return c.paths || (c.paths = jget(`/api/paths/${encodeURIComponent(sample)}`).then((d) => d.features).catch(() => null));
}
// Boundary of one glomerulus: the "All gloms" polygon whose region id equals the glom's (glom id "All gloms#59 403-019" -> region id "All gloms#59").
function glomOutline(features, label) {
  if (!features) return null; const rid = label.replace(/ \d{3}-\d{3}$/, "");
  const fs = features.filter((f) => f.properties.annotation_name === "All gloms" && f.properties.region_id === rid && !f.properties.negative_roa); if (!fs.length) return null;
  return { rings: fs.flatMap((f) => f.geometry.coordinates), color: regionColor("All gloms", fs[0].properties.line_color_hex) };
}
function galPoly(sample) {
  const c = galData[sample] || (galData[sample] = { cats: {} }), q = `/api/sample/${encodeURIComponent(sample)}/poly/`;
  return c.poly || (c.poly = Promise.all([bget(q + "offsets", Int32Array), bget(q + "xy", Float32Array)]).then(([off, xy]) => ({ off, xy })));
}
function openGlom() {
  closeAllWindows();
  const subs = sortedSubtypes();
  $("glomscope").innerHTML = `<option value="cur">This donor (${dlabel(S.sample)})</option><option value="all">All donors</option>` + subs.map((t) => `<option value="sub:${t}">Disease subtype: ${t}</option>`).join("");
  if (!$("glom").dataset.init) { $("glomscope").value = "cur"; $("glom").dataset.init = "1"; } else $("glomscope").value = S.glomScope || "cur";
  $("glom").hidden = false; renderGlom();
}
function closeGlom() { $("glom").hidden = true; if ($("glomzoom")) $("glomzoom").hidden = true; }
let glomToken = 0;
async function renderGlom() {
  if ($("glom").hidden) return;
  const tok = ++glomToken, scope = $("glomscope").value, qcSel = $("glomqc").value; S.glomScope = scope;
  const only = $("glomonly").checked; $("glomn").textContent = S.glomSel.size; $("glomscope").disabled = $("glomqc").disabled = only; $("glomzoomctl").hidden = !only; $("glomgrid").classList.toggle("cmpgrid", only); if (!only) $("glomgrid").removeAttribute("style");
  if (only) return renderGlomOnly(tok);
  const rank = (t) => sortedSubtypes().indexOf(t);
  const samples = S.meta.samples.filter((s) => scope === "all" || (scope === "cur" ? s.id === S.sample : subtypeOf(s) === scope.slice(4)))
    .sort((x, y) => rank(subtypeOf(x)) - rank(subtypeOf(y)));  // stable: keeps sample order within a subtype
  const totals = { passed: 0, failed: 0, drawn: 0 }, qcLv = catMeta("glom_quality").levels, gidLv = catMeta("glom_id").levels;
  const info = (done) => { $("glominfo").textContent = `${totals.passed + totals.failed} glomeruli (${totals.passed} passed qc, ${totals.failed} failed qc) · outline = pathologist boundary, passed-qc gloms only (drawn for ${totals.drawn}) · colour: ${colorName()}${winNote()}${done ? "" : " · loading…"} · click one to open it`; };
  $("glomgrid").innerHTML = samples.map((s) => `<div class="glomsec" data-s="${s.id}" hidden><div class="sechead"></div><div class="secgrid"></div></div>`).join("");
  info(false);
  const bp = samples.map((s) => winBundle(s.id)), gp = samples.map((s) => Promise.all([galGet(s.id, "glom_id"), galGet(s.id, "glom_quality"), galPoly(s.id), galPaths(s.id)]));
  const sc = await winScale(bp, () => glomToken, tok).catch((e) => { console.error(e); return null; });
  if (!sc || tok !== glomToken) return;
  $("glomleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  const shown = S.ctxs.glom = { sc, shown: [] };
  let pending = samples.length;
  await Promise.all(samples.map(async (s, k) => {
    try {
      const [b, [gid, gq, poly, feats]] = await Promise.all([bp[k], gp[k]]);
      if (tok !== glomToken) return;
      const gl = new Map();  // glom code -> cell indices
      for (let i = 0; i < b.n; i++) if (gid[i] !== MISSING) (gl.get(gid[i]) || gl.set(gid[i], []).get(gid[i])).push(i);
      const rows = [...gl.entries()].filter(([, cells]) => cells.some((i) => b.pass[i])).map(([g, cells]) => ({ g, cells, qc: qcLv[gq[cells[0]]] }))
        .filter((r) => qcSel === "all" || r.qc === qcSel)
        .sort((x, y) => +gidLv[x.g].match(/#(\d+)/)[1] - +gidLv[y.g].match(/#(\d+)/)[1]);
      const sec = $("glomgrid").querySelector(`.glomsec[data-s="${CSS.escape(s.id)}"]`); if (!sec || !rows.length) return;
      sec.hidden = false;
      sec.querySelector(".sechead").textContent = `${dlabel(s.id)} · ${s.disease_type} · ${subtypeOf(s)} · ${rows.length} glomeruli`;
      const grid = sec.querySelector(".secgrid");
      for (const r of rows) {
        r.qc === "passed qc" ? totals.passed++ : totals.failed++;
        const t = document.createElement("div"), gkey = `${s.id}|${gidLv[r.g]}`; t.className = "gtile" + (S.glomSel.has(gkey) ? " sel" : ""); t.dataset.key = gkey; t.title = `${gidLv[r.g]} — click to zoom; tick the box to compare`;
        t.innerHTML = `<input type="checkbox" class="gsel" title="tick to compare this glom" ${S.glomSel.has(gkey) ? "checked" : ""}><canvas width="328" height="328"></canvas><div class="gid">${gidLv[r.g]}</div><div class="gm">${subtypeOf(s)}</div><div class="gm"><span class="chip ${r.qc === "passed qc" ? "pass" : "fail"}">${r.qc === "passed qc" ? "✓" : "✕"} ${r.qc}</span>${fmt(r.cells.length)} cells</div>`;
        const outline = r.qc === "passed qc" ? glomOutline(feats, gidLv[r.g]) : null; totals.drawn += outline ? 1 : 0;   // boundary only for gloms that passed qc
        t.addEventListener("click", () => openGlomZoom({ src: "glom", s, b, sc, poly, cells: r.cells, label: gidLv[r.g], qc: r.qc, outline, note: outline ? "outline = pathologist \"All gloms\" boundary" : undefined }));
        const cb = t.querySelector(".gsel"); cb.addEventListener("click", (e) => e.stopPropagation());
        cb.addEventListener("change", () => { glomTick(gkey, { sample: s.id, label: gidLv[r.g] }, cb.checked); });
        grid.appendChild(t);
        drawGlom(Object.assign(t.querySelector("canvas"), { _hw: "glom" }), b, sc, poly, r.cells, true, -1, outline); shown.shown.push({ b, cells: r.cells, el: t, key: gkey, rec: { sample: s.id, label: gidLv[r.g] } });
      }
    } catch (e) { console.error(e); }
    if (tok === glomToken) { pending--; info(pending === 0); }
  }));
  if (tok === glomToken) { info(true); updateGlomSelUI(); if (!totals.passed && !totals.failed) $("glomgrid").innerHTML = `<div class="muted">No glomeruli match the current filters.</div>`; }
}
// "only ticked": the window is split into one equal chunk per ticked glom, each zoomable (wheel), pannable (drag), double-click = reset; "link zoom" shares one view.
const glomCfg = { views: () => S.glomView || (S.glomView = {}), sync: () => $("glomsync").checked, tiles: () => $("glomgrid").querySelectorAll(".gtile") };
async function renderGlomOnly(tok) {
  const grid = $("glomgrid"), order = [...S.glomSel.keys()], bySample = new Map();
  for (const [k, r] of S.glomSel) (bySample.get(r.sample) || bySample.set(r.sample, []).get(r.sample)).push([k, r.label]);
  const samples = S.meta.samples.filter((s) => bySample.has(s.id));
  $("glomleg").innerHTML = ""; $("glominfo").textContent = "";
  if (!samples.length) { grid.innerHTML = `<div class="muted">No glomerulus is ticked. Untick “only ticked” and tick the gloms to compare (or use “Add to comparison” in Glom heterogeneity).</div>`; return; }
  const gidLv = catMeta("glom_id").levels, qcLv = catMeta("glom_quality").levels, bp = samples.map((s) => winBundle(s.id));
  const sc = await winScale(bp, () => glomToken, tok).catch((e) => { console.error(e); return null; }); if (!sc || tok !== glomToken) return;
  const items = [];
  await Promise.all(samples.map(async (s, q) => {
    const [b, gid, gq, poly, feats] = await Promise.all([bp[q], galGet(s.id, "glom_id"), galGet(s.id, "glom_quality"), galPoly(s.id), galPaths(s.id)]);
    for (const [gkey, label] of bySample.get(s.id)) {
      const code = gidLv.indexOf(label), cells = []; if (code < 0) continue;
      for (let i = 0; i < b.n; i++) if (gid[i] === code) cells.push(i); if (!cells.length) continue;
      const qc = qcLv[gq[cells[0]]]; items.push({ s, b, poly, cells, gkey, label, qc, outline: qc === "passed qc" ? glomOutline(feats, label) : null });
    }
  }));
  if (tok !== glomToken) return; items.sort((x, y) => order.indexOf(x.gkey) - order.indexOf(y.gkey));
  $("glomleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  $("glominfo").textContent = `${items.length} ticked glomeruli · colour: ${colorName()}${winNote()} · window split per glom · wheel = zoom, drag = pan, double-click = reset, click a name to open the glom`;
  grid.innerHTML = items.map((r) => `<div class="gtile sel" data-id="${r.gkey}" data-key="${r.gkey}"><input type="checkbox" class="gsel" title="untick to remove from this view" checked><span class="zb"></span><canvas width="300" height="300"></canvas>
    <div class="cap"><b class="gid">${r.label}</b><div class="gm"><span class="chip ${r.qc === "passed qc" ? "pass" : "fail"}">${r.qc === "passed qc" ? "✓" : "✕"} ${r.qc}</span>${dlabel(r.s.id)} · ${subtypeOf(r.s)} · ${fmt(r.cells.length)} cells</div></div></div>`).join("");
  const dims = cmpLayout("glom", "glomgrid", Math.max(1, items.length), 70);
  items.forEach((r) => {
    const t = grid.querySelector(`.gtile[data-id="${CSS.escape(r.gkey)}"]`), cv = t.querySelector("canvas"), B = glomBox(r.cells, r.poly, 12, r.outline); cv._hw = "glom";
    const cv2 = hePane(t, cv, "glom", r.s.id), dm = cv2 ? hePairDims(dims, "glom") : dims;
    const rec = t._dt = { bb: { x0: B.x0, y0: B.y0, w: B.x1 - B.x0, h: B.y1 - B.y0 }, last: null, last2: null, draw: (vw, vh) => { rec.last = vw; rec.last2 = vh; drawGlom(cv, r.b, sc, r.poly, r.cells, true, -1, r.outline, () => heRedo(t), vw, dm); if (cv2) drawGlom(cv2, r.b, sc, r.poly, r.cells, true, -1, r.outline, () => heRedo(t), HEW.glom.link ? vw : vh, dm); } };
    const gk = glomCfg.sync() ? "__all" : r.gkey, v = glomCfg.views()[gk] || null; rec.draw(v, glomCfg.views()[gk + "#h"] || null); attachZoom(t, rec, glomCfg); if (v && v.z > 1.01) t.querySelector(".zb").textContent = "×" + v.z.toFixed(1);
    t.querySelector(".gsel").addEventListener("change", (e) => { glomTick(r.gkey, { sample: r.s.id, label: r.label }, e.target.checked); renderGlom(); });
    t.querySelector(".gid").addEventListener("click", () => openGlomZoom({ src: "glom", s: r.s, b: r.b, sc, poly: r.poly, cells: r.cells, label: r.label, qc: r.qc, outline: r.outline, note: r.outline ? "outline = pathologist \"All gloms\" boundary" : undefined }));
  });
}
function glomBox(cells, poly, margin = 12, outline = null) {
  let x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18;
  for (const i of cells) for (let k = poly.off[i] * 2; k < poly.off[i + 1] * 2; k += 2) {
    const x = poly.xy[k], y = poly.xy[k + 1]; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
  }
  if (outline) for (const ring of outline.rings) for (const [x, y] of ring) { if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }   // the whole boundary stays in view
  return { x0: x0 - margin, x1: x1 + margin, y0: y0 - margin, y1: y1 + margin };
}
// view = { z, rx, ry } (zoom and centre as a fraction of the glom's box) and dims = { W, H } give the zoomable compare tile; without them: the square, fitted tile.
function drawGlom(cv, b, sc, poly, cells, context = true, focus = -1, outline = null, onHe = null, view = null, dims = null) {  // returns pixels per µm; `focus` = cell outlined in white; `outline` = pathologist boundary
  const B = glomBox(cells, poly, 12, outline), bw = B.x1 - B.x0, bh = B.y1 - B.y0;
  if (dims) { cv.width = dims.W; cv.height = dims.H; }
  const W = cv.width, H = cv.height, k = Math.min((W - 8) / bw, (H - 8) / bh) * (view ? view.z : 1);
  const ox = W / 2 - (view ? view.rx : 0.5) * bw * k, oy = H / 2 - (view ? view.ry : 0.5) * bh * k;
  const vx0 = B.x0 - ox / k, vx1 = vx0 + W / k, vy0 = B.y0 - oy / k, vy1 = vy0 + H / k;   // visible window (um)
  const ctx = cv.getContext("2d", { willReadFrequently: true }); ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, W, H); ctx.lineWidth = Math.max(1, W / 328); ctx.lineJoin = "round";   // CPU canvas: hundreds of small tiles would otherwise each get a GPU surface
  cv._redraw = () => drawGlom(cv, b, sc, poly, cells, context, focus, outline, onHe, view, dims);
  const st = hw(cv._hw), hs = cv._heOnly ? (HE.data?.[b.sample] ? b.sample : null) : heFor(b.sample, cv._hw), al = hs ? st.a : 1, ol = st.outl, lw = Math.max(1, Math.min(W, H) / 328);
  if (hs) {
    const redo = onHe || (() => heRedoCv(cv)), cr = (cv._heOnly ? st.sides : st.imgs).has("he") ? heCrop(hs, cells) : null, sl = 1;
    if (cr && vx0 >= cr.x0 - sl && vx1 <= cr.x0 + cr.w + sl && vy0 >= cr.y0 - sl && vy1 <= cr.y0 + cr.h + sl) {   // the glom's own crop covers the window: one image, one draw
      const im = heCropImg(hs, cr, redo, k <= 2.6); if (im) { ctx.imageSmoothingQuality = "low"; ctx.drawImage(im, ox + (cr.x0 - B.x0) * k, oy + (cr.y0 - B.y0) * k, cr.w * k, cr.h * k); }
    } else drawHE(ctx, { sample: hs, redo, w: cv._hw, side: !!cv._heOnly }, (x) => ox + (x - B.x0) * k, (y) => oy + (y - B.y0) * k, W, H, k);
  }
  cv._xf = { s: k, ox: ox - B.x0 * k, oy: oy - B.y0 * k };
  const paint = (i, fill, stroke) => { path(i); if (ol) { ctx.lineWidth = lw * 1.8; ctx.strokeStyle = stroke; ctx.stroke(); } else { ctx.globalAlpha = al; ctx.fillStyle = fill; ctx.fill(); ctx.lineWidth = lw; ctx.strokeStyle = "rgba(10,10,10,0.8)"; ctx.stroke(); ctx.globalAlpha = 1; } };
  const path = (i) => { ctx.beginPath(); for (let q = poly.off[i] * 2, f = true; q < poly.off[i + 1] * 2; q += 2, f = false) { const X = ox + (poly.xy[q] - B.x0) * k, Y = oy + (poly.xy[q + 1] - B.y0) * k; f ? ctx.moveTo(X, Y) : ctx.lineTo(X, Y); } ctx.closePath(); };
  const inGlom = new Set(cells);
  const ho = cv._heOnly || !st.cells;   // no cells: H&E beside the cells, or Xenium switched off
  for (let i = 0; context && !ho && i < b.n; i++) {   // neighbouring (non-glom) cells for context
    if (inGlom.has(i)) continue; const x = b.xy[i * 2], y = b.xy[i * 2 + 1];
    if (view || dims ? x >= vx0 - 10 && x <= vx1 + 10 && y >= vy0 - 10 && y <= vy1 + 10 : x >= B.x0 && x <= B.x1 && y >= B.y0 && y <= B.y1) paint(i, "rgba(105,110,122,0.22)", "rgba(150,155,168,0.45)");
  }
  for (const i of cells) if (!ho && !b.pass[i]) paint(i, "rgba(105,110,122,0.5)", "rgba(150,155,168,0.7)");
  for (const i of cells) if (!ho && b.pass[i]) { const c = rgbStr(sc.color(b, i)); paint(i, c, c); }
  if (outline) {   // pathologist boundary: dark halo, then the pen colour
    const ring = (r) => { ctx.beginPath(); r.forEach(([px, py], j) => { const X = ox + (px - B.x0) * k, Y = oy + (py - B.y0) * k; j ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); }); ctx.closePath(); };
    for (const [w, col] of [[Math.max(4.5, Math.min(W, H) / 75), "rgba(0,0,0,0.85)"], [Math.max(2, Math.min(W, H) / 160), outline.color]]) { ctx.lineWidth = w; ctx.strokeStyle = col; for (const r of outline.rings) { ring(r); ctx.stroke(); } }
  }
  if (focus >= 0) { path(focus); ctx.lineWidth = Math.max(2.5, Math.min(W, H) / 90); ctx.strokeStyle = "#ffffff"; ctx.stroke(); }
  if (dims) {   // scale bar (compare tiles)
    const target = 0.2 * (W - 8) / k, p = Math.pow(10, Math.floor(Math.log10(target))), nice = [1, 2, 5, 10].map((f) => f * p).reduce((a, c) => (Math.abs(c - target) < Math.abs(a - target) ? c : a));
    const lab = nice >= 1000 ? `${nice / 1000} mm` : `${nice} µm`; ctx.font = `11px ${FONT}`; ctx.textBaseline = "alphabetic";
    ctx.fillStyle = "rgba(15,18,24,0.75)"; ctx.fillRect(6, H - 32, Math.max(nice * k, ctx.measureText(lab).width) + 10, 28); ctx.fillStyle = "#fff"; ctx.fillRect(10, H - 12, nice * k, 3); ctx.fillText(lab, 10, H - 17);
  }
  return k;
}
// ---- zoomed single-glomerulus view (what you see is what "Save PNG" writes)
const FONT = "-apple-system, Segoe UI, Roboto, sans-serif";
// Legend drawn into a canvas context (for exported PNGs). `items` = [[name, css colour, count]] for categorical colouring.
// u = size multiplier; returns the height used so callers can size the figure (call once on a scratch context to measure).
function pngLegend(x, info, items, X, Y, W, u = 1, maxLines = 3) {
  const mode = S.colorBy; x.textBaseline = "alphabetic"; x.textAlign = "left";
  if (info.type === "cat") {
    x.font = `${14 * u}px ${FONT}`; let px = X, py = Y, lines = 1;
    for (const [name, col, n] of items) {
      const t = `${name} (${fmt(n)})`, w = x.measureText(t).width + 26 * u;
      if (px + w > X + W) {
        px = X; py += 22 * u; lines++;
        if (lines > maxLines) { x.fillStyle = "#8b95a5"; x.fillText("…", px, py); return lines * 22 * u; }
      }
      x.fillStyle = col; x.fillRect(px, py - 11 * u, 12 * u, 12 * u); x.fillStyle = "#d7dde6"; x.fillText(t, px + 18 * u, py); px += w;
    }
    return lines * 22 * u;
  }
  x.font = `${14 * u}px ${FONT}`; x.fillStyle = "#d7dde6";
  if (info.type === "dual") {
    const a = `${S.geneA} (to ${info.hiA.toFixed(2)})`, b = `${S.geneB} (to ${info.hiB.toFixed(2)})`;
    x.fillStyle = "#d33"; x.fillRect(X, Y - 11 * u, 12 * u, 12 * u); x.fillStyle = "#d7dde6"; x.fillText(a, X + 18 * u, Y);
    const w = 18 * u + x.measureText(a).width + 20 * u;
    x.fillStyle = "#3a3"; x.fillRect(X + w, Y - 11 * u, 12 * u, 12 * u); x.fillStyle = "#d7dde6"; x.fillText(b, X + w + 18 * u, Y);
    x.fillStyle = "#8b95a5"; x.fillText(`yellow = both · gray = neither detected${transformed() ? " · positive part of " + transformNote() : ""}`, X, Y + 22 * u);
    return 34 * u;
  }
  if (info.type === "sigdual") {
    let px = X; x.font = `${14 * u}px ${FONT}`;
    for (const [col, name] of [[SIG_A, info.labelA], [SIG_B, info.labelB || "(no signature B)"]]) { x.fillStyle = col; x.fillRect(px, Y - 11 * u, 12 * u, 12 * u); x.fillStyle = "#d7dde6"; x.fillText(name, px + 18 * u, Y); px += 18 * u + x.measureText(name).width + 26 * u; }
    x.fillStyle = "#8b95a5"; x.fillText(`near white = both high · intensity: score ${info.z ? "above the average" : "from the 1st percentile"} up to p${S.clip}, each signature scaled separately${transformed() ? " · " + transformNote() : ""}`, X, Y + 22 * u);
    return 34 * u;
  }
  if (info.type === "none") return 0;
  const L = info.type === "gene" && info.cont && info.div ? DIVL : LUT, bw = 300 * u, gr = x.createLinearGradient(X, 0, X + bw, 0);
  L.forEach((c, i) => { if (i % 16 === 0 || i === 255) gr.addColorStop(i / 255, `rgb(${c})`); });
  x.fillStyle = gr; x.fillRect(X, Y - 11 * u, bw, 12 * u); x.fillStyle = "#d7dde6";
  for (const k of scaleTicks(info)) { x.textAlign = k.t === 0 ? "left" : k.t === 1 ? "right" : "center"; x.fillText(fmtTick(k.v), X + bw * k.t, Y + 20 * u); }
  x.textAlign = "left";
  const lay = S.layer === "counts" ? "raw counts" : `log1p(CP${fmt(S.meta.norm_scale)})`;
  x.fillStyle = "#8b95a5";  // description goes under the bar (one or two lines) so it can never run off the figure
  x.fillText(info.type === "num" ? `${labelOf(mode)}${info.log ? " (log scale)" : ""} · ends: ${info.autoLo ? "p1" : "set by you"} / ${info.autoHi ? "p99" : "set by you"}` : `${info.label} — ${info.sig ? "mean of " : ""}${lay}${info.type === "gene" && !info.div ? " · 0 = gray: not detected" : ""}${info.type === "gene" ? " · top = p" + S.clip : ""}`, X, Y + 42 * u);
  if (info.type !== "num" && transformed()) { x.fillText(transformNote(), X, Y + 62 * u); return 70 * u; }
  return 50 * u;
}
let heKz = 0;
function heKickZoom() { if (heKz) return; heKz = requestAnimationFrame(() => { heKz = 0; if (S.gz) drawGlomZoom(); }); }   // an H&E tile arrived for the zoomed glom
function composeGlomPNG(v, context) {
  const G = 900, H = 180, out = document.createElement("canvas"); out.width = G; out.height = G + H;
  const g = document.createElement("canvas"); g.width = g.height = G; g._hw = "gz";
  const k = drawGlom(g, v.b, v.sc, v.poly, v.cells, context, v.focus ?? -1, v.outline, heKickZoom, v.view || null), x = out.getContext("2d"); v.xf = g._xf;
  x.fillStyle = "#0f1218"; x.fillRect(0, 0, G, G + H); x.drawImage(g, 0, 0);
  const target = 0.2 * (G - 8) / k, p = Math.pow(10, Math.floor(Math.log10(target)));
  const nice = [1, 2, 5, 10].map((m) => m * p).reduce((a, b) => (Math.abs(b - target) < Math.abs(a - target) ? b : a));
  x.fillStyle = "#fff"; x.fillRect(24, G - 30, nice * k, 4); x.font = `15px ${FONT}`; x.textBaseline = "alphabetic"; x.fillText(`${nice} µm`, 24, G - 38);
  x.fillStyle = "#d7dde6"; x.font = `600 24px ${FONT}`; x.fillText(v.label, 24, G + 36);
  const tw = x.measureText(v.label).width, chip = v.qc || ""; x.font = `600 14px ${FONT}`;
  if (chip) {
    const chipTxt = (chip === "passed qc" ? "✓ " : "✕ ") + chip;
    x.fillStyle = chip === "passed qc" ? STATUS.good : STATUS.critical; x.fillRect(24 + tw + 14, G + 18, x.measureText(chipTxt).width + 16, 24);
    x.fillStyle = chip === "passed qc" ? "#0b0b0b" : "#fff"; x.fillText(chipTxt, 24 + tw + 22, G + 35);
  }
  x.font = `14px ${FONT}`; x.fillStyle = "#8b95a5";
  x.fillText(`${dlabel(v.s.id)} · ${v.s.disease_type} · ${subtypeOf(v.s)} · ${v.cells.length} cells · colour: ${colorName()}${v.note ? " · " + v.note : ""}`, 24, G + 62);
  let items = null;
  if (v.sc.info.type === "cat") {
    const cnt = new Map(); let dimmed = 0;
    for (const i of v.cells) { if (!v.b.pass[i]) { dimmed++; continue; } cnt.set(v.b.codes[i], (cnt.get(v.b.codes[i]) || 0) + 1); }
    items = [...cnt.entries()].sort((a, b) => b[1] - a[1]).map(([c, n]) => [c === MISSING ? "missing" : display(S.colorBy, c), c === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, c), n]);
    if (dimmed) items.push(["filtered out", "rgba(105,110,122,0.6)", dimmed]);
  }
  pngLegend(x, v.sc.info, items, 24, G + 94, G - 48);
  return out;
}
function drawGlomZoom() {
  const v = S.gz; if (!v) return;
  const c = composeGlomPNG(v, $("gzctx").checked), cv = $("gzcanvas");
  cv.width = c.width; cv.height = c.height; cv.getContext("2d").drawImage(c, 0, 0);
}
function openGlomZoom(v) {
  const src = v.src && HEW[v.src]; if (src) { Object.assign(HEW.gz, { on: src.on, imgs: new Set(src.imgs), a: src.a, outl: src.outl, cells: src.cells }); heSync("gz"); }   // starts with the H&E settings of the window it was opened from; the popup's own controls can change them
  S.gz = v; $("gztitle").textContent = `${v.label}${v.qc ? " · " + v.qc : ""} · ${dlabel(v.s.id)}`; $("glomzoom").hidden = false; drawGlomZoom();
  if (v.focus != null) setHighlight(v.s.id, [v.focus], v.cells.filter((c) => c !== v.focus)); else setHighlight(v.s.id, v.cells);
}
function gzZoomBind() {   // wheel = zoom at the cursor, drag = pan, double-click = reset (the saved PNG is what you see)
  const cv = $("gzcanvas"), G = 900, clamp = (q) => Math.max(0, Math.min(1, q)); let raf = 0, drag = null;
  const geom = () => { const v = S.gz, B = glomBox(v.cells, v.poly, 12, v.outline), bw = B.x1 - B.x0, bh = B.y1 - B.y0; return { bw, bh, k: Math.min((G - 8) / bw, (G - 8) / bh) }; };
  const pos = (e) => { const r = cv.getBoundingClientRect(); return { X: (e.clientX - r.left) * cv.width / r.width, Y: (e.clientY - r.top) * cv.height / r.height, sc: cv.width / r.width }; };
  const apply = (nv) => { S.gz.view = nv; if (!raf) raf = requestAnimationFrame(() => { raf = 0; drawGlomZoom(); }); };
  cv.addEventListener("wheel", (e) => {
    if (!S.gz) return; const { X, Y } = pos(e); if (Y > G) return; e.preventDefault();
    const v = S.gz.view || ID_VIEW, g = geom(), z2 = Math.max(1, Math.min(60, v.z * Math.exp(-e.deltaY * 0.0018))), wx = (X - G / 2) / (g.k * v.z) + v.rx * g.bw, wy = (Y - G / 2) / (g.k * v.z) + v.ry * g.bh;
    apply({ z: z2, rx: z2 === 1 ? 0.5 : clamp((wx - (X - G / 2) / (g.k * z2)) / g.bw), ry: z2 === 1 ? 0.5 : clamp((wy - (Y - G / 2) / (g.k * z2)) / g.bh) });
  }, { passive: false });
  cv.addEventListener("pointerdown", (e) => { if (!S.gz || pos(e).Y > G) return; S.gzMoved = false; drag = { x: e.clientX, y: e.clientY, v: S.gz.view || ID_VIEW }; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener("pointermove", (e) => { if (!drag) return; if (Math.hypot(e.clientX - drag.x, e.clientY - drag.y) > 3) S.gzMoved = true; if (!S.gzMoved) return;
    const g = geom(), sc = pos(e).sc, v = drag.v; apply({ z: v.z, rx: clamp(v.rx - ((e.clientX - drag.x) * sc) / (g.k * v.z * g.bw)), ry: clamp(v.ry - ((e.clientY - drag.y) * sc) / (g.k * v.z * g.bh)) }); });
  cv.addEventListener("pointerup", () => { drag = null; });
  cv.addEventListener("dblclick", () => { if (S.gz) apply(null); });
}
function closeGlomZoom() { $("glomzoom").hidden = true; S.gz = null; }
function saveGlomPNG() {
  if (!S.gz) return;
  $("gzcanvas").toBlob((blob) => {
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob);
    a.download = `${S.gz.qc ? "glom" : "neighborhood"}_${S.gz.label.replace(/[^\w-]+/g, "_")}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }, "image/png");
}
async function openGlomIn(sample, cells, poly) {
  closeAllWindows();
  if (sample !== S.sample) { $("sample").value = sample; await loadSample(sample); }
  const B = glomBox(cells, poly, 30), el = $("deckwrap"), w = el.clientWidth || 800, h = el.clientHeight || 600;
  S.viewState = { target: [(B.x0 + B.x1) / 2, -(B.y0 + B.y1) / 2, 0], zoom: Math.max(1.2, Math.log2(Math.min(w / (B.x1 - B.x0), h / (B.y1 - B.y0)))) };
  S.selected = new Set(cells);
  setView(S.viewState); updateAll();
}

// Main-view figure: the WebGL canvas plus a scale bar on the image and a footer with sample / filters / legend.
function saveViewPNG() {
  deckgl.redraw(true);
  const src = deckgl.getCanvas(), side = HEW.main.side && HEX.deck ? (HEX.deck.redraw(true), HEX.deck.getCanvas()) : null;   // side = the image pane beside the cells (a second deck)
  const W1 = src.width, GAP = side ? Math.max(4, Math.round(W1 / 200)) : 0, W = W1 + GAP + (side ? side.width : 0), H = Math.max(src.height, side ? side.height : 0), u = Math.max(1, W1 / 1100), info = S.colorInfo, mode = S.colorBy;
  let items = null;
  if (info.type === "cat") {
    const { out, miss } = countsFor(mode);
    items = out.map((n, i) => [i, n]).filter(([, n]) => n > 0).sort((a, b) => b[1] - a[1]).map(([i, n]) => [display(mode, i), catColor(mode, i), n]);
    if (miss) items.push(["missing", rgbStr(GRAY), miss]);
    if (S.n - S.visIdx.length) items.push(["filtered out", "rgba(105,110,122,0.6)", S.n - S.visIdx.length]);
  }
  const X = 24 * u, LW = W - 2 * X, scratch = document.createElement("canvas").getContext("2d"), lh = pngLegend(scratch, info, items, X, 0, LW, u, 9);
  const pitems = S.showPath && S.pathInfo?.legend ? S.pathInfo.legend.map((e) => [e.name, e.color, e.n]) : null;   // pathologist regions (outline colours)
  const ph = pitems ? 30 * u + pngLegend(scratch, { type: "cat" }, pitems, X, 0, LW, u, 3) : 0;
  const FH = Math.round(124 * u + lh + ph + 14 * u), out = document.createElement("canvas"); out.width = W; out.height = H + FH;
  const x = out.getContext("2d"); x.fillStyle = "#0f1218"; x.fillRect(0, 0, W, H + FH); x.drawImage(src, 0, 0);
  if (side) { x.drawImage(side, W1 + GAP, 0); x.fillStyle = "#fff"; x.fillRect(W1, 0, GAP, H); }   // white line between the panels
  // scale bar (output px per µm = CSS px per µm × canvas scale); one per pane, each with its own zoom
  const bar = (cv, zoom, x0, pw) => {
    const ppu = Math.pow(2, zoom) * (cv.width / (cv.clientWidth || cv.width)), target = 0.14 * pw / ppu, p = Math.pow(10, Math.floor(Math.log10(target)));
    const nice = [1, 2, 5, 10].map((m) => m * p).reduce((a, b) => (Math.abs(b - target) < Math.abs(a - target) ? b : a)), lab = nice >= 1000 ? nice / 1000 + " mm" : nice + " µm", bx = x0 + X;
    x.fillStyle = "rgba(15,18,24,0.6)"; x.fillRect(bx - 8 * u, H - 56 * u, nice * ppu + 16 * u, 48 * u);
    x.fillStyle = "#fff"; x.fillRect(bx, H - 24 * u, nice * ppu, 5 * u); x.font = `${15 * u}px ${FONT}`; x.textBaseline = "alphabetic"; x.fillText(lab, bx, H - 32 * u);
  };
  bar(src, S.viewState.zoom, 0, W1);
  if (side) bar(side, (HEX.vs || S.viewState).zoom, W1 + GAP, side.width);
  // footer
  const inf = S.meta.samples.find((s) => s.id === S.sample);
  x.fillStyle = "#d7dde6"; x.font = `600 ${22 * u}px ${FONT}`; x.fillText(dlabel(S.sample), X, H + 34 * u);
  x.font = `${14 * u}px ${FONT}`; x.fillStyle = "#8b95a5";
  const flt = [...Object.entries(S.filters).filter(([k, set]) => set.size && S.meta.cats[k]).map(([k, set]) => `${labelOf(k)} = ${[...set].map((i) => display(k, i)).join(", ")}`), ...exprSpecs().map(exprDesc)].join(" · ");
  const clip = (t) => (t.length > 220 ? t.slice(0, 217) + "…" : t);
  x.fillText(`${inf.disease_type} · ${subtypeOf(inf)} · ${fmt(S.n)} cells · colour: ${colorName()}${S.visIdx.length < S.n ? ` · ${fmt(S.visIdx.length)} highlighted` : ""}`, X, H + 60 * u);
  if (flt) x.fillText(clip("filters: " + flt), X, H + 82 * u);
  pngLegend(x, info, items, X, H + 112 * u, LW, u, 9);
  if (pitems) { x.font = `600 ${14 * u}px ${FONT}`; x.fillStyle = "#aab4c3"; x.fillText("Pathologist region outlines (number of regions)", X, H + 112 * u + lh + 22 * u); pngLegend(x, { type: "cat" }, pitems, X, H + 112 * u + lh + 46 * u, LW, u, 3); }
  out.toBlob((blob) => { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = `view_${S.sample}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, "image/png");
}

// ---- Save PNG while the Gallery / Glomeruli window is open: one picture of the window (visible part or everything)
function onSavePNG() {
  if (!$("glomzoom").hidden) return saveGlomPNG();
  if (!$("glomcmp").hidden) return saveGlomCmpPNG();
  if (!$("ghet").hidden) return ghSavePNG();
  if (!$("vio").hidden) { alert("PNG export is not available in the Violins window yet; take a screenshot."); return; }
  if (!$("sigwin").hidden) { alert("PNG export is not available in the Signatures window yet. Use \"Colour the main viewer by this signature\" and save from the main viewer, or take a screenshot."); return; }
  if (!$("nbr").hidden) { $("nbrinfo").textContent = "PNG export is not available for the Neighbors window yet — use Export CSV (or a screenshot)."; return; }
  if (!$("slide").hidden) return $("slidesel").value === "__all" ? saveSlideZip() : saveSlidePNG();
  if (!$("gallery").hidden) return showPngMenu("gallery");
  if (!$("rgn").hidden) return showPngMenu("rgn");
  if (!$("nbh").hidden) return showPngMenu("nbh");
  if (!$("glom").hidden) return showPngMenu("glom");
  saveViewPNG();
}
function saveSlidePNG() {
  const a1 = $("slidecv"), a2 = $("slidecv2");
  let src = a1;
  if (!a2.hidden) {   // the cells and the images beside them in one picture
    const gap = Math.max(4, Math.round(a1.width / 200)); src = document.createElement("canvas"); src.width = a1.width + gap + a2.width; src.height = Math.max(a1.height, a2.height);
    const x = src.getContext("2d"); x.fillStyle = "#fff"; x.fillRect(0, 0, src.width, src.height); x.drawImage(a1, 0, 0); x.drawImage(a2, a1.width + gap, 0);   // white line between the panels
  }
  src.toBlob((blob) => { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = `slide_${($("slidesel").value || "slide").replace(/[^\w-]+/g, "_")}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, "image/png"); }
function showPngMenu(kind) {
  const m = $("pngmenu"); m.dataset.kind = kind; m.hidden = false;
  $("pngall").textContent = kind === "gallery" && !$("galonly").checked ? "Everything: sharp numbered pages (.zip of PNG files)" : "Everything (all tiles, full scroll height)";
}
async function saveWindowPNG(kind, whole) {
  if (kind === "gallery" && whole && !$("galonly").checked) return saveGalleryZip();
  const W_ = { gallery: ["gallery", "galhead", "galgrid"], glom: ["glom", "glomhead", "glomgrid"], nbh: ["nbh", "nbhhead", "nbhgrid"], rgn: ["rgn", "rgnhead", "rgngrid"] }[kind];
  const win = $(W_[0]), head = $(W_[1]), grid = $(W_[2]), ctxS = S.ctxs[kind];
  if (!ctxS) return;
  const wr = win.getBoundingClientRect(), W = win.clientWidth, top = (el) => el.getBoundingClientRect().top - wr.top + win.scrollTop;
  const gridTop = top(grid) - 4, y0 = whole ? gridTop : Math.max(win.scrollTop + head.offsetHeight, gridTop), y1 = whole ? win.scrollHeight : win.scrollTop + win.clientHeight;
  // tiles / headers that fall inside the exported range (computed first so the legend counts only what is in the picture)
  const els = [...grid.querySelectorAll(".sechead, .tile, .gtile")].filter((e) => e.getBoundingClientRect().width > 0 && !e.closest("[hidden]"))
    .filter((e) => { const r = e.getBoundingClientRect(), ey = r.top - wr.top + win.scrollTop; return !(ey + r.height < y0 || ey > y1); });
  const tiles = els.filter((e) => !e.classList.contains("sechead")), allTiles = grid.querySelectorAll(".tile, .gtile").length;
  const info = ctxS.sc.info; let items = null;
  if (info.type === "cat") {
    const cnt = new Map(); let dimmed = 0;
    const add = (b, idx) => { for (const i of idx) { if (!b.pass[i]) { dimmed++; continue; } cnt.set(b.codes[i], (cnt.get(b.codes[i]) || 0) + 1); } };
    if (kind === "gallery") {
      const ids = new Set(tiles.map((e) => e.dataset.id)), bs = await Promise.all(ctxS.bp);
      bs.forEach((b, k) => { if (ids.has(ctxS.ids[k])) add(b, Array.from({ length: b.n }, (_, i) => i)); });
    } else { const inc = new Set(tiles); for (const g of ctxS.shown) if (inc.has(g.el)) add(g.b, g.cells); }
    items = [...cnt.entries()].sort((a, b) => b[1] - a[1]).map(([c, n]) => [c === MISSING ? "missing" : display(S.colorBy, c), c === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, c), n]);
    if (dimmed) items.push(["filtered out", "rgba(105,110,122,0.6)", dimmed]);
  }
  const clean = (t) => t.replace(/ · (click [a-z ]+ to open it|scroll for more|loading…)/g, "");
  const title = kind === "gallery" ? `Gallery — ${$("galsub").selectedOptions[0].text}` : kind === "rgn" ? `Region — ${$("rgnsel").selectedOptions[0].text} (${$("rgnsub").selectedOptions[0].text})` : kind === "nbh" ? `Neighborhoods — ${$("nbhsample").selectedOptions[0].text}`
    : `Glomeruli — ${$("glomscope").selectedOptions[0].text}${$("glomqc").value !== "all" ? " · " + $("glomqc").value : ""}`;
  const sub = clean($(kind === "gallery" ? "galinfo" : kind === "nbh" ? "nbhinfo" : kind === "rgn" ? "rgninfo" : "glominfo").textContent) + (whole ? "" : ` · showing ${tiles.length} of ${allTiles} tiles (the part of the window in view)`);
  const cssH0 = 100 + 40 + (y1 - y0), U = Math.max(0.5, Math.min(2, 15500 / cssH0));  // keep the canvas under browser limits
  const X = 14 * U, scratch = document.createElement("canvas").getContext("2d");
  const lh = pngLegend(scratch, info, items, X, 0, W * U - 2 * X, U, 6);
  const HH = Math.round(84 * U + lh + 10 * U), contentH = Math.round((y1 - y0) * U);
  const out = document.createElement("canvas"); out.width = Math.round(W * U); out.height = HH + contentH;
  const x = out.getContext("2d"); x.fillStyle = "#0f1218"; x.fillRect(0, 0, out.width, out.height);
  x.textBaseline = "alphabetic"; x.fillStyle = "#e3e8ef"; x.font = `600 ${22 * U}px ${FONT}`; x.fillText(title, X, 32 * U);
  x.font = `${13 * U}px ${FONT}`; x.fillStyle = "#8b95a5";
  const subClip = sub.length > 200 ? sub.slice(0, 197) + "…" : sub; x.fillText(subClip, X, 54 * U);
  pngLegend(x, info, items, X, 82 * U, W * U - 2 * X, U, 6);
  // tiles + section headers, positioned exactly as laid out in the window
  x.save(); x.beginPath(); x.rect(0, HH, out.width, contentH); x.clip();
  const toX = (cx) => (cx - wr.left) * U, toY = (cy) => HH + (cy - wr.top + win.scrollTop - y0) * U;
  for (const e of els) {
    const r = e.getBoundingClientRect();
    if (e.classList.contains("sechead")) {
      x.strokeStyle = "#29303b"; x.lineWidth = U; x.beginPath(); x.moveTo(toX(r.left), toY(r.bottom)); x.lineTo(toX(r.right), toY(r.bottom)); x.stroke();
    } else {
      const bx = toX(r.left), by = toY(r.top), bw = r.width * U, bh = r.height * U, rad = 6 * U;
      x.beginPath(); x.moveTo(bx + rad, by); x.arcTo(bx + bw, by, bx + bw, by + bh, rad); x.arcTo(bx + bw, by + bh, bx, by + bh, rad); x.arcTo(bx, by + bh, bx, by, rad); x.arcTo(bx, by, bx + bw, by, rad); x.closePath();
      x.fillStyle = "#141922"; x.fill(); x.lineWidth = 2 * U; x.strokeStyle = e.classList.contains("cur") ? "#2a6fdb" : "#29303b"; x.stroke();
      const cv = e.querySelector("canvas"); if (cv) { const c = cv.getBoundingClientRect(); x.drawImage(cv, toX(c.left), toY(c.top), c.width * U, c.height * U); }
    }
    const tw = document.createTreeWalker(e, NodeFilter.SHOW_TEXT);
    for (let n = tw.nextNode(); n; n = tw.nextNode()) {
      const t = n.textContent.trim(); if (!t) continue;
      const par = n.parentElement, cs = getComputedStyle(par), rg = document.createRange(); rg.selectNodeContents(n);
      const q = rg.getClientRects()[0]; if (!q) continue;
      if (par.classList.contains("chip")) { const c = par.getBoundingClientRect(); x.fillStyle = cs.backgroundColor; x.fillRect(toX(c.left), toY(c.top), c.width * U, c.height * U); }
      x.font = `${cs.fontWeight} ${parseFloat(cs.fontSize) * U}px ${FONT}`; x.fillStyle = cs.color; x.textBaseline = "middle"; x.fillText(t, toX(q.left), toY(q.top + q.height / 2));
    }
  }
  x.restore();
  const base = kind === "rgn" ? `region_${$("rgnsel").selectedOptions[0].text.replace(/[^\w-]+/g, "_")}_${$("rgnsub").value === "__all" ? "all_samples" : $("rgnsub").value.replace(/[^\w-]+/g, "_")}` : kind === "gallery" ? `gallery_${$("galsub").value === "__all" ? "all_samples" : $("galsub").value.replace(/[^\w-]+/g, "_")}` : kind === "nbh" ? `neighborhoods_${$("nbhsample").value.replace(/[^\w-]+/g, "_")}` : `glomeruli_${$("glomscope").value.replace(/[^\w-]+/g, "_")}`;
  out.toBlob((blob) => { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = `${base}_${whole ? "all" : "visible"}.png`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, "image/png");
}


// ---- Gallery "Everything": every tile is drawn again at high resolution (not the small on-screen bitmap enlarged) onto numbered pages; the pages are
// delivered as ONE .zip of PNG files. The zip is written here without a library (files stored uncompressed: PNG is compressed already).
const CRC_T = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
const crc32 = (u8) => { let c = 0xffffffff; for (let i = 0; i < u8.length; i++) c = CRC_T[(c ^ u8[i]) & 255] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; };
function zipBlob(files) {   // files: [{ name, data: Uint8Array }]
  const enc = new TextEncoder(), parts = [], cen = [], now = new Date();
  const dt = ((now.getFullYear() - 1980) << 9) | ((now.getMonth() + 1) << 5) | now.getDate(), tm = (now.getHours() << 11) | (now.getMinutes() << 5) | (now.getSeconds() >> 1);
  let off = 0, csz = 0;
  for (const f of files) {
    const nm = enc.encode(f.name), crc = crc32(f.data), sz = f.data.length;
    const lh = new DataView(new ArrayBuffer(30));
    lh.setUint32(0, 0x04034b50, true); lh.setUint16(4, 20, true); lh.setUint16(6, 0x0800, true); lh.setUint16(10, tm, true); lh.setUint16(12, dt, true);
    lh.setUint32(14, crc, true); lh.setUint32(18, sz, true); lh.setUint32(22, sz, true); lh.setUint16(26, nm.length, true);
    parts.push(lh.buffer, nm, f.data);
    const ch = new DataView(new ArrayBuffer(46));
    ch.setUint32(0, 0x02014b50, true); ch.setUint16(4, 20, true); ch.setUint16(6, 20, true); ch.setUint16(8, 0x0800, true); ch.setUint16(12, tm, true); ch.setUint16(14, dt, true);
    ch.setUint32(16, crc, true); ch.setUint32(20, sz, true); ch.setUint32(24, sz, true); ch.setUint16(28, nm.length, true); ch.setUint32(42, off, true);
    cen.push(ch.buffer, nm); off += 30 + nm.length + sz; csz += 46 + nm.length;
  }
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true); end.setUint16(8, files.length, true); end.setUint16(10, files.length, true); end.setUint32(12, csz, true); end.setUint32(16, off, true);
  return new Blob([...parts, ...cen, end.buffer], { type: "application/zip" });
}
async function saveGalleryZip() {
  const ctxS = S.ctxs.gallery; if (!ctxS) return;
  const msg = $("galinfo"), keep = msg.textContent, info = S.grey ? { type: "none" } : ctxS.sc.info;   // Grey cells: no colour legend
  try {
    const PW = 5000, PH = 3500, M = 60, COLS = 5, GAP = 24, U = 2.2, box = +$("galsz").value || GAL_BOX;
    const CW = Math.floor((PW - 2 * M - GAP * (COLS - 1)) / COLS), IMH = Math.round(CW * 0.78), CAP = 110, SEC = 84, ROW = IMH + CAP + GAP, Q = CW / box;
    const bs = await Promise.all(ctxS.bp), byId = new Map(ctxS.ids.map((id, k) => [id, k])), smp = new Map(S.meta.samples.map((s) => [s.id, s]));
    let items = null;
    if (info.type === "cat") {
      const cnt = new Map(); let dimmed = 0;
      for (const b of bs) for (let i = 0; i < b.n; i++) { if (!b.pass[i]) { dimmed++; continue; } cnt.set(b.codes[i], (cnt.get(b.codes[i]) || 0) + 1); }
      items = [...cnt.entries()].sort((a, b) => b[1] - a[1]).map(([c, n]) => [c === MISSING ? "missing" : display(S.colorBy, c), c === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, c), n]);
      if (dimmed) items.push(["filtered out", "rgba(105,110,122,0.6)", dimmed]);
    }
    const clean = (t) => t.replace(/ · (click [a-z ]+ to open it|scroll for more|loading…)/g, "");
    const title = `Gallery — ${$("galsub").selectedOptions[0].text}`, sub = clean(keep) + (S.grey ? " · cells grey" : "");
    const scratch = document.createElement("canvas").getContext("2d"), LW = PW - 2 * M, lh = pngLegend(scratch, info, items, M, 0, LW, U, 6), HEAD = Math.round(84 * U + lh + 10 * U);
    // flow: section headers and tiles onto pages
    const flat = []; for (const sec of $("galgrid").querySelectorAll(".galsec")) { flat.push({ head: sec.querySelector(".sechead").textContent.trim() }); for (const t of sec.querySelectorAll(".tile")) flat.push({ id: t.dataset.id }); }
    const pages = []; let page = null, y = 0, col = 0;
    const fresh = () => { page = []; pages.push(page); y = HEAD; col = 0; };
    fresh();
    for (const it of flat) {
      if (it.head) { if (col) { y += ROW; col = 0; } if (y + SEC + ROW > PH - M && page.length) fresh(); page.push({ head: it.head, y }); y += SEC; }
      else { if (col === COLS) { y += ROW; col = 0; } if (y + ROW > PH - M) fresh(); page.push({ id: it.id, x: M + col * (CW + GAP), y }); col++; }
    }
    const wantPoly = HEW.gal.outl, files = [], base = `gallery_${$("galsub").value === "__all" ? "all_samples" : $("galsub").value.replace(/[^\w-]+/g, "_")}`;
    for (let pi = 0; pi < pages.length; pi++) {
      msg.textContent = `exporting page ${pi + 1} of ${pages.length}…`; await new Promise((r) => setTimeout(r, 0));
      const out = document.createElement("canvas"); out.width = PW; out.height = PH;
      const x = out.getContext("2d"); x.fillStyle = "#0f1218"; x.fillRect(0, 0, PW, PH); x.textBaseline = "alphabetic";
      x.fillStyle = "#e3e8ef"; x.font = `600 ${22 * U}px ${FONT}`; x.fillText(title, M, 32 * U);
      x.font = `${13 * U}px ${FONT}`; x.fillStyle = "#8b95a5"; x.fillText(sub.length > 260 ? sub.slice(0, 257) + "…" : sub, M, 54 * U);
      pngLegend(x, info, items, M, 82 * U, LW, U, 6);
      x.font = `${12 * U}px ${FONT}`; x.fillStyle = "#8b95a5"; x.textAlign = "right"; x.fillText(`page ${pi + 1} of ${pages.length}`, PW - M, PH - 20); x.textAlign = "left";
      const fit = (t, w) => { if (x.measureText(t).width <= w) return t; let u = t; while (u.length > 3 && x.measureText("…" + u).width > w) u = u.slice(1); return "…" + u; };
      for (const e of pages[pi]) {
        if (e.head) {
          x.font = `600 ${15 * U}px ${FONT}`; x.fillStyle = "#e3e8ef"; x.textBaseline = "middle"; x.fillText(e.head, M, e.y + SEC / 2 - 6); x.textBaseline = "alphabetic";
          x.strokeStyle = "#29303b"; x.lineWidth = 2; x.beginPath(); x.moveTo(M, e.y + SEC - 14); x.lineTo(PW - M, e.y + SEC - 14); x.stroke(); continue;
        }
        const k = byId.get(e.id), s = smp.get(e.id), b = bs[k], poly = wantPoly ? await galPoly(e.id).catch(() => null) : null;
        const tc = document.createElement("canvas"); drawTile(tc, s, b, ctxS.sc, box, null, null, { W: CW, H: IMH }, poly, null, "gal", Q);
        x.drawImage(tc, e.x, e.y); x.strokeStyle = "#29303b"; x.lineWidth = 2; x.strokeRect(e.x, e.y, CW, IMH); tc.width = 0;
        x.fillStyle = "#e3e8ef"; x.font = `600 ${11 * U}px ${FONT}`; x.fillText(fit(dlabel(e.id), CW), e.x, e.y + IMH + 15 * U);
        x.fillStyle = "#8b95a5"; x.font = `${9.5 * U}px ${FONT}`; x.fillText(fit(`${fmt(s.n)} cells · ${s.disease_type} · ${subtypeOf(s)}`, CW), e.x, e.y + IMH + 29 * U);
      }
      const blob = await new Promise((r) => out.toBlob(r, "image/png"));
      files.push({ name: `${base}_p${String(pi + 1).padStart(2, "0")}.png`, data: new Uint8Array(await blob.arrayBuffer()) }); out.width = 0;
    }
    const a = document.createElement("a"); a.href = URL.createObjectURL(zipBlob(files)); a.download = `${base}_pages.zip`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 4000);
  } catch (err) { console.error(err); alert("Export failed: " + err.message); }
  finally { msg.textContent = keep; }
}

// ---------------------------------------------------------------- signatures window
// Curated gene signatures (server: signatures.json) -> explorer (genes, change across disease subtypes, smoothed maps) and type I vs II localisation.
const SIG = { list: null, byId: {}, sel: "ifn1", sig2: "", group: "subtype", tab: "explore", sigma: 25, scope: "all", useFilters: false, view: "plot", a: "ifn1", b: "ifn2", locSigma: 25, locQ: 10, locSub: "", mapSample: null, tok: 0 };
const SIG_GROUPS = ["Interferon", "Fibrosis", "Immune & inflammation", "Kidney programmes"];
const mk = (tag, cls, txt) => { const e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; };
const svgEl = (tag, attrs = {}, parent) => { const e = document.createElementNS("http://www.w3.org/2000/svg", tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (parent) parent.append(e); return e; };
const fmtP = (p) => (p == null ? "n/a" : p < 0.001 ? "<0.001" : p.toFixed(3));
const pTxt = (p) => (p == null ? "p n/a" : p < 0.001 ? "p < 0.001" : "p = " + p.toFixed(3));
const fmtN = (v, d = 2) => (v == null || Number.isNaN(v) ? "–" : (v > 0 ? "+" : "") + v.toFixed(d));
const sigLabel = () => (S.sigName ? S.sigName : `Signature (${S.sigGenes.length} genes)`);
const sigLabel2 = () => (S.sigName2 ? S.sigName2 : `Signature B (${S.sigGenes2.length} genes)`);
async function postJSON(url, body) { const r = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) }); if (!r.ok) throw new Error((await r.text()).slice(0, 200)); return r.json(); }
async function loadSigList() {
  if (SIG.list) return; const d = await jget("/api/signatures"); SIG.list = d.signatures; SIG.order = d.subtype_order; SIG.byId = Object.fromEntries(d.signatures.map((x) => [x.id, x]));
  for (const sel of [$("sigpreset"), $("sigpreset2"), $("xfpreset")]) if (sel && sel.options.length < 2) for (const g of SIG_GROUPS) { const og = document.createElement("optgroup"); og.label = g;
    d.signatures.filter((x) => x.group === g).forEach((x) => { const o = document.createElement("option"); o.value = x.id; o.textContent = `${x.name} (${x.genes.length})`; og.append(o); }); sel.append(og); }
}
async function openSig() {
  closeAllWindows();
  try { await loadSigList(); } catch (e) { alert("Signatures not available: " + e.message); return; }
  if (!SIG.list.length) { alert("This project has no signature definitions (signatures.json), so the Signatures window has nothing to show. You can still filter by any gene or gene list under Filters → Gene / signature positive."); return; }
  $("sigwin").hidden = false; renderSigList(); renderSigPane();
}
function closeSig() { $("sigwin").hidden = true; }
function renderSigList() {
  const q = $("sigsearch").value.toLowerCase(), box = $("sigitems"); box.replaceChildren();
  for (const g of SIG_GROUPS) {
    const items = SIG.list.filter((x) => x.group === g && (x.name + " " + x.desc).toLowerCase().includes(q)); if (!items.length) continue;
    box.append(mk("div", "sgrp", g));
    for (const x of items) {
      const it = mk("div", "sitem" + (x.id === SIG.sel ? " sel" : "")), nm = mk("span", "nm", x.name); nm.title = x.desc;
      const weak = x.coherence_bin != null && x.coherence_bin < 0.10, ct = mk("span", "ct" + (weak ? " weak" : ""), `${x.genes.length} genes${weak ? " · weak" : ""}`);
      if (weak) ct.title = "genes barely co-vary at tissue scale; interpret with caution";
      it.append(nm, ct); it.addEventListener("click", () => { SIG.sel = x.id; renderSigList(); if (SIG.tab === "explore") renderSigPane(); }); box.append(it);
    }
  }
}
function renderSigPane() { $("sigtab-explore").classList.toggle("on", SIG.tab === "explore"); $("sigtab-loc").classList.toggle("on", SIG.tab === "loc"); $("siglist").style.display = SIG.tab === "explore" ? "" : "none"; SIG.tab === "explore" ? renderExplorer() : renderLoc(); }
const subtypeRank = (t) => { const o = SIG.order || []; return o.includes(t) ? o.indexOf(t) : 99; };
// ---- colour scale for score maps: diverging around 0 for z-scores, viridis otherwise; values come from all shown samples
function scoreScale(vecs, z) {
  const pool = []; for (const v of vecs) for (let i = 0; i < v.length; i += 3) pool.push(v[i]);
  const hi = quantile(pool, 0.995), lo = quantile(pool, z ? 0.005 : 0.01);
  if (z) { const m = Z_LIM; return { z, lo: -m, hi: m, color: (v) => DIVL[Math.round(Math.max(0, Math.min(1, (v + m) / (2 * m))) * 255)] }; }
  const sp = Math.max(hi - lo, 1e-9); return { z, lo, hi, color: (v) => LUT[Math.round(Math.max(0, Math.min(1, (v - lo) / sp)) * 255)] };
}
function drawScoreTile(cv, meta, xy, vec, sc, box = 260) {
  const w = meta.x[1] - meta.x[0] || 1, h = meta.y[1] - meta.y[0] || 1, k = Math.min((box - 8) / w, (box - 8) / h);
  cv.width = Math.max(40, Math.ceil(w * k) + 8); cv.height = Math.max(40, Math.ceil(h * k) + 8);
  const ctx = cv.getContext("2d"), n = vec.length, px = n > 30000 ? 1.5 : 2.2; ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, cv.width, cv.height);
  cv._redraw = () => drawScoreTile(cv, meta, xy, vec, sc, box);
  const hs = heFor(meta.id, "sig"); if (hs) { drawHE(ctx, { sample: hs, redo: () => heRedoCv(cv), w: "sig" }, (x) => 4 + (x - meta.x[0]) * k, (y) => 4 + (y - meta.y[0]) * k, cv.width, cv.height, k); ctx.globalAlpha = HEW.sig.a; }
  const o = new Uint32Array(n); for (let i = 0; i < n; i++) o[i] = i; o.sort((p, q) => vec[p] - vec[q]);   // low values first, hot cells on top
  if (HEW.sig.cells) for (const i of o) { const c = sc.color(vec[i]); ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`; ctx.fillRect(4 + (xy[2 * i] - meta.x[0]) * k, 4 + (xy[2 * i + 1] - meta.y[0]) * k, px, px); }
  ctx.globalAlpha = 1;
}
const scoreLegend = (sc, label) => { const w = mk("div"); const bar = mk("div", "legbar"); bar.style.background = `linear-gradient(90deg,${(sc.z ? DIVL : LUT).filter((_, i) => i % 32 === 0 || i === 255).map((c) => `rgb(${c})`).join(",")})`;
  const tk = mk("div"); tk.innerHTML = tickRowHTML(scaleTicks({ type: "num", lo: sc.lo, hi: sc.hi }));
  const t = mk("div", "snote", sc.z ? "z-score; 0 = average" : "smoothed log-normalised expression"); w.append(mk("div", "snote", label), bar, tk, t); return w; };
const exprUrl2 = (sid, genes, sigma, scope) => `/api/sample/${encodeURIComponent(sid)}/expr?genes=${encodeURIComponent(genes.join(","))}&layer=logcpx&sigma=${sigma}&scope=${scope}`;
const sigVecCache = {};
const sigVec = (sid, genes, sigma, scope) => { const k = [sid, sigma, scope, genes.join(",")].join("|"); return sigVecCache[k] || (sigVecCache[k] = bget(exprUrl2(sid, genes, sigma, scope), Float32Array)); };

// ---------------------------------------------------------------- explorer
const SIG_TYPES = ["Control", "No LN", "LN"], SIG_A = "#3987e5", SIG_B = "#d95926";   // series colours: validated slots 1 and 2 (blue / orange: colour-blind-safe pair)
const SCOPE_TXT = { none: "none (smoothed raw score)", sample: "within each donor", subtype: "across samples of the same disease subtype", type: "across samples of the same disease type", all: "across all samples", nocontrol: "against the non-control samples (Control excluded from the reference)" };
const groupLabel = (m) => (SIG.group === "type" ? m.disease_type || "unknown" : subtypeOf(m));
const groupRank = (g) => { const o = SIG.group === "type" ? SIG_TYPES : SIG.order || []; return o.includes(g) ? o.indexOf(g) : 99; };
function scopeWarning() {   // z-scoring inside a group removes the differences between the groups it pools over
  const sc = SIG.scope, gb = SIG.group;
  if (sc === "sample") return "Each donor is z-scored on its own (its tissue pieces pooled), so every donor's average is 0: differences between donors and groups are removed by design and only the spatial pattern inside a donor remains. The group tests below are not informative.";
  if (sc === "subtype") return gb === "subtype" ? "Samples are z-scored within their disease subtype, so each subtype's average is about 0: differences between subtypes are removed by design (the tests below are not informative); differences between samples of one subtype remain." : "Samples are z-scored within their disease subtype, which also removes most differences between disease types; spread between samples of one subtype remains.";
  if (sc === "type") return gb === "type" ? "Samples are z-scored within their disease type (Control / No LN / LN), so the three types are each centred on 0: differences between types are removed by design. The three LN subtypes can still be compared with each other under 'group by subtype'." : "Samples are z-scored within their disease type, so Control and No LN (one subtype each) are centred on 0, while the three LN subtypes (Class I/II, III, V) remain comparable with each other.";
  return "";
}
const scopeNote = () => (SIG.scope === "nocontrol" ? "The mean and SD for the z-score come from the non-control samples only (No LN and LN, all their cells); every sample, Control included, is then scored against that reference. Differences between groups are kept, and Control samples show up as below the disease average." : "");
async function renderExplorer() {
  const pane = $("sigpane"), sg = SIG.byId[SIG.sel], tok = ++SIG.tok; if (SIG.sig2 === SIG.sel) SIG.sig2 = ""; const sg2 = SIG.sig2 ? SIG.byId[SIG.sig2] : null; pane.replaceChildren();
  const info = (g, color, tag) => { const box = mk("div"); const h = mk("div", "sh2", (tag ? tag + " " : "") + g.name); h.style.marginTop = tag ? "14px" : "0"; if (color) h.style.borderBottomColor = color; box.append(h, mk("div", "sdesc", g.desc)); if (g.note) box.append(mk("div", "snote", "Note: " + g.note));
    const cov = mk("div", "snote", `${g.genes.length} of ${g.n_defined} candidate genes are on the panel. ` + (g.coherence_bin != null ? `Coherence: ${g.coherence_bin.toFixed(2)} across 100 µm tissue bins, ${g.coherence_cell.toFixed(2)} per cell after removing cell-type means. ` : "")
      + (g.coherence_bin != null && g.coherence_bin < 0.10 ? "The genes barely co-vary at tissue scale, so treat this score with caution. " : "") + "Dashed genes co-vary weakly with the rest (r < 0.05)."); cov.style.marginTop = "6px"; box.append(cov);
    const chips = mk("div", "chips"); chips.style.marginTop = "6px"; for (const gn of g.genes) { const r = g.gene_r_bin?.[gn], c = mk("span", "chip2" + (r != null && r < 0.05 ? " weak" : ""), gn); c.title = `mean correlation with the other genes (100 µm bins): ${r == null ? "n/a" : r.toFixed(2)}`; chips.append(c); } box.append(chips);
    if (g.missing.length) { const d = mk("details", "snote"); d.style.marginTop = "6px"; d.append(mk("summary", "", `${g.missing.length} candidate genes are not on the panel`)); const m = mk("div", "chips"); g.missing.forEach((x) => m.append(mk("span", "chip2 miss", x))); d.append(m); box.append(d); } return box; };
  pane.append(info(sg, sg2 ? SIG_A : null, sg2 ? "A ·" : "")); if (sg2) pane.append(info(sg2, SIG_B, "B ·"));
  // controls
  const ctl = mk("div", "sctl"), sl = mk("input"); sl.type = "range"; sl.min = 0; sl.max = 100; sl.step = 5; sl.value = SIG.sigma; const sv = mk("span", "", SIG.sigma ? SIG.sigma + " µm" : "off"); const l1 = mk("label"); l1.append("Smoothing σ ", sl, sv);
  const sc = mk("select"); sc.innerHTML = Object.entries(SCOPE_TXT).map(([k, t]) => `<option value="${k}">z-score ${k === "none" ? t : t}</option>`).join("").replace("z-score none (smoothed raw score)", "no z-score (smoothed raw score)"); sc.value = SIG.scope; const l2 = mk("label"); l2.append("Score ", sc);
  const gb = mk("select"); gb.innerHTML = `<option value="subtype">disease subtype (5 groups)</option><option value="type">disease type (Control / No LN / LN)</option>`; gb.value = SIG.group; const l4 = mk("label"); l4.append("Group by ", gb);
  const s2 = mk("select"); s2.innerHTML = `<option value="">— none —</option>` + SIG_GROUPS.map((g) => `<optgroup label="${g}">${SIG.list.filter((x) => x.group === g && x.id !== SIG.sel).map((x) => `<option value="${x.id}">${x.name}</option>`).join("")}</optgroup>`).join(""); s2.value = SIG.sig2 || ""; const l5 = mk("label"); l5.append("Add a second signature ", s2);
  const fl = mk("input"); fl.type = "checkbox"; fl.checked = SIG.useFilters; const l3 = mk("label"); l3.append(fl, "only cells passing the current filters");
  const tv = mk("button", "", SIG.view === "plot" ? "Table" : "Plot"), tovw = mk("button", "", sg2 ? "Show both signatures in the main viewer" : "Colour the main viewer by this signature");
  ctl.append(l1, l2, l4, l5, l3, tv, tovw); ctl.classList.add("stickyctl"); pane.prepend(ctl);   // controls stay at the top of the pane
  const warn = scopeWarning(), note = scopeNote();
  if (warn) { const w = mk("div", "snote", "⚠ " + warn); w.style.cssText = "color:#e0a458;max-width:900px;margin:4px 0 8px"; ctl.after(w); }
  else if (note) { const w = mk("div", "snote", "ℹ " + note); w.style.cssText = "max-width:900px;margin:4px 0 8px"; ctl.after(w); }
  sl.addEventListener("input", () => { sv.textContent = +sl.value ? sl.value + " µm" : "off"; }); sl.addEventListener("change", () => { SIG.sigma = +sl.value; renderExplorer(); });
  sc.addEventListener("change", () => { SIG.scope = sc.value; renderExplorer(); }); gb.addEventListener("change", () => { SIG.group = gb.value; renderExplorer(); }); s2.addEventListener("change", () => { SIG.sig2 = s2.value; renderExplorer(); });
  fl.addEventListener("change", () => { SIG.useFilters = fl.checked; renderExplorer(); }); tv.addEventListener("click", () => { SIG.view = SIG.view === "plot" ? "table" : "plot"; renderExplorer(); });
  tovw.addEventListener("click", () => { S.sigGenes = [...sg.genes]; S.sigName = sg.name; S.sigKey = ""; $("sigtext").value = sg.genes.join(" "); S.smooth = SIG.sigma; $("smooth").value = SIG.sigma; $("smoothlab").textContent = SIG.sigma ? SIG.sigma + " µm" : "off";
    S.zscope = SIG.scope; $("zscope").value = S.zscope; S.colorBy = "__sig";
    if (sg2) { S.sigGenes2 = [...sg2.genes]; S.sigName2 = sg2.name; S.sigKey2 = ""; $("sigtext2").value = sg2.genes.join(" "); S.colorBy = "__sigdual"; } $("colorby").value = S.colorBy; syncGeneUI(); closeSig(); refresh(); });
  const out = mk("div", "snote", "Computing…"); pane.append(out);
  const body = (id) => ({ sig: id, sigma: SIG.sigma, scope: SIG.scope, group_by: SIG.group, ...(SIG.useFilters ? nbrIndexSpec() : { filters: {}, hidden: {}, expr: [] }) });
  let dA, dB = null; try { [dA, dB] = await Promise.all([postJSON("/api/signature/summary", body(sg.id)), sg2 ? postJSON("/api/signature/summary", body(sg2.id)) : null]); } catch (e) { out.textContent = "Error: " + e.message; return; }
  if (tok !== SIG.tok) return; out.remove();
  const z = SIG.scope !== "none", gname = SIG.group === "type" ? "disease type" : "disease subtype";
  pane.append(mk("div", "sh2", `Change across ${gname} (one value per donor: mean score of its cells, tissue pieces pooled)`));
  pane.append(SIG.view === "plot" ? stripPlot(dA, dB, sg, sg2, z) : sampleTable(dA, dB, sg, sg2)); pane.append(testsBox(dA, dB, sg, sg2));
  pane.append(mk("div", "sh2", "Which genes drive it (mean expression per group, each gene scaled across groups)")); pane.append(geneHeat(dA, sg, sg2 ? "A · " : "")); if (sg2) pane.append(geneHeat(dB, sg2, "B · "));
  pane.append(mk("div", "sh2", sg2 ? "Both signatures in space (blue = A, orange = B, white = both high), by sample and group" : "Smoothed signature in space, by sample and group")); const maps = mk("div", "smaps"); pane.append(maps); await drawSigMaps(maps, sg, sg2, dA, dB, z, tok);
}
function stripPlot(dA, dB, sg, sg2, z) {
  const ser = [{ d: dA, name: sg.name, col: SIG_A }].concat(dB ? [{ d: dB, name: sg2.name, col: SIG_B }] : []), groups = dA.groups.map((g) => g.subtype);
  const rowsOf = (d) => d.samples.filter((s) => s.n >= d.tests.min_cells && s.mean != null), all = ser.flatMap((s) => rowsOf(s.d).map((r) => r.mean));
  const W = 700, H = 320, m = { l: 58, r: 18, t: dB ? 34 : 14, b: 54 }, svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` });
  let lo = Math.min(...all), hi = Math.max(...all); const pad = (hi - lo) * 0.12 || 0.5; lo -= pad; hi += pad;
  const ys = (v) => m.t + (1 - (v - lo) / (hi - lo)) * (H - m.t - m.b), bw = (W - m.l - m.r) / groups.length, xc = (i) => m.l + bw * (i + 0.5);
  const step = (() => { const r = (hi - lo) / 5, p = Math.pow(10, Math.floor(Math.log10(r))); return [1, 2, 5, 10].map((k) => k * p).find((k) => k >= r) || r; })();
  for (let v = Math.ceil(lo / step) * step; v <= hi; v += step) { svgEl("line", { x1: m.l, x2: W - m.r, y1: ys(v), y2: ys(v), stroke: Math.abs(v) < 1e-9 && z ? "#4b566a" : "#29303b", "stroke-width": 1 }, svg); svgEl("text", { x: m.l - 8, y: ys(v) + 4, "text-anchor": "end", fill: "#8b95a5", "font-size": 11 }, svg).textContent = v.toFixed(2); }
  svgEl("text", { x: 14, y: (H - m.b + m.t) / 2, fill: "#8b95a5", "font-size": 11, transform: `rotate(-90 14 ${(H - m.b + m.t) / 2})`, "text-anchor": "middle" }, svg).textContent = z ? `score (z-score ${SCOPE_TXT[SIG.scope].replace("across samples of the same ", "within ")})`.slice(0, 44) : "smoothed score";
  if (dB) ser.forEach((s, k) => { const x0 = m.l + k * 250; svgEl("circle", { cx: x0 + 6, cy: 14, r: 5, fill: s.col, stroke: "#0f1218", "stroke-width": 2 }, svg); svgEl("text", { x: x0 + 18, y: 18, fill: "#d7dde6", "font-size": 12 }, svg).textContent = s.name.length > 34 ? s.name.slice(0, 33) + "…" : s.name; });
  groups.forEach((g, i) => {
    svgEl("text", { x: xc(i), y: H - m.b + 18, "text-anchor": "middle", fill: "#d7dde6", "font-size": 12 }, svg).textContent = g;
    svgEl("text", { x: xc(i), y: H - m.b + 33, "text-anchor": "middle", fill: "#8b95a5", "font-size": 11 }, svg).textContent = donTxt(nDonors(rowsOf(dA).filter((r) => r.group === g).map((r) => r.sample)));
    ser.forEach((s, k) => {
      const off = ser.length === 1 ? 0 : (k - 0.5) * bw * 0.34, gs = rowsOf(s.d).filter((r) => r.group === g), gm = s.d.groups[i], cx = xc(i) + off, half = ser.length === 1 ? 24 : 14;
      if (gm.se != null) svgEl("line", { x1: cx, x2: cx, y1: ys(gm.mean - gm.se), y2: ys(gm.mean + gm.se), stroke: "#e3e8ef", "stroke-width": 2 }, svg); svgEl("line", { x1: cx - half, x2: cx + half, y1: ys(gm.mean), y2: ys(gm.mean), stroke: "#e3e8ef", "stroke-width": 3 }, svg);
      gs.forEach((r, j) => { const jit = (((j * 0.618) % 1) - 0.5) * bw * (ser.length === 1 ? 0.55 : 0.26); svgEl("circle", { cx: cx + jit, cy: ys(r.mean), r: 5, fill: s.col, stroke: "#0f1218", "stroke-width": 2 }, svg);
        const hit = svgEl("circle", { cx: cx + jit, cy: ys(r.mean), r: 11, fill: "transparent" }, svg);
        hit.addEventListener("pointermove", (ev) => nbrTip(ev, [{ tag: "b", t: r.mean.toFixed(3) }, { t: `${dlabel(r.sample)} · ${s.name}` }, { cls: "m", t: `${r.subtype} · ${fmt(r.n)} cells · median ${r.median.toFixed(2)}` }])); hit.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; }); });
    });
  });
  const wrap = mk("div"); wrap.append(svg, mk("div", "snote", "Dots = donors (the tissue pieces of a donor are pooled); bar = mean ± standard error over the donors.")); return wrap;
}
function sampleTable(dA, dB, sg, sg2) {
  const tb = mk("table", "st"), hr = mk("tr"); ["Group / sample", "Cells", `Mean${dB ? " A" : ""}`, ...(dB ? ["Mean B"] : []), "Median" + (dB ? " A" : ""), "IQR" + (dB ? " A" : "")].forEach((t) => hr.append(mk("th", "", t))); tb.append(hr);
  const bById = dB ? Object.fromEntries(dB.samples.map((x) => [x.sample, x])) : {};
  dA.groups.forEach((g, gi) => {
    const tr = mk("tr", "sub"); [g.subtype + ` (${donTxt(g.n_donors)})`, "", g.mean.toFixed(3), ...(dB ? [dB.groups[gi].mean.toFixed(3)] : []), "", g.se != null ? "SE " + g.se.toFixed(3) : ""].forEach((t) => tr.append(mk("td", "", t))); tb.append(tr);
    for (const r of dA.samples.filter((x) => x.group === g.subtype && x.mean != null).sort((a, b) => dlabel(a.sample).localeCompare(dlabel(b.sample)))) { const t2 = mk("tr"), rb = bById[r.sample];
      ["   " + dlabel(r.sample) + (r.n < dA.tests.min_cells ? " (too few cells, not tested)" : ""), fmt(r.n), r.mean.toFixed(3), ...(dB ? [rb?.mean != null ? rb.mean.toFixed(3) : "–"] : []), r.median.toFixed(3), `${r.q25.toFixed(2)} – ${r.q75.toFixed(2)}`].forEach((t) => t2.append(mk("td", "", t))); tb.append(t2); }
  }); return tb;
}
function testsBox(dA, dB, sg, sg2) {
  const b = mk("div", "sbox"); b.style.marginTop = "10px"; const ul = mk("ul"); ul.style.cssText = "margin:0;padding-left:18px"; const gname = SIG.group === "type" ? "disease types" : "disease subtypes";
  const one = (d, name) => { const t = d.tests, pre = name ? name + ": " : "";
    ul.append(Object.assign(mk("li"), { textContent: `${pre}Kruskal–Wallis across ${gname} (donors as the unit): ${pTxt(t.kruskal_p)}.` }));
    if (t.vs_control.length) ul.append(Object.assign(mk("li"), { textContent: `${pre}versus Control (Mann–Whitney, Holm-adjusted): ` + t.vs_control.map((x) => `${x.subtype} ${fmtN(x.diff)} (${pTxt(x.p_holm)})`).join("; ") + "." }));
    if (t.spearman_rho != null) ul.append(Object.assign(mk("li"), { textContent: `${pre}rank correlation with the order ${d.order.join(" → ")}: ρ = ${t.spearman_rho.toFixed(2)} (${pTxt(t.spearman_p)}); this is the order you asked for, not an established severity scale.` })); };
  one(dA, dB ? "A" : ""); if (dB) one(dB, "B");
  const dead = SIG.scope === "sample" || SIG.scope === "subtype" || (SIG.scope === "type" && SIG.group === "type");   // the pooled groups' differences are removed by construction
  if (dead) ul.append(Object.assign(mk("li"), { textContent: "With this z-score scope the group tests above are not informative: the z-scoring removed the differences between the groups (see the warning at the top).", style: "color:#e0a458" }));
  else if (scopeWarning()) ul.append(Object.assign(mk("li"), { textContent: "With this z-score scope only part of the group comparisons remain informative (see the warning at the top).", style: "color:#e0a458" }));
  ul.append(Object.assign(mk("li"), { textContent: `Samples with fewer than ${dA.tests.min_cells} cells (after filters) are left out. With 3–7 samples per group these tests are descriptive; they do not account for cell-type composition (use the cell-type filter to look within one cell type).`, style: "color:#aab4c3" }));
  b.append(ul); return b;
}
function geneHeat(d, sg, tag) {
  const rows = d.samples.filter((s) => s.n >= d.tests.min_cells), groups = d.groups.map((g) => g.subtype), genes = d.genes;
  const mean = genes.map((_, gi) => groups.map((g) => { const v = rows.filter((r) => r.group === g && r.gene_means[gi] != null).map((r) => r.gene_means[gi]); return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null; }));
  const z = mean.map((r) => { const v = r.filter((x) => x != null), mu = v.reduce((a, b) => a + b, 0) / v.length, sd = Math.sqrt(v.reduce((a, b) => a + (b - mu) ** 2, 0) / v.length) || 1; return r.map((x) => (x == null ? null : (x - mu) / sd)); });
  const order = genes.map((_, i) => i).sort((a, b) => (z[b][z[b].length - 1] ?? 0) - (z[a][z[a].length - 1] ?? 0));
  const wrap = mk("div"); wrap.style.cssText = "max-height:340px;overflow:auto;display:inline-block;border:1px solid #29303b;border-radius:6px;margin-right:14px;vertical-align:top"; const tb = mk("table", "st"), hr = mk("tr"); hr.append(mk("th", "", (tag || "") + "Gene")); groups.forEach((g) => hr.append(mk("th", "", g))); tb.append(hr);
  for (const gi of order) {
    const tr = mk("tr"); tr.append(mk("td", "", genes[gi]));
    z[gi].forEach((v, k) => { const td = mk("td", "", mean[gi][k] == null ? "–" : mean[gi][k].toFixed(2)); if (v != null) { const c = DIVL[Math.round(((Math.max(-1.6, Math.min(1.6, v)) + 1.6) / 3.2) * 255)]; td.style.background = `rgb(${c})`; td.style.color = inkOn(`rgb(${c})`); } td.title = `${genes[gi]} · ${groups[k]}: mean log-normalised expression ${mean[gi][k] == null ? "n/a" : mean[gi][k].toFixed(3)}`; tr.append(td); });
    tb.append(tr);
  }
  wrap.append(tb); const w2 = mk("div"); w2.style.display = "inline-block"; w2.append(wrap); if (!tag || tag.startsWith("A")) w2.append(mk("div", "snote", "Cell value = mean log-normalised expression of the gene (average of the samples in the group); colour = that value scaled across the groups (red = highest, blue = lowest). Genes are ordered by how high they are in the last column.")); return w2;
}
// two signatures in one map: blue = A, orange = B, both high -> near white (a colour-blind-safe pair; additive mixing)
const dualLimits = (vecs, z) => { const pool = []; for (const v of vecs) for (let i = 0; i < v.length; i += 3) pool.push(v[i]); const hi = z ? Z_LIM : quantile(pool, 0.995), lo = z ? 0 : quantile(pool, 0.01); return { lo, hi: Math.max(hi, lo + 1e-6) }; };
function drawDualTile(cv, meta, xy, vA, vB, lA, lB, box = 260) {
  const w = meta.x[1] - meta.x[0] || 1, h = meta.y[1] - meta.y[0] || 1, k = Math.min((box - 8) / w, (box - 8) / h); cv.width = Math.max(40, Math.ceil(w * k) + 8); cv.height = Math.max(40, Math.ceil(h * k) + 8);
  const ctx = cv.getContext("2d"), n = vA.length, px = n > 30000 ? 1.5 : 2.2, base = [45, 50, 60], cA = hex2rgb(SIG_A), cB = hex2rgb(SIG_B); ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, cv.width, cv.height);
  cv._redraw = () => drawDualTile(cv, meta, xy, vA, vB, lA, lB, box);
  const hs = heFor(meta.id, "sig"); if (hs) { drawHE(ctx, { sample: hs, redo: () => heRedoCv(cv), w: "sig" }, (x) => 4 + (x - meta.x[0]) * k, (y) => 4 + (y - meta.y[0]) * k, cv.width, cv.height, k); ctx.globalAlpha = HEW.sig.a; }
  const tA = new Float32Array(n), tB = new Float32Array(n), key = new Float32Array(n); for (let i = 0; i < n; i++) { tA[i] = Math.max(0, Math.min(1, (vA[i] - lA.lo) / (lA.hi - lA.lo))); tB[i] = Math.max(0, Math.min(1, (vB[i] - lB.lo) / (lB.hi - lB.lo))); key[i] = Math.max(tA[i], tB[i]); }
  const o = new Uint32Array(n); for (let i = 0; i < n; i++) o[i] = i; o.sort((p, q) => key[p] - key[q]);
  if (HEW.sig.cells) for (const i of o) { const m = key[i]; const c = [0, 1, 2].map((q) => Math.min(255, Math.round(base[q] * (1 - m) + tA[i] * cA[q] + tB[i] * cB[q]))); ctx.fillStyle = `rgb(${c})`; ctx.fillRect(4 + (xy[2 * i] - meta.x[0]) * k, 4 + (xy[2 * i + 1] - meta.y[0]) * k, px, px); }
  ctx.globalAlpha = 1;
}
async function drawSigMaps(host, sg, sg2, dA, dB, z, tok) {
  const metas = S.meta.samples.slice().sort((a, b) => groupRank(groupLabel(a)) - groupRank(groupLabel(b)) || subtypeRank(subtypeOf(a)) - subtypeRank(subtypeOf(b)) || a.id.localeCompare(b.id)); host.replaceChildren(mk("div", "snote", "Loading maps…"));
  let vA, vB = null, xys; try { [vA, vB, xys] = await Promise.all([Promise.all(metas.map((m) => sigVec(m.id, sg.genes, SIG.sigma, SIG.scope))), sg2 ? Promise.all(metas.map((m) => sigVec(m.id, sg2.genes, SIG.sigma, SIG.scope))) : null, Promise.all(metas.map((m) => galGet(m.id, "xy")))]); } catch (e) { host.textContent = "Error: " + e.message; return; }
  if (tok !== SIG.tok) return; const meanA = Object.fromEntries(dA.samples.map((s) => [s.sample, s.mean])), meanB = dB ? Object.fromEntries(dB.samples.map((s) => [s.sample, s.mean])) : {};
  let drawTile, sc, lA, lB;
  if (sg2) {
    lA = dualLimits(vA, z); lB = dualLimits(vB, z); const leg = mk("div"), row = mk("div", "cmpleg"); row.style.flexWrap = "wrap";
    for (const [c, t] of [[SIG_A, `A · ${sg.name}`], [SIG_B, `B · ${sg2.name}`], ["#ffe0ff", "both high"]]) { const sp = mk("span"), k = mk("span", "k"); k.style.background = c; sp.append(k, document.createTextNode(t)); row.append(sp); }
    leg.append(row, mk("div", "snote", `Intensity = score ${z ? "above the average (z > 0 up to the 99.5th percentile)" : "from the 1st to the 99.5th percentile"}, scaled separately for A and B across all samples; σ = ${SIG.sigma || 0} µm.`)); host.replaceChildren(leg);
    drawTile = (cv, m, xy, i) => drawDualTile(cv, m, xy, vA[i], vB[i], lA, lB);
  } else { sc = scoreScale(vA, z); host.replaceChildren(scoreLegend(sc, `${sg.name} (σ = ${SIG.sigma || 0} µm)`)); drawTile = (cv, m, xy, i) => drawScoreTile(cv, m, xy, vA[i], sc); }
  let cur = null, grid = null;
  metas.forEach((m, i) => {
    const st = groupLabel(m); if (st !== cur) { cur = st; const sec = mk("div", "glomsec"); sec.append(mk("div", "sechead", st)); grid = mk("div", "secgrid"); sec.append(grid); host.append(sec); }
    const t = mk("div", "tile"), cv = mk("canvas"); t.append(cv); const cap = mk("div", "cap"); cap.append(Object.assign(mk("b"), { textContent: dlabel(m.id) }), document.createElement("br"),
      mk("span", "muted", `${fmt(m.n)} cells${meanA[donorOfId(m.id)] != null ? " · donor mean " + meanA[donorOfId(m.id)].toFixed(2) + (sg2 && meanB[donorOfId(m.id)] != null ? " / " + meanB[donorOfId(m.id)].toFixed(2) : "") : ""}${SIG.group === "type" ? " · " + subtypeOf(m) : ""}`)); t.append(cap); grid.append(t);
    t.addEventListener("click", async () => { $("sample").value = m.id; closeSig(); await loadSample(m.id); }); drawTile(cv, m, xys[i], i);
  });
}

// ---------------------------------------------------------------- type I vs II localisation
async function renderLoc() {
  const pane = $("sigpane"), tok = ++SIG.tok; pane.replaceChildren();
  pane.append(Object.assign(mk("div", "sh2", "Do two signatures localise to different regions of the tissue?"), { style: "margin-top:0" }));
  pane.append(mk("div", "sdesc", "Both signatures are smoothed within each tissue piece and z-scored within each donor, so each has mean 0 and SD 1 and their difference is in SD units. The test compares the two maps (correlation, overlap of their hottest cells) and asks, region by region and cell type by cell type, whether one signature is systematically higher than the other. Samples are the unit of replication."));
  const ctl = mk("div", "sctl"), opts = (cur) => SIG_GROUPS.map((g) => `<optgroup label="${g}">${SIG.list.filter((x) => x.group === g).map((x) => `<option value="${x.id}" ${x.id === cur ? "selected" : ""}>${x.name}</option>`).join("")}</optgroup>`).join("");
  const sa = mk("select"), sb = mk("select"); sa.innerHTML = opts(SIG.a); sb.innerHTML = opts(SIG.b); const la = mk("label"), lb = mk("label"); la.append("A ", sa); lb.append("B ", sb);
  const sl = mk("input"); sl.type = "range"; sl.min = 0; sl.max = 100; sl.step = 5; sl.value = SIG.locSigma; const sv = mk("span", "", SIG.locSigma ? SIG.locSigma + " µm" : "off"); const l1 = mk("label"); l1.append("Smoothing σ ", sl, sv);
  const q = mk("input"); q.type = "number"; q.min = 5; q.max = 30; q.step = 5; q.value = SIG.locQ; q.style.width = "56px"; const lq = mk("label"); lq.append("Hotspot = top ", q, "%");
  const ss = mk("select"); ss.innerHTML = `<option value="">all donors</option>` + (SIG.order || []).map((t) => `<option value="${t}" ${t === SIG.locSub ? "selected" : ""}>${t} only</option>`).join(""); const ls = mk("label"); ls.append("Samples ", ss);
  const run = mk("button", "", "Run"); ctl.append(la, lb, l1, lq, ls, run); pane.append(ctl);
  sl.addEventListener("input", () => { sv.textContent = +sl.value ? sl.value + " µm" : "off"; });
  const go = () => { SIG.a = sa.value; SIG.b = sb.value; SIG.locSigma = +sl.value; SIG.locQ = +q.value || 10; SIG.locSub = ss.value; renderLoc(); };
  run.addEventListener("click", go); [sa, sb, ss].forEach((e) => e.addEventListener("change", go)); sl.addEventListener("change", go); q.addEventListener("change", go);
  const out = mk("div", "snote", "Computing…"); pane.append(out);
  let d; try { d = await postJSON("/api/signature/localization", { a: SIG.a, b: SIG.b, sigma: SIG.locSigma, q: SIG.locQ / 100, subtype: SIG.locSub || null }); } catch (e) { out.textContent = "Error: " + e.message; return; }
  if (tok !== SIG.tok) return; out.remove();
  const A = SIG.byId[SIG.a], B = SIG.byId[SIG.b], na = A.name.replace(/ \(.*/, "").replace("IFN ", "IFN "), nb = B.name.replace(/ \(.*/, "");
  const sig = (rows, key) => rows.filter((r) => r[key] != null && r[key] < 0.05), o = d.overall;
  const sum = mk("div", "sbox"), ul = mk("ul"); ul.style.cssText = "margin:0;padding-left:18px";
  const li = (t) => ul.append(Object.assign(mk("li"), { textContent: t }));
  li(`Across ${donTxt(d.n_donors)} the two maps are ${o.r_mean > 0.7 ? "strongly" : o.r_mean > 0.4 ? "moderately" : "weakly"} correlated: mean r = ${o.r_mean.toFixed(2)} (r = ${o.r_adj_mean.toFixed(2)} after removing each cell type's mean, so the shared part is not just cell-type identity).`);
  li(`The hottest ${SIG.locQ}% of cells overlap ${o.overlap_ratio_mean.toFixed(1)}× more than chance (1× = independent, ${(100 / SIG.locQ).toFixed(0)}× = identical maps).`);
  const rA = sig(d.regions, "d_adj_fdr"), cA = sig(d.cell_types, "d_fdr"), rR = sig(d.regions, "d_fdr");
  const fmtR = (rows, key, k) => rows.length ? rows.map((r) => `${r.name} (${fmtN(r[k])}, FDR ${fmtP(r[key])})`).join("; ") : "none";
  li(`Regions where ${na} is higher than ${nb} after adjusting for cell type (FDR < 0.05): ${fmtR(rA.filter((r) => r.d_adj > 0), "d_adj_fdr", "d_adj")}. Where ${nb} is higher: ${fmtR(rA.filter((r) => r.d_adj < 0), "d_adj_fdr", "d_adj")}.`);
  li(`Without the cell-type adjustment the significant regions are: ${fmtR(rR, "d_fdr", "d")}.`);
  li(`Cell types where one signature is significantly higher (raw contrast, FDR < 0.05; positive = ${na} higher): ${cA.length ? cA.sort((a, b) => b.d - a.d).map((r) => `${r.name} (${fmtN(r.d)})`).join("; ") : "none"}.`);
  const small = d.regions.filter((r) => r.n_donors < 6).map((r) => r.name); li(`Caveat: regions seen in fewer than 6 donors cannot reach significance (${small.join(", ") || "none here"}); nested regions are analysed independently; ${d.subtype ? "only " + d.subtype + " samples" : "all samples"} pooled.`);
  sum.append(ul); pane.append(sum);
  pane.append(mk("div", "sh2", "Where does each signature sit? Contrast d = mean(A) − mean(B) over the cells inside, in SD units")); pane.append(forest(d.regions, na, nb, "Pathologist regions")); pane.append(forest(d.cell_types, na, nb, "Cell types"));
  pane.append(mk("div", "sh2", "Per donor")); pane.append(locSampleTable(d));
  pane.append(mk("div", "sh2", "Maps of one donor")); const mp = mk("div"); pane.append(mp); await locMaps(mp, A, B, tok);
}
function forest(rows, na, nb, title) {
  const w = mk("div"); w.append(mk("div", "snote", title + (rows.length && rows[0].d_adj == null ? ` — raw contrast only (removing a cell type's own mean would make it zero);` : ` — hollow = raw contrast, filled = after removing cell-type means;`) + ` blue = ${na} higher, orange = ${nb} higher (FDR < 0.05), grey = not significant`));
  if (!rows.length) { w.append(mk("div", "snote", "Too few donors with enough cells.")); return w; }
  rows = rows.slice().sort((a, b) => (b.d_adj ?? b.d) - (a.d_adj ?? a.d)); const W = 760, rh = 28, m = { l: 200, r: 150, t: 24, b: 30 }, H = m.t + m.b + rows.length * rh, svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` });
  const vals = rows.flatMap((r) => [r.d, ...(r.d_ci || []), r.d_adj, ...(r.d_adj_ci || [])]).filter((v) => v != null), mx = Math.max(0.2, Math.max(...vals.map(Math.abs))) * 1.08, xs = (v) => m.l + ((v + mx) / (2 * mx)) * (W - m.l - m.r);
  for (const t of [-mx, -mx / 2, 0, mx / 2, mx]) { svgEl("line", { x1: xs(t), x2: xs(t), y1: m.t - 6, y2: H - m.b, stroke: t === 0 ? "#6b7785" : "#29303b", "stroke-width": 1 }, svg); svgEl("text", { x: xs(t), y: H - m.b + 16, "text-anchor": "middle", fill: "#8b95a5", "font-size": 11 }, svg).textContent = t.toFixed(2); }
  svgEl("text", { x: xs(0) - 8, y: 12, "text-anchor": "end", fill: PAL[1], "font-size": 11 }, svg).textContent = `← ${nb} higher`; svgEl("text", { x: xs(0) + 8, y: 12, fill: PAL[0], "font-size": 11 }, svg).textContent = `${na} higher →`;
  rows.forEach((r, i) => {
    const y = m.t + i * rh + rh / 2; svgEl("text", { x: m.l - 10, y: y + 4, "text-anchor": "end", fill: "#e3e8ef", "font-size": 12 }, svg).textContent = r.name.length > 30 ? r.name.slice(0, 29) + "…" : r.name;
    for (const [d, ci, fdr, dy, filled] of (r.d_adj == null ? [[r.d, r.d_ci, r.d_fdr, 0, true]] : [[r.d, r.d_ci, r.d_fdr, -5, false], [r.d_adj, r.d_adj_ci, r.d_adj_fdr, 5, true]])) {
      const col = fdr != null && fdr < 0.05 ? (d > 0 ? PAL[0] : PAL[1]) : "#8b95a5"; if (ci && ci[0] != null) svgEl("line", { x1: xs(ci[0]), x2: xs(ci[1]), y1: y + dy, y2: y + dy, stroke: col, "stroke-width": 1.6 }, svg);
      const c = svgEl("circle", { cx: xs(d), cy: y + dy, r: 4.5, fill: filled ? col : "#0f1218", stroke: col, "stroke-width": 2 }, svg);
      const hit = svgEl("circle", { cx: xs(d), cy: y + dy, r: 11, fill: "transparent" }, svg);
      hit.addEventListener("pointermove", (ev) => nbrTip(ev, [{ tag: "b", t: `d = ${fmtN(d)}` }, { t: `${r.name} · ${filled ? "cell-type-adjusted" : "raw"}` }, { cls: "m", t: `95% CI ${ci && ci[0] != null ? fmtN(ci[0]) + " to " + fmtN(ci[1]) : "n/a"} · FDR ${fmtP(fdr)}` }, { cls: "m", t: `${na} ${fmtN(r.mean_a)} · ${nb} ${fmtN(r.mean_b)} · ${donTxt(r.n_donors)}, ${fmt(r.n_cells)} cells` }])); hit.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
    }
    svgEl("text", { x: W - m.r + 12, y: y + 4, fill: "#8b95a5", "font-size": 11 }, svg).textContent = `${donTxt(r.n_donors)} · FDR ${fmtP(r.d_adj == null ? r.d_fdr : r.d_adj_fdr)}`;
  });
  w.append(svg); return w;
}
function locSampleTable(d) {
  const tb = mk("table", "st"), hr = mk("tr"); ["Sample", "Cells", "r (A, B)", "r, cell-type-adjusted", "Hotspot overlap (× chance)"].forEach((t) => hr.append(mk("th", "", t))); tb.append(hr);
  let cur = null; for (const r of d.samples.slice().sort((a, b) => subtypeRank(a.subtype) - subtypeRank(b.subtype) || a.sample.localeCompare(b.sample))) {
    if (r.subtype !== cur) { cur = r.subtype; const tr = mk("tr", "sub"); const t0 = mk("td", "", cur); t0.colSpan = 5; tr.append(t0); tb.append(tr); }
    const tr = mk("tr"); [dlabel(r.sample), fmt(r.n), r.r.toFixed(2), r.r_adj.toFixed(2), r.overlap_ratio != null ? r.overlap_ratio.toFixed(1) + "×" : "–"].forEach((t) => tr.append(mk("td", "", t))); tb.append(tr);
  } const w = mk("div"); w.style.cssText = "max-height:340px;overflow:auto;display:inline-block;border:1px solid #29303b;border-radius:6px"; w.append(tb); return w;
}
async function locMaps(host, A, B, tok) {
  const ids = S.meta.samples.slice().sort((a, b) => subtypeRank(subtypeOf(a)) - subtypeRank(subtypeOf(b)) || a.id.localeCompare(b.id));
  const pick = mk("select"); pick.innerHTML = ids.map((m) => `<option value="${m.id}">${dlabel(m.id)} · ${subtypeOf(m)}</option>`).join(""); pick.value = SIG.mapSample && ids.some((m) => m.id === SIG.mapSample) ? SIG.mapSample : S.sample;
  const row = mk("div", "sctl"), lab = mk("label"); lab.append("Donor ", pick); row.append(lab); host.append(row); const grid = mk("div", "mapgrid"); host.append(grid);
  const draw = async () => {
    SIG.mapSample = pick.value; grid.replaceChildren(mk("div", "snote", "Loading maps…")); const meta = S.meta.samples.find((m) => m.id === pick.value);
    let a, b, xy; try { [a, b, xy] = await Promise.all([sigVec(meta.id, A.genes, SIG.locSigma, "sample"), sigVec(meta.id, B.genes, SIG.locSigma, "sample"), galGet(meta.id, "xy")]); } catch (e) { grid.textContent = "Error: " + e.message; return; }
    if (tok !== SIG.tok) return; const diff = new Float32Array(a.length); for (let i = 0; i < a.length; i++) diff[i] = a[i] - b[i];
    const sAB = scoreScale([a, b], true), sD = scoreScale([diff], true); grid.replaceChildren();
    for (const [name, vec, sc] of [[A.name, a, sAB], [B.name, b, sAB], [`${A.name} − ${B.name}`, diff, sD]]) { const f = mk("figure"), cv = mk("canvas"); f.append(Object.assign(mk("figcaption"), { textContent: name }), cv, scoreLegend(sc, "")); grid.append(f); drawScoreTile(cv, meta, xy, vec, sc, 330); }
  };
  pick.addEventListener("change", draw); await draw();
}

// ---------------------------------------------------------------- violins: one violin per sample, by cell type or by cell state
// Every panel = one gene / signature in one cell type; samples sit side by side in disease-subtype order. In cell-state mode a panel holds all
// cell states of one cell type (one group of per-sample violins per state). Distributions are computed on the server (density, quartiles).
const VIO = { feats: [], cts: new Set([0]), level: "cell_type", unit: "sample", sigma: 0, scope: "none", nz: false, tok: 0, type: "violin", logy: false };
const vioKey = (f) => (f.sig ? "s:" + f.sig : "g:" + f.gene), vioName = (f) => (f.sig ? (SIG.byId[f.sig]?.name || f.sig) : f.gene);
async function openVio() {
  closeAllWindows();
  try { await loadSigList(); } catch (e) { /* signatures optional */ }
  const sel = $("viosig"); if (sel.options.length < 2 && SIG.list) for (const g of SIG_GROUPS) { const og = document.createElement("optgroup"); og.label = g;
    SIG.list.filter((x) => x.group === g).forEach((x) => { const o = document.createElement("option"); o.value = x.id; o.textContent = `${x.name} (${x.genes.length})`; og.append(o); }); sel.append(og); }
  if (!VIO.feats.length) { if (SIG.byId?.ifn1) VIO.feats.push({ sig: "ifn1" }); else { const l = await jget("/api/genes/search?q="); if (l.length) VIO.feats.push({ gene: l[0] }); } }   // projects without signatures.json: start from the first panel gene
  $("vio").hidden = false; vioControls(); renderVio();
}
function closeVio() { $("vio").hidden = true; }
function vioControls() {
  const chips = $("viochips"); chips.replaceChildren();
  VIO.feats.forEach((f, i) => { const c = mk("span", "chip2", vioName(f) + (f.sig ? " (signature)" : "")), x = mk("span", "x", "×"); x.title = "remove"; x.addEventListener("click", () => { VIO.feats.splice(i, 1); vioControls(); renderVio(); }); c.append(x); chips.append(c); });
  const names = S.meta.cats.cell_type.levels, list = $("vioctlist"); list.replaceChildren();
  const addOpt = (id, name, color, off) => { const l = mk("label"), cb = mk("input"); cb.type = "checkbox"; cb.checked = VIO.cts.has(id); cb.disabled = !!off; const sw = mk("span", "sw"); sw.style.cssText = `width:11px;height:11px;border-radius:2px;background:${color}`;
    cb.addEventListener("change", () => { cb.checked ? VIO.cts.add(id) : VIO.cts.delete(id); vioControls(); renderVio(); }); l.append(cb, sw, name + (off ? " (not in cell-state mode)" : "")); list.append(l); };
  names.forEach((n, i) => addOpt(i, n, catColor("cell_type", i), false)); addOpt(-1, "All cells (pooled)", "#6b7785", VIO.level === "cell_state");
  const bar = mk("div"); bar.style.marginTop = "6px"; for (const [t, fn] of [["all types", () => { VIO.cts = new Set(names.map((_, i) => i)); }], ["none", () => { VIO.cts = new Set(); }]]) { const b = mk("button", "sm", t); b.addEventListener("click", (e) => { e.preventDefault(); fn(); vioControls(); renderVio(); }); bar.append(b, " "); } list.append(bar);
  if (VIO.level === "cell_state") VIO.cts.delete(-1);
  const sel = [...VIO.cts].sort((a, b) => a - b).map((i) => (i < 0 ? "All cells" : names[i].replace(/ cells$/, "")));
  $("vioctssum").textContent = `Cell types: ${sel.length ? (sel.length > 3 ? sel.length + " selected" : sel.join(", ")) : "none"} ▾`;
  $("vioshow").value = VIO.level; $("vioscope").value = VIO.scope; $("vionz").checked = VIO.nz; $("vionit").value = VIO.unit; $("viotype").value = VIO.type; $("viology").checked = VIO.logy; $("viology").parentElement.hidden = VIO.type !== "hist"; $("viohead").querySelector(".muted").textContent = VIO.unit === "sample" ? "one violin per donor (tissue pieces pooled), grouped by disease subtype" : `one bar per disease ${VIO.unit} = mean of its donors' means, grouped by cell type`;
}
async function vioAddGene() {
  const v = $("viogene").value.trim(); if (!v) return; const hit = (await jget("/api/genes/search?q=" + encodeURIComponent(v) + "&limit=1"))[0];
  if (!hit || hit.toLowerCase() !== v.toLowerCase()) { $("viomsg").textContent = `“${v}” is not on the panel`; return; }
  $("viomsg").textContent = ""; if (!VIO.feats.some((f) => f.gene === hit)) VIO.feats.push({ gene: hit }); $("viogene").value = ""; vioControls(); renderVio();
}
async function renderVio() {
  if ($("vio").hidden) return; const tok = ++VIO.tok, body = $("viobody"), cts = [...VIO.cts].sort((a, b) => a - b);
  if (!VIO.feats.length || !cts.length) { body.replaceChildren(mk("div", "snote", "Add at least one gene or signature and choose at least one cell type.")); return; }
  if (VIO.feats.length * cts.length > 24) { body.replaceChildren(mk("div", "snote", `That is ${VIO.feats.length * cts.length} panels; please select fewer genes / signatures / cell types (at most 24).`)); return; }
  body.replaceChildren(mk("div", "snote", "Computing…"));
  if (VIO.type === "hist") return renderVioHist(tok, cts);
  let d; try { d = await postJSON("/api/violin", { features: VIO.feats.map((f) => (f.sig ? { sig: f.sig } : { gene: f.gene })), cell_types: cts, sigma: 0, layer: S.layer, scope: VIO.scope, nonzero_only: VIO.nz, level: VIO.level, unit: VIO.unit }); }
  catch (e) { body.replaceChildren(mk("div", "snote", "Error: " + e.message)); return; }
  if (tok !== VIO.tok) return; body.replaceChildren();
  let meta = Object.fromEntries(d.samples.map((s) => [s.id, s])), order = d.samples.slice().sort((a, b) => subtypeRank(a.subtype) - subtypeRank(b.subtype) || a.id.localeCompare(b.id)).map((s) => s.id), panels = d.panels;
  if (VIO.unit !== "sample") {   // pooled violins: one per disease subtype / type; all selected cell types go into ONE plot per feature
    const names = d.panels[0].groups[0].samples.map((r) => r.sample), rk = (n) => (VIO.unit === "type" ? NB_TYPE_ORDER.indexOf(n) : subtypeRank(n));
    order = names.slice().sort((a, b) => rk(a) - rk(b)); meta = Object.fromEntries(names.map((n) => [n, { subtype: n }]));
    if (VIO.level === "cell_type") panels = [...new Set(d.panels.map((p) => p.feature))].map((f) => { const ps = d.panels.filter((p) => p.feature === f);
      return { feature: f, cell_type: ps.length > 1 ? "all selected cell types" : ps[0].cell_type, cell_type_id: ps[0].cell_type_id, level: "cell_type", lo: ps[0].lo, hi: ps[0].hi, merged: ps.length > 1,
        groups: ps.map((p) => ({ ...p.groups[0], name: p.cell_type, id: p.cell_type_id })) }; });
  }
  for (const p of panels) body.append(VIO.unit !== "sample" ? barPanel(p, meta, order) : violinPanel(p, meta, order));
  const note = mk("div", "snote", (VIO.unit !== "sample" ? `Each bar is the mean over the donors of one disease ${VIO.unit} of the donor's own mean (every donor counts once, however many cells it has); error bar = standard error over donors; dots = donors. ` : "Each violin is one donor (its tissue pieces pooled). ") + "(Disease subtype order: Control, No LN, Class I/II, Class III, Class V). Dot = median, dark bar = interquartile range; every violin is scaled to the same maximum width, so read its cell count from the tooltip. A grey tick means fewer than 10 cells. The axis covers the 0.5–99.5th percentile of the plotted cells (values beyond are drawn at the edge). Many cells have a score of 0 for sparse genes: tick “expressing cells only” to see the distribution among expressing cells (the tooltip gives the percentage expressing).");
  note.style.maxWidth = "900px"; body.append(note);
}
// Histograms of the per-cell values behind the violins (all donors pooled): the x-axis is the same scale as the violins' y-axis, and for a single gene
// a second histogram shows its raw counts in the same cells, so the transformation can be checked by eye.
async function renderVioHist(tok, cts) {
  const body = $("viobody"); let d;
  try { d = await postJSON("/api/violin/hist", { features: VIO.feats.map((f) => (f.sig ? { sig: f.sig } : { gene: f.gene })), cell_types: cts, sigma: 0, layer: S.layer, scope: VIO.scope, nonzero_only: VIO.nz, level: "cell_type", unit: VIO.unit }); }
  catch (e) { body.replaceChildren(mk("div", "snote", "Error: " + e.message)); return; }
  if (tok !== VIO.tok) return; body.replaceChildren();
  const hist = (counts, labels, w, h, col, zeroFirst, xlab, o = {}) => {   // bars = cells (or, with o.frac, the donor-averaged share of cells); zeroFirst: the first bar (exact zero) is drawn in grey
    const L = 50, T = 8, B = 38, W = w, H = h, tot = counts.reduce((a, c) => a + c, 0) || 1, f = (c) => (VIO.logy ? Math.log10(1 + (o.frac ? c * 1e4 : c)) : c), mx = Math.max(o.ymax != null ? f(o.ymax) : 0, ...counts.slice(zeroFirst ? 1 : 0).map(f), 1e-9);   // the zero bar is clipped when it dwarfs the rest (its value is printed)
    const svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` }), bw = (W - L - 10) / counts.length;
    counts.forEach((c, i) => { const bh = Math.min(f(c) / mx, 1) * (H - T - B), clip = f(c) > mx * 1.0001; r = svgEl("rect", { x: L + i * bw + 0.5, y: H - B - bh, width: Math.max(1, bw - 1), height: bh, fill: zeroFirst && i === 0 ? "#6b7785" : col }, svg);
      if (clip) svgEl("text", { x: L + i * bw + bw + 4, y: T + 10, fill: "#aab4c3", "font-size": 10 }, svg).textContent = `← cut off at the top: ${o.frac ? (100 * c).toFixed(1) + "% of cells" : fmt(c) + " cells"}`;
      svgEl("title", {}, r).textContent = o.frac ? `${labels[i]}: ${(100 * c).toFixed(c < 0.01 ? 3 : 2)}% of a donor's cells (mean over donors)` : `${labels[i]}: ${fmt(c)} cells (${(100 * c / tot).toFixed(c / tot < 0.01 ? 2 : 1)}%)`; });
    svgEl("line", { x1: L, x2: W - 10, y1: H - B, y2: H - B, stroke: "#39424f" }, svg);
    const ticks = counts.length <= 12 ? counts.length - 1 : 5; for (let k = 0; k <= ticks; k++) { const i = Math.round((k / ticks) * (counts.length - 1)); svgEl("text", { x: L + i * bw + bw / 2, y: H - B + 14, "text-anchor": "middle", fill: "#8b95a5", "font-size": 10 }, svg).textContent = labels[i].split("–")[0]; }
    svgEl("text", { x: L + (W - L - 10) / 2, y: H - 6, "text-anchor": "middle", fill: "#8b95a5", "font-size": 10.5 }, svg).textContent = xlab;
    svgEl("text", { x: 12, y: T + (H - T - B) / 2, "text-anchor": "middle", fill: "#8b95a5", "font-size": 10.5, transform: `rotate(-90 12 ${T + (H - T - B) / 2})` }, svg).textContent = o.frac ? (VIO.logy ? "log10(1 + 10⁴ × share)" : "% of cells (donor mean)") : VIO.logy ? "log10(1 + cells)" : "cells";
    svgEl("text", { x: L - 6, y: T + 8, "text-anchor": "end", fill: "#8b95a5", "font-size": 10 }, svg).textContent = o.frac ? (100 * (VIO.logy ? (10 ** mx - 1) / 1e4 : mx)).toFixed(2) + "%" : VIO.logy ? fmt(Math.round(10 ** mx - 1)) : fmt(Math.round(mx));
    return svg;
  };
  for (const p of d.panels) {
    const box = mk("div", "vpanel"); box.style.cssText = "margin:0 0 18px;max-width:980px";
    const nzv = p.n - (p.n_zero || 0), fz = p.n ? (100 * (p.n_zero || 0) / p.n).toFixed(1) : "0";
    box.append(Object.assign(mk("div", "sh2"), { textContent: `${p.feature} · ${p.cell_type} · ${fmt(p.n)} cells${VIO.unit === "sample" ? ", all donors pooled" : ""}` }));
    if (!p.n) { box.append(mk("div", "snote", "No cells.")); body.append(box); continue; }
    box.append(mk("div", "snote", `${fz}% of these cells have exactly 0 (grey bar); ${fmt(nzv)} cells are > 0. Plotted values: mean ${p.mean.toFixed(3)}, median ${p.median.toFixed(3)}, max ${p.max.toFixed(3)}${p.min_nonzero != null ? `, smallest non-zero ${p.min_nonzero.toFixed(3)}` : ""}.`));
    const ecol = p.cell_type_id < 0 ? "#8b95a5" : catColor("cell_type", p.cell_type_id), lab = (i) => `${p.edges[i]}–${p.edges[i + 1]}`;
    const labels = ["0", ...p.counts.map((_, i) => lab(i))], xl = "value plotted on the violins' y-axis: " + vioYLabel() + " (bar 1 = exactly 0, the rest are 60 equal bins of the non-zero values)";
    const rl = p.raw ? p.raw.counts.map((_, i) => (i === p.raw.top && p.raw.max > p.raw.top ? `${i}+` : String(i))) : [];
    const rk = (n) => (VIO.unit === "type" ? NB_TYPE_ORDER.indexOf(n) : subtypeRank(n)), gs = p.groups ? p.groups.slice().sort((a, b) => rk(a.name) - rk(b.name)) : null;
    if (gs) {   // one histogram per disease subtype / type: the average over its donors of each donor's own histogram
      const okg = gs.filter((g) => g.share), ym = Math.max(...okg.flatMap((g) => g.share.slice(1)), 1e-9), rm = p.raw ? Math.max(...okg.flatMap((g) => (g.raw_share || []).slice(1)), 1e-9) : 0;
      box.append(mk("div", "snote", "One histogram per disease " + VIO.unit + ": the average over its donors of each donor's own histogram (share of that donor's cells per bin; every donor counts once; donors with < 10 cells are left out). Same y-scale in every histogram; the grey zero bar is cut off at the top when it is much taller than the rest (its value is printed)."));
      for (const g of gs) {
        box.append(Object.assign(mk("div", "vtitle"), { textContent: `${g.name} · ${g.n_donors_used} donor${g.n_donors_used === 1 ? "" : "s"}${g.n_donors_used < g.n_donors ? ` (of ${g.n_donors})` : ""} · ${fmt(g.n)} cells${g.mean != null ? ` · mean of donor means ${g.mean.toFixed(3)}` : ""}` }));
        if (!g.share) { box.append(mk("div", "snote", "Too few cells.")); continue; }
        box.append(hist(g.share, labels, 900, 130, ecol, true, xl, { frac: true, ymax: ym }));
        if (p.raw && g.raw_share) box.append(hist(g.raw_share, rl, 900, 100, "#d9a441", true, "raw transcript counts per cell", { frac: true, ymax: rm }));
      }
    } else box.append(hist([p.n_zero || 0, ...p.counts], labels, 900, 200, ecol, true, xl));
    if (p.raw) {
      if (!gs) { box.append(mk("div", "snote", `Same cells, raw counts of ${p.feature} (integers; max ${fmt(p.raw.max)}, mean ${p.raw.mean.toFixed(3)}, ${(100 * p.raw.frac_zero).toFixed(1)}% zero):`));
        box.append(hist(p.raw.counts, rl, 900, 150, "#d9a441", true, "raw transcript counts per cell")); }
      else box.append(mk("div", "snote", `Orange = the same cells' raw counts of ${p.feature} (integers; pooled max ${fmt(p.raw.max)}).`));
      box.append(mk("div", "snote", p.check_max_abs_diff != null ? `Check: for all ${fmt(p.n)} cells the plotted value equals log1p(raw count / total counts × ${fmt(d.norm_scale)}); largest absolute difference = ${p.check_max_abs_diff.toExponential(1)}.`
        : p.layer === "counts" ? "Layer = raw counts: the plotted values are these integer counts." : "Smoothing / z-scoring is applied, so the plotted values are no longer a direct function of the raw counts."));
    } else box.append(mk("div", "snote", "Signature score = mean of the genes' values on the chosen layer" + (S.layer === "counts" ? " (raw counts)" : ", each log1p(counts / total counts × " + fmt(d.norm_scale) + ")") + (VIO.scope !== "none" ? ", then z-scored" : "") + "."));
    body.append(box);
  }
  const note = mk("div", "snote", (VIO.unit === "sample" ? "Histograms pool all donors (choose “disease subtype” or “disease type” under “One violin per” for donor-averaged histograms per group)" : "Histograms are averaged over donors") + "; they ignore the cell-state option. Cell types, score, z-score and expressing-cells options are the same as for the violins."); note.style.maxWidth = "900px"; body.append(note);
}
function vioYLabel() {   // follows the main viewer's expression layer; z-scores say which scope they use
  const lay = S.layer === "counts" ? "raw counts" : `log1p(counts/total×${fmt(S.meta.norm_scale)})`;
  const zl = { sample: "within each donor", subtype: "within disease subtype", type: "within disease type", all: "across all samples", nocontrol: "vs non-control samples" }[VIO.scope];
  return VIO.scope === "none" ? lay : `z-score of ${lay} (${zl})`;
}
// Disease subtype / type view: one bar per subtype (or type) = mean over its donors of each donor's mean score; whiskers = SE over donors; dots = the donors.
function barPanel(p, meta, order) {
  const many = p.groups.length > 1, bw = many ? 16 : 40, gap = 8, ggap = 26, L = 56, T = 24, PH = 190, B = 46;
  const pos = []; let x = L; for (const g of p.groups) { const gx = []; order.forEach(() => { gx.push(x + bw / 2); x += bw + gap; }); pos.push(gx); x += ggap - gap; }
  const W = x - ggap + 16, H = T + PH + B, svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` });
  const vals = []; for (const g of p.groups) for (const r of g.samples) if (r.mean != null) { vals.push(r.mean + (r.se || 0), r.mean - (r.se || 0), ...r.donor_means); }
  let lo = Math.min(0, ...vals), hi = Math.max(0, ...vals); if (hi - lo < 1e-9) { hi = lo + 1; } const pad = 0.06 * (hi - lo); hi += pad; if (lo < 0) lo -= pad;
  const ys = (v) => T + (1 - (v - lo) / (hi - lo)) * PH, step = (() => { const r = (hi - lo) / 5, q = Math.pow(10, Math.floor(Math.log10(r))); return [1, 2, 5, 10].map((k) => k * q).find((k) => k >= r) || r; })();
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) { svgEl("line", { x1: L - 4, x2: W - 10, y1: ys(v), y2: ys(v), stroke: "#29303b", "stroke-width": 1 }, svg); svgEl("text", { x: L - 8, y: ys(v) + 4, "text-anchor": "end", fill: "#8b95a5", "font-size": 10 }, svg).textContent = +v.toFixed(3); }
  svgEl("line", { x1: L - 4, x2: W - 10, y1: ys(0), y2: ys(0), stroke: "#6b7785", "stroke-width": 1 }, svg);
  svgEl("text", { x: 12, y: T + PH / 2, fill: "#8b95a5", "font-size": 10.5, transform: `rotate(-90 12 ${T + PH / 2})`, "text-anchor": "middle" }, svg).textContent = "mean of per-donor means: " + vioYLabel();
  p.groups.forEach((g, gi) => {
    const col = p.merged ? (g.id < 0 ? "#8b95a5" : catColor("cell_type", g.id)) : many ? catColor("cell_state", g.id) : catColor("cell_type", p.cell_type_id < 0 ? 0 : p.cell_type_id), colr = p.cell_type_id < 0 && !many ? "#8b95a5" : col;
    g.samples.forEach((r) => { const i = order.indexOf(r.sample), cx = pos[gi][i];
      if (r.mean != null) {
        svgEl("rect", { x: cx - bw / 2, y: Math.min(ys(0), ys(r.mean)), width: bw, height: Math.max(1, Math.abs(ys(r.mean) - ys(0))), fill: colr, "fill-opacity": 0.85 }, svg);
        if (r.se != null) { for (const [a, b2] of [[r.mean - r.se, r.mean + r.se]]) svgEl("line", { x1: cx, x2: cx, y1: ys(a), y2: ys(b2), stroke: "#e3e8ef", "stroke-width": 1.5 }, svg); for (const v of [r.mean - r.se, r.mean + r.se]) svgEl("line", { x1: cx - 4, x2: cx + 4, y1: ys(v), y2: ys(v), stroke: "#e3e8ef", "stroke-width": 1.5 }, svg); }
        r.donor_means.forEach((v, j) => svgEl("circle", { cx: cx + ((j * 0.618 % 1) - 0.5) * bw * 0.7, cy: ys(v), r: 2.4, fill: "#0f1218", stroke: "#fff", "stroke-width": 0.9, "fill-opacity": 0.9 }, svg));
      } else svgEl("line", { x1: cx - bw / 2, x2: cx + bw / 2, y1: ys(0), y2: ys(0), stroke: "#4b566a", "stroke-width": 2 }, svg);
      const tip = (ev) => nbrTip(ev, [{ tag: "b", t: r.mean != null ? "mean " + r.mean.toFixed(3) + (r.se != null ? " ± " + r.se.toFixed(3) + " (SE)" : "") : "too few cells" }, { t: `${r.sample} · ${r.n_donors_used} of ${r.n_donors} donors used (each needs ≥ 10 cells)` }, { cls: "m", t: `${many ? g.name : p.cell_type} · ${fmt(r.n)} cells${r.frac_nonzero == null ? "" : " · " + (100 * r.frac_nonzero).toFixed(1) + "% expressing (mean over donors)"}` }, ...(r.donors ? [{ cls: "m", t: r.donors.map((d, j) => `${dlabel(d)}: ${r.donor_means[j].toFixed(3)}`).join(" · ") }] : [])]);
      const hit = svgEl("rect", { x: cx - bw / 2 - gap / 2, y: T, width: bw + gap, height: PH, fill: "transparent" }, svg); hit.addEventListener("pointermove", tip); hit.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
      if (!many) svgEl("text", { x: cx, y: T + PH + 16, "text-anchor": "middle", fill: "#8b95a5", "font-size": 10.5 }, svg).textContent = `${r.sample} (${r.n_donors})`;
    });
    if (many) { const gx = pos[gi]; svgEl("text", { x: (gx[0] + gx[gx.length - 1]) / 2, y: T + PH + 16, "text-anchor": "middle", fill: "#e3e8ef", "font-size": 11.5, "font-weight": 600 }, svg).textContent = p.merged ? g.name.replace(/ cells$/, "") : g.name.split(".")[0];
      order.forEach((nm, i) => svgEl("text", { x: gx[i], y: T + PH + 30, "text-anchor": "middle", fill: "#6b7785", "font-size": 8.5 }, svg).textContent = nm.replace("Class ", "").replace("Control", "Ctl").replace("No LN", "NoLN")); }
  });
  const w = mk("div", "vpanel"); w.append(mk("div", "vtitle", `${p.feature} — ${p.cell_type}${many && !p.merged ? " (cell states)" : ""}`));
  if (many) { const leg = mk("div", "vleg"); p.groups.forEach((g) => { const sp = mk("span"), k = mk("span", "k"); k.style.background = catColor(p.merged ? "cell_type" : "cell_state", g.id); sp.append(k, g.name); leg.append(sp); }); w.append(leg); }
  w.append(svg); return w;
}
function violinPanel(p, meta, order) {
  const many = p.groups.length > 1, pooled = VIO.unit !== "sample", vw = pooled ? 30 : many ? 9 : 24, hw = vw * 0.46, subs = order.map((id) => meta[id].subtype), gap = 5, ggap = 30, L = 56, T = 28, PH = 190, B = many ? 48 : 78;
  // x layout: groups -> samples in subtype order (small gap between subtypes)
  const pos = []; let x = L; for (const g of p.groups) { const gx = []; subs.forEach((sb, i) => { if (i && sb !== subs[i - 1]) x += gap; gx.push(x + vw / 2); x += vw; }); pos.push(gx); x += ggap; }
  const W = x - ggap + 16, H = T + PH + B, svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` }), ys = (v) => T + (1 - (v - p.lo) / (p.hi - p.lo)) * PH;
  const step = (() => { const r = (p.hi - p.lo) / 5, q = Math.pow(10, Math.floor(Math.log10(r))); return [1, 2, 5, 10].map((k) => k * q).find((k) => k >= r) || r; })();
  for (let v = Math.ceil(p.lo / step) * step; v <= p.hi + 1e-9; v += step) { svgEl("line", { x1: L - 4, x2: W - 10, y1: ys(v), y2: ys(v), stroke: "#29303b", "stroke-width": 1 }, svg); svgEl("text", { x: L - 8, y: ys(v) + 4, "text-anchor": "end", fill: "#8b95a5", "font-size": 10.5 }, svg).textContent = v.toFixed(step < 1 ? 1 : 0); }
  svgEl("text", { x: 12, y: T + PH / 2, fill: "#8b95a5", "font-size": 10.5, transform: `rotate(-90 12 ${T + PH / 2})`, "text-anchor": "middle" }, svg).textContent = vioYLabel();
  const ids = [...new Set(subs)];
  p.groups.forEach((g, gi) => {
    // alternating bands per disease subtype, with the subtype name on top
    ids.forEach((sb, k) => { const idx = subs.map((s, i) => (s === sb ? i : -1)).filter((i) => i >= 0), x0 = pos[gi][idx[0]] - vw / 2 - 1, x1 = pos[gi][idx[idx.length - 1]] + vw / 2 + 1;
      if (k % 2 === 0) svgEl("rect", { x: x0, y: T, width: x1 - x0, height: PH, fill: "rgba(255,255,255,0.035)" }, svg); if (many || gi === 0) svgEl("text", { x: (x0 + x1) / 2, y: T - 5, "text-anchor": "middle", fill: "#6b7785", "font-size": 9.5 }, svg).textContent = x1 - x0 < 46 ? sb.replace("Class ", "") : sb; });
    const col = p.merged ? (g.id < 0 ? "#8b95a5" : catColor("cell_type", g.id)) : many ? catColor("cell_state", g.id) : catColor("cell_type", p.cell_type_id < 0 ? 0 : p.cell_type_id), colr = p.cell_type_id < 0 && !many ? "#8b95a5" : col;
    g.samples.forEach((r0) => { const i = order.indexOf(r0.sample), cx = pos[gi][i], m = meta[r0.sample];
      const tip = (ev) => nbrTip(ev, [{ tag: "b", t: r0.median != null ? "median " + r0.median.toFixed(2) : "too few cells" }, { t: `${dlabel(r0.sample)} · ${m.subtype}` }, { cls: "m", t: `${many ? g.name : p.cell_type} · ${fmt(r0.n)} cells${r0.n_used !== r0.n ? " (" + fmt(r0.n_used) + " used)" : ""} · ${r0.frac_nonzero == null ? "–" : (100 * r0.frac_nonzero).toFixed(1) + "% expressing"}` }, ...(r0.mean != null ? [{ cls: "m", t: `mean ${r0.mean.toFixed(2)} · IQR ${r0.q25.toFixed(2)} to ${r0.q75.toFixed(2)}` }] : [])]);
      if (r0.density) {
        const n = r0.density.length, pts = []; for (let j = 0; j < n; j++) pts.push([cx - r0.density[j] * hw, ys(p.lo + ((p.hi - p.lo) * j) / (n - 1))]); for (let j = n - 1; j >= 0; j--) pts.push([cx + r0.density[j] * hw, ys(p.lo + ((p.hi - p.lo) * j) / (n - 1))]);
        svgEl("path", { d: "M" + pts.map((q) => q[0].toFixed(1) + "," + q[1].toFixed(1)).join("L") + "Z", fill: colr, "fill-opacity": 0.85, stroke: "none" }, svg);
        svgEl("line", { x1: cx, x2: cx, y1: ys(r0.q25), y2: ys(r0.q75), stroke: "rgba(15,18,24,0.9)", "stroke-width": many ? 2 : 3 }, svg); svgEl("circle", { cx, cy: ys(r0.median), r: many ? 1.8 : 2.6, fill: "#ffffff" }, svg);
      } else svgEl("line", { x1: cx - hw, x2: cx + hw, y1: ys(p.lo), y2: ys(p.lo), stroke: "#4b566a", "stroke-width": 2 }, svg);   // too few cells: a grey tick
      const hit = svgEl("rect", { x: cx - vw / 2, y: T, width: vw, height: PH, fill: "transparent" }, svg); hit.addEventListener("pointermove", tip); hit.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
      if (!many) { const t = svgEl("text", { x: cx, y: T + PH + 12, fill: "#8b95a5", "font-size": 9, transform: `rotate(55 ${cx} ${T + PH + 12})` }, svg); t.textContent = dlabel(r0.sample); }
    });
    if (many) svgEl("text", { x: (pos[gi][0] + pos[gi][pos[gi].length - 1]) / 2, y: T + PH + 16, "text-anchor": "middle", fill: "#e3e8ef", "font-size": 11.5, "font-weight": 600 }, svg).textContent = p.merged ? g.name.replace(/ cells$/, "") : g.name.split(".")[0];
  });
  const w = mk("div", "vpanel"); w.append(mk("div", "vtitle", `${p.feature} — ${p.cell_type}${many && !p.merged ? " (cell states)" : ""}`));
  if (many) { const leg = mk("div", "vleg"); p.groups.forEach((g) => { const sp = mk("span"), k = mk("span", "k"); k.style.background = catColor(p.merged ? "cell_type" : "cell_state", g.id); sp.append(k, g.name); leg.append(sp); }); w.append(leg); }
  w.append(svg); return w;
}

// ---------------------------------------------------------------- glom heterogeneity
// One point per passed-qc glomerulus (data/embeddings from gloms_for_claude_*.qs, see scripts/export_gloms.R + build_gloms.py).
// Main-window filters apply as in the Glomeruli window: a glom counts as passing when >= 1 of its cells passes. Heatmaps use passing gloms only.
const GH = { info: null, nPass: null, fkey: "", colCache: {}, sel: -1, tab: "gene", tok: 0, corr: {}, leg: null, drawn: null };
const ghPass = (i) => !GH.nPass || GH.nPass[i] > 0;
async function openGH() {
  closeAllWindows();
  if (!GH.info) {
    try { GH.info = await jget("/api/glomhet/info"); } catch (e) { alert("Glom data are not available: " + e.message + "\nRun scripts/export_gloms.R and scripts/build_gloms.py first."); return; }
    const opt = (v, t) => { const o = document.createElement("option"); o.value = v; o.textContent = t; return o; };
    for (const id of ["ghpx", "ghpy"]) for (let k = 1; k <= GH.info.n_pcs; k++) $(id).append(opt(k - 1, "PC" + k));
    $("ghpx").value = 0; $("ghpy").value = 1;
    for (let k = 1; k <= GH.info.n_pcs; k++) $("ghsort").append(opt(k - 1, "PC" + k));
    GH.info.type_levels.forEach((t, i) => $("ghct").append(opt(i, t))); GH.info.state_levels.forEach((t, i) => $("ghcs").append(opt(i, t)));
    GH.info.cluster_cols.forEach((c, i) => $("ghclr").append(opt(c, c.replace("cluster_res", "")))); $("ghclr").value = GH.info.cluster_cols[Math.min(2, GH.info.cluster_cols.length - 1)];
  }
  try { await loadSigList(); } catch (e) { /* optional */ }
  for (const [id, def] of [["ghsig", "ifn1"], ["ghsig2", "ifn2"]]) { const sel = $(id); if (!sel.options.length && SIG.list) { for (const g of SIG_GROUPS) { const og = document.createElement("optgroup"); og.label = g;
    SIG.list.filter((x) => x.group === g).forEach((x) => { const o = document.createElement("option"); o.value = x.id; o.textContent = `${x.name} (${x.genes.length})`; og.append(o); }); sel.append(og); } if (SIG.byId[def]) sel.value = def; } }
  syncTransformUI();
  $("ghet").hidden = false; renderGH();
}
function closeGH() { $("ghet").hidden = true; $("nbrtip").hidden = true; }
function ghCtlVisibility() {
  const m = $("ghcol").value, gene = m === "gene" || m === "dual", sig = m === "signature" || m === "sigdual";
  $("ghgenel").hidden = !gene; $("ghgene2l").hidden = m !== "dual"; $("ghsigl").hidden = !sig; $("ghsig2l").hidden = m !== "sigdual"; $("ghsiglab").textContent = m === "sigdual" ? "A (blue)" : "Signature";
  $("ghctl_").hidden = m !== "type"; $("ghcsl").hidden = m !== "state"; $("ghclrl").hidden = m !== "cluster";
  $("ghexp").hidden = !(gene || sig); $("ghsrcl").hidden = m === "dual" || m === "sigdual";
  if (gene) { if (S.geneA && document.activeElement !== $("ghgene")) $("ghgene").value = S.geneA; if (S.geneB && document.activeElement !== $("ghgene2")) $("ghgene2").value = S.geneB; }
}
async function ghFilter() {
  const spec = nbrIndexSpec(), key = JSON.stringify(spec);
  if (key !== GH.fkey) { GH.nPass = (await postJSON("/api/glomhet/filter", spec)).n_pass; GH.fkey = key; GH.corr = {}; }
  const hasF = Object.keys(spec.filters).length || Object.keys(spec.hidden).length || spec.expr.length;
  const n = GH.info.gloms.filter((_, i) => ghPass(i)).length;
  return { n, hasF };
}
// colour of each glom: {cat:{levels,codes,colors}} or {num:{vals,lo,hi,label}}
async function ghColor() {
  const g = GH.info.gloms, mode = $("ghcol").value, msg = $("ghmsg"); msg.textContent = "";
  const catOf = (key, levels, colorOf) => { const idx = new Map(levels.map((l, i) => [String(l), i])); return { cat: { levels: levels.map(String), codes: g.map((x) => idx.get(String(x[key]))), colors: levels.map((l, i) => colorOf(l, i)) } }; };
  const present = (key) => [...new Set(g.map((x) => x[key]))];
  if (mode === "disease_subtype") return catOf("disease_subtype", present("disease_subtype").sort((a, b) => subtypeRank(a) - subtypeRank(b)), (l) => PAL[Math.max(0, SUB_ORDER.indexOf(l)) % PAL.length]);
  if (mode === "disease_type") return catOf("disease_type", present("disease_type").sort((a, b) => NB_TYPE_ORDER.indexOf(a) - NB_TYPE_ORDER.indexOf(b)), (l, i) => PAL[i]);
  if (mode === "sample") { const sub = Object.fromEntries(g.map((x) => [x.sample, x.disease_subtype])), r = catOf("sample", present("sample").sort((a, b) => subtypeRank(sub[a]) - subtypeRank(sub[b]) || a.localeCompare(b)), (l, i) => PAL[i % PAL.length]); r.cat.levels = r.cat.levels.map(dlabel); return r; }
  if (mode === "activity" || mode === "chronicity") { const lv = present(mode).sort((a, b) => a - b); return catOf(mode, lv, (l, i) => rgbStr(LUT[Math.round((lv.length > 1 ? i / (lv.length - 1) : 1) * 255)])); }
  if (mode === "cluster") { const k = $("ghclr").value; return catOf(k, present(k).sort((a, b) => a - b), (l, i) => PAL[i % PAL.length]); }
  if (mode === "gene" || mode === "dual" || mode === "signature" || mode === "sigdual") return ghExprColor(mode, msg);
  let body, label;
  { const lv = +$(mode === "type" ? "ghct" : "ghcs").value; body = { kind: mode, level: lv }; label = (mode === "type" ? GH.info.type_levels : GH.info.state_levels)[lv] + " (fraction of the glom's cells)"; }
  const key = JSON.stringify(body);
  if (!GH.colCache[key]) { try { GH.colCache[key] = await postJSON("/api/glomhet/color", body); } catch (e) { msg.textContent = e.message.replace(/^.*?: /, "").slice(0, 160); return { none: true }; } }
  const n = ghNum(GH.colCache[key].values, label); return n ? { num: n } : { none: true };
}
// quantile-limited numeric colouring of the gloms (null / NaN = no value); z-scores get a diverging scale centred on 0
function ghNum(values, label, div = false) {
  const vals = Float32Array.from(values, (v) => (v == null ? NaN : v)), fin = vals.filter((v) => v === v), sorted = Float32Array.from(fin).sort(), q = (p) => sorted[Math.min(sorted.length - 1, Math.floor(p * (sorted.length - 1)))];
  if (!fin.length) return null;
  let lo = q(0.02), hi = q(0.98);
  if (div) { const m = Math.max(Math.abs(lo), Math.abs(hi), 1e-6); lo = -m; hi = m; } else if (hi - lo < 1e-9) { lo = sorted[0]; hi = Math.max(sorted[sorted.length - 1], lo + 1e-9); }
  return { vals, lo, hi, label, div, nMissing: vals.length - fin.length };
}
// gene / signature colouring (the main window's settings): from the glom object, or from the viewer cells with layer / smoothing / z-score / compute-over
async function ghExprColor(mode, msg) {
  const dual = mode === "dual" || mode === "sigdual", sigm = mode === "signature" || mode === "sigdual";
  const chan = sigm ? [[S.sigGenes, S.sigName || `Signature (${S.sigGenes.length} genes)`], [S.sigGenes2, S.sigName2 || `Signature B (${S.sigGenes2.length} genes)`]] : [[S.geneA ? [S.geneA] : [], S.geneA], [S.geneB ? [S.geneB] : [], S.geneB]];
  if (!chan[0][0].length || (dual && !chan[1][0].length)) { msg.textContent = sigm ? "choose a signature (or paste genes) for each channel" : "type a gene (and a second gene for dual colour)"; return { none: true }; }
  const cells = dual || transformed() || S.layer === "counts" || $("ghsrc").value === "cells";
  $("ghsrc").value = cells ? "cells" : "object";
  const get = async (genes, label) => {
    if (!cells) {   // glom-level object, as before
      const body = sigm ? { kind: "signature", genes } : { kind: "gene", gene: genes[0] }, key = JSON.stringify(body);
      if (!GH.colCache[key]) GH.colCache[key] = await postJSON("/api/glomhet/color", body);
      const r = GH.colCache[key]; if (r.genes_used != null && r.genes_used < r.genes_requested) msg.textContent = `${r.genes_used} of ${r.genes_requested} signature genes are in the region_heterogeneity object`;
      return ghNum(r.values, label + (sigm ? " (mean log-normalised, region_heterogeneity object)" : " (ROI-level log-normalised)"));
    }
    const body = { genes, layer: S.layer, sigma: S.smooth, scope: S.zscope, within: withinSpec() }, key = "cs|" + JSON.stringify(body);
    if (!GH.colCache[key]) GH.colCache[key] = await postJSON("/api/glomhet/cellscore", body);
    const r = GH.colCache[key], bits = [sigm ? "mean" : "", S.zscope !== "none" ? "z-score" : "", S.layer === "counts" ? "raw counts" : "log-normalised"].filter(Boolean).join(" ");
    return ghNum(r.values, `${label} (${bits} of the glom's cells${S.smooth ? `, σ ${S.smooth} µm` : ""}${withinSpec() !== "tissue" ? ", computed over " + (withinSpec() === "glom" ? "glomerulus cells" : "region " + rgName()) : ""})`, S.zscope !== "none");
  };
  try {
    const A = await get(chan[0][0], chan[0][1]); if (!A) { msg.textContent = "no glom has a value"; return { none: true }; }
    if (!dual) { if (A.nMissing) msg.textContent = `${A.nMissing} gloms have no cell inside the compute-over subset (dark grey)`; return { num: A }; }
    const B = await get(chan[1][0], chan[1][1]); if (!B) { msg.textContent = "no glom has a value"; return { none: true }; }
    if (!sigm || S.zscope !== "none") { A.lo = 0; B.lo = 0; for (const n of [A, B]) n.hi = Math.max(n.hi, 1e-6); }   // as in the main window: only the positive part is coloured (zero / below average = neutral grey)
    return { dual: { A, B, sig: sigm } };
  } catch (e) { msg.textContent = e.message.replace(/^.*?: /, "").slice(0, 160); return { none: true }; }
}
const ghT = (n, i) => Math.max(0, Math.min(1, (n.vals[i] - n.lo) / (n.hi - n.lo)));
const ghLut = (n) => (n.div ? DIVL : LUT);
function ghFill(c, i) {
  if (c.cat) return c.cat.colors[c.cat.codes[i]] || "#8b95a5";
  if (c.dual) {   // genes: A red, B green, both yellow, neither grey (as in the main window); signatures: A blue, B orange, both near white
    const { A, B, sig } = c.dual; if (A.vals[i] !== A.vals[i] || B.vals[i] !== B.vals[i]) return "#4a5260";
    const tA = ghT(A, i), tB = ghT(B, i);
    if (sig) { const cA = hex2rgb(SIG_A), cB = hex2rgb(SIG_B), m = Math.max(tA, tB); return rgbStr([0, 1, 2].map((q) => Math.min(255, Math.round(GRAY[q] * (1 - m) + tA * cA[q] + tB * cB[q])))); }
    let r = Math.round(255 * tA), g = Math.round(255 * tB); if (r < 40 && g < 40) { r = Math.max(r, 60); g = Math.max(g, 60); } return `rgb(${r},${g},0)`;
  }
  if (c.num.vals[i] !== c.num.vals[i]) return "#4a5260";   // no cell of this glom inside the compute-over subset
  return rgbStr(ghLut(c.num)[Math.round(ghT(c.num, i) * 255)]);
}
const ghRank = (c, i) => (c.cat ? 0 : c.dual ? Math.max(ghT(c.dual.A, i), ghT(c.dual.B, i)) : c.num.vals[i] === c.num.vals[i] ? c.num.vals[i] : -Infinity);
function ghScatter(svg, X, Y, col, xl, yl) {
  const W = 360, H = 340, L = 40, R = 10, T = 10, B = 34; svg.replaceChildren(); svg.setAttribute("viewBox", `0 0 ${W} ${H}`); svg.setAttribute("width", W); svg.setAttribute("height", H);
  const g = GH.info.gloms, hide = $("ghhide").checked, idx = g.map((_, i) => i).filter((i) => !hide || ghPass(i));
  const ext = (a) => { let lo = Infinity, hi = -Infinity; for (const i of idx.length ? idx : g.map((_, k) => k)) { lo = Math.min(lo, a[i]); hi = Math.max(hi, a[i]); } const p = (hi - lo) * 0.05 || 1; return [lo - p, hi + p]; };
  const [x0, x1] = ext(X), [y0, y1] = ext(Y), sx = (v) => L + ((v - x0) / (x1 - x0)) * (W - L - R), sy = (v) => T + (1 - (v - y0) / (y1 - y0)) * (H - T - B);
  svgEl("rect", { x: L, y: T, width: W - L - R, height: H - T - B, fill: "none", stroke: "#29303b" }, svg);
  const ticks = (lo, hi) => { const r = (hi - lo) / 4, q = Math.pow(10, Math.floor(Math.log10(r))), st = [1, 2, 5, 10].map((k) => k * q).find((k) => k >= r); const o = []; for (let v = Math.ceil(lo / st) * st; v <= hi; v += st) o.push(v); return o; };
  for (const v of ticks(x0, x1)) { svgEl("text", { x: sx(v), y: H - B + 13, "text-anchor": "middle", fill: "#6b7785", "font-size": 9.5, "font-family": "system-ui,sans-serif" }, svg).textContent = +v.toFixed(2); }
  for (const v of ticks(y0, y1)) { svgEl("text", { x: L - 4, y: sy(v) + 3, "text-anchor": "end", fill: "#6b7785", "font-size": 9.5, "font-family": "system-ui,sans-serif" }, svg).textContent = +v.toFixed(2); }
  svgEl("text", { x: (L + W - R) / 2, y: H - 6, "text-anchor": "middle", fill: "#8b95a5", "font-size": 11, "font-family": "system-ui,sans-serif" }, svg).textContent = xl;
  svgEl("text", { x: 11, y: (T + H - B) / 2, "text-anchor": "middle", fill: "#8b95a5", "font-size": 11, "font-family": "system-ui,sans-serif", transform: `rotate(-90 11 ${(T + H - B) / 2})` }, svg).textContent = yl;
  // draw order: failing filters first (grey), then passing gloms low -> high value so the high ones stay visible
  const order = idx.slice().sort((a, b) => (ghPass(a) - ghPass(b)) || (col.none ? 0 : ghRank(col, a) - ghRank(col, b) || 0));
  for (const i of order) {
    const pass = ghPass(i), c = svgEl("circle", { cx: sx(X[i]), cy: sy(Y[i]), r: pass ? +$("ghsize").value : Math.max(1.5, +$("ghsize").value * 0.7), fill: col.none ? "#8b95a5" : pass ? ghFill(col, i) : "#3a4250", "fill-opacity": pass ? +$("ghop").value : Math.min(0.55, +$("ghop").value), stroke: "rgba(10,12,16,0.9)", "stroke-width": 0.6, style: "cursor:pointer" }, svg);
    c.addEventListener("click", () => ghSelect(i));
    c.addEventListener("pointermove", (ev) => { const m = g[i]; nbrTip(ev, [{ tag: "b", t: m.id }, { t: [dlabel(m.sample), m.disease_type, m.disease_subtype].filter((x) => x != null && x !== "").join(" · ") }, { cls: "m", t: `${m.activity != null ? "activity " + m.activity + " · " : ""}${m.chronicity != null ? "chronicity " + m.chronicity + " · " : ""}${m.n_cells} cells` },
      ...(col.num ? [{ cls: "m", t: `${col.num.label.split(" (")[0]}: ${col.num.vals[i] === col.num.vals[i] ? col.num.vals[i].toFixed(3) : "no cell inside the subset"}` }] : []),
      ...(col.dual ? [col.dual.A, col.dual.B].map((n) => ({ cls: "m", t: `${n.label.split(" (")[0]}: ${n.vals[i] === n.vals[i] ? n.vals[i].toFixed(3) : "n/a"}` })) : []), ...(pass ? [] : [{ cls: "m", t: "does not pass the main-window filters" }])]); });
    c.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
  }
  if (GH.sel >= 0 && (!hide || ghPass(GH.sel))) svgEl("circle", { cx: sx(X[GH.sel]), cy: sy(Y[GH.sel]), r: +$("ghsize").value + 4.1, fill: "none", stroke: "#fff", "stroke-width": 2, "pointer-events": "none" }, svg);
}
function ghLegend(col) {
  const box = $("ghleg"); box.replaceChildren(); GH.leg = null;
  if (col.none) return;
  if (col.cat) { box.append(mk("b", "", $("ghcol").selectedOptions[0].text + ":")); col.cat.levels.forEach((l, i) => { const sp = mk("span"), k = mk("span", "k"); k.style.background = col.cat.colors[i]; sp.append(k, l); if (col.cat.levels.length <= 40) box.append(sp); });
    if (col.cat.levels.length > 40) box.append(mk("span", "", col.cat.levels.length + " groups (hover a dot)")); GH.leg = { title: $("ghcol").selectedOptions[0].text, items: col.cat.levels.length <= 40 ? col.cat.levels.map((l, i) => ({ l, c: col.cat.colors[i] })) : null }; return; }
  if (col.dual) {
    const { A, B, sig } = col.dual, cA = sig ? SIG_A : "#d33", cB = sig ? SIG_B : "#3a3";
    for (const [c, n] of [[cA, A], [cB, B]]) { const sp = mk("span"), k = mk("span", "k"); k.style.background = c; sp.append(k, `${n.label} (${n.lo.toFixed(2)} → ${n.hi.toFixed(2)}, 2nd–98th percentile)`); box.append(sp); }
    box.append(mk("span", "", sig ? "near white = both high · grey = low" : "yellow = both · grey = neither"));
    GH.leg = { title: sig ? "Two signatures" : "Two genes", items: [{ l: `${A.label.split(" (")[0]}`, c: cA }, { l: `${B.label.split(" (")[0]}`, c: cB }] }; return;
  }
  const cv = mk("canvas"); cv.width = 160; cv.height = 10; const ctx = cv.getContext("2d"), L = ghLut(col.num); for (let x = 0; x < 160; x++) { const c = L[Math.round((x / 159) * 255)]; ctx.fillStyle = rgbStr(c); ctx.fillRect(x, 0, 1, 10); }
  cv.style.cssText = "width:100%;height:10px;display:block"; const wrap = mk("div"); wrap.style.cssText = "display:inline-block;width:220px;vertical-align:middle;margin:0 8px"; const tkr = mk("div"); tkr.innerHTML = tickRowHTML(scaleTicks({ type: "num", lo: col.num.lo, hi: col.num.hi })); wrap.append(cv, tkr);
  box.append(mk("b", "", col.num.label + ":"), wrap, mk("span", "", col.num.div ? "z-score, symmetric around 0; dark grey = no value" : "ends: 2nd and 98th percentile of gloms")); GH.leg = { title: col.num.label, num: col.num };
}
async function renderGH() {
  if ($("ghet").hidden) return; const tok = ++GH.tok; ghCtlVisibility();
  let f, col; try { f = await ghFilter(); col = await ghColor(); } catch (e) { $("ghmsg").textContent = "Error: " + e.message.slice(0, 160); return; }
  if (tok !== GH.tok) return;
  const I = GH.info, px = +$("ghpx").value, py = +$("ghpy").value;
  $("ghsub").textContent = `${I.n} ROIs in the region_heterogeneity object · ${f.n} pass the main-window filters${f.hasF ? " (others grey)" : ""}`;
  const hasH = !!I.pca.harmony;   // the batch-corrected PCA is optional in the object
  $("ghpca").closest(".ghp").hidden = !hasH; $("ghpcaraw").closest(".ghp").querySelector(".ghpt").textContent = hasH ? "PCA (uncorrected)" : "PCA";
  $("ghlayer").querySelector("option[value=harmony]").hidden = !hasH; if (!hasH) $("ghlayer").value = "raw";
  for (const [id, layer, lab] of [["ghpca", "harmony", "Harmony PC"], ["ghpcaraw", "raw", hasH ? "PC" : "PC"]]) { const P = I.pca[layer]; if (!P) continue; ghScatter($(id), P.map((r) => r[px]), P.map((r) => r[py]), col, lab + (px + 1), lab + (py + 1)); }   // both PCA versions, same x / y PCs and colours
  const U = I.umap; ghScatter($("ghumap"), U.map((r) => r[0]), U.map((r) => r[1]), col, "UMAP1", "UMAP2");
  GH.col = col; ghLegend(col);
  if (GH.sel >= 0) ghShow(GH.sel);
  renderGHHeat();
}
function ghSelect(i) { GH.sel = i; renderGH(); }
async function ghShow(i) {
  const m = GH.info.gloms[i], tok = (GH.dtok = (GH.dtok || 0) + 1), smp = S.meta.samples.find((s) => s.id === m.sample), cv = $("ghcv");
  $("ghdett").textContent = m.id; $("ghdi").innerHTML = `${dlabel(m.sample)}<br>${[m.disease_type, m.disease_subtype].filter((x) => x != null && x !== "").join(" · ")}${m.activity != null || m.chronicity != null ? `<br>${m.activity != null ? "activity " + m.activity : ""}${m.activity != null && m.chronicity != null ? " · " : ""}${m.chronicity != null ? "chronicity " + m.chronicity : ""}` : ""}<br>${m.n_cells} cells in the viewer${m.n_file != null ? ` (${m.n_file} in the object)` : ""}${GH.info.cluster_cols.length ? `<br>clusters (res ${GH.info.cluster_cols.map((c) => c.replace("cluster_res", "") + ": " + m[c]).join(", ")})` : ""}<br><span class="muted">cells coloured by the main window's Colour by · hover a cell for its type and state</span>`;
  const gidLv = catMeta("glom_id").levels, code = gidLv.indexOf(m.id);
  const [gid, poly, feats, b] = await Promise.all([galGet(m.sample, "glom_id"), galPoly(m.sample), galPaths(m.sample), winBundle(m.sample)]);
  const sc = await winScale([b], () => tok, tok); if (!sc || tok !== GH.dtok) return;
  const cells = []; for (let k = 0; k < b.n; k++) if (gid[k] === code) cells.push(k);
  const outline = glomOutline(feats, m.id);
  const [ct, cs] = await Promise.all([galGet(m.sample, "cell_type"), galGet(m.sample, "cell_state")]); if (tok !== GH.dtok) return;
  cv.hidden = false; GH.det = { b, sc, poly, cells, outline, ct, cs, hover: -1, views: {} }; ghDraw();
  const ad = $("ghadd"); ad.hidden = false; ad.onclick = () => glomTick(`${m.sample}|${m.id}`, { sample: m.sample, label: m.id }, true);
  const z = $("ghzoom"); z.hidden = false; z.onclick = () => openGlomZoom({ src: "gh", s: smp, b, sc, poly, cells, label: m.id, qc: "passed qc", outline, note: outline ? "outline = pathologist \"All gloms\" boundary" : undefined });
}
// hover a cell of the glom drawn in the detail panel: its cell type and state (the cell is outlined in white)
function ghDraw() {   // the glom detail: cells, plus (H&E beside) a second canvas with the H&E only; both can be zoomed and panned
  const d = GH.det, cv = $("ghcv"), host = $("ghdet"); if (!d || cv.hidden) return;
  cv._hw = "gh"; const cv2 = hePane(host, cv, "gh", d.b.sample), B = glomBox(d.cells, d.poly, 12, d.outline);
  host.dataset.id = "gh"; if (!host.querySelector(".zb")) host.append(mk("span", "zb"));
  const rec = host._dt = { bb: { x0: B.x0, y0: B.y0, w: B.x1 - B.x0, h: B.y1 - B.y0 }, draw: (vw, vh) => { drawGlom(cv, d.b, d.sc, d.poly, d.cells, true, d.hover, d.outline, null, vw); if (cv2) drawGlom(cv2, d.b, d.sc, d.poly, d.cells, true, -1, d.outline, null, HEW.gh.link ? vw : vh); } };
  if (cv2) { cv2.width = cv2.height = cv.width; attachZoom(host, rec, { views: () => d.views, sync: () => false, tiles: () => [host] }); }
  rec.draw(d.views.gh || null, d.views["gh#h"] || null);
}
function ghHover(ev) {
  const d = GH.det, cv = $("ghcv"); if (!d || cv.hidden || !cv._xf) return;
  const r = cv.getBoundingClientRect(), px = ((ev.clientX - r.left) / r.width) * cv.width, py = ((ev.clientY - r.top) / r.height) * cv.height;
  const x = (px - cv._xf.ox) / cv._xf.s, y = (py - cv._xf.oy) / cv._xf.s;
  const inside = (i) => { const o = d.poly.off; let c = false; for (let a = o[i], q = o[i + 1] - 1; a < o[i + 1]; q = a++) { const x1 = d.poly.xy[a * 2], y1 = d.poly.xy[a * 2 + 1], x2 = d.poly.xy[q * 2], y2 = d.poly.xy[q * 2 + 1];
    if (y1 > y !== y2 > y && x < ((x2 - x1) * (y - y1)) / (y2 - y1) + x1) c = !c; } return c; };
  let hit = -1, best = 1e9;
  for (const i of d.cells) { if (inside(i)) { hit = i; break; } const dd = Math.hypot(d.b.xy[2 * i] - x, d.b.xy[2 * i + 1] - y); if (dd < best && dd < 5) { best = dd; hit = i; } }
  if (hit < 0) { $("nbrtip").hidden = true; if (d.hover >= 0) { d.hover = -1; ghDraw(); } return; }
  if (hit !== d.hover) { d.hover = hit; ghDraw(); }
  const tl = catMeta("cell_type").levels, sl = catMeta("cell_state").levels;
  nbrTip(ev, [{ tag: "b", t: d.ct[hit] === MISSING ? "no cell type" : tl[d.ct[hit]] }, { t: d.cs[hit] === MISSING ? "no cell state" : sl[d.cs[hit]] }, ...(d.b.pass[hit] ? [] : [{ cls: "m", t: "does not pass the main-window filters" }])]);
}
function ghHoverEnd() { const d = GH.det; $("nbrtip").hidden = true; if (d && d.hover >= 0) { d.hover = -1; ghDraw(); } }
// ---- heatmaps: Spearman correlation of every gene / cell type / cell state (rows) with PC1..PC10 across the shown gloms
// Rows of a correlation result ordered for display: sortPc = PC index (0-based) or -1 for the default order; dir = desc | asc | abs
function ghSortedRows(d, kind, sortPc, dir) {
  const ix = d.rows.map((_, i) => i);
  if (sortPc >= 0) { const v = (i) => d.r[i][sortPc]; ix.sort((x, y) => (dir === "asc" ? v(x) - v(y) : dir === "abs" ? Math.abs(v(y)) - Math.abs(v(x)) : v(y) - v(x)) || d.rows[x].localeCompare(d.rows[y])); }
  return { ...d, rows: ix.map((i) => d.rows[i]), r: ix.map((i) => d.r[i]), parents: d.parents ? ix.map((i) => d.parents[i]) : undefined, sorted: sortPc >= 0, sortPc };
}
function ghSetSort(j) {   // click on a PC header: first click sorts by it (most positive first), then flips to most negative, then back to the default order
  const cur = $("ghsort").value, dir = $("ghdir").value;
  if (cur !== String(j)) { $("ghsort").value = String(j); $("ghdir").value = "desc"; }
  else if (dir === "desc") $("ghdir").value = "asc"; else if (dir === "asc") $("ghdir").value = "abs"; else $("ghsort").value = "";
  renderGHHeat();
}
async function renderGHHeat() {
  const box = $("ghhmbox"), kind = GH.tab, layer = $("ghlayer").value, tok = GH.tok, view = $("ghview").value, sortPc = $("ghsort").value === "" ? -1 : +$("ghsort").value, dir = $("ghdir").value;
  const ntop = Math.max(1, Math.min(100, +$("ghtop").value || 10));
  for (const k of ["gene", "type", "state"]) $("ghtab-" + k).classList.toggle("on", k === kind);
  $("ghtop").parentElement.style.display = kind === "gene" ? "" : "none"; $("ghdir").parentElement.style.display = sortPc >= 0 ? "" : "none";
  $("ghtoplab").textContent = sortPc >= 0 ? `genes (top + bottom) for PC${sortPc + 1}` : "genes per PC (top + bottom)";
  const idx = GH.info.gloms.map((_, i) => i).filter(ghPass), spc = kind === "gene" && sortPc >= 0 ? sortPc : null, key = [kind, layer, ntop, spc, GH.fkey, idx.length].join("|");
  if (!GH.corr[key]) { $("ghhmnote").textContent = "Computing…"; try { GH.corr[key] = await postJSON("/api/glomhet/corr", { kind, layer, glom_idx: idx, n_top: ntop, sort_pc: spc }); } catch (e) { if (tok === GH.tok) { box.replaceChildren(); $("ghhmnote").textContent = e.message.replace(/^.*?: /, "").slice(0, 160); } return; } }
  if (tok !== GH.tok) return;
  const d0 = GH.corr[key], d = ghSortedRows(d0, kind, sortPc, dir); GH.hm = { d, kind };
  const dirTxt = { desc: "most positive first", asc: "most negative first", abs: "strongest |r| first" }[dir];
  $("ghhmnote").textContent = `Spearman correlation across ${d.n_gloms} gloms (${layer === "harmony" ? "Harmony-corrected" : "uncorrected"} PCs)` + (kind === "gene" ? (sortPc >= 0 ? `; ${d.rows.length} genes = the ${ntop} most positive and ${ntop} most negative for PC${sortPc + 1}` : `; ${d.rows.length} genes = the ${ntop} most positive and ${ntop} most negative per PC, unique`) : "")
    + (sortPc >= 0 ? `; sorted by PC${sortPc + 1}, ${dirTxt}` : "") + `; click a row to colour the plots by it, click a PC header to sort by it`;
  box.replaceChildren(view === "table" ? ghTable(d, kind) : ghHeatSVG(d, kind));
}
function ghTable(d, kind) {
  const t = mk("table", "st ghtab"), hr = t.createTHead().insertRow(), m = Math.max(0.3, Math.ceil(Math.max(...d.r.flat().map(Math.abs)) * 10) / 10);
  const th = (txt, cls) => { const e = mk("th", cls || "", txt); hr.append(e); return e; };
  if (kind === "state") th("Cell type"); th(kind === "gene" ? "Gene" : kind === "type" ? "Cell type" : "Cell state");
  for (let j = 0; j < d.npc; j++) { const e = th(`PC${j + 1}${d.sortPc === j ? ({ desc: " ▼", asc: " ▲", abs: " |r|" })[$("ghdir").value] : ""}`, d.sortPc === j ? "sorted" : ""); e.style.cursor = "pointer"; e.title = `sort by PC${j + 1}`; e.addEventListener("click", () => ghSetSort(j)); }
  const tb = t.createTBody();
  d.rows.forEach((name, i) => {
    const tr = tb.insertRow(); if (kind === "state") { const cp = tr.insertCell(); cp.textContent = d.parents[i]; cp.className = "lab"; }
    const c0 = tr.insertCell(); c0.className = "lab"; c0.textContent = name; c0.style.cursor = "pointer"; c0.addEventListener("click", () => ghColourBy(kind, name));
    for (let j = 0; j < d.npc; j++) {
      const v = d.r[i][j], c = DIVL[Math.round(Math.max(0, Math.min(1, (v + m) / (2 * m))) * 255)], td = tr.insertCell(); td.textContent = v.toFixed(2); td.style.background = rgbStr(c); td.style.color = Math.abs(v) > m * 0.6 ? "#fff" : "#1a1f27";
      if (d.sortPc === j) td.style.outline = "1px solid #fff";
    }
  });
  return t;
}
function ghExportCSV() {
  if (!GH.hm) return; const { d, kind } = GH.hm, q = (x) => `"${String(x).replace(/"/g, '""')}"`;
  const head = [...(kind === "state" ? ["cell_type"] : []), kind === "gene" ? "gene" : kind === "type" ? "cell_type" : "cell_state", ...Array.from({ length: d.npc }, (_, j) => "PC" + (j + 1))];
  download(`glom_pc_correlations_${kind}${d.sorted ? "_sorted_PC" + (d.sortPc + 1) : ""}.csv`, [head.join(","), ...d.rows.map((n, i) => [...(kind === "state" ? [q(d.parents[i])] : []), q(n), ...d.r[i]].join(","))].join("\n"));
}
function ghHeatSVG(d, kind) {
  const cw = 46, ch = kind === "type" ? 22 : kind === "state" ? 15 : 12, Lw = kind === "gene" ? 74 : kind === "type" ? 190 : 250, T = 40, nR = d.rows.length, nC = d.npc, W = Lw + nC * cw + 120, H = T + nR * ch + 30;
  const svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` }), m = Math.max(0.3, Math.ceil(Math.max(...d.r.flat().map(Math.abs)) * 10) / 10), ff = "system-ui,sans-serif";
  for (let j = 0; j < nC; j++) { const h = svgEl("text", { x: Lw + j * cw + cw / 2, y: T - 8, "text-anchor": "middle", fill: d.sortPc === j ? "#fff" : "#aab4c3", "font-weight": d.sortPc === j ? 700 : 400, "font-size": 11, "font-family": ff, style: "cursor:pointer" }, svg);
    h.textContent = "PC" + (j + 1) + (d.sortPc === j ? ({ desc: "▼", asc: "▲", abs: "|r|" })[$("ghdir").value] : ""); h.addEventListener("click", () => ghSetSort(j)); }
  d.rows.forEach((name, i) => {
    const y = T + i * ch;
    if (kind === "state" && !d.sorted && i && d.parents[i] !== d.parents[i - 1]) svgEl("line", { x1: 0, x2: Lw + nC * cw, y1: y, y2: y, stroke: "#0f1218", "stroke-width": 3 }, svg);
    const t = svgEl("text", { x: Lw - 6, y: y + ch / 2 + 3.5, "text-anchor": "end", fill: "#d7dde6", "font-size": Math.min(11, ch - 1), "font-family": ff, style: "cursor:pointer" }, svg); t.textContent = name;
    t.addEventListener("click", () => ghColourBy(kind, name));
    for (let j = 0; j < nC; j++) {
      const v = d.r[i][j], c = DIVL[Math.round(Math.max(0, Math.min(1, (v + m) / (2 * m))) * 255)], r = svgEl("rect", { x: Lw + j * cw + 0.5, y: y + 0.5, width: cw - 1, height: ch - 1, fill: rgbStr(c), style: "cursor:pointer" }, svg);
      r.addEventListener("pointermove", (ev) => nbrTip(ev, [{ tag: "b", t: `${name} · PC${j + 1}` }, { t: `Spearman r = ${v.toFixed(3)}` }, { cls: "m", t: `${d.n_gloms} gloms` }])); r.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
      r.addEventListener("click", () => ghColourBy(kind, name));
      if (ch >= 20) svgEl("text", { x: Lw + j * cw + cw / 2, y: y + ch / 2 + 3.5, "text-anchor": "middle", fill: Math.abs(v) > m * 0.6 ? "#fff" : "#1a1f27", "font-size": 9.5, "font-family": ff, "pointer-events": "none" }, svg).textContent = v.toFixed(2);
    }
  });
  if (d.sortPc >= 0) svgEl("rect", { x: Lw + d.sortPc * cw, y: T, width: cw, height: nR * ch, fill: "none", stroke: "#fff", "stroke-width": 1.5, "pointer-events": "none" }, svg);
  const x0 = Lw + nC * cw + 18;   // colour bar
  for (let k = 0; k < 100; k++) svgEl("rect", { x: x0, y: T + k * 1.2, width: 12, height: 1.3, fill: rgbStr(DIVL[Math.round(((99 - k) / 99) * 255)]) }, svg);
  for (const [v, y] of [[m, T + 4], [0, T + 60], [-m, T + 118]]) svgEl("text", { x: x0 + 16, y, fill: "#8b95a5", "font-size": 10, "font-family": ff }, svg).textContent = (v > 0 ? "+" : "") + v.toFixed(1);
  svgEl("text", { x: x0, y: T - 8, fill: "#8b95a5", "font-size": 10, "font-family": ff }, svg).textContent = "Spearman r";
  return svg;
}
function ghColourBy(kind, name) {
  const I = GH.info;
  if (kind === "gene") { $("ghcol").value = "gene"; $("ghgene").value = name; return ghPush(); }
  else if (kind === "type") { $("ghcol").value = "type"; $("ghct").value = I.type_levels.indexOf(name); }
  else { $("ghcol").value = "state"; $("ghcs").value = I.state_levels.indexOf(name); }
  renderGH();
}
async function ghSavePNG() {
  const svgs = [$("ghpca"), $("ghpcaraw"), $("ghumap"), $("ghhmbox").querySelector("svg")].filter(Boolean), S2 = 2, pad = 16;
  const toImg = (svg) => new Promise((res, rej) => { const c = svg.cloneNode(true); c.setAttribute("xmlns", "http://www.w3.org/2000/svg"); const im = new Image();
    im.onload = () => res(im); im.onerror = rej; im.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(new XMLSerializer().serializeToString(c)); });
  const imgs = await Promise.all(svgs.map(toImg)), w0 = +svgs[0].getAttribute("width"), h0 = +svgs[0].getAttribute("height");
  const hm = imgs[3], hmW = hm ? +svgs[3].getAttribute("width") : 0, hmH = hm ? +svgs[3].getAttribute("height") : 0;
  const W = Math.max(3 * w0 + 4 * pad, hmW + 2 * pad), H = pad + 22 + 20 + h0 + 26 + (hm ? hmH + pad : 0) + pad, cv = document.createElement("canvas"); cv.width = W * S2; cv.height = H * S2;
  const ctx = cv.getContext("2d"); ctx.scale(S2, S2); ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, W, H); ctx.fillStyle = "#d7dde6"; ctx.font = "600 14px system-ui,sans-serif";
  ctx.fillText(`ROI heterogeneity — ${$("ghsub").textContent}`, pad, pad + 14);
  ctx.font = "600 12px system-ui,sans-serif"; ctx.fillStyle = "#aab4c3"; if (GH.info.pca.harmony) ctx.fillText("PCA (Harmony-corrected)", pad, pad + 40); ctx.fillText(GH.info.pca.harmony ? "PCA (uncorrected)" : "PCA", 2 * pad + w0, pad + 40); ctx.fillText("UMAP", 3 * pad + 2 * w0, pad + 40);
  ctx.drawImage(imgs[0], pad, pad + 46, w0, h0); ctx.drawImage(imgs[1], 2 * pad + w0, pad + 46, w0, h0); ctx.drawImage(imgs[2], 3 * pad + 2 * w0, pad + 46, w0, h0);
  let y = pad + 46 + h0 + 16; ctx.font = "11px system-ui,sans-serif"; ctx.fillStyle = "#aab4c3";
  const L = GH.leg; if (L) { let x = pad; ctx.fillText(L.title + ":", x, y); x += ctx.measureText(L.title + ":").width + 8;
    if (L.num) { for (let k = 0; k < 120; k++) { ctx.fillStyle = rgbStr(ghLut(L.num)[Math.round((k / 119) * 255)]); ctx.fillRect(x + k, y - 9, 1, 10); } ctx.fillStyle = "#aab4c3"; const tk3 = scaleTicks({ type: "num", lo: L.num.lo, hi: L.num.hi }, 3); tk3.forEach((k) => { ctx.textAlign = k.t === 0 ? "left" : k.t === 1 ? "right" : "center"; ctx.fillText(fmtTick(k.v), x + 119 * k.t + (k.t === 1 ? 1 : 0), y + 12); }); ctx.textAlign = "left"; }
    else if (L.items) for (const it of L.items) { ctx.fillStyle = it.c; ctx.fillRect(x, y - 9, 10, 10); ctx.fillStyle = "#aab4c3"; ctx.fillText(it.l, x + 14, y); x += 24 + ctx.measureText(it.l).width; } }
  if (hm) ctx.drawImage(hm, pad, y + 22, hmW, hmH);
  const a = document.createElement("a"); a.download = "roi_heterogeneity.png"; a.href = cv.toDataURL("image/png"); a.click();
}

// ---------------------------------------------------------------- Xenium vs scFFPE: cell fractions of the current sample
// Xenium = every QC'ed Xenium cell of the sample in the viewer; scFFPE = the cells of the linked scFFPE sample (counts only, no coordinates).
const CMP_C = { xen: "#3987e5", sc: "#d95926" };   // validated categorical slots 1 and 2
function openCmp() { closeAllWindows(); $("cmp").hidden = false; S.cmpOverride = ""; renderCmp(); }
function closeCmp() { $("cmp").hidden = true; }
async function renderCmp() {
  if ($("cmp").hidden) return;
  const level = $("cmplevel").value, incl = $("cmpincl").checked, tok = (S.cmpTok = (S.cmpTok || 0) + 1), body = $("cmpbody"), info = $("cmpinfo");
  info.textContent = "Loading…";
  let d; try { d = await jget(`/api/compare/${encodeURIComponent(S.sample)}?level=${level}&include_excluded=${incl}${S.cmpOverride ? "&scffpe=" + encodeURIComponent(S.cmpOverride) : ""}`); } catch (e) { info.textContent = "Error: " + e.message; return; }
  if (tok !== S.cmpTok) return;
  const smp = S.meta.samples.find((s) => s.id === S.sample), el = (t, c, x) => { const e = document.createElement(t); if (c) e.className = c; if (x != null) e.textContent = x; return e; };
  body.replaceChildren();
  info.textContent = `Xenium donor ${dlabel(d.sample)} (${smp.disease_type} · ${subtypeOf(smp)}): ${fmt(d.xenium.n)} QC'ed Xenium cells.`;
  if (!d.scffpe) {   // no counterpart: say so, and why
    const box = el("div", "cmpmsg");
    const msg = d.status === "id_without_cells" ? `The metadata links this sample to scFFPE sample ${d.scffpe_id}, but no scFFPE cells carry that sample ID in the data (it is absent from the scFFPE cells).`
      : d.status === "no_scffpe_data" ? "No scFFPE data were loaded (scripts/build_scffpe.py has not been run)." : "This Xenium donor has no scFFPE counterpart: the metadata lists no scFFPE sample for it.";
    box.append(Object.assign(el("div"), { innerHTML: `<b>No scFFPE data for this donor.</b>` }), el("div", "", msg));
    if (d.candidates?.length) {
      box.append(el("div", "", `scFFPE samples with the same block ID exist but are not linked in the metadata: ${d.candidates.join(", ")}.`));
      const b = el("button", "", "Compare with these (unlinked, same block)"); b.addEventListener("click", () => { S.cmpOverride = d.candidates.join(","); renderCmp(); }); box.append(b);
    }
    body.append(box); return;
  }
  const sc = d.scffpe, linked = d.status === "ok";
  info.textContent += ` scFFPE ${sc.ids.length > 1 ? "samples" : "sample"} ${sc.ids.join(" + ")}${linked ? "" : " (NOT linked in the metadata; same block ID)"}: ${fmt(sc.n)} cells`
    + (Object.keys(sc.excluded).length ? ` after removing ${Object.entries(sc.excluded).map(([k, n]) => `${fmt(n)} “${k}”`).join(" and ")} (not present in the QC'ed Xenium cells; tick the box to include them)` : "") + "."
    + (d.shared_with.length ? ` This scFFPE sample is shared with ${d.shared_with.join(", ")}.` : "");
  if (!linked) { const b = el("button", "", "◂ back to the metadata link"); b.addEventListener("click", () => { S.cmpOverride = ""; renderCmp(); }); body.append(b); }
  const cats = [...new Set([...Object.keys(d.xenium.counts), ...Object.keys(sc.counts)])].map((k) => ({ k, x: d.xenium.counts[k] || 0, s: sc.counts[k] || 0 }));
  cats.forEach((c) => { c.fx = c.x / d.xenium.n; c.fs = c.s / sc.n; });
  cats.sort((a, b) => b.fx - a.fx || b.fs - a.fs);
  const pct = (f) => (100 * f).toFixed(f < 0.001 && f > 0 ? 2 : 1) + "%";
  if ($("cmpview").classList.contains("on")) {   // table view
    const tb = el("table", "ct"), hr = el("tr"); ["", "Xenium cells", "Xenium %", "scFFPE cells", "scFFPE %", "Xenium − scFFPE (pp)"].forEach((t) => hr.append(el("th", "", t))); tb.append(hr);
    for (const c of cats) { const tr = el("tr"); [c.k, fmt(c.x), pct(c.fx), fmt(c.s), pct(c.fs), ((c.fx - c.fs) * 100 > 0 ? "+" : "") + ((c.fx - c.fs) * 100).toFixed(1)].forEach((t) => tr.append(el("td", "", t))); tb.append(tr); }
    const tr = el("tr"); ["Total", fmt(d.xenium.n), "100%", fmt(sc.n), "100%", ""].forEach((t) => { const td = el("td", "", t); td.style.fontWeight = 700; tr.append(td); }); tb.append(tr); body.append(tb);
  } else {
    const leg = el("div", "cmpleg"); for (const [c, n] of [[CMP_C.xen, "Xenium"], [CMP_C.sc, "scFFPE"]]) { const sp = el("span"), k = el("span", "k"); k.style.background = c; sp.append(k, document.createTextNode(n)); leg.append(sp); } body.append(leg);
    const mx = Math.max(...cats.map((c) => Math.max(c.fx, c.fs)), 1e-9), W = 480;
    for (const c of cats) {
      const row = el("div", "crow"), lab = el("div", "clab", c.k); lab.title = c.k; const bars = el("div", "cbars");
      for (const [col, f, n, who] of [[CMP_C.xen, c.fx, c.x, "Xenium"], [CMP_C.sc, c.fs, c.s, "scFFPE"]]) {
        const r = el("div", "cbar"), b = el("div", "b"), v = el("div", "v", `${pct(f)} (${fmt(n)})`); b.style.width = Math.max(f > 0 ? 2 : 0, (f / mx) * W) + "px"; b.style.background = col;
        r.addEventListener("pointermove", (ev) => nbrTip(ev, [{ tag: "b", t: pct(f) }, { t: `${c.k} · ${who}` }, { cls: "m", t: `${fmt(n)} of ${fmt(who === "Xenium" ? d.xenium.n : sc.n)} cells` }, { cls: "m", t: `Xenium ${pct(c.fx)} · scFFPE ${pct(c.fs)} (${((c.fx - c.fs) * 100 > 0 ? "+" : "") + ((c.fx - c.fs) * 100).toFixed(1)} pp)` }]));
        r.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; }); r.append(b, v); bars.append(r);
      } row.append(lab, bars); body.append(row);
    }
  }
  body.append(Object.assign(el("div", "muted small"), { textContent: "Fractions of all cells of the sample in each platform. Descriptive only: the two platforms capture different cells and genes, so differences are not corrections or tests.", style: "margin-top:14px;max-width:780px" }));
}

// ---------------------------------------------------------------- neighbours window
// For the cells passing the current filters ("index cells"): composition of the cells within R µm, per disease subtype / type / sample.
const NB_TYPE_ORDER = ["Control", "No LN", "LN"], NB_OTHER = "#6b7785";
let nbrToken = 0;
const nbOff = new Set();   // donors left out of the Neighbors / Neighborhoods windows
const nbSampleOn = (s) => !nbOff.has(s.donor);
function nbrDonorPicker() {
  const pop = $("nbdonpop"), donors = [...new Set(S.meta.samples.map((s) => s.donor))], sub = sortedSubtypes(), sOf = (d) => subtypeOf(S.meta.samples.find((s) => s.donor === d));
  donors.sort((a, b) => sub.indexOf(sOf(a)) - sub.indexOf(sOf(b)) || a.localeCompare(b));
  const upd = () => { $("nbdonbtn").textContent = nbOff.size ? `Donors (${donors.length - nbOff.size}/${donors.length}) ▾` : "Donors ▾"; };
  const redraw = () => { upd(); if (S.nbr) drawNbr(); if (!$("nbh").hidden) renderNbh(); };
  pop.replaceChildren();
  const bar = document.createElement("div"); bar.className = "nbdbar";
  for (const [t, f] of [["All", () => nbOff.clear()], ["None", () => donors.forEach((d) => nbOff.add(d))]]) { const b = document.createElement("button"); b.className = "sm"; b.textContent = t; b.addEventListener("click", () => { f(); nbrDonorPicker(); redraw(); pop.hidden = false; }); bar.append(b); }
  pop.append(bar);
  for (const d of donors) {
    const l = document.createElement("label"), c = document.createElement("input"); c.type = "checkbox"; c.checked = !nbOff.has(d);
    c.addEventListener("change", () => { c.checked ? nbOff.delete(d) : nbOff.add(d); redraw(); });
    l.append(c, document.createTextNode(` ${d} · ${sOf(d)}`)); pop.append(l);
  }
  upd();
}
function nbrIndexSpec() {
  const filters = {}, hidden = {};
  for (const [k, set] of Object.entries(S.filters)) if (set.size && S.meta.cats[k]) filters[k] = [...set];
  const h = S.hiddenCats[S.colorBy]; if (S.meta.cats[S.colorBy] && h?.size) hidden[S.colorBy] = [...h];
  return { filters, hidden, expr: exprSpecs() };
}
function openNbr() {
  closeAllWindows();
  if (!$("nbfield").options.length) $("nbfield").innerHTML = ["cell_type", "cell_state", "region_lvl1", "region_lvl2", "region_lvl3", "glom_quality"].filter((k) => S.meta.cats[k])
    .map((k) => `<option value="${k}">${ROLE_LABEL[k] || labelOf(k)}</option>`).join("");
  $("nbr").hidden = false; nbrDonorPicker(); renderNbr();
}
function closeNbr() { $("nbr").hidden = true; $("nbrtip").hidden = true; }

// ---------------------------------------------------------------- neighbourhoods in space
// One tile per index cell (cells passing the filters) of one sample: the cell outlined in white, its neighbours (same definition / distance
// settings as the Neighbors window) coloured by the current Colour by, other nearby cells faint. Click a tile for the zoomed view / PNG.
let nbhToken = 0;
function nbrParams() { return { method: $("nbmethod").value, distance: +$("nbrad").value, k: +$("nbk").value || 6, steps: +$("nbsteps").value || 1 }; }
function openNbh() {
  if (!S.nbr) return;
  const subs = sortedSubtypes(), withIdx = (id) => (S.nbr.res.samples[id]?.n_index || 0) > 0;
  if (!S.meta.samples.some((s) => withIdx(s.id))) { $("nbrinfo").textContent = "No index cells: the current filters select no cells in any donor."; return; }
  $("nbhsample").innerHTML = `<option value="cur">This donor (${dlabel(S.sample)})</option><option value="all">All donors</option>` + subs.map((t) => `<option value="sub:${t}">Disease subtype: ${t}</option>`).join("");
  $("nbhsample").value = S.nbhScope && [...$("nbhsample").options].some((o) => o.value === S.nbhScope) ? S.nbhScope : (withIdx(S.sample) ? "cur" : "all");
  S.nbhScope = $("nbhsample").value; if (!S.nbhTouched) $("nbhn").value = S.nbhScope === "cur" ? 48 : 12;
  closeAllWindows(); $("nbh").hidden = false; renderNbh();
}
function closeNbh() { $("nbh").hidden = true; }
async function renderNbh() {
  if ($("nbh").hidden) return;
  const tok = ++nbhToken, scope = $("nbhsample").value; S.nbhScope = scope;
  const P = nbrParams(), spec = nbrIndexSpec(), limit = Math.max(4, Math.min(200, +$("nbhn").value || 12)), order = $("nbhorder").value, field = $("nbfield").value, nl = catMeta(field).levels;
  const sub = sortedSubtypes(), withIdx = (id) => !S.nbr || (S.nbr.res.samples[id]?.n_index || 0) > 0;
  const samples = S.meta.samples.filter((s) => nbSampleOn(s) && (scope === "all" || (scope === "cur" ? s.id === S.sample : subtypeOf(s) === scope.slice(4))) && withIdx(s.id))
    .sort((a, b) => sub.indexOf(subtypeOf(a)) - sub.indexOf(subtypeOf(b)) || a.id.localeCompare(b.id));
  const mtxt = nbrMethodText(P), noteTxt = P.method === "radius" ? `radius ${P.distance} µm` : P.method === "knn" ? `k=${P.k} ≤${P.distance} µm` : `Delaunay ${P.steps} step${P.steps > 1 ? "s" : ""} ≤${P.distance} µm`;
  const flt = fltText(spec);
  const orderTxt = order === "random" ? "random" : order === "most" ? "most neighbours first" : "fewest neighbours first";
  const grid = $("nbhgrid"); grid.style.opacity = 1;
  grid.innerHTML = samples.map((s) => `<div class="glomsec" data-s="${s.id}" hidden><div class="sechead"></div><div class="secgrid"></div></div>`).join("");
  const tot = { idx: 0, tiles: 0, done: 0 }, cnt = new Map(), ones = (b) => ({ ...b, pass: new Uint8Array(b.n).fill(1) });  // neighbours need not pass the filters
  const info = () => { $("nbhinfo").textContent = `${fmt(tot.idx)} index cells in ${donTxt(nDonors(samples))} (${flt || "no filter active"}); showing up to ${limit} per donor tissue piece (${orderTxt}), ${fmt(tot.tiles)} tiles. `
    + `White outline = the index cell; coloured cells = its neighbours (${mtxt}); colour: ${colorName()}; faint = other cells nearby.${tot.done < samples.length ? " Loading…" : ""} Click a tile to zoom.`; };
  info(); if (!samples.length) { grid.innerHTML = `<div class="muted">No index cells in the selected donors.</div>`; return; }
  // requests start together; continuous colour scales are pooled over all shown samples, categorical ones need no pooling
  const rp = samples.map((s) => fetch("/api/neighborhoods", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ sample: s.id, ...P, ...spec, limit, order, seed: S.nbhSeed || 1 }) })
    .then(async (r) => { if (!r.ok) throw new Error((await r.text()).slice(0, 200)); return r.json(); }));
  const bp = samples.map((s) => winBundle(s.id).then(ones)), pp = samples.map((s) => Promise.all([galPoly(s.id), galGet(s.id, field)]));
  const sc = await winScale(bp, () => nbhToken, tok).catch((e) => { console.error(e); return null; });
  if (!sc || tok !== nbhToken) return;
  S.ctxs.nbh = { sc, shown: [] }; const shown = S.ctxs.nbh.shown;
  const legend = () => {
    const leg = $("nbhleg");
    if (sc.info.type !== "cat") { leg.innerHTML = scaleLegendHTML(sc.info); return; }
    const row = document.createElement("div"); row.className = "nlegrow";
    [...cnt.entries()].filter(([c]) => c !== undefined).sort((x, y) => y[1] - x[1]).slice(0, 14).forEach(([c]) => {
      const sp = document.createElement("span"), k = document.createElement("span"); k.className = "k"; k.style.background = c === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, c);
      sp.append(k, document.createTextNode(c === MISSING ? "missing" : display(S.colorBy, c))); row.append(sp);
    }); leg.replaceChildren(row);
  };
  legend();
  await Promise.all(samples.map(async (s, k) => {
    try {
      const [res, b, [poly, fcodes]] = await Promise.all([rp[k], bp[k], pp[k]]);
      if (tok !== nbhToken) return;
      const sec = grid.querySelector(`.glomsec[data-s="${CSS.escape(s.id)}"]`); if (!sec || !res.index.length) return;
      sec.hidden = false; tot.idx += res.n_index; tot.tiles += res.index.length;
      sec.querySelector(".sechead").textContent = `${dlabel(s.id)} · ${s.disease_type} · ${subtypeOf(s)} · ${fmt(res.n_index)} index cells · showing ${res.index.length}`;
      const g = sec.querySelector(".secgrid");
      res.index.forEach((ci, r) => {
        const nbrs = res.neighbors[r], cells = [ci, ...nbrs], comp = new Map();
        for (const j of nbrs) comp.set(fcodes[j], (comp.get(fcodes[j]) || 0) + 1);
        const top = [...comp.entries()].sort((x, y) => y[1] - x[1]).slice(0, 2).map(([c, n]) => `${c === MISSING ? "missing" : nl[c]} ${Math.round(100 * n / Math.max(1, nbrs.length))}%`).join(", ");
        const t = document.createElement("div"); t.className = "gtile"; t.title = `${res.ids[r]} — open the zoomed view`;
        const cv = document.createElement("canvas"); cv.width = 328; cv.height = 328;
        const g1 = document.createElement("div"); g1.className = "gid"; g1.textContent = res.ids[r];
        const g2 = document.createElement("div"); g2.className = "gm"; g2.textContent = `${nbrs.length} neighbours${top ? " · " + top : ""}`;
        t.append(cv, g1, g2); g.append(t); shown.push({ b, cells, el: t });
        t.addEventListener("click", () => openGlomZoom({ src: "nbh", s, b, sc, poly, cells, focus: ci, label: `${res.ids[r]} · neighbourhood`, qc: null, note: noteTxt }));
        cv._hw = "nbh"; drawGlom(cv, b, sc, poly, cells, true, ci);
        if (b.codes) for (const i of cells) cnt.set(b.codes[i], (cnt.get(b.codes[i]) || 0) + 1);
      });
    } catch (e) { console.error(e); }
    if (tok === nbhToken) { tot.done++; info(); legend(); }
  }));
  if (tok === nbhToken) { tot.done = samples.length; info(); if (!tot.tiles) grid.innerHTML = `<div class="muted">No index cells in the selected donors.</div>`; }
}
async function renderNbr() {
  if ($("nbr").hidden) return;
  const tok = ++nbrToken, field = $("nbfield").value, method = $("nbmethod").value, distance = +$("nbrad").value, k = +$("nbk").value || 6, steps = +$("nbsteps").value || 1, spec = nbrIndexSpec();
  $("nbrinfo").textContent = "Computing neighbours…"; $("nbrbody").style.opacity = 0.45;
  let res;
  try {
    const r = await fetch("/api/neighbors", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ field, method, distance, k, steps, ...spec }) });
    if (!r.ok) throw new Error((await r.text()).slice(0, 200)); res = await r.json();
  } catch (e) { if (tok === nbrToken) $("nbrinfo").textContent = "Error: " + e.message; return; }
  if (tok !== nbrToken) return;
  S.nbr = { res, spec }; $("nbrbody").style.opacity = 1; drawNbr();
}
function nbrGroups(res) {  // aggregate per-sample counts into disease subtype / disease type / sample groups (skipping groups without index cells)
  const nl = res.levels.length, meta = Object.fromEntries(S.meta.samples.map((s) => [s.id, s]));
  const mk = (name) => ({ name, counts: new Array(nl).fill(0), bg: new Array(nl).fill(0), n_index: 0, n_pairs: 0, n_cells: 0 });
  const by = { subtype: new Map(), type: new Map(), sample: new Map() }, subOf = (dn) => subtypeOf(S.meta.samples.find((s) => s.donor === dn));   // "sample" groups are donors: the neighbours are found inside each tissue piece, the counts are then summed per donor
  for (const [sid, r] of Object.entries(res.samples)) {
    if (!r.n_index || !nbSampleOn(meta[sid])) continue; const m = meta[sid];
    for (const [lvl, name] of [["subtype", subtypeOf(m)], ["type", m.disease_type || "unknown"], ["sample", m.donor]]) {
      const g = by[lvl].get(name) || by[lvl].set(name, mk(name)).get(name);
      r.counts.forEach((c, i) => { g.counts[i] += c; g.bg[i] += r.bg[i]; }); g.n_index += r.n_index; g.n_pairs += r.n_pairs; g.n_cells += r.n_cells;
    }
  }
  const rank = (arr, t) => (arr.includes(t) ? arr.indexOf(t) : 99), sub = sortedSubtypes();
  return {
    subtype: [...by.subtype.values()].sort((a, b) => rank(sub, a.name) - rank(sub, b.name)),
    type: [...by.type.values()].sort((a, b) => rank(NB_TYPE_ORDER, a.name) - rank(NB_TYPE_ORDER, b.name)),
    sample: [...by.sample.values()].map((g) => ({ ...g, subtype: subOf(g.name) })).sort((a, b) => rank(sub, a.subtype) - rank(sub, b.subtype) || a.name.localeCompare(b.name)),
  };
}
const nbrColor = (field, i, nl) => (i === nl - 1 ? NB_OTHER : catColor(field, i));
const inkOn = (css) => { const m = css.match(/\d+/g).map(Number), L = (0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2]) / 255; return L > 0.55 ? "#0f1218" : "#ffffff"; };
const hexRgb = (c) => (c.startsWith("#") ? `rgb(${hex2rgb(c).join(",")})` : c);
function nbrEnrich(g, i) { const f = g.n_pairs ? g.counts[i] / g.n_pairs : 0, b = g.n_cells ? g.bg[i] / g.n_cells : 0; return { f, b, e: f > 0 && b > 0 ? Math.log2(f / b) : null }; }
function nbrTip(ev, lines) {
  const t = $("nbrtip"); t.replaceChildren(...lines.map((l) => { const d = document.createElement(l.tag || "div"); if (l.cls) d.className = l.cls; d.textContent = l.t; return d; }));
  t.hidden = false; t.style.left = Math.min(ev.clientX + 14, innerWidth - 350) + "px"; t.style.top = Math.min(ev.clientY + 14, innerHeight - 110) + "px";
}
function nbrMethodText(res) {
  if (res.method === "knn") return `the ${res.k} nearest cells of each index cell (only those within ${res.distance} µm)`;
  if (res.method === "delaunay") return `cells within ${res.steps} step${res.steps > 1 ? "s" : ""} of an index cell in the Delaunay triangulation of the cell centres (triangulation edges longer than ${res.distance} µm removed)`;
  return `all cells within ${res.distance} µm of an index cell`;
}
// Index cells per donor (tissue pieces summed, zeros included), grouped by disease subtype in the fixed order; reflects every current filter.
function nbrCountsTable(res) {
  const meta = Object.fromEntries(S.meta.samples.map((s) => [s.id, s])), sub = sortedSubtypes();
  const rowsBy = new Map();   // one row per donor: the index cells and cells of its tissue pieces are summed
  for (const s of S.meta.samples) { const r = rowsBy.get(s.donor) || rowsBy.set(s.donor, { id: s.donor, type: s.disease_type || "unknown", sub: subtypeOf(s), n: 0, tot: 0 }).get(s.donor); r.n += res.samples[s.id]?.n_index || 0; r.tot += res.samples[s.id]?.n_cells || s.n; }
  const rows = [...rowsBy.values()]
    .sort((a, b) => sub.indexOf(a.sub) - sub.indexOf(b.sub) || a.id.localeCompare(b.id));
  const el = (tag, cls, txt) => { const e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; };
  const pct = (n, t) => (t ? (100 * n / t).toFixed(1) + "%" : "–"), tb = el("table", "nc"), hr = el("tr");
  for (const [t, cls] of [["Disease subtype", ""], ["Donor", ""], ["Disease type", ""], ["Index cells", "num"], ["Cells in donor tissue", "num"], ["% of donor tissue", "num"]]) hr.append(el("th", cls, t));
  tb.append(hr);
  const add = (cells, cls) => { const tr = el("tr", cls); cells.forEach(([t, c]) => tr.append(el("td", c, t))); tb.append(tr); };
  for (const sb of sub) {
    const g = rows.filter((r) => r.sub === sb); if (!g.length) continue; const n = g.reduce((a, r) => a + r.n, 0), t = g.reduce((a, r) => a + r.tot, 0);
    add([[sb, ""], [donTxt(nDonors(g.map((r) => r.id))), ""], ["", ""], [fmt(n), "num"], [fmt(t), "num"], [pct(n, t), "num"]], "sub");
    for (const r of g) add([[r.sub, ""], [dlabel(r.id), ""], [r.type, ""], [fmt(r.n), "num" + (r.n ? "" : " zero")], [fmt(r.tot), "num"], [pct(r.n, r.tot), "num"]], "");
  }
  const N = rows.reduce((a, r) => a + r.n, 0), T = rows.reduce((a, r) => a + r.tot, 0);
  add([["All donors", ""], [donTxt(nDonors(rows.map((r) => r.id))), ""], ["", ""], [fmt(N), "num"], [fmt(T), "num"], [pct(N, T), "num"]], "tot");
  const box = el("div", "nsec"); box.append(el("h4", "", "Index cells per donor (cells passing the current filters, region filters included)"), tb); return box;
}
function drawNbr() {
  if (!S.nbr) return;
  const { res, spec } = S.nbr, nl = res.levels.length, view = $("nbview").value, field = res.field, groups = nbrGroups(res), body = $("nbrbody");
  const tot = groups.sample.reduce((a, g) => ({ idx: a.idx + g.n_index, pairs: a.pairs + g.n_pairs }), { idx: 0, pairs: 0 });
  const flt = fltText(spec);
  $("nbrinfo").textContent = `Xenium cells only (QC'ed; no scFFPE). ${fmt(tot.idx)} index cells in ${donTxt(nDonors(groups.sample.map((g) => g.name)))} (${flt || "no filter active: every cell is an index cell, so this is just the tissue composition"}${Object.keys(spec.hidden).length ? " · legend-hidden types excluded" : ""}). `
    + `Neighbours = ${nbrMethodText(res)}, counted per index cell–neighbour pair (${fmt(tot.pairs)} pairs). Descriptive only, no statistical test; groups with few pairs are noisy.`;
  const pooled = new Array(nl).fill(0); groups.sample.forEach((g) => g.counts.forEach((c, i) => { pooled[i] += c; }));
  const order = pooled.map((c, i) => [i, c]).filter(([, c]) => c > 0).sort((a, b) => b[1] - a[1]).map(([i]) => i);
  const secs = [["subtype", "By disease subtype"], ["type", "By disease type"], ["sample", "By donor"]];
  body.replaceChildren(); if (S.nbCounts) body.append(nbrCountsTable(res));
  const el = (tag, cls, txt) => { const e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; };
  $("nbgrouplab").hidden = !(view === "dist" || view === "forest");
  if (view === "dist" || view === "forest") { nbrDistView(view, groups, order); return; }
  if (view === "bars") {  // part-to-whole: top 8 neighbour types keep their entity colour, the rest fold into "Other"
    const top = order.slice(0, 8), rest = order.slice(8), leg = el("div", "nlegend");
    for (const i of top.concat(rest.length ? [-1] : [])) { const s = el("span"); const k = el("span", "k"); k.style.background = i < 0 ? NB_OTHER : nbrColor(field, i, nl); s.append(k, document.createTextNode(i < 0 ? `Other (${rest.length} types)` : res.levels[i])); leg.append(s); }
    body.append(leg);
    for (const [lvl, title] of secs) {
      if (!groups[lvl].length) continue; const sec = el("div", "nsec"); sec.append(el("h4", "", title));
      for (const g of groups[lvl]) {
        const row = el("div", "nrow"), bar = el("div", "nbar"), lab = el("div", "nlab", lvl === "sample" ? dlabel(g.name) : g.name); lab.title = (lvl === "sample" ? dlabel(g.name) : g.name) + (g.subtype ? ` · ${g.subtype}` : "");
        if (g.subtype) { const sm = el("small", "", g.subtype); lab.append(sm); }
        const parts = top.map((i) => ({ i, c: g.counts[i], col: nbrColor(field, i, nl) })); if (rest.length) parts.push({ i: -1, c: rest.reduce((a, i) => a + g.counts[i], 0), col: NB_OTHER });
        const px = Math.max(300, body.clientWidth - 340);
        for (const p of parts) {
          if (!p.c) continue; const frac = p.c / g.n_pairs, seg = el("div", "nseg"); seg.style.width = frac * 100 + "%"; seg.style.background = p.col;
          if (frac * px > 44) { seg.textContent = Math.round(frac * 100) + "%"; seg.style.color = inkOn(hexRgb(p.col)); }
          seg.addEventListener("pointermove", (ev) => {
            const names = p.i < 0 ? `Other (${rest.length} types)` : res.levels[p.i], bgc = p.i < 0 ? rest.reduce((a, i) => a + g.bg[i], 0) : g.bg[p.i], b = bgc / g.n_cells;
            nbrTip(ev, [{ tag: "b", t: (frac * 100).toFixed(1) + "%" }, { t: names }, { cls: "m", t: `${lvl === "sample" ? dlabel(g.name) : g.name}${g.subtype ? " (" + g.subtype + ")" : ""} · ${fmt(p.c)} of ${fmt(g.n_pairs)} neighbour pairs` },
              { cls: "m", t: `background ${(b * 100).toFixed(1)}% → ×${b ? (frac / b).toFixed(2) : "–"}` }]);
          });
          seg.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
          bar.append(seg);
        }
        row.append(lab, bar, el("div", "nn", `${fmt(g.n_index)} index · ${fmt(g.n_pairs)} pairs`)); sec.append(row);
      }
      body.append(sec);
    }
  } else {  // heatmap (log2 enrichment, diverging around 0) or table of percentages
    const heat = view === "heat", cols = order;
    if (heat) { const leg = el("div", "nlegend"); const bar = el("div"); bar.style.cssText = `width:260px;height:10px;border-radius:3px;background:linear-gradient(90deg,${DIVL.filter((_, i) => i % 32 === 0 || i === 255).map((c) => `rgb(${c})`).join(",")})`;
      leg.append(el("span", "", "−2 (depleted)"), bar, el("span", "", "+2 (enriched)"), el("span", "m", "log2( neighbour % ÷ background % ), background = all cells of the group; grey = not seen")); body.append(leg); }
    for (const [lvl, title] of secs) {
      if (!groups[lvl].length) continue; const sec = el("div", "nsec"); sec.append(el("h4", "", title));
      const tb = el("table", "nh" + (heat ? "" : " plain")), hr = el("tr"); hr.append(el("th", "r", ""));
      for (const i of cols) { const th = el("th", "c", res.levels[i]); th.title = res.levels[i]; hr.append(th); } hr.append(el("th", "", "pairs")); tb.append(hr);
      for (const g of groups[lvl]) {
        const tr = el("tr"), th = el("th", "r", g.subtype ? `${g.name} · ${g.subtype}` : g.name); th.title = th.textContent; tr.append(th);
        for (const i of cols) {
          const { f, b, e } = nbrEnrich(g, i), td = el("td");
          if (heat) {
            const v = e == null ? (g.counts[i] === 0 && b > 0 ? -2 : null) : Math.max(-2, Math.min(2, e));
            if (v == null) { td.style.background = "#141922"; td.textContent = "·"; td.style.color = "#5b6676"; }
            else { const c = DIVL[Math.round(((v + 2) / 4) * 255)]; td.style.background = `rgb(${c})`; td.style.color = inkOn(`rgb(${c})`); td.textContent = e == null ? "0" : (e > 0 ? "+" : "") + e.toFixed(1); }
          } else td.textContent = `${(f * 100).toFixed(1)}%`;
          td.addEventListener("pointermove", (ev) => nbrTip(ev, [{ tag: "b", t: `${(f * 100).toFixed(1)}% of neighbours` }, { t: res.levels[i] }, { cls: "m", t: `${g.name}${g.subtype ? " (" + g.subtype + ")" : ""} · ${fmt(g.counts[i])} of ${fmt(g.n_pairs)} pairs` },
            { cls: "m", t: `background ${(b * 100).toFixed(1)}%${e == null ? "" : " · log2 enrichment " + e.toFixed(2)}` }]));
          td.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; }); tr.append(td);
        }
        tr.append(el("td", "n", fmt(g.n_pairs))); tb.append(tr);
      }
      sec.append(tb); body.append(sec);
    }
  }
}
// Spread across index cells (violin + box per neighbour type) and pooled enrichment with a cluster-bootstrap 95% interval (forest plot). Both come from
// /api/neighbors/dist, which pools the index cells of each group server-side; the pooled stacked bars above hide how much single cells differ.
const NB_GC = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#B8A100", "#999999"];
const nbSvg = (tag, attrs = {}, txt) => { const e = document.createElementNS("http://www.w3.org/2000/svg", tag); for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v); if (txt != null) e.textContent = txt; return e; };
async function nbrDistData(res, spec) {
  const meta = S.meta.samples, groups = { subtype: {}, type: {}, sample: {} };
  for (const m of meta.filter(nbSampleOn)) { groups.subtype[m.id] = subtypeOf(m); groups.type[m.id] = m.disease_type || "unknown"; groups.sample[m.id] = m.donor; }
  const body = JSON.stringify({ field: res.field, method: res.method, distance: res.distance, k: res.k, steps: res.steps, ...spec, groups });
  if (S.nbr.dist?.key !== body) {
    const r = await fetch("/api/neighbors/dist", { method: "POST", headers: { "Content-Type": "application/json" }, body });
    if (!r.ok) throw new Error((await r.text()).slice(0, 200)); S.nbr.dist = { key: body, data: await r.json() };
  }
  return S.nbr.dist.data;
}
async function nbrDistView(view, groups, order) {
  const body = $("nbrbody"), tok = nbrToken, snap = S.nbr, lvl = $("nbgroup").value, names = groups[lvl].map((g) => g.name), nl = snap.res.levels.length;
  body.replaceChildren(Object.assign(document.createElement("div"), { className: "muted", textContent: "Computing per-cell distributions…" })); body.style.opacity = 0.45;
  let d; try { d = await nbrDistData(snap.res, snap.spec); } catch (e) { body.textContent = "Error: " + e.message; body.style.opacity = 1; return; }
  if (tok !== nbrToken || snap !== S.nbr || $("nbview").value !== view) return;
  body.style.opacity = 1; body.replaceChildren();
  const el = (tag, cls, txt) => { const e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; };
  if (S.nbCounts) body.append(nbrCountsTable(snap.res));
  const by = new Map((d.groups[lvl] || []).map((g) => [g.name, g])), gs = names.map((n) => by.get(n)).filter((g) => g && g.n_with_neighbours > 0), field = d.field;
  if (!gs.length) { body.append(el("div", "muted", "No index cell has a neighbour under these settings.")); return; }
  const types = order.filter((i) => i < nl).slice(0, view === "dist" ? 10 : 12), nm = (g) => (lvl === "sample" ? dlabel(g.name) : g.name);
  const pct = (x) => (x * 100).toFixed(1) + "%", fmtE = (x) => (x == null ? "–" : (x > 0 ? "+" : "") + x.toFixed(2));
  if (view === "dist") {
    body.append(el("div", "muted small", `Each violin = the distribution, over the index cells of a group, of the fraction of that cell's neighbours that are of this type (cells without neighbours left out). Box = 25–75th percentile with the median bar, whiskers = 5th–95th, dot = mean. `
      + `Shapes are scaled to equal width, so read the y-axis, not the width. A type at ~9% pooled can be ~0% for most cells and ~70% for a few.`));
    const wrap = el("div", "ndist"), bw = gs.length > 12 ? 26 : 46, H = 170, mt = 22, mb = 74, ml = 44;
    for (const i of types) {
      const W = ml + gs.length * bw + 10, ymax = Math.min(1, Math.max(0.05, Math.ceil(Math.max(...gs.map((g) => g.q[i]?.[4] ?? 0)) * 1.1 * 20) / 20));
      const svg = nbSvg("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` }), top = mt, bot = H - mb;
      svg.append(nbSvg("text", { x: 4, y: 13, fill: "#d7dde6", "font-size": 12, "font-weight": 600 }, d.levels[i]));
      for (const t of [0, ymax / 2, ymax]) { const y = bot - (bot - top) * (t / ymax); svg.append(nbSvg("line", { x1: ml, x2: W - 6, y1: y, y2: y, stroke: "#29303b" }), nbSvg("text", { x: ml - 4, y: y + 4, fill: "#8b96a6", "font-size": 10, "text-anchor": "end" }, Math.round(t * 100) + "%")); }
      const col = nbrColor(field, i, nl + 1), Yf = (f) => bot - (bot - top) * (Math.min(f, ymax) / ymax);
      gs.forEach((g, gi) => {
        const cx = ml + gi * bw + bw / 2, h = g.hist[i], nb = h.length, sm = h.map((_, b) => (h[Math.max(0, b - 1)] + 2 * h[b] + h[Math.min(nb - 1, b + 1)]) / 4), mx = Math.max(...sm, 1), hw = bw * 0.42;
        const pts = []; for (let b = 0; b < nb; b++) { const f = (b + 0.5) / nb; if (f - 0.5 / nb > ymax) break; pts.push([f, sm[b] / mx * hw]); }
        const left = [`${cx},${Yf(0)}`, ...pts.map(([f, w]) => `${cx - w},${Yf(f)}`)], right = pts.slice().reverse().map(([f, w]) => `${cx + w},${Yf(f)}`);
        const vi = nbSvg("polygon", { points: [...left, ...right].join(" "), fill: col, "fill-opacity": 0.35, stroke: col, "stroke-width": 0.8 });
        const [q5, q25, q50, q75, q95, mean] = g.q[i], tipL = [{ tag: "b", t: `${d.levels[i]} · ${nm(g)}` }, { t: `median ${pct(q50)} · mean ${pct(mean)}` }, { cls: "m", t: `IQR ${pct(q25)}–${pct(q75)} · 5–95%: ${pct(q5)}–${pct(q95)}` }, { cls: "m", t: `${fmt(g.n_with_neighbours)} index cells with neighbours (${fmt(g.n_index)} index) · background ${pct(g.bg[i])}` }];
        for (const e of [vi]) { e.addEventListener("pointermove", (ev) => nbrTip(ev, tipL)); e.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; }); }
        svg.append(vi, nbSvg("line", { x1: cx, x2: cx, y1: Yf(q5), y2: Yf(q95), stroke: "#e6ebf2", "stroke-width": 1 }), nbSvg("rect", { x: cx - 4, y: Yf(q75), width: 8, height: Math.max(1, Yf(q25) - Yf(q75)), fill: "#0f1218", stroke: "#e6ebf2" }),
          nbSvg("line", { x1: cx - 4, x2: cx + 4, y1: Yf(q50), y2: Yf(q50), stroke: "#ffffff", "stroke-width": 2 }), nbSvg("circle", { cx, cy: Yf(mean), r: 2.5, fill: "#ffd166" }));
        const lab = nbSvg("text", { x: cx, y: bot + 10, fill: "#aab4c3", "font-size": 10, "text-anchor": "end", transform: `rotate(-40 ${cx} ${bot + 10})` }, nm(g).length > 16 ? nm(g).slice(0, 15) + "…" : nm(g)); svg.append(lab);
      });
      wrap.append(svg);
    }
    body.append(wrap);
  } else {
    body.append(el("div", "muted small", `Pooled log2( neighbour fraction ÷ background fraction ) per group, background = all cells of the group's tissue. Whiskers = 95% interval of a bootstrap that resamples whole index cells (pairs sharing an index cell are not independent), so wide whiskers mean the enrichment rests on few or very uneven index cells. `
      + `× at the left edge = the type is never a neighbour. 0 = same as background; the interval does not guarantee significance and nothing here is multiple-test corrected.`));
    const lg = el("div", "nlegend"); gs.forEach((g, gi) => { const sp = el("span"), k = el("span", "k"); k.style.background = NB_GC[gi % NB_GC.length]; sp.append(k, document.createTextNode(`${nm(g)} (${fmt(g.n_with_neighbours)} cells)`)); lg.append(sp); }); body.append(lg);
    const LIM = 4, ml = 170, pw = 520, rowH = gs.length * 11 + 12, H = 28 + types.length * rowH + 22, X = (v) => ml + ((Math.max(-LIM, Math.min(LIM, v)) + LIM) / (2 * LIM)) * pw;
    const svg = nbSvg("svg", { width: ml + pw + 20, height: H, viewBox: `0 0 ${ml + pw + 20} ${H}` });
    for (let v = -LIM; v <= LIM; v += 2) { svg.append(nbSvg("line", { x1: X(v), x2: X(v), y1: 20, y2: H - 20, stroke: v === 0 ? "#8b96a6" : "#29303b", "stroke-dasharray": v === 0 ? "" : "3 3" }), nbSvg("text", { x: X(v), y: H - 6, fill: "#8b96a6", "font-size": 10, "text-anchor": "middle" }, (v > 0 ? "+" : "") + v)); }
    svg.append(nbSvg("text", { x: ml + pw / 2, y: 12, fill: "#8b96a6", "font-size": 10, "text-anchor": "middle" }, "← depleted      log2 enrichment vs background      enriched →"));
    types.forEach((i, ti) => {
      const y0 = 24 + ti * rowH; if (ti) svg.append(nbSvg("line", { x1: 0, x2: ml + pw, y1: y0 - 2, y2: y0 - 2, stroke: "#1d232d" }));
      svg.append(nbSvg("text", { x: ml - 8, y: y0 + rowH / 2, fill: "#d7dde6", "font-size": 12, "text-anchor": "end" }, d.levels[i]));
      gs.forEach((g, gi) => {
        const y = y0 + 6 + gi * 11, c = NB_GC[gi % NB_GC.length], e = g.enr[i], lo = g.lo[i], hi = g.hi[i];
        const tip = (ev) => nbrTip(ev, [{ tag: "b", t: `${d.levels[i]} · ${nm(g)}` }, { t: e == null ? "never a neighbour" : `log2 enrichment ${fmtE(e)} (×${Math.pow(2, e).toFixed(2)})` }, { cls: "m", t: lo == null ? "interval not available (type rarely a neighbour)" : `95% bootstrap interval ${fmtE(lo)} to ${fmtE(hi)}` }, { cls: "m", t: `${fmt(g.n_with_neighbours)} index cells · ${fmt(g.n_pairs)} pairs · background ${pct(g.bg[i])}` }]);
        const mark = e == null ? nbSvg("text", { x: X(-LIM), y: y + 4, fill: c, "font-size": 11, "text-anchor": "middle" }, "×") : nbSvg("circle", { cx: X(e), cy: y, r: 3.5, fill: c });
        const els = [mark]; if (lo != null && hi != null) els.unshift(nbSvg("line", { x1: X(lo), x2: X(hi), y1: y, y2: y, stroke: c, "stroke-width": 1.6 }));
        els.forEach((x) => { x.addEventListener("pointermove", tip); x.addEventListener("pointerleave", () => { $("nbrtip").hidden = true; }); svg.append(x); });
      });
    });
    body.append(svg);
  }
}
function exportNbrCSV() {
  if (!S.nbr) return; const { res, spec } = S.nbr, groups = nbrGroups(res), q = (x) => `"${String(x).replace(/"/g, '""')}"`;
  const lines = [["group_level", "group", "disease_subtype", "n_index_cells", "n_pairs", "neighbour_type", "neighbour_pairs", "neighbour_fraction", "background_cells", "background_fraction", "log2_enrichment", "method", "distance_um", "k", "steps", "neighbour_field", "index_filters"].join(",")];
  const flt = JSON.stringify({ ...Object.fromEntries(Object.entries(spec.filters).map(([k, ids]) => [k, ids.map((i) => display(k, i))])), ...((spec.expr || []).length ? { gene_signature_filters: spec.expr.map(exprDesc) } : {}) });
  for (const lvl of ["subtype", "type", "sample"]) for (const g of groups[lvl]) res.levels.forEach((name, i) => {
    const { f, b, e } = nbrEnrich(g, i); lines.push([lvl, q(g.name), q(g.subtype || (lvl === "subtype" ? g.name : "")), g.n_index, g.n_pairs, q(name), g.counts[i], f.toFixed(6), g.bg[i], b.toFixed(6), e == null ? "" : e.toFixed(4), res.method, res.distance, res.method === "knn" ? res.k : "", res.method === "delaunay" ? res.steps : "", q(res.field), q(flt)].join(","));
  });
  download(`neighbors_${res.field}_${res.method}.csv`, lines.join("\n"));
}

// ---------------------------------------------------------------- UMAP (floating window) + cell highlight
// All cells in one scatter (514k points). S.hl = { sample, primary: [...], secondary: [...] } (sample-local cell indices): primary = the clicked / selected
// cells (large, outlined), secondary = related cells (e.g. its neighbours); with a highlight every other cell is grey.
let umTok = 0; const UM = { emb: {}, cat: {}, deck: null, bounds: {} };
const umGlobal = (sample, i) => S.meta.samples.find((s) => s.id === sample).offset + i;
function umSample(g) {  // sample holding global cell index g (samples are contiguous slices)
  const ss = S.meta.samples; let lo = 0, hi = ss.length - 1;
  while (lo < hi) { const mid = (lo + hi + 1) >> 1; ss[mid].offset <= g ? (lo = mid) : (hi = mid - 1); } return ss[lo];
}
function setHighlight(sample, primary, secondary = []) { S.hl = { sample, primary: [...primary], secondary: [...secondary] }; updateUmap(); }
function clearHighlight() { S.hl = null; updateUmap(); }
function toggleUmap(on) {
  S.umOpen = on ?? !S.umOpen; $("umapwin").hidden = !S.umOpen; $("umapbtn").classList.toggle("on", S.umOpen);
  if (S.umOpen) { if (!UM.deck) initUmapDeck(); updateUmap(true); }
}
function initUmapDeck() {
  const el = $("umdeck");
  UM.deck = new deck.Deck({ parent: el, views: new deck.OrthographicView({ flipY: false }), initialViewState: { target: [0, 0, 0], zoom: 3 }, controller: true, layers: [], width: el.clientWidth || 440, height: el.clientHeight || 440,
    getTooltip: ({ index, layer }) => {
      if (!(index >= 0) || layer?.id !== "um-base" || !UM.last) return null; const g = index, s = umSample(g), c = UM.last.codes[g], lv = catMeta(UM.last.key).levels;
      return { text: `${dlabel(s.id)}\n${UM.last.label}: ${c === MISSING ? "missing" : lv[c]}`, style: { fontSize: "12px", background: "#1b212b", color: "#e3e8ef", border: "1px solid #39424f" } };
    } });
  if (window.ResizeObserver) new ResizeObserver(() => { if (S.umOpen && UM.deck) UM.deck.setProps({ width: el.clientWidth, height: el.clientHeight }); }).observe(el);
}
function umFit() {
  const b = UM.bounds[UM.last?.name], el = $("umdeck"); if (!b || !UM.deck) return;
  const w = el.clientWidth || 440, h = el.clientHeight || 440, dx = b[1] - b[0], dy = b[3] - b[2];
  UM.vn = (UM.vn || 0) + 1;
  UM.deck.setProps({ initialViewState: { target: [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, 0], zoom: Math.log2(Math.min(w / dx, h / dy) * 0.92), _n: UM.vn } });
}
async function updateUmap(fit = false) {
  if (!S.umOpen || !UM.deck) return;
  const tok = ++umTok, kind = $("umemb").value, px = +$("umpcx").value || 1, py = +$("umpcy").value || 2, mode = S.meta.cats[S.colorBy] ? S.colorBy : "cell_type";
  const name = kind === "pca" ? `pca-${px}-${py}` : "umap", url = kind === "pca" ? `/api/pca/${px}/${py}` : "/api/umap/harmony";
  const fkeys = Object.keys(S.filters).filter((k) => S.filters[k].size && S.meta.cats[k]), hid = S.hiddenCats[mode] || new Set();
  $("umnote").textContent = UM.emb[name] ? $("umnote").textContent : "Loading UMAP…";
  let xm = null;
  try {
    if (!UM.emb[name]) UM.emb[name] = await bget(url, Float32Array);
    for (const k of new Set([mode, ...fkeys])) if (!UM.cat[k]) UM.cat[k] = await bget(`/api/cat_all/${encodeURIComponent(k)}`, Uint16Array);
    if (S.exprFilters.length) { $("umnote").textContent = "Applying gene / signature filters…"; xm = await exprMaskAll(); }
  } catch (e) { $("umnote").textContent = "UMAP not available: " + e.message; return; }
  if (tok !== umTok) return;
  const emb = UM.emb[name], codes = UM.cat[mode], n = codes.length;
  if (!UM.bounds[name]) {  // fit to the 0.5-99.5 percentile range so a few outliers (common in PCA) don't squeeze the bulk
    const xs = new Float32Array(Math.ceil(n / 7)), ys = new Float32Array(xs.length); for (let i = 0, j = 0; i < n; i += 7, j++) { xs[j] = emb[2 * i]; ys[j] = emb[2 * i + 1]; } xs.sort(); ys.sort();
    const q = (a, f) => a[Math.min(a.length - 1, Math.floor(f * a.length))]; UM.bounds[name] = [q(xs, 0.005), q(xs, 0.995), q(ys, 0.005), q(ys, 0.995)]; fit = true; }
  const lut = catMeta(mode).levels.map((_, i) => hex2rgb(catColor(mode, i))), hl = S.hl, hasHl = !!hl && (hl.primary.length + hl.secondary.length) > 0;
  const colors = new Uint8Array(n * 4);
  const keep = new Uint8Array(n).fill(1);
  for (const k of fkeys) { const set = S.filters[k], c = UM.cat[k], mu = isMulti(k), sm = mu ? selMask(set) : 0; for (let i = 0; i < n; i++) if (keep[i] && !(mu ? (c[i] & sm) !== 0 : set.has(c[i]))) keep[i] = 0; }
  if (xm) for (let i = 0; i < n; i++) if (keep[i] && !xm[i]) keep[i] = 0;
  if (hid.size) for (let i = 0; i < n; i++) if (keep[i] && hid.has(codes[i])) keep[i] = 0;
  for (let i = 0; i < n; i++) {
    const o = i * 4;
    if (hasHl || !keep[i]) { colors[o] = 105; colors[o + 1] = 110; colors[o + 2] = 122; colors[o + 3] = hasHl ? 70 : 90; }   // everything else grey
    else { const c = codes[i] === MISSING ? GRAY : lut[codes[i]]; colors[o] = c[0]; colors[o + 1] = c[1]; colors[o + 2] = c[2]; colors[o + 3] = 215; }
  }
  const label = ROLE_LABEL[mode] || labelOf(mode); UM.last = { key: mode, label, codes, name };
  const gl = (a) => (hl ? a.map((i) => umGlobal(hl.sample, i)) : []), prim = gl(hl?.primary || []), sec = gl(hl?.secondary || []);
  const colOf = (g) => { const c = codes[g]; return c === MISSING ? GRAY : lut[c]; };
  const layers = [new deck.ScatterplotLayer({ id: "um-base", data: { length: n, attributes: { getPosition: { value: emb, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true, type: "unorm8" } } },
    radiusUnits: "pixels", getRadius: 1.3, radiusMinPixels: 0.8, pickable: true })];
  if (sec.length) layers.push(new deck.ScatterplotLayer({ id: "um-sec", data: sec, getPosition: (g) => [emb[2 * g], emb[2 * g + 1]], getFillColor: (g) => [...colOf(g), 255], radiusUnits: "pixels", getRadius: 3.2, pickable: false }));
  if (prim.length) layers.push(new deck.ScatterplotLayer({ id: "um-prim", data: prim, getPosition: (g) => [emb[2 * g], emb[2 * g + 1]], getFillColor: (g) => [...colOf(g), 255], stroked: true, getLineColor: [255, 255, 255, 255],
    lineWidthUnits: "pixels", getLineWidth: 1.6, radiusUnits: "pixels", getRadius: 6, pickable: false }));
  UM.deck.setProps({ layers });
  $("umtitle").textContent = "Embeddings"; $("umtitle").title = kind === "pca" ? `${S.meta.embedding_label ? S.meta.embedding_label + " PCA" : "Harmony-corrected PCA"}: PC${px} vs PC${py}` : (S.meta.embedding_label ? `${S.meta.embedding_label} UMAP` : "Harmony-integrated UMAP");
  $("umnote").textContent = (hasHl ? `${fmt(prim.length)} highlighted${sec.length ? ` + ${fmt(sec.length)} related` : ""} in ${hl.sample} · others grey` : "Click a cell in the viewer, a tile or a zoom view to highlight it here")
    + ` · coloured by ${label}${S.colorBy !== mode ? " (gene/numeric colouring isn't shown here)" : ""}`;
  if (fit) umFit();
}
function polyContains(poly, i, x, y) {  // ray casting over cell i's boundary
  let inside = false; const a = poly.off[i] * 2, b = poly.off[i + 1] * 2;
  for (let q = a, p = b - 2; q < b; p = q, q += 2) {
    const xi = poly.xy[q], yi = poly.xy[q + 1], xj = poly.xy[p], yj = poly.xy[p + 1];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  } return inside;
}
// Cell under the pointer in the zoomed glomerulus / neighbourhood image (-1 = none).
function glomZoomCellAt(ev) {
  const v = S.gz; if (!v || !v.xf) return -1; const cv = $("gzcanvas"), r = cv.getBoundingClientRect(), G = 900;
  const cx = ((ev.clientX - r.left) * cv.width) / r.width, cy = ((ev.clientY - r.top) * cv.height) / r.height; if (cy > G) return -1;
  const px = (cx - v.xf.ox) / v.xf.s, py = (cy - v.xf.oy) / v.xf.s; let hit = v.cells.find((i) => polyContains(v.poly, i, px, py));
  for (let i = 0; hit == null && i < v.b.n; i++) { const x = v.b.xy[2 * i], y = v.b.xy[2 * i + 1]; if (Math.abs(x - px) < 40 && Math.abs(y - py) < 40 && polyContains(v.poly, i, px, py)) hit = i; }
  return hit == null ? -1 : hit;
}
// Click on a cell inside the zoomed glomerulus / neighbourhood image: outline it there and highlight it in the UMAP.
function glomZoomPick(ev) {
  if (S.gzMoved) { S.gzMoved = false; return; }
  const v = S.gz, hit = glomZoomCellAt(ev); if (!v || hit < 0) return;
  v.focus = hit; drawGlomZoom(); setHighlight(v.s.id, [hit], v.cells.filter((c) => c !== hit));
}
// Hover: cell type and state of the cell under the pointer (also for the neighbouring cells drawn around the glom).
async function glomZoomHover(ev) {
  const v = S.gz; if (!v) return;
  if (!v.ctp) v.ctp = Promise.all([galGet(v.s.id, "cell_type"), galGet(v.s.id, "cell_state")]);
  const hit = glomZoomCellAt(ev); if (hit < 0) { $("nbrtip").hidden = true; return; }
  const [ct, cs] = await v.ctp, tl = catMeta("cell_type").levels, sl = catMeta("cell_state").levels, inGlom = v.cells.includes(hit);
  nbrTip(ev, [{ tag: "b", t: ct[hit] === MISSING ? "no cell type" : tl[ct[hit]] }, { t: cs[hit] === MISSING ? "no cell state" : sl[cs[hit]] },
    { cls: "m", t: (inGlom ? "in this glomerulus" : "neighbouring cell") + (v.b.pass[hit] ? "" : " · does not pass the filters") }]);
}
// The zoomed glomerulus follows the main window: filters / colour-by / gene / smoothing changes redraw it (scale = this sample's cells).
async function refreshGlomZoom() {
  const v = S.gz; if ($("glomzoom").hidden || !v) return; const tok = (S.gzTok = (S.gzTok || 0) + 1);
  try { const b = await winBundle(v.s.id), sc = await winScale([b], () => S.gzTok, tok); if (!sc || tok !== S.gzTok || S.gz !== v) return; v.b = b; v.sc = sc; drawGlomZoom(); } catch (e) { console.error(e); }
}

// ---------------------------------------------------------------- export
function download(name, text) { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type: "text/csv" })); a.download = name; a.click(); }
function exportCSV(indices, name) {
  const keys = ROLE_ORDER.filter((k) => S.cats[k]);
  const head = ["cell_id", "sample", "donor_id_blinded", "x", "y", ...keys.map(labelOf)];
  const q = (s) => `"${String(s).replace(/"/g, '""')}"`;
  const lines = [head.join(",")];
  for (const i of indices) lines.push([q(S.ids[i]), q(S.sample), q(donorOfId(S.sample)), S.xy[i * 2], S.xy[i * 2 + 1], ...keys.map((k) => { const c = S.cats[k][i]; return c === MISSING ? "" : q(catMeta(k).levels[c]); })].join(","));
  download(name, lines.join("\n"));
}

// ---------------------------------------------------------------- UI wiring
function syncGeneUI() {
  const m = S.colorBy, g = m === "__gene" || m === "__dual" || m === "__sig" || m === "__sigdual";
  $("genebox").style.display = g ? "block" : "none"; $("numbox").style.display = S.meta.numeric[m] ? "block" : "none";
  $("rowA").style.display = m === "__gene" || m === "__dual" ? "flex" : "none";
  $("rowB").style.display = m === "__dual" ? "flex" : "none";
  $("sigbox").style.display = m === "__sig" || m === "__sigdual" ? "block" : "none"; $("sigbox2").style.display = m === "__sigdual" ? "block" : "none";
  $("siglabA").textContent = m === "__sigdual" ? "Signature A (blue)" : "";
}
function bindUI() {
  $("sample").addEventListener("change", (e) => loadSample(e.target.value));
  $("colorby").addEventListener("change", async (e) => {
    S.colorBy = e.target.value;
    if (S.colorBy === "__sigdual" && !S.sigGenes.length && !S.sigGenes2.length) {   // nothing chosen yet: start with the two interferon signatures
      try { await loadSigList(); const a = SIG.byId.ifn1, b = SIG.byId.ifn2;
        S.sigGenes = [...a.genes]; S.sigName = a.name; S.sigKey = ""; $("sigtext").value = a.genes.join(" "); $("sigstatus").textContent = `${a.genes.length} genes of "${a.name}" loaded`;
        S.sigGenes2 = [...b.genes]; S.sigName2 = b.name; S.sigKey2 = ""; $("sigtext2").value = b.genes.join(" "); $("sigstatus2").textContent = `${b.genes.length} genes of "${b.name}" loaded`; } catch (err) { /* signatures not available */ }
    }
    syncGeneUI(); refresh();
  });
  syncGeneUI();
  const searchGenes = debounce(async (e) => { const l = await jget("/api/genes/search?q=" + encodeURIComponent(e.target.value)); $("genelist").innerHTML = l.map((g) => `<option value="${g}">`).join(""); });
  const setGene = (key) => async (e) => {
    const v = e.target.value.trim(); if (v && !S.meta) return;
    if (!v) { S[key] = ""; return refresh(); }
    const hit = (await jget("/api/genes/search?q=" + encodeURIComponent(v) + "&limit=1"))[0];
    if (hit && hit.toLowerCase() === v.toLowerCase()) {
      S[key] = hit; e.target.value = hit;
      const want = key === "geneB" ? "__dual" : (S.colorBy === "__dual" ? "__dual" : "__gene");
      S.colorBy = want; $("colorby").value = want; syncGeneUI(); refresh();
    }
  };
  for (const id of ["geneA", "geneB"]) { $(id).addEventListener("input", searchGenes); $(id).addEventListener("change", setGene(id === "geneA" ? "geneA" : "geneB")); }
  $("applysig").addEventListener("click", async () => {
    const names = $("sigtext").value.split(/[\s,;]+/).filter(Boolean);
    if (!names.length) { S.sigGenes = []; $("sigstatus").textContent = ""; }
    else {
      const r = await jget("/api/genes/resolve?names=" + encodeURIComponent(names.join(",")));
      S.sigGenes = r.found;
      $("sigstatus").innerHTML = `${r.found.length} of ${new Set(names.map((n) => n.toLowerCase())).size} genes found` + (r.missing.length ? `<br>not in panel: ${r.missing.join(", ")}` : "") + (r.excluded?.length ? `<br>excluded from all analyses: ${r.excluded.join(", ")}` : "");
    }
    S.sigName = ""; S.sigKey = ""; S.colorBy = S.colorBy === "__sigdual" ? "__sigdual" : "__sig"; $("colorby").value = S.colorBy; syncGeneUI(); refresh();
  });
  $("applysig2").addEventListener("click", async () => {   // signature B (typed genes)
    const names = $("sigtext2").value.split(/[\s,;]+/).filter(Boolean);
    if (!names.length) { S.sigGenes2 = []; $("sigstatus2").textContent = ""; }
    else {
      const r = await jget("/api/genes/resolve?names=" + encodeURIComponent(names.join(","))); S.sigGenes2 = r.found;
      $("sigstatus2").innerHTML = `${r.found.length} of ${new Set(names.map((n) => n.toLowerCase())).size} genes found` + (r.missing.length ? `<br>not in panel: ${r.missing.join(", ")}` : "") + (r.excluded?.length ? `<br>excluded from all analyses: ${r.excluded.join(", ")}` : "");
    }
    S.sigName2 = ""; S.sigKey2 = ""; S.colorBy = "__sigdual"; $("colorby").value = "__sigdual"; syncGeneUI(); refresh();
  });
  $("sigpreset2").addEventListener("change", (e) => { const sg = SIG.byId[e.target.value]; if (!sg) return; S.sigGenes2 = [...sg.genes]; S.sigName2 = sg.name; S.sigKey2 = ""; $("sigtext2").value = sg.genes.join(" ");
    $("sigstatus2").textContent = `${sg.genes.length} genes of "${sg.name}" loaded`; S.colorBy = "__sigdual"; $("colorby").value = "__sigdual"; syncGeneUI(); refresh(); });
  $("smooth").addEventListener("input", (e) => { $("smoothlab").textContent = +e.target.value ? e.target.value + " µm" : "off"; });
  $("smooth").addEventListener("change", (e) => { S.smooth = +e.target.value; refresh(); });
  $("zscope").addEventListener("change", (e) => { S.zscope = e.target.value; refresh(); });
  $("layer").addEventListener("change", (e) => { S.layer = e.target.value; refresh(); });
  const numIn = () => { const v = (id) => { const t = $(id).value.trim(); return t === "" || !Number.isFinite(+t) ? null : +t; }; S.numScale = { log: $("numlog").checked, cmin: v("numcmin"), cmax: v("numcmax"), fmin: v("numfmin"), fmax: v("numfmax") }; updateAll(); };
  for (const id of ["numlog", "numcmin", "numcmax", "numfmin", "numfmax"]) $(id).addEventListener("change", numIn);
  $("numreset").addEventListener("click", () => { for (const id of ["numcmin", "numcmax", "numfmin", "numfmax"]) $(id).value = ""; $("numlog").checked = false; numIn(); });
  $("clip").addEventListener("input", (e) => { S.clip = +e.target.value; $("cliplab").textContent = "p" + S.clip; updateAll(); });
  $("psize").addEventListener("input", (e) => { S.psize = +e.target.value; $("psizelab").textContent = S.psize + " µm"; render(); });
  $("shape").addEventListener("change", (e) => { S.shape = e.target.value; render(); });
  $("opacity").addEventListener("input", (e) => { S.opacity = +e.target.value; render(); });
  bindHE(); ifInit(); scanInit(); ptInit();
  $("cellsearch").addEventListener("input", debounce((e) => { S.cellSearch = e.target.value.trim(); refresh(); }));
  $("resetf").addEventListener("click", () => { S.filters = {}; S.exprFilters = []; S.hiddenCats = {}; S.cellSearch = ""; $("cellsearch").value = ""; refresh(); });
  $("fit").addEventListener("click", () => { fitView(); });
  const toggleFocus = () => {  // viewer-only mode: hide header + side panels, then let deck.gl re-measure its container
    const on = document.body.classList.toggle("focus"); $("focusbtn").textContent = on ? "Show panels" : "Hide panels"; $("focusbtn").classList.toggle("on", on);
    window.dispatchEvent(new Event("resize")); requestAnimationFrame(() => { window.dispatchEvent(new Event("resize")); updateScaleBar(); if (S.colors && deckgl) deckgl.setProps({ layers: layers(S.colors) }); });
  };
  $("focusbtn").addEventListener("click", toggleFocus);
  initOvLegend();
  const toggleLegend = () => { S.ovLegend = !S.ovLegend; $("legbtn").classList.toggle("on", S.ovLegend); updateOvLegend(); };
  $("legbtn").addEventListener("click", toggleLegend);
  window.addEventListener("keydown", (e) => {
    if (e.ctrlKey || e.metaKey || e.altKey || /^(INPUT|TEXTAREA|SELECT)$/.test(document.activeElement?.tagName || "")) return;
    if (e.key === "f" || e.key === "F") toggleFocus(); else if (e.key === "l" || e.key === "L") toggleLegend();
  });
  $("galbtn").addEventListener("click", () => { closeGlom(); openGallery(); });
  slideBind();
  $("glombtn").style.display = S.meta.cats.glom_id && S.meta.cats.glom_quality ? "" : "none";
  $("umapbtn").style.display = S.meta.umaps?.length ? "" : "none";
  if (S.meta.has_gloms === false) $("ghbtn").style.display = "none";   // no region_heterogeneity object in this project
  txInit(); hesInit(); greyInit();
  if (S.meta.embedding_label) $("umemb").querySelector("option[value=pca]").textContent = "PCA (10x, not batch-corrected)";   // projects whose embeddings are not Harmony-corrected (Xenium import)
  $("glombtn").addEventListener("click", openGlom);
  $("glomclose").addEventListener("click", closeGlom);
  $("umapbtn").addEventListener("click", () => toggleUmap()); $("umx").addEventListener("click", () => toggleUmap(false)); $("umclear").addEventListener("click", clearHighlight);
  const npc = S.meta.n_pcs || 0; if (!npc) [...$("umemb").options].forEach((o) => { if (o.value === "pca") o.remove(); });
  for (const [id, def] of [["umpcx", 1], ["umpcy", 2]]) { $(id).innerHTML = Array.from({ length: npc }, (_, i) => `<option value="${i + 1}">PC${i + 1}</option>`).join(""); $(id).value = def; }
  const umChange = () => { $("umpcs").hidden = $("umemb").value !== "pca"; updateUmap(true); };
  $("umfit").addEventListener("click", umFit); $("umemb").addEventListener("change", umChange); $("umpcx").addEventListener("change", umChange); $("umpcy").addEventListener("change", umChange); $("gzcanvas").addEventListener("click", glomZoomPick);
  initFloating($("umapwin"), $("umhead"), "xb_umap");
  $("ghbtn").addEventListener("click", openGH); $("ghclose").addEventListener("click", closeGH);
  for (const id of ["ghct", "ghcs", "ghclr", "ghlayer", "ghpx", "ghpy", "ghhide"]) $(id).addEventListener("change", renderGH);
  for (const id of ["ghsize", "ghop"]) $(id).addEventListener("input", () => { clearTimeout(GH.sT); GH.sT = setTimeout(renderGH, 60); });
  $("ghtop").addEventListener("change", renderGHHeat); for (const id of ["ghview", "ghsort", "ghdir"]) $(id).addEventListener("change", renderGHHeat); $("ghcsv").addEventListener("click", ghExportCSV); $("ghcv").addEventListener("pointermove", ghHover); $("ghcv").addEventListener("pointerleave", ghHoverEnd);
  for (const k of ["gene", "type", "state"]) $("ghtab-" + k).addEventListener("click", () => { GH.tab = k; renderGHHeat(); });
  $("viobtn").addEventListener("click", openVio); $("vioclose").addEventListener("click", closeVio); $("vioaddg").addEventListener("click", vioAddGene);
  $("viogene").addEventListener("keydown", (e) => { if (e.key === "Enter") vioAddGene(); }); $("viogene").addEventListener("input", debounce(async (e) => { const l = await jget("/api/genes/search?q=" + encodeURIComponent(e.target.value)); $("genelist").innerHTML = l.map((g) => `<option value="${g}">`).join(""); }));
  $("vioadds").addEventListener("click", () => { const v = $("viosig").value; if (v && !VIO.feats.some((f) => f.sig === v)) { VIO.feats.push({ sig: v }); vioControls(); renderVio(); } });
  $("vioshow").addEventListener("change", (e) => { VIO.level = e.target.value; if (VIO.level === "cell_state") VIO.cts.delete(-1); vioControls(); renderVio(); });
  $("vioscope").addEventListener("change", (e) => { VIO.scope = e.target.value; renderVio(); });
  $("viotype").addEventListener("change", (e) => { VIO.type = e.target.value; $("viology").parentElement.hidden = VIO.type !== "hist"; renderVio(); }); $("viology").addEventListener("change", (e) => { VIO.logy = e.target.checked; renderVio(); }); $("vionz").addEventListener("change", (e) => { VIO.nz = e.target.checked; renderVio(); });
  $("layer").addEventListener("change", () => { if (!$("vio").hidden) renderVio(); });
  $("vionit").addEventListener("change", (e) => { VIO.unit = e.target.value; renderVio(); });
  $("sigbtn").addEventListener("click", openSig); $("sigclose").addEventListener("click", closeSig); $("sigsearch").addEventListener("input", renderSigList);
  $("sigtab-explore").addEventListener("click", () => { SIG.tab = "explore"; renderSigPane(); }); $("sigtab-loc").addEventListener("click", () => { SIG.tab = "loc"; renderSigPane(); });
  loadSigList().catch(() => {}); $("sigpreset").addEventListener("change", (e) => { const sg = SIG.byId[e.target.value]; if (!sg) return; S.sigGenes = [...sg.genes]; S.sigName = sg.name; S.sigKey = ""; $("sigtext").value = sg.genes.join(" ");
    $("sigstatus").textContent = `${sg.genes.length} genes of "${sg.name}" loaded`; S.colorBy = "__sig"; $("colorby").value = "__sig"; syncGeneUI(); refresh(); });
  $("cmpbtn").addEventListener("click", openCmp); $("cmpclose").addEventListener("click", closeCmp);
  $("cmplevel").addEventListener("change", renderCmp); $("cmpincl").addEventListener("change", renderCmp);
  $("cmpview").addEventListener("click", () => { $("cmpview").classList.toggle("on"); $("cmpview").textContent = $("cmpview").classList.contains("on") ? "Bars" : "Table"; renderCmp(); });
  $("nbrbtn").addEventListener("click", openNbr); $("nbclose").addEventListener("click", closeNbr); $("nbcsv").addEventListener("click", exportNbrCSV);
  $("nbfield").addEventListener("change", renderNbr); $("nbview").addEventListener("change", drawNbr); $("nbgroup").addEventListener("change", drawNbr);
  $("nbdonbtn").addEventListener("click", (e) => { e.stopPropagation(); $("nbdonpop").hidden = !$("nbdonpop").hidden; });
  document.addEventListener("click", (e) => { if (!e.target.closest("#nbdonwrap")) $("nbdonpop").hidden = true; });
  $("nbrad").addEventListener("input", (e) => { $("nbradlab").textContent = e.target.value + " µm"; }); $("nbrad").addEventListener("change", renderNbr);
  const syncNbrMethod = () => {  // show only the inputs the chosen definition uses
    const m = $("nbmethod").value; $("nbkwrap").hidden = m !== "knn"; $("nbstwrap").hidden = m !== "delaunay";
    $("nbdistlab").textContent = m === "radius" ? "Radius" : m === "knn" ? "Max distance" : "Max edge length";
    const def = { radius: 25, knn: 50, delaunay: 50 }[m]; if (S.nbMethod !== m) { $("nbrad").value = def; $("nbradlab").textContent = def + " µm"; S.nbMethod = m; }
  };
  $("nbmethod").addEventListener("change", () => { syncNbrMethod(); renderNbr(); }); syncNbrMethod();
  $("nbk").addEventListener("change", renderNbr); $("nbsteps").addEventListener("change", renderNbr);
  $("nbhbtn").addEventListener("click", openNbh); $("nbhclose").addEventListener("click", closeNbh);
  $("nbhback").addEventListener("click", () => { closeNbh(); openNbr(); });
  $("nbhsample").addEventListener("change", () => { $("nbhn").value = $("nbhsample").value === "cur" ? 48 : 12; renderNbh(); }); $("nbhn").addEventListener("change", () => { S.nbhTouched = true; renderNbh(); }); $("nbhorder").addEventListener("change", renderNbh);
  $("nbhshuffle").addEventListener("click", () => { S.nbhSeed = 1 + Math.floor(Math.random() * 1e6); $("nbhorder").value = "random"; renderNbh(); });
  $("nbcounts").addEventListener("click", () => { S.nbCounts = !S.nbCounts; $("nbcounts").classList.toggle("on", S.nbCounts); drawNbr(); });
  $("gzclose").addEventListener("click", closeGlomZoom); gzZoomBind();
  $("gzsave").addEventListener("click", saveGlomPNG);
  $("gzctx").addEventListener("change", drawGlomZoom); $("gzcanvas").addEventListener("pointermove", glomZoomHover); $("gzcanvas").addEventListener("pointerleave", () => { $("nbrtip").hidden = true; });
  $("gzopen").addEventListener("click", () => { const v = S.gz; closeGlomZoom(); if (v) openGlomIn(v.s.id, v.cells, v.poly); });
  $("glomscope").addEventListener("change", renderGlom);
  $("glomonly").addEventListener("change", renderGlom); $("glomrz").addEventListener("click", () => resetZoom(glomCfg)); $("glomsync").addEventListener("change", () => resetZoom(glomCfg));
  let glT; window.addEventListener("resize", () => { clearTimeout(glT); glT = setTimeout(() => { if (!$("glom").hidden && $("glomonly").checked) renderGlom(); }, 250); });
  $("glomqc").addEventListener("change", renderGlom);
  $("galclose").addEventListener("click", closeGallery);
  $("galsub").addEventListener("change", renderGallery); $("galsz").addEventListener("change", renderGallery); $("galonly").addEventListener("change", () => { const on = $("galonly").checked; if (on) { S.galPrevSz = $("galsz").value; $("galsz").value = "0"; } else if (S.galPrevSz) $("galsz").value = S.galPrevSz; renderGallery(); });
  let gT; window.addEventListener("resize", () => { clearTimeout(gT); gT = setTimeout(() => { if (!$("gallery").hidden && $("galonly").checked) renderGallery(); }, 250); });
  $("galrz").addEventListener("click", () => resetZoom(galCfg)); $("galsync").addEventListener("change", () => resetZoom(galCfg));
  $("galclr").addEventListener("click", () => { (S.galSel || (S.galSel = new Set())).clear(); renderGallery(); });
  window.addEventListener("keydown", (e) => { if (e.key === "Escape") { if (!$("vio").hidden && $("pngmenu").hidden && $("glomzoom").hidden) closeVio(); else if (!$("sigwin").hidden && $("pngmenu").hidden && $("glomzoom").hidden) closeSig(); else if (!$("cmp").hidden && $("pngmenu").hidden && $("glomzoom").hidden) closeCmp(); else if (!$("nbh").hidden && $("pngmenu").hidden && $("glomzoom").hidden) closeNbh(); else if (!$("nbr").hidden && $("pngmenu").hidden && $("glomzoom").hidden) closeNbr(); else if (!$("pngmenu").hidden) $("pngmenu").hidden = true; else if (!$("glomzoom").hidden) closeGlomZoom(); else if (!$("glomcmp").hidden) closeGlomCmp(); else if (!$("rgn").hidden) closeRgn(); else { closeSlide(); closeGallery(); closeGlom(); } } });
  $("clearsel").addEventListener("click", () => { S.selected = new Set(); $("cell").innerHTML = "Click a cell."; clearHighlight(); updateAll(); });
  $("boxsel").addEventListener("click", () => { S.boxMode = !S.boxMode; $("boxsel").classList.toggle("on", S.boxMode); deckgl.setProps({ controller: S.boxMode ? { dragPan: false } : { dragPan: true, scrollZoom: true, doubleClickZoom: true } }); });
  $("compfield").addEventListener("change", (e) => { S.compField = e.target.value; refresh(); });
  $("legsearch").addEventListener("input", updateLegend);
  $("expsel").addEventListener("click", () => exportCSV([...S.selected].filter((i) => S.vis[i]), `selected_${S.sample}.csv`));
  $("expfilt").addEventListener("click", () => exportCSV(Array.from(S.visIdx), `filtered_${S.sample}.csv`));
  $("png").addEventListener("click", (e) => { e.stopPropagation(); onSavePNG(); });
  const pm = $("pngmenu");
  $("pngvis").addEventListener("click", () => { pm.hidden = true; saveWindowPNG(pm.dataset.kind, false); });
  $("pngall").addEventListener("click", () => { pm.hidden = true; saveWindowPNG(pm.dataset.kind, true); });
  document.addEventListener("click", (e) => { if (!pm.hidden && !pm.contains(e.target)) pm.hidden = true; });
  $("showpath").addEventListener("change", async (e) => { S.showPath = e.target.checked; S.pathInfo = null; syncPathUI(); if (S.showPath) await loadPaths(); syncPathUI(); render(); });
  bindBox(); bindNewUI();
  window.addEventListener("resize", () => updateScaleBar());
}

// Pathologist region outlines. Coordinates are in the cells' system, which the viewer draws with y flipped; each annotation keeps the pen colour
// the pathologist used (line_color_hex); "negative" regions (excluded areas) are grey.
async function loadPaths() {
  S.pathLayer = null; S.pathInfo = null;
  const smp = S.meta.samples.find((x) => x.id === S.sample);
  if (!smp?.has_paths) { S.pathInfo = { none: true }; return; }
  try {
    const g = await jget(`/api/paths/${encodeURIComponent(S.sample)}`);
    const flip = (c) => (typeof c[0] === "number" ? [c[0], -c[1]] : c.map(flip));
    const feats = g.features.map((f) => ({ ...f, geometry: { ...f.geometry, coordinates: flip(f.geometry.coordinates) } }));
    const rgb = (f) => (f.properties.negative_roa ? [150, 160, 175] : hex2rgb(regionColor(f.properties.annotation_name, f.properties.line_color_hex)));
    const legend = new Map();
    for (const f of g.features) { const k = (f.properties.annotation_name || "").trim() + (f.properties.negative_roa ? " (excluded area)" : ""), e = legend.get(k) || { name: k, color: f.properties.negative_roa ? "#96a0af" : regionColor(f.properties.annotation_name, f.properties.line_color_hex), n: 0 }; e.n++; legend.set(k, e); }
    S.pathInfo = { n: g.features.length, tissues: g.tissues, legend: [...legend.values()].sort((a, b) => b.n - a.n) };
    const data = { type: "FeatureCollection", features: feats };
    const selNames = () => { const st = S.filters.path_region; return st?.size ? new Set([...st].map((i) => catMeta("path_region").levels[i])) : null; };
    S.pathLayer = () => { const sel = selNames(); return new deck.GeoJsonLayer({ id: "paths", data, stroked: true, filled: false, pickable: false, lineWidthUnits: "pixels", getLineWidth: 2, lineWidthMinPixels: 1.5, getLineColor: (f) => [...rgb(f), sel && !sel.has(f.properties.annotation_name) ? 55 : 235], updateTriggers: { getLineColor: [S.filters.path_region && [...S.filters.path_region].join()] } }); };
  } catch (e) { S.pathInfo = { error: e.message }; }
}
const pathLegendHTML = () => (S.showPath && S.pathInfo?.legend ? S.pathInfo.legend.map((e) => `<div class="opt"><span class="sw" style="background:${e.color}"></span><span class="l">${e.name}</span><span class="n">${e.n}</span></div>`).join("") : "");
function syncPathUI() {  // status line + legend under the checkbox, and the same legend in the floating legend
  const i = S.pathInfo, el = $("pathstat");
  el.textContent = !S.showPath ? "" : !i ? "loading…" : i.none ? "no pathologist annotation was matched to this sample" : i.error ? "could not load: " + i.error : `${i.n} regions from ${i.tissues.join(", ")}`;
  updateOvLegend();
}

// ---------------------------------------------------------------- Glom heterogeneity: gene / signature settings are the main window's
function syncMainGeneUI() {   // mirror S into the main panel after the Glom heterogeneity window changed it
  $("colorby").value = S.colorBy; $("geneA").value = S.geneA; $("geneB").value = S.geneB; $("sigtext").value = S.sigGenes.join(" "); $("sigtext2").value = S.sigGenes2.join(" ");
  $("sigstatus").textContent = S.sigGenes.length ? `${S.sigGenes.length} genes${S.sigName ? ` of "${S.sigName}"` : ""} loaded` : "";
  $("sigstatus2").textContent = S.sigGenes2.length ? `${S.sigGenes2.length} genes${S.sigName2 ? ` of "${S.sigName2}"` : ""} loaded` : "";
  syncGeneUI();
}
async function ghPush() {
  const mode = $("ghcol").value, MAP = { gene: "__gene", dual: "__dual", signature: "__sig", sigdual: "__sigdual" }, msg = $("ghmsg");
  if (!MAP[mode]) return renderGH();
  msg.textContent = "";
  const geneOf = async (id) => { const v = $(id).value.trim(); if (!v) return ""; const hit = (await jget("/api/genes/search?q=" + encodeURIComponent(v) + "&limit=1"))[0];
    if (hit && hit.toLowerCase() === v.toLowerCase()) { $(id).value = hit; return hit; } msg.textContent = `“${v}” is not on the panel`; return ""; };
  const sigOf = async (selId, txtId) => { const t = $(txtId).value.split(/[\s,;]+/).filter(Boolean);
    if (t.length) { const r = await jget("/api/genes/resolve?names=" + encodeURIComponent(t.join(","))); if (r.missing.length) msg.textContent = `not on the panel: ${r.missing.join(", ")}`; return { genes: r.found, name: "" }; }
    const sg = SIG.byId[$(selId).value]; return sg ? { genes: [...sg.genes], name: sg.name } : { genes: [], name: "" }; };
  try {
    if (mode === "gene" || mode === "dual") { S.geneA = await geneOf("ghgene"); if (mode === "dual") S.geneB = await geneOf("ghgene2"); }
    else { const a = await sigOf("ghsig", "ghsigtxt"); S.sigGenes = a.genes; S.sigName = a.name; S.sigKey = "";
      if (mode === "sigdual") { const b = await sigOf("ghsig2", "ghsigtxt2"); S.sigGenes2 = b.genes; S.sigName2 = b.name; S.sigKey2 = ""; } }
  } catch (e) { msg.textContent = e.message.replace(/^.*?: /, "").slice(0, 160); return; }
  S.colorBy = MAP[mode]; syncMainGeneUI(); refresh();
}

// ---------------------------------------------------------------- tick glomeruli (Glomeruli / Glom heterogeneity windows) and compare them
S.glomSel = new Map();   // "sample|glom id" -> { sample, label }
function updateGlomSelUI() {
  const n = S.glomSel.size; $("glomn").textContent = n; $("glomcmpbtn").textContent = `Compare selected (${n})`; $("ghcmp").textContent = `Compare (${n})`; $("glomselclr").disabled = !n;
  document.querySelectorAll("#glomgrid .gtile").forEach((t) => { const on = S.glomSel.has(t.dataset.key); t.classList.toggle("sel", on); const c = t.querySelector(".gsel"); if (c) c.checked = on; });
}
function glomTick(key, rec, on) { on ? S.glomSel.set(key, rec) : S.glomSel.delete(key); updateGlomSelUI(); if (!$("glomcmp").hidden) renderGlomCmp(); }
function glomUntickAll() { S.glomSel.clear(); updateGlomSelUI(); if (!$("glomcmp").hidden) renderGlomCmp(); if (!$("glom").hidden && $("glomonly").checked) renderGlom(); }
function openGlomCmp() {
  if (!S.glomSel.size) { alert("Tick one or more glomeruli first: the checkbox on each tile in the Glomeruli window, or “Add to comparison” in ROI heterogeneity."); return; }
  $("glomcmp").hidden = false; renderGlomCmp();
}
function closeGlomCmp() { $("glomcmp").hidden = true; }
let gcToken = 0;
async function renderGlomCmp() {
  if ($("glomcmp").hidden) return;
  const tok = ++gcToken, sel = [...S.glomSel.values()], size = +$("gcsize").value, ctxOn = $("gcctx").checked, own = $("gcown").checked, grid = $("gcgrid");
  S.gcCards = [];
  if (!sel.length) { grid.innerHTML = `<div class="muted">No glomeruli ticked. Go back and tick some.</div>`; $("gcleg").innerHTML = ""; $("gcinfo").textContent = ""; return; }
  $("gcinfo").textContent = `${sel.length} glomeruli · colour: ${colorName()}${winNote()} · click a tile to zoom`;
  const gidLv = catMeta("glom_id").levels, qcLv = catMeta("glom_quality").levels, samples = [...new Set(sel.map((x) => x.sample))], data = {};
  try {
    await Promise.all(samples.map(async (sid) => { const [b, gid, gq, poly, feats, ct] = await Promise.all([winBundle(sid), galGet(sid, "glom_id"), galGet(sid, "glom_quality"), galPoly(sid), galPaths(sid), galGet(sid, "cell_type")]); data[sid] = { b, gid, gq, poly, feats, ct }; }));
  } catch (e) { grid.textContent = "Error: " + e.message; return; }
  if (tok !== gcToken) return;
  const cards = sel.map((x) => { const d = data[x.sample], code = gidLv.indexOf(x.label), cells = []; for (let i = 0; i < d.b.n; i++) if (d.gid[i] === code) cells.push(i);
    const qc = cells.length ? qcLv[d.gq[cells[0]]] : ""; return { ...x, key: `${x.sample}|${x.label}`, d, cells, qc, s: S.meta.samples.find((z) => z.id === x.sample), outline: qc === "passed qc" ? glomOutline(d.feats, x.label) : null }; }).filter((c) => c.cells.length);
  const restrict = (c) => { const p = new Uint8Array(c.d.b.n); for (const i of c.cells) if (c.d.b.pass[i]) p[i] = 1; return { ...c.d.b, pass: p }; };
  const sc = S.meta.cats[S.colorBy] ? makeScale([]) : makeScale(own ? cards.map(restrict) : samples.map((sid) => data[sid].b));
  S.gcCards = cards; S.gcSc = sc;
  const types = new Map();
  grid.replaceChildren();
  for (const c of cards) {
    const el = mk("div", "gcard"), cv = mk("canvas"), comp = new Map(); cv.width = cv.height = size; cv.style.width = cv.style.height = size + "px"; el.style.width = size + 20 + "px";
    for (const i of c.cells) { const t = c.d.ct[i]; comp.set(t, (comp.get(t) || 0) + 1); types.set(t, (types.get(t) || 0) + 1); }
    c.comp = [...comp.entries()].sort((a, b) => b[1] - a[1]);
    const nPass = c.cells.filter((i) => c.d.b.pass[i]).length, A = c.d.b.A; let mean = "";
    if (A && sc.info.type !== "cat") {
      const mn = (V) => { let sum = 0, n = 0; for (const i of c.cells) if (c.d.b.pass[i] && V[i] === V[i]) { sum += V[i]; n++; } return n ? (sum / n).toFixed(3) : "n/a"; }, B = c.d.b.B;
      mean = ` · mean ${B ? `A ${mn(A)} · B ${mn(B)}` : "score " + mn(A)}${S.zscope !== "none" ? " (z)" : ""}`;
    }
    const rm = mk("span", "rm", "✕"); rm.title = "remove from the comparison"; rm.addEventListener("click", () => glomTick(c.key, { sample: c.sample, label: c.label }, false));
    const gm = mk("div", "gm"), chip = mk("span", "chip " + (c.qc === "passed qc" ? "pass" : "fail"), (c.qc === "passed qc" ? "✓ " : "✕ ") + c.qc);
    gm.append(chip, `${dlabel(c.s.id)} · ${subtypeOf(c.s)}`, mk("br"), `${fmt(c.cells.length)} cells${nPass < c.cells.length ? ` (${fmt(nPass)} pass the filters)` : ""}${mean}`);
    const bar = mk("div", "cbar2"); for (const [t, n] of c.comp) { const seg = mk("i"); seg.style.width = (100 * n / c.cells.length) + "%"; seg.style.background = t === MISSING ? rgbStr(GRAY) : catColor("cell_type", t); seg.title = `${t === MISSING ? "no cell type" : catMeta("cell_type").levels[t]}: ${n} (${(100 * n / c.cells.length).toFixed(0)}%)`; bar.append(seg); }
    el.append(rm, cv, mk("div", "gid", c.label), gm, bar); grid.append(el); c.cv = cv; cv._hw = "gcmp";
    const cv2 = hePane(el, cv, "gcmp", c.s.id);   // H&E beside the cells: a second canvas, wheel = zoom, drag = pan
    if (cv2) {
      cv2.width = cv2.height = size; cv2.style.width = cv2.style.height = size + "px"; el.style.width = 2 * size + 26 + "px"; el.dataset.id = c.key; el.append(mk("span", "zb")); c.views = {};
      const B = glomBox(c.cells, c.d.poly, 12, c.outline), rec = el._dt = { bb: { x0: B.x0, y0: B.y0, w: B.x1 - B.x0, h: B.y1 - B.y0 }, draw: (vw, vh) => { drawGlom(cv, c.d.b, sc, c.d.poly, c.cells, ctxOn, -1, c.outline, null, vw); drawGlom(cv2, c.d.b, sc, c.d.poly, c.cells, ctxOn, -1, c.outline, null, HEW.gcmp.link ? vw : vh); } };
      attachZoom(el, rec, { views: () => c.views, sync: () => false, tiles: () => [el] }); rec.draw(null, null);
    } else drawGlom(cv, c.d.b, sc, c.d.poly, c.cells, ctxOn, -1, c.outline);
    cv.addEventListener("click", () => el._moved ? (el._moved = false) : openGlomZoom({ src: "gcmp", s: c.s, b: c.d.b, sc, poly: c.d.poly, cells: c.cells, label: c.label, qc: c.qc, outline: c.outline, note: c.outline ? "outline = pathologist \"All gloms\" boundary" : undefined }));
  }
  const tl = catMeta("cell_type").levels;
  $("gcleg").innerHTML = (sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info) + (own ? '<div class="small muted">colour scale fitted to the ticked gloms</div>' : ""))
    + `<div class="small muted" style="margin-top:6px">Bar under each tile = cell-type composition of the glom:</div><div class="nlegrow">` + [...types.entries()].sort((a, b) => b[1] - a[1]).map(([t]) => `<span><span class="k" style="background:${t === MISSING ? rgbStr(GRAY) : catColor("cell_type", t)}"></span>${t === MISSING ? "no cell type" : tl[t]}</span>`).join("") + `</div>`;
}
function saveGlomCmpPNG() {
  const cards = S.gcCards || []; if (!cards.length) return;
  const size = +$("gcsize").value, pad = 16, cw = size + 16, ch = size + 8 + 80, cols = Math.max(1, Math.min(cards.length, Math.floor(2200 / (cw + pad)))), rows = Math.ceil(cards.length / cols), W = pad + cols * (cw + pad), info = S.gcSc.info;
  let items = null;
  if (info.type === "cat") { const cnt = new Map(); let dim = 0; for (const c of cards) for (const i of c.cells) { if (!c.d.b.pass[i]) { dim++; continue; } const k = c.d.b.codes[i]; cnt.set(k, (cnt.get(k) || 0) + 1); }
    items = [...cnt.entries()].sort((a, b) => b[1] - a[1]).map(([k, n]) => [k === MISSING ? "missing" : display(S.colorBy, k), k === MISSING ? rgbStr(GRAY) : catColor(S.colorBy, k), n]); if (dim) items.push(["filtered out", "rgba(105,110,122,0.6)", dim]); }
  const lh = pngLegend(document.createElement("canvas").getContext("2d"), info, items, pad, 0, W - 2 * pad, 1, 6), HH = 74 + lh + 14, H = HH + rows * (ch + pad) + pad;
  const out = document.createElement("canvas"); out.width = W; out.height = H; const x = out.getContext("2d");
  x.fillStyle = "#0f1218"; x.fillRect(0, 0, W, H); x.textBaseline = "alphabetic"; x.fillStyle = "#e3e8ef"; x.font = `600 22px ${FONT}`; x.fillText(`Glomeruli comparison (${cards.length})`, pad, 32);
  x.font = `13px ${FONT}`; x.fillStyle = "#8b95a5"; x.fillText(`colour: ${colorName()}${winNote()}${$("gcown").checked && info.type !== "cat" ? " · colour scale fitted to these gloms" : ""}`, pad, 54);
  pngLegend(x, info, items, pad, 80, W - 2 * pad, 1, 6);
  cards.forEach((c, k) => {
    const cx = pad + (k % cols) * (cw + pad), cy = HH + Math.floor(k / cols) * (ch + pad);
    x.fillStyle = "#141922"; x.fillRect(cx, cy, cw, ch); x.strokeStyle = "#29303b"; x.lineWidth = 2; x.strokeRect(cx, cy, cw, ch);
    x.drawImage(c.cv, cx + 8, cy + 8, size, size);
    x.font = `600 13px ${FONT}`; x.fillStyle = "#d7dde6"; x.fillText(c.label, cx + 8, cy + size + 26);
    const chipTxt = (c.qc === "passed qc" ? "✓ " : "✕ ") + c.qc; x.font = `600 11px ${FONT}`; const cwid = x.measureText(chipTxt).width + 12;
    x.fillStyle = c.qc === "passed qc" ? STATUS.good : STATUS.critical; x.fillRect(cx + 8, cy + size + 33, cwid, 16); x.fillStyle = c.qc === "passed qc" ? "#0b0b0b" : "#fff"; x.fillText(chipTxt, cx + 14, cy + size + 45);
    x.font = `11.5px ${FONT}`; x.fillStyle = "#8b95a5"; x.fillText(`${dlabel(c.s.id)} · ${subtypeOf(c.s)} · ${fmt(c.cells.length)} cells`, cx + 8 + cwid + 8, cy + size + 45);
    let bx = cx + 8; for (const [t, n] of c.comp) { const w = (size * n) / c.cells.length; x.fillStyle = t === MISSING ? rgbStr(GRAY) : catColor("cell_type", t); x.fillRect(bx, cy + size + 56, w, 12); bx += w; }
  });
  out.toBlob((blob) => { const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = "glomeruli_comparison.png"; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }, "image/png");
}


// ---------------------------------------------------------------- Region neighbourhoods window
// One tile per region object (e.g. each sclerotic glomerulus): its pathologist outline, the cells inside it and the cells within a distance of its outline,
// all coloured by the current colour-by. Cells that lie inside another glomerulus are not part of the ring.
function pointInRings(x, y, rings) {   // even-odd rule over all rings (holes included)
  let c = false;
  for (const r of rings) for (let i = 0, j = r.length - 1; i < r.length; j = i++) { const x1 = r[i][0], y1 = r[i][1], x2 = r[j][0], y2 = r[j][1]; if (y1 > y !== y2 > y && x < ((x2 - x1) * (y - y1)) / (y2 - y1) + x1) c = !c; }
  return c;
}
function distToRings(x, y, rings) {   // distance (um) from a point to the nearest outline segment
  let best = Infinity;
  for (const r of rings) for (let i = 0, j = r.length - 1; i < r.length; j = i++) {
    const ax = r[j][0], ay = r[j][1], dx = r[i][0] - ax, dy = r[i][1] - ay, l2 = dx * dx + dy * dy; let t = l2 ? ((x - ax) * dx + (y - ay) * dy) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
    const d = Math.hypot(x - (ax + t * dx), y - (ay + t * dy)); if (d < best) best = d;
  }
  return best;
}
function rnbShape(f) {
  const rings = ringsOf(f); let x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18, cx = 0, cy = 0, n = 0;
  for (const [x, y] of rings[0]) { if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; cx += x; cy += y; n++; }
  return { f, rings, bb: { x0, x1, y0, y1 }, cen: [cx / n, cy / n] };
}
const rnbInBB = (bb, x, y) => x >= bb.x0 && x <= bb.x1 && y >= bb.y0 && y <= bb.y1;
function rnbCells(b, inst, others, d) {   // -> { inside, ring } cell indices
  const bb = inst.bb, inside = [], ring = [];
  for (let i = 0; i < b.n; i++) {
    const x = b.xy[i * 2], y = b.xy[i * 2 + 1]; if (x < bb.x0 - d || x > bb.x1 + d || y < bb.y0 - d || y > bb.y1 + d) continue;
    if (rnbInBB(bb, x, y) && pointInRings(x, y, inst.rings)) { inside.push(i); continue; }
    if (distToRings(x, y, inst.rings) > d) continue;
    let oth = false; for (const o of others) if (rnbInBB(o.bb, x, y) && pointInRings(x, y, o.rings)) { oth = true; break; }
    if (!oth) ring.push(i);
  }
  return { inside, ring };
}
function openRnb() {
  closeAllWindows();
  const names = catMeta("path_region").levels.slice(0, -1), subs = sortedSubtypes();
  $("rnbsel").innerHTML = names.map((n) => `<option value="${n}">${n === "All gloms" ? "Non-sclerotic gloms (pathologist “All gloms”)" : n}</option>`).join("");
  $("rnbsel").value = S.rnbName && names.includes(S.rnbName) ? S.rnbName : names.includes("Sclerotic Gloms") ? "Sclerotic Gloms" : names[0];
  $("rnbscope").innerHTML = `<option value="all">All donors</option><option value="cur">This donor (${dlabel(S.sample)})</option>` + subs.map((t) => `<option value="sub:${t}">Disease subtype: ${t}</option>`).join("");
  $("rnbscope").value = S.rnbScope || "all"; if ($("rnbscope").value !== (S.rnbScope || "all")) $("rnbscope").value = "all";
  $("rnb").hidden = false; renderRnb();
}
let rnbToken = 0;
const rnbCfg = { views: () => S.rnbView || (S.rnbView = {}), sync: () => $("rnbsync").checked, tiles: () => $("rnbgrid").querySelectorAll(".gtile") };
async function rnbTiles(s, name, d) {   // the region objects of one tissue with their cells -> { poly, tiles: [{ key, label, it, inside, ring, cells, outline }] } (null when the tissue has no outlines)
  const [b, poly, feats] = await Promise.all([winBundle(s.id), galPoly(s.id), galPaths(s.id)]); if (!feats) return null;
  const live = feats.filter((f) => !f.properties.negative_roa), gl = live.filter((f) => /glom/i.test(f.properties.annotation_name)).map(rnbShape), num = (x) => +x.f.properties.region_id.match(/#(\d+)/)?.[1] || 0;
  const inst = live.filter((f) => f.properties.annotation_name === name).map(rnbShape).sort((p, q) => num(p) - num(q)), seen = new Map(), tiles = [];
  for (const it of inst) {
    const others = gl.filter((o) => o.f !== it.f && !pointInRings(o.cen[0], o.cen[1], it.rings) && !pointInRings(it.cen[0], it.cen[1], o.rings));   // other glomeruli (not this one under a second outline)
    const { inside, ring } = rnbCells(b, it, others, d), cells = inside.concat(ring); if (!cells.some((i) => b.pass[i])) continue;
    const rid = it.f.properties.region_id, nth = (seen.get(rid) || 0) + 1; seen.set(rid, nth); const label = nth > 1 ? `${rid} (${nth})` : rid;
    tiles.push({ key: `${s.id}|${label}`, label, it, inside, ring, cells, outline: { rings: it.rings, color: regionColor(name === "All gloms" ? "Non-sclerotic Gloms" : name, it.f.properties.line_color_hex) } });
  }
  return { b, poly, tiles };
}
async function renderRnb() {
  if ($("rnb").hidden) return;
  const tok = ++rnbToken, name = $("rnbsel").value, scope = $("rnbscope").value, d = Math.max(1, +$("rnbd").value || 50), grid = $("rnbgrid"), sel = S.rnbSel || (S.rnbSel = new Map()), only = $("rnbonly").checked, size = +$("rnbsz").value;
  S.rnbName = name; S.rnbScope = scope; $("rnbn").textContent = sel.size; $("rnbclr").disabled = !sel.size; $("rnbscope").disabled = $("rnbsz").disabled = only; $("rnbzoomctl").hidden = !only; grid.classList.toggle("cmpgrid", only); if (!only) grid.removeAttribute("style");
  const rank = (t) => sortedSubtypes().indexOf(t), pool = S.meta.samples.filter((s) => s.has_paths);
  const samples = (only ? pool.filter((s) => [...sel.values()].some((r) => r.sample === s.id)) : pool.filter((s) => scope === "all" || (scope === "cur" ? s.id === S.sample : subtypeOf(s) === scope.slice(4)))).sort((x, y) => rank(subtypeOf(x)) - rank(subtypeOf(y)));
  const totals = { n: 0, donors: new Set() }, subInfo = new Map();
  const info = (done) => ($("rnbinfo").textContent = only ? `${totals.n} ticked “${name}” region${totals.n === 1 ? "" : "s"} · cells inside plus within ${d} µm · colour: ${colorName()}${winNote()} · window split per region · wheel = zoom, drag = pan, double-click = reset, click a name to open it${done ? "" : " · loading…"}`
    : `${totals.n} “${name}” region${totals.n === 1 ? "" : "s"} in ${donTxt(totals.donors.size)} · cells inside the outline plus cells within ${d} µm of it (cells inside other glomeruli left out) · white lines separate donors · colour: ${colorName()}${winNote()}${done ? "" : " · loading…"}`);
  $("rnbleg").innerHTML = ""; info(false);
  if (only && !samples.length) { grid.innerHTML = `<div class="muted">No region is ticked. Untick “only ticked” and tick the tiles to compare.</div>`; $("rnbinfo").textContent = ""; return; }
  const bp = samples.map((s) => winBundle(s.id)), sc = await winScale(bp, () => rnbToken, tok).catch((e) => { console.error(e); return null; });
  if (!sc || tok !== rnbToken) return;
  $("rnbleg").innerHTML = sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info);
  const tick = (r, on, s) => { on ? sel.set(r.key, { sample: s.id, label: r.label }) : sel.delete(r.key); $("rnbn").textContent = sel.size; $("rnbclr").disabled = !sel.size; if (only) renderRnb(); };
  if (only) {   // compare view: one equal chunk per ticked region, wheel = zoom, drag = pan
    const items = [];
    await Promise.all(samples.map(async (s) => { const r = await rnbTiles(s, name, d); if (!r) return; for (const t of r.tiles) if (sel.has(t.key)) items.push({ ...t, s, b: r.b, poly: r.poly }); }));
    if (tok !== rnbToken) return; const order = [...sel.keys()]; items.sort((x, y) => order.indexOf(x.key) - order.indexOf(y.key)); totals.n = items.length; info(true);
    grid.innerHTML = items.map((r) => `<div class="gtile sel" data-id="${r.key}"><input type="checkbox" class="gsel" title="untick to remove from this view" checked><span class="zb"></span><canvas width="300" height="300"></canvas>
      <div class="cap"><b class="gid">${r.label}</b><div class="gm">${dlabel(r.s.id)} · ${subtypeOf(r.s)}<br>${fmt(r.inside.length)} cells inside · ${fmt(r.ring.length)} within ${d} µm</div></div></div>`).join("");
    const dims = cmpLayout("rnb", "rnbgrid", Math.max(1, items.length), 70);
    items.forEach((r) => {
      const t = grid.querySelector(`.gtile[data-id="${CSS.escape(r.key)}"]`), cv = t.querySelector("canvas"), B = glomBox(r.cells, r.poly, 12, r.outline); cv._hw = "rnb";
      const cv2 = hePane(t, cv, "rnb", r.s.id), dm = cv2 ? hePairDims(dims, "rnb") : dims;
      const rec = t._dt = { bb: { x0: B.x0, y0: B.y0, w: B.x1 - B.x0, h: B.y1 - B.y0 }, last: null, last2: null, draw: (vw, vh) => { rec.last = vw; rec.last2 = vh; drawGlom(cv, r.b, sc, r.poly, r.cells, false, -1, r.outline, () => heRedo(t), vw, dm); if (cv2) drawGlom(cv2, r.b, sc, r.poly, r.cells, false, -1, r.outline, () => heRedo(t), HEW.rnb.link ? vw : vh, dm); } };
      const gk = rnbCfg.sync() ? "__all" : r.key, v = rnbCfg.views()[gk] || null; rec.draw(v, rnbCfg.views()[gk + "#h"] || null); attachZoom(t, rec, rnbCfg); if (v && v.z > 1.01) t.querySelector(".zb").textContent = "×" + v.z.toFixed(1);
      t.querySelector(".gsel").addEventListener("change", (e) => tick(r, e.target.checked, r.s));
      t.querySelector(".gid").addEventListener("click", () => openGlomZoom({ src: "rnb", s: r.s, b: r.b, sc, poly: r.poly, cells: r.cells, label: r.label, qc: null, outline: r.outline, note: `${name} outline · cells inside + within ${d} µm` }));
    });
    return;
  }
  // one row per disease subtype; tissues of the same donor sit together and a white line separates donors
  const subs = [...new Set(samples.map(subtypeOf))]; let prevDonor = null;
  grid.innerHTML = subs.map((t) => `<div class="rnbsec" data-t="${t}" hidden><div class="sechead"></div><div class="rnbrow">` + samples.filter((s) => subtypeOf(s) === t).map((s) => { const nw = s.donor !== prevDonor; prevDonor = s.donor; return `<div class="dgrp${nw ? " dnew" : ""}" data-s="${s.id}"></div>`; }).join("") + `</div></div>`).join("");
  let pending = samples.length;
  await Promise.all(samples.map(async (s, k) => {
    try {
      const r = await rnbTiles(s, name, d); if (tok !== rnbToken || !r || !r.tiles.length) return;
      const slot = grid.querySelector(`.dgrp[data-s="${CSS.escape(s.id)}"]`), sec = slot.closest(".rnbsec"), t0 = subtypeOf(s); sec.hidden = false; totals.donors.add(s.donor);
      const si = subInfo.get(t0) || subInfo.set(t0, { n: 0, donors: new Set() }).get(t0); si.donors.add(s.donor);
      for (const t of r.tiles) {
        totals.n++; si.n++;
        const el = document.createElement("div"); el.className = "gtile" + (sel.has(t.key) ? " sel" : ""); el.style.width = size + 12 + "px"; el.title = `${t.label} — click to zoom; tick the box to compare`;
        el.innerHTML = `<input type="checkbox" class="gsel" title="tick to compare this region" ${sel.has(t.key) ? "checked" : ""}><canvas width="${size}" height="${size}" style="width:${size}px;height:${size}px"></canvas><div class="gid">${t.label}</div><div class="gm">${dlabel(s.id)}<br>${subtypeOf(s)}<br>${fmt(t.inside.length)} cells inside · ${fmt(t.ring.length)} within ${d} µm</div>`;
        el.addEventListener("click", () => openGlomZoom({ src: "rnb", s, b: r.b, sc, poly: r.poly, cells: t.cells, label: t.label, qc: null, outline: t.outline, note: `${name} outline · cells inside + within ${d} µm` }));
        const cb = el.querySelector(".gsel"); cb.addEventListener("click", (e) => e.stopPropagation()); cb.addEventListener("change", () => { el.classList.toggle("sel", cb.checked); tick(t, cb.checked, s); });
        slot.appendChild(el); const cv = el.querySelector("canvas"); cv._hw = "rnb"; drawGlom(cv, r.b, sc, r.poly, t.cells, false, -1, t.outline);
      }
      sec.querySelector(".sechead").textContent = `${t0} · ${si.n} region${si.n === 1 ? "" : "s"} · ${donTxt(si.donors.size)}`;
    } catch (e) { console.error(e); }
    if (tok === rnbToken) { pending--; info(pending === 0); }
  }));
  if (tok === rnbToken) { info(true); if (!totals.n) grid.innerHTML = `<div class="muted">No “${name}” outlines in this selection (or none with cells passing the filters).</div>`; }
}

// ---------------------------------------------------------------- Regions window: one region (not a glomerulus) in every sample, with its boundary
const GLOMISH = /^(glomerulus|glom tessera|all gloms|poorly preserved gloms|intraglom crescents|gloms w\/)/i;
const rgnKeys = () => ["path_region", "region_lvl1", "region_lvl2", "region_lvl3"].filter((k) => S.meta.cats[k]);
function rgnLevels(k) { const L = catMeta(k).levels; return L.map((_, i) => i).filter((i) => !(k === "path_region" && i === L.length - 1) && !GLOMISH.test(String(L[i]).trim())); }
function rgnDefault(k) { const L = catMeta(k).levels, ls = rgnLevels(k); return ls.find((i) => /^cortex$/i.test(String(L[i]).trim())) ?? ls[0] ?? -1; }
function rgnInit() { if (!S.meta.cats[RG.src]) RG.src = rgnKeys()[0]; if (RG.src && RG.level < 0) RG.level = rgnDefault(RG.src); }
function rgnFillSel() { const k = RG.src, L = catMeta(k).levels, ls = rgnLevels(k); $("rgnsel").innerHTML = ls.map((i) => `<option value="${i}">${L[i]}</option>`).join(""); if (!ls.includes(RG.level)) RG.level = rgnDefault(k); $("rgnsel").value = RG.level; }
function openRgn() {
  closeAllWindows(); rgnInit(); if (!RG.src) { alert("This dataset has no region fields."); return; }
  $("rgnsrc").innerHTML = rgnKeys().map((k) => `<option value="${k}">${k === "path_region" ? "Pathologist outlines (cells inside the outline)" : ROLE_LABEL[k] + " (cell labels)"}</option>`).join(""); $("rgnsrc").value = RG.src;
  const subs = sortedSubtypes(); $("rgnsub").innerHTML = `<option value="__all">All donors (${nDonors(S.meta.samples)})</option>` + subs.map((t) => `<option value="${t}">${t} (${nDonors(S.meta.samples.filter((s) => subtypeOf(s) === t))})</option>`).join("");
  $("rgnsub").value = S.rgnSub && (S.rgnSub === "__all" || subs.includes(S.rgnSub)) ? S.rgnSub : "__all";
  rgnFillSel(); syncTransformUI(); $("rgn").hidden = false; renderRgn();
}
function closeRgn() { $("rgn").hidden = true; }
function rgnChanged() { syncTransformUI(); S.within === "region" ? refresh() : renderRgn(); }
const ringsOf = (f) => (f.geometry.type === "MultiPolygon" ? f.geometry.coordinates.flat() : f.geometry.coordinates);
let rgnToken = 0;
async function renderRgn() {
  if ($("rgn").hidden || RG.level < 0) return;
  const tok = ++rgnToken, key = RG.src, lvl = RG.level, multi = isMulti(key), name = String(display(key, lvl)).split("›").pop().trim(), grid = $("rgngrid");
  const sub = $("rgnsub").value, box = +$("rgnsz").value, wantOut = $("rgnout").checked, ctxOn = $("rgnctx").checked; S.rgnSub = sub;
  const sel = S.rgnSel || (S.rgnSel = new Set()), only = $("rgnonly").checked;
  const groups = sortedSubtypes().filter((t) => sub === "__all" || t === sub).map((t) => ({ t, samples: S.meta.samples.filter((s) => subtypeOf(s) === t && (!only || sel.has(s.id))) })).filter((g) => g.samples.length), samples = groups.flatMap((g) => g.samples);
  $("rgnn").textContent = sel.size; $("rgnsz").disabled = only; $("rgnzoomctl").hidden = !only; grid.classList.toggle("cmpgrid", only); if (!only) grid.removeAttribute("style");
  if (only && !samples.length) { grid.innerHTML = `<div class="muted">No ticked tissue${sub === "__all" ? "" : " in this group"}. Untick “only ticked” and tick the tissues to compare.</div>`; $("rgnleg").innerHTML = ""; $("rgninfo").textContent = ""; return; }
  const what = `Region “${name}” (${key === "path_region" ? "cells inside the pathologist outline" : ROLE_LABEL[key] + " label"})`;
  $("rgninfo").textContent = `${what} · ${donTxt(nDonors(samples))} · colour: ${colorName()}${winNote()} · loading…`; $("rgnleg").innerHTML = "";
  grid.innerHTML = groups.map((g) => (only ? "" : `<div class="galsec"><div class="sechead">${g.t} · ${donTxt(nDonors(g.samples))}</div><div class="secgrid">`) + g.samples.map((s) => `<div class="tile" data-id="${s.id}" title="Open ${dlabel(s.id)} in the main viewer"><input type="checkbox" class="gsel" title="tick to compare this tissue" ${sel.has(s.id) ? "checked" : ""}><span class="zb"></span><canvas width="${box}" height="60"></canvas><div class="cap"><b>${dlabel(s.id)}</b><br><span class="muted">loading…</span></div></div>`).join("") + (only ? "" : `</div></div>`)).join("");
  grid.querySelectorAll(".tile").forEach((t) => {
    t.classList.toggle("sel", sel.has(t.dataset.id));
    t.addEventListener("click", async (e) => { if (only && e.target.closest("canvas")) return; $("sample").value = t.dataset.id; closeRgn(); await loadSample(t.dataset.id); });   // compare view: the picture is for zooming; click the caption to open the tissue
    const cb = t.querySelector(".gsel"); cb.addEventListener("click", (e) => e.stopPropagation());
    cb.addEventListener("change", () => { cb.checked ? sel.add(t.dataset.id) : sel.delete(t.dataset.id); t.classList.toggle("sel", cb.checked); $("rgnn").textContent = sel.size; if (only) renderRgn(); });
  });
  let rows;
  try {
    rows = await Promise.all(samples.map(async (s) => { const [b, codes, feats, poly] = await Promise.all([winBundle(s.id), galGet(s.id, key), galPaths(s.id), only || HEW.rgn.outl ? galPoly(s.id).catch(() => null) : null]); const inR = new Uint8Array(b.n), cells = [];
      for (let i = 0; i < b.n; i++) if (multi ? (codes[i] >> lvl) & 1 : codes[i] === lvl) { inR[i] = 1; cells.push(i); } return { s, b, inR, cells, feats, poly }; }));
  } catch (e) { $("rgninfo").textContent = "Error: " + e.message; return; }
  const heo = heOpt("rgn");
  if (tok !== rgnToken) return;
  const sc = S.meta.cats[S.colorBy] ? makeScale([]) : makeScale(rows.map((r) => ({ ...r.b, pass: r.b.pass.map((p, i) => p & r.inR[i]) })));   // scale pooled over the region's cells only
  $("rgnleg").innerHTML = (wantOut ? `<div class="opt"><span class="sw" style="background:${regionColor(name, null)}"></span><span class="l">${name} boundary</span></div>` : "") + (sc.info.type === "cat" ? "" : scaleLegendHTML(sc.info));
  const shown = S.ctxs.rgn = { sc, shown: [] }; let withOut = 0, drawn = 0; const empty = [], drawnIds = [], outIds = [];
  const dims = only ? cmpLayout("rgn", "rgngrid", Math.max(1, rows.filter((r) => r.cells.length).length), 104) : null;   // compare view: one equal chunk per tissue that has cells in the region
  for (const r of rows) {
    const el = grid.querySelector(`.tile[data-id="${CSS.escape(r.s.id)}"]`); if (!el) continue;
    const oname = name === "Non-sclerotic Gloms" ? "All gloms" : name,   // the non-sclerotic gloms are the pathologist's "All gloms" outlines (sclerotic ones are annotated separately)
      fs = (r.feats || []).filter((f) => f.properties.annotation_name === oname && !f.properties.negative_roa), rings = wantOut ? fs.flatMap(ringsOf) : [];
    if (!r.cells.length) { empty.push(r.s.id); el.style.display = "none"; continue; }
    drawn++; drawnIds.push(r.s.id); if (fs.length) { withOut++; outIds.push(r.s.id); }
    const pen = regionColor(name, fs[0]?.properties.line_color_hex), cv = el.querySelector("canvas"), bb = regionBBox(r, rings), cache = {};
    const he = heo.on && HE.data?.[r.s.id] ? { a: heo.a, sample: r.s.id, redo: () => heRedo(el) } : null;
    const cv2 = only ? hePane(el, cv, "rgn", r.s.id) : null, dm = cv2 ? hePairDims(dims, "rgn") : dims, he2 = cv2 ? { a: 1, sample: r.s.id, redo: () => heRedo(el) } : null;
    const rec = el._dt = { he, last: null, last2: null, bb, draw: (vw, vh) => { rec.last = vw; rec.last2 = vh; drawRegionTile(cv, r, sc, rings, pen, box, ctxOn, bb, vw, dm, cache, r.poly, he, "rgn"); if (cv2) drawRegionTile(cv2, r, sc, rings, pen, box, ctxOn, bb, HEW.rgn.link ? vw : vh, dm, null, null, he2, "rgn"); } }, v0 = only ? rgnCfg.views()[rgnCfg.sync() ? "__all" : r.s.id] || null : null;
    rec.draw(v0, only ? rgnCfg.views()[(rgnCfg.sync() ? "__all" : r.s.id) + "#h"] || null : null); if (only) { attachZoom(el, rec, rgnCfg); if (v0 && v0.z > 1.01) el.querySelector(".zb").textContent = "×" + v0.z.toFixed(1); }
    const area = fs.reduce((a, f) => a + (f.properties.area_um2 || 0), 0) / 1e6, nPass = r.cells.filter((i) => r.b.pass[i]).length;
    el.querySelector(".cap").innerHTML = `<b>${dlabel(r.s.id)}</b><br><span class="muted">${fmt(r.cells.length)} cells in the region (${(100 * r.cells.length / r.b.n).toFixed(0)}% of the donor tissue)${nPass < r.cells.length ? ` · ${fmt(nPass)} pass the filters` : ""}${area ? ` · outline area ${area.toFixed(2)} mm²` : ""}<br>${r.s.disease_type} · ${subtypeOf(r.s)}</span>`;
    shown.shown.push({ el, b: r.b, cells: r.cells });
  }
  grid.querySelectorAll(".galsec").forEach((g) => { if (![...g.querySelectorAll(".tile")].some((t) => t.style.display !== "none")) g.style.display = "none"; });
  $("rgninfo").title = empty.length ? `No cells of this region in: ${empty.join(", ")}` : "";
  $("rgninfo").textContent = `${what} · ${nDonors(drawnIds)} of ${nDonors(samples)} donors have cells in it${empty.length ? ` (${nDonors(empty)} have none)` : ""} · colour: ${colorName()}${winNote()} · ` +
    (!wantOut ? "boundary hidden" : key === "path_region" || withOut ? `boundary drawn for ${donTxt(nDonors(outIds))}` : `no pathologist outline is named “${name}”, so no boundary is drawn`) + " · only cells of the region are shown · " + (only ? "window split per tissue · wheel = zoom (cell outlines appear when zoomed in), drag = pan, double-click = reset, click a name to open the tissue" : "click a tissue to open it · tick tissues to compare them");
}
function regionBBox(r, rings) {   // bounding box of the region's cells and outline (with a small margin)
  const { b, cells } = r; let x0 = 1e18, x1 = -1e18, y0 = 1e18, y1 = -1e18;
  const grow = (x, y) => { if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; };
  for (const i of cells) grow(b.xy[i * 2], b.xy[i * 2 + 1]);
  for (const ring of rings) for (const [x, y] of ring) grow(x, y);
  const m = 0.02 * Math.max(x1 - x0, y1 - y0, 1); x0 -= m; x1 += m; y0 -= m; y1 += m;
  return { x0, y0, w: x1 - x0, h: y1 - y0 };
}
// view = { z, rx, ry } (null = whole region); dims = { W, H } fills a fixed rectangle (compare view).
function drawRegionTile(cv, r, sc, rings, pen, box, ctxOn, bb = regionBBox(r, rings), view = null, dims = null, cache = null, poly = null, he = null, w8 = "rgn") {
  const { b, cells, inR } = r, { x0, y0, w, h } = bb, k = dims ? Math.min((dims.W - 8) / w, (dims.H - 8) / h) : (box - 8) / Math.max(w, h);
  cv.width = dims ? dims.W : Math.max(40, Math.ceil(w * k) + 8); cv.height = dims ? dims.H : Math.max(40, Math.ceil(h * k) + 8);
  const z = view ? view.z : 1, W = cv.width, H = cv.height, kz = k * z, cx = x0 + (view ? view.rx : 0.5) * w, cy = y0 + (view ? view.ry : 0.5) * h;
  const ctx = cv.getContext("2d"), X = (x) => W / 2 + (x - cx) * kz, Y = (y) => H / 2 + (y - cy) * kz, sp = Math.sqrt((w * k * h * k) / Math.max(1, cells.length)), px = Math.max(1.2, Math.min(4, sp * 0.8)) * Math.min(z, 5);
  ctx.fillStyle = "#0f1218"; ctx.fillRect(0, 0, W, H);
  cv._xf = { s: kz, ox: W / 2 - cx * kz, oy: H / 2 - cy * kz };
  if (he) drawHE(ctx, { ...(he.w ? he : { ...he, w: w8 }), side: !!cv._heOnly }, X, Y, W, H, kz);
  const ho = cv._heOnly || !hw(w8).cells;   // H&E beside the cells, or Xenium switched off: no cells, but the region outline and scale bar stay
  ctx.globalAlpha = he ? he.a : 1;
  const inView = (x, y) => { const sx = X(x), sy = Y(y); return sx >= -px && sy >= -px && sx <= W && sy <= H; };
  if (ctxOn && !ho) { ctx.fillStyle = "rgba(105,110,122,0.16)"; for (let i = 0; i < b.n; i++) if (!inR[i]) { const x = b.xy[i * 2], y = b.xy[i * 2 + 1]; if (x >= x0 && x <= x0 + w && y >= y0 && y <= y0 + h && inView(x, y)) ctx.fillRect(X(x), Y(y), px, px); } }
  let order = cache?.order;
  if (!order && !ho) {
    const groups = new Map();
    for (const i of cells) { const col = b.pass[i] ? (S.grey ? rgbStr(GREY) : rgbStr(sc.color(b, i))) : "rgba(105,110,122,0.45)"; (groups.get(col) || groups.set(col, []).get(col)).push(i); }
    order = [...groups.entries()].sort((p, q) => (p[0].startsWith("rgba") ? -1 : 0) - (q[0].startsWith("rgba") ? -1 : 0)); if (cache) cache.order = order;
  }
  if (!ho) drawCells(ctx, order, b, X, Y, W, H, px, poly, he ? he.a : 1, hw(w8).outl);
  ctx.globalAlpha = 1; ctx.lineJoin = "round";
  for (const [lw, col] of [[4, "rgba(0,0,0,0.85)"], [2, pen]]) { ctx.lineWidth = lw; ctx.strokeStyle = col; for (const ring of rings) { ctx.beginPath(); ring.forEach(([x, y], j) => (j ? ctx.lineTo(X(x), Y(y)) : ctx.moveTo(X(x), Y(y)))); ctx.closePath(); ctx.stroke(); } }
  const target = 0.2 * (W - 8) / kz, p = Math.pow(10, Math.floor(Math.log10(target))), nice = [1, 2, 5, 10].map((f) => f * p).reduce((a, c) => (Math.abs(c - target) < Math.abs(a - target) ? c : a));
  const lab = nice >= 1000 ? `${nice / 1000} mm` : `${nice} µm`; ctx.font = `11px ${FONT}`; ctx.textBaseline = "alphabetic";
  ctx.fillStyle = "rgba(15,18,24,0.75)"; ctx.fillRect(6, H - 32, Math.max(nice * kz, ctx.measureText(lab).width) + 10, 28);   // keeps the scale bar legible over cells
  ctx.fillStyle = "#fff"; ctx.fillRect(10, H - 12, nice * kz, 3); ctx.fillText(lab, 10, H - 17);
}

function bindNewUI() {
  $("zregion").addEventListener("change", (e) => { S.within = e.target.value; refresh(); });
  // Regions window
  $("rnbbtn").style.display = S.meta.cats.path_region ? "" : "none"; $("rnbbtn").addEventListener("click", openRnb); $("rnbclose").addEventListener("click", () => ($("rnb").hidden = true));
  $("rnbsel").addEventListener("change", () => { S.rnbSel = new Map(); renderRnb(); }); ["rnbscope", "rnbsz", "rnbonly"].forEach((id) => $(id).addEventListener("change", renderRnb));
  $("rnbrz").addEventListener("click", () => resetZoom(rnbCfg)); $("rnbsync").addEventListener("change", () => resetZoom(rnbCfg)); $("rnbclr").addEventListener("click", () => { S.rnbSel = new Map(); renderRnb(); }); $("rnbd").addEventListener("input", () => { clearTimeout(S.rnbT); S.rnbT = setTimeout(renderRnb, 350); });
  $("rgnbtn").addEventListener("click", openRgn); $("rgnclose").addEventListener("click", closeRgn);
  $("rgnsrc").addEventListener("change", (e) => { RG.src = e.target.value; RG.level = rgnDefault(RG.src); rgnFillSel(); rgnChanged(); });
  $("rgnsel").addEventListener("change", (e) => { RG.level = +e.target.value; rgnChanged(); });
  for (const id of ["rgnsub", "rgnsz", "rgnout", "rgnctx"]) $(id).addEventListener("change", renderRgn);
  $("rgnonly").addEventListener("change", () => { const on = $("rgnonly").checked; if (on) { S.rgnPrevSz = $("rgnsz").value; } else if (S.rgnPrevSz) $("rgnsz").value = S.rgnPrevSz; renderRgn(); });
  $("rgnclr").addEventListener("click", () => { (S.rgnSel || (S.rgnSel = new Set())).clear(); renderRgn(); });
  $("rgnrz").addEventListener("click", () => resetZoom(rgnCfg)); $("rgnsync").addEventListener("change", () => resetZoom(rgnCfg));
  let rT; window.addEventListener("resize", () => { clearTimeout(rT); rT = setTimeout(() => { if (!$("rgn").hidden && $("rgnonly").checked) renderRgn(); }, 250); });
  $("rgnfilt").addEventListener("click", () => { if (RG.level < 0) return; S.filters[RG.src] = new Set([RG.level]); closeRgn(); refresh(); });
  // ticking / comparing glomeruli
  $("glomcmpbtn").addEventListener("click", openGlomCmp); $("ghcmp").addEventListener("click", openGlomCmp);
  $("glomselclr").addEventListener("click", glomUntickAll); $("gcclr").addEventListener("click", glomUntickAll);
  $("gcclose").addEventListener("click", closeGlomCmp); $("gcsave").addEventListener("click", saveGlomCmpPNG);
  for (const id of ["gcsize", "gcctx", "gcown"]) $(id).addEventListener("change", renderGlomCmp);
  updateGlomSelUI();
  // Glom heterogeneity: gene / signature controls act on the main window's state
  const searchG = debounce(async (e) => { const l = await jget("/api/genes/search?q=" + encodeURIComponent(e.target.value)); $("genelist").innerHTML = l.map((g) => `<option value="${g}">`).join(""); });
  for (const id of ["ghgene", "ghgene2"]) $(id).addEventListener("input", searchG);
  $("ghcol").addEventListener("change", () => { ghCtlVisibility(); ghPush(); });
  for (const id of ["ghgene", "ghgene2", "ghsig", "ghsig2", "ghsigtxt", "ghsigtxt2"]) $(id).addEventListener("change", ghPush);
  const ghSet = (fn) => { fn(); $("ghsrc").value = "cells"; refresh(); };
  $("ghsm").addEventListener("input", (e) => { $("ghsmlab").textContent = +e.target.value ? e.target.value + " µm" : "off"; });
  $("ghsm").addEventListener("change", (e) => ghSet(() => { S.smooth = +e.target.value; })); $("ghzs").addEventListener("change", (e) => ghSet(() => { S.zscope = e.target.value; }));
  $("ghexl").addEventListener("change", (e) => ghSet(() => { S.layer = e.target.value; })); $("ghwi").addEventListener("change", (e) => ghSet(() => { S.within = e.target.value; }));
  $("ghsrc").addEventListener("change", renderGH);
}

// arrows on the canvases of every side-by-side pair (gallery, regions, glomeruli, compare, neighbourhoods, ROI detail)
drawTile = ptWrap(drawTile, (a) => a[2].sample); drawGlom = ptWrap(drawGlom, (a) => a[1].sample); drawRegionTile = ptWrap(drawRegionTile, (a) => a[1].s.id);

init().catch((e) => { loading(false); document.body.insertAdjacentHTML("afterbegin", `<pre style="color:#b00;padding:12px">Failed to start: ${e.message}</pre>`); });
