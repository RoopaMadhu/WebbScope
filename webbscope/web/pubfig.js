// Publication figure of the main view: Nature-style layout (sans-serif Arial/Helvetica, 5-7 pt text, 8 pt bold panel letter, 89 / 136 / 183 mm widths),
// drawn at a physical size (mm) and resolution (dpi). One layout (list of ops in points) feeds two writers: a PNG with its dpi recorded, and a PDF
// whose text, scale bar and colour bar are vector (Helvetica) while the cells are one raster image at the chosen dpi.
const PUB = { FONT: "Arial, Helvetica, 'Liberation Sans', sans-serif", M: 2, BG: [255, 255, 255], DARK: [15, 18, 24] };
const mm2pt = (v) => (v * 72) / 25.4;
const pubMeasure = document.createElement("canvas").getContext("2d");
const pubTW = (s, size, bold) => { pubMeasure.font = `${bold ? "bold " : ""}${size}px ${PUB.FONT}`; return pubMeasure.measureText(s).width; };
const pubT = (s, x, y, size, bold, c) => ({ t: "t", s, x, y, size, bold, c });
const pubR = (x, y, w, h, c) => ({ t: "r", x, y, w, h, c });


function pubPoly() {   // the viewer keeps polygon y flipped (deck.gl); the figure uses raw tissue coordinates like S.xy
  if (!S.poly) return null;
  if (S._pubPoly?.src !== S.poly) { const xy = new Float32Array(S.poly.xy); for (let k = 1; k < xy.length; k += 2) xy[k] = -xy[k]; S._pubPoly = { src: S.poly, val: { off: S.poly.off, xy } }; }
  return S._pubPoly.val;
}
function pubLegend(o, W, ink, info, tickFor, catItems) {   // catItems: [[name, rgb]] overrides the main tissue's categories (window figures)   // ops relative to y = 0; returns { ops, h }
  const M = PUB.M * o.F, B = o.size, sm = Math.max(5, B - 1), ops = [], mode = S.colorBy; let y = 0;
  const cap = (t) => { ops.push(pubT(t, M, y + B, B, true, ink)); y += B * 1.5; };
  if (info.type === "none") return { ops, h: 0 };
  if (info.type === "cat") {
    cap(colorName());
    let items;
    if (catItems) items = catItems;
    else { const { out, miss } = countsFor(mode); items = out.map((n, i) => [i, n]).filter(([, n]) => n > 0).sort((a, b) => b[1] - a[1]).map(([i]) => [display(mode, i), hex2rgb(catColor(mode, i))]); if (miss) items.push(["missing", GRAY]); }
    const SW = sm * 0.85, row = sm * 1.7, cap_ = 60; let px = M;
    for (const [name, col] of items.slice(0, cap_)) {
      const w = SW + 3 + pubTW(name, sm) + 9; if (px + w > W - M && px > M) { px = M; y += row; }
      ops.push(pubR(px, y + sm - SW + 0.3, SW, SW, col), pubT(name, px + SW + 3, y + sm, sm, false, ink)); px += w;
    }
    if (items.length > cap_) { const t = `+ ${items.length - cap_} more`; if (px + pubTW(t, sm) > W - M) { px = M; y += row; } ops.push(pubT(t, px, y + sm, sm, false, ink)); }
    return { ops, h: y + row };
  }
  if (info.type === "dual" || info.type === "sigdual") {
    const names = info.type === "dual" ? [S.geneA, S.geneB] : [info.labelA, info.labelB || "(no signature B)"], cols = info.type === "dual" ? [[211, 51, 51], [51, 170, 51]] : [hex2rgb(SIG_A), hex2rgb(SIG_B)], SW = sm * 0.85;
    let px = M; names.forEach((n, k) => { ops.push(pubR(px, y + sm - SW + 0.3, SW, SW, cols[k]), pubT(n, px + SW + 3, y + sm, sm, false, ink)); px += SW + 3 + pubTW(n, sm) + 12; });
    y += sm * 1.7; ops.push(pubT(info.type === "dual" ? "yellow = both; grey = neither detected" : "near white = both high; each signature scaled separately", M, y + sm, sm, false, ink));
    return { ops, h: y + sm * 1.7 };
  }
  // continuous bar (gene / signature / z-score / numeric column)
  const lay = S.layer === "counts" ? "raw counts" : S.meta.ln_given ? "log-normalised" : `log1p(CP${fmt(S.meta.norm_scale)})`;
  cap(info.type === "num" ? `${labelOf(mode)}${info.log ? " (log scale)" : ""}` : `${info.label} (${lay})`);
  const L = info.type === "gene" && info.cont && info.div ? DIVL : LUT, g = document.createElement("canvas"); g.width = 256; g.height = 1;
  const gx = g.getContext("2d"); L.forEach((c, i) => { gx.fillStyle = rgbStr(c); gx.fillRect(i, 0, 1, 1); });
  const bw = Math.min(W - 2 * M - 24 * o.F, mm2pt(50 * o.F)), bx = M + 12 * o.F;   // 12 pt of room each side so the end labels are not clipped
  ops.push({ t: "i", c: g, x: bx, y, w: bw, h: B * 0.8, smooth: true }); y += B * 0.8;
  for (const k of scaleTicks(info)) {
    const s = fmtTick(k.v), w = pubTW(s, sm); ops.push(pubT(s, bx + bw * k.t - w / 2, y + sm + 1.5 * o.F, sm, false, ink));
    ops.push(pubR(bx + bw * k.t - 0.125 * o.F, y, 0.25 * o.F, 1.5 * o.F, ink));
  }
  return { ops, h: y + sm * 1.8 };
}

// ---- Window figures: the tissues of the open Gallery, or the ticked glomeruli of the Glomeruli window, as a grid of tiles with a scale bar and a name each.
// An entry = { name, meta, short, wu, hu (extent in um), view, draw(cv, view) -> Promise<codes passing the filters | null> }.
const pubSideOn = (w, id, on = true) => on && hw(w).side && !!(HE.data?.[id] || imcSide(w));   // this tile shows an image pane beside the cells
const pubSideView = (w, cfg, key, v) => (hw(w).link ? v : cfg.views()[(cfg.sync() ? "__all" : key) + "#h"] || null);
function pubWinSpec() {
  if (!$("gallery").hidden && S.ctxs.gallery?.ids.length) {
    const G = S.ctxs.gallery, only = $("galonly").checked;
    return { key: G, sc: G.sc, title: only ? `Gallery — ${G.ids.length} ticked` : `Gallery — ${$("galsub").selectedOptions[0].text}`, entries: G.ids.map((id, k) => {
      const smp = S.meta.samples.find((z) => z.id === id), v = only ? galCfg.views()[galCfg.sync() ? "__all" : id] || null : null;
      return { sample: id, name: dlabel(id), meta: `${fmt(smp.n)} cells${gTag(smp)}`, short: `${fmt(smp.n)} cells`, wu: smp.x[1] - smp.x[0] || 1, hu: smp.y[1] - smp.y[0] || 1, view: v, z: v ? v.z : 1,
        hw: "gal", side: pubSideOn("gal", id, only), viewSide: only ? pubSideView("gal", galCfg, id, v) : null,
        drawSide: async (cv, view) => { const b = await G.bp[k], he = { a: 1, sample: id, redo: () => {}, miss: { n: 0 } }; cv._heOnly = true; for (let t = 0; t < 40; t++) { he.miss.n = 0; drawTile(cv, smp, b, G.sc, 0, view, null, { W: cv.width, H: cv.height }, null, he, "gal", cv._ps || 1); if (!he.miss.n) break; await new Promise((r) => setTimeout(r, 250)); } },
        draw: async (cv, view) => {
          const b = await G.bp[k], heo = heOpt("gal"), poly = only || HEW.gal.outl ? await galPoly(id).catch(() => null) : null;
          const he = heo.on && HE.data?.[id] ? { a: heo.a, sample: id, redo: () => {}, miss: { n: 0 } } : null;
          for (let t = 0; t < 40; t++) { if (he) he.miss.n = 0; drawTile(cv, smp, b, G.sc, 0, view, null, { W: cv.width, H: cv.height }, poly, he, "gal", cv._ps || 1); if (!he || !he.miss.n) break; await new Promise((r) => setTimeout(r, 250)); }
          const out = []; for (let j = 0; j < b.n; j++) if (b.pass[j]) out.push(b.codes[j]); return out;
        } };
    }) };
  }
  // tiles drawn by drawGlom (glomeruli, region neighbourhoods): one entry per glom / region object
  const glomEntries = (items, sc, hwk, ctxFlag, cfg, nameOf, metaOf, shortOf) => items.map((r) => {
    const B = glomBox(r.cells, r.poly, 12, r.outline), key = r.gkey || r.key, v = cfg.views()[cfg.sync() ? "__all" : key] || null;
    const paint = async (cv, view, ctx) => {   // waits until every image tile it needs has arrived (at most ~10 s)
      cv._hw = hwk; cv._noSB = true; cv._heMiss = { n: 0 };
      for (let t = 0; t < 40; t++) { cv._heMiss.n = 0; drawGlom(cv, r.b, sc, r.poly, r.cells, ctx, -1, r.outline, () => {}, view, { W: cv.width, H: cv.height }); if (!cv._heMiss.n) break; await new Promise((q) => setTimeout(q, 250)); }
    };
    return { sample: r.s.id, hw: hwk, name: nameOf(r), meta: metaOf(r), short: shortOf(r), wu: B.x1 - B.x0, hu: B.y1 - B.y0, view: v, z: v ? v.z : 1,
      side: pubSideOn(hwk, r.s.id, hwk !== "rnb" || !!S.ctxs.rnbFig?.only), viewSide: pubSideView(hwk, cfg, key, v),
      draw: async (cv, view) => { await paint(cv, view, ctxFlag); return r.cells.filter((i) => r.b.pass[i]).map((i) => r.b.codes[i]); },
      drawSide: async (cv, view) => { cv._heOnly = true; await paint(cv, view, ctxFlag); } };
  });
  if (!$("glomcmp").hidden && S.gcCards?.length && S.gcSc) {   // Compare glomeruli: every ticked glom, with the same colour scale as the window
    const items = S.gcCards.map((c) => ({ ...c, b: c.d.b, poly: c.d.poly })), cfg = { views: () => ({}), sync: () => false };
    return { key: S.gcCards, sc: S.gcSc, title: `Glomeruli comparison (${items.length})`, entries: glomEntries(items, S.gcSc, "gcmp", $("gcctx").checked, cfg, (r) => r.label,
      (r) => `${dlabel(r.s.id)}${gTag(r.s)} · ${fmt(r.cells.length)} cells · ${r.qc}`, (r) => `${dlabel(r.s.id)} · ${fmt(r.cells.length)} cells`) };
  }
  const GO = S.ctxs.glomOnly;
  if (!$("glom").hidden && $("glomonly").checked && GO?.items.length)
    return { key: GO, sc: GO.sc, title: `Glomeruli — ${GO.items.length} ticked`, entries: glomEntries(GO.items, GO.sc, "glom", true, glomCfg, (r) => r.label,
      (r) => `${dlabel(r.s.id)}${gTag(r.s)} · ${fmt(r.cells.length)} cells · ${r.qc}`, (r) => `${dlabel(r.s.id)} · ${fmt(r.cells.length)} cells`) };
  const RF = S.ctxs.rnbFig;
  if (!$("rnb").hidden && RF) {   // ticked regions; with none ticked, every region listed (at most 40)
    const items = RF.only ? RF.items : (() => { const t = RF.all.filter((x) => RF.sel.has(x.key)); return t.length ? t : RF.all.slice(0, 40); })();
    if (items.length) return { key: RF, sc: RF.sc, title: `Region neighbourhoods — ${RF.name} (${items.length} ${items.length === 1 ? "region" : "regions"}, cells inside + within ${RF.d} µm)`,
      entries: glomEntries(items, RF.sc, "rnb", false, rnbCfg, (r) => r.label, (r) => `${dlabel(r.s.id)}${gTag(r.s)} · ${fmt(r.inside.length)} inside · ${fmt(r.ring.length)} within ${RF.d} µm`, (r) => `${dlabel(r.s.id)} · ${fmt(r.inside.length)} cells inside`) };
  }
  const RG = S.ctxs.rgnFig;
  if (!$("rgn").hidden && RG?.tiles.length) {
    return { key: RG, sc: RG.sc, title: `Region — ${RG.name}`, entries: RG.tiles.map((T) => {
      const r = T.r, id = r.s.id, v = RG.only ? rgnCfg.views()[rgnCfg.sync() ? "__all" : id] || null : null;
      return { sample: id, hw: "rgn", name: dlabel(id), meta: `${fmt(r.cells.length)} cells in region${gTag(r.s)}`, short: `${fmt(r.cells.length)} cells`, wu: T.bb.w, hu: T.bb.h, view: v, z: v ? v.z : 1,
        side: pubSideOn("rgn", id, RG.only), viewSide: RG.only ? pubSideView("rgn", rgnCfg, id, v) : null,
        drawSide: async (cv, view) => { const he = { a: 1, sample: id, redo: () => {}, miss: { n: 0 } }; cv._heOnly = true; cv._noSB = true; for (let t = 0; t < 40; t++) { he.miss.n = 0; drawRegionTile(cv, r, RG.sc, T.rings, T.pen, 0, T.ctxOn, T.bb, view, { W: cv.width, H: cv.height }, null, null, he, "rgn"); if (!he.miss.n) break; await new Promise((q) => setTimeout(q, 250)); } },
        draw: async (cv, view) => {
          const heo = heOpt("rgn"), he = heo.on && HE.data?.[id] ? { a: heo.a, sample: id, redo: () => {}, miss: { n: 0 } } : null; cv._noSB = true;
          for (let t = 0; t < 40; t++) { if (he) he.miss.n = 0; drawRegionTile(cv, r, RG.sc, T.rings, T.pen, 0, T.ctxOn, T.bb, view, { W: cv.width, H: cv.height }, null, r.poly, he, "rgn"); if (!he || !he.miss.n) break; await new Promise((q) => setTimeout(q, 250)); }
          return r.cells.filter((i) => r.b.pass[i]).map((i) => r.b.codes[i]);
        } };
    }) };
  }
  if (!$("slide").hidden) {
    const A = $("slidesel").value === "__all" ? S.slideAll : null, one = S.ctxs.slideFig, cur = A && A.tok === slideTok ? A : null;
    const list = cur ? cur.bySlide.map((q, i) => ({ g: q.g, samples: q.samples, geom: cur.geoms[i] })) : !A && one && one.g === $("slidesel").value ? [{ g: one.g, samples: one.samples, geom: slideGeom(one.samples) }] : [];
    if (list.length) { const sc = cur ? cur.sc : one.sc, getBs = (q) => (cur ? Promise.all(q.samples.map((z) => cur.getB(z.id))) : Promise.resolve(one.bs));
      const fit = (q, W, H, p) => Math.min(W / p / q.geom.w, (H / p - SLIDE_TOP - SLIDE_BOT) / q.geom.h);   // CSS px per um of drawSlideInto
      const paint = async (cv, q, imgOnly) => {
        const p = cv._ps || 1, k = fit(q, cv.width, cv.height, p), bs = await getBs(q), tmp = document.createElement("canvas"), miss = { n: 0 }, opt = { imgs: imgOnly || HEW.slide.on, imgOnly, redo: () => {}, miss, fig: true, bg: cv._bg, light: cv._bg === "#ffffff" };
        for (let t = 0; t < 40; t++) { miss.n = 0; drawSlideInto(tmp, q.g, q.samples, bs, sc, k, q.geom, p, 1, opt); if (!miss.n) break; await new Promise((r) => setTimeout(r, 250)); }
        const c = cv.getContext("2d"); c.fillStyle = cv._bg; c.fillRect(0, 0, cv.width, cv.height); c.drawImage(tmp, Math.round((cv.width - tmp.width) / 2), Math.round((cv.height - tmp.height) / 2)); cv._xf = { s: k * p, ox: 0, oy: 0 };
      };
      return { key: cur || one, sc, title: `Slide — ${list.length > 1 ? `all ${list.length} slides` : list[0].g}`, entries: list.map((q) => {
        const n = q.samples.reduce((a, z) => a + z.n, 0);
        return { name: q.g, meta: `${q.samples.length} tissue${q.samples.length === 1 ? "" : "s"} · ${fmt(n)} cells`, short: `${fmt(n)} cells`, wu: q.geom.w, hu: q.geom.h, view: null, z: 1, hw: "slide", side: !cur && !!HEW.slide.side, viewSide: null, scale: (W, H, p) => fit(q, W, H, p) * p,
          draw: async (cv) => { await paint(cv, q, false); const bs = await getBs(q), out = []; for (const b of bs) for (let j = 0; j < b.n; j++) if (b.pass[j]) out.push(b.codes[j]); return out; },
          drawSide: async (cv) => { await paint(cv, q, true); } };
      }) }; }
  }
  return null;
}
// Panel letters of a multi-panel figure: the first letter is the user's choice, the rest follow (a, b, c ... z, aa, ab ...; A, B, C ...; 1, 2, 3 ...).
function pubLetter(first, i) {
  const f = String(first).trim(); if (!f) return "";
  if (/^\d+$/.test(f)) return String(+f + i);
  const up = /^[A-Z]/.test(f), al = up ? "ABCDEFGHIJKLMNOPQRSTUVWXYZ" : "abcdefghijklmnopqrstuvwxyz", k = al.indexOf(f[0]);
  if (f.length === 1 && k >= 0) { const n = k + i; return n < 26 ? al[n] : al[Math.floor(n / 26) - 1] + al[n % 26]; }
  return f + (i + 1);   // anything else: the text followed by a number
}
function pubBuildWin(o, dpi, spec) {
  const M = PUB.M * o.F, W = mm2pt(o.wmm), dark = o.dark, ink = dark ? [255, 255, 255] : [0, 0, 0], bg = dark ? PUB.DARK : PUB.BG, B = o.size, sm = Math.max(5, B - 1), ops = [], s = dpi / 72, pd = o.panel === "dark";
  PUB.panelBg = pd ? "#0f1218" : "#ffffff"; const panelInk = pd ? [255, 255, 255] : [0, 0, 0];
  const E = spec.entries, n = E.length, info = spec.sc.info, pair = E.some((e) => e.side);   // pair: each tile is the cells and the image beside them
  const LH = o.letter ? (B + 1) + 2 * o.F : 0, gap = mm2pt(1.5), g2 = pair ? mm2pt(0.8) : 0, cols = pair ? (n <= 2 ? n : Math.min(4, Math.ceil(Math.sqrt(n)))) : n <= 3 ? n : Math.min(8, Math.ceil(Math.sqrt(n))), rows = Math.ceil(n / cols), tw = (W - gap * (cols - 1)) / cols, pw = pair ? (tw - g2) / 2 : tw, capH = sm * 2.9 + LH;
  const asp = Math.max(0.5, Math.min(1.6, Math.max(...E.map((e) => e.hu / e.wu)))), hdr = o.title ? (B + 1) + 4 * o.F : 0;
  const catItems = info.type === "cat" ? (S.pubLegKey === spec.key ? S.pubLegItems : []) : null;   // categories seen in the tiles: known after the first render
  const leg = pubLegend(o, W, ink, info, null, catItems), availH = Math.max(40, mm2pt(o.maxh) - hdr - 4 - leg.h - 2), th = Math.max(20, Math.min(pw * asp, (availH - rows * capH - (rows - 1) * gap) / rows));
  if (o.title) ops.push(pubT(`${spec.title} · colour: ${colorName()}`, M, B + 1, B, false, ink));
  const cwp = Math.max(8, Math.round(pw * s)), chp = Math.max(8, Math.round(th * s)), top = hdr, cvs = [], cvs2 = [], ps = Math.max(1, cwp / 800);
  const newCv = () => { const cv = document.createElement("canvas"); cv.width = cwp; cv.height = chp; const c = cv.getContext("2d"); c.fillStyle = PUB.panelBg; c.fillRect(0, 0, cwp, chp); cv._bg = PUB.panelBg; cv._ps = ps; cv._ptFree = true; return cv; };
  const ifOnly = (e, side) => { if (!e.sample || !HE.data?.[e.sample] || (!side && !hw(e.hw).on)) return false; const l = heSel(e.sample, e.hw, !!side); return l.length > 0 && l.every((q) => q.key.startsWith("if:")); };   // IF is drawn on black: its scale bar is white
  const bar = (e, z, x, yi, side) => {
    const kpt = (e.scale ? e.scale(cwp, chp, ps) : Math.min((cwp - 8) / e.wu, (chp - 8) / e.hu) * z) / s;   // points per um in this tile (same formula as the tile renderers)
    const target = 0.2 * pw / kpt, p = 10 ** Math.floor(Math.log10(target)), nice = [1, 2, 5, 10].map((q) => q * p).reduce((a, c2) => (Math.abs(c2 - target) < Math.abs(a - target) ? c2 : a));
    const lab = nice >= 1000 ? nice / 1000 + " mm" : nice + " µm", WH = ifOnly(e, side) ? [255, 255, 255] : panelInk;
    ops.push(pubR(x + 4 * o.F, yi + th - 5 * o.F, nice * kpt, 1.6 * o.F, WH), pubT(lab, x + 4 * o.F, yi + th - 7 * o.F, sm, false, WH));
  };
  E.forEach((e, i) => {
    const x = (i % cols) * (tw + gap), y = top + Math.floor(i / cols) * (th + capH + gap), cv = newCv(); cvs.push(cv);
    if (LH) ops.push(pubT(pubLetter(o.letter, i), x + 1, y + B + 1, B + 1, true, ink));   // panel letter: the chosen first letter, then the next ones
    const yi = y + LH;
    ops.push({ t: "i", c: cv, x, y: yi, w: pw, h: th, smooth: false }); bar(e, e.z, x, yi, false);
    if (e.side) {
      const c2 = newCv(), x2 = x + pw + g2; cvs2[i] = c2; ops.push({ t: "i", c: c2, x: x2, y: yi, w: pw, h: th, smooth: true });
      if (e.hw !== "slide" && imcSide(e.hw)) ops.push(pubT("IMC not aligned", x2 + 4 * o.F, yi + th - 4 * o.F, sm, false, pd ? [255, 210, 122] : [130, 85, 0]));   // not in the cells' frame: no scale bar from them
      else bar(e, e.viewSide ? e.viewSide.z : 1, x2, yi, true);
    }
    const nw = pubTW(e.name, sm, true); let nm = e.name;
    while (nm.length > 4 && pubTW(nm, sm, true) > tw - 2) nm = nm.slice(0, -2);
    if (nm !== e.name) nm += "…";
    ops.push(pubT(nm, x + 1, yi + th + sm + 1, sm, true, ink));
    if (nm === e.name && nw + 4 + pubTW(e.meta, sm) <= tw) ops.push(pubT(e.meta, x + 1 + nw + 4, yi + th + sm + 1, sm, false, ink));
    else { let t = e.short; while (t.length > 4 && pubTW(t, sm) > tw - 2) t = t.slice(0, -2) + "…"; ops.push(pubT(t, x + 1, yi + th + sm * 2 + 2, sm, false, ink)); }
  });
  const ly = top + rows * (th + capH + gap) - gap + 4; for (const op of leg.ops) { op.y += ly; ops.push(op); }
  const H = ly + leg.h + 2, plan = { ops, W, H, bg, dpi, px: [Math.round(W * s), Math.round(H * s)], what: spec.title, file: spec.title.split(/[ —]/)[0].toLowerCase() };
  plan.fill = async () => {
    const seen = new Map();
    for (let i = 0; i < E.length; i++) { cvs[i]._ptSample = E[i].sample; const codes = await E[i].draw(cvs[i], E[i].view); if (codes) for (const cd of codes) seen.set(cd, (seen.get(cd) || 0) + 1); if (cvs2[i]) { cvs2[i]._ptSample = E[i].sample; cvs2[i]._hw = E[i].hw; await E[i].drawSide(cvs2[i], E[i].viewSide); } }
    if (info.type === "cat" && S.pubLegKey !== spec.key) {   // first render: legend of the categories actually drawn, then lay out again around it
      S.pubLegItems = [...seen.entries()].sort((a, c) => c[1] - a[1]).map(([cd]) => [cd === MISSING ? "missing" : display(S.colorBy, cd), cd === MISSING ? GRAY : hex2rgb(catColor(S.colorBy, cd))]); S.pubLegKey = spec.key;
      const again = pubBuildWin(o, dpi, spec); Object.assign(plan, again); await again.fill();
    }
  };
  return plan;
}

// ---- Figures of the plot windows (Expression plots, Signatures, Neighbours, Neighbourhoods): the plots as they are on screen, redrawn as layout ops,
// so the PDF keeps text, lines, bars and dots as vectors. No screenshot is involved. A light panel swaps the dark theme's neutral inks; data colours stay.
const pubCx = document.createElement("canvas").getContext("2d");
function pubColor(v) {   // { c: [r, g, b], a } of a CSS colour, or null when it paints nothing
  if (!v || v === "none" || v === "transparent") return null;
  pubCx.fillStyle = "#010203"; pubCx.fillStyle = v; const t = pubCx.fillStyle; if (t === "#010203" && !/010203/.test(v)) return null;
  if (t[0] === "#") return { c: [1, 3, 5].map((i) => parseInt(t.slice(i, i + 2), 16)), a: 1 };
  const m = t.match(/[\d.]+/g).map(Number); return m[3] === 0 ? null : { c: m.slice(0, 3), a: m[3] ?? 1 };
}
function pubTheme(light, pbg) {
  const neutral = (c) => Math.max(...c) - Math.min(...c) < 30, lum = (c) => 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2], blend = (p) => p.c.map((v, i) => Math.round(v * p.a + pbg[i] * (1 - p.a)));
  const surface = (c) => c[0] === 15 && c[1] === 18 && c[2] === 24;    // the windows' background colour
  return {
    text: (p, onColour) => { const c = blend(p); return !light || onColour || !neutral(c) ? c : lum(c) > 165 ? [0, 0, 0] : [70, 74, 80]; },
    fill: (p) => { if (light && p.a === 1 && surface(p.c)) return pbg; return blend(p); },
    line: (p, w) => { if (light && p.a === 1 && surface(p.c)) return pbg; const c = blend(p); if (!light || !neutral(c)) return c; if (lum(c) > 150) return [34, 36, 40]; return w <= 1.2 ? [217, 217, 217] : c; },
    box: (p) => { const c = blend(p); return light && neutral(c) && lum(c) < 75 ? null : c; },    // dark panel fills drop out on a light page
  };
}
function pubCapture(root, th) {   // everything drawn inside `root`, in CSS px from its top left
  const R = root.getBoundingClientRect(), ops = []; let maxR = 0, maxB = 0;
  const grow = (x, y) => { if (x > maxR) maxR = x; if (y > maxB) maxB = y; };
  const svg = (el, onCol) => {
    const r = el.getBoundingClientRect(), vb = el.viewBox?.baseVal, sx = vb && vb.width ? r.width / vb.width : 1, sy = vb && vb.height ? r.height / vb.height : 1, ox = r.left - R.left - (vb?.x || 0) * sx, oy = r.top - R.top - (vb?.y || 0) * sy;
    const X = (v) => ox + v * sx, Y = (v) => oy + v * sy, num = (e, a, d = 0) => { const v = parseFloat(e.getAttribute(a)); return Number.isFinite(v) ? v : d; };
    const paint = (e, a) => { const p = pubColor(e.getAttribute(a) ?? (a === "fill" ? "#000" : null)); if (!p) return null; const q = e.hasAttribute(a + "-opacity") ? +e.getAttribute(a + "-opacity") : 1; return q > 0 ? { c: p.c, a: p.a * q } : null; };
    const shape = (e, pts, close) => { const f = close ? paint(e, "fill") : null, st = paint(e, "stroke"), lw = st ? num(e, "stroke-width", 1) : 0; if (!f && !st) return; ops.push({ t: "p", pts, close, fill: f ? th.fill(f) : null, stroke: st ? th.line(st, lw) : null, lw: lw * sx }); for (const q of pts) grow(q[0], q[1]); };
    const walk = (e) => {
      const tag = e.tagName.toLowerCase();
      if (tag === "g" || tag === "svg") { for (const c of e.children) walk(c); return; }
      if (tag === "line") shape(e, [[X(num(e, "x1")), Y(num(e, "y1"))], [X(num(e, "x2")), Y(num(e, "y2"))]], false);
      else if (tag === "rect") { const x = X(num(e, "x")), y = Y(num(e, "y")), w = num(e, "width") * sx, h = num(e, "height") * sy; if (w > 0 && h > 0) shape(e, [[x, y], [x + w, y], [x + w, y + h], [x, y + h]], true); }
      else if (tag === "circle") { const cx = X(num(e, "cx")), cy = Y(num(e, "cy")), rr = num(e, "r") * sx, pts = []; if (rr <= 0) return; for (let i = 0; i < 28; i++) pts.push([cx + rr * Math.cos(i * Math.PI / 14), cy + rr * Math.sin(i * Math.PI / 14)]); shape(e, pts, true); }
      else if (tag === "path" || tag === "polygon") {
        const d = tag === "path" ? e.getAttribute("d") || "" : "M" + (e.getAttribute("points") || "");
        if (/[^MLZmlz\d\s,.eE+-]/.test(d)) { const n = 80, L = e.getTotalLength(), pts = []; for (let i = 0; i <= n; i++) { const q = e.getPointAtLength((i / n) * L); pts.push([X(q.x), Y(q.y)]); } shape(e, pts, /z\s*$/i.test(d)); return; }
        for (const sub of d.split(/(?=M)/i)) { const v = sub.match(/-?[\d.]+(?:e[-+]?\d+)?/gi); if (!v || v.length < 4) continue; const pts = []; for (let i = 0; i + 1 < v.length; i += 2) pts.push([X(+v[i]), Y(+v[i + 1])]); shape(e, pts, tag === "polygon" || /z\s*$/i.test(sub)); }
      } else if (tag === "text") {
        const str = (e.textContent || "").trim(); if (!str) return; const col = paint(e, "fill"); if (!col) return;
        const size = num(e, "font-size", 11) * sx, bold = +(e.getAttribute("font-weight") || 400) >= 600, w = pubTW(str, size, bold), an = e.getAttribute("text-anchor"), off = an === "middle" ? -w / 2 : an === "end" ? -w : 0;
        let x = num(e, "x"), y = num(e, "y"), rot = 0; const m = /rotate\(\s*(-?[\d.]+)(?:[ ,]+(-?[\d.]+)[ ,]+(-?[\d.]+))?\s*\)/.exec(e.getAttribute("transform") || "");
        if (m) { rot = +m[1]; const px = m[2] != null ? +m[2] : 0, py = m[3] != null ? +m[3] : 0, t = rot * Math.PI / 180, dx = x - px, dy = y - py; x = px + dx * Math.cos(t) - dy * Math.sin(t); y = py + dx * Math.sin(t) + dy * Math.cos(t); }
        const t = rot * Math.PI / 180, tx = X(x) + off * Math.cos(t), ty = Y(y) + off * Math.sin(t); ops.push({ t: "t", s: str, x: tx, y: ty, size, bold, c: th.text(col, onCol), rot }); grow(tx + w * Math.cos(t), ty + w * Math.sin(t)); grow(tx, ty);
      }
    };
    for (const c of el.children) walk(c);
  };
  const textNode = (n, cs, onCol) => {
    const raw = n.nodeValue; if (!raw || !raw.trim()) return; const p = pubColor(cs.color); if (!p) return;
    const rg = document.createRange(), lines = [], re = /\S+/g; let m;
    while ((m = re.exec(raw))) { rg.setStart(n, m.index); rg.setEnd(n, m.index + m[0].length); const rs = rg.getClientRects(); if (!rs.length) continue; const a = rs[0], L = lines[lines.length - 1]; if (L && Math.abs(L.top - a.top) < 3) L.s += " " + m[0]; else lines.push({ s: m[0], left: a.left, top: a.top, h: a.height }); }
    const size = parseFloat(cs.fontSize) || 12, bold = +cs.fontWeight >= 600, par = n.parentElement, clip = cs.textOverflow === "ellipsis" ? par.getBoundingClientRect().right - R.left : Infinity;
    for (const L of lines) {
      let t = cs.textTransform === "uppercase" ? L.s.toUpperCase() : L.s; const x = L.left - R.left, room = clip - x;
      if (pubTW(t, size, bold) > room) { while (t.length > 1 && pubTW(t + "…", size, bold) > room) t = t.slice(0, -1); t += "…"; }
      const y = L.top - R.top + L.h / 2 + size * 0.35; ops.push({ t: "t", s: t, x, y, size, bold, c: th.text(p, onCol), rot: 0 }); grow(x + pubTW(t, size, bold), y + size * 0.3);
    }
  };
  const walk = (el, onCol) => {
    if (el.nodeType === 3) return textNode(el, getComputedStyle(el.parentElement), onCol);
    if (el.nodeType !== 1 || ["INPUT", "SELECT", "BUTTON", "TEXTAREA", "SCRIPT", "STYLE", "DATALIST", "TITLE"].includes(el.tagName)) return;
    if (el.matches(".sctl, .sdesc, .gsel, .zb") || (el.matches(".snote") && !el.parentElement.querySelector(":scope > .legbar"))) return;
    const cs = getComputedStyle(el); if (cs.display === "none" || cs.visibility === "hidden") return;
    const r = el.getBoundingClientRect(), x = r.left - R.left, y = r.top - R.top;
    if (el.tagName.toLowerCase() === "svg") { svg(el, onCol); grow(x + r.width, y + r.height); return; }
    if (el.tagName === "CANVAS") { if (r.width && r.height) { ops.push({ t: "i", c: el, x, y, w: r.width, h: r.height, smooth: true, fig: el._figDraw || null }); grow(x + r.width, y + r.height); } return; }
    if (r.width > 0 && r.height > 0) {
      const b = pubColor(cs.backgroundColor), box = b ? th.box(b) : null; if (box) { ops.push(pubR(x, y, r.width, r.height, box)); onCol = true; grow(x + r.width, y + r.height); }
      if (cs.backgroundImage.includes("gradient")) { const st = cs.backgroundImage.match(/rgba?\([^)]+\)/g) || [], g = document.createElement("canvas"); g.width = 256; g.height = 1; const gx = g.getContext("2d"), lg = gx.createLinearGradient(0, 0, 256, 0); st.forEach((q, i) => lg.addColorStop(st.length > 1 ? i / (st.length - 1) : 0, q)); gx.fillStyle = lg; gx.fillRect(0, 0, 256, 1); ops.push({ t: "i", c: g, x, y, w: r.width, h: r.height, smooth: true }); grow(x + r.width, y + r.height); }
      if (el.tagName === "TH" && parseFloat(cs.borderBottomWidth) > 0) ops.push(pubR(x, y + r.height - 0.6, r.width, 0.6, th.line({ c: [41, 48, 59], a: 1 }, 1)));
    }
    for (const c of el.childNodes) walk(c, onCol);
  };
  if (root.tagName.toLowerCase() === "svg" || root.tagName === "CANVAS") walk(root, false); else for (const c of root.childNodes) walk(c, false);   // a block can be a plot itself
  return { ops, w: maxR, h: maxB };
}
// The panels of the open plot window; { title, file, blocks: [{ el, title, key }] } (key = a colour key: no panel letter), or null with no such window open.
function pubDomSpec() {
  const live = (e) => e.offsetParent !== null && e.getBoundingClientRect().height > 0;
  if (!$("vio").hidden) return { title: "Expression plots", file: "expression", blocks: [...$("viobody").children].filter((e) => live(e) && e.matches(".vpanel, .dotkey")).map((e) => ({ el: e, key: e.matches(".dotkey") })) };
  if (!$("sigwin").hidden) { const blocks = []; let head = ""; for (const e of $("sigpane").children) { if (e.matches(".sh2")) { head = e.textContent; continue; } if (!live(e) || e.matches(".sctl, .sdesc, .snote, .sbox")) continue; if (e.matches(".smaps")) { [...e.children].filter(live).forEach((q, i) => blocks.push(q.matches(".glomsec") ? { el: q, title: i === 1 ? head : "" } : { el: q, key: true })); head = ""; continue; } if (e.matches("svg, table, canvas") || e.querySelector("svg, table, canvas")) { blocks.push({ el: e, title: head }); head = ""; } } return { title: "Signatures", file: "signatures", blocks }; }
  if (!$("nbr").hidden) return { title: "Neighbours", file: "neighbours", blocks: [...$("nbrbody").children].filter((e) => live(e) && !e.matches(".nbempty, .muted")).map((e) => ({ el: e, key: e.matches(".nlegend") })) };
  if (!$("nbh").hidden) return { title: "Neighbourhoods in space", file: "neighbourhoods", blocks: [...($("nbhleg").children.length && live($("nbhleg")) ? [{ el: $("nbhleg"), key: true }] : []), ...[...$("nbhgrid").querySelectorAll(".glomsec")].filter(live).map((e) => ({ el: e }))] };
  return null;
}
function pubBuildDom(o, dpi, spec) {
  const M = PUB.M * o.F, W = mm2pt(o.wmm), ink = o.dark ? [255, 255, 255] : [0, 0, 0], bg = o.dark ? PUB.DARK : PUB.BG, pbg = o.panel === "dark" ? PUB.DARK : [255, 255, 255], B = o.size, s = dpi / 72, ops = [], th = pubTheme(o.panel !== "dark", pbg);
  const caps = spec.blocks.map((b) => ({ ...pubCapture(b.el, th), title: b.title || "", key: !!b.key })).filter((c) => c.w > 0 && c.h > 0), hdr = o.title ? (B + 1) + 4 * o.F : 0, gap = mm2pt(2), avail = Math.max(40, mm2pt(o.maxh) - hdr - 2);
  const head = (c) => (c.key ? 0 : o.letter || c.title ? (B + 1) + 2 * o.F : 0), maxW = Math.max(...caps.map((c) => c.w));
  const pack = (k) => { let x = 0, y = 0, rh = 0; const pos = []; for (const c of caps) { const w = c.w * k, h = c.h * k + head(c); if (x > 0 && x + w > W) { x = 0; y += rh + gap; rh = 0; } pos.push({ x, y }); x += w + gap; rh = Math.max(rh, h); } return { pos, H: y + rh }; };
  let k = Math.min(W / maxW, o.size / 10.5); const floor = Math.min(k, 5 / 10.5); let P = pack(k);   // on-screen text of about 10.5 px becomes the chosen text size; shrinks to fit the height, down to about 5 pt
  while (P.H > avail && k > floor) { k = Math.max(floor, k * 0.94); P = pack(k); }
  if (o.title) ops.push(pubT(spec.title, M, B + 1, B, false, ink));
  let letter = 0;
  caps.forEach((c, i) => {
    const x = P.pos[i].x, y = hdr + P.pos[i].y, hh = head(c);
    if (hh) { let tx = x + 1; if (o.letter) { const l = pubLetter(o.letter, letter++); ops.push(pubT(l, tx, y + B + 1, B + 1, true, ink)); tx += pubTW(l, B + 1, true) + 4 * o.F; } if (c.title) { let t = c.title; const room = c.w * k - (tx - x); while (t.length > 4 && pubTW(t, B) > room) t = t.slice(0, -2); if (t !== c.title) t += "…"; ops.push(pubT(t, tx, y + B + 1, B, false, ink)); } }
    const yb = y + hh; ops.push(pubR(x, yb, c.w * k, c.h * k, pbg));
    for (const q of c.ops) {
      if (q.t === "t") ops.push({ ...q, x: x + q.x * k, y: yb + q.y * k, size: q.size * k });
      else if (q.t === "r") ops.push({ ...q, x: x + q.x * k, y: yb + q.y * k, w: q.w * k, h: q.h * k });
      else if (q.t === "p") ops.push({ ...q, pts: q.pts.map(([a, b2]) => [x + a * k, yb + b2 * k]), lw: q.lw * k });
      else ops.push({ t: "i", c: q.fig ? q.fig(Math.max(40, Math.round(Math.max(q.w, q.h) * k * s))) : q.c, x: x + q.x * k, y: yb + q.y * k, w: q.w * k, h: q.h * k, smooth: q.smooth });
    }
  });
  const H = hdr + P.H + 2;
  return { ops, W, H, bg, dpi, px: [Math.round(W * s), Math.round(H * s)], what: spec.title, file: spec.file, dom: true };
}

// The embedding panel of a figure: all cells as small dots (the viewer's colours, filters and hidden categories; filtered-out cells faint grey). The legend names the categories.
async function pubEmbedding(cve, ex, ey, Wp, imgH, o, ink, pbg, s, ops, B) {
  const kind = $("umemb").value, px = +$("umpcx").value || 1, py = +$("umpcy").value || 2, mode = S.meta.cats[S.colorBy] ? S.colorBy : "cell_type";
  const name = kind === "pca" ? `pca-${px}-${py}` : "umap", url = kind === "pca" ? `/api/pca/${px}/${py}` : "/api/umap/harmony";
  const fkeys = Object.keys(S.filters).filter((k) => S.filters[k].size && S.meta.cats[k]), hid = S.hiddenCats[mode] || new Set();
  if (!UM.emb[name]) UM.emb[name] = await bget(url, Float32Array);
  for (const k of new Set([mode, ...fkeys])) if (!UM.cat[k]) UM.cat[k] = await bget(`/api/cat_all/${encodeURIComponent(k)}`, Uint16Array);
  const xm = S.exprFilters.length ? await exprMaskAll() : null, emb = UM.emb[name], codes = UM.cat[mode], n = codes.length;
  const keep = new Uint8Array(n).fill(1); applyCatFilters(keep, n, fkeys, (k) => UM.cat[k]);
  if (xm) for (let i = 0; i < n; i++) if (keep[i] && !xm[i]) keep[i] = 0;
  if (hid.size) for (let i = 0; i < n; i++) if (keep[i] && hid.has(codes[i])) keep[i] = 0;
  const R = kind === "pca" ? null : umRun(emb, n), bk = R ? `${name}|${R.run}` : name;
  let b = R ? umBounds(emb, n, bk, R) : UM.bounds[name];
  if (R && !b) return;                                        // this run has no UMAP: the panel stays empty
  if (!b) { const xs = new Float32Array(Math.ceil(n / 7)), ys = new Float32Array(xs.length); for (let i = 0, j = 0; i < n; i += 7, j++) { xs[j] = emb[2 * i]; ys[j] = emb[2 * i + 1]; } xs.sort(); ys.sort(); const q = (a, f) => a[Math.min(a.length - 1, Math.floor(f * a.length))]; b = UM.bounds[name] = [q(xs, 0.005), q(xs, 0.995), q(ys, 0.005), q(ys, 0.995)]; }
  const W = cve.width, H = cve.height, k = Math.min(W / (b[1] - b[0]), H / (b[3] - b[2])) * 0.94, ox = W / 2 - k * (b[0] + b[1]) / 2, oy = H / 2 + k * (b[2] + b[3]) / 2;   // y up
  const X = (x) => ox + k * x, Y = (y) => oy - k * y, ctx = cve.getContext("2d"), dot = Math.max(1, 0.5 * s * Math.sqrt(o.F));
  const smp = S.meta.samples.find((z) => z.id === S.sample), lo = smp.offset, hi = smp.offset + smp.n, tissue = o.emb === "tissue";
  const groups = new Map(), GREYC = "rgba(150,155,165,0.25)";
  for (let i = 0; i < n; i++) {
    if (R && !R.mask[i]) continue;                            // another run: its UMAP has its own coordinates
    const c = codes[i], faint = !keep[i] || (tissue && (i < lo || i >= hi));
    const col = faint ? GREYC : c === MISSING ? rgbStr(GRAY) : catColor(mode, c);
    (groups.get(col) || groups.set(col, []).get(col)).push(i);
  }
  for (const [col, list] of [...groups.entries()].sort((p, q) => (p[0] === GREYC ? -1 : 0) - (q[0] === GREYC ? -1 : 0))) { ctx.fillStyle = col; for (const i of list) { const x = X(emb[2 * i]), y = Y(emb[2 * i + 1]); if (x > -dot && y > -dot && x < W && y < H) ctx.fillRect(x - dot / 2, y - dot / 2, dot, dot); } }
}

// Layout for output resolution `dpi`. Returns { ops, W, H } in points (1/72 in).
function pubBuild(o, dpi) {
  const ds = pubDomSpec(); if (ds?.blocks.length) return pubBuildDom(o, dpi, ds);
  { const spec = pubWinSpec(); if (spec) return pubBuildWin(o, dpi, spec); }
  const M = PUB.M * o.F, W = mm2pt(o.wmm), dark = o.dark, ink = dark ? [255, 255, 255] : [0, 0, 0], bg = dark ? PUB.DARK : PUB.BG, pbg = o.panel === "dark" ? PUB.DARK : [255, 255, 255], B = o.size, ops = [];
  const b = { ...mainBundle(), xy: S.xy }, sc = makeScale([b]), info = sc.info, s = dpi / 72;
  const two = !!HEW.main.side, snap = two && scanMode(), wEmb = !!o.emb, np = 1 + (two ? 1 : 0) + (wEmb ? 1 : 0), gap = np > 1 ? mm2pt(2) : 0, Wp = (W - gap * (np - 1)) / np;   // images beside the cells: a second panel of the same size
  const hdr = o.letter || o.title ? (B + 1) + 4 * o.F : 0, leg = pubLegend(o, W, ink, info), top = hdr, availH = Math.max(40, mm2pt(o.maxh) - top - 4 - leg.h - 2);
  if (o.letter) for (let q = 0; q < np; q++) ops.push(pubT(pubLetter(o.letter, q), q * (Wp + gap) + M, B + 1, B + 1, true, ink));   // one letter per panel
  if (o.title) { ops.push(pubT(`${dlabel(S.sample)} · ${fmt(S.n)} cells`, M + (o.letter ? pubTW(o.letter, B + 1, true) + 6 : 0), B + 1, B, false, ink)); }
  // region (µm window): what the viewer shows now, or the whole tissue
  let cx, cy, k, imgH; const el = $("deckwrap"), vs = S.viewState;
  if (o.region === "view" && vs) { const vw = el.clientWidth / 2 ** vs.zoom; cx = vs.target[0]; cy = -vs.target[1]; k = Wp / vw; imgH = Math.min(Wp * el.clientHeight / el.clientWidth, availH); }
  else { const bd = bounds(), w = (bd.x[1] - bd.x[0]) * 1.04, h = (bd.y[1] - bd.y[0]) * 1.04; cx = (bd.x[0] + bd.x[1]) / 2; cy = (bd.y[0] + bd.y[1]) / 2; k = Math.min(Wp / w, availH / h); imgH = Math.min(availH, h * k); }
  const cw = Math.round(Wp * s), ch = Math.max(8, Math.round(imgH * s)), ks = k * s, x0 = cx - cw / 2 / ks, x1 = cx + cw / 2 / ks, y0 = cy - ch / 2 / ks, y1 = cy + ch / 2 / ks;
  // the window of the panel with the images: the cells' own, or (when the view is not linked) what the side view shows now
  let ks2 = ks, X0b = x0, Y0b = y0;
  if (two && !snap && o.region === "view" && HEX.vs && !HEW.main.link) { const el2 = $("hewrap"), vw2 = el2.clientWidth / 2 ** HEX.vs.zoom; ks2 = (Wp / vw2) * s; X0b = HEX.vs.target[0] - cw / 2 / ks2; Y0b = -HEX.vs.target[1] - ch / 2 / ks2; }
  const cv = document.createElement("canvas"); cv.width = cw; cv.height = ch; const ctx = cv.getContext("2d");
  ctx.fillStyle = rgbStr(pbg); ctx.fillRect(0, 0, cw, ch);
  const ptDraw = (c, ax0, ay0, aks) => {   // the arrows placed on this tissue, in this panel's window
    const l = PT.list[S.sample]; if (!l || !l.length || !PT_IMG.complete) return;
    const sz = Math.max(24, 0.05 * cw), f = sz / 48;
    for (const p of l) c.drawImage(PT_IMG, (p.x - ax0) * aks - PT_ICON.anchorX * f, (-p.y - ay0) * aks - PT_ICON.anchorY * f, sz, sz);
  };
  const paintCells = (tc = ctx) => {
    const m = 15, groups = new Map(); let vis = 0;
    for (let i = 0; i < S.n; i++) {
      const x = S.xy[i * 2], y = S.xy[i * 2 + 1]; if (x < x0 - m || x > x1 + m || y < y0 - m || y > y1 + m) continue; vis++;
      const col = S.vis[i] ? (S.grey ? rgbStr(GREY) : rgbStr(sc.color(b, i))) : `rgba(105,110,122,${S.dimA})`;
      (groups.get(col) || groups.set(col, []).get(col)).push(i);
    }
    const order = [...groups.entries()].sort((p, q) => (p[0].startsWith("rgba") ? -1 : 0) - (q[0].startsWith("rgba") ? -1 : 0));
    const usePoly = !!S.poly && S.shape !== "centroid" && vis <= 400000 && 12 * ks >= 10, px = Math.max(1.2, 8 * ks), sh = usePoly ? 0 : px / 2;
    if (HEW.main.cells) drawCells(tc, order, b, (x) => (x - x0) * ks - sh, (y) => (y - y0) * ks - sh, cw, ch, px, usePoly ? pubPoly() : null, HEW.main.on ? HEW.main.a : 1, HEW.main.outl, 0.25 * s / 0.8, 400000);
    ptDraw(tc, x0, y0, ks);
  };
  const needImgs = (HEW.main.on || (two && !snap)) && HE.data && HE.data[S.sample];
  const imgsInto = async (c, ax0, ay0, aks, side) => {   // the ticked images of this tissue into c, once every tile has arrived (at most ~10 s)
    const he = { sample: S.sample, redo: () => {}, w: "main", side: !!side, miss: { n: 0 } };
    for (let t = 0; t < 40; t++) { c.fillStyle = rgbStr(pbg); c.fillRect(0, 0, cw, ch); he.miss.n = 0; drawHE(c, he, (x) => (x - ax0) * aks, (y) => (y - ay0) * aks, cw, ch, aks); if (!he.miss.n) break; await new Promise((r) => setTimeout(r, 250)); }
  };
  let cv2 = null, fill = null;
  if (two) {
    cv2 = document.createElement("canvas"); cv2.width = cw; cv2.height = ch; const c2 = cv2.getContext("2d"); c2.fillStyle = rgbStr(pbg); c2.fillRect(0, 0, cw, ch);
    if (snap) { try { HEX.deck.redraw(true); const sv = HEX.deck.getCanvas(), r = Math.min(cw / sv.width, ch / sv.height); c2.drawImage(sv, (cw - sv.width * r) / 2, (ch - sv.height * r) / 2, sv.width * r, sv.height * r); } catch (e) { /* no side pane */ } }
  }
  if (HEW.main.on && needImgs) fill = async () => { await imgsInto(ctx, x0, y0, ks); paintCells(); if (two && !snap) { await imgsInto(cv2.getContext("2d"), X0b, Y0b, ks2, true); ptDraw(cv2.getContext("2d"), X0b, Y0b, ks2); } };
  else { paintCells(); if (two && !snap && needImgs) fill = async () => { const c2 = cv2.getContext("2d"); await imgsInto(c2, X0b, Y0b, ks2, true); ptDraw(c2, X0b, Y0b, ks2); }; }
  if (wEmb) {   // the embedding as the last panel; drawn when the data are in (fill)
    const cve = document.createElement("canvas"); cve.width = cw; cve.height = ch; const ex = (two ? 2 : 1) * (Wp + gap), cc = cve.getContext("2d"); cc.fillStyle = rgbStr(pbg); cc.fillRect(0, 0, cw, ch);
    ops.push({ t: "i", c: cve, x: ex, y: top, w: Wp, h: imgH, smooth: true });
    const f0 = fill; fill = async () => { if (f0) await f0(); await pubEmbedding(cve, ex, top, Wp, imgH, o, ink, pbg, s, ops, B); };
  }
  ops.push({ t: "i", c: cv, x: 0, y: top, w: Wp, h: imgH, smooth: false });
  if (two) ops.push({ t: "i", c: cv2, x: Wp + gap, y: top, w: Wp, h: imgH, smooth: true });
  if (o.sb) {   // scale bar: about 20% of a panel's width, rounded to 1/2/5 x 10^n (one per panel, each with its own scale)
    const ifOnly = (side) => { const l = HE.data && HE.data[S.sample] ? heSel(S.sample, "main", side) : []; return l.length > 0 && l.every((e) => e.key.startsWith("if:")); };   // IF is drawn on black: the scale bar turns white there
    const bar = (kk, x, onImg, side) => {
      const target = 0.2 * Wp / kk, p = 10 ** Math.floor(Math.log10(target)), nice = [1, 2, 5, 10].map((q) => q * p).reduce((a, c) => (Math.abs(c - target) < Math.abs(a - target) ? c : a));
      const lab = nice >= 1000 ? nice / 1000 + " mm" : nice + " µm", sy = top + imgH - 6 * o.F, sc2 = o.panel === "dark" || (onImg && ifOnly(side)) ? [255, 255, 255] : [0, 0, 0];
      ops.push(pubR(x + M + 3 * o.F, sy, nice * kk, 2 * o.F, sc2), pubT(lab, x + M + 3 * o.F, sy - 2.5 * o.F, Math.max(5, B - 1), false, sc2));
    };
    bar(k, 0, HEW.main.on, false);
    if (two) bar(snap ? Wp / ($("hewrap").clientWidth / 2 ** (HEX.vs ? HEX.vs.zoom : 0)) : ks2 / s, Wp + gap, !snap, true);
  }
  const ly = top + imgH + 4; for (const op of leg.ops) { op.y += ly; ops.push(op); }
  const H = ly + leg.h + 2;
  return { ops, W, H, bg, dpi, px: [Math.round(W * s), Math.round(H * s)], fill, note: ds ? `The ${ds.title} window has nothing to draw yet, so this is the main view.` : "" };
}

function pubCanvas(plan) {   // rasterise the layout at plan.dpi
  const s = plan.dpi / 72, c = document.createElement("canvas"); c.width = Math.round(plan.W * s); c.height = Math.round(plan.H * s);
  const x = c.getContext("2d"); x.fillStyle = rgbStr(plan.bg); x.fillRect(0, 0, c.width, c.height); x.scale(s, s); x.textBaseline = "alphabetic";
  for (const op of plan.ops) {
    if (op.t === "r") { x.fillStyle = rgbStr(op.c); x.fillRect(op.x, op.y, op.w, op.h); }
    else if (op.t === "t") { x.fillStyle = rgbStr(op.c); x.font = `${op.bold ? "bold " : ""}${op.size}px ${PUB.FONT}`; if (op.rot) { x.save(); x.translate(op.x, op.y); x.rotate(op.rot * Math.PI / 180); x.fillText(op.s, 0, 0); x.restore(); } else x.fillText(op.s, op.x, op.y); }
    else if (op.t === "p") { x.beginPath(); op.pts.forEach(([a, b], i) => (i ? x.lineTo(a, b) : x.moveTo(a, b))); if (op.close) x.closePath(); if (op.fill) { x.fillStyle = rgbStr(op.fill); x.fill(); } if (op.stroke && op.lw > 0) { x.strokeStyle = rgbStr(op.stroke); x.lineWidth = op.lw; x.stroke(); } }
    else { x.imageSmoothingEnabled = op.smooth; x.drawImage(op.c, op.x, op.y, op.w, op.h); }
  }
  return c;
}
async function pubPNG(plan) {   // PNG bytes with a pHYs chunk so the file records its physical size (pixels per metre)
  const blob = await new Promise((r) => pubCanvas(plan).toBlob(r, "image/png")), u8 = new Uint8Array(await blob.arrayBuffer());
  const ppm = Math.round(plan.dpi / 0.0254), chunk = new Uint8Array(21), dv = new DataView(chunk.buffer);
  dv.setUint32(0, 9); chunk.set([0x70, 0x48, 0x59, 0x73], 4); dv.setUint32(8, ppm); dv.setUint32(12, ppm); chunk[16] = 1; dv.setUint32(17, crc32(chunk.subarray(4, 17)));
  const out = new Uint8Array(u8.length + 21); out.set(u8.subarray(0, 33)); out.set(chunk, 33); out.set(u8.subarray(33), 54);   // after the 33-byte signature + IHDR
  return new Blob([out], { type: "image/png" });
}
async function pubPDF(plan) {
  const l1 = (t) => Uint8Array.from(t, (c) => c.charCodeAt(0) & 255), parts = [], offs = []; let pos = 0;
  const put = (x) => { const u = typeof x === "string" ? l1(x) : x; parts.push(u); pos += u.length; };
  const n2 = (v) => (+v).toFixed(2), col = (c) => c.map((v) => (v / 255).toFixed(3)).join(" "), H = plan.H;
  const PDF_CH = { "–": "\x96", "—": "\x97", "…": "\x85", "−": "-", "σ": "sigma ", "ρ": "rho", "≤": "<=", "≥": ">=", "→": "->", "←": "<-", "⁴": "4" };   // Helvetica has no Greek or arrows
  const esc = (t) => t.replace(/[–—…−σρ≤≥→←⁴]/g, (c) => PDF_CH[c]).replace(/[^\x20-\xFF]/g, "-").replace(/([\\()])/g, "\\$1");
  const imgs = []; let content = `${col(plan.bg)} rg 0 0 ${n2(plan.W)} ${n2(H)} re f\n`;
  for (const op of plan.ops) {
    if (op.t === "r") content += `${col(op.c)} rg ${n2(op.x)} ${n2(H - op.y - op.h)} ${n2(op.w)} ${n2(op.h)} re f\n`;
    else if (op.t === "t") { const a = (op.rot || 0) * Math.PI / 180; content += `BT /F${op.bold ? 2 : 1} ${n2(op.size)} Tf ${col(op.c)} rg ${op.rot ? `${Math.cos(a).toFixed(4)} ${(-Math.sin(a)).toFixed(4)} ${Math.sin(a).toFixed(4)} ${Math.cos(a).toFixed(4)} ${n2(op.x)} ${n2(H - op.y)} Tm` : `${n2(op.x)} ${n2(H - op.y)} Td`} (${esc(op.s)}) Tj ET\n`; }
    else if (op.t === "p") { const st = op.stroke && op.lw > 0; content += `${op.fill ? col(op.fill) + " rg " : ""}${st ? `${col(op.stroke)} RG ${n2(op.lw)} w ` : ""}` + op.pts.map(([a, b], i) => `${n2(a)} ${n2(H - b)} ${i ? "l" : "m"}`).join(" ") + `${op.close ? " h" : ""} ${op.fill && st ? "B" : op.fill ? "f" : "S"}\n`; }
    else { imgs.push(op); content += `q ${n2(op.w)} 0 0 ${n2(op.h)} ${n2(op.x)} ${n2(H - op.y - op.h)} cm /Im${imgs.length} Do Q\n`; }
  }
  const obj = (n, body) => { offs[n] = pos; put(`${n} 0 obj\n`); put(body); put("\nendobj\n"); };
  const stream = (n, dict, data) => { offs[n] = pos; put(`${n} 0 obj\n<< ${dict} /Length ${data.length} >>\nstream\n`); put(data); put("\nendstream\nendobj\n"); };
  put("%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
  obj(1, "<< /Type /Catalog /Pages 2 0 R >>"); obj(2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
  obj(3, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${n2(plan.W)} ${n2(H)}] /Resources << /Font << /F1 5 0 R /F2 6 0 R >> /XObject << ${imgs.map((_, i) => `/Im${i + 1} ${7 + i} 0 R`).join(" ")} >> >> /Contents 4 0 R >>`);
  stream(4, "", l1(content));
  obj(5, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"); obj(6, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>");
  for (let i = 0; i < imgs.length; i++) {
    const op = imgs[i], w = op.c.width, h = op.c.height, d = op.c.getContext("2d").getImageData(0, 0, w, h).data, rgb = new Uint8Array(w * h * 3);
    for (let p = 0, q = 0; p < d.length; p += 4) { rgb[q++] = d[p]; rgb[q++] = d[p + 1]; rgb[q++] = d[p + 2]; }
    let data = rgb, filt = "";
    if (typeof CompressionStream !== "undefined") { data = new Uint8Array(await new Response(new Blob([rgb]).stream().pipeThrough(new CompressionStream("deflate"))).arrayBuffer()); filt = "/Filter /FlateDecode "; }
    stream(7 + i, `/Type /XObject /Subtype /Image /Width ${w} /Height ${h} /ColorSpace /DeviceRGB /BitsPerComponent 8 ${op.smooth ? "/Interpolate true " : ""}${filt}`, data);
  }
  const N = 7 + imgs.length, xref = pos; put(`xref\n0 ${N}\n0000000000 65535 f \n`);
  for (let n = 1; n < N; n++) put(`${String(offs[n]).padStart(10, "0")} 00000 n \n`);
  put(`trailer\n<< /Size ${N} /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF\n`);
  return new Blob(parts, { type: "application/pdf" });
}

// ---- dialog
function pubOpts() {
  const v = (id) => $(id).value, size = Math.max(5, Math.min(48, +v("pubsize") || 7));
  const hasEmb = !!S.meta?.umaps?.length;   // an open Embeddings window adds its panel to the figure (the Embedding panel menu still picks all cells / this tissue highlighted)
  return { wmm: +v("pubw"), maxh: Math.max(40, +v("pubmaxh") || 170), dpi: +v("pubdpi"), region: v("pubreg"), letter: v("publetter").trim(), title: $("pubtitle").checked, sb: $("pubsb").checked, dark: $("pubdark").checked, panel: $("pubpanel").value, emb: !hasEmb ? "" : $("pubemb").value || (S.umOpen && $("umapwin").clientWidth > 0 ? "all" : ""), size, F: Math.max(1, size / 7) };
}
let pubT_ = 0, pubTok = 0;
function pubPreview() {
  clearTimeout(pubT_);
  pubT_ = setTimeout(() => {
    if ($("pubwin").hidden || !S.meta) return;
    S._inFig = true; let plan; const o = pubOpts(); try { plan = pubBuild(o, 110); } finally { S._inFig = false; }
    const pv = $("pubprev"), show = () => { const c = pubCanvas(plan); pv.width = c.width; pv.height = c.height; pv.getContext("2d").drawImage(c, 0, 0); };
    show(); if (plan.fill) { const tk = ++pubTok; (async () => { S._inFig = true; try { await plan.fill(); } finally { S._inFig = false; } })().then(() => { if (tk === pubTok && !$("pubwin").hidden) show(); }); }   // the images beside / under the cells arrive after the first draw
    $("pubwhat").textContent = plan.what ? `Figure of the open window: ${plan.what}` : plan.note || ""; $("pubwhat").hidden = !plan.what && !plan.note; $("pubcols").hidden = !!plan.dom;   // plot windows keep the colours shown in them   // a window's figure: its own tiles, so the region and embedding choices do not apply
    $("pubreg").closest("label").hidden = !!plan.what; $("pubembl").hidden = !!plan.what || !S.meta?.umaps?.length;
    const f = o.dpi / 72, px = [Math.round(plan.W * f), Math.round(plan.H * f)];
    $("pubinfo").textContent = `${o.wmm} × ${(plan.H * 25.4 / 72).toFixed(0)} mm · ${px[0]} × ${px[1]} px at ${o.dpi} dpi · text ${o.size} pt (legend ${Math.max(5, o.size - 1)} pt), panel letter ${o.size + 1} pt bold`;
  }, 60);
}
async function pubSave(kind) {
  const o = pubOpts(), b = $(kind === "pdf" ? "pubpdf" : "pubpng"), t = b.textContent; b.disabled = true; b.textContent = "Rendering…";
  try {
    S._inFig = true; let plan; try { plan = pubBuild(o, o.dpi); if (plan.fill) await plan.fill(); } finally { S._inFig = false; }
    const blob = kind === "pdf" ? await pubPDF(plan) : await pubPNG(plan), a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = `figure_${plan.file || S.sample}_${o.wmm}mm_${o.dpi}dpi.${kind}`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 4000);
  } catch (e) { alert("Could not render the figure: " + e.message + "\nTry 300 dpi or a smaller width."); }
  b.disabled = false; b.textContent = t;
}
// Colours of the categories of the current colour-by, chosen for the figure (see ovHex in app.js)
const pubHex = (col) => { const x = document.createElement("canvas").getContext("2d"); x.fillStyle = col; return x.fillStyle; };
function pubColUI() {
  const body = $("pubcolbody"), k = S.colorBy; body.replaceChildren();
  if (!S.meta || !S.meta.cats[k]) { body.append(Object.assign(document.createElement("div"), { className: "small muted", textContent: "Colours can be chosen when the colour-by is a category (e.g. Cell type)." })); return; }
  const L = catMeta(k).levels, ov = (S.colOv ||= {})[k] || (S.colOv[k] = {});
  const row = document.createElement("div"); row.className = "pubcolhead small muted"; row.textContent = `${colorName()}: pick a colour per category`; body.append(row);
  const grid = document.createElement("div"); grid.className = "pubcolgrid";
  L.slice(0, 80).forEach((_, i) => {
    const lab = document.createElement("label"), inp = document.createElement("input"); inp.type = "color"; inp.value = pubHex(ov[i] || catColor0(k, i));
    inp.addEventListener("input", () => { ov[i] = inp.value; S.pubLegKey = null; pubAfterColour(); });
    lab.append(inp, " " + display(k, i)); grid.append(lab);
  });
  const rs = document.createElement("button"); rs.type = "button"; rs.textContent = "Reset colours"; rs.addEventListener("click", () => { S.colOv[k] = {}; S.pubLegKey = null; pubColUI(); pubAfterColour(); });
  const all = document.createElement("label"); all.className = "chk"; const ca = document.createElement("input"); ca.type = "checkbox"; ca.checked = !!S.colOvAll;
  ca.addEventListener("change", () => { S.colOvAll = ca.checked; pubAfterColour(); }); all.append(ca, " also use these colours in the viewer and windows");
  body.append(grid, rs, all);
}
function pubAfterColour() { pubPreview(); if (S.colOvAll) { render(); refreshWindows(); } }
// Presets: one click sets size, resolution, text size and background. Paper = journal 2-column width, 300 dpi, 7 pt text (Nature style);
// presentation = slide-wide (10 in) image, 18 pt text, dark, 200 dpi (screens / projectors); poster = 400 mm wide, 28 pt text (readable from ~1.5 m), 200 dpi.
const PUB_PRESETS = {
  paper: { pubw: "183", pubmaxh: 170, pubdpi: "300", pubsize: 7, publetter: "a", pubdark: false, pubpanel: "white", pubtitle: true },
  presentation: { pubw: "254", pubmaxh: 143, pubdpi: "200", pubsize: 18, publetter: "", pubdark: true, pubpanel: "dark", pubtitle: true },
  poster: { pubw: "400", pubmaxh: 500, pubdpi: "200", pubsize: 28, publetter: "", pubdark: false, pubpanel: "white", pubtitle: true },
};
function pubApplyPreset(k) { const p = PUB_PRESETS[k]; if (!p) return; for (const [id, v] of Object.entries(p)) { const e = $(id); if (e.type === "checkbox") e.checked = v; else e.value = v; } pubPreview(); }
const pubCustom = () => { $("pubpreset").value = ""; };   // a setting changed by hand: the figure no longer follows the preset
(function pubInit() {
  $("pubpreset").addEventListener("change", (e) => pubApplyPreset(e.target.value));
  $("pubbtn").addEventListener("click", () => { $("pubembl").hidden = !S.meta?.umaps?.length; $("pubwin").hidden = false; pubColUI(); pubPreview(); });
  for (const id of ["pubw", "pubmaxh", "pubdpi", "publetter", "pubtitle", "pubdark", "pubpanel", "pubsize"]) $(id).addEventListener("input", pubCustom);
  $("pubx").addEventListener("click", () => ($("pubwin").hidden = true));
  for (const id of ["pubw", "pubmaxh", "pubdpi", "pubreg", "publetter", "pubtitle", "pubsb", "pubdark", "pubpanel", "pubemb", "pubsize"]) $(id).addEventListener("input", pubPreview);
  $("pubpng").addEventListener("click", () => pubSave("png")); $("pubpdf").addEventListener("click", () => pubSave("pdf"));
  window.addEventListener("keydown", (e) => { if (e.key === "Escape" && !$("pubwin").hidden) { $("pubwin").hidden = true; e.stopPropagation(); } }, true);
})();
