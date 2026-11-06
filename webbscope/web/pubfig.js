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
function pubLegend(o, W, ink, info, tickFor) {   // ops relative to y = 0; returns { ops, h }
  const M = PUB.M, B = o.size, sm = Math.max(5, B - 1), ops = [], mode = S.colorBy; let y = 0;
  const cap = (t) => { ops.push(pubT(t, M, y + B, B, true, ink)); y += B * 1.5; };
  if (info.type === "none") return { ops, h: 0 };
  if (info.type === "cat") {
    cap(colorName());
    const { out, miss } = countsFor(mode), items = out.map((n, i) => [i, n]).filter(([, n]) => n > 0).sort((a, b) => b[1] - a[1]).map(([i]) => [display(mode, i), hex2rgb(catColor(mode, i))]);
    if (miss) items.push(["missing", GRAY]);
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
  const lay = S.layer === "counts" ? "raw counts" : `log1p(CP${fmt(S.meta.norm_scale)})`;
  cap(info.type === "num" ? `${labelOf(mode)}${info.log ? " (log scale)" : ""}` : `${info.label} (${lay})`);
  const L = info.type === "gene" && info.cont && info.div ? DIVL : LUT, g = document.createElement("canvas"); g.width = 256; g.height = 1;
  const gx = g.getContext("2d"); L.forEach((c, i) => { gx.fillStyle = rgbStr(c); gx.fillRect(i, 0, 1, 1); });
  const bw = Math.min(W - 2 * M - 24, mm2pt(50)), bx = M + 12;   // 12 pt of room each side so the end labels are not clipped
  ops.push({ t: "i", c: g, x: bx, y, w: bw, h: B * 0.8, smooth: true }); y += B * 0.8;
  for (const k of scaleTicks(info)) {
    const s = fmtTick(k.v), w = pubTW(s, sm); ops.push(pubT(s, bx + bw * k.t - w / 2, y + sm + 1.5, sm, false, ink));
    ops.push(pubR(bx + bw * k.t - 0.125, y, 0.25, 1.5, ink));
  }
  return { ops, h: y + sm * 1.8 };
}

// Layout for output resolution `dpi`. Returns { ops, W, H } in points (1/72 in).
function pubBuild(o, dpi) {
  const M = PUB.M, W = mm2pt(o.wmm), dark = o.dark, ink = dark ? [255, 255, 255] : [0, 0, 0], bg = dark ? PUB.DARK : PUB.BG, B = o.size, ops = [];
  const b = { ...mainBundle(), xy: S.xy }, sc = makeScale([b]), info = sc.info, s = dpi / 72;
  const two = !!HEW.main.side, snap = two && scanMode(), gap = two ? mm2pt(2) : 0, Wp = two ? (W - gap) / 2 : W;   // images beside the cells: a second panel of the same size
  const hdr = o.letter || o.title ? 8 + 4 : 0, leg = pubLegend(o, W, ink, info), top = hdr, availH = Math.max(40, mm2pt(o.maxh) - top - 4 - leg.h - 2);
  if (o.letter) ops.push(pubT(o.letter, M, 8, 8, true, ink));
  if (o.title) { ops.push(pubT(`${dlabel(S.sample)} · ${fmt(S.n)} cells`, M + (o.letter ? pubTW(o.letter, 8, true) + 6 : 0), 8, B, false, ink)); }
  // region (µm window): what the viewer shows now, or the whole tissue
  let cx, cy, k, imgH; const el = $("deckwrap"), vs = S.viewState;
  if (o.region === "view" && vs) { const vw = el.clientWidth / 2 ** vs.zoom; cx = vs.target[0]; cy = -vs.target[1]; k = Wp / vw; imgH = Math.min(Wp * el.clientHeight / el.clientWidth, availH); }
  else { const bd = bounds(), w = (bd.x[1] - bd.x[0]) * 1.04, h = (bd.y[1] - bd.y[0]) * 1.04; cx = (bd.x[0] + bd.x[1]) / 2; cy = (bd.y[0] + bd.y[1]) / 2; k = Math.min(Wp / w, availH / h); imgH = Math.min(availH, h * k); }
  const cw = Math.round(Wp * s), ch = Math.max(8, Math.round(imgH * s)), ks = k * s, x0 = cx - cw / 2 / ks, x1 = cx + cw / 2 / ks, y0 = cy - ch / 2 / ks, y1 = cy + ch / 2 / ks;
  // the window of the panel with the images: the cells' own, or (when the view is not linked) what the side view shows now
  let ks2 = ks, X0b = x0, Y0b = y0;
  if (two && !snap && o.region === "view" && HEX.vs && !HEW.main.link) { const el2 = $("hewrap"), vw2 = el2.clientWidth / 2 ** HEX.vs.zoom; ks2 = (Wp / vw2) * s; X0b = HEX.vs.target[0] - cw / 2 / ks2; Y0b = -HEX.vs.target[1] - ch / 2 / ks2; }
  const cv = document.createElement("canvas"); cv.width = cw; cv.height = ch; const ctx = cv.getContext("2d");
  ctx.fillStyle = rgbStr(bg); ctx.fillRect(0, 0, cw, ch);
  const ptDraw = (c, ax0, ay0, aks) => {   // the arrows placed on this tissue, in this panel's window
    const l = PT.list[S.sample]; if (!l || !l.length || !PT_IMG.complete) return;
    const sz = Math.max(24, 0.05 * cw), f = sz / 48;
    for (const p of l) c.drawImage(PT_IMG, (p.x - ax0) * aks - PT_ICON.anchorX * f, (-p.y - ay0) * aks - PT_ICON.anchorY * f, sz, sz);
  };
  const paintCells = () => {
    const m = 15, groups = new Map(); let vis = 0;
    for (let i = 0; i < S.n; i++) {
      const x = S.xy[i * 2], y = S.xy[i * 2 + 1]; if (x < x0 - m || x > x1 + m || y < y0 - m || y > y1 + m) continue; vis++;
      const col = S.vis[i] ? (S.grey ? rgbStr(GREY) : rgbStr(sc.color(b, i))) : "rgba(105,110,122,0.45)";
      (groups.get(col) || groups.set(col, []).get(col)).push(i);
    }
    const order = [...groups.entries()].sort((p, q) => (p[0].startsWith("rgba") ? -1 : 0) - (q[0].startsWith("rgba") ? -1 : 0));
    const usePoly = !!S.poly && S.shape !== "centroid" && vis <= 400000 && 12 * ks >= 10, px = Math.max(1.2, 8 * ks), sh = usePoly ? 0 : px / 2;
    if (HEW.main.cells) drawCells(ctx, order, b, (x) => (x - x0) * ks - sh, (y) => (y - y0) * ks - sh, cw, ch, px, usePoly ? pubPoly() : null, HEW.main.on ? HEW.main.a : 1, HEW.main.outl, 0.25 * s / 0.8, 400000);
    ptDraw(ctx, x0, y0, ks);
  };
  const needImgs = (HEW.main.on || (two && !snap)) && HE.data && HE.data[S.sample];
  const imgsInto = async (c, ax0, ay0, aks, side) => {   // the ticked images of this tissue into c, once every tile has arrived (at most ~10 s)
    const he = { sample: S.sample, redo: () => {}, w: "main", side: !!side, miss: { n: 0 } };
    for (let t = 0; t < 40; t++) { c.fillStyle = rgbStr(bg); c.fillRect(0, 0, cw, ch); he.miss.n = 0; drawHE(c, he, (x) => (x - ax0) * aks, (y) => (y - ay0) * aks, cw, ch, aks); if (!he.miss.n) break; await new Promise((r) => setTimeout(r, 250)); }
  };
  let cv2 = null, fill = null;
  if (two) {
    cv2 = document.createElement("canvas"); cv2.width = cw; cv2.height = ch; const c2 = cv2.getContext("2d"); c2.fillStyle = rgbStr(bg); c2.fillRect(0, 0, cw, ch);
    if (snap) { try { HEX.deck.redraw(true); const sv = HEX.deck.getCanvas(), r = Math.min(cw / sv.width, ch / sv.height); c2.drawImage(sv, (cw - sv.width * r) / 2, (ch - sv.height * r) / 2, sv.width * r, sv.height * r); } catch (e) { /* no side pane */ } }
  }
  if (HEW.main.on && needImgs) fill = async () => { await imgsInto(ctx, x0, y0, ks); paintCells(); if (two && !snap) { await imgsInto(cv2.getContext("2d"), X0b, Y0b, ks2, true); ptDraw(cv2.getContext("2d"), X0b, Y0b, ks2); } };
  else { paintCells(); if (two && !snap && needImgs) fill = async () => { const c2 = cv2.getContext("2d"); await imgsInto(c2, X0b, Y0b, ks2, true); ptDraw(c2, X0b, Y0b, ks2); }; }
  ops.push({ t: "i", c: cv, x: 0, y: top, w: Wp, h: imgH, smooth: false });
  if (two) ops.push({ t: "i", c: cv2, x: Wp + gap, y: top, w: Wp, h: imgH, smooth: true });
  if (o.sb) {   // scale bar: about 20% of a panel's width, rounded to 1/2/5 x 10^n (one per panel, each with its own scale)
    const ifOnly = (side) => { const l = HE.data && HE.data[S.sample] ? heSel(S.sample, "main", side) : []; return l.length > 0 && l.every((e) => e.key.startsWith("if:")); };   // IF is drawn on black: the scale bar turns white there
    const bar = (kk, x, onImg, side) => {
      const target = 0.2 * Wp / kk, p = 10 ** Math.floor(Math.log10(target)), nice = [1, 2, 5, 10].map((q) => q * p).reduce((a, c) => (Math.abs(c - target) < Math.abs(a - target) ? c : a));
      const lab = nice >= 1000 ? nice / 1000 + " mm" : nice + " µm", sy = top + imgH - 6, sc2 = dark || (onImg && ifOnly(side)) ? [255, 255, 255] : [0, 0, 0];
      ops.push(pubR(x + M + 3, sy, nice * kk, 2, sc2), pubT(lab, x + M + 3, sy - 2.5, Math.max(5, B - 1), false, sc2));
    };
    bar(k, 0, HEW.main.on, false);
    if (two) bar(snap ? Wp / ($("hewrap").clientWidth / 2 ** (HEX.vs ? HEX.vs.zoom : 0)) : ks2 / s, Wp + gap, !snap, true);
  }
  const ly = top + imgH + 4; for (const op of leg.ops) { op.y += ly; ops.push(op); }
  const H = ly + leg.h + 2;
  return { ops, W, H, bg, dpi, px: [Math.round(W * s), Math.round(H * s)], fill };
}

function pubCanvas(plan) {   // rasterise the layout at plan.dpi
  const s = plan.dpi / 72, c = document.createElement("canvas"); c.width = Math.round(plan.W * s); c.height = Math.round(plan.H * s);
  const x = c.getContext("2d"); x.fillStyle = rgbStr(plan.bg); x.fillRect(0, 0, c.width, c.height); x.scale(s, s); x.textBaseline = "alphabetic";
  for (const op of plan.ops) {
    if (op.t === "r") { x.fillStyle = rgbStr(op.c); x.fillRect(op.x, op.y, op.w, op.h); }
    else if (op.t === "t") { x.fillStyle = rgbStr(op.c); x.font = `${op.bold ? "bold " : ""}${op.size}px ${PUB.FONT}`; x.fillText(op.s, op.x, op.y); }
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
  const esc = (t) => t.replace(/[^\x20-\xFF]/g, "-").replace(/([\\()])/g, "\\$1");
  const imgs = []; let content = `${col(plan.bg)} rg 0 0 ${n2(plan.W)} ${n2(H)} re f\n`;
  for (const op of plan.ops) {
    if (op.t === "r") content += `${col(op.c)} rg ${n2(op.x)} ${n2(H - op.y - op.h)} ${n2(op.w)} ${n2(op.h)} re f\n`;
    else if (op.t === "t") content += `BT /F${op.bold ? 2 : 1} ${n2(op.size)} Tf ${col(op.c)} rg ${n2(op.x)} ${n2(H - op.y)} Td (${esc(op.s)}) Tj ET\n`;
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
  const v = (id) => $(id).value, size = Math.max(5, Math.min(7, +v("pubsize") || 7));
  return { wmm: +v("pubw"), maxh: Math.max(40, +v("pubmaxh") || 170), dpi: +v("pubdpi"), region: v("pubreg"), letter: v("publetter").trim(), title: $("pubtitle").checked, sb: $("pubsb").checked, dark: $("pubdark").checked, size };
}
let pubT_ = 0, pubTok = 0;
function pubPreview() {
  clearTimeout(pubT_);
  pubT_ = setTimeout(() => {
    if ($("pubwin").hidden || !S.meta) return;
    const o = pubOpts(), plan = pubBuild(o, 110), pv = $("pubprev"), show = () => { const c = pubCanvas(plan); pv.width = c.width; pv.height = c.height; pv.getContext("2d").drawImage(c, 0, 0); };
    show(); if (plan.fill) { const tk = ++pubTok; plan.fill().then(() => { if (tk === pubTok && !$("pubwin").hidden) show(); }); }   // the images beside / under the cells arrive after the first draw
    const f = o.dpi / 72, px = [Math.round(plan.W * f), Math.round(plan.H * f)];
    $("pubinfo").textContent = `${o.wmm} × ${(plan.H * 25.4 / 72).toFixed(0)} mm · ${px[0]} × ${px[1]} px at ${o.dpi} dpi · text ${o.size} pt (legend ${Math.max(5, o.size - 1)} pt), panel letter 8 pt bold`;
  }, 60);
}
async function pubSave(kind) {
  const o = pubOpts(), b = $(kind === "pdf" ? "pubpdf" : "pubpng"), t = b.textContent; b.disabled = true; b.textContent = "Rendering…";
  try {
    const plan = pubBuild(o, o.dpi); if (plan.fill) await plan.fill();
    const blob = kind === "pdf" ? await pubPDF(plan) : await pubPNG(plan), a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = `figure_${S.sample}_${o.wmm}mm_${o.dpi}dpi.${kind}`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 4000);
  } catch (e) { alert("Could not render the figure: " + e.message + "\nTry 300 dpi or a smaller width."); }
  b.disabled = false; b.textContent = t;
}
(function pubInit() {
  $("pubbtn").addEventListener("click", () => { $("pubwin").hidden = false; pubPreview(); });
  $("pubx").addEventListener("click", () => ($("pubwin").hidden = true));
  for (const id of ["pubw", "pubmaxh", "pubdpi", "pubreg", "publetter", "pubtitle", "pubsb", "pubdark", "pubsize"]) $(id).addEventListener("input", pubPreview);
  $("pubpng").addEventListener("click", () => pubSave("png")); $("pubpdf").addEventListener("click", () => pubSave("pdf"));
  window.addEventListener("keydown", (e) => { if (e.key === "Escape" && !$("pubwin").hidden) { $("pubwin").hidden = true; e.stopPropagation(); } }, true);
})();
