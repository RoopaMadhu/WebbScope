"""One-page QC report per tissue (self-contained HTML: inline SVG charts and embedded PNG maps; prints to one A4 page).

Everything here is computed from what a project already stores (counts matrix, cell table, cell boundaries, transcripts coordinates);
nothing is read from the original Xenium output. Quantities that depend on what was recorded at import (negative-control counts per cell,
the transcript quality histogram, the assigned share, negative-control positions) are shown when the project has them; for projects imported
before they were stored the report says "not recorded" and the run needs to be imported again. See docs/QC.md for how each number is calculated.
"""
import base64, datetime, html, io, json, math, re, zlib
from pathlib import Path
import numpy as np
from PIL import Image

TX_BIN_UM = 20.0                         # bin of the transcript-density map (µm); fixed so the map is comparable between tissues
BIN_UM = 25.0                          # grid for the occupied-area estimate
LOW_TX = 20                            # "low-count cell": fewer panel-gene transcripts than this (informational)
FLAG_Z, FLAG_MIN_TISSUES = 3.0, 6      # cohort flags: robust z beyond this, only with at least this many tissues
AREA_COLS, NUC_COLS = ("cell_area", "area_um", "area"), ("nucleus_area",)
_cache = {}                            # id(store) -> {"nnz": per-cell detected genes, "table": per-tissue summary rows}


# ------------------------------------------------------------------------------------------------------------------------------ numbers
def _col(S, names):
    return next((c for c in names if c in S.df.columns and S.df[c].dtype.kind in "fi"), None)


def _nnz_per_cell(S):
    """Genes detected (count > 0) per cell, all cells, cached: bincount of the CSR column indices in chunks (the matrix is genes x cells)."""
    c = _cache.setdefault(id(S), {})
    if "nnz" not in c:
        out = np.zeros(S.n, np.int64)
        for a in range(0, len(S.indices), 20_000_000):
            out += np.bincount(np.asarray(S.indices[a:a + 20_000_000]), minlength=S.n)
        c["nnz"] = out.astype(np.int32)
    return c["nnz"]


def _area_mm2(x, y):
    ix = np.floor((x - x.min()) / BIN_UM).astype(np.int64)
    iy = np.floor((y - y.min()) / BIN_UM).astype(np.int64)
    return len(np.unique(ix * (iy.max() + 1) + iy)) * BIN_UM * BIN_UM / 1e6


def _boundary_counts(S, a, b):
    if S.poly_off is None:
        return None
    n = np.diff(np.asarray(S.poly_off[a:b + 1]))
    return int((n == 0).sum()), int(((n > 0) & (n < 3)).sum())


def _manifest(S):
    c = _cache.setdefault(id(S), {})
    if "man" not in c:
        try:
            c["man"] = json.loads((Path(S.cache_dir) / "cache_info.json").read_text()).get("manifest", {})
        except Exception:
            c["man"] = {}
    return c["man"]


CTRL_LABEL = {"nc_probe": "negative control probes", "nc_codeword": "negative control codewords", "unassigned_codeword": "unassigned codewords",
              "deprecated_codeword": "deprecated codewords"}


def controls(S, sample, lib=None):
    """Control-feature numbers of a tissue from the per-cell columns written at import (None when the project has none):
    per type {n_features, mean per cell, % cells with >= 1}, the per-cell negative-control total (probes + codewords) and neg_ratio =
    negative-control counts per control feature as a percentage of the mean counts per gene (an estimate of the background level)."""
    a, b = S.slices[sample]
    run = str(S.df["run"].values[a]) if "run" in S.df.columns else ""
    feats = (_manifest(S).get("controls") or {}).get(run) or (_manifest(S).get("controls") or {}).get("") or {}
    out, neg, nneg = {}, np.zeros(b - a), 0
    for key in CTRL_LABEL:
        col = f"{key}_counts"
        if col in S.df.columns:
            v = S.df[col].values[a:b].astype(float)
            out[key] = {"n_features": feats.get(key), "mean": float(v.mean()), "pct_cells": float((v >= 1).mean() * 100)}
            if key in ("nc_probe", "nc_codeword"):
                neg += v; nneg += feats.get(key) or 0
    if not out:
        return None
    r = {"types": out, "neg_per_cell": neg if ("nc_probe" in out or "nc_codeword" in out) else None, "neg_ratio": None}
    if r["neg_per_cell"] is not None and nneg and lib is not None and lib.mean() > 0:
        r["neg_ratio"] = float(neg.mean() / nneg / (lib.mean() / len(S.genes)) * 100)
    return r


def tx_qc(S, TX, sample):
    """Assignment and quality numbers recorded when the transcripts were imported (None if absent)."""
    q = TX.qc(sample) if TX is not None else None
    if not q:
        return None
    out = {"assigned": q["assigned"], "unassigned_near": q["unassigned_near"], "unassigned_far": q.get("unassigned_far", 0), "run": q["run"], "qv_hist": q["qv_hist"], "gene_qv": q.get("gene_qv"), "tissue_rows": q.get("tissue_rows")}
    tot = q["assigned"] + q["unassigned_near"] + q.get("unassigned_far", 0)
    out["assigned_frac"] = q["assigned"] / tot * 100 if tot else None
    h = q["qv_hist"]
    out["q20_frac"] = (sum(h[20:]) / sum(h) * 100) if h and sum(h) else None
    r = q["run"]
    out["run_assigned_frac"] = r["q_filtered_panel_assigned"] / r["q_filtered_panel"] * 100 if r.get("q_filtered_panel") else None
    return out


def _row(S, sample, nnz, acol, TX=None):
    a, b = S.slices[sample]
    lib = np.asarray(S.lib[a:b]) if len(S.lib) else None
    r = {"sample": sample, "n_cells": b - a, "area_mm2": _area_mm2(S.x[a:b], S.y[a:b])}
    r["density"] = r["n_cells"] / r["area_mm2"] if r["area_mm2"] else float("nan")
    if lib is not None:
        r["median_tx"], r["low_frac"] = float(np.median(lib)), float((lib < LOW_TX).mean())
    r["median_genes"] = float(np.median(nnz[a:b]))
    if acol:
        r["median_area"] = float(np.nanmedian(S.df[acol].values[a:b].astype(float)))
    bc = _boundary_counts(S, a, b)
    if bc:
        r["no_boundary_frac"] = (bc[0] + bc[1]) / max(1, b - a)
    ct = controls(S, sample, lib)
    if ct and ct["neg_ratio"] is not None:
        r["neg_ratio"] = ct["neg_ratio"]
    if TX is not None and sample in TX.loc and r["area_mm2"]:
        r["transcripts_total"] = TX.total(sample)
        r["tx_density"] = r["transcripts_total"] / r["area_mm2"]            # stored (quality-filtered) transcripts per mm² of occupied tissue
        cp = TX.ctrl_points(sample)
        if cp is not None:
            r["nc_density"] = len(cp[0]) / r["area_mm2"]                    # negative-control probe transcripts (QV >= 20) per mm²
    tq = tx_qc(S, TX, sample)
    if tq:
        if tq["assigned_frac"] is not None:
            r["assigned_frac"] = tq["assigned_frac"]
        if tq["q20_frac"] is not None:
            r["q20_frac"] = tq["q20_frac"]
    return r


def table(S, TX=None):
    """One row per tissue of the project (cached; the transcripts-based columns need TX)."""
    c = _cache.setdefault(id(S), {})
    key = ("table", TX is not None)
    if key not in c:
        nnz, acol = _nnz_per_cell(S), _col(S, AREA_COLS)
        c[key] = [_row(S, s, nnz, acol, TX) for s in S.samples]
    return c[key]


METRICS = [("median_tx", "median transcripts / cell", 1), ("median_genes", "median genes / cell", 1), ("density", "cells per mm²", 1),
           ("low_frac", f"cells with < {LOW_TX} transcripts", -1), ("median_area", "median cell area (µm²)", 1),
           ("neg_ratio", "background vs genes (%)", -1), ("assigned_frac", "assigned to cells (%)", 1), ("q20_frac", "at Q20 or better (%)", 1)]


def _robust_z(vals, v):
    vals = np.asarray([x for x in vals if x == x], float)
    med = np.median(vals)
    sd = 1.4826 * np.median(np.abs(vals - med))
    return (v - med) / sd if sd > 0 else 0.0


def flags(S, sample, TX=None):
    """Cohort flags: this tissue's metric against the other tissues of the project (robust z). Not a pass / fail against any 10x specification."""
    rows = table(S, TX)
    me = next(r for r in rows if r["sample"] == sample)
    out = []
    for key, label, _ in METRICS:
        if key not in me:
            continue
        vals = [r[key] for r in rows if key in r]
        va = np.array(vals); pct = float((np.mean(va < me[key]) + 0.5 * np.mean(va == me[key])) * 100)   # mid-rank, so ties read 50
        z = _robust_z(vals, me[key]) if len(vals) >= FLAG_MIN_TISSUES else None
        state = "n/a" if z is None else ("high" if z > FLAG_Z else "low" if z < -FLAG_Z else "ok")
        out.append({"key": key, "label": label, "value": me[key], "percentile": pct, "z": z, "state": state})
    return out


def metrics(S, sample, TX=None):
    a, b = S.slices[sample]
    row = next(r for r in table(S, TX) if r["sample"] == sample)
    out = dict(row)
    out["flags"] = flags(S, sample, TX)
    return out


# ------------------------------------------------------------------------------------------------------------------------------ drawing
INK, MUTED, ACC, HI, GRID = "#222", "#777", "#2b5fb8", "#d9534f", "#e3e6ec"
VIRIDIS = np.array([[68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142], [31, 158, 137], [53, 183, 121], [110, 206, 88], [181, 222, 43], [253, 231, 37]], float)


def _lut(n=256):
    t = np.linspace(0, 1, len(VIRIDIS))
    u = np.linspace(0, 1, n)
    return np.stack([np.interp(u, t, VIRIDIS[:, k]) for k in range(3)], 1).astype(np.uint8)


def _e(s):
    return html.escape(str(s))


def _num(x, d=3):
    if x is None or x != x:
        return "–"
    if abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.{d}g}"


def _ticks(lo, hi, n=4):
    if hi <= lo:
        return [lo]
    step = 10 ** math.floor(math.log10((hi - lo) / n))
    for m in (1, 2, 5, 10):
        if (hi - lo) / (step * m) <= n:
            step *= m
            break
    t = math.ceil(lo / step) * step
    out = []
    while t <= hi + 1e-9:
        out.append(round(t, 10)); t += step
    return out


class Chart:
    """A small SVG chart with its own pixel box; sx / sy map data to pixels."""
    def __init__(self, w, h, title, xlabel="", ylabel="", ml=34, mr=8, mt=16, mb=26):
        self.w, self.h, self.ml, self.mr, self.mt, self.mb = w, h, ml, mr, mt, mb
        self.parts = [f'<text x="{ml}" y="10" font-size="9.5" font-weight="600" fill="{INK}">{_e(title)}</text>']
        if xlabel:
            self.parts.append(f'<text x="{(ml + w - mr) / 2}" y="{h - 2}" font-size="8.5" fill="{MUTED}" text-anchor="middle">{_e(xlabel)}</text>')
        if ylabel:
            self.parts.append(f'<text transform="translate(8 {(mt + h - mb) / 2}) rotate(-90)" font-size="8.5" fill="{MUTED}" text-anchor="middle">{_e(ylabel)}</text>')

    def scale(self, xlo, xhi, ylo, yhi):
        self.xlo, self.xhi, self.ylo, self.yhi = xlo, xhi, ylo, yhi
        self.sx = lambda v: self.ml + (v - xlo) / ((xhi - xlo) or 1) * (self.w - self.ml - self.mr)
        self.sy = lambda v: self.h - self.mb - (v - ylo) / ((yhi - ylo) or 1) * (self.h - self.mb - self.mt)

    def axes(self, xt=None, yt=None, xfmt=lambda v: _num(v), yfmt=lambda v: _num(v)):
        x0, x1, y0, y1 = self.ml, self.w - self.mr, self.h - self.mb, self.mt
        self.parts.append(f'<path d="M{x0} {y1}V{y0}H{x1}" fill="none" stroke="{MUTED}" stroke-width=".6"/>')
        for v in (xt if xt is not None else _ticks(self.xlo, self.xhi)):
            self.parts.append(f'<path d="M{self.sx(v):.1f} {y0}v3" stroke="{MUTED}" stroke-width=".6"/><text x="{self.sx(v):.1f}" y="{y0 + 11}" font-size="8" fill="{MUTED}" text-anchor="middle">{_e(xfmt(v))}</text>')
        for v in (yt if yt is not None else _ticks(self.ylo, self.yhi)):
            self.parts.append(f'<path d="M{x0 - 3} {self.sy(v):.1f}h3M{x0} {self.sy(v):.1f}H{x1}" stroke="{GRID}" stroke-width=".5"/><text x="{x0 - 5}" y="{self.sy(v) + 3:.1f}" font-size="8" fill="{MUTED}" text-anchor="end">{_e(yfmt(v))}</text>')

    def svg(self):
        return f'<svg viewBox="0 0 {self.w} {self.h}" xmlns="http://www.w3.org/2000/svg" font-family="Helvetica,Arial,sans-serif">{"".join(self.parts)}</svg>'


def _pow_fmt(v):
    return f"{10 ** v:,.0f}" if 10 ** v >= 1 else f"{10 ** v:.2g}"


def hist(values, title, xlabel, log=False, bins=40, mark=None, w=230, h=125):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if log:
        v = v[v > 0]
    c = Chart(w, h, title, xlabel, "cells")
    if not len(v):
        return c.svg()
    t = np.log10(v) if log else v
    lo, hi = float(np.percentile(t, 0.2)), float(np.percentile(t, 99.8))
    if hi <= lo:
        hi = lo + 1
    cnt, edges = np.histogram(np.clip(t, lo, hi), bins=bins, range=(lo, hi))
    c.scale(lo, hi, 0, cnt.max() * 1.05)
    ticks = [k for k in range(math.ceil(lo), math.floor(hi) + 1)] if log else None
    if log and len(ticks) < 2:
        ticks = None
    c.axes(xt=ticks, yt=[0, cnt.max()], xfmt=_pow_fmt if log else (lambda x: _num(x)), yfmt=lambda y: _num(y))
    for k in range(bins):
        x0, x1, y = c.sx(edges[k]), c.sx(edges[k + 1]), c.sy(cnt[k])
        c.parts.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(.5, x1 - x0 - .3):.1f}" height="{c.sy(0) - y:.1f}" fill="{ACC}"/>')
    if mark is not None and (not log or mark > 0):
        m = math.log10(mark) if log else mark
        m = min(max(m, lo), hi)
        c.parts.append(f'<path d="M{c.sx(m):.1f} {c.mt}V{c.sy(0):.1f}" stroke="{HI}" stroke-width="1" stroke-dasharray="3 2"/>'
                       f'<text x="{c.w - c.mr}" y="{c.mt + 8}" font-size="8" fill="{HI}" text-anchor="end">median {_num(mark)}</text>')
    return c.svg()


def density2d(x, y, title, xlabel, ylabel, logx=False, logy=False, nb=36, w=230, h=150):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & ((x > 0) if logx else True) & ((y > 0) if logy else True)
    c = Chart(w, h, title, xlabel, ylabel)
    if ok.sum() < 10:
        return c.svg()
    x, y = x[ok], y[ok]
    if logx: x = np.log10(x)
    if logy: y = np.log10(y)
    xr = (np.percentile(x, 0.5), np.percentile(x, 99.5)); yr = (np.percentile(y, 0.5), np.percentile(y, 99.5))
    if xr[1] <= xr[0] or yr[1] <= yr[0]:
        return c.svg()
    H, xe, ye = np.histogram2d(x, y, bins=nb, range=[xr, yr])
    c.scale(xr[0], xr[1], yr[0], yr[1])
    ticks = lambda r, lg: ([k for k in range(math.ceil(r[0]), math.floor(r[1]) + 1)] if lg and math.floor(r[1]) - math.ceil(r[0]) >= 1 else None)
    c.axes(xt=ticks(xr, logx), yt=ticks(yr, logy), xfmt=_pow_fmt if logx else (lambda v: _num(v)), yfmt=_pow_fmt if logy else (lambda v: _num(v)))
    L = _lut(64)
    mx = np.log1p(H.max())
    for i in range(nb):
        for j in range(nb):
            if H[i, j] > 0:
                col = L[int(np.log1p(H[i, j]) / mx * 63)]
                c.parts.append(f'<rect x="{c.sx(xe[i]):.1f}" y="{c.sy(ye[j + 1]):.1f}" width="{c.sx(xe[i + 1]) - c.sx(xe[i]) + .3:.1f}" height="{c.sy(ye[j]) - c.sy(ye[j + 1]) + .3:.1f}" fill="rgb({col[0]},{col[1]},{col[2]})"/>')
    return c.svg()


def hbars(items, title, w=230, h=150, unit="%"):
    c = Chart(w, h, title, ml=100, mb=8, mt=16)
    if not items:
        return c.svg()
    mx = max(v for _, v in items) or 1
    rh = (h - 24) / len(items)
    for k, (name, v) in enumerate(items):
        y = 18 + k * rh
        c.parts.append(f'<text x="{c.ml - 4}" y="{y + rh * .62:.1f}" font-size="8" fill="{INK}" text-anchor="end">{_e(name if len(name) < 20 else name[:19] + "…")}</text>'
                       f'<rect x="{c.ml}" y="{y + 1:.1f}" width="{v / mx * (w - c.ml - 36):.1f}" height="{max(2, rh - 3):.1f}" fill="{ACC}"/>'
                       f'<text x="{c.ml + v / mx * (w - c.ml - 36) + 3:.1f}" y="{y + rh * .62:.1f}" font-size="8" fill="{MUTED}">{v:.1f}{unit}</text>')
    return c.svg()


def qv_plot(hist41, title, w=230, h=125):
    c = Chart(w, h, title, "quality value (QV)", "transcripts")
    h41 = np.asarray(hist41, float)
    if h41.sum() == 0:
        return c.svg()
    c.scale(0, 41, 0, h41.max() * 1.05)
    c.axes(xt=[0, 10, 20, 30, 40], yt=[0, h41.max()])
    for q in range(41):
        c.parts.append(f'<rect x="{c.sx(q) + .2:.1f}" y="{c.sy(h41[q]):.1f}" width="{c.sx(1) - c.sx(0) - .4:.1f}" height="{c.sy(0) - c.sy(h41[q]):.1f}" fill="{ACC if q >= 20 else MUTED}"/>')
    c.parts.append(f'<path d="M{c.sx(20):.1f} {c.mt}V{c.sy(0):.1f}" stroke="{HI}" stroke-width="1" stroke-dasharray="3 2"/>'
                   f'<text x="{c.sx(20) + 3:.1f}" y="{c.mt + 8}" font-size="8" fill="{HI}">Q20</text>')
    return c.svg()


def rank_plot(totals, title, w=230, h=150):
    t = np.sort(np.asarray(totals, float))[::-1]
    pos = t[t > 0]
    c = Chart(w, h, title, "gene rank", "total counts", ml=44)
    if len(pos) < 2:
        return c.svg()
    ly = np.log10(pos)
    c.scale(0, len(t), float(ly.min()), float(ly.max()))
    c.axes(xt=_ticks(0, len(t), 3), yt=[k for k in range(math.ceil(ly.min()), math.floor(ly.max()) + 1)] or None, yfmt=_pow_fmt)
    pts = " ".join(f"{c.sx(i):.1f},{c.sy(v):.1f}" for i, v in enumerate(ly))
    c.parts.append(f'<polyline points="{pts}" fill="none" stroke="{ACC}" stroke-width="1.2"/>')
    c.parts.append(f'<text x="{c.w - c.mr}" y="{c.mt + 9}" font-size="8" fill="{MUTED}" text-anchor="end">{int((t == 0).sum())} of {len(t)} genes have no counts</text>')
    return c.svg()


def gene_scatter(mean_pc, det_pc, title, w=230, h=150):
    c = Chart(w, h, title, "mean counts per cell", "% cells detecting", ml=40)
    ok = mean_pc > 0
    if ok.sum() < 3:
        return c.svg()
    lx = np.log10(mean_pc[ok]); y = det_pc[ok]
    c.scale(float(lx.min()), float(lx.max()), 0, 100)
    c.axes(xt=[k for k in range(math.ceil(lx.min()), math.floor(lx.max()) + 1)] or None, yt=[0, 50, 100], xfmt=_pow_fmt)
    c.parts.append("".join(f'<circle cx="{c.sx(a):.1f}" cy="{c.sy(b):.1f}" r="1.1" fill="{ACC}" fill-opacity=".45"/>' for a, b in zip(lx, y)))
    return c.svg()


def gene_qv_plot(totals, mean_qv, genes, title, w=230, h=150):
    """One dot per gene: total counts in the tissue (log y) against the mean QV of the gene's transcripts (x); dashed red line at QV 20."""
    c = Chart(w, h, title, "average QV of the gene's transcripts", "total counts", ml=44)
    t = np.asarray(totals, float); q = np.asarray(mean_qv, float)
    ok = (t > 0) & np.isfinite(q)
    if ok.sum() < 2:
        return c.svg()
    ly = np.log10(t[ok])
    ylo, yhi = min(float(ly.min()), 2.0) - 0.2, max(float(ly.max()), 2.0)         # the axis always reaches the 100-count line, with a little room below it
    c.scale(0, 41, ylo, yhi)
    c.axes(xt=[0, 10, 20, 30, 40], yt=[k for k in range(math.ceil(ylo), math.floor(yhi) + 1)] or None, yfmt=_pow_fmt)
    gn = np.asarray(genes)[ok]
    c.parts.append("".join(f'<circle cx="{c.sx(a):.1f}" cy="{c.sy(b):.1f}" r="1.2" fill="{ACC}" fill-opacity=".5"><title>{_e(g)}: average QV {a:.1f}, {10 ** b:,.0f} counts</title></circle>' for a, b, g in zip(q[ok], ly, gn)))
    c.parts.append(f'<path d="M{c.sx(20):.1f} {c.mt}V{c.sy(c.ylo):.1f}" stroke="{HI}" stroke-width="1" stroke-dasharray="3 2"/>'
                   f'<text x="{c.sx(20) + 3:.1f}" y="{c.mt + 8}" font-size="8" fill="{HI}">QV 20</text>'
                   f'<path d="M{c.ml} {c.sy(2.0):.1f}H{c.w - c.mr}" stroke="{HI}" stroke-width="1" stroke-dasharray="3 2"/>')    # horizontal: 100 counts
    c.parts.append(f'<text x="{c.ml + 4}" y="{c.mt + 8}" font-size="8" fill="{MUTED}">{int(ok.sum())} of {len(t)} genes</text>')
    return c.svg()


def strip_plot(rows, sample, w=470, h=None):
    keys = [(k, lab) for k, lab, _ in METRICS if all(k in r for r in rows)]
    rh = 18
    h = h or 18 + rh * len(keys)
    parts = [f'<text x="2" y="10" font-size="9.5" font-weight="600" fill="{INK}">This tissue against the other {len(rows) - 1} tissue(s) of the project</text>']
    me = next(r for r in rows if r["sample"] == sample)
    for i, (k, lab) in enumerate(keys):
        v = np.array([r[k] for r in rows], float)
        lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
        x0, x1, y = 150, w - 70, 22 + i * rh + 8
        sx = lambda q: x0 + (q - lo) / ((hi - lo) or 1) * (x1 - x0)
        parts.append(f'<text x="{x0 - 6}" y="{y + 3}" font-size="8.5" fill="{INK}" text-anchor="end">{_e(lab)}</text><path d="M{x0} {y}H{x1}" stroke="{GRID}" stroke-width="1"/>')
        rng = np.random.default_rng(i)
        for q, r in zip(v, rows):
            if r["sample"] != sample and q == q:
                parts.append(f'<circle cx="{sx(q):.1f}" cy="{y + rng.uniform(-5, 5):.1f}" r="2" fill="{MUTED}" fill-opacity=".5"/>')
        parts.append(f'<circle cx="{sx(me[k]):.1f}" cy="{y}" r="3.4" fill="{HI}"/><text x="{x1 + 6}" y="{y + 3}" font-size="8.5" fill="{HI}">{_num(me[k])}</text>')
        parts.append(f'<text x="{x0}" y="{y + 12}" font-size="7" fill="{MUTED}">{_num(lo)}</text><text x="{x1}" y="{y + 12}" font-size="7" fill="{MUTED}" text-anchor="end">{_num(hi)}</text>')
    return f'<svg viewBox="0 0 {w} {h + 6}" xmlns="http://www.w3.org/2000/svg" font-family="Helvetica,Arial,sans-serif">{"".join(parts)}</svg>'


def _png_uri(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _nice_um(span):
    for v in (50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000):
        if v >= span / 6:
            return v
    return 10000


def _tick(v):
    a = abs(v)
    if a == 0:
        return "0"
    if a >= 1000:
        return f"{v:,.0f}"
    if a >= 100:
        return f"{v:.0f}"
    if a >= 10:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    return f"{v:.2g}"


def tissue_map(x, y, v=None, px=760, log_label="transcripts / cell", counts=False, bin_um=None, per_bin=False, tissue=None):
    """Map of cells (y downward, like the Xenium images): colour = mean of log10(v) per pixel (cells) or log count per pixel (counts=True,
    for transcript points). Returns HTML (image + colourbar + scale bar).
    per_bin (counts only): colour = the number of points in the bin (1, 2, 3 ... up to 8, the last label meaning "or more"); with tissue = (x, y) of
    points that mark the tissue, a bin that holds some of them but none of these points is drawn grey (empty tissue), distinct from the background."""
    x0, y0 = float(x.min()), float(y.min())
    x1, y1 = float(x.max()), float(y.max())
    if tissue is not None and len(tissue[0]):
        x0, y0 = min(x0, float(np.min(tissue[0]))), min(y0, float(np.min(tissue[1])))
        x1, y1 = max(x1, float(np.max(tissue[0]))), max(y1, float(np.max(tissue[1])))
    w_um, h_um = (x1 - x0) or 1.0, (y1 - y0) or 1.0
    k = (1.0 / bin_um) if (counts and bin_um) else px / max(w_um, h_um)    # fixed physical bin for the density map, so tissues are comparable
    W, H = max(8, int(w_um * k) + 1), max(8, int(h_um * k) + 1)
    cnt = np.zeros(W * H)
    for c0 in range(0, len(x), 4_000_000):             # chunks: every point is binned, without holding int64 copies of millions of points
        xs, ys = np.asarray(x[c0:c0 + 4_000_000]), np.asarray(y[c0:c0 + 4_000_000])
        flat = np.clip(((ys - y0) * k).astype(np.int64), 0, H - 1) * W + np.clip(((xs - x0) * k).astype(np.int64), 0, W - 1)
        cnt += np.bincount(flat, minlength=W * H)
    tm = None
    if counts and per_bin and tissue is not None:      # tissue bins: those holding any of the tissue marker points
        tm = np.zeros(W * H, bool)
        for c0 in range(0, len(tissue[0]), 4_000_000):
            xs, ys = np.asarray(tissue[0][c0:c0 + 4_000_000]), np.asarray(tissue[1][c0:c0 + 4_000_000])
            tm[np.clip(((ys - y0) * k).astype(np.int64), 0, H - 1) * W + np.clip(((xs - x0) * k).astype(np.int64), 0, W - 1)] = True
    if not counts:
        flat = np.clip(((y - y0) * k).astype(np.int64), 0, H - 1) * W + np.clip(((x - x0) * k).astype(np.int64), 0, W - 1)
    if counts and per_bin:
        hi_i = max(2, int(min(cnt.max(), 8))) if cnt.max() > 0 else 2
        val = cnt; lo, hi = 1.0, float(hi_i)              # 1 .. hi_i points per bin, linear
    elif counts:
        if bin_um:                                       # transcripts per mm²
            cnt = cnt / (bin_um * bin_um / 1e6)
        val = np.log1p(cnt); lo, hi = 0.0, float(np.percentile(val[cnt > 0], 99)) if (cnt > 0).any() else 1.0
    else:
        s = np.bincount(flat, weights=np.log10(np.maximum(v, 1)), minlength=W * H)
        val = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
        lo, hi = (float(np.nanpercentile(val, 1)), float(np.nanpercentile(val, 99))) if np.isfinite(val).any() else (0.0, 1.0)
    t = np.clip((val - lo) / ((hi - lo) or 1), 0, 1)
    img = np.full((H * W, 3), 245, np.uint8)
    if tm is not None:
        img[tm] = (196, 200, 207)                        # tissue bin without any of these points: grey, not the background colour
    on = cnt > 0
    img[on] = _lut()[(t[on] * 255).astype(int)]
    img = img.reshape(H, W, 3)
    if k < 2 and not (counts and bin_um):                                                # few pixels per cell: thicken by one pixel so isolated cells stay visible
        m = on.reshape(H, W)
        for dy, dx in ((0, 1), (1, 0), (1, 1)):
            sh = np.zeros_like(m); sh[dy:, dx:] = m[:H - dy, :W - dx]
            tgt = sh & ~m
            img[tgt] = np.roll(np.roll(img, dy, 0), dx, 1)[tgt]
    bar = _nice_um(w_um)
    lab = f"{bar / 1000:g} mm" if bar >= 1000 else f"{bar} µm"
    g = _lut(64)
    cb = "linear-gradient(90deg," + ",".join(f"rgb({r},{gg},{b})" for r, gg, b in g[::8]) + ")"
    if counts and per_bin:                             # colour = points per bin, linear: one label per whole number
        vals = list(range(1, hi_i + 1))
    elif counts:                                       # colour = log1p(count): the labels are log-spaced too
        vals = [math.expm1(hi * t) for t in (0, .25, .5, .75, 1)]
    else:                                              # colour = mean log10: labels at equal steps of it
        vals = [10 ** (lo + (hi - lo) * t) for t in (0, .25, .5, .75, 1)]
    tk = "".join(f"<span>{_tick(v)}</span>" for v in vals) if not per_bin else "".join(f"<span>{v}{'+' if (v == hi_i and cnt.max() > hi_i) else ''}</span>" for v in vals)
    return (f'<div class="map"><div class="im"><img src="{_png_uri(img)}" alt="map"><div class="sb" style="width:{bar * k / W * 100:.1f}%"><i></i><span>{lab}</span></div></div>'
            f'<div class="cb"><b style="background:{cb}"></b><div class="tk">{tk}</div><em>{_e(log_label)}{" (linear)" + ("; grey = tissue bin with none" if tm is not None else "") if per_bin else " (log scale)"}</em></div></div>')


# ------------------------------------------------------------------------------------------------------------------------------ report
CSS = """@page{size:A4;margin:9mm}*{box-sizing:border-box}body{font:11px/1.35 Helvetica,Arial,sans-serif;color:#222;margin:0;background:#fff}
#pg{width:718px;margin:0 auto;padding:6px 0}h1{font-size:17px;margin:0}.sub{color:#777;font-size:9.5px;margin:2px 0 6px}
.cards{display:grid;grid-template-columns:repeat(6,1fr);gap:6px;margin:6px 0}.card{border:1px solid #dde1e8;border-radius:4px;padding:4px 6px}.card b{display:block;font-size:15px}
.card span{color:#777;font-size:8.5px}.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:8px;margin-top:6px}.box{min-width:0;overflow:hidden}.span2{grid-column:span 2}.span3{grid-column:span 3}
.box{border:1px solid #eceff4;border-radius:4px;padding:3px}svg,img{width:100%;display:block}.map{position:relative}.map img{image-rendering:pixelated;max-height:215px;object-fit:contain;background:#f5f5f5}
.im{position:relative}.sb{position:absolute;left:6px;bottom:6px;color:#222;font-size:9px;text-align:center}.sb i{display:block;height:3px;background:#222}.cb{font-size:8.5px;color:#777;margin-top:3px}
.cb b{display:block;height:7px;max-width:190px}.cb .tk{display:flex;justify-content:space-between;max-width:190px;margin-top:1px}.cb em{display:block;font-style:normal}table.f{border-collapse:collapse;width:100%;font-size:9.5px}table.f td:not(:first-child),table.f th{white-space:nowrap}table.f td,table.f th{padding:1px 5px;border-top:1px solid #eceff4;text-align:left}
.ok{color:#2e7d32}.low,.high{color:#d9534f;font-weight:600}.na{color:#999}.note{color:#777;font-size:8.5px;margin-top:4px}.bignum{font-size:18px;font-weight:600;margin-top:4px}h2{font-size:10px;margin:0 0 2px;color:#222}
#bar{margin:8px auto;width:718px;display:flex;gap:8px}#bar button{font:inherit;padding:4px 10px}@media print{#bar{display:none}}"""


def report_html(S, sample, TX=None, project_name="", marker=None, nav=True):
    a, b = S.slices[sample]
    n = b - a
    x, y = np.asarray(S.x[a:b], float), np.asarray(S.y[a:b], float)
    m = metrics(S, sample, TX)
    lib = np.asarray(S.lib[a:b], float) if len(S.lib) else None
    genes_pc = _nnz_per_cell(S)[a:b]
    acol, ncol = _col(S, AREA_COLS), _col(S, NUC_COLS)
    area = S.df[acol].values[a:b].astype(float) if acol else None
    nuc = S.df[ncol].values[a:b].astype(float) if ncol else None
    # per-gene stats of this tissue: the matrix is genes x cells (CSR, sorted indices), so each gene row is cut to the tissue's column range
    G = len(S.genes)
    tot, det = np.zeros(G), np.zeros(G)
    for g in range(G):
        lo, hi = int(S.indptr[g]), int(S.indptr[g + 1])
        ind = S.indices[lo:hi]
        i0, i1 = int(np.searchsorted(ind, a)), int(np.searchsorted(ind, b))
        det[g] = i1 - i0
        tot[g] = float(np.asarray(S.data[lo + i0:lo + i1], np.float64).sum())
    top = np.argsort(-tot)[:10]
    man = {}
    try:
        man = json.loads((Path(S.cache_dir) / "cache_info.json").read_text()).get("manifest", {})
    except Exception:
        pass
    xi = man.get("xenium", {})
    bits = [project_name or man.get("name", "")]
    if "run" in S.df.columns:
        bits.append(f"run {S.df['run'].values[a]}")
    for k, lab in (("slide_id", "slide"), ("panel_name", "panel"), ("analysis_sw_version", "software")):
        if xi.get(k):
            bits.append(f"{lab} {xi[k]}")
    bits.append(f"{G:,} panel genes")
    unit = man.get("coordinate_unit", "um")
    cards = [("cells", f"{n:,}"), ("area (approx.)", f"{m['area_mm2']:.2f} mm²"), ("cells per mm²", f"{m['density']:,.0f}"),
             ("median transcripts / cell", _num(m.get("median_tx"))), ("median genes / cell", _num(m["median_genes"])),
             ("transcripts per mm²", f"{m['tx_density']:,.0f}") if "tx_density" in m else (f"cells with < {LOW_TX} transcripts", f"{m.get('low_frac', float('nan')) * 100:.1f}%")]
    _tq = tx_qc(S, TX, sample)
    if _tq and _tq.get("tissue_rows") is not None:                       # rows of the run's transcripts.parquet that belong to this tissue's cells (every quality, controls included), not the whole run
        cards.append(("transcripts in this tissue (all qualities)", f"{int(_tq['tissue_rows']):,}"))
    elif lib is not None:                                                # no recorded file total: the panel-gene counts summed over this tissue's cells
        cards.append(("total transcripts (from cells)", f"{int(lib.sum()):,}"))
    cell_map = tissue_map(x, y, lib if lib is not None else genes_pc.astype(float), log_label="transcripts / cell" if lib is not None else "genes / cell")
    box = lambda inner, cls="": f'<div class="box {cls}">{inner}</div>'
    row1 = box(f"<h2>Cells, coloured by transcripts per cell</h2>{cell_map}", "span2")
    row1 += '<div>' + box(hist(lib, "Transcripts per cell", "transcripts (panel genes)", log=True, mark=m.get("median_tx"))) + '<div style="height:8px"></div>' + \
        box(hist(genes_pc, "Genes detected per cell", "genes", mark=m["median_genes"])) + '</div>'
    seg = S.df["segmentation_method"].values[a:b] if "segmentation_method" in S.df.columns else None
    items = []
    if seg is not None:
        vc = {}
        for v in seg:
            vc[str(v)] = vc.get(str(v), 0) + 1
        items = sorted(((k.replace("Segmented by ", "").replace("segmented by ", ""), v / n * 100) for k, v in vc.items()), key=lambda t: -t[1])[:6]
    if nuc is not None:
        items.append(("cells with a nucleus", float(np.mean(np.nan_to_num(nuc) > 0) * 100)))
    row2 = box(hist(area, f"Cell area ({'µm²' if unit == 'um' else unit})", "area", mark=m.get("median_area"))) if area is not None else box("<h2>Cell area</h2><div class='note'>no cell-area column in this project</div>")
    row2 += box(hbars(items, "Segmentation / nuclei")) if items else box("<h2>Segmentation</h2><div class='note'>no segmentation columns in this project</div>")
    row2 += box(density2d(area, lib, "Transcripts vs cell area", "cell area", "transcripts", logx=True, logy=True)) if (area is not None and lib is not None) else box("")
    row3 = '<div>' + box(rank_plot(tot, "Genes ranked by counts")) + '<div style="height:8px"></div>' + box(gene_scatter(tot / n, det / n * 100, "Detection vs expression")) + '</div>'   # the two gene plots stacked in the first column
    nc_box = ""
    if TX is not None and sample in TX.loc:
        _tq2 = tx_qc(S, TX, sample)                                      # per-gene average QV recorded at import (every quality, transcripts assigned to this tissue's cells)
        if _tq2 and _tq2.get("gene_qv") is not None:
            row3 = row3[:-len('</div>')] + '<div style="height:8px"></div>' + box(gene_qv_plot(tot, _tq2["gene_qv"][0], S.genes, "Expression vs quality (QV)")) + '</div>'
        else:
            row3 = row3[:-len('</div>')] + '<div style="height:8px"></div>' + box("<h2>Expression vs quality (QV)</h2><div class='note'>not shown: this tissue was imported before the average QV per gene was recorded. Re-import the run to get this plot.</div>") + '</div>'
        px_, py_ = TX.sample_points(sample, None)                       # all stored transcripts of the tissue
        nt = TX.total(sample)
        q_ = tx_qc(S, TX, sample)
        tot_ = (q_['assigned'] + q_['unassigned_near'] + q_['unassigned_far']) if q_ else 0
        un = (f" {(q_['unassigned_near'] + q_['unassigned_far']) / tot_ * 100:.1f}% of them are not assigned to a cell: {q_['unassigned_near']:,} lie within 25 µm of one"
              + (f" and {q_['unassigned_far']:,} are farther (kept, in the tissue of the nearest cell)." if q_['unassigned_far'] else ".") if tot_ else "")
        cp_ = TX.ctrl_points(sample)
        if cp_ is None:
            nc_box = box("<h2>Negative-control density</h2><div class='note'>not shown: this tissue was imported before negative-control positions were stored. Re-import the run to get this map.</div>")
        else:
            ncn = len(cp_[0]); area_ = _area_mm2(S.x[a:b], S.y[a:b])
            cc_ = TX.ctrl_counts(sample); brk_ = ", ".join(f"{k.replace('_', ' ')} {v:,}" for k, v in cc_.items() if v) if cc_ else "classes not recorded"
            nc_box = box(f"<h2>Negative-control density</h2>{tissue_map(cp_[0], cp_[1], counts=True, bin_um=TX_BIN_UM, log_label=f'negative-control transcripts per {TX_BIN_UM:g} µm bin', per_bin=True, tissue=(px_, py_)) if ncn else ''}<div class='note'>{ncn:,} negative-control transcripts (QV ≥ 20; {brk_}) = {ncn / area_ if area_ else float('nan'):,.1f} per mm² of tissue; they should be rare and spread evenly. The pale background is outside the tissue.</div>")
        row3 += box(f"<h2>Transcript density</h2>{tissue_map(px_, py_, counts=True, bin_um=TX_BIN_UM, log_label=f'transcripts / mm² ({TX_BIN_UM:g} µm bins)')}<div class='note'>all {len(px_):,} stored transcripts used (gene probes, QV ≥ 20, control probes excluded).{un}{'' if len(px_) == nt else f' — WARNING: only {len(px_):,} of {nt:,} could be read'} Colour is log-scaled, top = 99th percentile of occupied bins; same bin size for every tissue. Not stored: QV &lt; 20 and control probes. Unassigned transcripts are kept whatever their distance from a cell (projects imported earlier dropped those more than 25 µm away).</div>") + nc_box   # the negative-control map beside it
    else:
        ctc = S.map.get("cell_type")
        ct = S.df[ctc].values[a:b] if ctc in S.df.columns else None
        if ct is not None:
            vc = {}
            for v in ct:
                k = str(v).split(": ", 1)[-1]
                vc[k] = vc.get(k, 0) + 1
            row3 += box(hbars(sorted(((k, v / n * 100) for k, v in vc.items()), key=lambda t: -t[1])[:8], "Largest 10x clusters"))
        row3 += box("<h2>Transcript density</h2><div class='note'>not shown: " + ("this project has no stored transcripts" if TX is None else "no transcripts are stored for this tissue (its run was imported without them)") + "</div>")
    # --- recorded at import: negative controls, transcript quality, assignment
    ct, tq = controls(S, sample, lib), tx_qc(S, TX, sample)
    row_b = ""
    if ct:
        trs = "".join(f"<tr><td>{CTRL_LABEL[k]}</td><td>{v['n_features'] if v['n_features'] is not None else '–'}</td><td>{_num(v['mean'])}</td><td>{v['pct_cells']:.1f}</td></tr>" for k, v in ct["types"].items())
        dist = ""
        if ct["neg_per_cell"] is not None:
            npc = ct["neg_per_cell"]
            dist = hbars([(lab, float(m_.mean() * 100)) for lab, m_ in (("0 per cell", npc == 0), ("1 per cell", npc == 1), ("2 per cell", npc == 2), ("3 or more", npc >= 3))], "Neg. controls per cell", h=78)
        row_b += box(f"<h2>Negative controls</h2><table class='f'><tr><th>feature type</th><th>n</th><th>mean/cell</th><th>% ≥ 1</th></tr>{trs}</table>{dist}"
                     + (f"<div class='note'>Background level: {ct['neg_ratio']:.2f}% of the mean gene (negative-control counts per control feature relative to counts per gene).</div>" if ct["neg_ratio"] is not None else ""))
    if tq and tq["qv_hist"]:
        row_b += box(qv_plot(tq["qv_hist"], "Transcript quality (QV)") + (f"<div class='note'>{tq['q20_frac']:.1f}% at Q20 or better (panel-gene transcripts assigned to cells of this tissue).</div>" if tq["q20_frac"] is not None else ""))
    if tq:
        af = f"{tq['assigned_frac']:.1f}%" if tq["assigned_frac"] is not None else "–"
        rf = f"{tq['run_assigned_frac']:.1f}%" if tq["run_assigned_frac"] is not None else "–"
        row_b += box(f"<h2>Transcript assignment (Q20+, panel genes)</h2><div class='bignum'>{af}</div><div class='note'>assigned to a cell, among transcripts assigned to a cell of this tissue or "
                     f"lying near it ({tq['assigned']:,} assigned, {tq['unassigned_near']:,} unassigned within 25 µm" + (f", {tq['unassigned_far']:,} unassigned farther away" if tq['unassigned_far'] else "") + f").</div><div class='bignum' style='font-size:12px'>{rf}</div>"
                     f"<div class='note'>of all such transcripts of the whole run (including those outside every tissue).</div>")
    missing = [lab for lab, ok in (("negative controls", bool(ct)), ("transcript quality and assignment", bool(tq))) if not ok]
    rows = table(S, TX)
    ctx = box(strip_plot(rows, sample), "span3") if len(rows) > 1 else box("<h2>Context</h2><div class='note'>only one tissue in this project, nothing to compare against</div>", "span3")
    fl = "".join(f"<tr><td>{_e(f['label'])}</td><td>{_num(f['value'])}</td><td>{f"{f['percentile']:.0f}" if len(rows) > 1 else "–"}</td><td class='{f['state']}'>{f['state']}</td></tr>" for f in m["flags"])
    bc = _boundary_counts(S, a, b)
    integ = [f"{(bc[0] + bc[1]):,} of {n:,} cells have no or a degenerate (&lt; 3 vertices) boundary ({(bc[0] + bc[1]) / n * 100:.2f}%)"] if bc else []
    if marker and marker.get("min_counts_per_cell") is not None:
        integ.append(f"cells with fewer than {marker['min_counts_per_cell']} panel-gene transcripts were left out at import")
    flags_box = box(f"<h2>Flags against the other tissues</h2><table class='f'><tr><th>metric</th><th>value</th><th>percentile</th><th>flag</th></tr>{fl}</table>"
                    f"<div class='note'>{'Flag = robust z-score beyond ±' + str(FLAG_Z) + ' versus the project&#39;s tissues' if len(rows) >= FLAG_MIN_TISSUES else 'Flags need at least ' + str(FLAG_MIN_TISSUES) + ' tissues in the project.'} "
                    f"This is a comparison within this project, not a check against a 10x specification.</div>"
                    + ("<div class='note'>" + "<br>".join(integ) + "</div>" if integ else ""), "span2")
    topg = ", ".join(f"{_e(S.genes[g])} ({tot[g]:,.0f})" for g in top)
    foot = ((f"Not recorded in this project: {' and '.join(missing)} (re-import the run{' with the transcripts' if 'transcript quality and assignment' in missing else ''} to add them). ") if missing else "") + (f"Top genes by counts: {topg}. Area is estimated from {BIN_UM:g} µm grid squares that hold at least one cell. Transcripts per cell are panel-gene counts of the cells kept at import. "
            f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by WebbScope.")
    nav_html = NAV if nav else ""
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>QC {_e(sample)}</title><style>{CSS}</style></head><body>{nav_html}<div id="pg">'
            f'<h1>QC report: {_e(sample)}</h1><div class="sub">{_e(" · ".join(str(b) for b in bits if b))}</div>'
            f'<div class="cards" style="grid-template-columns:repeat({len(cards)},1fr)">{"".join(f"<div class=card><b>{_e(v)}</b><span>{_e(k)}</span></div>" for k, v in cards)}</div>'
            f'<div class="grid">{row1}{row2}{row_b}{row3}{flags_box}{ctx}</div><div class="note">{foot}</div></div>{FIT}</body></html>')


FIT = """<script>(function(){var p=document.getElementById("pg");function fit(){p.style.zoom=1;var h=p.scrollHeight;p.style.zoom=Math.min(1,1010/h)}
window.addEventListener("beforeprint",fit);window.addEventListener("afterprint",function(){p.style.zoom=1});if(window.matchMedia){var m=window.matchMedia("print");if(m.addEventListener)m.addEventListener("change",function(e){e.matches?fit():(p.style.zoom=1)})}})()</script>"""
NAV = """<div id="bar"><button onclick="window.print()">Print / save as PDF</button><button onclick="location.href='/qc-summary'" title="one plot per metric across all tissues, outliers marked">All-tissue summary</button><button id="ex" title="write one report per tissue plus a summary table and an all-tissue summary page to a folder on the machine running WebbScope">Export all tissues…</button><span id="exs"></span></div>
<script>document.getElementById("ex").onclick=async()=>{let d0="";try{d0=(await(await fetch("/api/qc-folder")).json()).folder}catch(e){}
const f=prompt("Folder to write the reports to (full path, on the machine running WebbScope):",d0);if(!f)return;
const r=await fetch("/api/qc/export",{method:"POST",headers:{"content-type":"application/json"},body:JSON.stringify({folder:f})});const s=document.getElementById("exs");
if(!r.ok){s.textContent=(await r.json()).detail||"failed";return}
const t=setInterval(async()=>{const d=await(await fetch("/api/qc/export/status")).json();s.textContent=d.state==="running"?`writing ${d.done} of ${d.total}…`:d.state==="done"?`done: ${d.total} reports in ${d.folder}`:`error: ${d.error}`;if(d.state!=="running")clearInterval(t)},1000)}</script>"""


METRIC_HELP = {
    "median_tx": ("Median, over the tissue's cells, of the number of panel-gene transcripts counted in a cell (the sum of a cell's counts over the panel's genes; control probes are not included).", "Low: poorly captured or sparse cells."),
    "median_genes": ("Median number of different panel genes detected (count of at least 1) in a cell.", "Low: cells in which few genes were seen."),
    "density": ("Cells per mm² of occupied tissue. The area is that of the 25 µm grid squares which contain at least one cell.", "Low: sparse or damaged tissue, or cells missed by the segmentation."),
    "low_frac": (f"Fraction of the tissue's cells (0 to 1) that have fewer than {LOW_TX} panel-gene transcripts. Such cells carry little information.", "High: many nearly empty cells."),
    "median_area": ("Median cell area (µm²) from the cell segmentation.", "Low: possible over-segmentation (cells cut into pieces). Unusually large cells (merged cells) are not flagged here."),
    "neg_ratio": ("Background level: counts of the negative-control features (probes and codewords) per control feature, as a percentage of the average counts per gene. About 0 is clean.", "High: much of the signal could be background."),
    "assigned_frac": ("Of the quality-20+ panel-gene transcripts that belong to this tissue (assigned to one of its cells, or unassigned and given to the tissue of its nearest cell; projects imported earlier only kept unassigned ones within 25 µm), the percentage that are assigned to a cell.", "Low: many transcripts lie outside any segmented cell."),
    "q20_frac": ("Percentage of the panel-gene transcripts assigned to this tissue's cells that have a quality value (QV) of at least 20 (Q20), counted over all qualities.", "Low: noisier decoding."),
    "nc_density": ("Negative-control transcripts (all five classes; QV of at least 20, assigned to a cell or within 25 µm of one) per mm² of occupied tissue: the background that the gene counts compete with.", "High: much background signal per area (debris, bad hybridisation or decoding)."),
    "tx_density": ("Stored transcripts (panel genes, QV of at least 20, assigned to a cell or to the tissue of the nearest cell at any distance; projects imported earlier only kept those within 25 µm) per mm² of occupied tissue.", "Low: little signal per area (poor capture, or tissue holes inside the occupied area)."),
    "transcripts_total": ("Number of stored transcripts of the tissue (panel genes, QV of at least 20, assigned to a cell or to the tissue of the nearest cell at any distance; projects imported earlier only kept those within 25 µm).", "Low: a small or poorly captured tissue."),
}


def short_names(names):
    """Drop the leading words (split at _ - or a space) that every tissue name shares, always keeping the last two words.
    Returns (short names, removed prefix, removed suffix = ''); the full names come back unchanged if the short ones would not be unique."""
    names = list(names)
    if len(names) < 2:
        return names, "", ""
    toks = [re.split(r"([_\-\s]+)", n) for n in names]
    ws = [t[0::2] for t in toks]
    m = min(len(w) for w in ws)
    L = 0
    while L < m and all(w[L] == ws[0][L] for w in ws):
        L += 1
    L = max(0, min(L, m - 2))
    T = 0                                              # only a shared beginning is removed: a shared ending is usually the tissue number
    out = []
    for t, w in zip(toks, ws):
        a, b = L, len(w) - T
        out.append("".join(w[k] + (t[2 * k + 1] if k < b - 1 else "") for k in range(a, b)))
    if len(set(out)) != len(out):
        return names, "", ""
    t0, w0 = toks[0], ws[0]
    pre = "".join(w0[k] + t0[2 * k + 1] for k in range(L))
    suf = "".join(t0[2 * k - 1] + w0[k] for k in range(len(w0) - T, len(w0)))
    return out, pre, suf


COHORT_EXTRA = [("tx_density", "transcripts per mm²", 1), ("nc_density", "neg. controls per mm²", -1), ("transcripts_total", "transcripts stored (Q20+, panel genes)", 1)]


def _cohort_panel(key, label, direction, names, full, vals, w=560, h=76):
    """One strip plot: every tissue is a dot on a (log when it spans > 1.5 decades) axis, a line marks the 1st percentile on the bad side (lowest 1%
    for 'higher is better' metrics, highest 1% for 'lower is better'); dots beyond the line are red and named."""
    v = np.asarray(vals, float); ok = np.isfinite(v)
    lo, hi = float(v[ok].min()), float(v[ok].max())
    log = lo > 0 and (hi / lo > 30 or key in ("tx_density", "transcripts_total", "nc_density"))     # the transcript counts are always on a log axis
    f = (lambda q: np.log10(q)) if log else (lambda q: q)
    cut = float(np.percentile(v[ok], 1 if direction > 0 else 99))
    out = ok & ((v < cut) if direction > 0 else (v > cut))
    x0, x1, ya = 24, w - 24, 58
    a, b = f(lo), f(hi)
    sx = lambda q: x0 + (f(q) - a) / ((b - a) or 1) * (x1 - x0)
    parts = [f'<text x="{x0}" y="14" font-size="12" font-weight="600" fill="{INK}">{_e(label)}</text>',
             f'<text x="{x1}" y="14" font-size="10" fill="{MUTED}" text-anchor="end">{"log scale · " if log else ""}n = {int(ok.sum())} · line = {"1st" if direction > 0 else "99th"} percentile ({_num(cut)})</text>',
             f'<path d="M{x0} {ya}H{x1}" stroke="{GRID}" stroke-width="1"/>']
    if log:                                                   # decade ticks
        for e in range(int(np.ceil(a)), int(np.floor(b)) + 1):
            xx = x0 + (e - a) / ((b - a) or 1) * (x1 - x0)
            parts.append(f'<path d="M{xx:.1f} {ya - 4}V{ya + 4}" stroke="{MUTED}" stroke-width=".6"/><text x="{xx:.1f}" y="{ya + 24}" font-size="8" fill="{MUTED}" text-anchor="middle">{_num(10.0 ** e)}</text>')
    parts.append(f'<path d="M{sx(cut):.1f} {ya - 14}V{ya + 14}" stroke="{HI}" stroke-width="1.2" stroke-dasharray="3 2"/>')
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    for i in np.flatnonzero(ok & ~out):
        parts.append(f'<circle cx="{sx(v[i]):.1f}" cy="{ya + rng.uniform(-11, 11):.1f}" r="2.6" fill="{MUTED}" fill-opacity=".55"><title>{_e(full[i])}: {_num(v[i])}</title></circle>')
    for n_, i in enumerate(sorted(np.flatnonzero(out), key=lambda j: v[j] if direction > 0 else -v[j])):
        scr = full[i] if len(full[i]) <= 20 else "…" + full[i][-19:]                 # on screen: the full name (front-truncated); printed: without the shared part
        right = sx(v[i]) < (x0 + x1) / 2
        parts.append(f'<circle cx="{sx(v[i]):.1f}" cy="{ya}" r="3.8" fill="{HI}"><title>{_e(full[i])}: {_num(v[i])}</title></circle>'
                     f'<text x="{sx(v[i]) + (7 if right else -7):.1f}" y="{ya - 15 - (n_ % 3) * 9:.1f}" font-size="8.5" fill="{HI}" text-anchor="{"start" if right else "end"}"><tspan class="f">{_e(scr)}</tspan><tspan class="s">{_e(names[i])}</tspan> ({_num(v[i])})</text>')
    parts.append(f'<text x="{x0}" y="{ya + 34}" font-size="8.5" fill="{MUTED}">{_num(lo)}</text><text x="{x1}" y="{ya + 34}" font-size="8.5" fill="{MUTED}" text-anchor="end">{_num(hi)}</text>')
    return f'<svg viewBox="0 0 {w} {h + 22}" xmlns="http://www.w3.org/2000/svg" font-family="Helvetica,Arial,sans-serif">{"".join(parts)}</svg>', sorted(np.flatnonzero(out).tolist(), key=lambda j: v[j] if direction > 0 else -v[j])


def cohort_html(S, TX=None, project_name=""):
    """All-tissue summary: one strip plot per QC metric (every tissue a dot), a line at the 1st percentile on the bad side, outliers named."""
    rows = table(S, TX)
    full_all = [r["sample"] for r in rows]
    short_all, pre, suf = short_names(full_all)
    sh = dict(zip(full_all, short_all))
    nm = lambda n: f'<span class="f">{_e(n)}</span><span class="s">{_e(sh[n])}</span>'      # full name on screen, without the shared part when printed (PDF)
    panels, flagged = [], {}
    for key, label, direction in METRICS + COHORT_EXTRA:
        have = [r for r in rows if key in r and r[key] == r[key]]
        if len(have) < 3:
            continue
        full = [r["sample"] for r in have]
        vals = [r[key] for r in have]
        svg, bad = _cohort_panel(key, label, direction, [sh[n] for n in full], full, vals)
        what, badside = METRIC_HELP.get(key, ("", ""))
        trs = "".join(f'<tr><td title="{_e(full[i])}">{nm(full[i])}</td><td>{_num(vals[i])}</td></tr>' for i in bad)
        tab = (f'<table class="f"><tr><th>tissue beyond the line</th><th>value</th></tr>{trs}</table>' if bad else '<div class="note">no tissue beyond the line</div>')
        panels.append(f'<div class="box">{svg}<div class="note"><b>What it is:</b> {_e(what)} <b>Flagged when:</b> {_e(badside)}</div>{tab}</div>')
        for i in bad:
            flagged.setdefault(full[i], []).append(label)
    skipped = [lab for key, lab, _ in METRICS + COHORT_EXTRA if sum(1 for r in rows if key in r) < 3]
    ft = "".join(f'<tr><td title="{_e(n)}">{nm(n)}</td><td>{len(l)}</td><td>{_e("; ".join(l))}</td></tr>' for n, l in sorted(flagged.items(), key=lambda t: -len(t[1])))
    common = (f'<span class="f">Names are shown in full here; in the printed / PDF version the part shared by all names ({_e(pre)}…) is left out. </span>'
              f'<span class="s">The part shared by all tissue names is left out of the names ({_e(pre)}…). </span>' if pre else "")
    css = CSS + ".s{display:none}@media print{.f{display:none}.s{display:inline}}#pg{width:1180px}#bar{width:1180px}.pgrid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}.box{padding:6px}"
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>QC summary – {_e(project_name or "project")}</title><style>{css}</style></head><body><div id="bar"><button onclick="window.print()" title="in the print dialog choose Save as PDF; the part shared by all tissue names is left out of the printed names">Print / save as PDF</button></div><div id="pg">'
            f'<h1>QC summary across {len(rows)} tissues</h1><div class="sub">{_e(project_name)} · each dot is one tissue (hover for its name). The dashed red line is the 1st percentile on the bad side '
            f'(lowest 1% where higher is better; highest 1% for "cells with &lt; {LOW_TX} transcripts" and "background vs genes"); tissues beyond it are red and named. '
            f'{common}With fewer than 100 tissues the line catches at most the single most extreme tissue. A comparison within this project, not a 10x specification.</div>'
            f'<div class="pgrid">{"".join(panels)}</div>'
            f'<h2 style="margin-top:12px">Tissues beyond the line</h2>'
            + (f'<table class="f"><tr><th>tissue</th><th>metrics</th><th>which</th></tr>{ft}</table>' if ft else '<div class="note">none</div>')
            + (f'<div class="note">Not drawn (fewer than 3 tissues have the value): {_e(", ".join(skipped))}</div>' if skipped else "") + '</div></body></html>')


def summary_csv(S, TX=None):
    rows = table(S, TX)
    keys = ["sample", "n_cells", "area_mm2", "density", "transcripts_total", "tx_density", "nc_density", "median_tx", "median_genes", "low_frac", "median_area", "no_boundary_frac", "neg_ratio", "assigned_frac", "q20_frac"]
    lines = [",".join(keys + [f"flag_{k}" for k, _, _ in METRICS])]
    for r in rows:
        fl = {f["key"]: f["state"] for f in flags(S, r["sample"], TX)}
        lines.append(",".join([f'"{r["sample"]}"'] + [("" if k not in r else f"{r[k]:.6g}") for k in keys[1:]] + [fl.get(k, "") for k, _, _ in METRICS]))
    return "\n".join(lines) + "\n"
