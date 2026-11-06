"""Per-sample distributions of a gene / signature score within a cell type, for violin plots (one violin per sample).

For every (feature, cell type, sample): the number of cells, mean, median, quartiles, fraction of cells with a non-zero score and a density curve
on a grid shared by all violins of the feature, so violins can be compared. Scores are the same as everywhere else (mean log-normalised
expression of the genes; optional Gaussian smoothing and z-scoring, see smooth.py). The density is a Gaussian kernel estimate (Scott's rule,
bounded) scaled to a maximum of 1, i.e. every violin is drawn at the same maximum width and its cell count is reported separately.
"""
import numpy as np

GRID = 64


def _density(x, lo, hi):
    if len(x) < 2:
        return [0.0] * GRID
    h = float(np.clip(1.06 * x.std() * len(x) ** -0.2, (hi - lo) / 80, (hi - lo) / 6))
    step = (hi - lo) / (GRID - 1)
    cnt, _ = np.histogram(np.clip(x, lo, hi), bins=GRID, range=(lo - step / 2, hi + step / 2))
    k = np.exp(-0.5 * (np.arange(-GRID, GRID + 1) * step / h) ** 2)             # kernel centred at index GRID
    d = np.convolve(cnt.astype(float), k, mode="full")[GRID:2 * GRID]             # aligned back onto the GRID bins
    return [round(float(v), 3) for v in d / d.max()] if d.max() > 0 else [0.0] * GRID


def _aggregate(rows, label_of):
    by = {}
    for r in rows:
        by.setdefault(label_of[r["sample"]], []).append(r)
    out = []
    for lab, rs in by.items():
        ok = [r for r in rs if "mean" in r]
        fr = [r["frac_nonzero"] for r in rs if r["frac_nonzero"] is not None]
        o = {"sample": lab, "n_donors": len(rs), "n_donors_used": len(ok), "n": int(sum(r["n"] for r in rs)), "n_used": int(sum(r["n_used"] for r in rs)),
             "frac_nonzero": float(np.mean(fr)) if fr else None}
        if ok:
            m = np.array([r["mean"] for r in ok])
            o.update(mean=float(m.mean()), median=float(np.median(m)), sd=float(m.std(ddof=1)) if len(m) > 1 else None,
                     se=float(m.std(ddof=1) / np.sqrt(len(m))) if len(m) > 1 else None, donor_means=[float(v) for v in m], donors=[r["sample"] for r in ok])
        out.append(o)
    return out


def _natural(name):
    import re
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def violins(S, SM, feats, cell_types, sigma, scope, layer="logcpx", nonzero_only=False, min_cells=10, level="cell_type", unit="sample"):
    """feats: [{"label", "genes"}]; cell_types: level ids of cell_type (or -1 = all cells).
    level "cell_type": one panel per (feature, cell type), one group (the cell type) of per-sample violins.
    level "cell_state": one panel per (feature, cell type) with one group per cell state of that cell type (cell_types must be real types)."""
    ct_codes, ct_names = S.cats["cell_type"]["codes"], S.cats["cell_type"]["levels"]
    cs_codes, cs_names = S.cats["cell_state"]["codes"], S.cats["cell_state"]["levels"]
    parent = (S.state_parent or {}).get("parent", [])
    if level == "cell_state" and (not parent or any(c < 0 for c in cell_types)):
        raise ValueError("cell-state violins need real cell types (not 'all cells') and the cell-state parent mapping")

    def groups_of(ct):
        if level == "cell_type":
            return [(ct, "All cells" if ct < 0 else ct_names[ct])]
        return [(i, cs_names[i]) for i in sorted((i for i, p in enumerate(parent) if p == ct), key=lambda i: _natural(cs_names[i]))]

    def mask(ct, gid, a, b):
        if level == "cell_type":
            return np.ones(b - a, bool) if gid < 0 else ct_codes[a:b] == gid
        return cs_codes[a:b] == gid

    from . import sigstats
    units = sigstats.donors_of(S)             # per-donor values: the cells of a donor's tissue pieces (one biopsy) are pooled
    label_of = {d: (sigstats.subtype_of(SM, m[0]) if unit == "subtype" else SM.info[m[0]]["disease_type"]) for d, m in units.items()}   # used when unit is subtype / type

    def cells(vec, ct, gid, members):
        return np.concatenate([vec[s][mask(ct, gid, *S.slices[s])] for s in members])

    out = []
    for f in feats:
        vec = {s: SM.expression(s, tuple(f["genes"]), layer, sigma, scope) for s in S.samples}
        sel = {}
        for ct in cell_types:
            for gid, _ in groups_of(ct):
                for u, members in units.items():
                    x = cells(vec, ct, gid, members)
                    sel[(ct, gid, u)] = x[x > 0] if nonzero_only else x
        pool = np.concatenate([v for v in sel.values() if len(v)] or [np.zeros(1)])
        lo, hi = (float(np.percentile(pool, 0.5)), float(np.percentile(pool, 99.5))) if len(pool) else (0.0, 1.0)
        if hi - lo < 1e-6:
            lo, hi = lo - 0.5, hi + 0.5
        for ct in cell_types:
            groups = []
            for gid, gname in groups_of(ct):
                rows = []
                for u, members in units.items():
                    full = cells(vec, ct, gid, members)
                    x = sel[(ct, gid, u)]
                    r = {"sample": u, "n": int(len(full)), "frac_nonzero": float((full > 0).mean()) if len(full) else None, "n_used": int(len(x))}
                    if len(x) >= min_cells:
                        q = np.percentile(x, [25, 50, 75])
                        r.update(mean=float(x.mean()), median=float(q[1]), q25=float(q[0]), q75=float(q[2]), density=_density(x, lo, hi))
                    rows.append(r)
                if unit != "sample":      # bar plot data: the mean of the per-donor means of each disease subtype / type (each donor counts once)
                    rows = _aggregate(rows, label_of)
                groups.append({"id": int(gid), "name": gname, "samples": rows})
            out.append({"feature": f["label"], "cell_type": "All cells" if ct < 0 else ct_names[ct], "cell_type_id": int(ct), "level": level, "lo": lo, "hi": hi, "groups": groups})
    return out


def histograms(S, SM, feats, cell_types, sigma, scope, layer="logcpx", nonzero_only=False, bins=60, unit="sample", min_cells=10):
    """Histogram of the per-cell values that the violins are made from, one panel per (feature, cell type). Exact zeros get their own bar.
    unit "sample": all donors pooled in one histogram. unit "subtype" / "type": one histogram per disease subtype / type, each the AVERAGE over its
    donors of the donor's own histogram (share of the donor's cells per bin), so every donor counts once whatever its cell number.
    For a single gene there is also a histogram of its raw counts in the same cells, and a check that the plotted value equals
    log1p(raw count / total counts x scale) recomputed from the counts (only when no smoothing / z-score is applied)."""
    from . import sigstats
    ct_codes, ct_names = S.cats["cell_type"]["codes"], S.cats["cell_type"]["levels"]
    donors = sigstats.donors_of(S)
    label_of = {d: (sigstats.subtype_of(SM, m[0]) if unit == "subtype" else SM.info[m[0]]["disease_type"]) for d, m in donors.items()} if unit != "sample" else {}
    out = []
    for f in feats:
        genes = tuple(f["genes"]); single = len(genes) == 1
        vec = {s: SM.expression(s, genes, layer, sigma, scope) for s in S.samples}
        raw = {s: SM.expression(s, genes, "counts", 0.0, "none") for s in S.samples} if single else None
        for ct in cell_types:
            per, es = {}, []            # donor -> (values, raw counts)
            for d, members in donors.items():
                xs, rs = [], []
                for s in members:
                    a, b = S.slices[s]
                    m = np.ones(b - a, bool) if ct < 0 else ct_codes[a:b] == ct
                    v = vec[s]
                    keep = m & np.isfinite(v)
                    if nonzero_only:
                        keep &= v > 0
                    xs.append(v[keep])
                    if single:
                        c = raw[s][keep]
                        rs.append(c)
                        if layer == "logcpx" and sigma == 0 and scope == "none":
                            es.append(np.abs(v[keep] - np.log1p(c / np.maximum(S.lib[a:b][keep], 1.0) * S.norm_scale)))
                per[d] = (np.concatenate(xs), np.concatenate(rs) if single else None)
            x = np.concatenate([p[0] for p in per.values()])
            n = int(len(x))
            rec = {"feature": f["label"], "cell_type": "All cells" if ct < 0 else ct_names[ct], "cell_type_id": int(ct), "n": n, "layer": layer, "unit": unit}
            if n:
                nz = x[x != 0]
                rec.update(n_zero=int(n - len(nz)), mean=float(x.mean()), median=float(np.median(x)), max=float(x.max()),
                           min_nonzero=float(nz.min()) if len(nz) else None)
                lo = hi = None
                if len(nz):
                    lo, hi = float(nz.min()), float(np.percentile(nz, 99.9))
                    if hi - lo < 1e-9:
                        lo, hi = lo - 0.5, hi + 0.5
                    cnt, edges = np.histogram(np.clip(nz, lo, hi), bins=bins, range=(lo, hi))
                    rec.update(edges=[round(float(e), 4) for e in edges], counts=[int(c) for c in cnt])
                else:
                    rec.update(edges=[], counts=[])
                if single:
                    r = np.concatenate([p[1] for p in per.values()])
                    K = int(min(30, max(1, r.max())))
                    rc = np.bincount(np.minimum(r.astype(int), K), minlength=K + 1)
                    rec["raw"] = {"counts": [int(c) for c in rc], "top": K, "max": float(r.max()), "mean": float(r.mean()), "frac_zero": float((r == 0).mean())}
                    if es:
                        rec["check_max_abs_diff"] = float(np.concatenate(es).max())
                if unit != "sample":
                    def share(xd, rd):
                        z = int((xd == 0).sum()); nzd = xd[xd != 0]
                        h = np.histogram(np.clip(nzd, lo, hi), bins=bins, range=(lo, hi))[0] if lo is not None else np.zeros(0)
                        sh = [z / len(xd)] + [c / len(xd) for c in h]
                        rs = list(np.bincount(np.minimum(rd.astype(int), rec["raw"]["top"]), minlength=rec["raw"]["top"] + 1) / len(xd)) if rd is not None else None
                        return sh, rs, float(xd.mean())
                    by = {}
                    for d, (xd, rd) in per.items():
                        by.setdefault(label_of[d], []).append((d, xd, rd))
                    groups = []
                    for lab, items in by.items():
                        ok = [(d, xd, rd) for d, xd, rd in items if len(xd) >= min_cells]
                        g = {"name": lab, "n_donors": len(items), "n_donors_used": len(ok), "n": int(sum(len(xd) for _, xd, _ in items))}
                        if ok:
                            sh = [share(xd, rd) for _, xd, rd in ok]
                            g["share"] = [round(float(v), 5) for v in np.mean([t[0] for t in sh], axis=0)]
                            if sh[0][1] is not None:
                                g["raw_share"] = [round(float(v), 5) for v in np.mean([t[1] for t in sh], axis=0)]
                            g["mean"] = float(np.mean([t[2] for t in sh]))
                        groups.append(g)
                    rec["groups"] = groups
            out.append(rec)
    return out
