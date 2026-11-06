"""Signature statistics: change across disease subtypes, and spatial localisation of two signatures.

All tests use DONORS as the unit (cells within a tissue are spatially autocorrelated, so cell-level p-values would be meaningless). A donor's
tissue pieces (fragments) come from one biopsy, so their cells are pooled into one unit; the "sample" field of the rows holds the donor id.
Scores are the mean log-normalised expression of the signature genes, optionally Gaussian-smoothed over neighbouring cells and z-scored
(see smooth.py).
"""
import numpy as np
from scipy import stats

SUBTYPE_ORDER = ["Control", "No LN", "Class I/II", "Class III", "Class V"]
TYPE_ORDER = ["Control", "No LN", "LN"]
ORDERS = {"subtype": SUBTYPE_ORDER, "type": TYPE_ORDER}


def _holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for rank, i in enumerate(o):
        run = max(run, (len(p) - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def _bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if ok.sum():
        q = p[ok]
        o = np.argsort(q)
        r = q[o] * len(q) / (np.arange(len(q)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        adj = np.empty_like(q)
        adj[o] = np.minimum(1.0, r)
        out[ok] = adj
    return out


def donors_of(S):
    """donor id -> its samples (tissue pieces), in sample order."""
    d = {}
    for s in S.samples:
        d.setdefault(S.donor[s], []).append(s)
    return d


def subtype_of(SM, sample):
    t = SM.info[sample]["disease_subtype"]
    return t if t not in ("", "nan", "None") else SM.info[sample]["disease_type"]


def group_of(SM, sample, group_by):
    return subtype_of(SM, sample) if group_by == "subtype" else (SM.info[sample]["disease_type"] or "unknown")


def summary(S, SM, NB, genes, sigma, scope, filters, hidden, layer="logcpx", min_cells=30, group_by="subtype", expr=()):
    """Per-sample score statistics over the cells passing `filters`, aggregated by disease subtype (5 groups) or disease type (3 groups),
    with across-group tests. scope = how the score is z-scored: none | sample | subtype | type | all (see smooth.py)."""
    ORDER = ORDERS[group_by]
    genes = tuple(genes)
    rows = []
    for dn, smp in donors_of(S).items():
        s0 = smp[0]
        xs, gv, n_cells = [], [[] for _ in genes], 0
        for s in smp:                           # pool the donor's tissue pieces
            a, b = S.slices[s]
            m = NB.index_mask(s, filters, hidden, expr)
            xs.append(SM.expression(s, genes, layer, sigma, scope)[m])
            for k, g in enumerate(genes):
                gv[k].append(S.gene_vector(g, a, b, layer)[m])
            n_cells += b - a
        x = np.concatenate(xs)
        r = {"sample": dn, "subtype": subtype_of(SM, s0), "group": group_of(SM, s0, group_by), "disease_type": SM.info[s0]["disease_type"], "n": int(len(x)), "n_cells": int(n_cells)}
        if len(x):
            q = np.percentile(x, [25, 50, 75])
            r.update(mean=float(x.mean()), median=float(q[1]), q25=float(q[0]), q75=float(q[2]), sd=float(x.std()))
        r["gene_means"] = [float(np.concatenate(v).mean()) if len(x) else None for v in gv]
        rows.append(r)
    groups, gdon = {}, {}
    for r in rows:
        if r["n"] >= min_cells:
            groups.setdefault(r["group"], []).append(r["mean"])
            gdon.setdefault(r["group"], set()).add(r["sample"])
    order = [t for t in ORDER if t in groups] + sorted(t for t in groups if t not in ORDER)
    sub = []
    for t in order:
        v = np.array(groups[t])
        sub.append({"subtype": t, "n_samples": int(len(v)), "n_donors": len(gdon[t]), "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else None,
                    "se": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None})
    tests = {"min_cells": min_cells, "kruskal_p": None, "vs_control": [], "spearman_rho": None, "spearman_p": None}
    big = [t for t in order if len(groups[t]) >= 2]
    if len(big) >= 2:
        try:
            tests["kruskal_p"] = float(stats.kruskal(*[groups[t] for t in big]).pvalue)
        except ValueError:
            pass
    if "Control" in groups and len(groups["Control"]) >= 2:
        rest = [t for t in order if t != "Control" and len(groups[t]) >= 2]
        ps = [float(stats.mannwhitneyu(groups[t], groups["Control"], alternative="two-sided").pvalue) for t in rest]
        adj = _holm(ps) if ps else []
        tests["vs_control"] = [{"subtype": t, "p": p, "p_holm": float(h), "diff": float(np.mean(groups[t]) - np.mean(groups["Control"]))} for t, p, h in zip(rest, ps, adj)]
    ordinal = [(ORDER.index(t), x) for t in order if t in ORDER for x in groups[t]]
    if len(ordinal) >= 6 and len({o for o, _ in ordinal}) >= 3:
        rho = stats.spearmanr([o for o, _ in ordinal], [x for _, x in ordinal])
        tests["spearman_rho"], tests["spearman_p"] = float(rho.statistic), float(rho.pvalue)
    return {"genes": list(genes), "sigma": sigma, "scope": scope, "group_by": group_by, "order": ORDER, "samples": rows, "subtypes": sub, "groups": sub, "tests": tests}


def _boot_ci(vals, n=2000, seed=1):
    v = np.asarray(vals, float)
    if len(v) < 3:
        return [None, None]
    rng = np.random.default_rng(seed)
    m = v[rng.integers(0, len(v), (n, len(v)))].mean(1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def _wilcoxon(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) < 6 or np.allclose(v, 0):
        return None
    try:
        return float(stats.wilcoxon(v).pvalue)
    except ValueError:
        return None


def localization(S, SM, genes_a, genes_b, sigma, q=0.10, subtype=None, min_cells=50, min_samples=4):
    """Do two signatures sit in different places? Scores are smoothed and z-scored within each sample, so both have mean 0 / sd 1 per
    sample and their difference is in SD units. Per sample: correlation of the two maps (raw and after removing each cell type's mean),
    overlap of the top-q hotspots against chance. Per region / cell type: contrast d = mean(A) - mean(B) over the cells inside it (and the
    same after removing cell-type means), aggregated over samples with a bootstrap CI and a Wilcoxon signed-rank test (BH-corrected)."""
    ga, gb = tuple(genes_a), tuple(genes_b)
    pr = S.cats.get("path_region")
    reg_names = pr["levels"] if pr else []
    ct_names = S.cats["cell_type"]["levels"]
    don = {d: smp for d, smp in donors_of(S).items() if subtype is None or subtype_of(SM, smp[0]) == subtype}
    samples = list(don)
    per_sample, reg_rows, ct_rows = [], {}, {}
    pooled_a, pooled_b = [], []
    for s, smp in don.items():                  # s = donor id; a donor's tissue pieces are pooled (z-scores are already donor-wide)
        sl = [S.slices[x] for x in smp]
        za = np.concatenate([SM.expression(x, ga, "logcpx", sigma, "sample") for x in smp]).astype(np.float64)
        zb = np.concatenate([SM.expression(x, gb, "logcpx", sigma, "sample") for x in smp]).astype(np.float64)
        ct = np.concatenate([S.cats["cell_type"]["codes"][a:b] for a, b in sl]).astype(np.int64)
        ok = ct != 65535
        cnt = np.bincount(ct[ok], minlength=len(ct_names))
        ma = np.bincount(ct[ok], za[ok], minlength=len(ct_names)) / np.maximum(cnt, 1)
        mb = np.bincount(ct[ok], zb[ok], minlength=len(ct_names)) / np.maximum(cnt, 1)
        ra, rb = np.where(ok, za - ma[np.where(ok, ct, 0)], np.nan), np.where(ok, zb - mb[np.where(ok, ct, 0)], np.nan)
        n = len(za)
        r_all = float(np.corrcoef(za, zb)[0, 1])
        r_adj = float(np.corrcoef(ra[ok], rb[ok])[0, 1])
        ta, tb = za >= np.quantile(za, 1 - q), zb >= np.quantile(zb, 1 - q)
        expected = max(ta.sum(), tb.sum()) / n
        ov = float((ta & tb).sum() / max(1, min(ta.sum(), tb.sum())))
        per_sample.append({"sample": s, "subtype": subtype_of(SM, smp[0]), "n": int(n), "r": r_all, "r_adj": r_adj, "overlap": ov, "overlap_expected": float(expected),
                           "overlap_ratio": float(ov / expected) if expected else None, "jaccard": float((ta & tb).sum() / max(1, (ta | tb).sum()))})
        pooled_a.append(za)
        pooled_b.append(zb)
        # region contrasts (nested regions are analysed independently; the last level is "outside all outlines")
        if pr is not None:
            bits = np.concatenate([pr["codes"][a:b] for a, b in sl]).astype(np.int64)
            for i, name in enumerate(reg_names):
                m = ((bits >> i) & 1).astype(bool)
                if m.sum() >= min_cells:
                    reg_rows.setdefault(name, []).append((float(za[m].mean()), float(zb[m].mean()), float(np.nanmean(ra[m] - rb[m])), int(m.sum()), s))
        for i, name in enumerate(ct_names):
            m = ct == i
            if m.sum() >= min_cells:
                ct_rows.setdefault(name, []).append((float(za[m].mean()), float(zb[m].mean()), float(np.nanmean(ra[m] - rb[m])), int(m.sum()), s))

    def aggregate(d, adjust=True):   # adjust=False for cell types: removing each cell type's mean would make the adjusted contrast zero by construction
        out = []
        for name, v in d.items():
            if len(v) < min_samples:
                continue
            ndon = len({t[4] for t in v})
            arr = np.array([t[:4] for t in v])
            dd = arr[:, 0] - arr[:, 1]
            out.append({"name": name, "n_samples": int(len(arr)), "n_donors": ndon, "n_cells": int(arr[:, 3].sum()), "mean_a": float(arr[:, 0].mean()), "mean_b": float(arr[:, 1].mean()),
                        "d": float(dd.mean()), "d_ci": _boot_ci(dd), "d_p": _wilcoxon(dd),
                        "d_adj": float(arr[:, 2].mean()) if adjust else None, "d_adj_ci": _boot_ci(arr[:, 2]) if adjust else [None, None], "d_adj_p": _wilcoxon(arr[:, 2]) if adjust else None})
        for key in ("d_p", "d_adj_p"):
            f = _bh([r[key] if r[key] is not None else np.nan for r in out])
            for r, x in zip(out, f):
                r[key.replace("_p", "_fdr")] = None if np.isnan(x) else float(x)
        return sorted(out, key=lambda r: -r["d"])

    allr = np.array([[p["r"], p["r_adj"], p["overlap_ratio"]] for p in per_sample], float)
    pa, pb = np.concatenate(pooled_a), np.concatenate(pooled_b)
    return {"a": list(ga), "b": list(gb), "sigma": sigma, "q": q, "subtype": subtype, "min_cells": min_cells, "min_samples": min_samples, "n_samples": len(samples), "n_donors": len(samples),
            "samples": per_sample, "regions": aggregate(reg_rows), "cell_types": aggregate(ct_rows, adjust=False),
            "overall": {"r_mean": float(np.nanmean(allr[:, 0])), "r_adj_mean": float(np.nanmean(allr[:, 1])), "overlap_ratio_mean": float(np.nanmean(allr[:, 2])),
                        "r_pooled": float(np.corrcoef(pa, pb)[0, 1])}}
