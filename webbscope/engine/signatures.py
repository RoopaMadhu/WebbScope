"""Resolve the curated signature definitions against the panel and describe how well each one holds together in the data."""
import numpy as np
from .signature_defs import SIGNATURES, ALIASES, ALT_SYMBOLS, EXCLUDED_GENES


def _find(g, panel):
    """Panel symbol for candidate `g` (after aliases), trying the older / newer alternative symbol too; None if absent."""
    g2 = ALIASES.get(g, g)
    for c in (g2, ALT_SYMBOLS.get(g2)):
        if c and c in panel:
            return c
    return None


def resolve(genes_on_panel):
    """List of signatures with only the genes present on the panel (order kept, duplicates and EXCLUDED_GENES removed), plus what is missing."""
    panel = set(genes_on_panel)
    out = []
    for s in SIGNATURES:
        cand = list(dict.fromkeys(ALIASES.get(g, g) for g in s["genes"]))
        present, missing = [], []
        for g in cand:
            f = _find(g, panel)
            if g in EXCLUDED_GENES or f in EXCLUDED_GENES:
                continue
            (present if f else missing).append(f or g)
        present = list(dict.fromkeys(present))
        out.append({"id": s["id"], "group": s["group"], "name": s["name"], "desc": s["desc"], "note": s.get("note", ""),
                    "genes": present, "missing": missing, "n_defined": len(present) + len(missing)})
    return out


def _bins(store, size=100.0):
    """Global id of the (sample, size x size µm bin) of every cell, and the cell count per bin."""
    sid = np.zeros(store.n, np.int64)
    for k, smp in enumerate(store.samples):
        a, b = store.slices[smp]
        sid[a:b] = k
    key = sid * 10**9 + np.floor(store.x / size).astype(np.int64) * 10**4 + np.floor(store.y / size).astype(np.int64)
    _, gid = np.unique(key, return_inverse=True)
    return gid, np.bincount(gid)


def coherence(store, sigs, n_cells=40000, seed=0, bin_um=100.0):
    """How well a signature's genes move together, measured two ways (mean pairwise Pearson correlation of log-normalised expression):
      cell   over a random sample of cells AFTER subtracting each cell type's mean (so cell-type differences cannot create correlation);
      bin    over 100 µm tissue bins (mean per bin, bins with >= 5 cells), the scale the smoothed maps work at.
    Also returns each gene's mean correlation with the other genes of its signature, for both. Single-cell values are low for any imaging panel
    because counts per cell are sparse; the bin values are the informative ones."""
    rng = np.random.default_rng(seed)
    rows = np.sort(rng.choice(store.n, n_cells, replace=False))
    tc = store.cats["cell_type"]["codes"][rows].astype(np.int64)
    gid, cnt = _bins(store, bin_um)
    keep = cnt >= 5
    need = sorted({g for s in sigs for g in s["genes"]})
    cell, binned = {}, {}
    for g in need:
        v = store.gene_vector(g, 0, store.n)
        cell[g] = _resid(v[rows], tc)
        binned[g] = (np.bincount(gid, v, minlength=len(cnt)) / np.maximum(cnt, 1))[keep]
    for s in sigs:
        gs = s["genes"]
        if len(gs) < 2:
            s.update(coherence_cell=None, coherence_bin=None, gene_r_cell={g: None for g in gs}, gene_r_bin={g: None for g in gs})
            continue
        for tag, M in (("cell", cell), ("bin", binned)):
            C = np.corrcoef(np.array([M[g] for g in gs]))
            np.fill_diagonal(C, np.nan)
            s["coherence_" + tag] = float(np.nanmean(C))
            s["gene_r_" + tag] = {g: float(np.nanmean(C[i])) for i, g in enumerate(gs)}
    return rows, tc, cell


def _resid(v, codes):
    m = np.bincount(codes, v) / np.maximum(np.bincount(codes), 1)
    return v - m[codes]


def _rank(v):
    o = np.argsort(v, kind="stable")
    r = np.empty(len(v), np.float64)
    r[o] = np.arange(len(v))
    # average ties (most cells have 0 counts for a gene)
    u, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    first = np.cumsum(cnt) - cnt
    return (first + (cnt - 1) / 2.0)[inv]
