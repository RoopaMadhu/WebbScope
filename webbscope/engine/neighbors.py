"""Neighbourhood composition: for the cells that pass the viewer's filters ("index cells"), count the cell types
(or any categorical field) of their neighbours, per sample. Three neighbour definitions (all on cell centroids, µm):

  radius    every cell within `distance` of the index cell.
  knn       the `k` nearest cells, keeping only those within `distance`.
  delaunay  cells reachable within `steps` edges of the Delaunay triangulation of the sample's cell centres; triangulation
            edges longer than `distance` are dropped first (they only bridge gaps / the tissue boundary).

Counting is per index-cell/neighbour pair (a neighbour shared by two index cells is counted twice); the cell itself is never
its own neighbour. `bg` is the field's distribution over ALL cells of the sample, the reference for enrichment. The last bin
(index = number of levels) holds cells whose field value is missing.
"""
import threading
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree, Delaunay

MISSING = 65535
MAX_DISTANCE, MAX_K, MAX_STEPS = 150.0, 50, 8
METHODS = ("radius", "knn", "delaunay")


def match(S, key, codes, ids):
    """Boolean per cell: the cell's value of `key` is one of `ids`. Multi-membership fields (pathologist regions) store a bit mask per
    cell, so a cell matches when it lies in ANY of the selected regions."""
    if S.cats[key].get("multi"):
        sel = 0
        for i in ids:
            sel |= 1 << int(i)
        return (codes.astype(np.int64) & sel) != 0
    return np.isin(codes, ids)


class Neighbors:
    def __init__(self, store, max_cache=24):
        self.S, self.max_cache = store, max_cache
        self.SM = None      # Smoother; set by the server (needed by expression filters)
        self.lock = threading.Lock()
        self.trees = {}     # sample -> cKDTree over all cells of the sample
        self.graphs = {}    # (sample, distance) -> pruned Delaunay adjacency (CSR, symmetric)
        self.cache = {}

    def _tree(self, sample):
        with self.lock:
            if sample not in self.trees:
                a, b = self.S.slices[sample]
                self.trees[sample] = cKDTree(np.c_[self.S.x[a:b], self.S.y[a:b]].astype(np.float64))
            return self.trees[sample]

    def _delaunay(self, sample, distance):
        key = (sample, distance)
        with self.lock:
            if key in self.graphs:
                return self.graphs[key]
        pts = self._tree(sample).data
        n = len(pts)
        A = sp.csr_matrix((n, n), dtype=np.float32)
        if n >= 4:
            try:
                indptr, indices = Delaunay(pts).vertex_neighbor_vertices
                rows = np.repeat(np.arange(n), np.diff(indptr))
                keep = np.hypot(*(pts[rows] - pts[indices]).T) <= distance
                A = sp.csr_matrix((np.ones(keep.sum(), np.float32), (rows[keep], indices[keep])), shape=(n, n))
            except Exception:       # degenerate geometry (collinear points): no graph
                pass
        with self.lock:
            self.graphs[key] = A
            while len(self.graphs) > 6:
                self.graphs.pop(next(iter(self.graphs)))
        return A

    def neighbour_pairs(self, sample, idx, method, distance, k, steps):
        """(rows, neighbours): for every index-cell/neighbour pair, the position of the index cell in `idx` and the neighbour's
        sample-local cell index."""
        tree = self._tree(sample)
        pts = tree.data[idx]
        if method == "radius":
            e = cKDTree(pts).sparse_distance_matrix(tree, distance, output_type="ndarray")
            keep = idx[e["i"]] != e["j"]
            return e["i"][keep], e["j"][keep]
        if method == "knn":
            d, j = tree.query(pts, k=k + 1, distance_upper_bound=distance)
            ok = np.isfinite(d) & (j != idx[:, None])               # drop self and "no neighbour within distance"
            ok &= np.cumsum(ok, axis=1) <= k                          # at most k per index cell (self may sit anywhere among ties)
            return np.nonzero(ok)[0], j[ok]
        A = self._delaunay(sample, distance)                          # delaunay
        m = len(idx)
        R = sp.csr_matrix((np.ones(m, np.float32), (np.arange(m), idx)), shape=(m, A.shape[0]))
        M = A + sp.identity(A.shape[0], dtype=np.float32, format="csr")
        for _ in range(steps):
            R = R @ M
            R.data[:] = 1.0
        rows = np.repeat(np.arange(m), np.diff(R.indptr))
        keep = R.indices != idx[rows]
        return rows[keep], R.indices[keep]

    def neighbour_ids(self, sample, idx, method, distance, k, steps):
        """Neighbour cell indices (sample-local), one entry per index-cell/neighbour pair."""
        return self.neighbour_pairs(sample, idx, method, distance, k, steps)[1]

    def _norm(self, method, distance, k, steps):
        if method not in METHODS:
            raise KeyError(f"method must be one of {METHODS}")
        return float(min(max(distance, 1.0), MAX_DISTANCE)), int(min(max(k, 1), MAX_K)), int(min(max(steps, 1), MAX_STEPS))

    def expr_mask(self, sample, expr):
        """Boolean per cell of `sample` for gene / signature filters (all must hold). Each spec: genes, op ("gt": score > thr, "le": score <= thr),
        thr, layer, sigma, scope, within - the same score the viewer colours by (mean over genes -> smooth -> z-score). Cells without a value
        (NaN: outside the "compute over" subset) fail both operators."""
        a, b = self.S.slices[sample]
        m = np.ones(b - a, bool)
        for e in expr or ():
            v = self.SM.expression(sample, e["genes"], e["layer"], e["sigma"], e["scope"], e["within"])
            m &= (v > e["thr"]) if e["op"] == "gt" else (v <= e["thr"])
        return m

    def index_mask(self, sample, filters, hidden, expr=()):
        S = self.S
        a, b = S.slices[sample]
        mask = np.ones(b - a, bool)
        for kk, ids in filters.items():
            if kk not in S.cats:
                raise KeyError(f"unknown filter field {kk!r}")
            mask &= match(S, kk, S.cats[kk]["codes"][a:b], ids)
        for kk, ids in hidden.items():
            if kk not in S.cats:
                raise KeyError(f"unknown filter field {kk!r}")
            mask &= ~match(S, kk, S.cats[kk]["codes"][a:b], ids)
        if expr:
            mask &= self.expr_mask(sample, expr)
        return mask

    def neighbourhoods(self, sample, filters, hidden, method, distance, k, steps, limit, order, seed, expr=()):
        """Index cells of one sample with their neighbour lists (for drawing neighbourhoods in space).
        order: random (seeded) | most | fewest neighbours. Returns at most `limit` index cells."""
        S = self.S
        if sample not in S.slices:
            raise KeyError(f"unknown sample {sample!r}")
        distance, k, steps = self._norm(method, distance, k, steps)
        a, b = S.slices[sample]
        idx = np.flatnonzero(self.index_mask(sample, filters, hidden, expr))
        res = {"sample": sample, "n_index": int(len(idx)), "index": [], "ids": [], "neighbors": []}
        if not len(idx):
            return res
        rows, j = self.neighbour_pairs(sample, idx, method, distance, k, steps)
        cnt = np.bincount(rows, minlength=len(idx))
        if order == "most":
            sel = np.argsort(-cnt, kind="stable")[:limit]
        elif order == "fewest":
            sel = np.argsort(cnt, kind="stable")[:limit]
        else:
            sel = np.random.default_rng(seed).permutation(len(idx))[:limit]
        srt = np.argsort(rows, kind="stable")
        j_s, starts = j[srt], np.r_[0, np.cumsum(cnt)]
        cid = S.df[S.map["cell_id"]].values[a:b]
        res["index"] = idx[sel].tolist()
        res["ids"] = [str(cid[i]) for i in idx[sel]]
        res["neighbors"] = [j_s[starts[r]:starts[r + 1]].tolist() for r in sel]
        return res

    def compute(self, field, method, distance, k, steps, filters, hidden, expr=()):
        S = self.S
        if field not in S.cats:
            raise KeyError(f"unknown field {field!r}")
        if method not in METHODS:
            raise KeyError(f"method must be one of {METHODS}")
        distance = float(min(max(distance, 1.0), MAX_DISTANCE))
        k, steps = int(min(max(k, 1), MAX_K)), int(min(max(steps, 1), MAX_STEPS))
        norm = lambda d: {kk: sorted(int(i) for i in v) for kk, v in sorted(d.items())}
        filters, hidden = norm(filters), norm(hidden)
        for kk in list(filters) + list(hidden):
            if kk not in S.cats:
                raise KeyError(f"unknown filter field {kk!r}")
        key = (field, method, distance, k, steps, repr(filters), repr(hidden), repr(expr))
        if key in self.cache:
            return self.cache[key]
        nlev = len(S.cats[field]["levels"])
        fc = S.cats[field]["codes"].astype(np.int64)
        fc = np.where(fc == MISSING, nlev, fc)
        out = {}
        for sample in S.samples:
            a, b = S.slices[sample]
            mask = np.ones(b - a, bool)
            for kk, ids in filters.items():
                mask &= match(S, kk, S.cats[kk]["codes"][a:b], ids)
            for kk, ids in hidden.items():
                mask &= ~match(S, kk, S.cats[kk]["codes"][a:b], ids)
            if expr:
                mask &= self.expr_mask(sample, expr)
            idx = np.flatnonzero(mask)
            codes = fc[a:b]
            rec = {"n_index": int(len(idx)), "n_cells": int(b - a), "bg": np.bincount(codes, minlength=nlev + 1).tolist(),
                   "counts": [0] * (nlev + 1), "n_pairs": 0}
            if len(idx):
                j = self.neighbour_ids(sample, idx, method, distance, k, steps)
                rec["counts"] = np.bincount(codes[j], minlength=nlev + 1).tolist()
                rec["n_pairs"] = int(len(j))
            out[sample] = rec
        res = {"field": field, "levels": list(S.cats[field]["levels"]) + ["(missing)"], "method": method, "distance": distance,
               "k": k, "steps": steps, "samples": out}
        if len(self.cache) >= self.max_cache:
            self.cache.pop(next(iter(self.cache)))
        self.cache[key] = res
        return res

    def per_cell(self, field, method, distance, k, steps, filters, hidden, expr, groups, bins=40, boot=200, seed=1):
        """Spread of the neighbour composition ACROSS index cells (the pooled composition in compute() hides it).
        `groups`: {level: {sample: group name}}; the index cells of all samples in a group are pooled. Per group and neighbour type:
          q       [5, 25, 50, 75, 95]th percentile and mean of the per-cell neighbour fraction (cells with no neighbour are left out),
          hist    `bins` equal bins on [0, 1] of that fraction (cell counts),
          enr     pooled log2(neighbour fraction / background fraction) (null if the type is never a neighbour), lo / hi: 95% interval of
                  a cluster bootstrap that resamples whole index cells (pairs of one cell are not independent), background fixed."""
        S = self.S
        if field not in S.cats:
            raise KeyError(f"unknown field {field!r}")
        distance, k, steps = self._norm(method, distance, k, steps)
        norm = lambda d: {kk: sorted(int(i) for i in v) for kk, v in sorted(d.items())}
        filters, hidden = norm(filters), norm(hidden)
        gkey = repr(sorted((lv, sorted(m.items())) for lv, m in groups.items()))
        key = ("pc", field, method, distance, k, steps, repr(filters), repr(hidden), repr(expr), gkey)
        if key in self.cache:
            return self.cache[key]
        nlev = len(S.cats[field]["levels"])
        L = nlev + 1
        fc = S.cats[field]["codes"].astype(np.int64)
        fc = np.where(fc == MISSING, nlev, fc)
        per = {}                                    # sample -> (C: index cells x L neighbour counts, bg: counts over all cells)
        for sample in S.samples:
            a, b = S.slices[sample]
            idx = np.flatnonzero(self.index_mask(sample, filters, hidden, expr))
            codes = fc[a:b]
            bg = np.bincount(codes, minlength=L)
            if not len(idx):
                continue
            rows, j = self.neighbour_pairs(sample, idx, method, distance, k, steps)
            C = np.bincount(rows * L + codes[j], minlength=len(idx) * L).reshape(len(idx), L)
            per[sample] = (C, bg, b - a)
        rng = np.random.default_rng(seed)
        out = {}
        for lvl, mp in groups.items():
            res = []
            for name in dict.fromkeys(mp.values()):
                ss = [s for s in per if mp.get(s) == name]
                if not ss:
                    continue
                C = np.vstack([per[s][0] for s in ss])
                bgf = sum(per[s][1] for s in ss) / max(1, sum(per[s][2] for s in ss))
                T = C.sum(1)
                has = T > 0
                Ch, Th = C[has], T[has]
                rec = {"name": name, "n_index": int(len(C)), "n_with_neighbours": int(has.sum()), "n_pairs": int(T.sum()),
                       "bg": bgf.tolist(), "q": [], "hist": [], "enr": [], "lo": [], "hi": []}
                if not has.any():
                    res.append(rec)
                    continue
                F = Ch / Th[:, None]
                qs = np.percentile(F, [5, 25, 50, 75, 95], axis=0)
                for l in range(L):
                    rec["q"].append(qs[:, l].tolist() + [float(F[:, l].mean())])
                    rec["hist"].append(np.histogram(F[:, l], bins=bins, range=(0, 1))[0].tolist())
                n = len(Ch)
                pooled = Ch.sum(0) / Th.sum()
                E = np.full((boot, L), np.nan)
                for t in range(boot):
                    w = np.bincount(rng.integers(0, n, n), minlength=n)
                    f = (w @ Ch) / max(1, w @ Th)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        E[t] = np.where((f > 0) & (bgf > 0), np.log2(f / np.where(bgf > 0, bgf, 1)), np.nan)
                for l in range(L):
                    ok = np.isfinite(E[:, l])
                    lo, hi = np.percentile(E[ok, l], [2.5, 97.5]) if ok.sum() >= boot // 2 else (np.nan, np.nan)
                    e = float(np.log2(pooled[l] / bgf[l])) if pooled[l] > 0 and bgf[l] > 0 else None
                    rec["enr"].append(e)
                    rec["lo"].append(None if np.isnan(lo) else float(lo))
                    rec["hi"].append(None if np.isnan(hi) else float(hi))
                res.append(rec)
            out[lvl] = res
        res = {"field": field, "levels": list(S.cats[field]["levels"]) + ["(missing)"], "bins": bins, "groups": out}
        if len(self.cache) >= self.max_cache:
            self.cache.pop(next(iter(self.cache)))
        self.cache[key] = res
        return res
