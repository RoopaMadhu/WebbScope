"""On-the-fly spatial smoothing of expression over cell centroids, plus z-scoring within a chosen scope.

Per sample: score = mean over the requested genes of the chosen layer; if sigma > 0 the score is replaced by a
Gaussian-weighted mean over neighbouring cells (weights exp(-d^2 / 2 sigma^2), truncated at 3 sigma, the cell itself
included). Z-scoring then standardises the smoothed score with the mean/sd pooled over every cell of the samples
in the scope (this sample only, same disease subtype, same disease type, or all samples).

`within` restricts both steps to a subset of cells: "tissue" (every cell, the default), "glom" (cells inside a glomerulus: each glomerulus is
smoothed on its own, never mixing with the tissue around it or with another glomerulus, and the z-score mean/SD come from glomerulus cells only)
or "cat:<field>:<level>" (cells of one region level, e.g. the pathologist region Cortex). Cells outside the subset get NaN.
"""
import threading
from collections import OrderedDict
import numpy as np, scipy.sparse as sp
from .signature_defs import EXCLUDED_GENES
from .store import MISSING
from scipy.spatial import cKDTree

MAX_SIGMA = 100.0
WITHIN_TISSUE = "tissue"
SCOPES = ("none", "sample", "subtype", "type", "all", "nocontrol")   # nocontrol: mean/SD from the non-Control samples only, applied to every sample


class Smoother:
    def __init__(self, store, max_w_nnz=30_000_000, max_scores=300):
        self.S, self.max_w_nnz, self.max_scores = store, max_w_nnz, max_scores
        self.lock = threading.Lock()
        self.W = OrderedDict()       # (sample, sigma) -> row-normalised CSR weights
        self.scores = OrderedDict()  # (sample, genes, layer, sigma) -> (vector, n, sum, sumsq)
        self.stats = {}              # (genes, layer, sigma, scope, group) -> (mean, sd)
        df = store.df
        self.info = {s: {c: str(df[c].iloc[store.slices[s][0]]) if c in df.columns else ""
                         for c in ("disease_subtype", "disease_type")} for s in store.samples}

    # -- smoothing weights ----------------------------------------------------------
    def check_within(self, within):
        """Raise ValueError unless `within` is "tissue", "glom" or "cat:<categorical field>:<level index>"."""
        if within == WITHIN_TISSUE:
            return
        if within == "glom":
            if "glom_id" not in self.S.cats:
                raise ValueError("no glomerulus ids in the project (it needs a glom_id column in the cell table)")
            return
        if within.startswith("cat:"):
            key, _, lvl = within[4:].rpartition(":")
            c = self.S.cats.get(key)
            if c is not None and lvl.isdigit() and int(lvl) < len(c["levels"]):
                return
        raise ValueError(f"bad 'within' {within!r}: use tissue, glom or cat:<field>:<level>")

    def _subset(self, sample, within):
        """(bool mask over the sample's cells, group codes or None). Cells of one group are smoothed together only (glomerulus id); mask is None for the whole tissue."""
        if within == WITHIN_TISSUE:
            return None, None
        a, b = self.S.slices[sample]
        if within == "glom":
            g = self.S.cats["glom_id"]["codes"][a:b]
            return g != MISSING, g
        key, _, lvl = within[4:].rpartition(":")
        c, lvl = self.S.cats[key], int(lvl)
        codes = c["codes"][a:b]
        return ((codes.astype(np.int64) >> lvl) & 1).astype(bool) if c.get("multi") else codes == lvl, None

    def _weights(self, sample, sigma, within=WITHIN_TISSUE):
        key = (sample, sigma, within)
        with self.lock:
            if key in self.W:
                self.W.move_to_end(key)
                return self.W[key]
            a, b = self.S.slices[sample]
            mask, group = self._subset(sample, within)
            pts = np.c_[self.S.x[a:b], self.S.y[a:b]].astype(np.float64)
            if mask is not None:
                pts = pts[mask]
            tree = cKDTree(pts)
            e = tree.sparse_distance_matrix(tree, 3 * sigma, output_type="ndarray")
            keep = e["i"] != e["j"]                                   # self weight added explicitly below
            if group is not None:                                     # neighbours must belong to the same glomerulus
                gm = group[mask]
                keep &= gm[e["i"]] == gm[e["j"]]
            n = len(pts)
            i = np.r_[e["i"][keep], np.arange(n)]
            j = np.r_[e["j"][keep], np.arange(n)]
            w = np.r_[np.exp(-e["v"][keep].astype(np.float32) ** 2 / (2 * np.float32(sigma) ** 2)), np.ones(n, np.float32)]
            W = sp.csr_matrix((w, (i, j)), shape=(n, n))
            W = sp.diags((1.0 / np.asarray(W.sum(1)).ravel()).astype(np.float32)) @ W
            self.W[key] = W.astype(np.float32)
            while len(self.W) > 1 and sum(m.nnz for m in self.W.values()) > self.max_w_nnz:
                self.W.popitem(last=False)
            return self.W[key]

    def _score(self, sample, genes, layer, sigma, within=WITHIN_TISSUE):
        key = (sample, genes, layer, sigma, within)
        hit = self.scores.get(key)
        if hit is not None:
            return hit
        a, b = self.S.slices[sample]
        v = np.zeros(b - a, np.float32)
        for g in genes:
            v += self.S.gene_vector(g, a, b, layer)
        v /= len(genes)
        mask, _ = self._subset(sample, within)
        if mask is not None:                    # only the subset takes part; everything else is NaN
            vm = v[mask]
            if sigma > 0 and len(vm):
                vm = np.asarray(self._weights(sample, sigma, within) @ vm, dtype=np.float32)
            v = np.full(b - a, np.nan, np.float32)
            v[mask] = vm
            d = vm.astype(np.float64)
            out = (v, len(vm), float(d.sum()), float((d * d).sum()))
        else:
            if sigma > 0:
                v = np.asarray(self._weights(sample, sigma) @ v, dtype=np.float32)
            d = v.astype(np.float64)
            out = (v, len(v), float(d.sum()), float((d * d).sum()))
        self.scores[key] = out
        while len(self.scores) > self.max_scores:
            self.scores.popitem(last=False)
        return out

    # -- scopes ---------------------------------------------------------------------
    def members(self, sample, scope):
        if scope == "sample":                   # "sample" scope = the whole donor: its tissue pieces (fragments) come from one biopsy and share one mean / SD
            d = self.S.donor[sample]
            return [s for s in self.S.samples if self.S.donor[s] == d]
        if scope == "all":
            return list(self.S.samples)
        if scope == "nocontrol":
            return [s for s in self.S.samples if self.info[s]["disease_type"] != "Control"]
        col = "disease_subtype" if scope == "subtype" else "disease_type"
        return [s for s in self.S.samples if self.info[s][col] == self.info[sample][col]]

    def expression(self, sample, genes, layer="logcpx", sigma=0.0, scope="none", within=WITHIN_TISSUE):
        """Per-cell vector (float32) for `sample`: mean over genes, optionally smoothed, optionally z-scored; NaN outside `within`."""
        genes = tuple(genes)
        if len(genes) > 1:                      # signatures / multi-gene scores never use excluded genes (single genes can still be viewed)
            genes = tuple(g for g in genes if g not in EXCLUDED_GENES)
            if not genes:
                raise ValueError("no genes left after removing excluded genes")
        self.check_within(within)
        sigma = float(min(max(sigma, 0.0), MAX_SIGMA))
        v = self._score(sample, genes, layer, sigma, within)[0]
        if scope == "none":
            return v
        mem = self.members(sample, scope)
        skey = (genes, layer, sigma, scope, tuple(mem), within)
        if skey not in self.stats:
            parts = [self._score(s, genes, layer, sigma, within) for s in mem]
            n = sum(p[1] for p in parts)
            if n == 0:                          # no cell of the subset in any sample of the scope
                self.stats[skey] = (0.0, 1.0)
            else:
                mean = sum(p[2] for p in parts) / n
                sd = float(np.sqrt(max(sum(p[3] for p in parts) / n - mean * mean, 0.0)))
                self.stats[skey] = (mean, sd if sd > 1e-9 else 1.0)
        mean, sd = self.stats[skey]
        return ((v - mean) / sd).astype(np.float32)
