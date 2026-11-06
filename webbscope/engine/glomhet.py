"""ROI (glomerulus) heterogeneity: embeddings of the regions (one point per ROI), colouring, filtering by the viewer's cell filters, and
Spearman correlations of genes / cell types / cell states with the first PCs. Data: <cache>/gloms/ (written by the importer from gloms/ or from a
region_heterogeneity object, roi_object.py)."""
import json
import numpy as np
from scipy.stats import rankdata
from .smooth import EXCLUDED_GENES


def _rank_cols(a):
    return np.apply_along_axis(rankdata, 0, a) if a.ndim == 2 else rankdata(a)


class GlomHet:
    def __init__(self, store, nb, sm=None):
        d = store.cache_dir / "gloms"
        self.ok = (d / "gloms.json").exists()
        self.S, self.NB, self.SM = store, nb, sm
        if not self.ok:
            return
        j = json.loads((d / "gloms.json").read_text())
        z = np.load(d / "gloms.npz")
        self.genes, self.meta, self.cells = j["genes"], j["meta"], [np.asarray(c) for c in j["cells"]]
        self.gi = {g: i for i, g in enumerate(self.genes)}
        self.E, self.emb = z["expr"], {k: z[f] for k, f in (("harmony", "pca_harmony"), ("raw", "pca")) if f in z.files}   # the batch-corrected PCA is optional
        self.umap = z["umap"]
        self.Ct, self.Cs = z["type_counts"], z["state_counts"]
        self.type_levels, self.state_levels = j["type_levels"], j["state_levels"]
        self.parent = (store.state_parent or {}).get("parent", [])
        self.n = len(self.meta)
        self.sample = np.array([m["sample"] for m in self.meta])

    def info(self):
        S = self.S
        cols = [k for k in self.meta[0] if k.startswith("cluster_")]
        f = lambda v: None if v is None else float(v)
        return {"n": self.n, "gloms": [{"id": m["glom_id_final"], "sample": m["sample"], "disease_type": m.get("disease_type"), "disease_subtype": m.get("disease_subtype"),
                                        "activity": f(m.get("Activity")), "chronicity": f(m.get("Chronicity")), "n_file": m.get("N"), "n_cells": m["n_viewer_cells"],
                                        "logUMI": m.get("logUMI"), **{c: m[c] for c in cols}} for m in self.meta],
                "umap": np.round(self.umap, 4).tolist(), "pca": {k: np.round(v[:, :10], 4).tolist() for k, v in self.emb.items()},
                "n_pcs": 10, "cluster_cols": cols, "type_levels": self.type_levels, "state_levels": self.state_levels, "n_genes": len(self.genes)}

    def passing(self, filters, hidden, expr=()):
        """Number of each glom's viewer cells that pass the viewer's filters (a glom is shown when >= 1 passes)."""
        S = self.S
        out = np.zeros(self.n, int)
        masks = {}
        for i, c in enumerate(self.cells):
            s = self.meta[i]["sample"]
            if s not in masks:
                masks[s] = self.NB.index_mask(s, filters, hidden, expr)
            out[i] = int(masks[s][c - S.slices[s][0]].sum())
        return out

    def color(self, kind, genes=None, level=None):
        """One value per glom. kind: gene/signature (mean over genes of the glom-level log-normalised expression), type/state (fraction of the glom's cells)."""
        if kind in ("gene", "signature"):
            genes = list(genes)
            if len(genes) > 1:
                genes = [g for g in genes if g not in EXCLUDED_GENES]
            use = [g for g in genes if g in self.gi]
            if not use:
                raise ValueError(f"none of these genes is in the region_heterogeneity object ({len(self.genes):,} genes)")
            return self.E[[self.gi[g] for g in use]].mean(axis=0), {"genes_used": len(use), "genes_requested": len(genes), "missing": [g for g in genes if g not in self.gi][:20]}
        C, lv = (self.Ct, self.type_levels) if kind == "type" else (self.Cs, self.state_levels)
        tot = self.Ct.sum(1).clip(min=1)
        return C[:, level] / tot, {"label": lv[level]}

    def cell_score(self, genes, layer="logcpx", sigma=0.0, scope="none", within="tissue"):
        """One value per glom: mean over the glom's viewer cells of the cell-level score the main window shows (mean over genes of the layer,
        optionally Gaussian-smoothed, optionally z-scored; `within` = "glom" smooths each glom on its own and z-scores against glom cells only).
        NaN when the glom has no cell inside `within`."""
        out = np.full(self.n, np.nan)
        by_sample = {}
        for i, m in enumerate(self.meta):
            by_sample.setdefault(m["sample"], []).append(i)
        for s, ids in by_sample.items():
            v = self.SM.expression(s, genes, layer, sigma, scope, within)
            a = self.S.slices[s][0]
            for i in ids:
                x = v[self.cells[i] - a]
                x = x[~np.isnan(x)]
                if len(x):
                    out[i] = float(x.mean())
        return out

    def correlations(self, idx, layer, kind, npc=10, ntop=10, sort_pc=None):
        """Spearman correlation of gloms `idx` between PCs 1..npc and genes (top/bottom `ntop` per PC, unique), cell types or cell states.
        sort_pc (0-based, genes only): instead of the genes that stand out in any PC, return the `ntop` most positive and `ntop` most negative genes of that one PC
        (still with the correlations of every PC), ordered by their correlation with it."""
        idx = np.asarray(idx, int)
        if len(idx) < 8:
            raise ValueError("need at least 8 gloms for correlations")
        P = _rank_cols(self.emb[layer][idx, :npc])
        if kind == "gene":
            X = self.E[:, idx].T
        elif kind == "type":
            X = (self.Ct / self.Ct.sum(1, keepdims=True).clip(min=1))[idx]
        else:
            X = (self.Cs / self.Ct.sum(1, keepdims=True).clip(min=1))[idx]
        keep = X.std(axis=0) > 1e-9
        R = _rank_cols(X[:, keep])
        Z = lambda a: (a - a.mean(0)) / a.std(0).clip(min=1e-12)
        r = (Z(R).T @ Z(P)) / len(idx)                      # features x npc
        names = {"gene": self.genes, "type": self.type_levels, "state": self.state_levels}[kind]
        names = [n for n, k in zip(names, keep) if k]
        if kind == "gene" and sort_pc is not None:
            o = np.argsort(r[:, sort_pc])
            sel = list(dict.fromkeys(list(o[-ntop:][::-1]) + list(o[:ntop])))
            sel.sort(key=lambda f: -float(r[f, sort_pc]))
        elif kind == "gene":
            sel = []
            for j in range(npc):
                o = np.argsort(r[:, j]); sel += list(o[-ntop:][::-1]) + list(o[:ntop])
            sel = list(dict.fromkeys(sel))
            # order: by the PC where the gene's |r| is largest, then by r within that PC (positive first), so the heatmap is block-diagonal
            sel.sort(key=lambda f: (int(np.abs(r[f]).argmax()), -float(r[f][np.abs(r[f]).argmax()])))
        elif kind == "state":
            allp = self.parent
            kidx = [i for i, k in enumerate(keep) if k]
            sel = sorted(range(len(names)), key=lambda f: (allp[kidx[f]] if allp else 0, names[f]))
        else:
            sel = list(range(len(names)))
        out = {"rows": [names[f] for f in sel], "r": np.round(r[sel], 3).tolist(), "n_gloms": int(len(idx)), "npc": npc}
        if kind == "state":
            out["parents"] = [self.type_levels[self.parent[kidx[f]]] if self.parent else "" for f in sel]
        return out
