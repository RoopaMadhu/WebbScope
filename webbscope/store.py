"""Read access to an imported project: metadata in memory, sparse counts memory-mapped."""
import json
from pathlib import Path
import numpy as np, pandas as pd

MISSING = 65535


class Store:
    def __init__(self, cdir):
        c = Path(cdir)
        info = json.loads((c / "cache_info.json").read_text())
        self.manifest, self.report = info["manifest"], info["report"]
        self.name = self.manifest.get("name", Path(info["project"]).name)
        self.norm_scale = float(self.manifest.get("norm_scale", 10000.0))
        self.genes = json.loads((c / "genes.json").read_text())
        self.gene_index = {g: i for i, g in enumerate(self.genes)}
        self.df = pd.read_parquet(c / "cells.parquet")
        self.n = len(self.df)
        self.indptr = np.load(c / "csr_indptr.npy", mmap_mode="r")
        self.indices = np.load(c / "csr_indices.npy", mmap_mode="r")
        self.data = np.load(c / "csr_data.npy", mmap_mode="r")
        self.poly_off = np.load(c / "poly_offsets.npy")
        self.poly_xy = np.load(c / "poly_xy.npy")
        self.lib = np.load(c / "lib_size.npy")
        cols = self.manifest.get("columns", {})
        self.id_col, self.xcol, self.ycol = cols.get("cell_id", "cell_id"), cols.get("x", "x"), cols.get("y", "y")
        self.xy = np.ascontiguousarray(self.df[[self.xcol, self.ycol]].values, dtype=np.float32)
        self.cats, self.numeric = {}, {}
        for col in self.manifest.get("categorical", []):
            if col in self.df.columns:
                s = self.df[col].astype("string")
                levels = s.value_counts().index.tolist()
                lut = {l: i for i, l in enumerate(levels)}
                codes = s.map(lut).fillna(MISSING).astype(np.uint16).values
                self.cats[col] = {"levels": levels, "codes": codes, "n_missing": int((codes == MISSING).sum())}
        skip = {self.xcol, self.ycol, self.id_col}
        for col in self.df.columns:
            if col not in skip and self.df[col].dtype.kind in "fi":
                self.numeric[col] = self.df[col].values.astype(np.float32)

    def describe(self):
        return {"name": self.name, "n_cells": self.n, "n_genes": len(self.genes), "norm_scale": self.norm_scale,
                "bounds": [float(self.xy[:, 0].min()), float(self.xy[:, 1].min()), float(self.xy[:, 0].max()), float(self.xy[:, 1].max())],
                "cats": {k: {"levels": v["levels"], "n_missing": v["n_missing"]} for k, v in self.cats.items()},
                "numeric": list(self.numeric), "report": self.report, "unit": self.manifest.get("coordinate_unit", "")}

    def gene_vector(self, gene, layer="logcpx"):
        """Dense expression of one gene over all cells; only that gene's non-zeros are read.
        logcpx = log1p(counts / lib_size * norm_scale)."""
        g = self.gene_index[gene]
        a, b = int(self.indptr[g]), int(self.indptr[g + 1])
        idx, val = np.asarray(self.indices[a:b]), np.asarray(self.data[a:b])
        out = np.zeros(self.n, dtype=np.float32)
        if layer == "counts":
            out[idx] = val
        else:
            out[idx] = np.log1p(val / np.maximum(self.lib[idx], 1.0) * self.norm_scale)
        return out

    def resolve_gene(self, name):
        if name in self.gene_index:
            return name
        low = {g.lower(): g for g in self.genes}
        return low.get(name.lower())

    def search_genes(self, q, limit=30):
        q = q.lower()
        starts = [g for g in self.genes if g.lower().startswith(q)]
        rest = [g for g in self.genes if q in g.lower() and not g.lower().startswith(q)]
        return (starts + rest)[:limit]
