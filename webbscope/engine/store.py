"""In-memory metadata + memory-mapped sparse expression access."""
import re
import json
from pathlib import Path
import numpy as np, pandas as pd
from .config import load_mapping, CATEGORICAL_ROLES, REGION_ROLES

MISSING = 65535
SEP = " › "  # region path separator


class Store:
    def __init__(self, cache_dir, config_path=None):
        c = Path(cache_dir)
        self.meta = json.loads((c / "cache_meta.json").read_text())
        self.genes = json.loads((c / "genes.json").read_text())
        self.gene_index = {g: i for i, g in enumerate(self.genes)}
        self.norm_scale = float(self.meta.get("norm_scale", 10000.0))
        self.df = pd.read_parquet(c / "cells.parquet")
        self.n = len(self.df)
        # The viewer, its embeddings and every number in it describe QC'ed Xenium cells only (cells_with_regions_...qs). scFFPE cells have no
        # coordinates and are only used as per-sample counts in the comparison window (scffpe_composition.json), never as cells here.
        if "source" in self.df.columns and not (self.df["source"].astype(str) == "xenium").all():
            raise ValueError("cells.parquet contains non-Xenium cells: " + str(self.df["source"].value_counts().to_dict()))
        sg = c / "signatures.json"   # curated gene signatures resolved against the panel (scripts/build_signatures.py)
        self.signatures = json.loads(sg.read_text()) if sg.exists() else {"signatures": [], "isg_candidates": []}
        self.sig_by_id = {x["id"]: x for x in self.signatures["signatures"]}
        scf = c / "scffpe_composition.json"
        self.scffpe = json.loads(scf.read_text()) if scf.exists() else None
        self.indptr = np.load(c / "csr_indptr.npy", mmap_mode="r")
        self.indices = np.load(c / "csr_indices.npy", mmap_mode="r")
        self.data = np.load(c / "csr_data.npy", mmap_mode="r")
        self.poly_off = self.poly_xy = None
        if (c / "poly_offsets.npy").exists():
            self.poly_off = np.load(c / "poly_offsets.npy", mmap_mode="r")
            self.poly_xy = np.load(c / "poly_xy.npy", mmap_mode="r")
            assert len(self.poly_off) == self.n + 1, "polygon cache does not match cells table; rerun build_polygons.py"
        self.umap = {}    # name -> (n, 2) float32 memmap, rows in cache cell order (scripts/build_umap.py)
        for name in ("harmony",):
            f = c / f"umap_{name}.npy"
            if f.exists():
                e = np.load(f, mmap_mode="r")
                if len(e) == self.n:
                    self.umap[name] = e
        self.pca = None   # Harmony-corrected PCs, (n, 20) float32 memmap
        if (c / "pca_harmony.npy").exists():
            pc = np.load(c / "pca_harmony.npy", mmap_mode="r")
            if len(pc) == self.n:
                self.pca = pc
        self.paths = {}   # sample -> pathologist region polygons (scripts/build_paths.py), optional
        if (c / "paths_by_sample.json").exists():
            self.paths = json.loads((c / "paths_by_sample.json").read_text())["samples"]
            # the polygon export lists some outlines twice (same annotation, id and coordinates); keep one copy so they are not drawn, counted or summed twice
            self.paths_dropped = 0
            for p in self.paths.values():
                seen, keep = set(), []
                for f in p["features"]:
                    pr = f["properties"]
                    k = (pr.get("annotation_name"), pr.get("region_id"), bool(pr.get("negative_roa")), json.dumps(f["geometry"]["coordinates"]))
                    if k in seen:
                        self.paths_dropped += 1
                        continue
                    seen.add(k); keep.append(f)
                p["features"] = keep
        self.map = load_mapping(self.df.columns, config_path)
        self.cache_dir = c
        self._build()

    # -- metadata -----------------------------------------------------------------
    def _build(self):
        df, m = self.df, self.map
        self.sample_col = m["sample"]
        s = df[self.sample_col].astype(str).values
        self.samples = list(dict.fromkeys(s))  # sorted by build_cache, so contiguous
        starts = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
        ends = np.r_[starts[1:], self.n]
        self.slices = {s[a]: (int(a), int(b)) for a, b in zip(starts, ends)}
        self.donor = {k: self.donor_of(k, df.iloc[a]) for k, (a, _) in self.slices.items()}   # sample -> blinded donor id
        self.lib = df[m["lib_size"]].values.astype(np.float32) if m["lib_size"] else np.asarray(self.data[:0])
        self.x = df[m["x"]].values.astype(np.float32)
        self.y = df[m["y"]].values.astype(np.float32)

        # Categorical fields -> uint16 codes with stable level order (by global frequency).
        self.cats = {}
        for role in CATEGORICAL_ROLES:
            col = m.get(role)
            if col:
                self._add_cat(role, df[col], label=col)
        # Region levels are independent categorical fields (filters on different levels combine with AND).
        # Parent cell type of every cell state (the type holding most of its cells; purity = that share), for grouping in the UI.
        self.state_parent = None
        if "cell_state" in self.cats and "cell_type" in self.cats:
            sc, tc = self.cats["cell_state"]["codes"].astype(np.int64), self.cats["cell_type"]["codes"].astype(np.int64)
            nst, nty = len(self.cats["cell_state"]["levels"]), len(self.cats["cell_type"]["levels"])
            ok = (sc != MISSING) & (tc != MISSING)
            tab = np.bincount(sc[ok] * nty + tc[ok], minlength=nst * nty).reshape(nst, nty)
            tot = tab.sum(1)
            parent = np.where(tot > 0, tab.argmax(1), -1)
            self.state_parent = {"parent": parent.tolist(), "min_purity": float((tab.max(1)[tot > 0] / tot[tot > 0]).min())}
        # Pathologist regions (scripts/build_paths.py): multi-membership filter field. codes = bit mask per cell (bit i = inside annotation i),
        # so one cell can belong to several nested regions; the last level means "inside none of the outlines".
        pm, pr = self.cache_dir / "path_mask.npy", self.cache_dir / "path_regions.json"
        if pm.exists() and pr.exists():
            mask, reg = np.load(pm), json.loads(pr.read_text())
            if len(mask) == self.n:
                self.cats["path_region"] = {"label": "Pathologist region", "levels": reg["names"], "codes": mask, "n_missing": 0, "multi": True, "colors": reg["colors"]}
        # Glomerulus identity / QC label (the glom_id / glom_quality columns of the cell table); empty for cells outside gloms.
        for col, label in (("glom_quality", "Glom quality"), ("glom_id", "Glom id")):
            if col in df.columns:
                self._add_cat(col, df[col], label=label)
        # Other low-cardinality string columns become optional colour/filter fields.
        used = set(m.values())
        for col in df.columns:
            if col in used or col in self.cats or df[col].dtype.kind in "fiub":
                continue
            nu = df[col].nunique(dropna=True)
            if 1 < nu <= 300 and col != m["cell_id"]:
                self._add_cat("meta:" + col, df[col], label=col)
        self.numeric = {}
        for col in df.columns:
            if df[col].dtype.kind in "fi" and col not in (m["x"], m["y"]):
                self.numeric["meta:" + col] = col

    def _add_cat(self, key, series, label):
        s = series.astype("string")
        levels = s.value_counts().index.tolist()
        lut = {l: i for i, l in enumerate(levels)}
        codes = s.map(lut).fillna(MISSING).astype(np.uint16).values
        self.cats[key] = {"label": label, "levels": levels, "codes": codes,
                          "n_missing": int((codes == MISSING).sum())}

    def donor_of(self, sample, row):
        """Donor id: the donor column of the mapping (found by name or chosen when the project was opened), else the sample id without its trailing _fragment-N."""
        c = self.map.get("donor")
        d = row.get(c) if c else None
        if d is not None and str(d) not in ("", "nan", "None"):
            return str(d)
        return re.sub(r"_fragment[-_]?\d+$", "", str(sample))

    def describe(self):
        m = self.map
        samples = []
        for s in self.samples:
            a, b = self.slices[s]
            row = self.df.iloc[a]
            samples.append({"id": s, "n": b - a, "x": [float(self.x[a:b].min()), float(self.x[a:b].max())],
                            "y": [float(self.y[a:b].min()), float(self.y[a:b].max())],
                            "disease_type": str(row.get("disease_type", "")),
                            "disease_subtype": str(row.get("disease_subtype", "")),
                            "has_paths": s in self.paths, "offset": int(a),
                            "sample_type": str(row.get("sample_type", "")),
                            "donor": self.donor[s]})
        cats = {}
        for k, v in self.cats.items():
            lv = v["levels"]
            cats[k] = {"label": v["label"], "levels": lv, "n_missing": v["n_missing"]}
            if v.get("multi"):
                cats[k].update({"multi": True, "colors": v["colors"]})
            if k == "cell_state" and self.state_parent:
                cats[k].update(self.state_parent)
        return {"mapping": m, "n_cells": self.n, "n_genes": len(self.genes), "samples": samples,
                "cats": cats, "numeric": {k: v for k, v in self.numeric.items()},
                "matrix": self.meta, "norm_scale": self.norm_scale, "has_polygons": self.poly_off is not None, "umaps": list(self.umap), "n_pcs": int(self.pca.shape[1]) if self.pca is not None else 0, "coordinate_unit": "µm (assumed: area equals area_um in source)"}

    def validation(self):
        df, m = self.df, self.map
        rep = {"n_cells": self.n, "n_genes": len(self.genes), "orientation": self.meta["orientation"],
               "duplicate_cell_ids": int(df[m["cell_id"]].duplicated().sum()),
               "duplicate_genes": len(self.genes) - len(set(self.genes)),
               "missing_x": int(df[m["x"]].isna().sum()), "missing_y": int(df[m["y"]].isna().sum()),
               "missing_annotations": {r: self.cats[r]["n_missing"] for r in CATEGORICAL_ROLES if r in self.cats},
               "layers": self.meta["layers"], "mapping": m,
               "cell_sources": df["source"].astype(str).value_counts().to_dict() if "source" in df.columns else None,
               "n_samples": len(self.samples), "embeddings_rows": {k: int(len(v)) for k, v in self.umap.items()} | ({"pca_harmony": int(len(self.pca))} if self.pca is not None else {})}
        # parent/child consistency: a child path must have exactly one parent by construction;
        # report cells with a child label but missing parent label
        probs = {}
        cols = [m.get(r) for r in REGION_ROLES]
        for i in range(1, 3):
            if cols[i] and cols[i - 1]:
                probs[f"{REGION_ROLES[i]}_without_{REGION_ROLES[i-1]}"] = int((df[cols[i]].notna() & df[cols[i-1]].isna()).sum())
        rep["region_hierarchy_issues"] = probs
        return rep

    def polygons(self, lo, hi):
        """(offsets int32 rebased to 0, n+1; vertices float32 x,y pairs) for cells [lo, hi)."""
        off = np.asarray(self.poly_off[lo:hi + 1])
        return (off - off[0]).astype(np.int32), np.asarray(self.poly_xy[off[0]:off[-1]])

    # -- expression ---------------------------------------------------------------
    def gene_vector(self, gene, lo, hi, layer="logcpx", scale=None):
        """Dense expression for one gene over cells [lo, hi). Only that gene's nonzeros are read.
        logcpx = log1p(counts / lib_size * scale); scale defaults to the source pipeline's normalisation value (harmony_normval)."""
        scale = self.norm_scale if scale is None else scale
        g = self.gene_index[gene]
        a, b = int(self.indptr[g]), int(self.indptr[g + 1])
        idx = np.asarray(self.indices[a:b])
        val = np.asarray(self.data[a:b])
        l, r = np.searchsorted(idx, [lo, hi])
        idx, val = idx[l:r], val[l:r]
        out = np.zeros(hi - lo, dtype=np.float32)
        if layer == "counts":
            out[idx - lo] = val
        else:
            out[idx - lo] = np.log1p(val / np.maximum(self.lib[idx], 1.0) * scale)
        return out

    def resolve_genes(self, names):
        lut = getattr(self, "_gene_lower", None)
        if lut is None:
            lut = self._gene_lower = {g.lower(): g for g in self.genes}
        from .signature_defs import EXCLUDED_GENES, ALT_SYMBOLS
        found, missing, excluded = [], [], []
        for n in dict.fromkeys(x.strip() for x in names):
            g = lut.get(n.lower()) or lut.get(ALT_SYMBOLS.get(n.upper(), "").lower())    # also try the older / newer symbol
            if g is None:
                missing.append(n)
            elif g in EXCLUDED_GENES or n.upper() in EXCLUDED_GENES:
                excluded.append(g)
            elif g not in found:
                found.append(g)
        return {"found": found, "missing": missing, "excluded": excluded}

    def search_genes(self, q, limit=30):
        q = q.lower()
        starts = [g for g in self.genes if g.lower().startswith(q)]
        rest = [g for g in self.genes if q in g.lower() and not g.lower().startswith(q)]
        return (starts + rest)[:limit]
