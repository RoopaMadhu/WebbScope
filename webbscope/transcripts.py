"""Read access to a project's transcripts (transcripts/ folder written by the Xenium import): memory-mapped, one gene of one tissue = one slice.

Two layouts of transcripts/index.json:
  single : {"samples": [...], "n_genes": G, ...} with tx_x.npy / tx_y.npy / tx_ptr.npy next to it (one Xenium run)
  parts  : {"parts": [{"dir": <folder of a run's transcripts, relative to this folder>, "samples": [<project sample names, in that run's order>],
            "gene_index": null | [<index in that run's gene list of each project gene>]}, ...], "n_genes": G}   (merged batch projects: the points are read
            from each run's own folder, nothing is copied). A sample no part lists has no transcripts."""
import json
from pathlib import Path
import numpy as np


class Transcripts:
    def __init__(self, folder, samples, n_genes):
        d = Path(folder)
        info = json.loads((d / "index.json").read_text())
        if info["n_genes"] != n_genes:
            raise ValueError(f"transcripts were written for {info['n_genes']} genes, the project has {n_genes}")
        self.info, self.G = info, n_genes
        self.parts, self.loc = [], {}
        specs = info["parts"] if "parts" in info else [{"dir": ".", "samples": info["samples"], "gene_index": None}]
        for pi, sp in enumerate(specs):
            pd_ = (d / sp["dir"]).resolve()
            sub = json.loads((pd_ / "index.json").read_text())
            gi = sp.get("gene_index")
            P = {"x": np.load(pd_ / "tx_x.npy", mmap_mode="r"), "y": np.load(pd_ / "tx_y.npy", mmap_mode="r"), "ptr": np.load(pd_ / "tx_ptr.npy"),
                 "G": int(sub["n_genes"]), "gmap": None if gi is None else np.asarray(gi, np.int64), "qc": sub.get("qc")}
            if (pd_ / "ctl_ptr.npy").exists():   # negative-control probes (imports from this version on)
                P["cx"], P["cy"], P["cptr"] = np.load(pd_ / "ctl_x.npy", mmap_mode="r"), np.load(pd_ / "ctl_y.npy", mmap_mode="r"), np.load(pd_ / "ctl_ptr.npy")
                P["ccat"], P["cats"] = (np.load(pd_ / "ctl_cat.npy", mmap_mode="r") if (pd_ / "ctl_cat.npy").exists() else None), (sub.get("ctrl") or {}).get("categories")
            if int(P["ptr"][-1]) != len(P["x"]):
                raise ValueError(f"transcripts index does not match the coordinate files ({pd_})")
            if len(P["ptr"]) != len(sp["samples"]) * P["G"] + 1:
                raise ValueError(f"transcripts of {pd_} do not match the samples listed for them")
            if P["gmap"] is None and P["G"] != n_genes:
                raise ValueError(f"transcripts of {pd_} have {P['G']} genes, the project has {n_genes}")
            self.parts.append(P)
            for li, s in enumerate(sp["samples"]):
                self.loc[s] = (pi, li)
        missing = [s for s in samples if s not in self.loc]
        if missing and "parts" not in info:
            raise ValueError(f"transcripts do not cover the tissues {missing[:3]}")
        self.n_empty = len(missing)

    def _range(self, sample, g):
        l = self.loc.get(sample)
        if l is None:
            return None, 0, 0
        P = self.parts[l[0]]
        k = l[1] * P["G"] + (g if P["gmap"] is None else int(P["gmap"][g]))
        return P, int(P["ptr"][k]), int(P["ptr"][k + 1])

    def qc(self, sample):
        """Import-time QC numbers of the tissue (None for projects imported before they were recorded): assigned / unassigned_near / unassigned_far counts (far = farther than the import distance from every cell, kept since the import started keeping them; 0 before),
        qv_hist (41 bins, QV 0..40, or None without a quality column) and the run-level counts."""
        l = self.loc.get(sample)
        if l is None or not self.parts[l[0]].get("qc"):
            return None
        q, i = self.parts[l[0]]["qc"], l[1]
        out = {"assigned": q["assigned"][i], "unassigned_near": q["unassigned_near"][i], "unassigned_far": q["unassigned_far"][i] if q.get("unassigned_far") else 0, "qv_hist": q["qv_hist"][i] if q.get("qv_hist") else None, "run": q["run"], "gene_qv": None,
               "tissue_rows": q["tissue_rows"][i] if q.get("tissue_rows") else None}
        if q.get("gene_qv_mean"):                                   # mean QV and transcript number per gene of the project's gene list (NaN / 0 = no transcript)
            m_ = np.array([np.nan if v is None else v for v in q["gene_qv_mean"][i]], float)
            n_ = np.asarray(q["gene_qv_n"][i], float)
            g_ = self.parts[l[0]]["gmap"]
            out["gene_qv"] = (m_, n_) if g_ is None else (m_[g_], n_[g_])
        return out

    def ctrl_points(self, sample):
        """x, y of the tissue's negative-control transcripts (all five negative classes; QV >= the import threshold, same tissue rule as the genes), or None when
        the tissue's run was imported before they were recorded."""
        l = self.loc.get(sample)
        if l is None or "cx" not in self.parts[l[0]]:
            return None
        P = self.parts[l[0]]; a, b = int(P["cptr"][l[1]]), int(P["cptr"][l[1] + 1])
        return np.asarray(P["cx"][a:b]), np.asarray(P["cy"][a:b])

    def ctrl_counts(self, sample):
        """{category: transcripts} of the tissue's negative controls (None if the categories were not recorded)."""
        l = self.loc.get(sample)
        if l is None or self.parts[l[0]].get("ccat") is None:
            return None
        P = self.parts[l[0]]; c = np.asarray(P["ccat"][int(P["cptr"][l[1]]):int(P["cptr"][l[1] + 1])])
        return {n: int((c == i).sum()) for i, n in enumerate(P["cats"])}

    def total(self, sample):
        """Transcripts stored for the tissue (all genes): they are one contiguous range, sorted by gene."""
        l = self.loc.get(sample)
        if l is None:
            return 0
        P = self.parts[l[0]]
        return int(P["ptr"][(l[1] + 1) * P["G"]] - P["ptr"][l[1] * P["G"]])

    def sample_points(self, sample, max_points=400000):
        """x, y of the tissue's transcripts (all genes); max_points thins them evenly, None returns all of them."""
        l = self.loc.get(sample)
        if l is None:
            return np.zeros(0, np.float32), np.zeros(0, np.float32)
        P = self.parts[l[0]]
        a, b = int(P["ptr"][l[1] * P["G"]]), int(P["ptr"][(l[1] + 1) * P["G"]])
        if max_points is not None and b - a > max_points:
            i = a + np.linspace(0, b - a - 1, max_points).astype(np.int64)
            return np.asarray(P["x"][i]), np.asarray(P["y"][i])
        return np.asarray(P["x"][a:b]), np.asarray(P["y"][a:b])

    def count(self, sample, gene):
        return self._range(sample, gene)[2] - self._range(sample, gene)[1]

    def query(self, sample, genes, bbox=None, max_points=200000):
        """genes: project gene indices. Returns (float32 array (n, 3): x, y, position of the gene in `genes`; in-view count per gene; total per gene).
        When more than max_points fall in the box they are thinned evenly per gene (the share of each gene is kept)."""
        parts, in_view, total = [], [], []
        sl = []
        for g in genes:
            P, a, b = self._range(sample, g)
            total.append(b - a)
            x, y = (np.asarray(P["x"][a:b]), np.asarray(P["y"][a:b])) if P is not None else (np.zeros(0, np.float32), np.zeros(0, np.float32))
            if bbox is not None:
                m = (x >= bbox[0]) & (x <= bbox[2]) & (y >= bbox[1]) & (y <= bbox[3])
                x, y = x[m], y[m]
            in_view.append(int(len(x)))
            sl.append((x, y))
        n_in = sum(in_view)
        for j, (x, y) in enumerate(sl):
            cap = max(1, int(max_points * in_view[j] / n_in)) if n_in > max_points else len(x)
            if len(x) > cap:
                i = np.linspace(0, len(x) - 1, cap).astype(np.int64)
                x, y = x[i], y[i]
            parts.append(np.column_stack([x, y, np.full(len(x), j, np.float32)]).astype(np.float32))
        out = np.concatenate(parts) if parts else np.zeros((0, 3), np.float32)
        return out, in_view, total
