"""region_heterogeneity object (.rds, or .qs2 read through R) -> the project's gloms/ files (table, cell ids and embeddings of the ROI heterogeneity window).

A named list. Element names are matched ignoring case, with "." = "_"; the first name found of each row is used:
  expression  logcpx | exprs_norm          genes x ROIs log-normalised expression (matrix or Matrix::dgCMatrix), gene names as row names   [required]
  metadata    metadata | meta_data         one row per ROI (data frame): an id column (glom_id_final, roi_id or glom_id; else the row names) and
                                           optionally disease_type, disease_subtype, Activity, Chronicity, N, logUMI, annotation_clean ...    [required]
  pca         pca | pca_embeddings | pca_res$embeddings     ROIs x PCs                                                                    [required]
  pca_harmony pca_harmony | H              ROIs x PCs, batch-corrected                                                                    [optional]
  umap        umap | umap_embeddings | Humap$embedding      ROIs x 2                                                                       [required]
  clusters    Humap$clusters (Clust0.5 ...)  one column per resolution                                                                  [optional]
  counts      counts | counts_mat          only checked for the same size as the expression matrix, not stored                           [optional]
ROIs are matched across the elements by their ids (column names / row names) when present, else by order (the metadata rows). The cells of an ROI are the
cells of the project's cell table whose region column (roles.roi in project.json, else glom_id, else roi_id) equals the ROI id; the sample, and disease
type / subtype when the metadata has none, are taken from those cells. ROIs with no cell in the table are left out.
"""
import json, os, shutil, subprocess, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from .importer import ImportError_
from .rds import read_rds, Matrix

ID_COLS = ("glom_id_final", "roi_id", "glom_id", "id")
CELL_COLS = ("glom_id", "roi_id")
ALIASES = {"counts": ("counts", "counts_mat"), "metadata": ("metadata", "meta", "meta_data"), "logcpx": ("logcpx", "exprs_norm"),
           "pca": ("pca", "pca_embeddings", "pca_embedding", "pca_res/embeddings"), "pca_harmony": ("pca_harmony", "harmony", "harmony_embeddings", "h"),
           "umap": ("umap", "umap_embeddings", "umap_embedding", "humap/embedding"), "clusters": ("humap/clusters",)}


def _norm(s):
    return str(s).lower().replace(".", "_")


def _pick(obj, key):
    """The first element of the list found under one of the names of `key` (a name with "/" is looked up inside a nested list); also its name."""
    for a in ALIASES[key]:
        o, ok = obj, True
        for part in a.split("/"):
            low = {_norm(k): k for k in o} if isinstance(o, dict) else {}
            if part in low:
                o = o[low[part]]
            else:
                ok = False
                break
        if ok:
            return o, a
    return None, None


def _dense(m):
    x = m.x.toarray() if sp.issparse(m.x) else np.asarray(m.x)
    return x


def _table(v, what):
    """A matrix or numeric data frame -> (array, row names, column names)."""
    if isinstance(v, Matrix):
        return _dense(v).astype(np.float64), v.rownames, v.colnames
    if isinstance(v, pd.DataFrame):
        num = v.select_dtypes("number")
        return num.values.astype(np.float64), [str(i) for i in v.index], [str(c) for c in num.columns]
    raise ImportError_(f"{what} must be a matrix or a data frame, found {type(v).__name__}")


def _num(v):
    """A text value that is a number as that number ("3" -> 3), else unchanged."""
    try:
        return _js(float(v))
    except (TypeError, ValueError):
        return _js(v)


def _js(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else (int(v) if float(v).is_integer() else float(v))
    if v is None or (not isinstance(v, str) and pd.isna(v)):
        return None
    return v if isinstance(v, (str, int, bool)) else str(v)


def _as_rds(path, tmp, log):
    """The .rds to read: the file itself, or for a .qs2 a copy made with R (qs2::qs_read -> saveRDS) in `tmp`."""
    path = Path(path)
    if path.suffix.lower() == ".rds":
        return path
    if path.suffix.lower() != ".qs2":
        raise ImportError_(f"the region_heterogeneity object must be an .rds or .qs2 file, not {path.suffix}")
    r = os.environ.get("WEBBSCOPE_RSCRIPT") or shutil.which("Rscript") or next((c for c in ("/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript", "/usr/bin/Rscript") if os.path.exists(c)), None)
    cmd = f'The file can be converted in R: saveRDS(qs2::qs_read("{path}"), "{path.with_suffix(".rds").name}")'
    if not r:
        raise ImportError_(f"a .qs2 file is read through R and Rscript was not found. {cmd}, then give the .rds (or set WEBBSCOPE_RSCRIPT to the Rscript path).")
    out = Path(tmp) / (path.stem + ".rds")
    log(f"converting {path.name} to .rds with R (qs2)")
    p = subprocess.run([r, "-e", "args <- commandArgs(TRUE); saveRDS(qs2::qs_read(args[1]), args[2])", str(path), str(out)], capture_output=True, text=True)
    if p.returncode or not out.exists():
        raise ImportError_(f"R could not read {path.name} (is the qs2 package installed?): {(p.stderr or p.stdout).strip()[-300:]}. {cmd}")
    return out


def convert(rds_path, cells, ids, roles, out_dir, log=print):
    """Write <out_dir>/gloms/gloms.json + gloms.npz from the .rds / .qs2; `cells` is the cell table in the order of `ids`. Returns the number of ROIs."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out_dir) as tmp:
        obj = read_rds(_as_rds(rds_path, tmp, log))
    return _convert(obj, cells, ids, roles, Path(out_dir), log)


def _convert(obj, cells, ids, roles, out_dir, log):
    if not isinstance(obj, dict):
        raise ImportError_("the region_heterogeneity object must be a named list (counts, metadata, logcpx, pca, umap)")
    got = {k: _pick(obj, k) for k in ALIASES}
    parts = {k: v[0] for k, v in got.items()}
    missing = [f"{k} ({' / '.join(ALIASES[k])})" for k in ("metadata", "logcpx", "pca", "umap") if parts[k] is None]
    if missing:
        raise ImportError_(f"the region_heterogeneity object has no element for: {'; '.join(missing)}. It has: {', '.join(map(str, obj))}")
    log("elements used: " + ", ".join(f"{k} = {v[1]}" for k, v in got.items() if v[0] is not None))
    md = parts["metadata"]
    if not isinstance(md, pd.DataFrame):
        raise ImportError_("metadata must be a data frame")
    idc = next((c for c in ID_COLS if c in md.columns), None)
    mids = [str(x) for x in (md[idc] if idc else md.index)]
    if len(set(mids)) != len(mids):
        raise ImportError_(f"metadata has duplicate ROI ids ({idc or 'row names'})")
    E = parts["logcpx"]
    if not isinstance(E, Matrix):
        raise ImportError_("logcpx must be a matrix (genes x ROIs)")
    x, rn, cn = _dense(E).astype(np.float32), E.rownames, E.colnames
    if rn and cn and len(set(mids) & set(rn)) > len(set(mids) & set(cn)):
        x, rn, cn = x.T, cn, rn
        log("logcpx is ROIs x genes; transposing")
    if not rn:
        raise ImportError_("logcpx has no row names (gene names)")
    if not cn:
        if x.shape[1] != len(mids):
            raise ImportError_(f"the expression matrix has no column names and {x.shape[1]} columns, but metadata has {len(mids)} rows")
        log("the expression matrix has no column names: ROIs are taken in the order of the metadata rows")
        cn = mids
    rois = [str(c) for c in cn]
    if len(set(rois)) != len(rois):
        raise ImportError_("logcpx has duplicate ROI ids as column names")
    c = parts["counts"]
    if isinstance(c, Matrix):
        if c.shape != x.shape or (c.rownames and c.rownames != rn) or (c.colnames and [str(i) for i in c.colnames] != rois):
            raise ImportError_(f"the counts element ({c.shape[0]} x {c.shape[1]}) does not match the expression matrix ({x.shape[0]} x {x.shape[1]}) in size or names")
        log("counts: same size and genes as the expression matrix (not stored)")
    absent = [r for r in rois if r not in set(mids)]
    if absent:
        raise ImportError_(f"{len(absent)} ROIs of logcpx are not in metadata (first: {absent[0]}); metadata ids come from {idc or 'its row names'}")
    md = md.assign(__id=mids).set_index("__id").loc[rois]

    def emb(key, name):
        a, r, _ = _table(parts[key], name)
        if r is not None and set(r) == set(rois):
            a = pd.DataFrame(a, index=[str(i) for i in r]).loc[rois].values
        elif a.shape[0] != len(rois):
            raise ImportError_(f"{name} has {a.shape[0]} rows but there are {len(rois)} ROIs (it needs ROI ids as row names, or one row per ROI in the order of logcpx)")
        elif r is None:
            log(f"{name}: no row names, taken in the order of logcpx")
        return a.astype(np.float32)

    pca, um = emb("pca", "pca"), emb("umap", "umap")
    if um.shape[1] < 2:
        raise ImportError_("umap needs two columns")
    har = emb("pca_harmony", "pca_harmony") if parts["pca_harmony"] is not None else None

    ccol = roles.get("roi") if roles.get("roi") in cells.columns else next((k for k in CELL_COLS if k in cells.columns), None)
    if ccol is None:
        raise ImportError_("the cell table has no column with the ROI of each cell (glom_id or roi_id; or name it under roles.roi in project.json), "
                           "so the ROIs cannot be linked to cells")
    key = cells[ccol].astype(str).where(cells[ccol].notna())
    rows = key.groupby(key).indices
    samp = cells[roles["sample"]].astype(str).values
    keep = [i for i, r in enumerate(rois) if r in rows]
    if not keep:
        raise ImportError_(f"none of the {len(rois)} ROIs has a cell in the cell table (column {ccol})")
    if len(keep) < len(rois):
        log(f"{len(rois) - len(keep)} of {len(rois)} ROIs have no cell in the cell table (column {ccol}) and are left out")
    dis = {s: cells.loc[samp == s].iloc[0] for s in dict.fromkeys(samp)}
    meta, lists = [], []
    tot = None                                              # logUMI is only taken from the metadata
    cl = parts["clusters"] if isinstance(parts["clusters"], pd.DataFrame) and len(parts["clusters"]) == len(rois) else None
    if parts["clusters"] is not None and cl is None:
        log("clusters: not one row per ROI, ignored")
    num = {c: pd.to_numeric(md[c], errors="coerce") for c in md.columns if not pd.api.types.is_numeric_dtype(md[c])}
    for c_, v in num.items():                                # text columns that are numbers ("0", "3"): Activity, Chronicity ...
        if v.notna().sum() == md[c_].notna().sum() and c_ not in ID_COLS:
            md[c_] = v
    for i in keep:
        r = rois[i]
        idx = rows[r]
        s = pd.Series(samp[idx]).mode().iloc[0]
        m = {k: _js(v) for k, v in md.loc[r].items() if not k.startswith("__")}
        m.setdefault("glom_id_final", r); m["glom_id_final"] = r
        m["sample"] = s
        for k in ("disease_type", "disease_subtype"):
            if m.get(k) is None:
                m[k] = _js(dis[s][k]) if k in dis[s].index else None
        for lo, hi in (("activity", "Activity"), ("chronicity", "Chronicity")):
            if hi not in m:
                m[hi] = next((m[k] for k in m if k.lower() == lo), None)
        if "N" not in m:
            m["N"] = None
        if "logUMI" not in m:
            m["logUMI"] = float(np.log(tot[i])) if tot is not None and tot[i] > 0 else None
        if cl is not None:
            for c_ in cl.columns:
                m["cluster_" + (("res" + c_[5:]) if c_.startswith("Clust") else c_)] = _num(cl[c_].iloc[i])
        meta.append(m)
        lists.append([ids[j] for j in idx])
    keep = np.array(keep)
    d = out_dir / "gloms"
    d.mkdir(parents=True, exist_ok=True)
    (d / "gloms.json").write_text(json.dumps({"genes": list(rn), "meta": meta, "cells": lists}))
    z = {"expr": x[:, keep], "pca": pca[keep], "umap": um[keep][:, :2]}
    if har is not None:
        z["pca_harmony"] = har[keep]
    np.savez_compressed(d / "gloms.npz", **z)
    return len(keep)
