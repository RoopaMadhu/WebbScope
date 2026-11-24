"""Read an AnnData .h5ad file with h5py alone (anndata and scanpy are not needed).

Read: the cell table (obs, with its index as a column), the gene names (var index or a var column), one matrix (X, a layer, or raw/X; sparse
CSR / CSC or dense), the spatial coordinates in obsm (as extra columns <key>_x, <key>_y of the cell table), the UMAP / PCA in obsm and cell
outlines kept as WKT text in an obs column (found by its values; such columns are outlines, so they are not part of the cell table).
Files written by anndata 0.7 or later are supported (obs stored as a group); older ones (obs as one compound table) are refused with a message.
"""
import re
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
from .importer import ImportError_

UMAP_KEYS = ("X_umap_harmony", "X_umap")
PCA_KEYS = ("X_pca_harmony", "X_harmony", "X_pca")
WKT_POLYGON = re.compile(r"\s*(MULTI)?POLYGON\s*(Z\s*)?(\(|EMPTY)", re.I)
WKT_ANY = re.compile(r"\s*((MULTI)?(POINT|LINESTRING|POLYGON)|GEOMETRYCOLLECTION)\s*(Z\s*)?(\(|EMPTY)", re.I)
NO_TEXT = ("", "nan", "none", "na", "<na>")               # a missing value written as text
OUTLINE_NAMES = ("geometry", "shape", "geom")


def _attr(node, k, default=None):
    v = node.attrs.get(k, default)
    return v.decode() if isinstance(v, bytes) else v


def _values(ds):
    """A 1-d dataset as a numpy array; strings as str objects. Also a nullable string array (anndata 0.11+, a group of values and mask;
    the missing entries become None)."""
    if isinstance(ds, h5py.Group):
        if _attr(ds, "encoding-type", "") == "nullable-string-array" and "values" in ds and "mask" in ds:
            v = _values(ds["values"]).astype(object)
            v[ds["mask"][...].astype(bool)] = None
            return v
        raise ValueError("unsupported group")
    if h5py.check_string_dtype(ds.dtype) is not None:
        return np.asarray(ds.asstr()[...], dtype=object)
    if ds.dtype.kind == "O":
        raise ValueError("unsupported object column")
    return ds[...]


def _column(g, name, f):
    """One column of a dataframe group as something pandas takes; None when its encoding is not supported."""
    node = g[name]
    enc = _attr(node, "encoding-type", "")
    try:
        if isinstance(node, h5py.Group):
            if enc == "categorical":
                cats = _values(node["categories"])
                return pd.Categorical.from_codes(node["codes"][...].astype(np.int64), categories=pd.Index(cats.astype(str) if cats.dtype == object else cats))
            if enc == "nullable-string-array":
                return _values(node)
            if enc in ("nullable-integer", "nullable-boolean"):
                v, mask = _values(node["values"]), node["mask"][...].astype(bool)
                if enc == "nullable-integer":
                    return pd.arrays.IntegerArray(v.astype(np.int64), mask)
                return pd.arrays.BooleanArray(v.astype(bool), mask)
            return None
        v = _values(node)
        if v.ndim != 1:
            return None
        ref = node.attrs.get("categories")                  # anndata 0.7: codes here, the levels in obs/__categories/<name>
        if ref is not None and isinstance(ref, h5py.Reference):
            cats = _values(f[ref])
            return pd.Categorical.from_codes(v.astype(np.int64), categories=cats)
        return v
    except (ValueError, TypeError, KeyError):
        return None


def _order(g, idx):
    """The columns of a dataframe group, in order."""
    order = g.attrs.get("column-order")
    return [o.decode() if isinstance(o, bytes) else str(o) for o in order] if order is not None and len(order) else [k for k in g.keys() if k not in (idx, "__categories")]


def _frame(f, path, leave_out=()):
    """A dataframe group (obs, var, raw/var) -> (DataFrame of the readable columns, index values as str, names of the skipped columns)."""
    g = f[path]
    if not isinstance(g, h5py.Group):
        raise ImportError_(f"this .h5ad was written by an old anndata version (before 0.7: {path} is one table). Open it with anndata and save it again (adata.write_h5ad), then try again.")
    idx = _attr(g, "_index", "_index")
    if idx not in g:
        raise ImportError_(f"{path} has no index in this .h5ad")
    index = np.asarray(_values(g[idx])).astype(str)
    data, skipped = {}, []
    for c in _order(g, idx):
        if c not in g or c in leave_out:
            continue
        v = _column(g, c, f)
        if v is None or len(v) != len(index):
            skipped.append(c)
        else:
            data[c] = v
    return pd.DataFrame(data), index, skipped


def n_obs(f):
    """Number of cells (the length of the obs index; a nullable string index is a group of values and mask)."""
    node = f["obs"][_attr(f["obs"], "_index", "_index")]
    return len(node["values"]) if isinstance(node, h5py.Group) else len(node)


def _text_source(f, node):
    """The string dataset holding a column's text (its levels when categorical); None for a column that holds no text."""
    if isinstance(node, h5py.Group):
        node = node.get("categories") if _attr(node, "encoding-type", "") == "categorical" else node.get("values")
    elif isinstance(node.attrs.get("categories"), h5py.Reference):   # anndata 0.7: codes here, the levels elsewhere
        node = f[node.attrs["categories"]]
    return node if isinstance(node, h5py.Dataset) and node.ndim == 1 and h5py.check_string_dtype(node.dtype) is not None else None


def _holds_wkt(f, node, step=5000):
    """True when the first values that are not missing are all WKT and at least one is a polygon (a column of WKT points is not an outline)."""
    d = _text_source(f, node)
    if d is None:
        return False
    for a in range(0, len(d), step):
        v = [x for x in d.asstr()[a:a + step] if x.strip().lower() not in NO_TEXT]
        if v:
            return all(WKT_ANY.match(x) for x in v) and any(WKT_POLYGON.match(x) for x in v)
    return False


def outline_columns(f):
    """The obs columns holding WKT polygons: {"all": every one, "cell": the one to draw, or None}. A column named for nuclei is never drawn as
    the cell's outline; among the others a plain geometry / shape / geom name comes first, then a name with cell or bound in it."""
    g = f["obs"]
    if not isinstance(g, h5py.Group):
        return {"all": [], "cell": None}
    idx = _attr(g, "_index", "_index")
    found = [c for c in _order(g, idx) if c in g and _holds_wkt(f, g[c])]
    cell = sorted((c for c in found if "nuc" not in c.lower()), key=lambda c: (c.lower() not in OUTLINE_NAMES, not ("cell" in c.lower() or "bound" in c.lower())))
    return {"all": found, "cell": cell[0] if cell else None}


def outlines(f, column):
    """The WKT text of one obs column, one entry per cell, as objects (None = missing)."""
    v = _column(f["obs"], column, f)
    if v is None:
        raise ImportError_(f"the outline column {column!r} of this .h5ad cannot be read")
    return pd.Series(v).astype(object).where(lambda x: x.notna(), None).values


def obs(f):
    """(cell table, name of the column holding the obs index, skipped columns). The spatial arrays of obsm are added as <key>_x / <key>_y;
    columns of WKT outlines are left out."""
    df, index, skipped = _frame(f, "obs", outline_columns(f)["all"])
    idcol = "cell_id" if "cell_id" not in df.columns else "obs_names"
    df.insert(0, idcol, index)
    for key in spatial_keys(f):
        a = f["obsm"][key]
        for j, ax in enumerate("xy"):
            c = f"{key}_{ax}"
            if c not in df.columns:
                df[c] = a[:, j].astype(np.float64)
    return df, idcol, skipped


def obs_columns(f):
    """Column names of the cell table obs() would give, without reading the values (for the question asked when opening)."""
    g = f["obs"]
    if not isinstance(g, h5py.Group):
        raise ImportError_("this .h5ad was written by an old anndata version (before 0.7). Open it with anndata and save it again (adata.write_h5ad), then try again.")
    idx = _attr(g, "_index", "_index")
    wkt = outline_columns(f)["all"]
    cols = [c for c in _order(g, idx) if c not in wkt]
    idcol = "cell_id" if "cell_id" not in cols else "obs_names"
    return [idcol] + cols + [f"{k}_{ax}" for k in spatial_keys(f) for ax in "xy"], idcol


def spatial_keys(f):
    """obsm arrays that hold coordinates: 2-d, at least 2 columns, a name with spatial / centroid / coord in it."""
    if "obsm" not in f:
        return []
    n, out = n_obs(f), []
    for k, a in f["obsm"].items():
        if isinstance(a, h5py.Dataset) and a.ndim == 2 and a.shape[0] == n and a.shape[1] >= 2 and any(w in k.lower() for w in ("spatial", "centroid", "coord")):
            out.append(k)
    return out


def embeddings(f):
    """{"umap": (obsm key, (n, 2) array), "pca": (key, (n, k) array)} for what the file has."""
    out = {}
    if "obsm" not in f:
        return out
    n = n_obs(f)
    for kind, keys in (("umap", UMAP_KEYS), ("pca", PCA_KEYS)):
        for k in keys:
            a = f["obsm"].get(k)
            if isinstance(a, h5py.Dataset) and a.ndim == 2 and a.shape[0] == n and a.shape[1] >= 2:
                out[kind] = (k, a[:, :2] if kind == "umap" else a[...])
                break
    return out


def matrices(f):
    """The matrices that can be used as expression: [{key, kind (csr / csc / dense), shape, n_genes}]; key is "X", "layers/<name>" or "raw/X"."""
    out = []
    n = n_obs(f)
    cands = [("X", f.get("X"))] + [(f"layers/{k}", v) for k, v in (f["layers"].items() if "layers" in f else [])] + [("raw/X", f["raw"].get("X") if "raw" in f else None)]
    for key, node in cands:
        if node is None:
            continue
        if isinstance(node, h5py.Group):
            enc = _attr(node, "encoding-type", "") or _attr(node, "h5sparse_format", "")
            kind = {"csr_matrix": "csr", "csc_matrix": "csc", "csr": "csr", "csc": "csc"}.get(enc)
            shape = tuple(int(x) for x in node.attrs.get("shape", ()))
            if kind is None or len(shape) != 2 or not all(k in node for k in ("data", "indices", "indptr")):
                continue
        elif isinstance(node, h5py.Dataset) and node.ndim == 2:
            kind, shape = "dense", tuple(int(x) for x in node.shape)
        else:
            continue
        if shape[0] != n:
            continue
        out.append({"key": key, "kind": kind, "shape": list(shape), "n_genes": shape[1]})
    return out


def _node(f, key):
    if key not in f:
        raise ImportError_(f"this .h5ad has no matrix {key}")
    return f[key]


def peek(f, key, k=200_000):
    """A sample of the stored values of a matrix (the first k stored values): {integer, negative, max}."""
    node = _node(f, key)
    if isinstance(node, h5py.Group):
        v = node["data"][:k]
    else:
        rows = max(1, min(node.shape[0], k // max(1, node.shape[1])))
        v = node[:rows].ravel()
        v = v[v != 0]
    v = np.asarray(v, np.float64)
    if not len(v):
        return {"integer": True, "negative": False, "max": 0.0}
    return {"integer": bool(np.allclose(v, np.round(v))), "negative": bool((v < 0).any()), "max": float(v.max())}


def gene_columns(f, key):
    """String columns of the var table that belongs to the matrix (raw/var for raw/X), for choosing the gene names."""
    df, _, _ = _frame(f, "raw/var" if key == "raw/X" else "var")
    return [c for c in df.columns if pd.api.types.is_string_dtype(df[c].dtype) or df[c].dtype == object or isinstance(df[c].dtype, pd.CategoricalDtype)]


def gene_names(f, key, column=None):
    df, index, _ = _frame(f, "raw/var" if key == "raw/X" else "var")
    if not column:
        return [str(x) for x in index]
    if column not in df.columns:
        raise ImportError_(f"the gene table of this .h5ad has no column {column!r}")
    return [str(x) for x in df[column].astype(str)]


def read_matrix(f, key, log=print):
    """The matrix as genes x cells CSR (the cache layout), float64 values."""
    node = _node(f, key)
    m = next((x for x in matrices(f) if x["key"] == key), None)
    if m is None:
        raise ImportError_(f"{key} of this .h5ad is not a matrix of cells x genes")
    n, g = m["shape"]
    if m["kind"] == "dense":
        step, parts = max(1, int(5e7 // max(1, g))), []
        for a in range(0, n, step):
            parts.append(sp.csr_matrix(np.asarray(node[a:a + step], np.float64)))
        return sp.vstack(parts).T.tocsr() if parts else sp.csr_matrix((g, n))
    data, ind, ptr = np.asarray(node["data"][...], np.float64), node["indices"][...], node["indptr"][...]
    if m["kind"] == "csr":                                  # cells x genes rows = cells
        return sp.csr_matrix((data, ind, ptr), shape=(n, g)).T.tocsr()
    return sp.csc_matrix((data, ind, ptr), shape=(n, g)).T.tocsr()   # CSC of cells x genes = CSR of genes x cells


def guess(f, mats):
    """(matrix with the raw counts, matrix with log-normalised values), suggested from the stored values; either may be None.
    Raw counts = whole, non-negative numbers (looked for in layers/counts, raw/X, X, then the other layers); log-normalised = non-negative
    numbers that are not all whole (X first, then the layers). Negative values (scaled data) are neither."""
    pk = {m["key"]: (m if "integer" in m else {**m, **peek(f, m["key"])}) for m in mats}
    first = ("layers/counts", "raw/X", "X")
    order_c = [k for k in first if k in pk] + [k for k in pk if k not in first]
    order_l = (["X"] if "X" in pk else []) + [k for k in pk if k != "X"]
    ck = next((k for k in order_c if pk[k]["integer"] and not pk[k]["negative"]), None)
    lk = next((k for k in order_l if not pk[k]["integer"] and not pk[k]["negative"]), None)
    return ck, lk
