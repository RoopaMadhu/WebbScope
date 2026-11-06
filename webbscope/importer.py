"""Project folder -> fast internal cache (the layout the analysis engine reads).

Required inputs : project.json, matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz, cells.parquet, cell_boundaries.parquet
Optional inputs : regions.geojson            pathologist outlines (cell coordinates)  -> Regions, Region neighbourhoods, Glomeruli, region filters
                  embeddings.parquet         UMAP / PCA coordinates                    -> Embeddings
                  gloms/gloms.json + .npz    ROI-level object                          -> ROI heterogeneity (or give a region_heterogeneity .rds when opening: roi_object.py)
                  signatures.json            curated gene signatures                   -> Signatures, violins
                  scffpe_composition.json    scFFPE cell-type counts                   -> Xenium vs scFFPE
                  transcripts/               tx_x.npy tx_y.npy tx_ptr.npy index.json (sorted by tissue, gene)   -> transcript overlay (read in place)
                  he_slide.json (+ he/<key>/ tiles, he_slide_annotations.geojson)   H&E slide scan, not registered   -> H&E slide window
                  he/, he_gloms/             pre-tiled H&E pyramids and glom crops     -> H&E views (read from the project folder)

The cache lives in the user's cache directory (override with WEBBSCOPE_CACHE), never inside the project folder, so the input data is never
modified. It is rebuilt automatically when any input, or this format version, changes.
"""
import collections, gzip, hashlib, json, os, shutil, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.io import mmread

FORMAT = 3
REQUIRED = ["project.json", "matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "cells.parquet", "cell_boundaries.parquet"]
OPTIONAL = ["matrix.npz", "regions.geojson", "embeddings.parquet", "gloms/gloms.json", "gloms/gloms.npz", "signatures.json", "scffpe_composition.json",
            "he/index.json", "he_gloms/index.json", "transcripts/index.json", "transcripts/tx_x.npy", "transcripts/tx_ptr.npy", "he_slide.json", "he_slide_annotations.geojson"]
OUTSIDE = "(outside all outlines)"
ROLE_DEFAULTS = {"cell_id": "cell_id", "x": "x", "y": "y", "sample": "sample", "lib_size": "lib_size"}


class ImportError_(Exception):
    """Problem with the user's input files; the message is shown to the user."""


def cache_root():
    return Path(os.environ.get("WEBBSCOPE_CACHE") or Path.home() / ".webbscope" / "cache")


def cache_dir_for(project):
    project = Path(project).resolve()
    return cache_root() / f"{project.name}-{hashlib.sha1(str(project).encode()).hexdigest()[:10]}"


def _fingerprint(project, roi_object=None, cols=None):
    fp = {}
    if cols:
        fp["cols"] = {k: cols[k] for k in sorted(cols)}                      # sample / donor columns chosen when opening
    if roi_object:
        st = Path(roi_object).stat()
        fp["roi_object"] = [str(Path(roi_object).resolve()), st.st_size, int(st.st_mtime)]
    for f in REQUIRED + OPTIONAL:
        p = project / f
        if p.exists():
            st = p.stat()
            fp[f] = [st.st_size, int(st.st_mtime)]
    return fp


def is_current(project, cdir, roi_object=None, cols=None):
    info = cdir / "cache_info.json"
    if not info.exists():
        return False
    try:
        d = json.loads(info.read_text())
    except Exception:
        return False
    return d.get("format") == FORMAT and d.get("inputs") == _fingerprint(project, roi_object, cols)


def _lines(path):
    with gzip.open(path, "rt") as f:
        return [l.rstrip("\n").split("\t")[0] for l in f if l.strip()]


def _regions(project, cells, sample_col, x_col, y_col, cdir, log):
    """regions.geojson -> paths_by_sample.json (outlines per sample), path_mask.npy (bit i = cell inside annotation i; last bit = inside none)
    and path_regions.json. Same rules as the original browser: regions flagged negative_roa are subtracted from the same annotation's other
    polygons; a cell can be inside several nested annotations."""
    import shapely
    from shapely.geometry import shape
    from shapely.validation import make_valid
    gj = json.loads((project / "regions.geojson").read_text())
    feats = [f for f in gj["features"] if f.get("geometry")]
    if not feats:
        raise ImportError_("regions.geojson has no polygons")
    sid = cells[sample_col].astype(str).values
    xy_all = cells[[x_col, y_col]].values.astype(float)
    rows_of = {s: np.flatnonzero(sid == s) for s in dict.fromkeys(sid)}
    by_s = collections.defaultdict(list)
    for f in feats:
        p = f.get("properties") or {}
        s = str(p.get("sample"))
        if s not in rows_of:
            raise ImportError_(f"regions.geojson: sample '{s}' is not in the cell table")
        by_s[s].append(f)
    log(f"regions: {len(feats)} polygons in {len(by_s)} sample(s)")
    res = {"samples": {}}
    member, totals, colours = {}, collections.Counter(), collections.defaultdict(collections.Counter)
    for s, fs in by_s.items():
        fs = [{**f, "properties": {k: v for k, v in f["properties"].items() if k != "sample"}} for f in fs]       # `sample` is an input-only key
        res["samples"][s] = {"tissues": sorted({str((f["properties"] or {}).get("tissue_id")) for f in fs}), "iou": [], "features": fs}
        xy = xy_all[rows_of[s]]
        groups = collections.defaultdict(lambda: ([], []))
        for f in fs:
            p = f["properties"]
            neg = bool(p.get("negative_roa"))
            groups[p["annotation_name"]][1 if neg else 0].append(make_valid(shape(f["geometry"])))
            if not neg and p.get("line_color_hex"):
                colours[p["annotation_name"]][p["line_color_hex"]] += 1
        member[s] = {}
        for name, (pos, neg) in groups.items():
            if not pos:
                continue
            reg = shapely.union_all(pos)
            if neg:
                reg = reg.difference(shapely.union_all(neg))
            inside = np.zeros(len(xy), bool)
            x0, y0, x1, y1 = reg.bounds
            m = (xy[:, 0] >= x0) & (xy[:, 0] <= x1) & (xy[:, 1] >= y0) & (xy[:, 1] <= y1)
            if m.any():
                inside[m] = shapely.contains_xy(reg, xy[m, 0], xy[m, 1])
            member[s][name] = inside
            totals[name] += int(inside.sum())
    names = [n for n, _ in totals.most_common()]
    if len(names) + 1 > 16:
        raise ImportError_("regions.geojson has more than 15 annotation types (limit of the region mask)")
    mask = np.zeros(len(cells), np.uint16)
    for s, mem in member.items():
        r = rows_of[s]
        m = np.zeros(len(r), np.uint16)
        for i, n in enumerate(names):
            if n in mem:
                m |= (mem[n].astype(np.uint16) << i)
        mask[r] = m
    outside = mask == 0
    mask[outside] = np.uint16(1 << len(names))
    np.save(cdir / "path_mask.npy", mask)
    cols = [colours[n].most_common(1)[0][0] if colours[n] else "#888888" for n in names]
    (cdir / "path_regions.json").write_text(json.dumps({"names": names + [OUTSIDE], "colors": cols + ["#8b95a5"]}))
    (cdir / "paths_by_sample.json").write_text(json.dumps(res, separators=(",", ":")))
    return {"annotations": {n: int(totals[n]) for n in names}, "cells_inside_none": int(outside.sum())}


def _embeddings(project, ids, cdir, log):
    e = pd.read_parquet(project / "embeddings.parquet")
    if "cell_id" not in e.columns:
        raise ImportError_("embeddings.parquet needs a cell_id column")
    e = e.set_index("cell_id")
    absent = [i for i in ids if i not in e.index]
    if absent:
        raise ImportError_(f"embeddings.parquet is missing {len(absent)} cells (first: {absent[0]})")
    e = e.loc[ids]
    out = {}
    for name in sorted({c[5:-2] for c in e.columns if c.startswith("umap_") and c.endswith("_1")}):
        np.save(cdir / f"umap_{name}.npy", e[[f"umap_{name}_1", f"umap_{name}_2"]].values.astype(np.float32)); out[f"umap_{name}"] = 2
    for name in sorted({c[4:].rsplit("_", 1)[0] for c in e.columns if c.startswith("pca_")}):
        cols = sorted([c for c in e.columns if c.startswith(f"pca_{name}_")], key=lambda c: int(c.rsplit("_", 1)[1]))
        np.save(cdir / f"pca_{name}.npy", e[cols].values.astype(np.float32)); out[f"pca_{name}"] = len(cols)
    log("embeddings: " + ", ".join(f"{k} ({v}d)" for k, v in out.items()))
    return out


def _codes(series):
    """Level order and uint16 codes exactly as engine.store.Store._add_cat computes them (levels by frequency)."""
    s = series.astype("string")
    levels = s.value_counts().index.tolist()
    lut = {l: i for i, l in enumerate(levels)}
    return levels, s.map(lut).fillna(65535).astype(np.uint16).values


def _gloms(project, ids, cells, cdir, log):
    """gloms/gloms.json (+ gloms.npz): glom table, the cell ids of each glom, glom-level expression / PCA / UMAP.
    The per-glom cell-type and cell-state counts are NOT taken from the file: they are recomputed here from the cell table with the
    engine's own level order, so glom colouring by type/state can never disagree with the viewer's levels."""
    from .engine.config import load_mapping
    j = json.loads((project / "gloms" / "gloms.json").read_text())
    z = np.load(project / "gloms" / "gloms.npz")
    n = len(j["meta"])
    if len(j["cells"]) != n or z["expr"].shape[1] != n or z["umap"].shape[0] != n:
        raise ImportError_("gloms: meta, cells, expr and umap do not describe the same number of gloms")
    pos = pd.Index(ids).get_indexer
    rows = []
    for g, lst in zip(j["meta"], j["cells"]):
        r = pos(lst)
        if (r < 0).any():
            raise ImportError_(f"gloms.json: glom {g.get('glom_id_final')} lists cells that are not in the cell table")
        rows.append(r)
    mp = load_mapping(cells.columns, cdir / "mapping.json")
    lt, ct = _codes(cells[mp["cell_type"]])
    ls, cs = _codes(cells[mp["cell_state"]])
    Ct, Cs = np.zeros((n, len(lt)), np.float32), np.zeros((n, len(ls)), np.float32)
    for i, r in enumerate(rows):
        Ct[i] = np.bincount(ct[r][ct[r] < 65535], minlength=len(lt))
        Cs[i] = np.bincount(cs[r][cs[r] < 65535], minlength=len(ls))
    for g, r in zip(j["meta"], rows):
        g["n_viewer_cells"] = int(len(r))
    (cdir / "gloms").mkdir(exist_ok=True)
    (cdir / "gloms" / "gloms.json").write_text(json.dumps({"genes": j["genes"], "meta": j["meta"], "cells": [r.tolist() for r in rows],
                                                          "type_levels": lt, "state_levels": ls}))
    emb = {k: z[k] for k in ("pca_harmony", "pca") if k in z.files}                 # the batch-corrected PCA is optional
    np.savez_compressed(cdir / "gloms" / "gloms.npz", expr=z["expr"], umap=z["umap"], type_counts=Ct, state_counts=Cs, **emb)
    log(f"gloms: {n} glomeruli")
    return n


def import_project(project, cdir, log=print, roi_object=None, cols=None):
    t0 = time.time()
    project, cdir = Path(project), Path(cdir)
    has_npz = (project / "matrix.npz").exists()      # scipy.sparse.save_npz, genes x cells: a fast alternative to matrix.mtx.gz (used by the merged batch projects)
    missing = [f for f in REQUIRED if not (project / f).exists() and not (has_npz and f == "matrix.mtx.gz")]
    if missing:
        raise ImportError_("This folder is not a project: missing " + ", ".join(missing))
    man = json.loads((project / "project.json").read_text())
    roles = {**ROLE_DEFAULTS, **man.get("roles", {})}

    log("reading gene and cell lists")
    genes, ids = _lines(project / "features.tsv.gz"), _lines(project / "barcodes.tsv.gz")
    if len(set(ids)) != len(ids):
        raise ImportError_("barcodes.tsv.gz contains duplicate cell ids")
    if len(set(genes)) != len(genes):
        raise ImportError_("features.tsv.gz contains duplicate gene names")

    log("reading counts matrix")
    m = sp.load_npz(project / "matrix.npz") if has_npz else mmread(project / "matrix.mtx.gz")
    m = m.tocsr() if sp.issparse(m) else sp.csr_matrix(m)
    if m.shape == (len(ids), len(genes)) and m.shape != (len(genes), len(ids)):
        log("matrix is cells x genes; transposing")
        m = m.T.tocsr()
    if m.shape != (len(genes), len(ids)):
        raise ImportError_(f"matrix is {m.shape[0]} x {m.shape[1]} but there are {len(genes)} genes and {len(ids)} cells")
    if m.nnz and not np.allclose(m.data, np.round(m.data)):
        raise ImportError_("counts matrix has non-integer values (expected raw counts)")
    if m.nnz and m.data.min() < 0:
        raise ImportError_("counts matrix has negative values")

    log("reading cell table")
    cells = pd.read_parquet(project / "cells.parquet")
    from .engine.config import find_columns, suggest
    cols = cols or {}
    if cols.get("sample"):
        roles["sample"] = cols["sample"]                                          # chosen by the user when opening the project
    elif roles["sample"] not in cells.columns and suggest(cells.columns)["sample"]:
        roles["sample"] = suggest(cells.columns)["sample"]                        # not set in project.json: the column whose name says sample
        log(f"sample column: {roles['sample']} (found by name)")
    if "donor" in cols:
        roles["donor"] = cols["donor"] or None                                    # "" = the user said there is no donor column
    elif roles.get("donor") not in cells.columns:
        roles["donor"] = suggest(cells.columns)["donor"]
    for r in ("cell_id", "x", "y", "sample"):
        if roles[r] not in cells.columns:
            raise ImportError_(f"cells.parquet has no column '{roles[r]}' (role: {r}); set it under 'roles' in project.json" + (f" or choose it when opening the project. Columns with 'sample' in the name: {find_columns(cells.columns, 'sample')[:6] or 'none'}" if r == "sample" else ""))
    if cells[roles["cell_id"]].duplicated().any():
        raise ImportError_("cells.parquet has duplicate cell ids")
    cells = cells.set_index(roles["cell_id"], drop=False)
    absent = [i for i in ids if i not in cells.index]
    if absent:
        raise ImportError_(f"{len(absent)} cells in the matrix are not in cells.parquet (first: {absent[0]})")
    cells = cells.loc[ids].reset_index(drop=True)
    if cells[[roles["x"], roles["y"]]].isna().any().any():
        raise ImportError_("cells.parquet has missing x/y coordinates")
    # every sample must be one contiguous block of rows (the analyses address a sample as a row range): stable sort by first appearance
    samp = cells[roles["sample"]].astype(str).values
    rank = {s: i for i, s in enumerate(dict.fromkeys(samp))}
    order = np.argsort(np.array([rank[s] for s in samp]), kind="stable")
    if not (order == np.arange(len(order))).all():
        log("cells of one sample were not contiguous; reordering by sample")
        cells = cells.iloc[order].reset_index(drop=True)
        ids = [ids[i] for i in order]
        m = m[:, order]
    m = m.tocsr()
    m.sum_duplicates(); m.sort_indices()
    m.data = np.round(m.data).astype(np.float32)           # float32 like the original cache: same arithmetic, same numbers (counts < 2^24 are exact)

    log("reading cell boundaries")
    bd = pd.read_parquet(project / "cell_boundaries.parquet", columns=["cell_id", "vertex_x", "vertex_y"])
    pos = pd.Index(ids).get_indexer(bd["cell_id"].astype(str))
    ok = pos >= 0
    bd, pos = bd[ok], pos[ok]
    o = np.argsort(pos, kind="stable")                      # keeps each cell's vertex order
    xy = np.column_stack([bd["vertex_x"].values[o], bd["vertex_y"].values[o]]).astype(np.float32)
    cnt = np.bincount(pos, minlength=len(ids))
    off = np.r_[0, np.cumsum(cnt)].astype(np.int64)

    colsum = np.asarray(m.sum(0)).ravel().astype(np.float32)
    lib_mismatch = None
    if roles["lib_size"] in cells.columns:
        lib_mismatch = int((~np.isclose(cells[roles["lib_size"]].values.astype(np.float32), colsum)).sum())
    else:
        cells[roles["lib_size"]] = colsum                   # library size = total counts over the panel genes

    cdir.mkdir(parents=True, exist_ok=True)
    for stale in ("path_mask.npy", "path_regions.json", "paths_by_sample.json", "signatures.json", "scffpe_composition.json", "gloms",
                  "umap_harmony.npy", "pca_harmony.npy"):
        p = cdir / stale
        shutil.rmtree(p) if p.is_dir() else p.unlink(missing_ok=True)
    log("writing cache")
    np.save(cdir / "csr_indptr.npy", m.indptr.astype(np.int64))
    np.save(cdir / "csr_indices.npy", m.indices.astype(np.int32))
    np.save(cdir / "csr_data.npy", m.data.astype(np.float32))
    np.save(cdir / "poly_offsets.npy", off)
    np.save(cdir / "poly_xy.npy", xy)
    cells.to_parquet(cdir / "cells.parquet", index=False)
    (cdir / "genes.json").write_text(json.dumps(genes))
    (cdir / "mapping.json").write_text(json.dumps({k: v for k, v in man.get("roles", {}).items() if v in cells.columns} | {"lib_size": roles["lib_size"], "sample": roles["sample"]}
                                                  | ({"donor": roles["donor"]} if roles.get("donor") in cells.columns else {})))
    (cdir / "cache_meta.json").write_text(json.dumps({
        "n_genes": len(genes), "n_cells": len(ids), "nnz": int(m.nnz), "sample_col": roles["sample"],
        "orientation": "genes x cells (rows=genes), stored CSR", "layers": ["counts"], "norm_scale": float(man.get("norm_scale", 10000.0)),
        "note": "'logcpx' = log1p(counts / lib_size * norm_scale), lib_size = total counts over the panel genes."}))

    extras = {}
    if (project / "regions.geojson").exists():
        extras["regions"] = _regions(project, cells, roles["sample"], roles["x"], roles["y"], cdir, log)
    if (project / "embeddings.parquet").exists():
        extras["embeddings"] = _embeddings(project, ids, cdir, log)
    if roi_object:
        from . import roi_object as ro
        src = cdir / "_roi_src"
        shutil.rmtree(src, ignore_errors=True)
        log(f"reading the region_heterogeneity object {Path(roi_object).name}")
        try:
            ro.convert(roi_object, cells, ids, roles, src, log)
            extras["gloms"] = _gloms(src, ids, cells, cdir, log)
        finally:
            shutil.rmtree(src, ignore_errors=True)
    elif (project / "gloms" / "gloms.json").exists():
        if not (project / "gloms" / "gloms.npz").exists():
            raise ImportError_("gloms/gloms.json needs gloms/gloms.npz")
        extras["gloms"] = _gloms(project, ids, cells, cdir, log)
    for f in ("signatures.json", "scffpe_composition.json"):
        if (project / f).exists():
            shutil.copyfile(project / f, cdir / f); extras[f] = True
    if (project / "transcripts" / "index.json").exists():
        extras["transcripts"] = True
    if (project / "he_slide.json").exists():
        extras["he_slide"] = True
    for d in ("he", "he_gloms"):
        if (project / d / "index.json").exists():
            extras[d] = True
    report = {"n_cells": len(ids), "n_genes": len(genes), "n_samples": len(rank), "nnz": int(m.nnz),
              "cells_without_boundary": int((cnt == 0).sum()), "cells_with_degenerate_boundary": int(((cnt > 0) & (cnt < 3)).sum()),
              "boundary_rows_for_unknown_cells": int((~ok).sum()), "lib_size_from_table_differs_from_matrix_column_sum": lib_mismatch,
              "genes_with_zero_counts": int((np.diff(m.indptr) == 0).sum()), "optional_inputs_found": extras, "seconds": round(time.time() - t0, 2)}
    (cdir / "cache_info.json").write_text(json.dumps({"format": FORMAT, "inputs": _fingerprint(project, roi_object, cols), "manifest": man,
                                                      "project": str(project.resolve()), "report": report}, indent=1))
    log(f"imported {len(ids):,} cells x {len(genes):,} genes in {report['seconds']} s")
    return report


def ensure_cache(project, log=print, force=False, roi_object=None, cols=None):
    cdir = cache_dir_for(project)
    if not force and is_current(Path(project), cdir, roi_object, cols):
        log("using cached import")
        return cdir
    import_project(project, cdir, log, roi_object, cols)
    return cdir
