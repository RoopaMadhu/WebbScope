"""Project folder -> fast internal cache (the layout the analysis engine reads).

Required inputs : either a matrix + features.tsv.gz, barcodes.tsv.gz and the cell metadata (cells_metadata.parquet or cells.parquet),
                  or one AnnData .h5ad file (read by h5ad.py). The matrix: counts.mtx or matrix.mtx = raw counts (matrix.mtx as in 10x's
                  cell_feature_matrix), counts_logcpx.mtx = log-normalised values; one or both. Each also .mtx.gz or .npz (a sparse matrix
                  written by scipy.sparse.save_npz). The names are not bound to 10x's: both spellings are accepted, one of each kind per folder.
Optional inputs : project.json               which column is which (roles), the normalisation scale; without it the columns are chosen when opening
                  cell_boundaries.parquet    cell outlines (without it the cells are drawn as dots): one row per outline point (cell_id,
                                             vertex_x, vertex_y, as Xenium writes it) or one row per cell with its outline as GeoParquet (WKB)
                                             or WKT text in a geometry / shape / geom column (read by read_boundaries); an .h5ad
                                             may instead keep WKT outlines in an obs column (read by boundaries_from_wkt)
                  regions.geojson            pathologist outlines (cell coordinates)  -> Regions, Region neighbourhoods, Glomeruli, region filters
                  embeddings.parquet         UMAP / PCA coordinates                    -> Embeddings
                  gloms/gloms.json + .npz    ROI-level object                          -> ROI heterogeneity (or give a region_heterogeneity .rds when opening: roi_object.py)
                  signatures.json            curated gene signatures                   -> Signatures, violins
                  scffpe_composition.json    scFFPE cell-type counts                   -> Xenium vs scFFPE
                  transcripts/               tx_x.npy tx_y.npy tx_ptr.npy index.json (sorted by tissue, gene)   -> transcript overlay (read in place)
                  he_slide.json (+ he/<key>/ tiles, he_slide_annotations.geojson)   H&E slide scan, not registered   -> H&E slide window
                  he/, he_gloms/             pre-tiled H&E pyramids and glom crops     -> H&E views (read from the project folder)

Raw counts are log-normalised here as log1p(counts / total counts x scale) (scale: project.json, else typed when opening, else the median
total counts per cell). Log-normalised values that are given are stored and shown as they are; without raw counts the raw-counts layer
is not offered. With both, the given values are shown and the counts feed the raw-counts layer and the QC numbers.

The cache lives in the user's cache directory (override with WEBBSCOPE_CACHE), never inside the project folder, so the input data is never
modified. It is rebuilt automatically when any input, or this format version, changes. The columns chosen when opening are remembered there.
"""
import collections, gzip, hashlib, json, os, shutil, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.io import mmread

FORMAT = 3                               # older caches stay valid: values = raw counts when cache_meta.json does not say otherwise
TABLES = ["features.tsv.gz", "barcodes.tsv.gz"]                       # + the cell metadata (CELLS) and a matrix (MATRICES)
CELLS = ["cells_metadata.parquet", "cells.parquet"]                   # the cell metadata: either name (cells.parquet as in Xenium output), not both
MATRICES = {"counts": ["counts.mtx.gz", "counts.mtx", "counts.npz", "matrix.mtx.gz", "matrix.mtx", "matrix.npz"],
            "lognorm": ["counts_logcpx.mtx.gz", "counts_logcpx.mtx", "counts_logcpx.npz"]}     # raw counts / log-normalised values: at most one file of each
REQUIRED = ["project.json", "matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz", "cells.parquet", "cell_boundaries.parquet"]   # the files whose changes redo the import
OPTIONAL = [f for names in MATRICES.values() for f in names if f != "matrix.mtx.gz"] + ["cells_metadata.parquet"] + [ "regions.geojson", "embeddings.parquet", "gloms/gloms.json", "gloms/gloms.npz", "signatures.json", "scffpe_composition.json",
            "he/index.json", "he_gloms/index.json", "transcripts/index.json", "transcripts/tx_x.npy", "transcripts/tx_ptr.npy", "he_slide.json", "he_slide_annotations.geojson"]
OUTSIDE = "(outside all outlines)"
ROLE_DEFAULTS = {"cell_id": "cell_id", "x": "x", "y": "y", "sample": "sample", "lib_size": "lib_size"}
ROLES = ["cell_id", "x", "y", "sample", "donor", "cell_type", "cell_state", "region_lvl1", "region_lvl2", "region_lvl3", "lib_size"]
CHOICES = "open_choices.json"            # in the cache folder: the columns and options chosen when the project was last opened


class ImportError_(Exception):
    """Problem with the user's input files; the message is shown to the user."""


def cache_root():
    return Path(os.environ.get("WEBBSCOPE_CACHE") or Path.home() / ".webbscope" / "cache")


def cache_dir_for(project, h5ad=None):
    """The cache folder of a project folder (of one .h5ad in it, when the folder holds several)."""
    project = Path(project).resolve()
    if h5ad:
        f = project / h5ad
        return cache_root() / f"{project.name}-{Path(h5ad).stem}-{hashlib.sha1(str(f).encode()).hexdigest()[:10]}"
    return cache_root() / f"{project.name}-{hashlib.sha1(str(project).encode()).hexdigest()[:10]}"


def _one_of(project, names, what):
    """The file of `names` present in the folder, None if none; an error when more than one is (which to use would be a guess)."""
    found = [Path(project) / n for n in names if (Path(project) / n).is_file()]
    if len(found) > 1:
        raise ImportError_(f"this folder has {len(found)} files with the {what} ({', '.join(f.name for f in found)}); keep one of them")
    return found[0] if found else None


def matrix_files(project):
    """{"counts": path, "lognorm": path} of the matrix files present (None where absent)."""
    return {"counts": _one_of(project, MATRICES["counts"], "raw counts"), "lognorm": _one_of(project, MATRICES["lognorm"], "log-normalised values")}


def cells_file(project):
    """The cell-metadata table of a project folder (cells_metadata.parquet or cells.parquet), None when there is none."""
    return _one_of(project, CELLS, "cell metadata")


def has_tables(project):
    project = Path(project)
    return all((project / f).is_file() for f in TABLES) and cells_file(project) is not None and any(matrix_files(project).values())


def find_h5ad(project, name=None):
    """The .h5ad file of a project folder: the one named, else the only one in the folder; None when there is none."""
    project = Path(project)
    if name:
        f = project / name
        if not f.is_file():
            raise ImportError_(f"not found: {f}")
        return f
    hs = sorted(p for p in project.glob("*.h5ad") if p.is_file() and not p.name.startswith("."))
    if len(hs) > 1:
        raise ImportError_(f"this folder holds {len(hs)} .h5ad files ({', '.join(h.name for h in hs[:5])}): give the path to one of them "
                           "(only that file is read), or put each in its own folder with its cell outlines and other files")
    return hs[0] if hs else None


def source_of(project, cols=None):
    """The .h5ad to read, or None for the table files (the tables win when a folder has both, unless an .h5ad was named)."""
    name = (cols or {}).get("h5ad")
    if name:
        return find_h5ad(project, name)
    if has_tables(project):
        return None
    return find_h5ad(project)


def remembered(project, h5ad=None):
    """The choices made when this project was last opened ({} when none)."""
    f = cache_dir_for(project, h5ad) / CHOICES
    try:
        return json.loads(f.read_text()) if f.is_file() else {}
    except ValueError:
        return {}


def saved_choices(project, h5ad=None):
    """The choices saved in the project's own project.json ("choices", written from the open page); {} when none, or when they are for another .h5ad."""
    try:
        ch = json.loads((Path(project) / "project.json").read_text()).get("choices") or {}
    except (OSError, ValueError):
        return {}
    if h5ad and ch.get("h5ad") and ch["h5ad"] != h5ad:
        return {}
    return dict(ch)


def known_choices(project, h5ad=None, file_only=False):
    """(name of the .h5ad, the choices known for the project: project.json's, overridden by those remembered in the cache).
    file_only: an .h5ad file was given (not its folder), so the folder's project.json is not used."""
    pj = {} if file_only else saved_choices(project, h5ad)
    name = h5ad or pj.get("h5ad")
    mem = {k: v for k, v in remembered(project, name).items() if k != "file_only"}   # file or folder is decided by what is opened now, not by the last time
    return name, {**pj, **mem}


def _h5ad_outline(h5):
    """The obs column of the .h5ad with the cells' WKT outlines, or None."""
    import h5py
    from . import h5ad
    try:
        with h5py.File(h5, "r") as f:
            return h5ad.outline_columns(f)["cell"]
    except (OSError, KeyError, ValueError, ImportError_):
        return None


def _fingerprint(project, roi_object=None, cols=None):
    fp = {}
    if cols:
        fp["cols"] = {k: cols[k] for k in sorted(cols)}                      # the columns and options chosen when opening
    h5 = source_of(project, cols)
    if h5 is not None:
        st = h5.stat()
        fp["h5ad"] = [h5.name, st.st_size, int(st.st_mtime)]
        wkt = _h5ad_outline(h5)
        if wkt:
            fp["h5ad_outlines"] = wkt                                        # a cache made before these were read is redone once
    if roi_object:
        st = Path(roi_object).stat()
        fp["roi_object"] = [str(Path(roi_object).resolve()), st.st_size, int(st.st_mtime)]
    for f in ([] if (cols or {}).get("file_only") else REQUIRED + OPTIONAL):    # an .h5ad given as a file: the folder's files do not matter
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


def _scale(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        raise ImportError_(f"the normalisation scale must be a number, not {x!r}")
    if not v > 0:
        raise ImportError_("the normalisation scale must be above 0")
    return v


def import_project(project, cdir, log=print, roi_object=None, cols=None):
    t0 = time.time()
    project, cdir = Path(project), Path(cdir)
    cols = dict(cols or {})
    asked = dict(cols)                                                        # as given: the fingerprint is taken from these
    h5 = source_of(project, cols)
    if h5 is None and not has_tables(project):
        raise ImportError_("This folder is not a project: it needs one .h5ad file, or a matrix (counts.mtx or matrix.mtx, and/or counts_logcpx.mtx; also .mtx.gz or .npz) with features.tsv.gz, barcodes.tsv.gz and the cell metadata (cells_metadata.parquet or cells.parquet)")
    only = bool(cols.get("file_only")) and h5 is not None              # an .h5ad file was given: read that file alone, nothing else in its folder
    if only:
        log(f"the file {h5.name} was given: only it is read (the files in its folder are not; give the folder to use them)")
    here = (lambda f: False) if only else (lambda f: (project / f).exists())   # a file of the project folder, when the folder is used
    has_pj = here("project.json")
    man = json.loads((project / "project.json").read_text()) if has_pj else {"name": h5.stem if h5 else project.name, "format_version": 2}
    roles = {**ROLE_DEFAULTS, **man.get("roles", {})}
    emb, C, L, src = {}, None, None, {}           # C = raw counts, L = log-normalised values (genes x cells CSR; either may be absent)
    wkt, wkt_col = None, None                     # an .h5ad's WKT outlines, in the order of the cell table

    if h5 is not None:
        import h5py
        from . import h5ad
        log(f"reading {h5.name}")
        try:
            f = h5py.File(h5, "r")
        except OSError as e:
            raise ImportError_(f"cannot read {h5.name} as an .h5ad (HDF5) file: {e}")
        with f:
            cells, idcol, skipped = h5ad.obs(f)
            if skipped:
                log(f"cell table: {len(skipped)} column(s) of a kind that cannot be read were left out: {', '.join(skipped[:8])}")
            avail = h5ad.matrices(f)
            if not avail:
                raise ImportError_(f"{h5.name} has no matrix of cells x genes (X, layers or raw/X)")
            keys = [x["key"] for x in avail]
            if "counts_matrix" in cols or "lognorm_matrix" in cols:
                ck, lk = cols.get("counts_matrix") or None, cols.get("lognorm_matrix") or None
                for k in (ck, lk):
                    if k and k not in keys:
                        raise ImportError_(f"{h5.name} has no matrix {k} (it has {', '.join(keys)})")
            else:
                ck, lk = h5ad.guess(f, avail)
            if not ck and not lk:
                raise ImportError_("choose the matrix with the raw counts, the one with the log-normalised values, or both")
            genes = h5ad.gene_names(f, lk or ck, cols.get("genes") or None)
            if lk:
                log(f"log-normalised values: {lk}"); L = h5ad.read_matrix(f, lk); src["lognorm"] = lk
            if ck:
                log(f"raw counts: {ck}"); C = h5ad.read_matrix(f, ck); src["counts"] = ck
                if lk:                                                       # raw/X may hold more genes than X: match the counts to X's genes
                    li, ci = h5ad.gene_names(f, lk), h5ad.gene_names(f, ck)
                    if li != ci:
                        pos = pd.Index(ci).get_indexer(li)
                        if (pos < 0).any():
                            raise ImportError_(f"{(pos < 0).sum()} genes of {lk} are not in {ck} (for example {[g for g, q in zip(li, pos) if q < 0][:3]})")
                        log(f"{ck}: the {len(li):,} genes of {lk} are used (of {len(ci):,})")
                        C = C[pos]
            emb = h5ad.embeddings(f)
            oc = h5ad.outline_columns(f)
            if here("cell_boundaries.parquet"):
                if oc["all"]:
                    log(f"cell_boundaries.parquet is used; the WKT in obs ({', '.join(oc['all'])}) is not")
            elif oc["cell"]:
                wkt_col, wkt = oc["cell"], h5ad.outlines(f, oc["cell"])
            elif oc["all"]:
                log(f"obs holds only nucleus outlines ({', '.join(oc['all'])}): the cells are drawn as dots")
        roles["cell_id"] = idcol
        src["h5ad"] = h5.name
    else:
        log("reading gene and cell lists")
        genes, ids = _lines(project / "features.tsv.gz"), _lines(project / "barcodes.tsv.gz")
        if len(set(ids)) != len(ids):
            raise ImportError_("barcodes.tsv.gz contains duplicate cell ids")
        fs = matrix_files(project)
        if fs["counts"] is not None:
            C = _read_matrix_file(fs["counts"], genes, ids, log); src["counts"] = fs["counts"].name
        if fs["lognorm"] is not None:
            L = _read_matrix_file(fs["lognorm"], genes, ids, log); src["lognorm"] = fs["lognorm"].name
        cf = cells_file(project)
        log(f"reading the cell metadata ({cf.name})")
        cells = pd.read_parquet(cf)
    if len(set(genes)) != len(genes):
        dup = pd.Index(genes)[pd.Index(genes).duplicated()][:3].tolist()
        raise ImportError_(f"the gene names are not unique (for example {dup})" + ("; choose another gene-name column when opening" if h5 is not None else ""))
    for M, what in ((C, "raw counts"), (L, "log-normalised values")):
        if M is not None and M.nnz and not np.isfinite(M.data).all():
            raise ImportError_(f"the {what} have missing or infinite values")
    if C is not None and C.nnz:
        if not np.allclose(C.data, np.round(C.data)):
            raise ImportError_(f"the raw counts ({src['counts']}) have values that are not whole numbers. If they are log-normalised, "
                               + ("choose them as the log-normalised values when opening." if h5 is not None else "name the file counts_logcpx.mtx."))
        if C.data.min() < 0:
            raise ImportError_(f"the raw counts ({src['counts']}) have negative values")
    if L is not None and L.nnz:
        if L.data.min() < 0:
            raise ImportError_(f"the log-normalised values ({src['lognorm']}) have negative values: give log-normalised values (log1p of normalised counts), not scaled or z-scored ones")
        if np.allclose(L.data[:200_000], np.round(L.data[:200_000])):
            log(f"note: the log-normalised values ({src['lognorm']}) look like whole numbers; are they raw counts?")
    values = "both" if C is not None and L is not None else "counts" if C is not None else "lognorm"
    log({"both": "raw counts and log-normalised values: the given values are shown, the counts feed the raw-counts layer and QC",
         "counts": "raw counts: log-normalised here", "lognorm": "log-normalised values only: used as they are, no raw-counts layer"}[values])

    from .engine.config import find_columns, suggest
    for r in ROLES:                                                           # the columns chosen when opening ("" = none)
        if r in cols:
            roles[r] = cols[r] or None
    sug = suggest(cells.columns)
    if cols.get("sample"):
        pass
    elif roles["sample"] not in cells.columns and sug["sample"]:
        roles["sample"] = sug["sample"]                                       # not set: the column whose name says sample
        log(f"sample column: {roles['sample']} (found by name)")
    if "donor" not in cols and roles.get("donor") not in cells.columns:
        roles["donor"] = sug["donor"]
    for r in ("cell_id", "x", "y"):
        if roles.get(r) not in cells.columns and r not in cols:
            guess = sug.get(r) or next((c for c in (f"spatial_{r}", f"X_spatial_{r}", f"{r}_centroid", f"centroid_{r}", f"center_{r}") if c in cells.columns), None)
            if guess:
                roles[r] = guess
    for r in ("cell_id", "x", "y", "sample"):
        if roles.get(r) not in cells.columns:
            raise ImportError_(f"the cell table has no column '{roles.get(r)}' (role: {r}); choose it when opening the project" + (" or set it under 'roles' in project.json" if has_pj else "")
                               + (f". Columns with 'sample' in the name: {find_columns(cells.columns, 'sample')[:6] or 'none'}" if r == "sample" else ""))
    made = set()                                                              # roles filled in here, not columns of the input
    if cols and ("cell_type" in cols or "cell_state" in cols):              # the viewer needs a cell type and a cell state (violins, signature statistics)
        if roles.get("cell_type") not in cells.columns:
            roles["cell_type"] = "cell_type" if "cell_type" not in cells.columns else "cell_type_none"
            cells[roles["cell_type"]] = "unlabelled"
            made.add("cell_type")
            log("no cell-type column chosen: every cell is 'unlabelled'")
        if roles.get("cell_state") not in cells.columns or not cols.get("cell_state"):
            roles["cell_state"] = roles["cell_type"]
            made.add("cell_state")
            log(f"no cell-state column chosen: the cell type ({roles['cell_type']}) is used as the cell state")
        cols = {k: v for k, v in cols.items() if k not in ("cell_type", "cell_state") or v}
    if cells[roles["cell_id"]].duplicated().any():
        raise ImportError_(f"the cell id column '{roles['cell_id']}' has duplicate values")
    if h5 is not None:
        ids = cells[roles["cell_id"]].astype(str).tolist()                   # the matrix rows are the cell table's rows
        cells[roles["cell_id"]] = ids
    else:
        cells = cells.set_index(roles["cell_id"], drop=False)
        absent = [i for i in ids if i not in cells.index]
        if absent:
            raise ImportError_(f"{len(absent)} cells in the matrix are not in {cells_file(project).name} (first: {absent[0]})")
        cells = cells.loc[ids].reset_index(drop=True)
    for r in ("x", "y"):
        if cells[roles[r]].dtype.kind not in "fiu":
            raise ImportError_(f"the {r} column '{roles[r]}' is not numeric")
    if cells[[roles["x"], roles["y"]]].isna().any().any():
        raise ImportError_("the cell table has missing x/y coordinates")
    # every sample must be one contiguous block of rows (the analyses address a sample as a row range): stable sort by first appearance
    samp = cells[roles["sample"]].astype(str).values
    rank = {s: i for i, s in enumerate(dict.fromkeys(samp))}
    order = np.argsort(np.array([rank[s] for s in samp]), kind="stable")
    if not (order == np.arange(len(order))).all():
        log("cells of one sample were not contiguous; reordering by sample")
        cells = cells.iloc[order].reset_index(drop=True)
        ids = [ids[i] for i in order]
        wkt = wkt[order] if wkt is not None else None
        C = C[:, order] if C is not None else None
        L = L[:, order] if L is not None else None
        emb = {k: (n, a[order]) for k, (n, a) in emb.items()}
    for M in (C, L):
        if M is not None:
            M.sum_duplicates(); M.sort_indices()
            if values == "both":
                M.eliminate_zeros()                              # so that the two can share their positions when they have the same non-zeros
    if C is not None:
        C = C.tocsr(); C.data = np.round(C.data).astype(np.float32)   # float32 like the original cache (counts < 2^24 are exact)
    if L is not None:
        L = L.tocsr(); L.data = L.data.astype(np.float32)
    m = C if C is not None else L                                # the main matrix of the cache: the counts when there are counts
    shared = values == "both" and np.array_equal(C.indptr, L.indptr) and np.array_equal(C.indices, L.indices)

    cnt, binfo = None, None
    if here("cell_boundaries.parquet"):
        log("reading cell boundaries")
        xy, cnt, binfo = read_boundaries(project / "cell_boundaries.parquet", ids, roles["cell_id"], log)
        off = np.r_[0, np.cumsum(cnt)].astype(np.int64)
    elif wkt is not None:
        log(f"reading cell outlines from the obs column {wkt_col}")
        xy, cnt, binfo = boundaries_from_wkt(wkt, ids, wkt_col, log)
        off = np.r_[0, np.cumsum(cnt)].astype(np.int64)
    elif not only:
        log("no cell_boundaries.parquet: the cells are drawn as dots")

    lib_mismatch = None
    if roles.get("lib_size") in cells.columns:
        if C is not None:
            colsum = np.asarray(C.sum(0)).ravel().astype(np.float32)
            lib_mismatch = int((~np.isclose(cells[roles["lib_size"]].values.astype(np.float32), colsum)).sum())
    elif C is not None:
        roles["lib_size"] = "lib_size" if "lib_size" not in cells.columns else "lib_size_matrix"
        cells[roles["lib_size"]] = np.asarray(C.sum(0)).ravel().astype(np.float32)    # library size = total counts over the panel genes
        made.add("lib_size")
    else:
        roles["lib_size"] = None                            # log-normalised values only: no library size unless a column gives it
    norm_scale = None
    if values == "counts":                                  # only then are the values computed here
        if cols.get("norm_scale") not in (None, ""):
            norm_scale = _scale(cols["norm_scale"])
        elif "norm_scale" in man:
            norm_scale = _scale(man["norm_scale"])
        elif has_pj:
            norm_scale = 10000.0
        else:
            lib = cells[roles["lib_size"]].values.astype(np.float64)
            norm_scale = round(max(1.0, float(np.median(lib))), 1)
            log(f"normalisation scale = median total counts per cell = {norm_scale:g}")

    cdir.mkdir(parents=True, exist_ok=True)
    for stale in ("path_mask.npy", "path_regions.json", "paths_by_sample.json", "signatures.json", "scffpe_composition.json", "gloms",
                  "umap_harmony.npy", "pca_harmony.npy", "poly_offsets.npy", "poly_xy.npy", "csr_data_ln.npy", "ln_indptr.npy", "ln_indices.npy", "ln_data.npy"):
        p = cdir / stale
        shutil.rmtree(p) if p.is_dir() else p.unlink(missing_ok=True)
    log("writing cache")
    np.save(cdir / "csr_indptr.npy", m.indptr.astype(np.int64))
    np.save(cdir / "csr_indices.npy", m.indices.astype(np.int32))
    np.save(cdir / "csr_data.npy", m.data.astype(np.float32))
    if values == "both":                                    # the given log-normalised values next to the counts
        if shared:
            np.save(cdir / "csr_data_ln.npy", L.data)
        else:
            np.save(cdir / "ln_indptr.npy", L.indptr.astype(np.int64)); np.save(cdir / "ln_indices.npy", L.indices.astype(np.int32)); np.save(cdir / "ln_data.npy", L.data)
    if cnt is not None:
        np.save(cdir / "poly_offsets.npy", off)
        np.save(cdir / "poly_xy.npy", xy)
    cells.to_parquet(cdir / "cells.parquet", index=False)
    (cdir / "genes.json").write_text(json.dumps(genes))
    mp = {k: v for k, v in man.get("roles", {}).items() if v in cells.columns}
    mp |= {r: roles.get(r) for r in ROLES if roles.get(r) in cells.columns}
    mp |= {r: None for r in ROLES if r in cols and not cols[r]}                # "none" chosen when opening: not guessed from the names later
    mp |= {"lib_size": roles["lib_size"], "sample": roles["sample"]}
    (cdir / "mapping.json").write_text(json.dumps(mp))
    (cdir / "cache_meta.json").write_text(json.dumps({
        "n_genes": len(genes), "n_cells": len(ids), "nnz": int(m.nnz), "sample_col": roles["sample"],
        "orientation": "genes x cells (rows=genes), stored CSR", "values": values, "layers": {"counts": ["counts"], "lognorm": ["lognorm"], "both": ["counts", "lognorm"]}[values],
        "norm_scale": norm_scale, "lognorm_storage": ("shares the counts' positions" if shared else "own positions") if values == "both" else None,
        "note": "'logcpx' = log1p(counts / lib_size * norm_scale), lib_size = total counts over the panel genes." if values == "counts"
        else "the log-normalised values are used as they are given."}))

    extras = {}
    if here("regions.geojson"):
        extras["regions"] = _regions(project, cells, roles["sample"], roles["x"], roles["y"], cdir, log)
    if here("embeddings.parquet"):
        extras["embeddings"] = _embeddings(project, ids, cdir, log)
    elif emb:
        for kind, (k, a) in emb.items():
            np.save(cdir / f"{kind}_harmony.npy", np.asarray(a, np.float32))
        extras["embeddings"] = {kind: f"{k} ({a.shape[1]}d)" for kind, (k, a) in emb.items()}
        if not all("harmony" in k.lower() for k, _ in emb.values()):
            man["embedding_label"] = "h5ad " + " / ".join(k for k, _ in emb.values())
        log("embeddings from the h5ad: " + ", ".join(extras["embeddings"].values()))
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
    elif here("gloms/gloms.json"):
        if not here("gloms/gloms.npz"):
            raise ImportError_("gloms/gloms.json needs gloms/gloms.npz")
        extras["gloms"] = _gloms(project, ids, cells, cdir, log)
    for f in ("signatures.json", "scffpe_composition.json"):
        if here(f):
            shutil.copyfile(project / f, cdir / f); extras[f] = True
    if here("transcripts/index.json"):
        extras["transcripts"] = True
    if here("he_slide.json"):
        extras["he_slide"] = True
    for d in ("he", "he_gloms"):
        if here(f"{d}/index.json"):
            extras[d] = True
    report = {"n_cells": len(ids), "n_genes": len(genes), "n_samples": len(rank), "nnz": int(m.nnz),
              "cells_without_boundary": int((cnt == 0).sum()) if cnt is not None else len(ids), "cells_with_degenerate_boundary": int(((cnt > 0) & (cnt < 3)).sum()) if cnt is not None else 0,
              "boundaries": binfo, "source": h5.name if h5 is not None else "tables", "values": values, "lib_size_from_table_differs_from_matrix_column_sum": lib_mismatch,
              "matrices": {k: str(v) for k, v in src.items()}, "genes_with_zero_counts": int((np.diff(m.indptr) == 0).sum()), "optional_inputs_found": extras, "seconds": round(time.time() - t0, 2)}
    resolved = {"roles": {r: (None if r in made else roles.get(r)) for r in ROLES}, "values": values, "norm_scale": norm_scale, "matrices": {k: str(v) for k, v in src.items()}}
    (cdir / "cache_info.json").write_text(json.dumps({"format": FORMAT, "inputs": _fingerprint(project, roi_object, asked or None), "manifest": man,
                                                      "project": str(project.resolve()), "report": report, "resolved": resolved}, indent=1))
    log(f"imported {len(ids):,} cells x {len(genes):,} genes in {report['seconds']} s")
    return report


def ensure_cache(project, log=print, force=False, roi_object=None, cols=None):
    """Import the project into its cache unless the cache is current. cols = the columns and options chosen when opening; what is not given
    is taken from the choices made the last time this project was opened (so opening again from the command line asks nothing)."""
    project = Path(project)
    h5name, known = known_choices(project, (cols or {}).get("h5ad"), bool((cols or {}).get("file_only")))
    cols = {**known, **(cols or {})} or None
    cdir = cache_dir_for(project, h5name)
    if not force and is_current(project, cdir, roi_object, cols):
        log("using cached import")
        return cdir
    import_project(project, cdir, log, roi_object, cols)
    if cols:
        (cdir / CHOICES).write_text(json.dumps(cols, indent=1))
    return cdir


def save_project_json(project, cdir, roi_object=None):
    """(Not for an .h5ad given as a file: that reads nothing of its folder, so nothing is written there either.)"""
    """Write project.json into the project folder from the choices of the import that was just made (the roles as resolved, the kind of
    values, the scale, and the choices themselves so the next opening asks nothing), then mark the cache current with it."""
    project, cdir = Path(project), Path(cdir)
    info = json.loads((cdir / "cache_info.json").read_text())
    ch = json.loads((cdir / CHOICES).read_text()) if (cdir / CHOICES).is_file() else {}
    res = info["resolved"]
    out = {"name": info["manifest"].get("name") or project.name, "format_version": 2,
           "roles": {k: v for k, v in res["roles"].items() if v}, "values": res["values"]}
    if res.get("norm_scale"):
        out["norm_scale"] = res["norm_scale"]
    out["matrices"] = res.get("matrices", {})
    out["choices"] = ch
    out["saved_by"] = f"WebbScope open page, {time.strftime('%Y-%m-%d %H:%M')}"
    f, tmp = project / "project.json", project / ".project.json.tmp"
    tmp.write_text(json.dumps(out, indent=2))
    os.replace(tmp, f)
    info["inputs"] = _fingerprint(project, roi_object, ch or None)       # the same import: the new project.json does not make the cache stale
    (cdir / "cache_info.json").write_text(json.dumps(info, indent=1))
    return f


GEOM_COLUMNS = ("geometry", "shape", "geom")


def read_boundaries(path, ids, id_role=None, log=print):
    """cell_boundaries.parquet -> (xy float32 (n_vertices, 2) ordered by cell, then as given; vertices per cell (len(ids)); a summary dict).
    Two layouts:
      one row per outline point: cell_id, vertex_x, vertex_y (as Xenium writes it), the points of a cell in drawing order;
      one row per cell: a cell-id column and the outline in a geometry, shape or geom column, as WKB (GeoParquet; the "geo" metadata names
      the column and must say WKB) or as WKT text ("POLYGON ((x y, ...))").
    The viewer draws one outline per cell: of a multipolygon the largest part is used, holes are left out, and other kinds of geometry
    (points, lines) and empty ones give no outline (the cell is drawn as a dot). Rows for cells that are not in the project are ignored."""
    import pyarrow as pa, pyarrow.parquet as pq
    sch = pq.read_schema(path)
    names = sch.names
    n = len(ids)
    if "vertex_x" in names and "vertex_y" in names:
        if "cell_id" not in names:
            raise ImportError_("cell_boundaries.parquet has vertex_x / vertex_y but no cell_id column")
        bd = pd.read_parquet(path, columns=["cell_id", "vertex_x", "vertex_y"])
        pos = pd.Index(ids).get_indexer(bd["cell_id"].astype(str))
        ok = pos >= 0
        bd, pos = bd[ok], pos[ok]
        o = np.argsort(pos, kind="stable")                      # keeps each cell's vertex order
        xy = np.column_stack([bd["vertex_x"].values[o], bd["vertex_y"].values[o]]).astype(np.float32)
        return xy, np.bincount(pos, minlength=n), {"format": "one row per outline point", "rows_for_unknown_cells": int((~ok).sum())}

    import shapely
    geo = None
    if sch.metadata and b"geo" in sch.metadata:
        try:
            geo = json.loads(sch.metadata[b"geo"])
        except ValueError:
            geo = None
    gcol = geo.get("primary_column") if geo and geo.get("primary_column") in names else next((c for c in GEOM_COLUMNS if c in names), None)
    if gcol is None:
        raise ImportError_("cell_boundaries.parquet needs vertex_x and vertex_y columns (one row per outline point) or a geometry, shape or geom "
                           f"column (one row per cell, WKB or WKT); it has {', '.join(names[:8])}")
    enc = str(((geo or {}).get("columns") or {}).get(gcol, {}).get("encoding", "WKB")) if geo else None
    if enc is not None and enc.upper() != "WKB":
        raise ImportError_(f"cell_boundaries.parquet: the geometry is stored as {enc}; WKB (GeoParquet's usual encoding) or WKT text can be read")
    gtype = sch.field(gcol).type
    if pa.types.is_binary(gtype) or pa.types.is_large_binary(gtype):
        kind = "GeoParquet (WKB)" if geo else "WKB"
    elif pa.types.is_string(gtype) or pa.types.is_large_string(gtype):
        kind = "WKT text"
    else:
        raise ImportError_(f"cell_boundaries.parquet: the {gcol} column is {gtype}; WKB (binary) or WKT (text) can be read")
    # the cell-id column: the one of the text columns whose values match the project's cell ids best (cell_id, the id role's name, cell, ... first)
    first = [c for c in dict.fromkeys(["cell_id", id_role, "cell", "barcode", "cell_ID", "CellID"]) if c and c in names and c != gcol]
    cands = first + [c for c in names if c not in first and c != gcol and (pa.types.is_string(sch.field(c).type) or pa.types.is_large_string(sch.field(c).type) or pa.types.is_integer(sch.field(c).type))]
    if not cands:
        raise ImportError_("cell_boundaries.parquet has a geometry column but no cell-id column")
    idx = pd.Index(ids)
    best, best_n = None, -1
    for c in cands[:12]:
        v = pq.read_table(path, columns=[c]).column(0).to_pandas().astype(str).values
        k = int((idx.get_indexer(v) >= 0).sum())             # the whole column: a file may hold many more cells than the project
        if k > best_n:
            best, best_n, cid = c, k, v
        if k == n:
            break
    if best_n <= 0:
        raise ImportError_(f"cell_boundaries.parquet: none of its columns ({', '.join(cands[:6])}) holds the project's cell ids (for example {ids[0]})")
    raw = pq.read_table(path, columns=[gcol]).column(0).to_numpy(zero_copy_only=False)
    present = np.array([r is not None for r in raw], bool)
    geoms = np.full(len(raw), None, dtype=object)
    if kind == "WKT text":
        geoms[present] = shapely.from_wkt(raw[present], on_invalid="ignore")
    else:
        geoms[present] = shapely.from_wkb(raw[present], on_invalid="ignore")
    bad = present & shapely.is_missing(geoms)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise ImportError_(f"cell_boundaries.parquet: {int(bad.sum()):,} outline(s) could not be read as {kind} (first: cell {cid[i]}: {str(raw[i])[:60]!r})")
    pos = idx.get_indexer(cid)
    known = pos >= 0
    dup = int(pd.Series(pos[known]).duplicated().sum())
    if dup:
        raise ImportError_(f"cell_boundaries.parquet: {dup:,} cells have more than one row")
    geoms, pos = geoms[known], pos[known]
    xy, cnt, st = _outline_points(geoms, pos, n)
    info = {"format": f"one row per cell, {kind}", "id_column": best, "geometry_column": gcol, "rows_for_unknown_cells": int((~known).sum()), **st}
    log(f"cell boundaries: {kind}, cell ids from {best}" + _outline_notes(st))
    return xy, cnt, info


def _outline_points(geoms, pos, n):
    """Shapely geometries of the cells at positions pos -> (xy float32 ordered by cell, vertices per cell, counts of what was simplified or left out).
    One outline per cell: of a multipolygon the largest part; holes are left out; points, lines and empty geometries give none."""
    import shapely
    tid = shapely.get_type_id(geoms)                          # 3 polygon, 6 multipolygon, -1 missing
    empty = shapely.is_missing(geoms) | shapely.is_empty(geoms)
    multi = (tid == 6) & ~empty
    other = ~empty & (tid != 3) & (tid != 6)
    if multi.any():                                           # the largest part of each multipolygon
        parts, owner = shapely.get_parts(geoms[multi], return_index=True)
        o = np.lexsort((-shapely.area(parts), owner))
        firsts = o[np.r_[True, owner[o][1:] != owner[o][:-1]]]
        geoms[np.flatnonzero(multi)[owner[firsts]]] = parts[firsts]
    use = ~empty & ~other
    polys, ppos = geoms[use], pos[use]
    holes = int((shapely.get_num_interior_rings(polys) > 0).sum())
    coords, which = shapely.get_coordinates(shapely.get_exterior_ring(polys), return_index=True)
    cellof = ppos[which]
    o = np.argsort(cellof, kind="stable")
    return coords[o].astype(np.float32), np.bincount(cellof, minlength=n), {"empty": int(empty.sum()), "multipolygons_largest_part_used": int(multi.sum()),
                                                                             "with_holes_left_out": holes, "not_polygons": int(other.sum())}


def _outline_notes(st):
    notes = [f"{v:,} {k.replace('_', ' ')}" for k, v in st.items() if v]
    return f" ({'; '.join(notes)})" if notes else ""


def boundaries_from_wkt(text, ids, column, log=print):
    """WKT outlines of an .h5ad's obs column, one entry per cell in the order of ids -> (xy, vertices per cell, a summary dict), as read_boundaries."""
    import shapely
    from .h5ad import NO_TEXT
    n = len(ids)
    present = np.array([isinstance(t, str) and t.strip().lower() not in NO_TEXT for t in text], bool)
    geoms = np.full(n, None, dtype=object)
    geoms[present] = shapely.from_wkt(text[present].astype(object), on_invalid="ignore")
    bad = present & shapely.is_missing(geoms)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise ImportError_(f"the obs column {column}: {int(bad.sum()):,} outline(s) could not be read as WKT (first: cell {ids[i]}: {str(text[i])[:60]!r})")
    xy, cnt, st = _outline_points(geoms, np.arange(n), n)
    log(f"cell boundaries: WKT text in obs column {column}" + _outline_notes(st))
    return xy, cnt, {"format": "one row per cell, WKT text in the .h5ad", "geometry_column": column, **st}


def _read_matrix_file(path, genes, ids, log):
    """A matrix file (.mtx, .mtx.gz, or an .npz written by scipy.sparse.save_npz) as genes x cells CSR; cells x genes is transposed."""
    log(f"reading {path.name}")
    if path.suffix == ".npz":
        try:
            m = sp.load_npz(path)
        except (ValueError, KeyError, OSError) as e:
            raise ImportError_(f"{path.name} is not a sparse matrix written by scipy.sparse.save_npz ({e}); a dense array is not accepted")
    else:
        m = mmread(path)
    m = m.tocsr() if sp.issparse(m) else sp.csr_matrix(m)
    if m.shape == (len(ids), len(genes)) and m.shape != (len(genes), len(ids)):
        log(f"{path.name} is cells x genes; transposing")
        m = m.T.tocsr()
    if m.shape != (len(genes), len(ids)):
        raise ImportError_(f"{path.name} is {m.shape[0]} x {m.shape[1]} but there are {len(genes)} genes and {len(ids)} cells")
    return sp.csr_matrix(m)
