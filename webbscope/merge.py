"""Merge several processed datasets (project folders or .h5ad files) into one project folder, then open it.

Each dataset is first imported on its own, exactly as when it is opened alone (importer.ensure_cache: h5ad, raw counts and / or
log-normalised values, either metadata name, any outline layout), with its own saved choices (project.json or the cache) or the answers of
the merge page's column form. The merged project is then written from those caches:
  genes            the genes every dataset has (in the first dataset's order)
  values           raw counts when every dataset has them (log-normalised together, one scale); the given log-normalised values when every
                   dataset has those too (both kept), or when not every dataset has counts (then they were normalised separately: a warning)
  cell ids         "<dataset>|<cell id>" (the original id in cell_id_in_dataset)
  samples          kept; a sample name used by two datasets becomes "<dataset>_<sample>" in both
  cell types ...   kept as they are (the datasets' labels are taken to be comparable)
  role columns     written under the role's name (cell_id, x, y, sample, donor, cell_type, cell_state, region_lvl1-3, lib_size); every other
                   metadata column is kept, empty for the cells of datasets without it; merged_from = the dataset of each cell
  outlines         every dataset's outlines, one row per outline point (datasets without outlines: dots)
  regions          the datasets' regions.geojson, with the renamed samples
Embeddings (computed per dataset, not comparable), H&E, transcripts and ROI objects are not merged.
App-written Xenium projects are not merged here (their labels are 10x clusters per run): the Xenium batch page merges those.
"""
import gzip, json, re, shutil, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from . import importer, xenium
from .importer import ImportError_, ROLES

MARKER = ".webbscope_merge.json"          # in a merged project folder: what was merged (also tells it from other folders)


def _is_merged(p, finished=True):
    """A merge written here (finished=False: also one that was stopped half-way, which may be replaced)."""
    return xenium.finished(p, MARKER) if finished else (Path(p) / MARKER).is_file()


def _known(folder, h5name=None, file_only=False):
    """Whether a dataset's columns are known without the form: choices saved in its project.json or remembered in the cache, or a
    project.json that names its columns (roles)."""
    if importer.known_choices(folder, h5name, file_only)[1]:
        return True
    if file_only:
        return False
    try:
        return bool((json.loads((Path(folder) / "project.json").read_text()).get("roles") or {}).get("sample"))
    except (OSError, ValueError):
        return False


def discover(parent):
    """The datasets directly inside `parent`: project folders (tables or one .h5ad) and .h5ad files. Returns {parent, datasets: [{name, path,
    kind, known}], skipped: [{name, reason}], merged: [paths of earlier merges]}. known = its columns are known (saved or remembered)."""
    parent = Path(parent).expanduser()
    if not parent.is_dir():
        raise ImportError_(f"not a folder: {parent}")
    out, skipped, merged = [], [], []
    for p in sorted(parent.iterdir(), key=lambda q: q.name.lower()):
        if p.name.startswith("."):
            continue
        if p.is_file() and p.suffix.lower() == ".h5ad":
            name = p.stem
            out.append({"name": name, "path": str(p), "kind": ".h5ad file", "known": _known(parent, p.name, True)})
            continue
        if not p.is_dir():
            continue
        if _is_merged(p):
            merged.append(str(p)); continue
        if _is_merged(p, finished=False):
            skipped.append({"name": p.name, "reason": "a merge that was stopped before it finished: merge again into this folder to replace it"}); continue
        if (p / xenium.MARKER).is_file():
            skipped.append({"name": p.name, "reason": "written by WebbScope's Xenium import (10x cluster labels per run): merge it on the raw Xenium batch page"}); continue
        try:
            tables = importer.has_tables(p)
            h5 = None if tables else importer.find_h5ad(p)
        except ImportError_ as e:
            skipped.append({"name": p.name, "reason": str(e)}); continue
        if not tables and h5 is None:
            skipped.append({"name": p.name, "reason": "not a dataset (no matrix with cell metadata, no .h5ad)"}); continue
        out.append({"name": p.name, "path": str(p), "kind": "tables" if tables else ".h5ad folder", "known": _known(p)})
    seen = {}
    for d in out:                                               # names prefix the cell ids: they must be unique
        k = d["name"]
        seen[k] = seen.get(k, 0) + 1
        if seen[k] > 1:
            d["name"] = f"{k}_{seen[k]}"
    return {"parent": str(parent), "datasets": out, "skipped": skipped, "merged": merged}


def _where(path):
    """(project folder, cols for a file alone) of a dataset path."""
    p = Path(path).expanduser()
    if p.is_file():
        return p.parent, {"h5ad": p.name, "file_only": True}
    return p, {}


FORMAT_ROLES = ("cell_id", "x", "y")      # columns that depend on the format (an .h5ad: the obs index, obsm["spatial"]), not on what they mean


def _columns(folder, extra):
    """Column names of a dataset's cell table (for an .h5ad: obs, its index and the obsm coordinates), without reading the values."""
    try:
        h5 = importer.source_of(folder, extra or None)
        if h5 is not None:
            import h5py
            from . import h5ad
            with h5py.File(h5, "r") as f:
                return set(h5ad.obs_columns(f)[0])
        import pyarrow.parquet as pq
        return set(pq.read_schema(importer.cells_file(folder)).names)
    except Exception:
        return None


def _answers_for(answers, folder, extra, log, name):
    """The form's answers for one dataset. A cell-id or coordinate column the dataset does not have is left to the dataset's own guess (an
    .h5ad keeps its coordinates as spatial_x / spatial_y, a table as x / y); the other columns must exist (the import says which do not)."""
    cols = dict(answers or {})
    have = _columns(folder, extra)
    if have is not None:
        for r in FORMAT_ROLES:
            if cols.get(r) and cols[r] not in have:
                log(f"  {name}: has no column {cols[r]!r} ({r}); its own is guessed")
                cols.pop(r)
    return cols


def _load(name, cdir):
    """What the merge needs from one dataset's cache."""
    c = Path(cdir)
    meta = json.loads((c / "cache_meta.json").read_text())
    genes = json.loads((c / "genes.json").read_text())
    mp = json.loads((c / "mapping.json").read_text())
    cells = pd.read_parquet(c / "cells.parquet")
    n = len(cells)
    ptr, ind, dat = (np.load(c / f"csr_{k}.npy") for k in ("indptr", "indices", "data"))
    main = sp.csr_matrix((dat, ind, ptr), shape=(len(genes), n))
    values = meta.get("values", "counts")
    counts = main if values in ("counts", "both") else None
    if values == "lognorm":
        ln = main
    elif values == "both":
        ln = sp.csr_matrix((np.load(c / "csr_data_ln.npy"), ind, ptr), shape=(len(genes), n)) if (c / "csr_data_ln.npy").is_file() else \
            sp.csr_matrix((np.load(c / "ln_data.npy"), np.load(c / "ln_indices.npy"), np.load(c / "ln_indptr.npy")), shape=(len(genes), n))
    else:
        ln = None
    poly = (np.load(c / "poly_offsets.npy"), np.load(c / "poly_xy.npy")) if (c / "poly_offsets.npy").is_file() else None
    return {"name": name, "genes": genes, "map": mp, "cells": cells, "counts": counts, "ln": ln, "poly": poly,
            "norm_scale": meta.get("norm_scale"), "values": values}


def _table(d, rename):
    """One dataset's cell table with the role columns under the role names, the other columns as they are (one clashing with a role name gets
    _in_data), cell ids prefixed, samples renamed."""
    t, mp = d["cells"], d["map"]
    src = {r: mp.get(r) for r in ROLES if mp.get(r) in t.columns}
    out = {}
    used = set(src.values())
    for c in t.columns:
        if c in used:
            continue
        out[c + "_in_data" if c in ROLES or c in ("merged_from", "cell_id_in_dataset") else c] = t[c].values
    for r, c in src.items():
        out[r] = t[c].values
    out = pd.DataFrame(out)
    out["cell_id_in_dataset"] = out["cell_id"].astype(str)
    out["cell_id"] = d["name"] + "|" + out["cell_id_in_dataset"]
    out["sample"] = out["sample"].astype(str).map(lambda s: rename.get(s, s))
    out["merged_from"] = d["name"]
    return out


def _clean(df):
    """Columns that pandas could not keep in one type after concatenating (text in one dataset, numbers in another) become text."""
    for c in df.columns:
        if df[c].dtype == object:
            v = df[c]
            kinds = {type(x) for x in v.dropna().values[:20000]}
            if len(kinds) > 1 or (kinds and not kinds <= {str}):
                df[c] = v.map(lambda x: x if x is None or (isinstance(x, float) and np.isnan(x)) else str(x))
    return df


def merge(datasets, answers, out, log=print):
    """datasets: [{name, path}] (from discover); answers: the column form's cols, for the datasets whose columns are not known; out: the folder
    to write (an earlier merge there is replaced). Returns out."""
    t0 = time.time()
    out = Path(out).expanduser()
    if len(datasets) < 2:
        raise ImportError_("choose at least two datasets to merge")
    if out.exists() and any(out.iterdir()) and not _is_merged(out, finished=False):
        raise ImportError_(f"{out} exists and is not an earlier merge; choose a new or empty output folder. {xenium.STOPPED}")
    names = [d["name"] for d in datasets]
    if len(set(names)) != len(names):
        raise ImportError_("two datasets have the same name")
    loaded, errors = [], []
    for i, d in enumerate(datasets):
        folder, extra = _where(d["path"])
        known = _known(folder, extra.get("h5ad"), bool(extra.get("file_only")))
        cols = {**extra} if known else {**_answers_for(answers, folder, extra, log, d["name"]), **extra}
        log(f"[{i + 1}/{len(datasets)}] {d['name']}: importing" + (" (its saved choices)" if known else " (the form's choices)"))
        try:
            cdir = importer.ensure_cache(folder, lambda m, n=d["name"]: log(f"  {n}: {m}"), cols=cols or None)
            loaded.append(_load(d["name"], cdir))
        except ImportError_ as e:
            errors.append(f"{d['name']}: {e}")
    if errors:
        raise ImportError_("these datasets could not be read, so nothing was merged (untick them, or open each alone first to choose its "
                           "columns, which are then remembered):\n" + "\n".join(errors))

    # genes and values
    sets = [set(d["genes"]) for d in loaded]
    common = [g for g in loaded[0]["genes"] if all(g in s for s in sets)]
    if not common:
        raise ImportError_("the datasets have no gene in common")
    dropped = {d["name"]: len(d["genes"]) - len(common) for d in loaded if len(d["genes"]) != len(common)}
    log(f"{len(common):,} genes shared by all datasets" + (f"; left out per dataset: {dropped}" if dropped else " (the gene lists are the same)"))
    all_c, all_l = all(d["counts"] is not None for d in loaded), all(d["ln"] is not None for d in loaded)
    if all_c and all_l:
        values = "both"
    elif all_c:
        values = "counts"
        if any(d["ln"] is not None for d in loaded):
            log("only some datasets have log-normalised values: they are not used; the raw counts of all datasets are log-normalised together")
    elif all_l:
        values = "lognorm"
        none = not any(d["counts"] is not None for d in loaded)
        log(f"WARNING: {'none of the datasets has' if none else 'not every dataset has'} raw counts, so the given log-normalised values are merged; "
            "each dataset was normalised on its own, so the values may not be comparable between datasets")
    else:
        lack = [d["name"] for d in loaded if d["counts"] is None]
        raise ImportError_(f"cannot merge: {', '.join(lack)} {'has' if len(lack) == 1 else 'have'} no raw counts and some other datasets have no "
                           "log-normalised values, so there is no kind of value that all of them have")

    def stack(key):
        return sp.hstack([d[key][pd.Index(d["genes"]).get_indexer(common)] if d["genes"] != common else d[key] for d in loaded], format="csr")

    # samples: a name used by two datasets is prefixed with the dataset in both
    owners = {}
    for d in loaded:
        for s in dict.fromkeys(d["cells"][d["map"]["sample"]].astype(str)):
            owners.setdefault(s, []).append(d["name"])
    clash = {s for s, o in owners.items() if len(o) > 1}
    if clash:
        log(f"{len(clash)} sample name(s) are used by more than one dataset and get the dataset's name in front: {sorted(clash)[:6]}")
    tabs = [_table(d, {s: f"{d['name']}_{s}" for s in clash}) for d in loaded]
    cells = _clean(pd.concat(tabs, ignore_index=True, join="outer"))
    for r in ("x", "y"):
        cells[r] = cells[r].astype(np.float64)
    ids = cells["cell_id"].tolist()
    if len(set(ids)) != len(ids):
        raise ImportError_("duplicate cell ids after merging")
    roles = {r: r for r in ROLES if r in cells.columns and cells[r].notna().any()}
    for r in ("cell_type", "cell_state"):                         # filled in each dataset's import, so always present; keep it that way
        roles.setdefault(r, r)
    if values in ("counts", "both"):
        roles.pop("lib_size", None)                              # recomputed from the merged counts on import
        cells = cells.drop(columns=[c for c in ("lib_size",) if c in cells.columns])

    # write
    out = xenium.claim_folder(out, MARKER, "merge")
    log("writing the merged matrix")
    if values in ("counts", "both"):
        sp.save_npz(out / "counts.npz", stack("counts"), compressed=False)
    if values in ("lognorm", "both"):
        sp.save_npz(out / "counts_logcpx.npz", stack("ln"), compressed=False)
    with gzip.open(out / "features.tsv.gz", "wt") as f:
        f.write("\n".join(common) + "\n")
    with gzip.open(out / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(ids) + "\n")
    log("writing the merged cell metadata and outlines")
    cells.to_parquet(out / "cells_metadata.parquet", index=False)
    parts, start = [], 0
    for d, t in zip(loaded, tabs):
        if d["poly"] is not None:
            off, xy = d["poly"]
            cnt = np.diff(off)
            parts.append(pd.DataFrame({"cell_id": np.repeat(t["cell_id"].values, cnt), "vertex_x": xy[:, 0], "vertex_y": xy[:, 1]}))
    if parts:
        pd.concat(parts, ignore_index=True).to_parquet(out / "cell_boundaries.parquet", index=False)
    n_reg = _regions(datasets, loaded, clash, out, log)
    norm = None
    if values == "counts":
        scales = {d["norm_scale"] for d in loaded}
        if len(scales) == 1 and None not in scales:
            norm = float(scales.pop())                             # the datasets used the same scale: keep it
        else:
            lib = np.asarray(stack("counts").sum(0)).ravel()
            norm = round(max(1.0, float(np.median(lib))), 1)
            log(f"normalisation scale = median total counts per cell over all datasets = {norm:g}")
    choices = {r: roles.get(r, "") for r in ROLES if r != "lib_size"} | {"lib_size": roles.get("lib_size", "")}
    if norm is not None:
        choices["norm_scale"] = str(norm)
    info = [{"name": d["name"], "path": str(x["path"]), "cells": int(len(d["cells"])), "values": d["values"]} for d, x in zip(loaded, datasets)]
    man = {"name": f"{len(loaded)} datasets merged", "format_version": 2, "source": "merged-datasets", "roles": roles, "values": values,
           **({"norm_scale": norm} if norm is not None else {}), "coordinate_unit": "um", "datasets": info, "choices": choices,
           "saved_by": f"WebbScope merge page, {time.strftime('%Y-%m-%d %H:%M')}"}
    (out / "project.json").write_text(json.dumps(man, indent=2))
    (out / MARKER).write_text(json.dumps({"datasets": info, "n_cells": len(ids), "n_genes": len(common), "genes_left_out": dropped,
                                          "values": values, "samples_renamed": sorted(clash), "regions": n_reg}, indent=1))
    log(f"merged {len(loaded)} datasets: {len(ids):,} cells, {len(common):,} genes, {cells['sample'].nunique()} samples, values: {values} "
        f"({time.time() - t0:.0f} s)")
    return out


def _regions(datasets, loaded, clash, out, log):
    """The datasets' regions.geojson (folders only: a file alone has none) with the renamed samples, into one file. Returns the feature count."""
    feats = []
    for x, d in zip(datasets, loaded):
        p = Path(x["path"])
        f = p / "regions.geojson" if p.is_dir() else None
        if f is None or not f.is_file():
            continue
        try:
            for ft in json.loads(f.read_text()).get("features", []):
                pr = dict(ft.get("properties") or {})
                s = str(pr.get("sample"))
                pr["sample"] = f"{d['name']}_{s}" if s in clash else s
                feats.append({**ft, "properties": pr})
        except (OSError, ValueError) as e:
            log(f"{d['name']}: regions.geojson ignored ({e})")
    if not feats:
        return 0
    kinds = {(ft.get("properties") or {}).get("annotation_name") for ft in feats}
    if len(kinds) > 15:
        log(f"regions not merged: together the datasets have {len(kinds)} annotation types (at most 15)")
        return 0
    (out / "regions.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    log(f"regions: {len(feats)} outlines from the datasets' regions.geojson")
    return len(feats)
