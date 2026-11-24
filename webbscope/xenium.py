"""Xenium Onboard Analysis output -> WebbScope project folder.

Input  : the .tar.gz or .zip of a Xenium output bundle, or its (extracted) folder, or the individual files. A folder may keep cell_feature_matrix and analysis
         packed as cell_feature_matrix.tar.gz / analysis.tar.gz (older output): the needed files are then extracted to the staging folder.
Needed : cell_feature_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}, cells.parquet (or cells.csv.gz),
         cell_boundaries.parquet (or cell_boundaries.csv.gz)
Optional: transcripts.parquet (decoded transcripts, only if asked for), experiment.xenium (run / panel metadata), analysis/clustering/gene_expression_graphclust/clusters.csv (cell labels),
         analysis/umap/gene_expression_2_components/projection.csv, analysis/pca/gene_expression_*_components/projection.csv
Optional: morphology_focus.ome.tif (else morphology_mip.ome.tif): the IF image, only if asked for; it is not copied, the project records where it is (if_image.py)
Never read: the z-stack morphology.ome.tif, per-cycle and background images, the *.zarr.zip copies, nucleus boundaries, differential expression.

From a .tar.gz or .zip only the needed members are written (to a staging folder, never into the tarball's folder) and reading stops as soon as all of
them were seen. One slide can hold several tissue sections; they are split by cell density (editable in the preview) and become the
project's samples. Nothing is written to the input; the project folder goes to the output folder the user chooses.
"""
import gzip, hashlib, json, math, os, re, shutil, tarfile, time, zipfile
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import ndimage
from scipy.io import mmread, mmwrite
from . import importer, if_image

MARKER = ".webbscope_xenium_import.json"
STOPPED = ("If it is a WebbScope import that was stopped half-way by an earlier version (one that did not mark its folders until the end), "
           "delete that folder or choose another.")


def claim_folder(out, marker=MARKER, label="WebbScope import"):
    """Empty `out` (when it holds an earlier import of this kind) and mark it as WebbScope's before anything is written, so that a job stopped
    half-way leaves a folder the next import may replace. The marker says "incomplete" until the import rewrites it with its details at the end.
    Refuses a folder that is not empty and has no marker. The marker is set to "incomplete" before anything is deleted and is never deleted itself,
    so a job stopped while emptying leaves a folder that is known to be unfinished."""
    out = Path(out).expanduser()
    if out.exists() and not (out / marker).exists() and any(out.iterdir()):
        raise importer.ImportError_(f"{out} exists and is not an earlier {label}; choose a new or empty folder. {STOPPED}")
    out.mkdir(parents=True, exist_ok=True)
    (out / marker).write_text(json.dumps({"incomplete": True, "started": time.strftime("%Y-%m-%d %H:%M:%S")}))
    for f in out.iterdir():
        if f.name != marker:
            shutil.rmtree(f) if f.is_dir() and not f.is_symlink() else f.unlink()
    return out


def finished(folder, marker=MARKER):
    """True when `folder` holds an import of this kind that ran to the end."""
    f = Path(folder) / marker
    try:
        return f.is_file() and not json.loads(f.read_text()).get("incomplete")
    except (OSError, ValueError, AttributeError):
        return False
IF_KEY = "run"                                  # name of the IF image in a single run's if.json
CTRL_KEYS = {"Negative Control Probe": "nc_probe", "Negative Control Codeword": "nc_codeword", "Unassigned Codeword": "unassigned_codeword", "Deprecated Codeword": "deprecated_codeword"}
BIN_UM = 50.0                                   # density grid for tissue detection from cells
TX_BIN_UM = 100.0                               # ... and from transcripts
CLUSTERS_REL = "analysis/clustering/gene_expression_graphclust/clusters.csv"
UMAP_REL = "analysis/umap/gene_expression_2_components/projection.csv"
PCA_RE = re.compile(r"^analysis/pca/gene_expression_\d+_components/projection\.csv$")
FIXED = {"experiment.xenium": "experiment", "cells.parquet": "cells_parquet", "cells.csv.gz": "cells_csv",
         "cell_boundaries.parquet": "bnd_parquet", "cell_boundaries.csv.gz": "bnd_csv",
         "cell_feature_matrix/matrix.mtx.gz": "mtx", "cell_feature_matrix/features.tsv.gz": "features",
         "cell_feature_matrix/barcodes.tsv.gz": "barcodes", CLUSTERS_REL: "clusters", UMAP_REL: "umap", "transcripts.parquet": "tx"}
REQUIRED_GROUPS = [("experiment",), ("cells_parquet", "cells_csv"), ("bnd_parquet", "bnd_csv"), ("mtx",), ("features",), ("barcodes",)]
OPTIONAL_ROLES = ("clusters", "umap", "pca")


def _classify(name):
    """Role of a tar member / relative path, or None. Accepts a leading run folder (output-XETG...)."""
    cands = [name] + ([name.split("/", 1)[1]] if "/" in name else [])
    for rel in cands:
        if rel in FIXED:
            return FIXED[rel], rel
        if PCA_RE.match(rel):
            return "pca", rel
        r = if_image.role_of(rel)
        if r:
            return r, rel
    return None, None


def staging_dir(tar_path):
    st = Path(tar_path).stat()
    h = hashlib.sha1(f"{Path(tar_path).resolve()}|{st.st_size}|{int(st.st_mtime)}".encode()).hexdigest()[:12]
    return importer.cache_root().parent / "staging" / h


def is_archive(name):
    """A Xenium bundle archive by its name. The *.zarr.zip files inside an output are 10x's zarr copies, not bundles."""
    n = str(name).lower()
    return n.endswith((".tar.gz", ".tgz", ".tar")) or (n.endswith(".zip") and not n.endswith(".zarr.zip"))


NESTED = {"cell_feature_matrix": ("mtx", "features", "barcodes"), "analysis": OPTIONAL_ROLES}      # parts that older output keeps packed in their own archive


def _sniff(path):
    """What the file really is, from its first bytes: 'zip', 'gzip', 'tar' or None. The extension can be wrong (a .tar.gz saved as .zip, a download that did not finish)."""
    with open(path, "rb") as f:
        head = f.read(512)
    if head[:4] in (b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08"):
        return "zip"
    if head[:2] == b"\x1f\x8b":
        return "gzip"
    if head[257:262] == b"ustar":
        return "tar"
    return None


def extract_archive(path, log, want_tx=False, want_if=False):
    """The needed files of a .tar.gz / .tgz / .tar / .zip bundle -> {role: Path} in the staging folder. The format is taken from the file's content, not its name."""
    path = Path(path)
    kind = _sniff(path)
    by_name = "zip" if path.name.lower().endswith(".zip") else "gzip"
    if kind is None:
        with open(path, "rb") as f:
            first = f.read(8)
        raise importer.ImportError_(f"{path.name} ({path.stat().st_size:,} bytes) is not a zip, tar or tar.gz archive (its first bytes are {first!r}). "
                                    "A download that did not finish, or a file saved under the wrong name, are the usual causes: compare its size with the download page.")
    if kind != by_name and not (kind == "tar" and by_name == "gzip"):
        log(f"{path.name} is named like a {by_name} file but is a {kind} archive; reading it as one")
    return _extract_zip(path, log, want_tx, want_if) if kind == "zip" else extract_tar(path, log, want_tx, want_if)


def _extract_zip(zip_path, log, want_tx=False, want_if=False):
    """Like extract_tar for a .zip (members are read directly, nothing is streamed through). A cell_feature_matrix.tar.gz / analysis.tar.gz inside
    the zip (older 10x output keeps them packed) is unpacked too, when the zip has no unpacked copy of those files."""
    zip_path = Path(zip_path)
    stage = staging_dir(zip_path)
    done = stage / "complete.json"
    if done.exists():
        found = {k: (stage / v).resolve() for k, v in json.loads(done.read_text()).items() if (stage / v).exists()}
        if found and (not want_tx or "tx" in found) and (not want_if or any(k.startswith("if_") for k in found)):
            log(f"using files already extracted from this archive ({stage})")
            return found
    stage.mkdir(parents=True, exist_ok=True)
    found, nested = {}, {}
    log(f"reading {zip_path.name} ({zip_path.stat().st_size / 1e9:.1f} GB): only the needed files are extracted")
    try:
        zf_cm = zipfile.ZipFile(zip_path)
    except zipfile.BadZipFile as e:
        raise importer.ImportError_(f"{zip_path.name} ({zip_path.stat().st_size:,} bytes) starts like a zip file but cannot be read as one ({e}): its end-of-archive record is missing or unreadable. "
                                    "Usual causes: the download did not finish (compare the size with the download page, or run `unzip -t` on it), or a zip over 4 GB made by macOS Finder "
                                    "(extract it with `ditto -x -k FILE.zip FOLDER` and give the folder instead).")
    with zf_cm as zf:
        for m in zf.infolist():
            if m.is_dir():
                continue
            role, rel = _classify(m.filename)
            if role is None:
                base = m.filename.rsplit("/", 1)[-1]
                for key in NESTED:
                    if base in (key + ".tar.gz", key + ".tgz", key + ".tar", key + ".zip"):
                        nested[key] = m
                continue
            if role in found or (role == "tx" and not want_tx) or (role.startswith("if_") and not want_if):
                continue
            dest = stage / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            log(f"  extracting {rel} ({m.file_size / 1e6:.0f} MB)")
            with zf.open(m) as src, open(str(dest) + ".part", "wb") as out:
                shutil.copyfileobj(src, out, 1 << 20)
            Path(str(dest) + ".part").rename(dest)
            found[role] = dest
        for key, m in nested.items():                                       # a packed part, only when its files were not found unpacked
            if all(r in found for r in NESTED[key]) or (key == "analysis" and any(r in found for r in NESTED[key])):
                continue
            dest = stage / m.filename.rsplit("/", 1)[-1]
            with zf.open(m) as src, open(str(dest) + ".part", "wb") as out:
                shutil.copyfileobj(src, out, 1 << 20)
            Path(str(dest) + ".part").rename(dest)
            try:
                for r, q in extract_archive(dest, lambda _m: None).items():
                    found.setdefault(r, q)
            except (tarfile.TarError, zipfile.BadZipFile, OSError):
                pass
    done.write_text(json.dumps({k: os.path.relpath(v, stage) for k, v in found.items()}))
    log("extracted: " + ", ".join(sorted(found)))
    return found


def extract_tar(tar_path, log, want_tx=False, want_if=False):
    """Stream the archive once and write only the wanted members; stop early when all were seen. Returns {role: Path}."""
    tar_path = Path(tar_path)
    stage = staging_dir(tar_path)
    done = stage / "complete.json"
    if done.exists():
        found = {k: stage / v for k, v in json.loads(done.read_text()).items() if (stage / v).exists()}
        if found and (not want_tx or "tx" in found) and (not want_if or any(k.startswith("if_") for k in found)):
            log(f"using files already extracted from this archive ({stage})")
            return found
    stage.mkdir(parents=True, exist_ok=True)
    total = tar_path.stat().st_size
    found, last = {}, 0.0
    log(f"reading {tar_path.name} ({total / 1e9:.1f} GB) — only the needed files are extracted; this takes a few minutes for a full run")
    with open(tar_path, "rb") as raw, tarfile.open(fileobj=raw, mode="r|*") as tf:            # gzip, bzip2, xz or none: whatever the file is
        for m in tf:
            if time.time() - last > 5:
                log(f"  {100 * raw.tell() / total:.0f}% of the archive read; extracted: {', '.join(sorted(found)) or 'nothing yet'}")
                last = time.time()
            if not m.isfile():
                continue
            role, rel = _classify(m.name)
            if role is None or role in found or (role == "tx" and not want_tx) or (role.startswith("if_") and not want_if):
                continue
            dest = stage / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            with tf.extractfile(m) as src, open(str(dest) + ".part", "wb") as out:
                shutil.copyfileobj(src, out, 1 << 20)
            Path(str(dest) + ".part").rename(dest)
            found[role] = dest
            if all(any(r in found for r in g) for g in REQUIRED_GROUPS) and all(r in found for r in OPTIONAL_ROLES) and (not want_tx or "tx" in found) and (not want_if or "if_focus" in found):
                log("all needed files found; skipping the rest of the archive")
                break
    done.write_text(json.dumps({k: str(v.relative_to(stage)) for k, v in found.items()}))
    log("extracted: " + ", ".join(sorted(found)))
    return found


def find_in_folder(folder):
    """Locate the Xenium files in an extracted output folder."""
    folder = Path(folder)
    found = {}
    for rel, role in FIXED.items():
        if (folder / rel).exists():
            found[role] = folder / rel
    for rel, role in if_image.ROLE_FILES.items():
        if (folder / rel).exists():
            found[role] = folder / rel
    for p in sorted(list((folder / "morphology_focus").glob("*.ome.tif")) + list(folder.glob("morphology_focus_*.ome.tif"))):
        r = if_image.role_of(p.relative_to(folder).as_posix())
        if r:
            found[r] = p
    d = folder / "analysis" / "pca"
    if d.is_dir():
        for p in sorted(d.glob("gene_expression_*_components/projection.csv")):
            if PCA_RE.match(p.relative_to(folder).as_posix()):
                found["pca"] = p
    for name, role in (("matrix.mtx.gz", "mtx"), ("features.tsv.gz", "features"), ("barcodes.tsv.gz", "barcodes")):     # flat layout
        if role not in found and (folder / name).exists():
            found[role] = folder / name
    for key, roles in NESTED.items():                                          # older output keeps these packed (.tar.gz, or .zip)
        for suffix in (".tar.gz", ".tgz", ".tar", ".zip"):
            t = folder / (key + suffix)
            if t.is_file() and not all(r in found for r in roles):
                try:
                    for r, q in extract_archive(t, lambda m: None).items():      # only the needed members are written, to the staging folder
                        found.setdefault(r, q)
                except (tarfile.TarError, zipfile.BadZipFile, OSError):
                    pass
    return found


def resolve_sources(path, overrides, log, want_tx=False, want_if=False):
    """{role: Path} from a tar.gz / folder (`path`) and explicit per-file overrides (matrix may be a folder)."""
    src = {}
    if path:
        p = Path(path).expanduser()
        if p.is_dir():
            src = find_in_folder(p)
        elif p.is_file() and is_archive(p.name):
            src = extract_archive(p, log, want_tx, want_if)
        elif p.exists():
            raise importer.ImportError_(f"{p} is neither a folder nor a .tar.gz / .zip archive")
        else:
            raise importer.ImportError_(f"not found: {p}")
    if not want_tx:
        tx = src.pop("tx", None)
        if tx is not None:
            src["tx_probe"] = tx                              # not imported, but readable for the transcript-density tissue detection
    if not want_if:
        for k in [k for k in src if k.startswith("if_")]:
            del src[k]
    for k, v in (overrides or {}).items():
        v = (v or "").strip()
        if not v:
            continue
        p = Path(v).expanduser()
        if not p.exists():
            raise importer.ImportError_(f"not found: {p}")
        if k == "matrix":
            d = p if p.is_dir() else p.parent
            for name, role in (("matrix.mtx.gz", "mtx"), ("features.tsv.gz", "features"), ("barcodes.tsv.gz", "barcodes")):
                if (d / name).exists():
                    src[role] = d / name
            if p.is_file():
                src["mtx"] = p
        elif k == "cells":
            src["cells_csv" if p.name.endswith(".csv.gz") or p.name.endswith(".csv") else "cells_parquet"] = p
        elif k == "boundaries":
            src["bnd_csv" if p.name.endswith(".csv.gz") or p.name.endswith(".csv") else "bnd_parquet"] = p
        elif k in ("clusters", "umap", "pca", "experiment"):
            src[k] = p
        elif k == "transcripts":
            src["tx"] = p
        elif k == "if_image":                                 # an .ome.tif, or a Xenium output folder that holds one
            if p.is_dir():
                found = {r: q for r, q in find_in_folder(p).items() if r.startswith("if_")}
                if not found:
                    raise importer.ImportError_(f"no morphology_focus / morphology_mip .ome.tif in {p}")
                src.update(found)
            else:
                src["if_focus"] = p
        else:
            raise importer.ImportError_(f"unknown file role {k!r}")
    missing = []
    for g in REQUIRED_GROUPS:
        if g == ("experiment",):
            continue                                        # metadata only
        if not any(r in src for r in g):
            missing.append({"mtx": "cell_feature_matrix/matrix.mtx.gz", "features": "cell_feature_matrix/features.tsv.gz",
                            "barcodes": "cell_feature_matrix/barcodes.tsv.gz", "cells_parquet": "cells.parquet (or cells.csv.gz)",
                            "bnd_parquet": "cell_boundaries.parquet (or cell_boundaries.csv.gz)"}[g[0]])
    if missing:
        raise importer.ImportError_("Xenium output is missing: " + ", ".join(missing))
    return src


def _read_table(p):
    p = Path(p)
    if p.name.endswith(".parquet"):
        return pd.read_parquet(p)
    return pd.read_csv(p)


def _read_lines(p, col=0):
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt") as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def tx_density(path, bin_um=TX_BIN_UM, min_qv=20.0):
    """Gene transcripts (quality >= min_qv) per square bin of bin_um micrometres, streamed one row group at a time (about 0.4 GB for 70 million
    transcripts, 4 s). Returns (x origin, y origin, 2-D counts [x bin, y bin])."""
    import pyarrow.parquet as pq
    pf = pq.ParquetFile(path)
    names, md = pf.schema_arrow.names, pf.metadata
    for c in ("x_location", "y_location"):
        if c not in names:
            raise importer.ImportError_(f"{Path(path).name} has no column {c}")
    cols = ["x_location", "y_location"] + [c for c in ("qv", "is_gene") if c in names]
    ci = [names.index("x_location"), names.index("y_location")]
    lo, hi = [np.inf, np.inf], [-np.inf, -np.inf]
    for g in range(md.num_row_groups):                        # extent from the row-group statistics (no data read)
        for k, c in enumerate(ci):
            st = md.row_group(g).column(c).statistics
            if st is None or not st.has_min_max:
                lo = None
                break
            lo[k], hi[k] = min(lo[k], st.min), max(hi[k], st.max)
        if lo is None:
            break
    if lo is None:                                            # no statistics: one pass to find the extent
        lo, hi = [np.inf, np.inf], [-np.inf, -np.inf]
        for g in range(md.num_row_groups):
            t = pf.read_row_group(g, columns=cols[:2])
            for k in (0, 1):
                v = t.column(k).to_numpy(zero_copy_only=False)
                if len(v):
                    lo[k], hi[k] = min(lo[k], float(np.nanmin(v))), max(hi[k], float(np.nanmax(v)))
    if not (np.isfinite(lo).all() and np.isfinite(hi).all()):
        raise importer.ImportError_(f"{Path(path).name} has no transcript positions")
    nx, ny = int((hi[0] - lo[0]) // bin_um) + 1, int((hi[1] - lo[1]) // bin_um) + 1
    if nx * ny > 5e7:
        raise importer.ImportError_(f"the transcript extent ({nx} x {ny} bins) is implausibly large")
    H = np.zeros(nx * ny, np.int64)
    for g in range(md.num_row_groups):
        t = pf.read_row_group(g, columns=cols)
        x, y = t.column("x_location").to_numpy(zero_copy_only=False), t.column("y_location").to_numpy(zero_copy_only=False)
        m = np.isfinite(x) & np.isfinite(y)
        if "qv" in cols:
            m &= t.column("qv").to_numpy(zero_copy_only=False) >= min_qv
        if "is_gene" in cols:
            m &= t.column("is_gene").to_numpy(zero_copy_only=False).astype(bool)
        H += np.bincount(((x[m] - lo[0]) // bin_um).astype(np.int64) * ny + ((y[m] - lo[1]) // bin_um).astype(np.int64), minlength=nx * ny)
    return float(lo[0]), float(lo[1]), H.reshape(nx, ny)


def auto_tx_threshold(H):
    """Minimum transcripts per bin that counts as tissue when the user gives none: 0.5% of the median count of the busiest quarter of the occupied bins
    (about 120 per 100 µm bin on a 5000-gene panel, where tissue bins hold thousands and cell-free background a handful), never below 5."""
    nz = H[H > 0]
    if not len(nz):
        return 5
    top = nz[nz >= np.percentile(nz, 75)]
    return int(max(5, round(0.005 * float(np.median(top)))))


class Session:
    """A scanned Xenium output: counts of the panel genes, the cell table and optional analysis, all in memory."""

    def __init__(self, src, log):
        self.src, self.log = src, log
        self.meta = {}
        if "experiment" in src:
            try:
                self.meta = json.loads(Path(src["experiment"]).read_text())
            except Exception:
                self.meta = {}
        log("reading gene and cell lists")
        feats = _read_lines(src["features"])
        ids = [r[0] for r in _read_lines(src["barcodes"])]
        if len(set(ids)) != len(ids):
            raise importer.ImportError_("barcodes.tsv.gz contains duplicate cell ids")
        has_type = all(len(r) >= 3 for r in feats)
        keep = np.array([(r[2] == "Gene Expression") if has_type else True for r in feats])
        self.n_features, self.n_genes = len(feats), int(keep.sum())
        names = [r[1] if len(r) > 1 else r[0] for r in feats]
        self.feature_types = pd.Series([r[2] if has_type else "Gene Expression" for r in feats]).value_counts().to_dict()
        log(f"features: {len(feats)} rows, {self.n_genes} are genes ({', '.join(f'{k}: {v}' for k, v in self.feature_types.items())})")
        log("reading counts matrix")
        m = mmread(src["mtx"])
        m = m.tocsr() if sp.issparse(m) else sp.csr_matrix(m)
        if m.shape == (len(ids), len(feats)) and m.shape != (len(feats), len(ids)):
            m = m.T.tocsr()
        if m.shape != (len(feats), len(ids)):
            raise importer.ImportError_(f"matrix is {m.shape[0]} x {m.shape[1]} but there are {len(feats)} features and {len(ids)} cells")
        if m.nnz and not np.allclose(m.data, np.round(m.data)):
            raise importer.ImportError_("counts matrix has non-integer values (expected raw counts)")
        self.ctrl, self.ctrl_n = {}, {}                      # non-gene features (negative control probes / codewords, unassigned, deprecated): counts per cell, kept for the QC report
        if has_type:
            types = np.array([r[2] for r in feats])
            for t in pd.unique(types[~keep]):
                rows = np.flatnonzero(types == t)
                self.ctrl[str(t)] = np.asarray(m[rows].sum(0)).ravel().astype(np.int32)
                self.ctrl_n[str(t)] = int(len(rows))
        gi = np.flatnonzero(keep)
        m = m[gi].tocsr()
        genes = [names[i] for i in gi]
        ens = [feats[i][0] for i in gi]
        dup = pd.Series(genes).duplicated(keep=False).values
        if dup.any():
            genes = [f"{g}_{e}" if d else g for g, e, d in zip(genes, ens, dup)]
            log(f"{int(dup.sum())} features share a gene symbol; their Ensembl id was appended")
        self.genes, self.m, self.ids = genes, m, ids
        self.gene_counts = np.asarray(m.sum(0)).ravel().astype(np.int64)

        log("reading cell table")
        key = "cells_parquet" if "cells_parquet" in src else "cells_csv"
        c = _read_table(src[key])
        c["cell_id"] = c["cell_id"].astype(str)
        for col in ("x_centroid", "y_centroid"):
            if col not in c.columns:
                raise importer.ImportError_(f"cells table has no column {col}")
        c = c.drop_duplicates("cell_id").set_index("cell_id", drop=False)
        absent = [i for i in ids if i not in c.index]
        if absent:
            raise importer.ImportError_(f"{len(absent)} cells of the matrix are not in the cell table (first: {absent[0]})")
        self.cells = c.loc[ids].reset_index(drop=True)
        self.x = self.cells["x_centroid"].values.astype(np.float64)
        self.y = self.cells["y_centroid"].values.astype(np.float64)
        self.n = len(ids)
        self.clusters = None
        if "clusters" in src:
            cl = pd.read_csv(src["clusters"])
            self.clusters = pd.Series(cl.iloc[:, 1].values, index=cl.iloc[:, 0].astype(str).values)
            self.clusters = self.clusters[~self.clusters.index.duplicated()].reindex(ids)
            log(f"graph clusters: {int(self.clusters.notna().sum()):,} cells in {self.clusters.nunique()} clusters")
        self.emb = {}
        if "umap" in src:
            u = pd.read_csv(src["umap"]); u = u.set_index(u.columns[0])
            self.emb["umap"] = u.iloc[:, :2].reindex(ids)
        if "pca" in src:
            p = pd.read_csv(src["pca"]); p = p.set_index(p.columns[0])
            self.emb["pca"] = p.reindex(ids)
        self._auto, self._txview = {}, {}
        self.xy32 = np.column_stack([self.x, self.y]).astype(np.float32)
        log(f"{self.n:,} cells x {self.n_genes:,} genes ready")

    # ---- tissue detection --------------------------------------------------------------------------------------------------------
    def tx_grid(self):
        """Transcript counts per bin of the run (read once, then kept): (x origin, y origin, counts)."""
        if getattr(self, "_txg", None) is None:
            p = self.src.get("tx") or self.src.get("tx_probe")
            if not p:
                raise importer.ImportError_("tissue detection from transcripts needs transcripts.parquet, and this run has none that could be found. For a .tar.gz or .zip "
                                            "bundle the transcripts have to be extracted: tick the transcripts option and read the run again.")
            self._txg = tx_density(p)
        return self._txg

    def piece_boxes(self, gap, min_cells, min_counts, detect, tx_min):
        """{piece id: [x0, y0, x1, y1]} of the detected tissues: the extent of the piece's cells, or (transcript detection) of its transcript bins."""
        comp, comps = self.auto(gap, min_cells, min_counts, detect, tx_min)
        tv = self._txview.get((round(float(gap), 3), int(min_cells), int(min_counts), "transcripts", int(tx_min))) if detect == "transcripts" else None
        out = {}
        for c in comps:
            if c["id"] <= 0:
                continue
            if tv is not None and c["id"] in tv["boxes"]:
                out[c["id"]] = [round(float(v), 1) for v in tv["boxes"][c["id"]]]
            else:
                m = comp == c["id"]
                if m.any():
                    out[c["id"]] = [round(float(self.x[m].min()), 1), round(float(self.y[m].min()), 1), round(float(self.x[m].max()), 1), round(float(self.y[m].max()), 1)]
        return out

    def tx_payload(self, gap, min_cells, min_counts, detect, tx_min):
        """Transcript density for display (None in cell mode): the bins that hold transcripts as flat [x bin, y bin, brightness 0-255 (log scale), tissue id] quadruples."""
        self.auto(gap, min_cells, min_counts, detect, tx_min)
        v = self._txview.get((round(float(gap), 3), int(min_cells), int(min_counts), "transcripts", int(tx_min))) if detect == "transcripts" else None
        if v is None:
            return None
        H = v["counts"]
        xs, ys = np.nonzero(H)
        c = H[xs, ys].astype(np.float64)
        top = max(float(np.percentile(c, 99.5)), 1.0)
        bright = np.clip(np.round(255 * np.log1p(c) / np.log1p(top)), 1, 255).astype(np.int64)
        return {"ox": v["ox"], "oy": v["oy"], "bin": v["bin"], "thr": v["thr"], "q": np.column_stack([xs, ys, bright, v["lab"][xs, ys]]).astype(np.int64).ravel().tolist()}

    def auto(self, gap, min_cells, min_counts, detect="cells", tx_min=0):
        """Connected components of the cell-density map. Returns (comp id per cell: -1 = below min_counts, 0 = isolated cell,
        1..K = pieces in the order they would be named, sorted by x), and a list of piece descriptions."""
        detect = "transcripts" if detect == "transcripts" else "cells"
        key = (round(float(gap), 3), int(min_cells), int(min_counts), detect, int(tx_min) if detect == "transcripts" else 0)
        if key in self._auto:
            return self._auto[key]
        ok = self.gene_counts >= min_counts
        if detect == "transcripts":                          # tissue = bins with enough transcripts, wherever the cells are
            ox, oy, Ht = self.tx_grid()
            nx, ny = Ht.shape
            ix = np.clip(((self.x - ox) // TX_BIN_UM).astype(np.int64), 0, nx - 1)
            iy = np.clip(((self.y - oy) // TX_BIN_UM).astype(np.int64), 0, ny - 1)
            thr = int(tx_min) if int(tx_min) > 0 else auto_tx_threshold(Ht)
            self.tx_threshold = thr
            occ = Ht >= thr
        else:
            ix = ((self.x - self.x.min()) // BIN_UM).astype(np.int64)
            iy = ((self.y - self.y.min()) // BIN_UM).astype(np.int64)
            nx, ny = ix.max() + 1, iy.max() + 1
            H = np.bincount((ix * ny + iy)[ok], minlength=nx * ny).reshape(nx, ny)
            occ = H >= 2
        B = TX_BIN_UM if detect == "transcripts" else BIN_UM
        it = max(0, int(math.ceil(gap / (2 * B))))
        dil = ndimage.binary_dilation(occ, iterations=it) if it else occ
        lab, k = ndimage.label(dil, structure=np.ones((3, 3), int))
        cell_lab = lab[ix, iy]
        cell_lab = np.where(occ[ix, iy] | (cell_lab > 0), cell_lab, 0)
        cnt = np.bincount(cell_lab[ok], minlength=k + 1)
        cx = np.array([self.x[ok & (cell_lab == i)].mean() if cnt[i] else 0.0 for i in range(k + 1)])
        big = [i for i in range(1, k + 1) if cnt[i] >= min_cells]
        small = [i for i in range(1, k + 1) if 0 < cnt[i] < min_cells]
        order = sorted(big, key=lambda i: cx[i]) + sorted(small, key=lambda i: -cnt[i])
        remap = np.zeros(k + 1, np.int64)
        comps = []
        for new, old in enumerate(order, 1):
            remap[old] = new
            is_big = old in big
            nb = sum(1 for o in order[:new] if o in big)
            comps.append({"id": new, "n": int(cnt[old]), "name": f"T{nb - 1}" if is_big else f"S{new - len(big) - 1}", "keep": bool(is_big)})
        if detect == "transcripts":                              # what a viewer needs: the density, the tissue of every bin, and each piece's box around its bins
            bl = remap[lab]
            boxes = {}
            for c in comps:
                m = occ & (bl == c["id"])
                if m.any():
                    xs, ys = np.nonzero(m)
                    boxes[c["id"]] = [ox + xs.min() * B, oy + ys.min() * B, ox + (xs.max() + 1) * B, oy + (ys.max() + 1) * B]
            self._txview[key] = {"ox": ox, "oy": oy, "bin": B, "counts": Ht, "lab": bl, "boxes": boxes, "thr": thr}
        comp = np.where(ok, remap[cell_lab], -1)
        n0 = int(((comp == 0)).sum())
        if n0:
            comps.insert(0, {"id": 0, "n": n0, "name": "isolated_cells", "keep": False})
        self._auto[key] = (comp, comps)
        return self._auto[key]

    def labels(self, spec):
        """Apply the user's edits. Returns (names, idx per cell (-1 = dropped, -2 = below min_counts), pieces, per-name counts)."""
        gap, mc, mn = float(spec.get("gap", 200)), int(spec.get("min_cells", 500)), int(spec.get("min_counts", 1))
        comp, comps = self.auto(gap, mc, mn, spec.get("detect", "cells"), int(spec.get("tx_min", 0) or 0))
        edits = spec.get("comps", {}) or {}
        auto_off = bool(spec.get("auto_off"))                       # tissue boxes loaded from a GeoJSON define the tissues: detected pieces start switched off
        if auto_off:
            comps = [{**c, "keep": False} for c in comps]
        nm, kp = {}, {}
        for c in comps:
            e = edits.get(str(c["id"]), {})
            nm[c["id"]] = (str(e.get("name") or c["name"]).strip() or c["name"])
            kp[c["id"]] = bool(e.get("keep", c["keep"]))
        ok = comp >= 0
        names_of = {i: nm[i] for i in nm}
        lab = np.empty(self.n, dtype=object)
        keep = np.zeros(self.n, bool)
        for i in nm:
            sel = comp == i
            lab[sel] = nm[i]
            keep[sel] = kp[i]
        for r in spec.get("rects", []) or []:
            sel = self._inside(r) & ok
            lab[sel] = (str(r.get("name") or "region").strip() or "region")
            keep[sel] = bool(r.get("keep", True))
        order = []
        for c in comps:
            if kp[c["id"]] and nm[c["id"]] not in order:
                order.append(nm[c["id"]])
        for r in spec.get("rects", []) or []:
            n = (str(r.get("name") or "region").strip() or "region")
            if r.get("keep", True) and n not in order:
                order.append(n)
        pos = {n: i for i, n in enumerate(order)}
        idx = np.full(self.n, -1, np.int32)
        for n, i in pos.items():
            idx[keep & ok & (lab == n)] = i
        idx[~ok] = -2
        counts = {n: int((idx == i).sum()) for n, i in pos.items()}
        return order, idx, comps, counts

    def _inside(self, r):
        """Cells inside a drawn shape: a rectangle {x0,y0,x1,y1} or a free-hand polygon {pts: [[x, y], ...], holes?: [[[x, y], ...], ...]}."""
        if r.get("pts"):
            import shapely
            g = shape_geom(r)
            x0, y0, x1, y1 = g.bounds
            sel = np.zeros(self.n, bool)
            m = (self.x >= x0) & (self.x <= x1) & (self.y >= y0) & (self.y <= y1)
            if m.any():
                sel[m] = shapely.contains_xy(g, self.x[m], self.y[m])
            return sel
        x0, x1 = sorted((float(r["x0"]), float(r["x1"]))); y0, y1 = sorted((float(r["y0"]), float(r["y1"])))
        return (self.x >= x0) & (self.x <= x1) & (self.y >= y0) & (self.y <= y1)

    # ---- cell labels -------------------------------------------------------------------------------------------------------------
    def inspect_annotation(self, path):
        p = Path(path).expanduser()
        if not p.is_file():
            raise importer.ImportError_(f"not found: {p}")
        d = _read_annotation(p)
        ids = set(self.ids)
        rate = {c: float(d[c].astype(str).isin(ids).mean()) for c in d.columns}
        best = max(rate, key=rate.get)
        if rate[best] == 0:
            raise importer.ImportError_("no column of the annotation file contains any of this run's cell ids")
        from .engine.config import CANDIDATES
        ctype = next((c for c in CANDIDATES["cell_type"] if c in d.columns), None)
        cstate = next((c for c in CANDIDATES["cell_state"] if c in d.columns), None)
        return {"columns": list(d.columns), "n_rows": int(len(d)), "id_col": best, "match": {c: round(r, 4) for c, r in rate.items() if r > 0},
                "type_col": ctype, "state_col": cstate,
                "n_matched": int(d[best].astype(str).isin(ids).sum()),
                "cells_covered": {c: int(len(ids & set(d[c].astype(str)))) for c in d.columns if rate[c] > 0}, "n_cells": self.n}

    def cell_labels(self, ann):
        """(cell_type, cell_state, extra columns) for every cell."""
        clus = None
        if self.clusters is not None:
            clus = self.clusters.map(lambda v: "unclustered" if pd.isna(v) else f"Cluster {int(v)}").values
        extra = pd.DataFrame(index=range(self.n))
        ctype, cstate = clus, clus
        if ann and ann.get("path"):
            d = _read_annotation(Path(ann["path"]).expanduser())
            d[ann["id_col"]] = d[ann["id_col"]].astype(str)
            d = d[~d[ann["id_col"]].duplicated()].set_index(ann["id_col"]).reindex(self.ids)
            for c in d.columns:
                if c not in self.cells.columns:
                    extra[c] = d[c].values
            if ann.get("type_col"):
                ctype = d[ann["type_col"]].astype("string").fillna("unlabelled").values
                cstate = clus if ann.get("state_col") in (None, "") else d[ann["state_col"]].astype("string").fillna("unlabelled").values
            elif ann.get("state_col"):
                cstate = d[ann["state_col"]].astype("string").fillna("unlabelled").values
        if ctype is None:
            ctype = cstate = np.array(["all cells"] * self.n, dtype=object)
        if cstate is None:
            cstate = ctype
        return np.asarray(ctype, dtype=object), np.asarray(cstate, dtype=object), extra


class Light(Session):
    """A run without its counts matrix: cell centroids and per-cell transcript counts only. Enough to detect tissues, draw shapes on them and preview
    the result (the slide gallery of the batch import). The per-cell counts stand in for the panel-gene counts the full import uses, so with a
    'minimum transcripts per cell' above 1 the pieces can differ slightly from the converted ones; the shapes are plain coordinates and do not."""

    def __init__(self, src, log=lambda m: None):
        self.src, self.log, self.meta = src, log, {}
        key = "cells_parquet" if "cells_parquet" in src else "cells_csv"
        p = Path(src[key])
        want = ("x_centroid", "y_centroid", "transcript_counts")
        if p.name.endswith(".parquet"):
            import pyarrow.parquet as pq
            have = set(pq.ParquetFile(p).schema.names)
            c = pd.read_parquet(p, columns=[k for k in want if k in have])
        else:
            c = pd.read_csv(p, usecols=lambda k: k in want)
        for col in want[:2]:
            if col not in c.columns:
                raise importer.ImportError_(f"cells table has no column {col}")
        c = c.dropna(subset=list(want[:2]))
        self.x = c["x_centroid"].values.astype(np.float64)
        self.y = c["y_centroid"].values.astype(np.float64)
        self.n = len(self.x)
        if not self.n:
            raise importer.ImportError_("the cells table is empty")
        self.gene_counts = c["transcript_counts"].fillna(0).values.astype(np.int64) if "transcript_counts" in c.columns else np.ones(self.n, np.int64)
        self._auto, self._txview = {}, {}
        self.xy32 = np.column_stack([self.x, self.y]).astype(np.float32)

    def thumb_polys(self, gap, min_cells, min_counts, detect, tx_min, comp, comps):
        """The carved box of every kept piece (default names), as {piece id: polygons}, for the gallery thumbnails."""
        spec = {"gap": gap, "min_cells": min_cells, "min_counts": min_counts, "detect": detect, "tx_min": tx_min}
        names, idx, _, _ = self.labels(spec)
        by_name = {c["name"]: c["id"] for c in comps}
        return {str(by_name[p["name"]]): p["polys"] for p in polys_payload(tissue_boxes(self, spec, names, idx)) if p["name"] in by_name}

    def thumb(self, gap, min_cells, min_counts, max_pts=12000, detect="cells", tx_min=0):
        """What a thumbnail needs: a thinned set of cells (whole micrometres) with the piece each belongs to, and every piece's name and bounding box."""
        comp, comps = self.auto(gap, min_cells, min_counts, detect, tx_min)
        sel = np.arange(0, self.n, max(1, self.n // max_pts))
        pieces = []
        tv = self._txview.get((round(float(gap), 3), int(min_cells), int(min_counts), "transcripts", int(tx_min))) if detect == "transcripts" else None
        for c in comps:
            m = comp == c["id"]
            if tv is not None and c["id"] in tv["boxes"]:
                pieces.append({**c, "bbox": [round(v) for v in tv["boxes"][c["id"]]]})              # the box around the piece's transcript bins
            elif c["id"] > 0 and m.any():
                pieces.append({**c, "bbox": [round(float(self.x[m].min())), round(float(self.y[m].min())), round(float(self.x[m].max())), round(float(self.y[m].max()))]})
            else:
                pieces.append({**c, "bbox": None})
        return {"n": int(self.n), "bounds": view_bounds(self),
                "xy": np.column_stack([np.round(self.x[sel]), np.round(self.y[sel])]).astype(np.int32).ravel().tolist(), "comp": comp[sel].astype(int).tolist(), "comps": pieces,
                "grid": self.tx_payload(gap, min_cells, min_counts, detect, tx_min), "polys": self.thumb_polys(gap, min_cells, min_counts, detect, tx_min, comp, comps)}


def _read_annotation(p):
    """A cell table: .parquet, or text (.csv / .tsv / .txt, also .gz) with commas or tabs (the separator is taken from the first line)."""
    p = Path(p)
    if p.name.endswith(".parquet"):
        return pd.read_parquet(p)
    import gzip
    with (gzip.open(p, "rt") if p.name.endswith(".gz") else open(p, encoding="utf-8", errors="replace")) as f:
        first = f.readline()
    return pd.read_csv(p, sep="\t" if first.count("\t") > first.count(",") else ",", low_memory=False)


NEG_CATS = ["genomic_control_probe", "deprecated_codeword", "negative_control_probe", "negative_control_codeword", "unassigned_codeword"]   # codeword_category values that are negative controls (predesigned_gene / custom_gene are the positive probes)
NEG_PREFIX = {"GenomicControlProbe": 0, "DeprecatedCodeword": 1, "NegControlProbe": 2, "NegControlCodeword": 3, "UnassignedCodeword": 4}   # fallback by feature name when the file has no codeword_category


def _write_tissue_table(path, tcell, near_rows, names, prefix, out, log):
    """transcript_tissue.parquet from the tissue index of every row of the transcripts file (see write_transcripts)."""
    import pyarrow as pa, pyarrow.parquet as pq
    rows = np.flatnonzero(tcell >= 0)
    cols = {"row": pa.array(rows.astype(np.int64))}
    if "transcript_id" in pq.ParquetFile(path).schema_arrow.names:
        cols["transcript_id"] = pq.read_table(path, columns=["transcript_id"])["transcript_id"].combine_chunks().take(pa.array(rows))
    ti = tcell[rows].astype(np.int32)
    nm = [str(n) for n in names]
    cols["tissue"] = pa.DictionaryArray.from_arrays(pa.array(ti), pa.array(nm))
    cols["sample"] = pa.DictionaryArray.from_arrays(pa.array(ti), pa.array([f"{prefix}_{n}" if prefix else n for n in nm]))
    near = np.zeros(len(tcell), bool)
    near[near_rows] = True
    cols["assignment"] = pa.DictionaryArray.from_arrays(pa.array(near[rows].astype(np.int8)), pa.array(["cell", "nearest_cell"]))
    pq.write_table(pa.table(cols), Path(out) / "transcript_tissue.parquet", compression="zstd")
    log(f"transcript_tissue.parquet: {len(rows):,} of {len(tcell):,} transcripts assigned to a tissue ({int(near[rows].sum()):,} by the nearest cell)")
    return int(len(rows))


def _txt(v):
    """A value of a transcripts.parquet text column as str: older Xenium output stores these columns as binary (bytes)."""
    return v.decode("utf-8", "replace") if isinstance(v, (bytes, bytearray)) else str(v)


def write_transcripts(sess, names, idx, ids, path, out, log, min_qv=20.0, max_assign_um=25.0, store=True, table_prefix=None, threads=1):
    """Decoded transcripts (quality >= min_qv, panel genes only) of the kept cells, stored sorted by (tissue, gene) so that one gene of one tissue is a
    contiguous slice: transcripts/{tx_x,tx_y}.npy (float32, um), tx_ptr.npy (offsets, length n_tissues * n_genes + 1), index.json.
    A transcript belongs to the tissue of the cell it was assigned to; unassigned panel-gene ones go to the tissue of the nearest kept cell, at any
    distance (all quality >= min_qv panel-gene transcripts are kept; those farther than max_assign_um from every kept cell are counted separately as
    "unassigned_far" in the QC block). Unassigned negative-control transcripts still need a kept cell within max_assign_um, else they are dropped.
    Transcripts of cells that were left out are dropped.
    index.json also carries a "qc" block for the QC report: per tissue the quality-value histogram (QV 0..40, panel-gene transcripts assigned to a
    cell of the tissue, any quality) and the counts of quality >= min_qv transcripts assigned to a cell / unassigned but within max_assign_um of the
    tissue; and for the whole run the same quality-filtered counts over all panel-gene transcripts (assigned to any cell or not).
    table_prefix (a string, "" for none) also writes transcript_tissue.parquet: for every row of transcripts.parquet that belongs to a tissue, its row number in
    the file (0-based), transcript_id (when the file has it), tissue, sample (= table_prefix + "_" + tissue, the name in a merged batch project) and
    assignment ("cell" = its own cell is kept; "nearest_cell" = unassigned, within max_assign_um of a kept cell; that rule only for quality >= min_qv
    panel-gene / negative-control rows). Rows of no tissue are not listed. store=False skips the .npy files and index.json (table only).
    threads = cores the nearest-cell search may use.
    Returns (number of points stored or None, number of table rows or None)."""
    import pyarrow.parquet as pq, pyarrow.compute as pc
    from scipy.spatial import cKDTree
    t0 = time.time()
    has_qv = "qv" in pq.ParquetFile(path).schema_arrow.names
    names_ = pq.ParquetFile(path).schema_arrow.names
    has_cat = "codeword_category" in names_
    cols = ["x_location", "y_location", "feature_name", "cell_id"] + (["qv"] if has_qv else []) + (["codeword_category"] if has_cat else [])
    tb = pq.read_table(path, columns=cols)
    n_all = tb.num_rows
    G, nT = len(sess.genes), len(names)
    gidx = {g: i for i, g in enumerate(sess.genes)}
    fn = pc.dictionary_encode(tb["feature_name"]).combine_chunks()
    lut = np.array([gidx.get(_txt(n), -1) for n in fn.dictionary.to_pylist()], np.int64)
    gene = lut[fn.indices.to_numpy(zero_copy_only=False)]
    if has_cat:   # negative-control class of every row (index in NEG_CATS, -1 = not a negative control), from the codeword_category column
        cc = pc.dictionary_encode(tb["codeword_category"]).combine_chunks()
        cl = np.array([NEG_CATS.index(_txt(n)) if _txt(n) in NEG_CATS else -1 for n in cc.dictionary.to_pylist()], np.int8)
        negcat = cl[cc.indices.to_numpy(zero_copy_only=False)]
        del cc
    else:
        cl = np.array([next((k for p, k in NEG_PREFIX.items() if _txt(n).startswith(p)), -1) for n in fn.dictionary.to_pylist()], np.int8)
        negcat = cl[fn.indices.to_numpy(zero_copy_only=False)]
    cd = pc.dictionary_encode(tb["cell_id"]).combine_chunks()
    pos_all = pd.Index(sess.ids).get_indexer([_txt(v) for v in cd.dictionary.to_pylist()])   # row in the scanned cells, -1 = unassigned (ids as text: the cell table's are text; older runs have integer ids in both files)
    cpos = pos_all[cd.indices.to_numpy(zero_copy_only=False)]
    if n_all and (cpos >= 0).mean() < 0.01:
        ex = [v for v in cd.dictionary.to_pylist() if v is not None][:3]
        raise importer.ImportError_(f"fewer than 1% of the {n_all:,} transcripts carry a cell id that is in the cell table (transcripts: {ex}, cell table: {list(sess.ids[:3])}); "
                                    "the two files do not use the same cell ids")
    if table_prefix is not None:                                                  # tissue of every row through its own cell (any quality, any gene), by row number
        tcell = np.full(n_all, -1, np.int16)
        a0 = cpos >= 0
        tcell[a0] = np.maximum(idx[cpos[a0]], -1)
        del a0
        rid = np.arange(n_all, dtype=np.int64 if n_all >= 2 ** 31 else np.int32)   # row number of every row still in x / y / gene / cpos below
    qv = tb["qv"].to_numpy(zero_copy_only=False).astype(np.float32) if has_qv else None
    x = tb["x_location"].to_numpy(zero_copy_only=False).astype(np.float32)
    y = tb["y_location"].to_numpy(zero_copy_only=False).astype(np.float32)
    del tb
    tr_ = idx[cpos[cpos >= 0]] if len(idx) else np.zeros(0, np.int64)                 # tissue of every row through its own cell (any quality, any feature)
    tissue_rows = np.bincount(tr_[tr_ >= 0], minlength=nT)
    del tr_
    qvh = np.zeros((nT, 41), np.int64)
    gq_mean = gq_n = None                                                        # per tissue and gene: mean QV and number of transcripts (every quality, assigned to a cell of the tissue)
    if has_qv:
        ta = np.full(n_all, -1, np.int64)
        asg = cpos >= 0
        ta[asg] = idx[cpos[asg]]
        sel = (gene >= 0) & (ta >= 0) & np.isfinite(qv)
        qvh = np.bincount(ta[sel] * 41 + np.clip(np.floor(qv[sel]), 0, 40).astype(np.int64), minlength=nT * 41).reshape(nT, 41)
        gk = ta[sel] * G + gene[sel]
        gq_n = np.bincount(gk, minlength=nT * G).astype(np.int64)
        gq_s = np.bincount(gk, weights=qv[sel].astype(np.float64), minlength=nT * G)
        gq_mean = np.where(gq_n > 0, gq_s / np.maximum(gq_n, 1), np.nan).reshape(nT, G)
        gq_n = gq_n.reshape(nT, G)
        del ta, asg, sel, gk, gq_s
        ok = qv >= min_qv
        x, y, gene, cpos, negcat = x[ok], y[ok], gene[ok], cpos[ok], negcat[ok]
        if table_prefix is not None:
            rid = rid[ok]
    del qv
    log(f"transcripts: {n_all:,} rows, {len(x):,} with quality >= {min_qv:g}")
    tis = np.full(len(x), -1, np.int64)
    assigned = cpos >= 0
    tis[assigned] = idx[cpos[assigned]]                                          # -1 / -2 = that cell was left out
    run_panel, run_assigned = int((gene >= 0).sum()), int(((gene >= 0) & assigned).sum())
    neg = negcat >= 0
    un = (~assigned) & ((gene >= 0) | neg)
    farm = np.zeros(len(x), bool)                                                # panel-gene rows kept although no kept cell is within max_assign_um
    if un.any() and len(ids):
        kept = np.flatnonzero(idx >= 0)
        tree = cKDTree(np.column_stack([sess.x[kept], sess.y[kept]]))
        pts = np.column_stack([x[un], y[un]])
        d, j = tree.query(pts, distance_upper_bound=max_assign_um, workers=threads)
        t = np.full(len(d), -1, np.int64); okk = np.isfinite(d); t[okk] = idx[kept[j[okk]]]
        fg = (~okk) & (gene[un] >= 0)
        if fg.any():                                                             # no distance limit for panel genes: the tissue of the nearest kept cell
            _, j2 = tree.query(pts[fg], workers=threads)
            t[fg] = idx[kept[j2]]
            farm[np.flatnonzero(un)[fg]] = True
        del pts
        tis[un] = t
        if table_prefix is not None:
            tcell[rid[un][okk]] = t[okk]
            near = rid[un][okk]
    n_table = None
    if table_prefix is not None:
        n_table = _write_tissue_table(path, tcell, near if (un.any() and len(ids)) else np.zeros(0, np.int64), names, table_prefix, out, log)
        del tcell, rid
    if not store:
        return None, n_table
    cnt_a = np.bincount(tis[(gene >= 0) & assigned & (tis >= 0)], minlength=nT)
    cnt_far = np.bincount(tis[(gene >= 0) & farm & (tis >= 0)], minlength=nT)
    cnt_u = np.bincount(tis[(gene >= 0) & ~assigned & (tis >= 0)], minlength=nT) - cnt_far     # unassigned within max_assign_um of a kept cell
    kn = neg & (tis >= 0)                                                        # negative-control probes (same quality filter and tissue rule), sorted by tissue
    nx, ny, nt_, ncat = x[kn], y[kn], tis[kn], negcat[kn]
    on = np.argsort(nt_, kind="stable")
    nptr = np.r_[0, np.cumsum(np.bincount(nt_, minlength=nT))].astype(np.int64)
    keep = (gene >= 0) & (tis >= 0)
    x, y, gene, tis = x[keep], y[keep], gene[keep], tis[keep]
    key = tis * G + gene
    o = np.argsort(key, kind="stable")
    ptr = np.r_[0, np.cumsum(np.bincount(key, minlength=nT * G))].astype(np.int64)
    d = Path(out) / "transcripts"
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "tx_x.npy", x[o]); np.save(d / "tx_y.npy", y[o]); np.save(d / "tx_ptr.npy", ptr)
    np.save(d / "ctl_x.npy", nx[on]); np.save(d / "ctl_y.npy", ny[on]); np.save(d / "ctl_ptr.npy", nptr); np.save(d / "ctl_cat.npy", ncat[on].astype(np.int8))
    qc = {"assigned": cnt_a.tolist(), "unassigned_near": cnt_u.tolist(), "unassigned_far": cnt_far.tolist(), "tissue_rows": tissue_rows.tolist(), "qv_hist": qvh.tolist() if has_qv else None,
          "gene_qv_mean": [[None if v != v else round(float(v), 2) for v in row] for row in gq_mean] if gq_mean is not None else None,
          "gene_qv_n": gq_n.tolist() if gq_n is not None else None,
          "run": {"rows": int(n_all), "q_filtered_panel": run_panel, "q_filtered_panel_assigned": run_assigned}}
    (d / "index.json").write_text(json.dumps({"samples": list(names), "n_genes": G, "n": int(len(x)), "min_qv": min_qv, "n_rows_in_file": int(n_all),
                                              "unassigned_max_distance_um": max_assign_um, "genes_beyond_max_distance": True, "qc": qc,
                                              "ctrl": {"categories": NEG_CATS, "source": "codeword_category" if has_cat else "feature_name", "n": int(len(nx)), "min_qv": min_qv if has_qv else None}}))
    log(f"negative controls ({'codeword_category' if has_cat else 'feature names'}): {len(nx):,} stored, " + ", ".join(f"{c} {int((ncat == i).sum()):,}" for i, c in enumerate(NEG_CATS)))
    log(f"transcripts: {len(x):,} stored ({(x.nbytes * 2) / 1e6:.0f} MB) in {time.time() - t0:.0f} s")
    return int(len(x)), n_table


def write_signatures(out, genes, log=print):
    """signatures.json for the Signatures window / violins: the curated definitions (engine/signature_defs.py, kidney / interferon / lupus
    nephritis) resolved against this project's gene panel. No coherence statistics (those need cell types that are comparable between
    samples), so the window shows no "weak" flag. Returns the number of signatures with at least one gene on the panel."""
    from .engine.signatures import resolve
    sigs = resolve(list(genes))
    (Path(out) / "signatures.json").write_text(json.dumps({"signatures": sigs, "isg_candidates": [], "n_cells_tested": 0}))
    n = sum(1 for x in sigs if x["genes"])
    log(f"signatures: {n} of {len(sigs)} curated signatures have genes on this panel")
    return n


def _fmt_int(x):
    return int(round(float(x)))


def shape_geom(r):
    """The polygon of a drawn shape: a rectangle {x0,y0,x1,y1}, or a free-hand outline {pts} with optional holes {holes: [ring, ...]} (a tissue box
    loaded from a GeoJSON can have holes, e.g. a tissue box with another tissue's box cut out of it). Made valid."""
    from shapely.geometry import box, Polygon
    from shapely.validation import make_valid
    if not r.get("pts"):
        x0, x1 = sorted((float(r["x0"]), float(r["x1"]))); y0, y1 = sorted((float(r["y0"]), float(r["y1"])))
        return box(x0, y0, x1, y1)
    rings = [r["pts"]] + list(r.get("holes") or [])
    if len(rings) > 200:
        raise importer.ImportError_("a drawn shape has too many holes")
    arr = []
    for q in rings:
        a = np.asarray(q, float)
        if a.ndim != 2 or a.shape[1] != 2 or len(a) < 3 or len(a) > 5000 or not np.isfinite(a).all():
            raise importer.ImportError_("a drawn shape needs 3 to 5000 valid points per outline")
        arr.append(a)
    return make_valid(Polygon(arr[0], arr[1:]))


def _ring_pts(coords):
    """Shapely ring coordinates -> [[x, y], ...] without the closing point, rounded to 0.1 um."""
    c = [[round(float(x), 1), round(float(y), 1)] for x, y in list(coords)[:-1]]
    return c


def _poly_shape(q, name, keep):
    """A shapely Polygon -> a drawn shape: a rectangle when it is an axis-aligned box, else a free-hand outline (with holes), simplified to at most 2000
    points per outline (tolerance 0.5 um, doubled until it fits)."""
    x0, y0, x1, y1 = q.bounds
    if not q.interiors and len(q.exterior.coords) <= 6 and abs(q.area - (x1 - x0) * (y1 - y0)) <= 1e-6 * max(1.0, q.area):
        return {"x0": round(x0, 1), "y0": round(y0, 1), "x1": round(x1, 1), "y1": round(y1, 1), "name": name, "keep": keep}
    npts = lambda g: max([len(g.exterior.coords)] + [len(i.coords) for i in g.interiors])
    g, tol = q, 0.5
    while npts(g) > 2001:
        g = q.simplify(tol, preserve_topology=True)
        tol *= 2
    r = {"pts": _ring_pts(g.exterior.coords), "name": name, "keep": keep}
    holes = [_ring_pts(i.coords) for i in g.interiors if len(i.coords) >= 4]
    if holes:
        r["holes"] = holes
    return r


def load_tissue_boxes(path, name_col=None):
    """Read a tissue-box GeoJSON (the <name>_regions.geojson written by write_regions, or any FeatureCollection of polygons in Xenium um) back as the
    shapes that separate the tissues on the slide. Only for telling tissues apart: nothing here becomes the project's anatomical regions
    (regions.geojson, the outlines within a tissue). Returns
    {shapes, tissues: [{name, polys}], name_col, n_tissues, n_removed, n_skipped, source}. Each exported tissue box becomes a kept shape named after the
    tissue (a tissue of several pieces becomes several shapes of one name); an exported drawn region that was switched off becomes a removal shape;
    kept drawn regions are skipped (their tissue box already holds them, cut where it met another tissue). A file without tissue features (another
    tool's) has every polygon loaded as a kept shape, or a removal shape when its "keep" property is false. The import then switches the detected
    tissues off (spec auto_off), so these shapes define the tissues."""
    from shapely.geometry import shape as sshape
    from shapely.validation import make_valid
    p = Path(path).expanduser()
    if not p.is_file():
        raise importer.ImportError_(f"not found: {p}")
    try:
        gj = json.loads(p.read_text())
    except (OSError, ValueError) as e:
        raise importer.ImportError_(f"{p.name} is not a readable GeoJSON file: {e}")
    feats = gj.get("features") if isinstance(gj, dict) and gj.get("type") == "FeatureCollection" else [gj] if isinstance(gj, dict) and gj.get("type") == "Feature" else None
    if not feats:
        raise importer.ImportError_(f"{p.name} is not a GeoJSON FeatureCollection with features")
    col = (name_col or "").strip() or gj.get("name_column") or "sample_id_unblinded"

    def name_of(props, i):
        for k in (col, gj.get("name_column"), "sample_id_unblinded", "name", "Name", "label", "id"):
            if k and props.get(k) not in (None, ""):
                return str(props[k]).strip()
        return f"region_{i + 1}"

    def polys(f):
        try:
            g = make_valid(sshape(f.get("geometry")))
        except Exception:
            return []
        return [q for q in ([g] if g.geom_type == "Polygon" else getattr(g, "geoms", [])) if q.geom_type == "Polygon" and not q.is_empty and q.area > 0]

    has_tissue = any((f.get("properties") or {}).get("role") == "tissue" for f in feats)
    keep_shapes, removal, skipped, tissues = [], [], 0, {}
    for i, f in enumerate(feats):
        pr = f.get("properties") or {}
        role, keep = pr.get("role"), pr.get("keep", True) is not False
        qs = polys(f)
        if not qs:
            skipped += 1
            continue
        nm = name_of(pr, i)
        if has_tissue and role == "drawn region" and keep:
            skipped += 1                                                       # held by its tissue box
            continue
        if not keep:
            removal += [_poly_shape(q, nm, False) for q in qs]
            continue
        keep_shapes += [_poly_shape(q, nm, True) for q in qs]
        tissues.setdefault(nm, []).extend(qs)
    if not keep_shapes:
        raise importer.ImportError_(f"{p.name} has no tissue boxes or drawn shapes to load")
    shapes = keep_shapes + removal                                             # removal shapes last, so they take their cells out of any tissue
    if len(shapes) > 200:
        raise importer.ImportError_(f"{p.name} has {len(shapes)} polygons; at most 200 shapes per run can be loaded")
    tl = [{"name": n, "polys": [[r["pts"]] + r.get("holes", []) if r.get("pts") else [[[r["x0"], r["y0"]], [r["x1"], r["y0"]], [r["x1"], r["y1"]], [r["x0"], r["y1"]]]]
                                for r in keep_shapes if r["name"] == n]} for n in tissues]
    return {"shapes": shapes, "tissues": tl, "name_col": col, "n_tissues": len(tissues), "n_removed": len(removal), "n_skipped": skipped, "source": str(p)}


def find_tissue_file(folder, run):
    """The tissue-box file <run>_regions.geojson of a batch run: in `folder`, in folder/<run>/, or one more level down (an output folder of an earlier batch import). None if absent."""
    d, fn = Path(folder).expanduser(), regions_filename(run)
    for c in (d / fn, d / run / fn):
        if c.is_file():
            return c
    hits = sorted(d.glob(f"*/{fn}")) + sorted(d.glob(f"*/*/{fn}"))
    return hits[0] if hits else None


def regions_filename(slide):
    """<slide name>_regions.geojson, the slide name made safe for a file name."""
    return (re.sub(r"[^A-Za-z0-9._-]+", "_", str(slide)).strip("_") or "slide") + "_regions.geojson"


ERODE_UM = 20.0                                 # boxes of neighbouring tissues keep this distance from the cut between them (each is pulled back by it)


def _split_overlap(o, xa, ya, xb, yb):
    """The parts of the overlap polygon `o` that lie nearer to the cells of A and nearer to the cells of B: (part of A, part of B). Done on a pixel grid
    (20 um or coarser) by the nearest cell, so a cell of either tissue is always on its own side."""
    import shapely
    from scipy.spatial import cKDTree
    from shapely.geometry import box as sbox
    from shapely.ops import unary_union
    x0, y0, x1, y1 = o.bounds
    px = max(20.0, max(x1 - x0, y1 - y0) / 300.0)
    nx, ny = int(math.ceil((x1 - x0) / px)), int(math.ceil((y1 - y0) / px))
    gx, gy = np.meshgrid(x0 + (np.arange(nx) + 0.5) * px, y0 + (np.arange(ny) + 0.5) * px, indexing="ij")
    inside = shapely.contains_xy(o, gx.ravel(), gy.ravel()).reshape(nx, ny)
    ka, kb = max(1, len(xa) // 60000), max(1, len(xb) // 60000)
    da = cKDTree(np.column_stack([xa[::ka], ya[::ka]])).query(np.column_stack([gx.ravel(), gy.ravel()]))[0].reshape(nx, ny)
    db = cKDTree(np.column_stack([xb[::kb], yb[::kb]])).query(np.column_stack([gx.ravel(), gy.ravel()]))[0].reshape(nx, ny)
    near_a = inside & (da <= db)
    rects = []
    for iy in range(ny):                                     # runs of A pixels in every row, as rectangles
        row = near_a[:, iy]
        ix = 0
        while ix < nx:
            if row[ix]:
                j = ix
                while j + 1 < nx and row[j + 1]:
                    j += 1
                rects.append(sbox(x0 + ix * px, y0 + iy * px, x0 + (j + 1) * px, y0 + (iy + 1) * px))
                ix = j + 1
            else:
                ix += 1
    oa = unary_union(rects).intersection(o) if rects else o.difference(o)
    return oa, o.difference(oa)


def _disjoint_boxes(items, erode_um=ERODE_UM):
    """Make tissue bounding boxes non-overlapping and apart. `items` = [(name, box polygon, x, y)] with the coordinates of the tissue's own cells;
    returns the geometries in the same order. Where two boxes overlap, the overlap is divided between the two tissues by which one's cells are nearer
    (each box loses the part that went to the other), and then each box is pulled back by `erode_um` from the other, so the two never touch: a
    transcript can fall in at most one box. A box may become a rectangle with a notch, or several pieces."""
    g = [it[1] for it in items]
    pairs = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            o = g[i].intersection(g[j])
            if o.is_empty or o.area <= 0:
                continue
            oa, ob = _split_overlap(o, items[i][2], items[i][3], items[j][2], items[j][3])
            g[i], g[j] = g[i].difference(ob), g[j].difference(oa)
            pairs.append((i, j))
    pre = list(g)
    for i, j in pairs:
        g[i] = g[i].difference(pre[j].buffer(erode_um))
        g[j] = g[j].difference(pre[i].buffer(erode_um))
    return g



def _data_bounds(sess):
    """The extent [x0, y0, x1, y1] of the run's data: all its cells, and its transcript bins when those were read (transcript detection)."""
    b = [float(sess.x.min()), float(sess.y.min()), float(sess.x.max()), float(sess.y.max())]
    g = getattr(sess, "_txg", None)
    if g is not None:
        ox, oy, H = g
        b = [min(b[0], float(ox)), min(b[1], float(oy)), max(b[2], float(ox) + H.shape[0] * TX_BIN_UM), max(b[3], float(oy) + H.shape[1] * TX_BIN_UM)]
    return b


# The full Xenium slide, in the frame of the cell coordinates: the largest imaged area seen on 11 real runs (10 kidney runs, Xenium software 3.0 and
# 4.0, plus one more), 10,867 um wide (51,138 px at 0.2125 um/px, run slidearray0073) by 23,817 um long (112,078 px, run slidearray0140). Each run's own
# morphology image covers only its scanned region (from 7,862 x 20,208 um to that size), so the image is not read: the same slide is used for every
# run. Coordinates start at the corner of the area that was exported; some outputs hold only half a slide (the lung preview run: data within
# 10,816 x 10,128 um), which the coordinates show: all data within half the length (or width).
XENIUM_SLIDE_UM = (10867.0, 23817.0)


def slide_bounds(sess):
    """The slide's extent [x0, y0, x1, y1] in micrometres, which tissue boxes and drawn shapes are capped at: the Xenium imaged area (XENIUM_SLIDE_UM),
    each dimension halved when all the data lies within its first half (an output that holds half a slide), and widened if any data lies outside it."""
    d = _data_bounds(sess)
    ext = []
    for k, full in enumerate(XENIUM_SLIDE_UM):
        ext.append(full / 2 if d[k + 2] <= full / 2 else full)
    return [min(d[0], 0.0), min(d[1], 0.0), max(d[2], ext[0]), max(d[3], ext[1])]


def view_bounds(sess):
    """What the review canvases show: the data's extent plus room for the largest box growth (200 um), within the slide."""
    d, sl, m = _data_bounds(sess), slide_bounds(sess), max(DILATE_STEPS_UM)
    return [max(sl[0], d[0] - m), max(sl[1], d[1] - m), min(sl[2], d[2] + m), min(sl[3], d[3] + m)]


DILATE_STEPS_UM = (200.0, 100.0, 50.0, 20.0)    # a detected tissue box grows by the first of these that keeps it clear of every other box


def _dilate_boxes(boxes, steps=DILATE_STEPS_UM, keep_um=ERODE_UM, limit=None, fixed=None):
    """Grow every tissue box on both sides of a dimension, first along its length (longer side), then along its breadth. For each dimension the
    first of `steps` (200, 100, 50, 20 um) is used whose grown box comes no nearer than `keep_um` to any other box (detected, drawn rectangle or
    free-hand) than it already was, i.e. makes no new overlap or near-contact; if none does, that dimension keeps its original extent. Growth stops at
    `limit` = [x0, y0, x1, y1] (the slide). Boxes with fixed[i] true (drawn shapes) are not grown but are obstacles. Boxes are processed in order
    against the current (possibly already grown) others. Returns new geometries."""
    from shapely.geometry import box
    out = list(boxes)
    for i, b in enumerate(boxes):
        if fixed and fixed[i]:
            continue
        x0, y0, x1, y1 = b.bounds
        dims = ("x", "y") if (x1 - x0) >= (y1 - y0) else ("y", "x")
        for d in dims:
            cur = box(x0, y0, x1, y1).buffer(keep_um)
            others = [(o, cur.intersection(o).area) for j, o in enumerate(out) if j != i]
            for st in steps:
                c = [x0 - st, y0, x1 + st, y1] if d == "x" else [x0, y0 - st, x1, y1 + st]
                if limit is not None:                       # capped at the slide, but never shrunk below the original box
                    c = [max(c[0], min(limit[0], x0)), max(c[1], min(limit[1], y0)), min(c[2], max(limit[2], x1)), min(c[3], max(limit[3], y1))]
                grown = box(*c).buffer(keep_um)
                if all(grown.intersection(o).area <= a + 1e-6 for o, a in others):
                    x0, y0, x1, y1 = c
                    break
        out[i] = box(x0, y0, x1, y1)
    return out


def _drawn_geoms(spec, limit):
    """{tissue name: union of its kept drawn shapes, capped at the slide `limit`}."""
    from shapely.geometry import box
    from shapely.ops import unary_union
    by = {}
    for r in spec.get("rects", []) or []:
        if not r.get("keep", True):
            continue
        n = str(r.get("name") or "region").strip() or "region"
        g = shape_geom(r).intersection(box(*limit))
        if not g.is_empty and g.area > 0:
            by.setdefault(n, []).append(g)
    return {n: unary_union(gs) for n, gs in by.items()}


def tissue_boxes(sess, spec, order, idx):
    """The bounding box of every tissue after the cell labels in `idx` (index into `order`): [{name, box, geom, n, x, y, cut, outside}], a tissue made by
    drawn shapes has the drawn outline as its box (not grown); a detected tissue's box is grown by up to 200 um (capped at the slide) where no other box
    is near (see _dilate_boxes); boxes never overlap (see _disjoint_boxes); with spec cut_off (the default) the switched-off shapes (drawn with "use"
    unticked) are cut out of the boxes, so nothing in such an area belongs to a tissue box. With transcript detection a box also grows to the transcript bins of its piece. Used by the export and by the review windows, so
    what is shown is what is written."""
    import shapely
    from shapely.geometry import box
    tis, fixed, tx_boxes, tx_comp = [], [], None, None
    limit = slide_bounds(sess)
    drawn = _drawn_geoms(spec, limit)
    if spec.get("detect") == "transcripts":
        k = (float(spec.get("gap", 200)), int(spec.get("min_cells", 500)), int(spec.get("min_counts", 1)), "transcripts", int(spec.get("tx_min", 0) or 0))
        tx_comp, _ = sess.auto(*k)
        v = sess._txview.get((round(k[0], 3),) + k[1:])
        tx_boxes = v["boxes"] if v else None
    for i, n in enumerate(order):
        m = idx == i
        if m.any():
            x, y = sess.x[m], sess.y[m]
            x0, y0, x1, y1 = float(x.min()), float(y.min()), float(x.max()), float(y.max())
            if tx_boxes:                                                       # grow the box to the piece's transcript bins (a piece the drawn shapes did not split)
                for c in np.unique(tx_comp[m]):
                    if c > 0 and c in tx_boxes and ((tx_comp == c) & m).sum() >= 0.95 * (tx_comp == c).sum():
                        b = tx_boxes[c]
                        x0, y0, x1, y1 = min(x0, b[0]), min(y0, b[1]), max(x1, b[2]), max(y1, b[3])
            if n in drawn:                                                     # a drawn tissue: its box is the outline that was drawn (plus any of its cells outside it)
                g = drawn[n]
                rest = ~shapely.intersects_xy(g.buffer(1.0), x, y)
                if rest.any():
                    g = g.union(box(float(x[rest].min()), float(y[rest].min()), float(x[rest].max()), float(y[rest].max())))
                tis.append((n, g, x, y)); fixed.append(True)
            else:
                tis.append((n, box(x0, y0, x1, y1), x, y)); fixed.append(False)
    # grow the detected boxes (200 -> 100 -> 50 -> 20 um, capped at the slide) where that touches no other box, then carve any overlap
    tis = [(n, b, x, y) for (n, _, x, y), b in zip(tis, _dilate_boxes([t[1] for t in tis], limit=limit, fixed=fixed))]
    geoms = _disjoint_boxes(tis)
    raw = list(geoms)                                                      # before the switched-off shapes are cut out (the batch page's geometry-only export cuts them itself)
    if spec.get("cut_off", True):                                          # areas drawn with "use" unticked are cut out of every tissue box
        from shapely.ops import unary_union
        off = [shape_geom(r) for r in spec.get("rects", []) or [] if not r.get("keep", True)]
        if off:
            off = unary_union(off)
            geoms = [g.difference(off) for g in geoms]
    out = []
    for (n, b, x, y), g, g0 in zip(tis, geoms, raw):
        cut = not g.equals(b)
        out.append({"name": n, "box": b, "geom": g, "raw": g0, "n": int(len(x)), "x": x, "y": y, "cut": cut,
                    "outside": int((~shapely.intersects_xy(g, x, y)).sum()) if cut else 0})
    return out


def polys_payload(geoms):
    """For the page: [{name, polys: [[ring, ...], ...]}] with rings as [[x, y], ...] (a tissue may be several polygons)."""
    out = []
    for t in geoms:
        g = t["geom"]
        parts = [g] if g.geom_type == "Polygon" else [q for q in getattr(g, "geoms", []) if q.geom_type == "Polygon"]
        out.append({"name": t["name"], "polys": [[[[round(float(a), 1), round(float(b), 1)] for a, b in [*q.exterior.coords]]] + [[[round(float(a), 1), round(float(b), 1)] for a, b in r.coords] for r in q.interiors] for q in parts if not q.is_empty]})
    return out


def write_regions(sess, spec, order, idx, out, log, name_col="sample_id_unblinded", slide="slide", extra_dir=None):
    """Export the import-time selection as GeoJSON (Xenium micrometre frame, y down, same as the cell table): one Feature per drawn
    region (switched-off ones too, keep false, so the file can be loaded again) and one bounding box per kept tissue (with spec cut_off, the default,
    the switched-off shapes are cut out of the boxes). Written to <out>/<slide>_regions.geojson, and also to `extra_dir` when one is given. A rectangle is exported as its box; a free-hand region as the outline that was traced
    (the same valid polygon that selects the cells), not its convex hull. Tissue boxes never overlap or touch each other (see _disjoint_boxes: the overlap is divided between the two tissues by which one's cells are nearer, and each box is pulled back 20 um from the other). Returns the number of features."""
    from shapely.geometry import box, mapping, Polygon
    from shapely.validation import make_valid
    name_col = (name_col or "").strip() or "sample_id_unblinded"
    feats = []
    ok = sess.gene_counts >= int(spec.get("min_counts", 1))
    slide_box = box(*slide_bounds(sess))
    for r in spec.get("rects", []) or []:
        name = str(r.get("name") or "region").strip() or "region"
        g, kind = shape_geom(r), ("freehand" if r.get("pts") else "rectangle")
        g = g.intersection(slide_box)                                      # a drawn shape is capped at the slide
        if g.is_empty:
            continue
        feats.append({"type": "Feature", "geometry": mapping(g),
                      "properties": {name_col: name, "kind": kind, "role": "drawn region", "keep": bool(r.get("keep", True)),
                                     "n_cells": int((sess._inside(r) & ok).sum()), "area_um2": round(float(g.area), 1)}})
    for t in tissue_boxes(sess, spec, order, idx):
        if t["outside"]:
            log(f"tissue {t['name']}: {t['outside']} of its {t['n']:,} cells lie outside its box after the overlap with a neighbouring tissue was removed")
        feats.append({"type": "Feature", "geometry": mapping(t["geom"]),
                      "properties": {name_col: t["name"], "kind": "bounding box", "role": "tissue", "keep": True, "n_cells": t["n"], "switched_off_cut": bool(spec.get("cut_off", True)),
                                     "area_um2": round(float(t["geom"].area), 1), "overlap_removed": bool(t["cut"]), "cells_outside": t["outside"]}})
    fc = {"type": "FeatureCollection", "coordinate_unit": "um (Xenium frame, y down)", "name_column": name_col, "features": feats}
    fn, txt = regions_filename(slide), json.dumps(fc)
    (out / fn).write_text(txt)
    where = [str(out / fn)]
    if extra_dir:
        d = Path(extra_dir).expanduser()
        d.mkdir(parents=True, exist_ok=True)
        (d / fn).write_text(txt)
        where.append(str(d / fn))
    log(f"exported {len(feats)} outline(s) (tissue boxes and drawn shapes) to " + " and ".join(where))
    return len(feats)


def _polygonal(g):
    """Only the polygon parts of a geometry (a cut or make_valid can leave lines or points behind)."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    if g.geom_type in ("Polygon", "MultiPolygon"):
        return g
    parts = [q for q in getattr(g, "geoms", []) if q.geom_type in ("Polygon", "MultiPolygon")]
    return unary_union(parts) if parts else Polygon()


def _polys_geom(polys):
    """A tissue outline as the page holds it, [[outer ring, hole, ...], ...] with rings as [[x, y], ...], -> one valid geometry."""
    from shapely.ops import unary_union
    if not isinstance(polys, list) or not polys or len(polys) > 200:
        raise importer.ImportError_("a tissue outline needs 1 to 200 parts")
    return _polygonal(unary_union([shape_geom({"pts": p[0], "holes": p[1:]}) for p in polys if isinstance(p, list) and p]))


def write_tissue_geometry(folder, slide, shapes, tissues=None, cut_off=True, name_col="sample_id_unblinded"):
    """Write <slide>_regions.geojson into `folder` from geometry alone: no cells or transcripts are read. One Feature per tissue: its outline
    (`tissues` = [{name, polys}] as already built on the batch page, or, when None, the union of the kept shapes of each name, e.g. tissue boxes
    loaded from a GeoJSON), with the switched-off shapes cut out when `cut_off`; and one Feature per switched-off shape (keep false, so the file
    loads back with it switched off). Properties: the tissue name (under `name_col`) and the slide. Tissues are not made disjoint here: overlapping
    ones are written as they are and reported. Returns (path, [[name, name], ...] of overlapping tissues, [names left empty by the cut])."""
    from shapely.geometry import mapping
    from shapely.ops import unary_union
    nm = lambda r: str(r.get("name") or "region").strip() or "region"
    shapes = [r for r in (shapes or []) if isinstance(r, dict)]
    by = {}
    if tissues is None:
        for r in shapes:
            if r.get("keep", True):
                by.setdefault(nm(r), []).append(shape_geom(r))
    else:
        for t in tissues:
            if isinstance(t, dict):
                by.setdefault(nm(t), []).append(_polys_geom(t.get("polys")))
    geoms = {n: _polygonal(unary_union(gs)) for n, gs in by.items()}
    off = [r for r in shapes if not r.get("keep", True)]
    if cut_off and off:
        cut = unary_union([shape_geom(r) for r in off])
        geoms = {n: _polygonal(g.difference(cut)) for n, g in geoms.items()}
    empty = [n for n, g in geoms.items() if g.is_empty or g.area <= 0]
    geoms = {n: g for n, g in geoms.items() if n not in empty}
    if not geoms:
        raise importer.ImportError_("no tissue is switched on" if not empty else "nothing is left of the tissues after the switched-off shapes were cut out")
    names = list(geoms)
    overlaps = [[a, b] for i, a in enumerate(names) for b in names[i + 1:] if geoms[a].intersection(geoms[b]).area > 1e-6]
    feats = [{"type": "Feature", "geometry": mapping(g), "properties": {name_col: n, "slide": slide}} for n, g in geoms.items()]
    feats += [{"type": "Feature", "geometry": mapping(_polygonal(shape_geom(r))), "properties": {name_col: nm(r), "slide": slide, "keep": False}} for r in off]
    fc = {"type": "FeatureCollection", "coordinate_unit": "um (Xenium frame, y down)", "name_column": name_col, "features": feats}
    f = Path(folder).expanduser() / regions_filename(slide)
    f.write_text(json.dumps(fc))
    return str(f), overlaps, empty


def convert(sess, spec, out, project_name, ann, log, transcripts=False, he=None, tx_table=None, threads=1, export_regions=False, regions_name_col="sample_id_unblinded", regions_dir=None):
    """Write the project folder `out`; returns a short report. Refuses to touch a folder that is not a previous import of this kind."""
    out = claim_folder(out, MARKER, "WebbScope Xenium import")
    order, idx, comps, counts = sess.labels(spec)
    kept = np.flatnonzero(idx >= 0)
    if not len(kept):
        raise importer.ImportError_("no cells are selected: every tissue is switched off")
    for nme, c in counts.items():
        log(f"tissue {nme}: {c:,} cells")
    log(f"{len(kept):,} of {sess.n:,} cells kept; {int((idx == -1).sum()):,} dropped by tissue selection, {int((idx == -2).sum()):,} below the minimum transcript count")
    ids = np.array(sess.ids, dtype=object)[kept]

    slide = (project_name or "").strip() or f"{sess.meta.get('run_name', 'Xenium run')} / {sess.meta.get('region_name', '')}".strip(" /")
    n_reg = write_regions(sess, spec, order, idx, out, log, regions_name_col, slide, regions_dir) if export_regions else None
    log("writing counts")
    m = sess.m[:, kept].tocoo()
    m = sp.coo_matrix((np.round(m.data).astype(np.int64), (m.row, m.col)), shape=m.shape)
    with gzip.open(out / "matrix.mtx.gz", "wb", compresslevel=3) as f:
        mmwrite(f, m, field="integer")
    with gzip.open(out / "features.tsv.gz", "wt") as f:
        f.write("\n".join(sess.genes) + "\n")
    with gzip.open(out / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(ids) + "\n")

    log("writing cell table")
    ctype, cstate, extra = sess.cell_labels(ann)
    c = sess.cells.iloc[kept].reset_index(drop=True)
    keep_cols = [k for k in ("cell_id", "x_centroid", "y_centroid", "transcript_counts", "cell_area", "nucleus_area", "nucleus_count", "segmentation_method")
                 if k in c.columns]
    t = c[keep_cols].copy()
    for typ, key in CTRL_KEYS.items():                     # control-feature counts per cell (not part of the gene matrix)
        if typ in sess.ctrl:
            t[f"{key}_counts"] = sess.ctrl[typ][kept]
    t["tissue"] = np.array(order, dtype=object)[idx[kept]]
    t["cell_type"] = ctype[kept]
    t["cell_state"] = cstate[kept]
    if sess.clusters is not None:
        t["cluster_10x"] = sess.clusters.map(lambda v: np.nan if pd.isna(v) else f"Cluster {int(v)}").values[kept]
    for col in extra.columns:
        t[col] = extra[col].values[kept]
    t.to_parquet(out / "cells.parquet", index=False)

    log("writing cell boundaries")
    key = "bnd_parquet" if "bnd_parquet" in sess.src else "bnd_csv"
    bd = pd.read_parquet(sess.src[key], columns=["cell_id", "vertex_x", "vertex_y"]) if key == "bnd_parquet" else \
        pd.read_csv(sess.src[key], usecols=["cell_id", "vertex_x", "vertex_y"])
    bd["cell_id"] = bd["cell_id"].astype(str)
    bd = bd[bd["cell_id"].isin(set(ids))]
    bd.to_parquet(out / "cell_boundaries.parquet", index=False)
    n_nobnd = int(len(ids) - bd["cell_id"].nunique())

    emb_note = None
    if sess.emb:
        e = {"cell_id": ids}
        okall = True
        if "umap" in sess.emb:
            u = sess.emb["umap"].iloc[kept]
            okall &= not u.isna().any().any()
            e["umap_harmony_1"], e["umap_harmony_2"] = u.iloc[:, 0].values, u.iloc[:, 1].values
        if "pca" in sess.emb:
            p = sess.emb["pca"].iloc[kept]
            okall &= not p.isna().any().any()
            for k in range(p.shape[1]):
                e[f"pca_harmony_{k + 1}"] = p.iloc[:, k].values
        if okall:
            pd.DataFrame(e).to_parquet(out / "embeddings.parquet", index=False)
            emb_note = "10x Xenium analysis (not batch-corrected)"
            log("embeddings: 10x UMAP" + (f" and {sess.emb['pca'].shape[1]} PCs" if "pca" in sess.emb else ""))
        else:
            log("the 10x embeddings do not cover all kept cells; embeddings were not written (raise the minimum transcripts per cell to match the 10x analysis)")

    n_tx = n_tab = None
    if (transcripts or tx_table is not None) and "tx" in sess.src:
        n_tx, n_tab = write_transcripts(sess, order, idx, ids, sess.src["tx"], out, log, store=transcripts, table_prefix=tx_table, threads=threads)
    if he:
        from . import he_slide
        he_slide.build(he["czi"], out, he.get("annotations"), log, workers=int(he.get("threads") or 4))
    n_if = None
    ifps = if_image.pick_all(sess.src)
    if ifps:
        try:
            ent = if_image.describe_many(ifps, out)
            if_image.write_index(out, {IF_KEY: ent}, {t: IF_KEY for t in counts})
            n_if = ", ".join(c["name"] for c in ent["channels"])
            log(f"IF image: {len(ent['channels'])} channel(s) ({n_if}; {ent['w']:,} x {ent['h']:,} px, {ent['px']} um/px) are read in place, not copied; kept in {ifps[0].parent}")
        except if_image.IFError as e:
            log(f"IF image ignored: {e}")
    meta = sess.meta
    name = slide
    manifest = {"name": name, "format_version": 2,
                "roles": {"cell_id": "cell_id", "x": "x_centroid", "y": "y_centroid", "sample": "tissue", "cell_type": "cell_type", "cell_state": "cell_state"},
                "norm_scale": 10000.0, "coordinate_unit": "um", "source": "xenium",
                "xenium": {k: meta[k] for k in ("run_name", "region_name", "slide_id", "panel_name", "panel_organism", "panel_tissue_type",
                                                 "preservation_method", "pixel_size", "analysis_sw_version") if k in meta},
                "controls": {"": {CTRL_KEYS[k]: n for k, n in sess.ctrl_n.items() if k in CTRL_KEYS}},      # number of control features per type (run "" = this single run)
                "cell_labels": "annotation file" if (ann and ann.get("path")) else ("10x graph clusters (stand-in for cell types)" if sess.clusters is not None else "none")}
    if emb_note:
        manifest["embedding_label"] = emb_note
    (out / "project.json").write_text(json.dumps(manifest, indent=2))
    write_signatures(out, sess.genes, log)
    (out / MARKER).write_text(json.dumps({"spec": spec, "tissues": counts, "n_cells": int(len(kept)), "n_genes": sess.n_genes,
                                          "min_counts_per_cell": int(spec.get("min_counts", 1)), "cells_without_boundary": n_nobnd}, indent=1))
    return {"n_cells": int(len(kept)), "n_genes": sess.n_genes, "tissues": counts, "cells_without_boundary": n_nobnd, "transcripts": n_tx, "tissue_table": n_tab, "if": n_if, "regions_exported": n_reg}
