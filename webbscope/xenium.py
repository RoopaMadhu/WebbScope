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
IF_KEY = "run"                                  # name of the IF image in a single run's if.json
CTRL_KEYS = {"Negative Control Probe": "nc_probe", "Negative Control Codeword": "nc_codeword", "Unassigned Codeword": "unassigned_codeword", "Deprecated Codeword": "deprecated_codeword"}
BIN_UM = 50.0                                   # density grid for tissue detection
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
        src.pop("tx", None)
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
        self._auto = {}
        self.xy32 = np.column_stack([self.x, self.y]).astype(np.float32)
        log(f"{self.n:,} cells x {self.n_genes:,} genes ready")

    # ---- tissue detection --------------------------------------------------------------------------------------------------------
    def auto(self, gap, min_cells, min_counts):
        """Connected components of the cell-density map. Returns (comp id per cell: -1 = below min_counts, 0 = isolated cell,
        1..K = pieces in the order they would be named, sorted by x), and a list of piece descriptions."""
        key = (round(float(gap), 3), int(min_cells), int(min_counts))
        if key in self._auto:
            return self._auto[key]
        ok = self.gene_counts >= min_counts
        ix = ((self.x - self.x.min()) // BIN_UM).astype(np.int64)
        iy = ((self.y - self.y.min()) // BIN_UM).astype(np.int64)
        nx, ny = ix.max() + 1, iy.max() + 1
        H = np.bincount((ix * ny + iy)[ok], minlength=nx * ny).reshape(nx, ny)
        occ = H >= 2
        it = max(0, int(math.ceil(gap / (2 * BIN_UM))))
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
            comps.append({"id": new, "n": int(cnt[old]), "name": f"tissue_{nb}" if is_big else f"small_piece_{new - len(big)}", "keep": bool(is_big)})
        comp = np.where(ok, remap[cell_lab], -1)
        n0 = int(((comp == 0)).sum())
        if n0:
            comps.insert(0, {"id": 0, "n": n0, "name": "isolated_cells", "keep": False})
        self._auto[key] = (comp, comps)
        return self._auto[key]

    def labels(self, spec):
        """Apply the user's edits. Returns (names, idx per cell (-1 = dropped, -2 = below min_counts), pieces, per-name counts)."""
        gap, mc, mn = float(spec.get("gap", 200)), int(spec.get("min_cells", 500)), int(spec.get("min_counts", 1))
        comp, comps = self.auto(gap, mc, mn)
        edits = spec.get("comps", {}) or {}
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
        """Cells inside a drawn shape: a rectangle {x0,y0,x1,y1} or a free-hand polygon {pts: [[x, y], ...]}."""
        if r.get("pts"):
            import shapely
            from shapely.geometry import Polygon
            from shapely.validation import make_valid
            pts = np.asarray(r["pts"], float)
            if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 3 or len(pts) > 5000 or not np.isfinite(pts).all():
                raise importer.ImportError_("a drawn shape needs 3 to 5000 valid points")
            g = make_valid(Polygon(pts))
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


def _read_annotation(p):
    n = p.name
    if n.endswith(".parquet"):
        return pd.read_parquet(p)
    return pd.read_csv(p, sep=None, engine="python")           # comma or tab


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


def convert(sess, spec, out, project_name, ann, log, transcripts=False, he=None, tx_table=None, threads=1):
    """Write the project folder `out`; returns a short report. Refuses to touch a folder that is not a previous import of this kind."""
    out = Path(out).expanduser()
    if out.exists():
        if not (out / MARKER).exists() and any(out.iterdir()):
            raise importer.ImportError_(f"{out} exists and is not an earlier WebbScope Xenium import; choose a new or empty folder")
        for f in out.iterdir():
            shutil.rmtree(f) if f.is_dir() else f.unlink()
    out.mkdir(parents=True, exist_ok=True)
    order, idx, comps, counts = sess.labels(spec)
    kept = np.flatnonzero(idx >= 0)
    if not len(kept):
        raise importer.ImportError_("no cells are selected: every tissue is switched off")
    for nme, c in counts.items():
        log(f"tissue {nme}: {c:,} cells")
    log(f"{len(kept):,} of {sess.n:,} cells kept; {int((idx == -1).sum()):,} dropped by tissue selection, {int((idx == -2).sum()):,} below the minimum transcript count")
    ids = np.array(sess.ids, dtype=object)[kept]

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
        he_slide.build(he["czi"], out, he.get("annotations"), log)
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
    name = (project_name or "").strip() or f"{meta.get('run_name', 'Xenium run')} / {meta.get('region_name', '')}".strip(" /")
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
    return {"n_cells": int(len(kept)), "n_genes": sess.n_genes, "tissues": counts, "cells_without_boundary": n_nobnd, "transcripts": n_tx, "tissue_table": n_tab, "if": n_if}
