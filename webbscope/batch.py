"""(A folder of already-converted project folders is merged without converting again: kind 'project'.)
A folder of Xenium runs -> one project per run, plus one merged project holding every tissue of every run.

Each run goes through the same code as the single-run import (xenium.Session / xenium.convert) with automatic tissue detection and the
default settings, one run at a time, so only one run's counts are in memory while converting. The merged project is a normal project folder
(`<out>/_merged`) that the viewer opens like any other:
  sample         = "<run folder>_<tissue>"         (tissues are T0 .. Tn-1, ordered left to right on the slide)
  run            = the run folder name            (the Slide window groups the tissues of one run; not offered as a group for statistics)
  cell_type/state= "<run folder>_<tissue>: Cluster N"   (10x graph clusters of that run; not comparable between runs, hence the prefix)
  genes          = those present in every run (reported if any were dropped)
  embeddings     = each run's own 10x UMAP (embedding_per_run in project.json): the runs share no coordinates, so the viewer shows one run
                   at a time; no PCA. Cells are drawn grey until a colouring is chosen.
"""
import gc, json, os, re, shutil, sys, time, traceback
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from . import importer, xenium, he_slide, if_image

TAR_SUFFIXES = (".tar.gz", ".tgz", ".tar", ".zip")                       # bundle archives (.zarr.zip files are not bundles: xenium.is_archive)
MERGED = "_merged"
HE_DIR = "_he"                                                  # H&E tiles of a batch, one sub-folder per scan file (shared by the runs that use it)
HE_EXT = (".czi", ".tif", ".tiff")                                       # H&E scans: Zeiss CZI, or OME-TIFF (.ome.tif / .ome.tiff)
SPEC = {"gap": 200, "min_cells": 500, "min_counts": 1, "comps": {}, "rects": []}


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_") or "run"


def _peak_gb():
    import resource
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3 if sys.platform == "darwin" else 1024 ** 2)      # bytes on macOS, KiB on Linux


def cpu_count():
    try:
        return len(os.sched_getaffinity(0))                       # the cores this process may use (respects a scheduler's allocation)
    except AttributeError:
        return os.cpu_count() or 1


def mem_available_gb():
    """Memory this process can still use: the smaller of the container / job limit (cgroup) and what the system reports free."""
    lim = []
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                lim.append(int(line.split()[1]) / 1024 ** 2)
    except OSError:
        pass
    for lf, uf in (("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
                   ("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes")):
        try:
            v = open(lf).read().strip()
            if v.isdigit() and int(v) < 1 << 50:
                lim.append((int(v) - int(open(uf).read().strip())) / 1024 ** 3)
        except (OSError, ValueError):
            pass
    if not lim:
        try:
            lim.append(os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 1024 ** 3 * 0.6)      # no "available" figure (macOS): 60% of RAM
        except (ValueError, OSError, AttributeError):
            return None
    return min(lim)


GB_PER_RUN = {True: 20.0, False: 6.0}                             # peak memory of one run being converted, with / without transcripts (seen: 11-19 GB with)
AUTO_CORES = 4


def plan_workers(requested, n_runs, with_tx, log):
    """(runs converted at once, cores for each). requested 0 = automatic (up to AUTO_CORES cores). Never more than the runs, the cores or what
    the memory allows (about GB_PER_RUN per run at once)."""
    cores = cpu_count()
    want = int(requested) if requested and int(requested) > 0 else min(AUTO_CORES, cores)
    w = max(1, min(want, cores, max(n_runs, 1)))
    avail = mem_available_gb()
    if avail is not None and w > 1:
        fit = max(1, int(avail * 0.85 / GB_PER_RUN[bool(with_tx)]))
        if fit < w:
            log(f"parallel runs limited to {fit} (asked for {w}): about {GB_PER_RUN[bool(with_tx)]:.0f} GB per run and {avail:.0f} GB of memory available")
            w = fit
    log(f"{cores} core(s) available: converting {w} run(s) at a time, {max(1, cores // w)} core(s) for each run's nearest-cell search")
    return w, max(1, cores // w)


def _complete(found):
    return all(any(r in found for r in g) for g in xenium.REQUIRED_GROUPS if g != ("experiment",))


def existing_merged(folder):
    """The merged project an earlier batch import wrote to `folder`/_merged, or None: {path, runs, n_tissues, n_cells, n_genes, transcripts}."""
    d = Path(folder).expanduser() / MERGED
    mk = d / xenium.MARKER
    if not ((d / "project.json").is_file() and mk.is_file()):
        return None
    try:
        m = json.loads(mk.read_text())
    except Exception:
        return None
    if "runs" not in m:
        return None
    return {"path": str(d), "runs": list(m["runs"]), "n_tissues": len(m.get("tissues") or {}), "n_cells": m.get("n_cells"), "n_genes": m.get("n_genes"), "transcripts": m.get("transcripts")}


def _is_project(p):
    """A project folder written by this app's Xenium import (manifest source 'xenium'): returns True / False / None (not a project)."""
    m = Path(p) / "project.json"
    if not (m.is_file() and (Path(p) / "cells.parquet").is_file()):
        return None
    if (Path(p) / xenium.MARKER).is_file() and not xenium.finished(p):
        return None                                             # an import that was stopped half-way: not a project
    try:
        return json.loads(m.read_text()).get("source") == "xenium"
    except Exception:
        return None


def discover(parent):
    """Runs found directly inside `parent`: extracted Xenium output folders, .tar.gz bundles, and project folders already written by this app's
    import (kind 'project': merged without converting again). Returns {runs: [{name, path, kind}], skipped: [...]}."""
    parent = Path(parent).expanduser()
    if not parent.is_dir():
        raise importer.ImportError_(f"not a folder: {parent}")
    dirs, tars, projs, skipped = [], [], [], []
    for p in sorted(parent.iterdir(), key=lambda q: q.name):
        if p.name.startswith("."):
            continue
        if p.is_dir() and p.name in (MERGED, HE_DIR):           # an earlier merge: offered for opening (existing_merged); _he = tiles written by this import; neither is a run
            continue
        if p.is_dir():
            isp = _is_project(p)
            if isp:
                projs.append(p)
            elif isp is False or p.name == MERGED:               # a merged project, or one of another source: not merged again
                skipped.append(p)
            else:
                (dirs if _complete(xenium.find_in_folder(p)) else skipped).append(p)
        elif p.is_file() and xenium.is_archive(p.name):
            tars.append(p)
    runs = []
    for p in dirs:
        runs.append({"name": p.name, "path": str(p), "kind": "folder"})
    for p in projs:
        runs.append({"name": p.name, "path": str(p), "kind": "project"})
    for p in tars:
        runs.append({"name": re.sub(r"(\.tar\.gz|\.tgz|\.tar|\.zip)$", "", p.name, flags=re.I), "path": str(p), "kind": "zip" if p.name.lower().endswith(".zip") else "tar.gz"})
    seen = {}
    for r in runs:                                              # run names become file names and sample prefixes: they must be unique
        k = safe_name(r["name"])
        seen[k] = seen.get(k, 0) + 1
        if seen[k] > 1:
            r["name"] = f"{r['name']}_{seen[k]}"
    return {"parent": str(parent), "runs": runs, "skipped": [p.name for p in skipped if p.is_dir()], "merged": existing_merged(parent)}


# ---------------------------------------------------------------- H&E folder: match scans to runs by name, build the tiles once per scan
def he_discover(folder, depth=2):
    """H&E scans (.czi, .ome.tif, .ome.tiff) and ImageScope annotation files (.annotations) in `folder` and its sub-folders (up to `depth` levels)."""
    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise importer.ImportError_(f"not a folder: {folder}")
    czi, ann = [], []

    def walk(d, lvl):
        try:
            items = sorted(d.iterdir(), key=lambda q: q.name)
        except OSError:
            return
        for p in items:
            if p.name.startswith("."):
                continue
            if p.is_dir():
                if lvl < depth:
                    walk(p, lvl + 1)
            elif p.suffix.lower() in HE_EXT:
                czi.append({"name": p.name, "path": str(p), "size": p.stat().st_size})
            elif p.suffix.lower() == ".annotations":
                ann.append({"name": p.name, "path": str(p)})
    walk(folder, 0)
    return {"folder": str(folder), "czi": czi, "annotations": ann}


def _name_tokens(name):
    """Identifying parts of a file / run name: every run of digits (one-digit runs, as in Region_1, are ignored) with the letters just before it,
    separators dropped, so 'chip0073' and 'Chip-0073' agree. List of (letters, digits without leading zeros, digits as written)."""
    stem = re.sub(r"(\.tar\.gz|\.tgz|\.tar|\.zip|\.czi|\.ome\.tiff?|\.tiff?|\.annotations)$", "", name, flags=re.I).lower()
    out = []
    for m in re.finditer(r"\d+", stem):
        if len(m.group()) < 2:
            continue
        before = re.sub(r"[^a-z]", "", re.sub(r"\d+", "|", stem[:m.start()]).split("|")[-1])          # letters since the previous digit run
        out.append((before, m.group().lstrip("0") or "0", m.group()))
    return out


def _suffix(a, b):
    n = 0
    while n < min(len(a), len(b)) and a[-1 - n] == b[-1 - n]:
        n += 1
    return n


def _shared(a, b):
    """Tokens two names share: the same digits, and letter contexts that end alike (equal, one a tail of the other, or sharing a tail of 4+ letters); a bare number only matches
    a bare number or one with at least 4 digits. Returns (weight, the shared tokens as written in a, their digit strings)."""
    w, shown, digs = 0, [], set()
    for la, da, oa in a:
        for lb, db, _ in b:
            if da != db:
                continue
            k = _suffix(la, lb)
            if (la == lb) or (k and (k == min(len(la), len(lb)) or k >= 4)) or (not (la and lb) and len(da) >= 4):
                w += len(da) + k
                shown.append(la[-k:] + oa if k else oa)
                digs.add(da)
                break
    return w, sorted(shown), digs


def he_match(run_names, czi, ann):
    """Suggest an H&E scan (and annotation file) for each run from the file names. Returns {run: {czi, annotations, why, tied}} where czi / annotations
    are paths or None. A run gets a scan when exactly one scan shares the most identifying numbers with it (numbers common to all runs are ignored);
    a tie or no shared number leaves it empty for the user to choose. Several runs may share one scan (one slide holding two regions). The annotation
    file must contain one of the numbers that tied the run to its scan (a shared date alone is not enough)."""
    rt = {r: _name_tokens(r) for r in run_names}
    if len(rt) > 1:
        common = set.intersection(*[{t[1] for t in v} for v in rt.values()])
        rt = {r: [t for t in v if t[1] not in common] for r, v in rt.items()}
    ct = {c["path"]: _name_tokens(c["name"]) for c in czi}
    at = {a["path"]: _name_tokens(a["name"]) for a in ann}
    out = {}
    for r in run_names:
        sc = sorted(((_shared(rt[r], t)[0], pth) for pth, t in ct.items()), reverse=True)
        best = [pth for v, pth in sc if v and v == sc[0][0]]
        rec = {"czi": None, "annotations": None, "why": "", "tied": []}
        if len(best) == 1:
            rec["czi"] = best[0]
            _, shown, digs = _shared(rt[r], ct[best[0]])
            rec["why"] = "shares " + ", ".join(shown)
            cand = [(_shared(ct[best[0]], t)[0], pth) for pth, t in at.items() if _shared(rt[r], t)[2] & digs]
            cand.sort(reverse=True)
            top = [pth for v, pth in cand if v and v == cand[0][0]] if cand else []
            if len(top) == 1:
                rec["annotations"] = top[0]
        elif len(best) > 1:
            rec["tied"] = best
            rec["why"] = f"{len(best)} scans match equally, choose one"
        else:
            rec["why"] = "no scan name shares a number with this run"
        out[r] = rec
    return out


# ---- cell labels from the user's own annotation file (one table for all runs), applied when the runs are merged
ID_PREFIX = "__id_prefix__"                                     # run column choice: the part of the cell id before its last _ | : (as Seurat / Scanpy prefix merged cells)
XENIUM_ID = re.compile(r"^[a-p]{8}-\d+$|^\d+$")                 # a Xenium cell id (current and early formats)
_TAIL = re.compile(r"^.*[_|:]")


def _id_parts(v):
    """(cell id as written, its part after the last _ | :, the part before it or "")."""
    v = pd.Series(v, dtype="string").fillna("")
    tail = v.str.replace(_TAIL, "", regex=True)
    head = pd.Series([a[: len(a) - len(b) - 1] if a != b else "" for a, b in zip(v, tail)], dtype="string")
    return v, tail, head


def ann_inspect(path, run_names):
    """The columns of an annotation file and the suggested roles: {columns, n_rows, id_col, run_col, type_col, state_col, donor_col, id_prefix}.
    id_prefix = the cell ids carry a prefix (sample_cellid), which can stand for the run when there is no run column."""
    from .engine.config import CANDIDATES, find_columns
    p = Path(path).expanduser()
    if not p.is_file():
        raise importer.ImportError_(f"not found: {p}")
    d = xenium._read_annotation(p)
    if not len(d.columns):
        raise importer.ImportError_("the annotation file has no columns")
    cols = [str(c) for c in d.columns]
    smp = d.head(5000)
    def xen(c):                                                 # share of values that are Xenium cell ids, as written or after a prefix
        v, tail, _ = _id_parts(smp[c].astype(str))
        return max(v.str.match(XENIUM_ID).mean(), tail.str.match(XENIUM_ID).mean())
    rate = {c: float(xen(c)) for c in d.columns if d[c].dtype == object or str(d[c].dtype).startswith(("string", "category"))}
    id_col = max(rate, key=rate.get) if rate and max(rate.values()) > 0.5 else next((c for c in CANDIDATES["cell_id"] if c in cols), cols[0])
    v, tail, head = _id_parts(smp[id_col].astype(str))
    id_prefix = bool((head != "").mean() > 0.9 and tail.str.match(XENIUM_ID).mean() > 0.5)
    run_col, best = None, 0
    for c in d.columns:
        if c == id_col:
            continue
        vals = d[c].dropna().astype(str).unique()
        if not 1 <= len(vals) <= 2000:
            continue
        m = ann_match(vals, run_names)
        k = sum(1 for x in m.values() if x["value"] is not None)
        if k > best:
            run_col, best = c, k
    if run_col is None and id_prefix:
        run_col = ID_PREFIX
    pick = lambda role: next((c for c in CANDIDATES[role] if c in cols), None)
    donor = pick("donor") or next(iter(find_columns(cols, "donor")), None)
    return {"columns": cols, "n_rows": int(len(d)), "id_col": id_col, "run_col": run_col, "type_col": pick("cell_type"), "state_col": pick("cell_state"),
            "donor_col": donor, "id_prefix": id_prefix}


def _ann_values(d, run_col, id_col):
    """The value saying which run each row is from (the cell id's prefix for ID_PREFIX)."""
    if run_col == ID_PREFIX:
        return _id_parts(d[id_col].astype(str))[2]
    return d[run_col].astype("string")


def ann_values(path, run_col, id_col):
    """{value: number of rows} of the run column (at most 2000 values)."""
    d = xenium._read_annotation(Path(path).expanduser())
    for c in (id_col, run_col):
        if c and c != ID_PREFIX and c not in d.columns:
            raise importer.ImportError_(f"the annotation file has no column {c!r}")
    vc = _ann_values(d, run_col, id_col).dropna().value_counts()
    if len(vc) > 2000:
        raise importer.ImportError_(f"the column {run_col} has {len(vc):,} different values: it does not say which run a cell is from")
    return {str(k): int(n) for k, n in vc.items()}


def ann_match(values, run_names):
    """Suggest the run-column value of each run: {run: {value, why}}. The value equal to the run name, else the one value sharing the most
    identifying numbers with it (as for the H&E scans: 'slidearray077' matches '..._slidearray0077_...'); a tie or no shared number leaves it empty."""
    values = [str(v) for v in values]
    rt = {r: _name_tokens(r) for r in run_names}
    if len(rt) > 1:
        common = set.intersection(*[{t[1] for t in v} for v in rt.values()])
        rt = {r: [t for t in v if t[1] not in common] for r, v in rt.items()}
    vt = {v: _name_tokens(v) for v in values}
    out = {}
    for r in run_names:
        if r in vt:
            out[r] = {"value": r, "why": "same name"}; continue
        sc = sorted(((_shared(rt[r], t)[0], v) for v, t in vt.items()), reverse=True)
        best = [v for w, v in sc if w and w == sc[0][0]]
        if len(best) == 1:
            out[r] = {"value": best[0], "why": "shares " + ", ".join(_shared(rt[r], vt[best[0]])[1])}
        else:
            out[r] = {"value": None, "why": f"{len(best)} values match equally, choose one" if best else "no value shares a number with this run"}
    return out


def ann_prepare(ann, log):
    """Read the annotation file once for the merge: {spec, rows: {value: rows}} (rows for one run) or {spec, rows: None} (one run, no run column)."""
    d = xenium._read_annotation(Path(ann["path"]).expanduser())
    for k in ("id_col", "run_col", "type_col", "state_col", "donor_col"):
        c = ann.get(k)
        if c and c != ID_PREFIX and c not in d.columns:
            raise importer.ImportError_(f"the annotation file has no column {c!r} (chosen as {k.replace('_col', '').replace('_', ' ')})")
    if not (ann.get("type_col") or ann.get("state_col")):
        raise importer.ImportError_("choose the cell type column (or the cell state column) of the annotation file")
    log(f"annotation file {Path(ann['path']).name}: {len(d):,} rows")
    if not ann.get("run_col"):
        return {"spec": ann, "rows": None, "df": d}
    key = _ann_values(d, ann["run_col"], ann["id_col"])
    return {"spec": ann, "rows": {str(k): g for k, g in d.groupby(key.values, sort=False)}}


def _annotate(t, run, A, log):
    """Cell type, cell state and the other columns of the annotation rows for this run, matched by cell id (as written, or after a prefix).
    Returns (table, True) or (t unchanged, False) when the run has no rows or no cell matches."""
    s = A["spec"]
    if A["rows"] is None:
        rows = A["df"]
    else:
        val = (s.get("runs") or {}).get(run)
        rows = A["rows"].get(str(val)) if val is not None else None
        if rows is None:
            log(f"{run}: no annotation rows chosen; its 10x clusters are kept")
            return t, False
    ids = t["cell_id_in_run"].astype(str)
    v, tail, _ = _id_parts(rows[s["id_col"]].astype(str))
    key = tail if tail.isin(ids).sum() > v.isin(ids).sum() else v
    r = rows.assign(_k=key.values)
    dup = int(r["_k"].duplicated().sum())
    r = r.drop_duplicates("_k").set_index("_k")
    hit = ids.isin(r.index)
    if not hit.any():
        log(f"{run}: none of its cell ids is in the annotation rows chosen for it; its 10x clusters are kept")
        return t, False
    r = r.reindex(ids.values)
    t = t.copy()
    tc, sc = s.get("type_col") or s.get("state_col"), s.get("state_col") or s.get("type_col")
    t["cell_type"] = r[tc].astype("string").fillna("unlabelled").values
    t["cell_state"] = r[sc].astype("string").fillna("unlabelled").values
    skip = {s["id_col"], s.get("run_col"), s.get("type_col"), s.get("state_col")}
    clash = [c for c in r.columns if c not in skip and c in t.columns]
    for c in r.columns:
        if c not in skip and c not in t.columns:
            t[c] = r[c].values
    log(f"{run}: {int(hit.sum()):,} of {len(t):,} cells labelled from the annotation file"
        + (f"; {dup:,} repeated cell ids (first kept)" if dup else "") + (f"; columns already in the cell table not taken: {', '.join(map(str, clash[:5]))}" if clash else ""))
    return t, True


def _he_key(czi_path):
    return safe_name(Path(czi_path).stem)


def _fp(p):
    st = Path(p).stat()
    return [Path(p).name, st.st_size, int(st.st_mtime)]


def build_he(czi_path, ann_path, out_root, log, workers=4):
    """Tiles of one scan into <out_root>/_he/<key>/ (he/, he_slide.json, annotations). Skipped when that folder already holds this exact file."""
    key = _he_key(czi_path)
    d = Path(out_root) / HE_DIR / key
    mk = d / ".source.json"
    want = {"czi": _fp(czi_path), "ann": _fp(ann_path) if ann_path else None}
    if not he_slide.is_ome(czi_path):
        want["levels"] = "halved from full resolution"           # tiles of a CZI made before this read the coarser levels from the file's own copies: made again
    if mk.is_file() and (d / he_slide.INDEX).is_file():
        try:
            if json.loads(mk.read_text()) == want:
                log(f"H&E {Path(czi_path).name}: tiles from an earlier import are current, not rebuilt")
                return key, json.loads((d / he_slide.INDEX).read_text())
        except Exception:
            pass
    d.mkdir(parents=True, exist_ok=True)
    if mk.exists():
        mk.unlink()
    he_slide.build(czi_path, d, ann_path, log=log, workers=workers)
    mk.write_text(json.dumps(want))
    return key, json.loads((d / he_slide.INDEX).read_text())


def attach_he(out_root, merged, he_map, run_names, log, cores=1, threads=4):
    """Build the tiles of every scan that a converted run uses and write the merged project's he_slide.json: one entry per scan with the runs that use it.
    he_map = {run: {czi, annotations}}. A scan that cannot be read is logged and skipped (its runs get no H&E).
    Each scan is tiled with `threads` threads; scans are built side by side when the cores allow (cores // threads at once)."""
    out_root, merged = Path(out_root), Path(merged)
    scans, runs = {}, {}
    todo = {}
    for r in run_names:
        m = (he_map or {}).get(r) or {}
        if m.get("czi"):
            todo.setdefault((m["czi"], m.get("annotations") or None), []).append(r)
    if not todo:
        return None
    log(f"=== H&E: {len(todo)} scan(s) for {sum(len(v) for v in todo.values())} run(s)")
    items = list(todo.items())
    threads = max(1, int(threads))
    n_par = max(1, min(len(items), cores // threads))

    def one(i):
        (czi_path, ann_path), rs = items[i]
        tag = f"[{Path(czi_path).name}] " if n_par > 1 else ""
        log(f"=== H&E [{i + 1}/{len(items)}] {Path(czi_path).name} -> {', '.join(rs)}")
        try:
            return build_he(czi_path, ann_path, out_root, (lambda m: log(tag + m)), workers=threads)
        except Exception as e:
            if not isinstance(e, importer.ImportError_):
                traceback.print_exc()
            log(f"=== H&E {Path(czi_path).name} FAILED, skipped: {e}")
            return None
    log(f"H&E: {n_par} scan(s) at a time, {threads} thread(s) each ({cores} core(s))")
    if n_par > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(n_par) as ex:
            res = list(ex.map(one, range(len(items))))
    else:
        res = [one(i) for i in range(len(items))]
    for ((czi_path, ann_path), rs), r in zip(items, res):
        if r is None:
            continue
        key, idx = r
        if key in scans:                                       # two entries with the same scan but another annotation file: the first wins
            scans[key]["runs"] += rs
        else:
            scans[key] = {"dir": os.path.relpath((out_root / HE_DIR / key).resolve(), merged.resolve()), "source": idx["source"],
                          "n_scenes": len(idx["scenes"]), "annotations": bool(idx.get("annotations")), "runs": list(rs)}
        for r_ in rs:
            runs[r_] = key
    if not scans:
        return None
    (merged / he_slide.INDEX).write_text(json.dumps({"multi": True, "scans": scans, "runs": runs}, indent=1))
    log(f"H&E: {len(scans)} scan(s) attached to {len(runs)} run(s) of the merged project (tiles stay in {out_root / HE_DIR}; keep the output folder together)")
    return scans


def _convert_one(run, out, spec, log, transcripts=False, tx_table=False, threads=1, ifimg=False, reg=None):
    """Write the run's own project folder; return what the merge needs (counts of the kept cells, genes, cell table, boundaries)."""
    src = xenium.resolve_sources(run["path"], {}, log, want_tx=transcripts or tx_table or spec.get("detect") == "transcripts", want_if=ifimg)
    if ifimg and not any(k.startswith("if_") for k in src):
        log(f"{run['name']}: no morphology_focus / morphology_mip .ome.tif found; this run will have no IF image")
    if (transcripts or tx_table) and "tx" not in src:
        log(f"{run['name']}: no transcripts.parquet found; this run will have no transcript points" + (" and no transcript_tissue.parquet" if tx_table else ""))
    reg = reg or {}
    # shapes drawn / tissues renamed for this run (none = automatic tissues); auto_off: tissue boxes loaded from a GeoJSON, the detected pieces start switched off
    spec = {**spec, "rects": reg.get("rects", {}).get(run["name"], []), "comps": reg.get("comps", {}).get(run["name"], {}),
            "auto_off": bool(reg.get("auto_off", {}).get(run["name"])), "cut_off": bool(reg.get("cut_off", True))}
    sess = xenium.Session(src, log)
    _, idx, _, _ = sess.labels(spec)
    kept = np.flatnonzero(idx >= 0)
    rep = xenium.convert(sess, spec, out, run["name"], None, log, transcripts, None, run["name"] if tx_table else None, threads,
                        export_regions=bool(reg.get("export")), regions_name_col=reg.get("col") or "sample_id_unblinded", regions_dir=reg.get("dir"))
    m = sess.m[:, kept].tocsr()
    m.data = np.round(m.data).astype(np.int32)
    sess_ctrl_n = dict(sess.ctrl_n)
    piece = {"name": run["name"], "dir": Path(out), "m": m, "genes": list(sess.genes), "tissues": rep["tissues"], "n_cells": rep["n_cells"],
             "controls": {xenium.CTRL_KEYS[k]: n for k, n in sess_ctrl_n.items() if k in xenium.CTRL_KEYS},
             "tx_dir": (out / "transcripts") if rep.get("transcripts") is not None else None,
             "if": if_image.read_index(out)}
    del sess
    piece["cells"] = pd.read_parquet(out / "cells.parquet")                     # same cells, same order as the columns of m
    piece["bnd"] = pd.read_parquet(out / "cell_boundaries.parquet")
    return piece


def _load_project(run, log, transcripts=False):
    """The merge input of a project folder this app already wrote (no conversion): counts of its cells, genes, cell table, boundaries."""
    d = Path(run["path"])
    man = json.loads((d / "project.json").read_text())
    genes, ids = importer._lines(d / "features.tsv.gz"), importer._lines(d / "barcodes.tsv.gz")
    m = sp.load_npz(d / "matrix.npz") if (d / "matrix.npz").is_file() else __import__("scipy.io", fromlist=["mmread"]).mmread(d / "matrix.mtx.gz")
    m = m.tocsr() if sp.issparse(m) else sp.csr_matrix(m)
    if m.shape == (len(ids), len(genes)) and m.shape != (len(genes), len(ids)):
        m = m.T.tocsr()
    if m.shape != (len(genes), len(ids)):
        raise importer.ImportError_(f"{run['name']}: matrix is {m.shape[0]} x {m.shape[1]} but there are {len(genes)} genes and {len(ids)} cells")
    m.data = np.round(m.data).astype(np.int32)
    cells = pd.read_parquet(d / "cells.parquet")
    if len(cells) != len(ids) or not (cells["cell_id"].astype(str).values == np.asarray(ids)).all():
        raise importer.ImportError_(f"{run['name']}: cells.parquet and barcodes.tsv.gz do not list the same cells in the same order")
    mk = d / xenium.MARKER
    tissues = json.loads(mk.read_text()).get("tissues") if mk.is_file() else None
    if not tissues:
        tissues = {str(k): int(v) for k, v in cells["tissue"].value_counts(sort=False).items()}
    ctrl = (man.get("controls") or {}).get("") or {}
    tx = d / "transcripts"
    has_tx = (tx / "index.json").is_file()
    if has_tx and not transcripts:
        log(f"{run['name']}: has stored transcripts, not included (tick the transcripts box to include them)")
    log(f"{run['name']}: existing project read ({len(ids):,} cells, {len(genes):,} genes, {len(tissues)} tissue(s))")
    return {"name": run["name"], "dir": d, "m": m, "genes": list(genes), "tissues": tissues, "n_cells": len(ids), "controls": ctrl,
            "tx_dir": tx if (has_tx and transcripts) else None, "if": if_image.read_index(d), "cells": cells, "bnd": pd.read_parquet(d / "cell_boundaries.parquet")}


def merge_transcripts(pieces, common, out, log):
    """Write out/transcripts/index.json pointing at each run's own transcripts/ folder (nothing is copied; see transcripts.py). Runs without
    transcripts are simply not listed. The runs' folders must stay where they are (relative paths, so the whole output folder can be moved)."""
    d = Path(out) / "transcripts"
    d.mkdir(parents=True, exist_ok=True)
    parts, n = [], 0
    for p in pieces:
        if not p["tx_dir"]:
            continue
        tis = list(p["tissues"])
        sub = json.loads((Path(p["tx_dir"]) / "index.json").read_text())
        if list(sub["samples"]) != tis:
            raise importer.ImportError_(f"{p['name']}: transcript tissues do not match the converted tissues")
        gi = None if p["genes"] == common else pd.Index(p["genes"]).get_indexer(common).tolist()
        parts.append({"dir": os.path.relpath(Path(p["tx_dir"]).resolve(), d.resolve()), "samples": [f"{p['name']}_{t}" for t in tis], "gene_index": gi})
        n += int(sub["n"])
    if not parts:
        return None
    (d / "index.json").write_text(json.dumps({"n_genes": len(common), "n": n, "parts": parts}))
    log(f"transcripts: {n:,} points in {sum(len(q['samples']) for q in parts)} tissues, read from the {len(parts)} runs' own folders (not copied)")
    return n


def merge_if(pieces, out, log):
    """Write out/if.json: each run's IF image (read where it is, recorded relative to the merged folder as well) and the merged sample names (<run>_<tissue>)."""
    images, samples = {}, {}
    for p in pieces:
        ix = p.get("if")
        if not ix:
            continue
        key = p["name"]
        ent = dict(next(iter(ix["images"].values())))
        ent["channels"] = [dict(c) for c in ent["channels"]]
        for e in [ent] + ent["channels"]:                                # the entry's file and each channel's own file, relative to the merged folder
            if not e.get("file"):
                continue
            try:
                e["rel"] = os.path.relpath(e["file"], Path(out).resolve())
            except ValueError:
                e.pop("rel", None)
        images[key] = ent
        samples.update({f"{p['name']}_{t}": key for t in ix["samples"]})
    if not images:
        return 0
    if_image.write_index(out, images, samples)
    log(f"IF images: {len(images)} run(s) have one (read in place, not copied)")
    return len(images)


def _run_umap(p, ids):
    """The run's own 10x UMAP for its cells ids (in that order) as an (n, 2) array, or None when the run folder has none."""
    f = Path(p["dir"]) / "embeddings.parquet" if p.get("dir") else None
    if f is None or not f.is_file():
        return None
    import pyarrow.parquet as pq
    if not {"cell_id", "umap_harmony_1", "umap_harmony_2"} <= set(pq.read_schema(f).names):
        return None
    e = pd.read_parquet(f, columns=["cell_id", "umap_harmony_1", "umap_harmony_2"])
    e = e.set_index(e["cell_id"].astype(str)).reindex(pd.Index(ids, dtype=str))
    return e[["umap_harmony_1", "umap_harmony_2"]].to_numpy(np.float32)


def merge(pieces, out, log, ann=None):
    """Write the merged project folder `out` from the per-run pieces. ann: the user's annotation file and its columns (cell labels, see _annotate)."""
    if ann and ann.get("path") and not ann.get("run_col") and len(pieces) > 1:
        raise importer.ImportError_("choose the column of the annotation file that says which run each cell is from: Xenium cell ids repeat between runs")
    A = ann_prepare(ann, log) if ann and ann.get("path") else None
    out = xenium.claim_folder(out, xenium.MARKER, "merged import")
    sets = [set(p["genes"]) for p in pieces]
    common = [g for g in pieces[0]["genes"] if all(g in s for s in sets)]
    if not common:
        raise importer.ImportError_("the runs have no gene in common")
    dropped = {p["name"]: len(p["genes"]) - len(common) for p in pieces if len(p["genes"]) != len(common)}
    log(f"merging {len(pieces)} runs: {len(common):,} genes shared by all runs" + (f"; genes left out per run: {dropped}" if dropped else " (the panels are identical)"))
    mats, tabs, bnds, umaps, own = [], [], [], [], []
    for p in pieces:
        run = p["name"]
        if p["genes"] != common:
            m = p["m"][pd.Index(p["genes"]).get_indexer(common)]
        else:
            m = p["m"]
        mats.append(m)
        t = p["cells"].copy()
        t["cell_id_in_run"] = t["cell_id"].astype(str)
        t["cell_id"] = run + "|" + t["cell_id_in_run"]
        sample = run + "_" + t["tissue"].astype(str)
        t["tissue"] = sample
        t["run"] = run
        lab = False
        if A is not None:
            t, lab = _annotate(t, run, A, log)
        if not lab:                                             # 10x clusters: not comparable between runs, so they carry the tissue's name
            for col in ("cell_type", "cell_state"):
                t[col] = sample + ": " + t[col].astype(str)
        own.append(lab)
        tabs.append(t)
        b = p["bnd"].copy()
        b["cell_id"] = run + "|" + b["cell_id"].astype(str)
        bnds.append(b)
        umaps.append(_run_umap(p, t["cell_id_in_run"].tolist()))
    log("writing merged counts")
    M = sp.hstack(mats, format="csr")
    sp.save_npz(out / "matrix.npz", M, compressed=False)
    cells = pd.concat(tabs, ignore_index=True)
    for c in cells.columns:                                     # annotation columns that only some runs have: object columns of mixed values become text
        if cells[c].dtype == object:
            cells[c] = cells[c].where(cells[c].isna(), cells[c].astype(str))
    ids = cells["cell_id"].tolist()
    if len(set(ids)) != len(ids):
        raise importer.ImportError_("duplicate cell ids after merging")
    import gzip
    with gzip.open(out / "features.tsv.gz", "wt") as f:
        f.write("\n".join(common) + "\n")
    with gzip.open(out / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(ids) + "\n")
    log("writing merged cell table and boundaries")
    cells.to_parquet(out / "cells.parquet", index=False)
    pd.concat(bnds, ignore_index=True).to_parquet(out / "cell_boundaries.parquet", index=False)
    have = [u is not None for u in umaps]
    if any(have):                                              # each run's own UMAP; cells of a run without one get NaN
        u = np.concatenate([x if x is not None else np.full((len(t), 2), np.nan, np.float32) for x, t in zip(umaps, tabs)])
        pd.DataFrame({"cell_id": ids, "umap_harmony_1": u[:, 0], "umap_harmony_2": u[:, 1]}).to_parquet(out / "embeddings.parquet", index=False)
        log(f"embeddings: the 10x UMAP of {sum(have)} of {len(pieces)} runs, each in its own coordinates (shown one run at a time)"
            + ("" if all(have) else "; without one: " + ", ".join(p["name"] for p, h in zip(pieces, have) if not h)))
    n_tx = merge_transcripts(pieces, common, out, log) if any(p["tx_dir"] for p in pieces) else None
    merge_if(pieces, out, log)
    tissues = {f"{p['name']}_{k}": v for p in pieces for k, v in p["tissues"].items()}
    manifest = {"name": f"{len(pieces)} Xenium run" + ("s" if len(pieces) != 1 else ""), "format_version": 2,
                "roles": {"cell_id": "cell_id", "x": "x_centroid", "y": "y_centroid", "sample": "tissue", "cell_type": "cell_type", "cell_state": "cell_state"},
                "norm_scale": 10000.0, "coordinate_unit": "um", "source": "xenium-batch", "default_grey": True,
                "runs": [p["name"] for p in pieces], "controls": {p["name"]: p["controls"] for p in pieces},
                "cell_labels": "10x graph clusters per tissue ('<run>_<tissue>: Cluster N'); not comparable between runs"}
    if any(own):
        donor = (ann or {}).get("donor_col")
        if donor and donor in cells.columns:
            manifest["roles"]["donor"] = donor
            d = cells[donor].astype("string")
            if d.isna().any():                                  # a tissue's cells without a donor (not in the file, or a run without labels): the tissue's own donor, else the tissue
                top = d.groupby(cells["tissue"].values).agg(lambda s: s.mode().iloc[0] if s.notna().any() else pd.NA)
                d = d.fillna(pd.Series(cells["tissue"].map(top).values, dtype="string")).fillna(cells["tissue"].astype("string"))
                cells[donor] = d.astype(str)
                cells.to_parquet(out / "cells.parquet", index=False)
                log(f"donor: cells without one take their tissue's donor (or the tissue's name when none of its cells has one)")
        manifest["cell_labels"] = (f"from {Path(ann['path']).name} ({sum(own)} of {len(pieces)} runs)"
                                   + ("; the other runs keep their 10x graph clusters ('<run>_<tissue>: Cluster N')" if not all(own) else ""))
        manifest["default_grey"] = not all(own)
        log(f"cell labels: from the annotation file for {sum(own)} of {len(pieces)} runs" + ("" if all(own) else "; the others keep their 10x clusters"))
    if any(have):
        manifest.update(embedding_label="10x Xenium analysis of each run (not comparable between runs)", embedding_per_run=True)
    (out / "project.json").write_text(json.dumps(manifest, indent=2))
    xenium.write_signatures(out, common, log)
    (out / xenium.MARKER).write_text(json.dumps({"runs": [p["name"] for p in pieces], "tissues": tissues, "n_cells": len(ids), "n_genes": len(common),
                                                 "genes_left_out_per_run": dropped, "transcripts": n_tx}, indent=1))
    log(f"merged project: {len(tissues)} tissues, {len(ids):,} cells, {len(common):,} genes (peak memory so far {_peak_gb():.1f} GB)")
    return out


_Q = None


def _init_worker(q):
    global _Q
    _Q = q


def _worker(i, run, out, spec, transcripts, tx_table, threads, ifimg, reg=None):
    """One run converted in a worker process: (index, piece or None, error text or None). Log lines go to the parent through the queue."""
    def log(m):
        _Q.put((run["name"], m))
    try:
        return i, _convert_one(run, out, spec, log, transcripts, tx_table, threads, ifimg, reg), None
    except Exception as e:
        if not isinstance(e, importer.ImportError_):
            traceback.print_exc()
        return i, None, str(e)


def run_batch(runs, out_root, log, spec=None, do_merge=True, transcripts=False, he=None, tx_table=False, workers=0, ifimg=False, he_reg=None, reg=None, he_threads=4, ann=None):
    """Convert every run (a failing run is logged and skipped), then merge. `workers` = runs converted at once (0 = automatic, 1 = one by one).
    Returns {done, failed, merged, first}."""
    spec = {**SPEC, **(spec or {})}
    out_root = Path(out_root).expanduser()
    out_root.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_conv = sum(1 for r in runs if r.get("kind") != "project")
    w, threads = plan_workers(workers, n_conv, transcripts or tx_table, log)
    got, failed = {}, {}

    def done(i, p):
        got[i] = p
        log(f"=== {runs[i]['name']} done: {p['n_cells']:,} cells in {len(p['tissues'])} tissue(s) ({time.time() - t0:.0f} s elapsed, peak memory {_peak_gb():.1f} GB)")

    def fail(i, e):
        failed[i] = {"name": runs[i]["name"], "error": str(e)}
        log(f"=== {runs[i]['name']} FAILED, skipped: {e}")

    def serial(i):
        run = runs[i]
        log(f"=== [{i + 1}/{len(runs)}] {run['name']}")
        out = out_root / safe_name(run["name"])
        try:
            done(i, _load_project(run, log, transcripts) if run.get("kind") == "project" else _convert_one(run, out, spec, log, transcripts, tx_table, threads, ifimg, reg))
        except Exception as e:
            if not isinstance(e, importer.ImportError_):
                traceback.print_exc()
            fail(i, e)
        gc.collect()

    if w <= 1:
        for i in range(len(runs)):
            serial(i)
    else:
        for i, run in enumerate(runs):
            if run.get("kind") == "project":
                serial(i)
        import multiprocessing as mp, threading
        from concurrent.futures import ProcessPoolExecutor, as_completed
        ctx = mp.get_context("fork")                    # the workers inherit the loaded modules; the log queue is handed over the same way
        q = ctx.Queue()

        def drain():
            while True:
                m = q.get()
                if m is None:
                    return
                log(f"[{m[0]}] {m[1]}")
        th = threading.Thread(target=drain, daemon=True)
        th.start()
        try:
            with ProcessPoolExecutor(max_workers=w, mp_context=ctx, initializer=_init_worker, initargs=(q,)) as ex:
                futs = {}
                for i, run in enumerate(runs):
                    if run.get("kind") != "project":
                        log(f"=== [{i + 1}/{len(runs)}] {run['name']} queued")
                        futs[ex.submit(_worker, i, run, out_root / safe_name(run["name"]), spec, transcripts, tx_table, threads, ifimg, reg)] = i
                for f in as_completed(futs):
                    i = futs[f]
                    try:
                        _, p, err = f.result()
                    except Exception as e:                      # the worker process died (for example out of memory)
                        p, err = None, f"worker process failed: {e}"
                    done(i, p) if p is not None else fail(i, err)
        finally:
            q.put(None)
            th.join(5)
    pieces = [got[i] for i in sorted(got)]
    first = next((Path(runs[i]["path"]) if runs[i].get("kind") == "project" else out_root / safe_name(runs[i]["name"]) for i in sorted(got)), None)
    failed = list(failed.values())
    merged = None
    if do_merge and pieces:
        merged = merge(pieces, out_root / MERGED, log, ann)
        if he:
            attach_he(out_root, merged, he, [p["name"] for p in pieces], log, cores=cpu_count() if not workers else min(int(workers), cpu_count()), threads=he_threads)
    elif he and pieces:
        log("H&E is only attached to the merged project: tick the merge option to use it")
    if ann and ann.get("path") and not (do_merge and pieces):
        log("the annotation file is only used for the merged project: tick the merge option to use it")
    if he_reg and merged:                                       # registered H&E OME-TIFFs for the merged samples (<run>_<tissue>)
        try:
            from . import he_registered
            he_registered.run(he_reg["folder"], he_reg["index"], merged, log)
        except (importer.ImportError_, OSError, ValueError) as e:
            log(f"registered H&E ignored: {e}")
    elif he_reg:
        log("registered H&E is only attached to the merged project: tick the merge option to use it")
    log(f"batch finished: {len(pieces)} run(s) converted, {len(failed)} failed" + (": " + ", ".join(f["name"] for f in failed) if failed else ""))
    return {"done": [p["name"] for p in pieces], "failed": failed, "merged": merged, "first": first}
