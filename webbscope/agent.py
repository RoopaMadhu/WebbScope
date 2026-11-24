"""Claude's access to the open project (the MCP connector, webbscope/mcp.py), under the user's control.

The connector holds the agent key (access.agent_key(): derived from the access token, so it changes with it), which reaches only /api/agent/...;
the access token itself (the user's browser) reaches everything else, including the sharing level of each project, which only the user can set.
Levels, each including the ones before it:
  off         nothing but this status (the default for every project)
  structure   names and sizes: genes, columns, categories and their values, tissues, donors, groups; drive the viewer (show a tissue, a colouring)
  summaries   numbers aggregated over at least MIN_CELLS cells: cell counts and composition, mean expression and the share of expressing cells
  full        cell-level rows (metadata, coordinates, expression) and pictures of the viewer
  code        run code with the user's own Python or R on the project's data, here on the machine running WebbScope
Every call is written to the activity log (~/.webbscope/agent_log/<project>.jsonl), shown with the level in the viewer (About > Claude access).
Settings live in ~/.webbscope/sharing.json (the home folder, not the cache, which may be a shared disk).
"""
import base64, json, os, shlex, shutil, signal, subprocess, sys, threading, time
from pathlib import Path
import numpy as np
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import ORJSONResponse
from pydantic import BaseModel

LEVELS = ["off", "structure", "summaries", "full", "code"]
LEVEL_TEXT = {
    "off": "Off: Claude sees nothing of this project.",
    "structure": "Structure: names and sizes (genes, columns and their values, tissues, donors, groups), no measurements. Claude can also show a tissue or a colouring in the viewer.",
    "summaries": "Summaries: also numbers over groups of at least 10 cells (cell counts, composition, mean expression, share of expressing cells).",
    "full": "Full: also single cells (metadata, coordinates, expression) and pictures of the viewer.",
    "code": "Run code: also runs Python or R code that Claude writes, with your account's rights, on this machine, with the project's data.",
}
MIN_CELLS = 10                  # a summary row of fewer cells shows its count only
MAX_GENES = 200                 # genes per summary / cell request
MAX_ROWS = 50_000               # cells per request
router = APIRouter()
_lock = threading.Lock()


def _home():
    return Path.home() / ".webbscope"


def _srv():
    from . import server
    return server


def project_key():
    """The open project, as the cache folder's name (one per project, and per .h5ad of a folder)."""
    s = _srv()
    return Path(s.S.cache_dir).name if s.S is not None else None


def _settings():
    f = _home() / "sharing.json"
    try:
        return json.loads(f.read_text())
    except (OSError, ValueError):
        return {}


def _save_settings(d):
    f = _home() / "sharing.json"
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_name(f".sharing.{os.getpid()}.tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as fh:
        json.dump(d, fh, indent=1)
    os.replace(tmp, f)


def settings_of(key=None):
    key = key or project_key()
    p = _settings().get("projects", {}).get(key or "", {})
    return {"level": p.get("level", "off") if p.get("level") in LEVELS else "off", "python": p.get("python") or _settings().get("python") or "",
            "r": p.get("r") or _settings().get("r") or ""}


def level():
    return settings_of()["level"]


def _log_file(key=None):
    return _home() / "agent_log" / f"{key or project_key()}.jsonl"


def log(tool, params, outcome, **extra):
    """One line of the activity log of the open project."""
    key = project_key()
    if key is None:
        return
    f = _log_file(key)
    rec = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "tool": tool, "level": level(), "outcome": outcome,
           "params": {k: (v if len(json.dumps(v, default=str)) < 300 else str(v)[:300] + "…") for k, v in (params or {}).items()}, **extra}
    with _lock:
        f.parent.mkdir(parents=True, exist_ok=True)
        new = not f.exists()
        with open(f, "a") as fh:
            fh.write(json.dumps(rec, default=str) + "\n")
        if new:
            os.chmod(f, 0o600)


def need(tool, minimum, params=None):
    """Refuse (and log) a call the project's sharing level does not allow."""
    s = _srv()
    if s.S is None:
        raise HTTPException(409, "no project is open in WebbScope")
    lv = level()
    if LEVELS.index(lv) < LEVELS.index(minimum):
        log(tool, params, "refused")
        raise HTTPException(403, f"WebbScope's sharing level for this project is '{lv}'; {tool} needs '{minimum}' or above. "
                                 "Only the user can change it, in the viewer: About > Claude access.")


# ---------------------------------------------------------------- user side (access token only; the agent key does not reach these)
class SharingReq(BaseModel):
    level: str
    python: str = ""
    r: str = ""


def connector_command():
    """How to start the MCP connector of this installation (it runs where Claude runs; on this machine it finds the server by itself)."""
    if getattr(sys, "frozen", False):
        exe = Path(sys.executable)
        return f"{shlex.quote(str(exe))} --mcp"
    return f"{shlex.quote(sys.executable)} {shlex.quote(str(Path(__file__).resolve().parent.parent / 'run_app.py'))} --mcp"


@router.get("/api/sharing", response_class=ORJSONResponse)
def sharing_get():
    from . import access
    key = project_key()
    if key is None:
        raise HTTPException(409, "no project is open")
    f, rows = _log_file(key), []
    if f.exists():
        rows = [json.loads(l) for l in f.read_text().splitlines()[-200:] if l.strip()]
    return {**settings_of(key), "levels": [{"key": k, "text": LEVEL_TEXT[k]} for k in LEVELS], "log": rows[::-1],
            "agent_key": access.agent_key(), "project": str(_srv().PROJ), "connector": connector_command()}


@router.post("/api/sharing", response_class=ORJSONResponse)
def sharing_set(r: SharingReq):
    key = project_key()
    if key is None:
        raise HTTPException(409, "no project is open")
    if r.level not in LEVELS:
        raise HTTPException(400, f"level must be one of {', '.join(LEVELS)}")
    d = _settings()
    old = settings_of(key)["level"]
    d.setdefault("projects", {})[key] = {"level": r.level, "python": r.python.strip(), "r": r.r.strip(), "path": str(_srv().PROJ), "changed": time.strftime("%Y-%m-%d %H:%M:%S")}
    _save_settings(d)
    if old != r.level:
        log("sharing level", {"from": old, "to": r.level}, "set by the user")
    return {"level": r.level}


@router.delete("/api/sharing/log", response_class=ORJSONResponse)
def sharing_clear_log():
    f = _log_file()
    if f.exists():
        f.unlink()
    return {"ok": True}


# ---------------------------------------------------------------- agent side
@router.get("/api/agent/status", response_class=ORJSONResponse)
def agent_status():
    """Always allowed: which project is open and what the user lets Claude do with it."""
    s = _srv()
    if s.S is None:
        return {"project": None, "level": "off", "note": "no project is open in WebbScope"}
    lv = level()
    log("status", {}, "ok")
    return {"project": s.PROJ.name, "level": lv, "allowed": LEVEL_TEXT[lv], "levels": LEVEL_TEXT,
            "note": "The level is set by the user in the viewer (About > Claude access); ask them if you need more."}


def _grouping(S):
    """Per-cell keys the summaries can group or filter by: {name: (codes int64, levels)}."""
    n = S.n
    samp = np.empty(n, np.int64)
    for i, smp in enumerate(S.samples):
        a, b = S.slices[smp]
        samp[a:b] = i
    donors = sorted(set(S.donor.values()))
    di = {d: i for i, d in enumerate(donors)}
    sdon = np.array([di[S.donor[x]] for x in S.samples], np.int64)
    out = {"sample": (samp, list(S.samples)), "donor": (sdon[samp], donors)}
    if S.groups.column:
        gl = list(S.groups.order) or sorted({S.groups.of(x) for x in S.samples})
        gi = {g: i for i, g in enumerate(gl)}
        out["group"] = (np.array([gi.get(S.groups.of(x), -1) for x in S.samples], np.int64)[samp], gl)
    for k, v in S.cats.items():
        if v.get("multi"):
            continue
        c = v["codes"].astype(np.int64)
        out[k] = (np.where(c == 65535, -1, c), list(v["levels"]))
    return out


def _fields(S):
    return {k: (S.cats[k]["label"] if k in S.cats else {"sample": "tissue (sample)", "donor": "donor", "group": f"group ({S.groups.column})"}[k]) for k in _grouping(S)}


@router.get("/api/agent/describe", response_class=ORJSONResponse)
def agent_describe():
    """Structure: what the project holds. Summaries and above add the number of cells per tissue and per value."""
    need("describe", "structure")
    s = _srv(); S = s.S; lv = level(); summ = LEVELS.index(lv) >= LEVELS.index("summaries")
    G = _grouping(S)
    fields = {}
    for k, label in _fields(S).items():
        codes, levels = G[k]
        f = {"label": label, "n_values": len(levels), "values": levels[:300]}
        if summ:
            cnt = np.bincount(codes[codes >= 0], minlength=len(levels))
            f["cells_per_value"] = {str(l): int(c) for l, c in zip(levels[:300], cnt[:300])}
        fields[k] = f
    m = S.map
    out = {"project": s.PROJ.name, "n_cells": int(S.n), "n_genes": len(S.genes), "genes": list(S.genes),
           "values": {"counts": "raw counts; the viewer's expression = log1p(counts / total counts of the cell * scale)", "lognorm": "log-normalised values as given",
                      "both": "raw counts and the given log-normalised values"}[S.values],
           "norm_scale": S.norm_scale, "layers": ["logcpx"] + (["counts"] if S.has_counts else []),
           "roles": {k: v for k, v in m.items() if v}, "fields": fields,
           "groups": {"column": S.groups.column, "reference": S.groups.reference, "order": S.groups.order, "ordered": S.groups.ordered,
                      "note": "donors are the unit of the viewer's statistics; the group of a donor is set by the user"},
           "numeric_columns": sorted(v for v in S.numeric.values()), "embeddings": sorted(S.umap) + (["pca"] if S.pca is not None else []),
           "per_run_embedding": bool(getattr(s, "EMB_PER_RUN", False)), "cell_outlines": S.poly_off is not None,
           "coordinates": "micrometres, per tissue section", "level": lv}
    log("describe", {}, "ok")
    return out


class SummaryReq(BaseModel):
    kind: str = "counts"            # counts | expression
    by: list = []                   # fields to group by: sample, donor, group, cell_type, cell_state, meta:<column> ...
    genes: list = []                # expression: genes (at most MAX_GENES)
    layer: str = "logcpx"
    filters: dict = {}              # {field: [values]}: only these cells


def _mask(S, G, filters):
    keep = np.ones(S.n, bool)
    for k, vals in (filters or {}).items():
        if k not in G:
            raise HTTPException(400, f"unknown field {k!r}; fields: {', '.join(G)}")
        codes, levels = G[k]
        want = {i for i, l in enumerate(levels) if str(l) in {str(v) for v in vals}}
        if not want:
            raise HTTPException(400, f"none of {list(vals)[:5]} is a value of {k}")
        keep &= np.isin(codes, list(want))
    return keep


@router.post("/api/agent/summary", response_class=ORJSONResponse)
def agent_summary(r: SummaryReq):
    """Cell counts, or mean expression and the share of cells with the gene (> 0), per combination of the `by` fields. Rows of fewer than MIN_CELLS
    cells keep their count only."""
    params = r.model_dump()
    need("summary", "summaries", params)
    S = _srv().S
    G = _grouping(S)
    for b in r.by:
        if b not in G:
            raise HTTPException(400, f"unknown field {b!r}; fields: {', '.join(G)}")
    if len(r.by) > 4:
        raise HTTPException(400, "group by at most 4 fields")
    keep = _mask(S, G, r.filters)
    for b in r.by:
        keep &= G[b][0] >= 0                                    # cells without a value of that field
    idx = np.flatnonzero(keep)
    if r.by:
        codes = np.stack([G[b][0][idx] for b in r.by])
        dims = [len(G[b][1]) for b in r.by]
        key = np.ravel_multi_index(codes, dims)
        uk, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
        combos = np.array(np.unravel_index(uk, dims)).T
    else:
        inv, cnt, combos = np.zeros(len(idx), np.int64), np.array([len(idx)]), np.zeros((1, 0), np.int64)
    head = list(r.by) + ["n_cells"]
    rows = [[G[b][1][c] for b, c in zip(r.by, combo)] + [int(n)] for combo, n in zip(combos, cnt)]
    hidden = int((cnt < MIN_CELLS).sum())
    if r.kind == "expression":
        found = S.resolve_genes([str(g) for g in r.genes])
        genes = found["found"] + found["excluded"]
        if not genes:
            raise HTTPException(400, f"none of these genes is on the panel: {r.genes[:8]}")
        if len(genes) > MAX_GENES:
            raise HTTPException(400, f"at most {MAX_GENES} genes per request")
        if r.layer not in ("logcpx", "counts"):
            raise HTTPException(400, "layer must be logcpx or counts")
        for g in genes:
            v = S.gene_vector(g, 0, S.n, r.layer)[idx]
            mean = np.bincount(inv, weights=v, minlength=len(cnt)) / cnt
            pct = np.bincount(inv, weights=(v > 0).astype(np.float64), minlength=len(cnt)) / cnt
            head += [f"{g} mean", f"{g} % expressing"]
            for row, m_, p_, n_ in zip(rows, mean, pct, cnt):
                row += [round(float(m_), 5), round(100 * float(p_), 2)] if n_ >= MIN_CELLS else [None, None]
        missing = found["missing"]
    else:
        missing = []
    out = {"columns": head, "rows": rows, "cells_used": int(len(idx)),
           "note": (f"{hidden} row(s) of fewer than {MIN_CELLS} cells show their count only. " if hidden else "")
           + ("Expression: " + ("log1p(counts / total counts * scale)" if r.layer == "logcpx" and S.values == "counts" else r.layer) + "." if r.kind == "expression" else ""),
           **({"genes_not_on_panel": missing} if missing else {})}
    log("summary", params, "ok", rows=len(rows))
    return out


class CellsReq(BaseModel):
    columns: list = []              # fields or numeric columns of the cell table (sample, donor, group, cell_type, meta:<column>, a numeric column ...)
    genes: list = []
    layer: str = "logcpx"
    filters: dict = {}
    limit: int = 10_000
    offset: int = 0


@router.post("/api/agent/cells", response_class=ORJSONResponse)
def agent_cells(r: CellsReq):
    """Single cells: id, tissue, x, y, the chosen fields and genes, for the cells passing the filters (in the project's order, `limit` from `offset`)."""
    params = r.model_dump()
    need("cells", "full", params)
    S = _srv().S
    G = _grouping(S)
    if not 1 <= r.limit <= MAX_ROWS or r.offset < 0:
        raise HTTPException(400, f"limit must be 1 to {MAX_ROWS}")
    idx = np.flatnonzero(_mask(S, G, r.filters))
    total = len(idx)
    idx = idx[r.offset:r.offset + r.limit]
    m = S.map
    cols = {"cell_id": S.df[m["cell_id"]].astype(str).values[idx].tolist(), "sample": [S.samples[i] for i in G["sample"][0][idx]],
            "x": np.round(S.x[idx], 2).tolist(), "y": np.round(S.y[idx], 2).tolist()}
    for c in r.columns:
        if c in G:
            codes, levels = G[c]
            cols[c] = [levels[k] if k >= 0 else None for k in codes[idx]]
        elif c in S.df.columns and S.df[c].dtype.kind in "fiub":
            cols[c] = S.df[c].values[idx].tolist()
        else:
            raise HTTPException(400, f"unknown column {c!r}")
    found = S.resolve_genes([str(g) for g in r.genes])
    genes = found["found"] + found["excluded"]
    if len(genes) > MAX_GENES:
        raise HTTPException(400, f"at most {MAX_GENES} genes per request")
    for g in genes:
        cols[g] = np.round(S.gene_vector(g, 0, S.n, r.layer)[idx], 5).tolist()
    log("cells", params, "ok", rows=len(idx))
    return {"total_matching": int(total), "offset": r.offset, "n": int(len(idx)), "columns": cols,
            **({"genes_not_on_panel": found["missing"]} if found["missing"] else {})}


# ---------------------------------------------------------------- run code (the user's own Python or R, here)
def export_dir():
    """The project's data for code, made once per project in its cache: cells.parquet (the cell table plus tissue, donor, group), genes.txt and the
    expression matrix as links to the cache's arrays (genes x cells CSR; no copy), with a README saying how to load them."""
    S = _srv().S
    d = Path(S.cache_dir) / "agent_export"
    stamp = d / ".made"
    src = Path(S.cache_dir) / "cache_info.json"
    if stamp.exists() and stamp.stat().st_mtime >= src.stat().st_mtime and (d / "groups.json").exists() \
            and json.loads((d / "groups.json").read_text()) == S.groups.describe():
        return d
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    G = _grouping(S)
    t = S.df.copy()
    t.insert(0, "webbscope_sample", [S.samples[i] for i in G["sample"][0]])
    t.insert(1, "webbscope_donor", [G["donor"][1][i] for i in G["donor"][0]])
    if "group" in G:
        t.insert(2, "webbscope_group", [G["group"][1][i] if i >= 0 else None for i in G["group"][0]])
    t.to_parquet(d / "cells.parquet", index=False)
    (d / "genes.txt").write_text("\n".join(S.genes) + "\n")
    for f in ("csr_indptr.npy", "csr_indices.npy", "csr_data.npy"):
        os.symlink(Path(S.cache_dir) / f, d / f)
    (d / "groups.json").write_text(json.dumps(S.groups.describe()))
    stored = {"counts": "raw counts", "lognorm": "log-normalised values (as given)", "both": "raw counts"}[S.values]
    (d / "README.txt").write_text(f"""WebbScope project data for code (read only; do not change these files).
cells.parquet  one row per cell, in the order of the matrix columns: the project's cell table plus webbscope_sample (tissue), webbscope_donor
               and webbscope_group (the group the user chose in the viewer: {S.groups.column or 'none'}; reference {S.groups.reference or 'none'})
genes.txt      the genes, in the order of the matrix rows
csr_indptr.npy, csr_indices.npy, csr_data.npy   the expression matrix, genes x cells, CSR, values = {stored}
Python:
    import numpy as np, pandas as pd, scipy.sparse as sp
    X = sp.csr_matrix((np.load("csr_data.npy"), np.load("csr_indices.npy"), np.load("csr_indptr.npy")), shape=({len(S.genes)}, {S.n}))   # genes x cells
    cells = pd.read_parquet("cells.parquet"); genes = open("genes.txt").read().split()
""" + (f"""    The viewer's expression (logcpx) = log1p(counts / total counts of the cell * {S.norm_scale}); total counts = column {S.map.get('lib_size')} of cells.parquet.
""" if S.values == "counts" else "") + """R: arrow::read_parquet("cells.parquet"); the matrix through reticulate (scipy) or RcppCNPy.
""")
    stamp.write_text(time.strftime("%Y-%m-%d %H:%M:%S"))
    return d


class RunReq(BaseModel):
    language: str = "python"        # python | r
    code: str
    timeout: int = 300              # seconds, at most 1800


def _interpreter(lang):
    st = settings_of()
    if lang == "python":
        exe = st["python"] or shutil.which("python3") or shutil.which("python")
        if not exe:
            raise HTTPException(400, "no Python found on this machine; the user can set one in About > Claude access")
        return shlex.split(exe), "script.py"
    exe = st["r"] or shutil.which("Rscript")
    if not exe:
        raise HTTPException(400, "no Rscript found on this machine; the user can set one in About > Claude access")
    return shlex.split(exe), "script.R"


@router.post("/api/agent/run", response_class=ORJSONResponse)
def agent_run(r: RunReq):
    """Run the code with the user's interpreter in a new folder of the project's cache, the data folder in WEBBSCOPE_DATA. Returns the output,
    the files made and up to 4 PNG pictures (base64)."""
    params = {"language": r.language, "timeout": r.timeout, "code": r.code}
    need("run code", "code", params)
    if r.language not in ("python", "r"):
        raise HTTPException(400, "language must be python or r")
    if not 1 <= r.timeout <= 1800:
        raise HTTPException(400, "timeout must be 1 to 1800 seconds")
    cmd, name = _interpreter(r.language)
    data = export_dir()
    S = _srv().S
    with _lock:
        VIEW["seq"] += 1
        run = Path(S.cache_dir) / "agent_runs" / (time.strftime("%Y%m%d-%H%M%S-") + str(VIEW["seq"]))
    run.mkdir(parents=True)
    (run / name).write_text(r.code)
    env = {**os.environ, "WEBBSCOPE_DATA": str(data), "MPLBACKEND": "Agg"}
    env.pop("WEBBSCOPE_TOKEN", None)
    t0 = time.time()
    p = subprocess.Popen(cmd + [name], cwd=run, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        out, err = p.communicate(timeout=r.timeout)
        code, note = p.returncode, ""
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, err = p.communicate()
        code, note = None, f"stopped after {r.timeout} s"
    files = sorted(str(f.relative_to(run)) for f in run.rglob("*") if f.is_file() and f.name != name)
    pics = []
    for f in files:
        q = run / f
        if q.suffix.lower() == ".png" and q.stat().st_size <= 3_000_000 and len(pics) < 4:
            pics.append({"file": f, "png_base64": base64.b64encode(q.read_bytes()).decode()})
    cut = lambda s: s if len(s) <= 60_000 else s[:20_000] + f"\n… {len(s) - 50_000:,} characters left out …\n" + s[-30_000:]
    log("run code", params, "ok" if code == 0 else f"exit {code}" if code is not None else note, seconds=round(time.time() - t0, 1), folder=str(run))
    return {"exit_code": code, "note": note, "seconds": round(time.time() - t0, 1), "stdout": cut(out), "stderr": cut(err), "folder": str(run),
            "files": files, "pictures": pics, "data": str(data)}


# ---------------------------------------------------------------- the viewer: what it shows, and showing something else
VIEW = {"state": None, "commands": [], "snap": None, "snap_event": threading.Event(), "seq": 0}


class ViewCmd(BaseModel):
    sample: str = ""                # a tissue to show
    color_by: str = ""              # a field (cell_type, cell_state, meta:<column> ...) or a gene
    filters: dict = {}              # {field: [values]} shown, the rest greyed out; {} leaves the filters as they are


@router.get("/api/agent/view", response_class=ORJSONResponse)
def agent_view_get():
    need("view", "structure")
    log("view", {}, "ok")
    return {"viewer": VIEW["state"] or "no viewer page has reported yet (the user may have no WebbScope page open)"}


@router.post("/api/agent/view", response_class=ORJSONResponse)
def agent_view_set(r: ViewCmd):
    params = r.model_dump()
    need("show in viewer", "structure", params)
    S = _srv().S
    if r.sample and r.sample not in S.slices:
        raise HTTPException(400, f"no tissue {r.sample!r}")
    if r.color_by and r.color_by not in S.cats and r.color_by not in S.gene_index:
        raise HTTPException(400, f"{r.color_by!r} is neither a field nor a gene of this project")
    for k, vals in (r.filters or {}).items():                   # the viewer filters by its own fields only (sample, donor and group are chosen with the tissue)
        if k not in S.cats or S.cats[k].get("multi"):
            raise HTTPException(400, f"the viewer cannot filter by {k!r}; fields it can filter by: {', '.join(k2 for k2, v in S.cats.items() if not v.get('multi'))}")
        bad = [v for v in vals if str(v) not in set(map(str, S.cats[k]["levels"]))]
        if bad:
            raise HTTPException(400, f"not values of {k}: {bad[:5]}")
    with _lock:
        VIEW["seq"] += 1
        VIEW["commands"].append({"seq": VIEW["seq"], **params})
        VIEW["commands"] = VIEW["commands"][-20:]
    log("show in viewer", params, "sent")
    return {"sent": True, "note": "the open viewer page applies it within a few seconds"}


@router.post("/api/agent/snapshot", response_class=ORJSONResponse)
def agent_snapshot():
    """A picture of the viewer as the user sees it (PNG, base64)."""
    need("snapshot", "full")
    with _lock:
        VIEW["seq"] += 1
        VIEW["commands"].append({"seq": VIEW["seq"], "snapshot": True})
        VIEW["snap"] = None
        VIEW["snap_event"].clear()
    if not VIEW["snap_event"].wait(25) or not VIEW["snap"]:
        log("snapshot", {}, "no viewer answered")
        raise HTTPException(504, "no open viewer page answered (the user needs a WebbScope page open on this project)")
    log("snapshot", {}, "ok")
    return {"png_base64": VIEW["snap"]}


# the viewer page's side (access token): report its state, fetch the commands, deliver a snapshot
@router.post("/api/viewer/state", response_class=ORJSONResponse)
async def viewer_state(request: Request):
    d = await request.json()
    VIEW["state"] = {k: d.get(k) for k in ("sample", "color_by", "gene", "filters", "selected", "windows", "zoom")}
    return {"ok": True}


@router.get("/api/viewer/commands", response_class=ORJSONResponse)
def viewer_commands(after: int = 0):
    return {"level": level() if _srv().S is not None else "off", "commands": [c for c in VIEW["commands"] if c["seq"] > after], "seq": VIEW["seq"]}


@router.post("/api/viewer/snapshot", response_class=ORJSONResponse)
async def viewer_snapshot(request: Request):
    d = await request.json()
    png = str(d.get("png_base64") or "")
    if len(png) > 12_000_000:
        raise HTTPException(413, "picture too large")
    VIEW["snap"] = png
    VIEW["snap_event"].set()
    return {"ok": True}
