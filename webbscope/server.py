"""FastAPI app: the browser's analysis endpoints (unchanged) + opening a project folder.

Derived from xbrowser/server.py. Differences: engine modules live in webbscope.engine; the data come from an imported project
(importer.py) opened at run time instead of a fixed --cache; H&E folders are read from the project folder; / shows an Open-project
page until a project is loaded.
"""
import json, os, re, threading, traceback
from pathlib import Path
import numpy as np
from fastapi import FastAPI, HTTPException, Response
from fastapi.responses import FileResponse, ORJSONResponse
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.staticfiles import StaticFiles
from .engine.store import Store
from .engine.smooth import Smoother, SCOPES, MAX_SIGMA
from .engine.signature_defs import EXCLUDED_GENES
from .engine.neighbors import Neighbors
from .engine import sigstats, violin
from .engine.glomhet import GlomHet
from . import importer, xenium, batch, qc
from .profiler import PROF
import time
from .transcripts import Transcripts
from . import he_slide, he_registered, if_image, imc
from pydantic import BaseModel

WEB = Path(__file__).resolve().parent / "web"
app = FastAPI(title="WebbScope")
app.add_middleware(GZipMiddleware, minimum_size=2000)
S: Store = None
HE_EXTRA = None         # folder holding registered H&E tiles (he/<key>/…, he/index.json) added when a project was opened, in the cache; else None
IMX = None              # IMC slide of the open project: {store: Path, index: dict (imc.read_store)}, else None
IFX = None              # IF images of the open project: {images: {key: entry}, samples: {sample: key}} (if_image.py), else None
SM: Smoother = None
NB: Neighbors = None
GH: GlomHet = None
PROJ: Path = None          # the open project folder (H&E pyramids are served from it)
_gene_cache: dict = {}
TX: Transcripts = None    # transcript overlay data of the open project, if it has any
DEFAULT_GREY = False   # project.json "default_grey": start with the Grey cells display on (merged batch projects: the cell types are per tissue)
EMB_LABEL = None       # set when the project's embeddings are not Harmony-corrected (Xenium import)


def raw(arr):
    return Response(np.ascontiguousarray(arr).tobytes(), media_type="application/octet-stream")


class ExprReq(BaseModel):
    """Gene / signature filter: keep cells whose score (mean over `genes`, smoothed / z-scored like the viewer's colouring) is > thr ("gt") or <= thr ("le")."""
    genes: list[str]
    op: str = "gt"
    thr: float = 0.0
    layer: str = "logcpx"
    sigma: float = 0.0
    scope: str = "none"
    within: str = "tissue"


def _expr(specs):
    """Validate gene / signature filters and return them as plain dicts (hashable once repr'd)."""
    out = []
    if len(specs) > 6:
        raise HTTPException(400, "at most 6 gene / signature filters")
    for e in specs:
        genes = [g for g in e.genes if g]
        if not genes or len(genes) > 500:
            raise HTTPException(400, "each gene filter needs 1-500 genes")
        bad = [g for g in genes if g not in S.gene_index]
        if bad:
            raise HTTPException(404, f"genes not in panel: {bad[:5]}")
        if e.op not in ("gt", "le") or e.layer not in ("logcpx", "counts") or e.scope not in SCOPES or not 0 <= e.sigma <= MAX_SIGMA:
            raise HTTPException(400, "bad op, layer, scope or sigma")
        try:
            SM.check_within(e.within)
        except ValueError as ex:
            raise HTTPException(400, str(ex))
        out.append({"genes": genes, "op": e.op, "thr": float(e.thr), "layer": e.layer, "sigma": float(e.sigma), "scope": e.scope, "within": e.within})
    return out


def sl(sample):
    if sample not in S.slices:
        raise HTTPException(404, f"unknown sample {sample!r}")
    return S.slices[sample]


@app.middleware("http")
async def no_stale_ui(request, call_next):
    """Make browsers revalidate the UI files so a restart always serves the latest JS/CSS/HTML."""
    p = request.url.path
    if S is None and p.startswith("/api/") and p not in ("/api/status", "/api/open", "/api/pick", "/api/columns") and not p.startswith("/api/xenium/") and not p.startswith("/api/profile"):
        return Response("no project is open", status_code=409)
    t0 = time.perf_counter()
    r = await call_next(request)
    if p.startswith("/api/") and not p.startswith("/api/profile") and p != "/api/status":       # not the profile page's own polling
        route = request.scope.get("route")
        PROF.request(route.path if route is not None else p, time.perf_counter() - t0, int(r.headers.get("content-length") or 0))
    if p in ("/", "/welcome") or p.startswith("/static/"):
        r.headers["Cache-Control"] = "no-cache"
    return r


@app.get("/api/meta", response_class=ORJSONResponse)
def meta():
    d = S.describe()
    if EMB_LABEL:
        d["embedding_label"] = EMB_LABEL
    d["default_grey"] = bool(DEFAULT_GREY)
    d["has_gloms"] = (Path(S.cache_dir) / "gloms" / "gloms.json").exists()
    d["has_transcripts"] = TX is not None
    d["has_he_slide"] = (PROJ / "he_slide.json").exists()
    d["has_if"] = bool(IFX and IFX["samples"])
    d["has_imc"] = bool(IMX)
    if TX is not None:                                           # transcripts in the whole run (= slide) each tissue came from: the rows of its transcripts.parquet
        for smp in d["samples"]:
            try:
                q = TX.qc(smp["id"])
                smp["tx_slide_total"] = int(q["run"]["rows"]) if q and q.get("run") and q["run"].get("rows") else None
                smp["tx_tissue_total"] = int(q["tissue_rows"]) if q and q.get("tissue_rows") is not None else None   # rows of the transcripts file that belong to this tissue
            except Exception:
                smp["tx_slide_total"] = smp["tx_tissue_total"] = None
    return d


@app.get("/api/validation", response_class=ORJSONResponse)
def validation():
    return S.validation()


@app.get("/api/sample/{sample}/xy")
def xy(sample: str):
    a, b = sl(sample)
    return raw(np.stack([S.x[a:b], S.y[a:b]], 1).astype(np.float32))


@app.get("/api/paths/{sample}", response_class=ORJSONResponse)
def paths(sample: str):
    """Pathologist region polygons of one sample as a GeoJSON FeatureCollection (same coordinates as the cells)."""
    sl(sample)
    p = S.paths.get(sample)
    if p is None:
        raise HTTPException(404, "no pathologist annotation matched this sample")
    return {"type": "FeatureCollection", "features": p["features"], "tissues": p["tissues"]}


# Coarse labels the cell-level QC removed from the Xenium cells (present in scFFPE); dropped from the scFFPE side by default for a like-for-like comparison.
SC_EXCLUDED = ("Low quality cells", "Proliferating cells")


@app.get("/api/compare/{sample}", response_class=ORJSONResponse)
def compare(sample: str, level: str = "cell_type", include_excluded: bool = False, scffpe: str = ""):
    """Cell fractions of one Xenium sample vs its scFFPE sample(s). `scffpe` (comma-separated scFFPE sample ids) overrides the metadata link."""
    a, b = sl(sample)
    if level not in ("cell_type", "cell_state") or level not in S.cats:
        raise HTTPException(400, "level must be cell_type or cell_state")
    lv, codes = S.cats[level]["levels"], S.cats[level]["codes"][a:b]
    xc = np.bincount(codes[codes != 65535], minlength=len(lv))
    out = {"sample": sample, "level": level, "xenium": {"n": int(b - a), "counts": {lv[i]: int(c) for i, c in enumerate(xc) if c}},
           "excluded_labels": list(SC_EXCLUDED), "include_excluded": include_excluded, "scffpe": None, "status": "no_scffpe_data", "scffpe_id": None,
           "shared_with": [], "candidates": []}
    SCF = S.scffpe
    if not SCF:
        return out
    link = SCF["xenium"].get(sample, {"status": "no_scffpe_id", "scffpe": None, "shared_with": []})
    ids = [x for x in scffpe.split(",") if x]
    if ids and any(i not in SCF["scffpe"] for i in ids):
        raise HTTPException(404, "unknown scFFPE sample")
    out.update(status=("override" if ids else link["status"]), scffpe_id=link["scffpe"], shared_with=link.get("shared_with", []))
    if link["status"] != "ok":
        blk = sample.split("-")[0]
        out["candidates"] = sorted(x for x in SCF["scffpe"] if x.split("-")[0] == blk)     # same block id, NOT linked in the metadata
    ids = ids or ([link["scffpe"]] if link["status"] == "ok" else [])
    if ids:
        counts, excluded, n_all = {}, {}, 0
        for i in ids:
            for coarse, fine, n in SCF["scffpe"][i]["rows"]:
                n_all += n
                if coarse in SC_EXCLUDED and not include_excluded:
                    excluded[coarse] = excluded.get(coarse, 0) + n
                    continue
                k = coarse if level == "cell_type" else fine
                counts[k] = counts.get(k, 0) + n
        out["scffpe"] = {"ids": ids, "n": int(sum(counts.values())), "n_all": int(n_all), "counts": counts, "excluded": excluded}
    return out


_sig_cache: dict = {}


def _genes_of(sig_id, genes):
    if sig_id:
        if sig_id not in S.sig_by_id:
            raise HTTPException(404, f"unknown signature {sig_id!r}")
        return list(S.sig_by_id[sig_id]["genes"])
    bad = [g for g in genes or [] if g not in S.gene_index]
    if not genes or bad:
        raise HTTPException(400, f"give a signature id or panel genes (not on the panel: {bad[:5]})")
    return [g for g in genes if g not in EXCLUDED_GENES]


@app.get("/api/signatures", response_class=ORJSONResponse)
def signatures():
    """Curated signatures (genes present on the panel, genes missing, coherence checks) + subtype order used by the analyses."""
    return {**S.signatures, "subtype_order": sigstats.SUBTYPE_ORDER}


class SigSummaryReq(BaseModel):
    sig: str | None = None
    genes: list[str] | None = None
    sigma: float = 25.0
    scope: str = "all"                      # z-scoring: none | sample | subtype | type | all (pooled over the samples of that group)
    group_by: str = "subtype"               # subtype (5 groups) | type (Control / No LN / LN)
    filters: dict[str, list[int]] = {}
    hidden: dict[str, list[int]] = {}
    expr: list[ExprReq] = []


@app.post("/api/signature/summary", response_class=ORJSONResponse)
def signature_summary(req: SigSummaryReq):
    """Per-sample and per-disease-subtype statistics of a (smoothed) signature over the cells passing `filters`, with across-subtype tests."""
    genes = _genes_of(req.sig, req.genes)
    if req.scope not in SCOPES or req.group_by not in ("subtype", "type") or not 0 <= req.sigma <= MAX_SIGMA:
        raise HTTPException(400, "bad scope, group_by or sigma")
    ex = _expr(req.expr)
    key = ("sum", tuple(genes), req.sigma, req.scope, req.group_by, repr(sorted((k, sorted(v)) for k, v in req.filters.items())), repr(sorted((k, sorted(v)) for k, v in req.hidden.items())), repr(ex))
    if key not in _sig_cache:
        _sig_cache[key] = sigstats.summary(S, SM, NB, genes, req.sigma, req.scope, req.filters, req.hidden, group_by=req.group_by, expr=ex)
        _sig_cache[key]["sig"] = req.sig
    return _sig_cache[key]


class SigLocReq(BaseModel):
    a: str = "ifn1"
    b: str = "ifn2"
    sigma: float = 25.0
    q: float = 0.10
    subtype: str | None = None


@app.post("/api/signature/localization", response_class=ORJSONResponse)
def signature_localization(req: SigLocReq):
    """Do two signatures localise to different regions / cell types? (see sigstats.localization)"""
    ga, gb = _genes_of(req.a, None), _genes_of(req.b, None)
    if not 0 <= req.sigma <= MAX_SIGMA or not 0.01 <= req.q <= 0.5:
        raise HTTPException(400, "bad sigma or q")
    key = ("loc", req.a, req.b, req.sigma, req.q, req.subtype)
    if key not in _sig_cache:
        _sig_cache[key] = sigstats.localization(S, SM, ga, gb, req.sigma, req.q, req.subtype)
        _sig_cache[key].update(a_id=req.a, b_id=req.b)
    return _sig_cache[key]


class ViolinFeature(BaseModel):
    sig: str | None = None
    gene: str | None = None


class ViolinReq(BaseModel):
    features: list[ViolinFeature]
    cell_types: list[int] = [0]             # cell_type level ids; -1 = all cells
    sigma: float = 0.0
    scope: str = "none"
    nonzero_only: bool = False
    layer: str = "logcpx"                   # logcpx | counts (same as the main viewer)
    unit: str = "sample"                    # sample | subtype | type: one violin per sample, or cells pooled per disease subtype / type
    level: str = "cell_type"                # cell_type | cell_state (one panel per cell type, one group per cell state of it)


@app.post("/api/violin", response_class=ORJSONResponse)
def violin_api(req: ViolinReq):
    """Per-sample score distributions (one violin per sample) for each selected gene / signature and cell type."""
    if not req.features or len(req.features) > 8 or not req.cell_types or len(req.cell_types) > 12:
        raise HTTPException(400, "choose 1-8 features and 1-12 cell types")
    if req.scope not in SCOPES or not 0 <= req.sigma <= MAX_SIGMA:
        raise HTTPException(400, "bad scope or sigma")
    if req.layer not in ("logcpx", "counts"):
        raise HTTPException(400, "layer must be logcpx or counts")
    if req.unit not in ("sample", "subtype", "type"):
        raise HTTPException(400, "unit must be sample, subtype or type")
    if req.level not in ("cell_type", "cell_state"):
        raise HTTPException(400, "level must be cell_type or cell_state")
    nct = len(S.cats["cell_type"]["levels"])
    if any(c < -1 or c >= nct for c in req.cell_types):
        raise HTTPException(400, "unknown cell type")
    feats = []
    for f in req.features:
        if f.sig:
            feats.append({"label": S.sig_by_id[f.sig]["name"] if f.sig in S.sig_by_id else f.sig, "genes": _genes_of(f.sig, None)})
        elif f.gene:
            feats.append({"label": f.gene, "genes": _genes_of(None, [f.gene])})
        else:
            raise HTTPException(400, "each feature needs a gene or a signature id")
    key = ("vio", tuple((x["label"], tuple(x["genes"])) for x in feats), tuple(req.cell_types), req.sigma, req.scope, req.nonzero_only, req.level, req.unit, req.layer)
    if key not in _sig_cache:
        try:
            _sig_cache[key] = {"panels": violin.violins(S, SM, feats, req.cell_types, req.sigma, req.scope, layer=req.layer, nonzero_only=req.nonzero_only, level=req.level, unit=req.unit),
                               "samples": [{"id": d, "subtype": sigstats.subtype_of(SM, smp[0]), "disease_type": SM.info[smp[0]]["disease_type"]} for d, smp in sigstats.donors_of(S).items()]}
        except ValueError as e:
            raise HTTPException(400, str(e))
    return _sig_cache[key]


@app.post("/api/violin/hist", response_class=ORJSONResponse)
def violin_hist(req: ViolinReq):
    """Histograms of the per-cell values behind the violins (all donors pooled), so the plotted scale can be checked."""
    if not req.features or len(req.features) > 8 or not req.cell_types or len(req.cell_types) > 12:
        raise HTTPException(400, "choose 1-8 features and 1-12 cell types")
    if req.scope not in SCOPES or not 0 <= req.sigma <= MAX_SIGMA or req.layer not in ("logcpx", "counts"):
        raise HTTPException(400, "bad scope, sigma or layer")
    nct = len(S.cats["cell_type"]["levels"])
    if any(c < -1 or c >= nct for c in req.cell_types) or req.unit not in ("sample", "subtype", "type"):
        raise HTTPException(400, "unknown cell type or unit")
    feats = [{"label": S.sig_by_id[f.sig]["name"] if f.sig in S.sig_by_id else f.sig, "genes": _genes_of(f.sig, None)} if f.sig else {"label": f.gene, "genes": _genes_of(None, [f.gene])}
             for f in req.features if f.sig or f.gene]
    key = ("viohist", tuple((x["label"], tuple(x["genes"])) for x in feats), tuple(req.cell_types), req.sigma, req.scope, req.nonzero_only, req.layer, req.unit)
    if key not in _sig_cache:
        _sig_cache[key] = {"panels": violin.histograms(S, SM, feats, req.cell_types, req.sigma, req.scope, layer=req.layer, nonzero_only=req.nonzero_only, unit=req.unit), "norm_scale": S.norm_scale}
    return _sig_cache[key]


# ---------------------------------------------------------------- glom heterogeneity
def _gh():
    if GH is None or not GH.ok:
        raise HTTPException(404, "this project has no ROI heterogeneity data (a region_heterogeneity object, or gloms/gloms.json + gloms.npz)")
    return GH


@app.get("/api/glomhet/info", response_class=ORJSONResponse)
def glomhet_info():
    """Glom table (ids, sample, disease, activity / chronicity, cell counts, clusters) and the PC / UMAP embeddings (first 10 PCs)."""
    return _gh().info()


class GlomFilterReq(BaseModel):
    filters: dict[str, list[int]] = {}
    hidden: dict[str, list[int]] = {}
    expr: list[ExprReq] = []


@app.post("/api/glomhet/filter", response_class=ORJSONResponse)
def glomhet_filter(req: GlomFilterReq):
    """Cells of each glom (viewer cells) that pass the viewer's filters, in the order of /api/glomhet/info."""
    try:
        return {"n_pass": _gh().passing(req.filters, req.hidden, _expr(req.expr)).tolist()}
    except KeyError as e:
        raise HTTPException(400, str(e))


class GlomColorReq(BaseModel):
    kind: str                      # gene | signature | type | state
    gene: str | None = None
    sig: str | None = None
    genes: list[str] | None = None  # custom signature (panel genes) instead of a curated `sig`
    level: int | None = None       # cell type / state level id


@app.post("/api/glomhet/color", response_class=ORJSONResponse)
def glomhet_color(req: GlomColorReq):
    """One value per glom: glom-level log-normalised expression (mean over genes) or the fraction of the glom's cells of a type / state."""
    g = _gh()
    try:
        if req.kind == "signature":
            v, info = g.color("signature", _genes_of(req.sig, req.genes))
        elif req.kind == "gene":
            v, info = g.color("gene", _genes_of(None, [req.gene]) or [req.gene])
        elif req.kind in ("type", "state"):
            n = len(g.type_levels if req.kind == "type" else g.state_levels)
            if req.level is None or not 0 <= req.level < n:
                raise HTTPException(400, "bad level")
            v, info = g.color(req.kind, level=req.level)
        else:
            raise HTTPException(400, "kind must be gene, signature, type or state")
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"values": np.round(v, 4).tolist(), **info}


class GlomCellScoreReq(BaseModel):
    genes: list[str]               # panel genes (a signature = its genes); mean over genes, like the main window
    layer: str = "logcpx"
    sigma: float = 0.0
    scope: str = "none"
    within: str = "tissue"         # tissue | glom: smoothing and z-score over the whole tissue or over glomerulus cells only


@app.post("/api/glomhet/cellscore", response_class=ORJSONResponse)
def glomhet_cellscore(req: GlomCellScoreReq):
    """One value per glom (order of /api/glomhet/info): mean of the cell-level, optionally smoothed / z-scored score over the glom's cells."""
    g = _gh()
    bad = [x for x in req.genes if x not in S.gene_index]
    if not req.genes or bad:
        raise HTTPException(400, f"give panel genes (not on the panel: {bad[:5]})")
    if req.layer not in ("logcpx", "counts") or req.scope not in SCOPES or not 0 <= req.sigma <= MAX_SIGMA:
        raise HTTPException(400, "bad layer, scope or sigma")
    try:
        v = g.cell_score(req.genes, req.layer, req.sigma, req.scope, req.within)
    except ValueError as e:
        raise HTTPException(400, str(e))
    return {"values": np.round(v, 4).tolist(), "n_missing": int(np.isnan(v).sum())}


class GlomCorrReq(BaseModel):
    kind: str = "gene"             # gene | type | state
    layer: str = "harmony"         # harmony | raw (uncorrected PCA); only the ones the object has
    glom_idx: list[int]            # gloms to include (those currently shown)
    n_top: int = 10
    sort_pc: int | None = None     # genes: only the genes most correlated (+ / -) with this PC (0-based)


@app.post("/api/glomhet/corr", response_class=ORJSONResponse)
def glomhet_corr(req: GlomCorrReq):
    g = _gh()
    if req.kind not in ("gene", "type", "state") or req.layer not in g.emb or not 1 <= req.n_top <= 100 or (req.sort_pc is not None and not 0 <= req.sort_pc < 10):
        raise HTTPException(400, "bad kind, layer or n_top")
    if any(i < 0 or i >= g.n for i in req.glom_idx):
        raise HTTPException(400, "bad glom index")
    try:
        return g.correlations(req.glom_idx, req.layer, req.kind, ntop=req.n_top, sort_pc=req.sort_pc)
    except ValueError as e:
        raise HTTPException(400, str(e))



@app.get("/api/umap/{name}")
def umap(name: str):
    """UMAP coordinates of ALL cells (float32, n x 2, cache cell order; a sample's cells start at its `offset` in /api/meta)."""
    if name not in S.umap:
        raise HTTPException(404, f"no UMAP {name!r}")
    return raw(np.asarray(S.umap[name], dtype=np.float32))


@app.get("/api/pca/{x}/{y}")
def pca(x: int, y: int):
    """Two Harmony-corrected PCs (1-based) of ALL cells as float32, n x 2, cache cell order."""
    if S.pca is None or not (1 <= x <= S.pca.shape[1] and 1 <= y <= S.pca.shape[1]):
        raise HTTPException(404, "no such PC")
    return raw(np.stack([S.pca[:, x - 1], S.pca[:, y - 1]], 1).astype(np.float32))


@app.get("/api/cat_all/{key:path}")
def cat_all(key: str):
    """Categorical codes (uint16) of ALL cells in cache order, to colour / filter the UMAP."""
    if key not in S.cats:
        raise HTTPException(404, f"unknown categorical field {key!r}")
    return raw(S.cats[key]["codes"])


class ExprMaskReq(BaseModel):
    expr: list[ExprReq]


_mask_cache: dict = {}


@app.post("/api/expr_mask")
def expr_mask(req: ExprMaskReq):
    """uint8 per cell (1 = passes every gene / signature filter) for ALL cells in cache order, to filter the UMAP."""
    ex = _expr(req.expr)
    key = repr(ex)
    if key not in _mask_cache:
        if len(_mask_cache) >= 8:
            _mask_cache.pop(next(iter(_mask_cache)))
        out = np.zeros(len(S.df), np.uint8)
        for s, (a, b) in S.slices.items():
            out[a:b] = NB.expr_mask(s, ex)
        _mask_cache[key] = out
    return raw(_mask_cache[key])


@app.get("/api/sample/{sample}/poly/{part}")
def poly(sample: str, part: str):
    """Cell boundary polygons: part=offsets (Int32, n+1 vertex offsets) or xy (Float32 interleaved)."""
    a, b = sl(sample)
    if S.poly_off is None or part not in ("offsets", "xy"):
        raise HTTPException(404, "no polygon data" if S.poly_off is None else f"unknown part {part!r}")
    off, xy = S.polygons(a, b)
    return raw(off if part == "offsets" else xy)


@app.get("/api/sample/{sample}/cat/{key:path}")
def cat(sample: str, key: str):
    a, b = sl(sample)
    if key not in S.cats:
        raise HTTPException(404, f"unknown categorical field {key!r}")
    return raw(S.cats[key]["codes"][a:b])


@app.get("/api/sample/{sample}/num/{key:path}")
def num(sample: str, key: str):
    a, b = sl(sample)
    if key not in S.numeric:
        raise HTTPException(404, f"unknown numeric field {key!r}")
    return raw(S.df[S.numeric[key]].values[a:b].astype(np.float32))


@app.get("/api/sample/{sample}/ids", response_class=ORJSONResponse)
def ids(sample: str):
    a, b = sl(sample)
    return S.df[S.map["cell_id"]].values[a:b].tolist()


class NeighborReq(BaseModel):
    field: str = "cell_type"
    method: str = "radius"                 # radius | knn | delaunay
    distance: float = 25.0                 # µm: radius / max neighbour distance (knn) / max Delaunay edge length
    k: int = 6                             # knn: number of nearest cells
    steps: int = 1                         # delaunay: graph steps
    filters: dict[str, list[int]] = {}     # field -> allowed level ids (AND across fields)
    hidden: dict[str, list[int]] = {}      # field -> level ids to exclude
    expr: list[ExprReq] = []               # gene / signature filters (AND with the above)


@app.post("/api/neighbors", response_class=ORJSONResponse)
def neighbors(req: NeighborReq):
    """Neighbour composition of the cells passing `filters`, per sample (see neighbors.py)."""
    try:
        return NB.compute(req.field, req.method, req.distance, req.k, req.steps, req.filters, req.hidden, _expr(req.expr))
    except KeyError as e:
        raise HTTPException(400, str(e))


class NeighborDistReq(NeighborReq):
    groups: dict[str, dict[str, str]] = {}  # grouping level -> {sample: group name}


@app.post("/api/neighbors/dist", response_class=ORJSONResponse)
def neighbors_dist(req: NeighborDistReq):
    """Per-index-cell spread of the neighbour composition and bootstrap enrichment intervals, per group (see Neighbors.per_cell)."""
    try:
        return NB.per_cell(req.field, req.method, req.distance, req.k, req.steps, req.filters, req.hidden, _expr(req.expr), req.groups)
    except KeyError as e:
        raise HTTPException(400, str(e))


class NeighborhoodReq(BaseModel):
    sample: str
    method: str = "radius"
    distance: float = 25.0
    k: int = 6
    steps: int = 1
    filters: dict[str, list[int]] = {}
    hidden: dict[str, list[int]] = {}
    expr: list[ExprReq] = []
    limit: int = 48
    order: str = "random"                  # random | most | fewest (neighbours)
    seed: int = 1


@app.post("/api/neighborhoods", response_class=ORJSONResponse)
def neighborhoods(req: NeighborhoodReq):
    """Index cells of one sample (those passing `filters`) with their neighbour cell indices, for drawing in space."""
    try:
        return NB.neighbourhoods(req.sample, req.filters, req.hidden, req.method, req.distance, req.k, req.steps,
                                 max(1, min(req.limit, 200)), req.order, req.seed, _expr(req.expr))
    except KeyError as e:
        raise HTTPException(400, str(e))


@app.get("/api/sample/{sample}/expr")
def expr(sample: str, genes: str, layer: str = "logcpx", sigma: float = 0.0, scope: str = "none", within: str = "tissue"):
    """Mean expression of `genes` per cell, optionally Gaussian-smoothed (sigma in µm) and z-scored within `scope`.
    `within` limits smoothing and the z-score reference to the whole tissue, the glomeruli (each smoothed alone) or one region level
    (`cat:<field>:<level index>`); cells outside it are NaN."""
    sl(sample)
    gl = [g for g in genes.split(",") if g]
    if not gl or len(gl) > 500:
        raise HTTPException(400, "give 1-500 comma-separated genes")
    bad = [g for g in gl if g not in S.gene_index]
    if bad:
        raise HTTPException(404, f"genes not in panel: {bad[:5]}")
    if layer not in ("logcpx", "counts"):
        raise HTTPException(400, "layer must be 'logcpx' or 'counts'")
    if scope not in SCOPES:
        raise HTTPException(400, f"scope must be one of {SCOPES}")
    if not 0 <= sigma <= MAX_SIGMA:
        raise HTTPException(400, f"sigma must be in [0, {MAX_SIGMA}] µm")
    try:
        return raw(SM.expression(sample, gl, layer, sigma, scope, within))
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.get("/api/sample/{sample}/gene/{gene}")
def gene(sample: str, gene: str, layer: str = "logcpx"):
    a, b = sl(sample)
    if gene not in S.gene_index:
        raise HTTPException(404, f"gene {gene!r} not in panel")
    if layer not in ("logcpx", "counts"):
        raise HTTPException(400, "layer must be 'logcpx' or 'counts'")
    key = (sample, gene, layer)
    if key not in _gene_cache:
        if len(_gene_cache) > 400:
            _gene_cache.pop(next(iter(_gene_cache)))
        _gene_cache[key] = S.gene_vector(gene, a, b, layer)
    return raw(_gene_cache[key])


@app.get("/api/sample/{sample}/tx")
def transcripts(sample: str, genes: str, x0: float = -1e12, y0: float = -1e12, x1: float = 1e12, y1: float = 1e12, max: int = 200000):
    """Transcript positions of up to 8 genes in a box (cell coordinates, y not flipped). Body: float32 (n, 3) = x, y, index of the gene in `genes`;
    header X-Counts = {"genes": [...], "in_view": [...], "total": [...]}."""
    sl(sample)
    if TX is None:
        raise HTTPException(404, "this project has no transcripts")
    names = [g for g in genes.split(",") if g]
    if not names or len(names) > 8:
        raise HTTPException(400, "give 1-8 genes")
    bad = [g for g in names if g not in S.gene_index]
    if bad:
        raise HTTPException(404, f"genes not in panel: {bad[:5]}")
    if not 1000 <= max <= 1_000_000:
        raise HTTPException(400, "max must be 1000-1000000")
    arr, in_view, total = TX.query(sample, [S.gene_index[g] for g in names], (x0, y0, x1, y1), max)
    return Response(np.ascontiguousarray(arr).tobytes(), media_type="application/octet-stream",
                    headers={"X-Counts": json.dumps({"genes": names, "in_view": in_view, "total": total}), "Access-Control-Expose-Headers": "X-Counts"})


@app.get("/api/genes/search", response_class=ORJSONResponse)
def genes(q: str = "", limit: int = 30):
    return S.search_genes(q, limit) if q else S.genes[:limit]


@app.get("/api/genes/resolve", response_class=ORJSONResponse)
def resolve_genes(names: str = ""):
    """Case-insensitive match of comma-separated gene symbols to panel genes: {found: [...], missing: [...]}."""
    return S.resolve_genes([n for n in names.split(",") if n.strip()])


@app.get("/api/cell/{sample}/{i}", response_class=ORJSONResponse)
def cell(sample: str, i: int, top: int = 40):
    """All metadata + the cell's expressed genes (counts and logcpx)."""
    a, b = sl(sample)
    if not 0 <= i < b - a:
        raise HTTPException(404, "cell index out of range")
    gi = a + i
    row = S.df.iloc[gi]
    md = {k: (None if (isinstance(v, float) and np.isnan(v)) else (v.item() if hasattr(v, "item") else v))
          for k, v in row.items()}
    # expression of one cell needs a column of the gene-major CSR: scan genes (4k) with searchsorted
    out = []
    for g, name in enumerate(S.genes):
        p0, p1 = int(S.indptr[g]), int(S.indptr[g + 1])
        j = np.searchsorted(S.indices[p0:p1], gi)
        if j < p1 - p0 and S.indices[p0 + j] == gi:
            c = float(S.data[p0 + j])
            out.append({"gene": name, "counts": c, "logcpx": float(np.log1p(c / max(S.lib[gi], 1) * S.norm_scale))})
    out.sort(key=lambda r: -r["counts"])
    return {"metadata": md, "expression": out, "n_genes_detected": len(out)}


@app.get("/api/he", response_class=ORJSONResponse)
def he_index():
    """Registered H&E images per sample (pre-tiled pyramids in the project's he/ folder): placement in um, size, pyramid levels. {} when none."""
    f = PROJ / "he" / "index.json"
    out = json.loads(f.read_text()) if f.exists() else {}
    if HE_EXTRA is not None:                                     # registered images given when the project was opened (kept in the cache, not in the project folder)
        out.update(json.loads((HE_EXTRA / "he" / "index.json").read_text()))
    return out


_HE_KEY = re.compile(r"^[A-Za-z0-9_.\-]+$")


_HE_CACHE = {}


def _he_multi():
    """The merged project's he_slide.json when it lists several scans ({multi, scans: {key: {dir, source, runs}}, runs: {run: key}}), else None."""
    f = PROJ / "he_slide.json"
    if not f.exists():
        return None
    k = (str(f), f.stat().st_mtime_ns)
    if _HE_CACHE.get("k") != k:
        d = json.loads(f.read_text())
        _HE_CACHE.clear()
        _HE_CACHE.update({"k": k, "m": d if d.get("multi") else None, "tiles": None})
    return _HE_CACHE["m"]


def _he_tile_dirs():
    """{pyramid key: folder holding its he/<key>/ tiles} over all scans of a multi-scan project (read once)."""
    m = _he_multi()
    if m and _HE_CACHE.get("tiles") is None:
        t = {}
        for sk in m["scans"]:
            d = _he_scan_dir(sk)
            try:
                for sc in json.loads((d / he_slide.INDEX).read_text())["scenes"]:
                    t[sc["key"]] = d
            except Exception:
                pass
        _HE_CACHE["tiles"] = t
    return (_HE_CACHE.get("tiles") or {}) if m else {}


def _he_scan_dir(key):
    m = _he_multi()
    if not m or key not in m["scans"] or not _HE_KEY.match(key):
        raise HTTPException(404, "no such H&E scan in this project")
    return (PROJ / m["scans"][key]["dir"]).resolve()


@app.get("/api/heslide", response_class=ORJSONResponse)
def heslide_index(scan: str = ""):
    """H&E slide scan (unregistered): scenes with their pyramids (placement in um of the slide frame) and the annotation summary. {} when none.
    A project that holds several scans (batch import with an H&E folder) answers with the list of scans and the run -> scan map; ?scan=<key> gives one."""
    m = _he_multi()
    if m is None:
        f = PROJ / "he_slide.json"
        return json.loads(f.read_text()) if f.exists() else {}
    if not scan:
        return {"multi": True, "scans": [{"key": k, **{q: v[q] for q in ("source", "n_scenes", "annotations", "runs")}} for k, v in m["scans"].items()], "runs": m["runs"]}
    return json.loads((_he_scan_dir(scan) / he_slide.INDEX).read_text())


@app.get("/api/heslide/annotations")
def heslide_annotations(scan: str = ""):
    f = (_he_scan_dir(scan) if _he_multi() is not None else PROJ) / he_slide.ANN
    if not f.is_file():
        raise HTTPException(404, "this project has no H&E annotations")
    return FileResponse(f, media_type="application/json")


@app.get("/api/he/{key}/{z}/{x}_{y}.jpg")
def he_tile(key: str, z: int, x: int, y: int):
    if not _HE_KEY.match(key) or min(z, x, y) < 0:
        raise HTTPException(404)
    f = PROJ / "he" / key / str(z) / f"{x}_{y}.jpg"
    if not f.is_file() and HE_EXTRA is not None:
        f = HE_EXTRA / "he" / key / str(z) / f"{x}_{y}.jpg"
    if not f.is_file():                                       # a batch project keeps each scan's tiles in <out>/_he/<scan>/he/<key>/
        d = _he_tile_dirs().get(key)
        f = d / "he" / key / str(z) / f"{x}_{y}.jpg" if d else None
        if f is None or not f.is_file():
            raise HTTPException(404)
    return FileResponse(f, media_type="image/jpeg", headers={"Cache-Control": "public, max-age=86400"})


@app.get("/api/if", response_class=ORJSONResponse)
def if_index():
    """IF image channels per sample, in the shape of /api/he (key, w, h, px, x0, y0, tsize, nz) plus name and the default intensity window (lo, hi). {} when none."""
    if not IFX:
        return {}
    out = {}
    for smp in IFX["samples"]:
        lst = if_image.sample_list(IFX, smp)
        if lst:
            out[smp] = lst
    return out


@app.get("/api/if/{key}/{ch}/{z}/{x}_{y}.jpg")
def if_tile(key: str, ch: int, z: int, x: int, y: int, c: str = "white", lo: float = 0, hi: float = 0, g: float = 1.0):
    ent = IFX["images"].get(key) if IFX else None
    if ent is None or not 0 <= ch < len(ent["channels"]) or min(z, x, y) < 0:
        raise HTTPException(404)
    chn = ent["channels"][ch]
    lo, hi = (lo, hi) if hi > lo else (chn["lo"], chn["hi"])
    try:
        b = if_image.render(chn.get("file") or ent["file"], chn["ifd"], z, x, y, c, lo, hi, min(max(g, 0.2), 5.0))
    except (OSError, if_image.IFError) as e:
        raise HTTPException(404, f"IF image not readable: {e}")
    if b is None:
        raise HTTPException(404)
    return Response(b, media_type="image/jpeg", headers={"Cache-Control": "public, max-age=3600"})


# ---------------------------------------------------------------- IMC (.mcd) beside the cells
class IMCState:
    state = "idle"           # idle | converting | error
    progress = 0.0
    msg = ""
    error = None
    mcd = None


IMS = IMCState()


def _imc_links():
    f = importer.cache_root() / "imc" / "_links.json"
    try:
        return f, json.loads(f.read_text())
    except (OSError, ValueError):
        return f, {}


def _imc_link(project, store=None):
    """Remember (store given) or look up which converted IMC store belongs to a project folder."""
    f, d = _imc_links()
    k = str(Path(project).resolve())
    if store is None:
        return d.get(k)
    d[k] = str(store)
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(d))


def _imc_open_store(project):
    global IMX
    IMX = None
    s = _imc_link(project)
    ix = imc.read_store(s) if s else None
    if ix:
        IMX = {"store": Path(s), "index": ix}


def _imc_map_file():
    return IMX["store"] / ("map-" + __import__("hashlib").sha1(str(PROJ.resolve()).encode()).hexdigest()[:10] + ".json")


def _imc_map():
    try:
        return json.loads(_imc_map_file().read_text())
    except (OSError, ValueError):
        return {}


def _imc_start(path):
    """Convert an .mcd file into the cache (background thread) and attach it to the open project when it is done."""
    global IMX
    if IMS.state == "converting":
        raise HTTPException(409, "an IMC conversion is already running")
    try:
        info = imc.scan(path)
    except imc.IMCError as e:
        raise HTTPException(400, str(e))
    proj = PROJ
    IMS.state, IMS.progress, IMS.msg, IMS.error, IMS.mcd = "converting", 0.0, "starting", None, info["file"]

    def work():
        global IMX
        try:
            d = imc.convert(info["file"], importer.cache_root(), _log, lambda f, m="": (setattr(IMS, "progress", float(f)), setattr(IMS, "msg", m)), info)
            if PROJ is proj:
                _imc_link(proj, d)
                _imc_open_store(proj)
            IMS.state = "idle"
            _log("IMC ready")
        except Exception as e:
            traceback.print_exc()
            IMS.state, IMS.error = "error", str(e)
            _log(f"IMC error: {e}")
    threading.Thread(target=work, daemon=True).start()


class IMCAddReq(BaseModel):
    path: str


@app.post("/api/imc/add", response_class=ORJSONResponse)
def imc_add(r: IMCAddReq):
    p = Path(r.path.strip()).expanduser()
    if not p.is_file() or p.suffix.lower() != ".mcd":
        raise HTTPException(400, f"not an .mcd file: {p}")
    _imc_start(p)
    return {"ok": True}


@app.get("/api/imc", response_class=ORJSONResponse)
def imc_info():
    """The IMC slide of the open project: ROIs (acquisitions) with size, slide position and whether they are complete, the channels, the ROI chosen
    for each tissue, and the conversion state."""
    st = {"state": IMS.state, "progress": IMS.progress, "msg": IMS.msg, "error": IMS.error}
    if not IMX:
        return {"has": False, "status": st}
    ix = IMX["index"]
    mp = _imc_map()
    first = complete = None
    for r in ix["rois"]:
        try:
            m = json.loads((IMX["store"] / r["key"] / "meta.json").read_text())
        except (OSError, ValueError):
            continue
        r["hi"] = m["hi"][0]
        first = first or r["key"]
        if r["complete"] and not complete:
            complete = r["key"]
    return {"has": True, "status": st, "file": ix["file"], "rois": ix["rois"], "channels": ix["channels"], "map": mp, "default": complete or first}


class IMCMapReq(BaseModel):
    sample: str
    roi: str


@app.post("/api/imc/map", response_class=ORJSONResponse)
def imc_map(r: IMCMapReq):
    if not IMX or r.roi not in {x["key"] for x in IMX["index"]["rois"]}:
        raise HTTPException(404, "unknown ROI")
    m = _imc_map()
    m[r.sample] = r.roi
    _imc_map_file().write_text(json.dumps(m))
    return {"ok": True}


def _imc_roi(roi):
    if not IMX:
        raise HTTPException(404, "no IMC data")
    for r in IMX["index"]["rois"]:
        if r["key"] == roi and (IMX["store"] / roi / "meta.json").is_file():
            return r
    raise HTTPException(404, "unknown ROI")


@app.get("/api/imc/index", response_class=ORJSONResponse)
def imc_scene(roi: str):
    """One ROI in the shape of the H&E slide scan (/api/heslide?scan=...), so the side pane shows it like a scan."""
    r = _imc_roi(roi)
    m = json.loads((IMX["store"] / roi / "meta.json").read_text())
    note = "" if m["complete"] else f" (incomplete: {m['rows']:,} of {m['h']:,} rows acquired)"
    return {"source": f"IMC {r['name']}{note}", "bounds_um": [0, 0, m["w"], m["h"]],
            "scenes": [{"key": f"imc:{roi}", "tsize": m["tile"], "px": 1.0, "nz": m["nz"], "x0": 0, "y0": 0, "w": m["w"], "h": m["h"]}]}


def _imc_args(ch, c, b):
    try:
        chs = [int(x) for x in ch.split(",") if x != ""]
        bs = [float(x) for x in b.split(",")] if b else []
    except ValueError:
        raise HTTPException(400, "bad channel list")
    cols = c.split(",") if c else []
    n = len(IMX["index"]["channels"])
    if not chs or len(chs) > 16 or any(not 0 <= x < n for x in chs):
        raise HTTPException(400, "choose 1 to 16 channels")
    cols = (cols + ["white"] * len(chs))[:len(chs)]
    bs = (bs + [1.0] * len(chs))[:len(chs)]
    return chs, cols, bs


@app.get("/api/imc/tile/{roi}/{z}/{x}_{y}.jpg")
def imc_tile(roi: str, z: int, x: int, y: int, ch: str = "", c: str = "", b: str = "", g: float = 1.0):
    _imc_roi(roi)
    chs, cols, bs = _imc_args(ch, c, b)
    out = imc.tile(IMX["store"], roi, z, x, y, chs, cols, bs, min(max(g, 0.2), 5.0))
    if out is None:
        raise HTTPException(404)
    return Response(out, media_type="image/jpeg", headers={"Cache-Control": "public, max-age=3600"})


@app.get("/api/imc/fit/{roi}.jpg")
def imc_fit(roi: str, w: int = 300, h: int = 300, ch: str = "", c: str = "", b: str = "", g: float = 1.0):
    _imc_roi(roi)
    chs, cols, bs = _imc_args(ch, c, b)
    out = imc.fit(IMX["store"], roi, min(max(w, 32), 1600), min(max(h, 32), 1600), chs, cols, bs, min(max(g, 0.2), 5.0))
    return Response(out, media_type="image/jpeg", headers={"Cache-Control": "public, max-age=3600"})


# ---------------------------------------------------------------- opening a project
class State:
    state = "idle"          # idle | importing | ready | error
    log: list = []
    error = None
    project = None


ST = State()


def _log(msg):
    ST.log.append(msg)
    PROF.mark(msg)
    print("[webbscope]", msg, flush=True)


def _load(cache_dir, project):
    """Build the engine objects for an imported project."""
    global S, SM, NB, GH, PROJ, EMB_LABEL, TX, DEFAULT_GREY, IFX, HE_EXTRA
    with PROF.step("load: Store (counts, cells, polygons from the cache)"):
        st = Store(cache_dir, Path(cache_dir) / "mapping.json")
    with PROF.step("load: Smoother, Neighbors, GlomHet"):
        sm = Smoother(st)
        nb = Neighbors(st)
        nb.SM = sm
        gh = GlomHet(st, nb, sm)
    _gene_cache.clear(); _sig_cache.clear(); _mask_cache.clear()
    S, SM, NB, GH, PROJ = st, sm, nb, gh, Path(project)
    _man = json.loads((Path(cache_dir) / "cache_info.json").read_text())["manifest"]
    EMB_LABEL = _man.get("embedding_label")
    DEFAULT_GREY = bool(_man.get("default_grey"))
    HE_EXTRA = None
    _d = Path(cache_dir) / "he_registered"
    if (_d / "he" / "index.json").is_file():                    # registered H&E added earlier when this project was opened
        HE_EXTRA = _d
    IFX = None
    try:
        IFX = if_image.read_index(project)
    except Exception as e:                           # a damaged if.json must not stop the project from opening
        _log(f"IF images ignored: {e}")
    _imc_open_store(project)
    TX = None
    if (Path(project) / "transcripts" / "index.json").exists():
        try:
            with PROF.step("load: transcripts index"):
                TX = Transcripts(Path(project) / "transcripts", list(st.slices), len(st.genes))
        except Exception as e:                       # a damaged / mismatching transcript folder must not stop the project from opening
            _log(f"transcripts ignored: {e}")


def open_project(path, force=False, prepare=None, roi_object=None, if_path=None, cols=None, he_reg=None, imc_path=None):
    """Import (if needed) and load a project folder, in a background thread."""
    global S
    if ST.state == "importing":
        raise HTTPException(409, "an import is already running")
    ST.state, ST.log, ST.error, ST.project = "importing", [], None, str(path)
    PROF.start_open(path)

    def work():
        global S
        try:
            if prepare:
                with PROF.step("prepare (convert / merge before opening)"):
                    p = Path(prepare()).expanduser()
            else:
                p = Path(path).expanduser()
            if not p.is_dir():
                raise importer.ImportError_(f"not a folder: {p}")
            with PROF.step("import into the cache (skipped when current)"):
                cdir = importer.ensure_cache(p, _log, force, roi_object, cols)
            S = None                                  # until the new project is fully loaded
            _log("loading")
            _load(cdir, p)
            if if_path:
                _user_if(if_path)
            if he_reg:
                _register_he(he_reg, cdir)
            if imc_path:
                try:
                    _imc_start(imc_path)
                except HTTPException as e:
                    _log(f"IMC ignored: {e.detail}")
            ST.state = "ready"
            _log("ready")
        except importer.ImportError_ as e:
            ST.state, ST.error = "error", str(e); _log("error: " + str(e))
        except Exception as e:                        # unexpected: full traceback on the console, short message in the page
            traceback.print_exc()
            ST.state, ST.error = "error", f"unexpected error: {e}"
    threading.Thread(target=work, daemon=True).start()


def _register_he(he_reg, cdir):
    """Registered H&E images (a folder of OME-TIFFs, optional index.json) for the samples of the project that was just opened. Tiles go to the cache, not into the project."""
    global HE_EXTRA
    try:
        out = Path(cdir) / "he_registered"
        _log("registered H&E: tiling the images (once; later opens reuse the tiles)")
        he_registered.run(he_reg["folder"], he_reg.get("index") or "", PROJ, _log, out_root=out, samples=set(S.samples))
        if (out / "he" / "index.json").is_file():
            HE_EXTRA = out
    except (importer.ImportError_, OSError, ValueError, KeyError) as e:
        _log(f"registered H&E ignored: {e}")


def _he_reg_args(folder, index):
    """Validated (folder, index) of the registered-H&E fields, or None when the folder is empty; HTTP 400 when a path is wrong."""
    folder, index = (folder or "").strip(), (index or "").strip()
    if not folder:
        return None
    f = Path(folder).expanduser()
    if not f.is_dir():
        raise HTTPException(400, f"registered H&E folder not found: {f}")
    i = Path(index).expanduser() if index else None
    if i is not None and not i.is_file():
        raise HTTPException(400, f"index.json not found: {i}")
    return {"folder": str(f), "index": str(i) if i else ""}


def _user_if(path):
    """An IF image given on the open page (counts + metadata import): one image for every sample. It has to be in the cells' coordinates (the
    Xenium morphology image is); a quick check refuses an image that is smaller than the cells' extent."""
    global IFX
    try:
        p = Path(path)
        ent = if_image.describe_many([p] if p.is_file() else _pick_in_folder(p), None)
        W, Hh = ent["w"] * ent["px"], ent["h"] * ent["px"]
        mx, my = float(np.nanmax(S.x)), float(np.nanmax(S.y))
        if mx > W * 1.02 or my > Hh * 1.02:
            raise if_image.IFError(f"the image covers {W:.0f} x {Hh:.0f} um but the cells extend to {mx:.0f} x {my:.0f} um: it is not in the cells' coordinates")
        IFX = {"images": {"user": ent}, "samples": {smp: "user" for smp in S.slices}}
        _log(f"IF image: {len(ent['channels'])} channel(s) ({', '.join(c['name'] for c in ent['channels'])}; {ent['w']:,} x {ent['h']:,} px, {ent['px']} um/px), read in place; assumed to be in the cells' coordinates (only an image smaller than the cells is refused)")
    except Exception as e:
        _log(f"IF image ignored: {e}")


def _pick_in_folder(p):
    q = if_image.pick_all(xenium.find_in_folder(p))
    if not q:
        raise if_image.IFError(f"no morphology_focus / morphology_mip .ome.tif in {p}")
    return q


class OpenReq(BaseModel):
    path: str
    force: bool = False
    roi_object: str = ""           # region_heterogeneity object (.rds), optional
    if_image: str = ""             # IF image (.ome.tif or a Xenium output folder), optional
    he_reg_folder: str = ""        # folder of registered H&E OME-TIFFs (one per tissue), optional
    he_reg_index: str = ""         # index.json mapping image names to sample ids, optional
    imc_file: str = ""             # IMC slide (.mcd), optional: converted into the cache in the background
    sample_col: str | None = None  # cell-table column that identifies the samples (None = project.json / found by name)
    donor_col: str | None = None   # cell-table column that identifies the donors ("" = there is none, None = found by name)


@app.get("/api/status", response_class=ORJSONResponse)
def status():
    return {"state": ST.state, "log": ST.log[-40:], "error": ST.error, "project": ST.project, "open": S is not None}


@app.get("/api/columns", response_class=ORJSONResponse)
def project_columns(path: str):
    """Columns of a project folder's cell table, with the ones that look like the sample id and the donor id (by name), for the question asked when opening.
    ask = False for projects this app wrote itself (their columns are known)."""
    import pyarrow.parquet as pq, pyarrow.compute as pc
    from .engine.config import find_columns, suggest
    p = Path(path).expanduser()
    f = p / "cells.parquet"
    if not f.is_file():
        raise HTTPException(400, "this folder has no cells.parquet: it is not a project folder")
    names = pq.ParquetFile(f).schema_arrow.names
    roles = {}
    if (p / "project.json").is_file():
        try:
            roles = json.loads((p / "project.json").read_text()).get("roles") or {}
        except ValueError:
            pass
    sug = suggest(names)
    out = {"columns": names, "ask": not (p / xenium.MARKER).is_file()}
    tb = None
    for role in ("sample", "donor"):
        cand = find_columns(names, role)
        dflt = roles.get(role) if roles.get(role) in names else sug.get(role)
        if dflt and dflt not in cand:
            cand = [dflt] + cand
        n = {}
        for c in cand[:8]:
            try:
                n[c] = int(pc.count_distinct(pq.read_table(f, columns=[c])[c]).as_py())
            except Exception:
                n[c] = None
        out[role] = {"default": dflt, "candidates": cand, "n_unique": n}
    return out


@app.post("/api/open", response_class=ORJSONResponse)
def open_(r: OpenReq):
    roi = None
    if r.roi_object.strip():
        roi = Path(r.roi_object.strip()).expanduser()
        if not roi.is_file():
            raise HTTPException(400, f"region_heterogeneity object not found: {roi}")
        if roi.suffix.lower() not in (".rds", ".qs2"):
            raise HTTPException(400, "the region_heterogeneity object must be an .rds file (or .qs2, read through R)")
    ifp = None
    if r.if_image.strip():
        ifp = Path(r.if_image.strip()).expanduser()
        if not ifp.exists():
            raise HTTPException(400, f"IF image not found: {ifp}")
    cols = {k: v.strip() for k, v in (("sample", r.sample_col), ("donor", r.donor_col)) if v is not None and (k == "donor" or v.strip())}
    imp = None
    if r.imc_file.strip():
        imp = Path(r.imc_file.strip()).expanduser()
        if not imp.is_file() or imp.suffix.lower() != ".mcd":
            raise HTTPException(400, f"IMC file not found or not an .mcd file: {imp}")
    open_project(r.path, r.force, roi_object=roi, if_path=ifp, cols=cols or None, he_reg=_he_reg_args(r.he_reg_folder, r.he_reg_index), imc_path=imp)
    return {"ok": True}


# ---------------------------------------------------------------- Xenium output -> project
class XState:
    state = "idle"          # idle | scanning | ready | error
    log: list = []
    error = None
    sess = None             # xenium.Session
    he = None               # he_slide.inspect() of the H&E given at scan time


XST = XState()


def _xlog(msg):
    XST.log.append(msg)
    print("[webbscope xenium]", msg, flush=True)


class XScanReq(BaseModel):
    path: str = ""
    files: dict = {}
    transcripts: bool = False
    ifimg: bool = False
    he_czi: str = ""
    he_ann: str = ""


@app.post("/api/xenium/scan", response_class=ORJSONResponse)
def xenium_scan(r: XScanReq):
    if XST.state == "scanning" or ST.state == "importing":
        raise HTTPException(409, "another import is running")
    XST.state, XST.log, XST.error, XST.sess, XST.he = "scanning", [], None, None, None

    def work():
        try:
            if r.he_czi.strip():
                _xlog("reading the H&E header")
                XST.he = he_slide.inspect(Path(r.he_czi.strip()).expanduser(), Path(r.he_ann.strip()).expanduser() if r.he_ann.strip() else None)
                XST.he["paths"] = {"czi": str(Path(r.he_czi.strip()).expanduser()), "annotations": str(Path(r.he_ann.strip()).expanduser()) if r.he_ann.strip() else None}
                _xlog(f"H&E: {len(XST.he['czi']['scenes'])} scene(s)" + (f", {XST.he['annotations']['n_regions']} annotation polygons" if XST.he["annotations"] else ""))
            src = xenium.resolve_sources(r.path, r.files, _xlog, r.transcripts, r.ifimg)
            XST.sess = xenium.Session(src, _xlog)
            XST.state = "ready"
        except importer.ImportError_ as e:
            XST.state, XST.error = "error", str(e); _xlog("error: " + str(e))
        except Exception as e:
            traceback.print_exc()
            XST.state, XST.error = "error", f"unexpected error: {e}"
    threading.Thread(target=work, daemon=True).start()
    return {"ok": True}


@app.get("/api/xenium/status", response_class=ORJSONResponse)
def xenium_status():
    return {"state": XST.state, "log": XST.log[-40:], "error": XST.error}


def _xs():
    if XST.sess is None:
        raise HTTPException(409, "no Xenium output has been scanned")
    return XST.sess


@app.get("/api/xenium/info", response_class=ORJSONResponse)
def xenium_info():
    s = _xs()
    m = s.meta
    name = f"{m.get('run_name', 'Xenium run')} / {m.get('region_name', '')}".strip(" /")
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "xenium_project"
    return {"name": name, "n_cells": s.n, "n_genes": s.n_genes, "feature_types": s.feature_types,
            "bounds": [float(s.x.min()), float(s.y.min()), float(s.x.max()), float(s.y.max())],
            "has_clusters": s.clusters is not None, "n_clusters": int(s.clusters.nunique()) if s.clusters is not None else 0,
            "has_umap": "umap" in s.emb, "has_pca": "pca" in s.emb,
            "panel": {k: m.get(k) for k in ("panel_name", "panel_organism", "panel_tissue_type", "preservation_method", "slide_id")},
            "zero_count_cells": int((s.gene_counts == 0).sum()), "has_transcripts": "tx" in s.src, "he": XST.he,
            "default_output": str(Path.home() / ".webbscope" / "projects" / safe)}


@app.get("/api/xenium/xy")
def xenium_xy():
    return raw(_xs().xy32)


class XSpec(BaseModel):
    gap: float = 200
    min_cells: int = 500
    min_counts: int = 1
    comps: dict = {}
    rects: list = []


@app.post("/api/xenium/preview", response_class=ORJSONResponse)
def xenium_preview(sp: XSpec):
    s = _xs()
    if not (0 <= sp.gap <= 5000 and 1 <= sp.min_cells <= 10**7 and 0 <= sp.min_counts <= 10**6 and len(sp.rects) <= 200):
        raise HTTPException(400, "parameter out of range")
    try:
        names, idx, comps, counts = s.labels(sp.model_dump())
    except importer.ImportError_ as e:
        raise HTTPException(400, str(e))
    return {"names": names, "labels": idx.tolist(), "comps": comps, "counts": counts,
            "n_kept": int((idx >= 0).sum()), "n_below_min": int((idx == -2).sum())}


class XAnnReq(BaseModel):
    path: str


@app.post("/api/xenium/annotation", response_class=ORJSONResponse)
def xenium_annotation(r: XAnnReq):
    try:
        return _xs().inspect_annotation(r.path)
    except importer.ImportError_ as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(400, f"cannot read the annotation file: {e}")


class XConvertReq(BaseModel):
    spec: XSpec
    output: str
    name: str = ""
    annotation: dict = {}
    transcripts: bool = False
    he: bool = False
    he_reg_folder: str = ""
    he_reg_index: str = ""
    imc_file: str = ""              # IMC slide (.mcd), optional: converted into the cache after the project opens


@app.post("/api/xenium/convert", response_class=ORJSONResponse)
def xenium_convert(r: XConvertReq):
    s = _xs()
    if ST.state == "importing":
        raise HTTPException(409, "an import is already running")
    out = Path(r.output).expanduser()
    if not r.output.strip():
        raise HTTPException(400, "choose an output folder")

    hreg = _he_reg_args(r.he_reg_folder, r.he_reg_index)
    imp = None
    if r.imc_file.strip():
        imp = Path(r.imc_file.strip()).expanduser()
        if not imp.is_file() or imp.suffix.lower() != ".mcd":
            raise HTTPException(400, f"IMC file not found or not an .mcd file: {imp}")

    def prepare():
        _log("converting the Xenium output to a project folder")
        rep = xenium.convert(s, r.spec.model_dump(), out, r.name, r.annotation or None, _log, r.transcripts, XST.he["paths"] if (r.he and XST.he) else None)
        if hreg:
            try:
                he_registered.run(hreg["folder"], hreg["index"], out, _log)
            except (importer.ImportError_, OSError, ValueError) as e:
                _log(f"registered H&E ignored: {e}")
        _log(f"project written to {out}")
        return out
    open_project(str(out), False, prepare, imc_path=imp)
    return {"ok": True}


class BScanReq(BaseModel):
    path: str


@app.post("/api/xenium/batch/scan", response_class=ORJSONResponse)
def xenium_batch_scan(r: BScanReq):
    """Xenium runs inside a parent folder (extracted folders, else .tar.gz / .zip bundles)."""
    try:
        d = batch.discover(r.path)
    except importer.ImportError_ as e:
        raise HTTPException(400, str(e))
    d["default_output"] = str(Path.home() / ".webbscope" / "projects" / batch.safe_name(Path(d["parent"]).name))
    return d


@app.get("/api/xenium/batch/merged", response_class=ORJSONResponse)
def xenium_batch_merged(folder: str):
    """The merged project written to <folder>/_merged by an earlier batch import (null if there is none)."""
    return batch.existing_merged(folder)


class BHeReq(BaseModel):
    folder: str
    runs: list[str]                  # run names to match


@app.post("/api/xenium/batch/he", response_class=ORJSONResponse)
def xenium_batch_he(r: BHeReq):
    """H&E scans (.czi, .ome.tif, .ome.tiff) and annotation files in a folder, with a suggested scan for each run name (matched by the numbers in the file names)."""
    try:
        d = batch.he_discover(r.folder)
    except importer.ImportError_ as e:
        raise HTTPException(400, str(e))
    d["match"] = batch.he_match(r.runs, d["czi"], d["annotations"])
    return d


class BRunReq(BaseModel):
    runs: list                      # [{name, path, kind}] as returned by the scan, minus the ones the user switched off
    output: str
    merge: bool = True
    transcripts: bool = False
    tx_table: bool = False          # write transcript_tissue.parquet per run
    workers: int = 0                # runs converted at once; 0 = automatic
    ifimg: bool = False             # record each run's morphology (IF) image
    he_reg_folder: str = ""         # registered H&E OME-TIFFs for the merged project's samples (<run>_<tissue>), optional
    he_reg_index: str = ""
    gap: float = 200
    min_cells: int = 500
    min_counts: int = 1
    he: dict = {}                   # {run name: {czi, annotations}} from the H&E folder step


@app.post("/api/xenium/batch/run", response_class=ORJSONResponse)
def xenium_batch_run(r: BRunReq):
    if ST.state == "importing" or XST.state == "scanning":
        raise HTTPException(409, "another import is running")
    if not r.output.strip():
        raise HTTPException(400, "choose an output folder")
    if not r.runs or len(r.runs) > 500 or not all(isinstance(x, dict) and x.get("name") and x.get("path") for x in r.runs):
        raise HTTPException(400, "no runs selected")
    if not (0 <= r.workers <= 256 and 0 <= r.gap <= 5000 and 1 <= r.min_cells <= 10**7 and 0 <= r.min_counts <= 10**6):
        raise HTTPException(400, "parameter out of range")
    out = Path(r.output).expanduser()
    spec = {"gap": r.gap, "min_cells": r.min_cells, "min_counts": r.min_counts}
    runs = [{"name": str(x["name"]), "path": str(x["path"]), "kind": str(x.get("kind", ""))} for x in r.runs]
    he = {str(k): {"czi": str(Path(v.get("czi") or "").expanduser()), "annotations": str(Path(v["annotations"]).expanduser()) if v.get("annotations") else ""} for k, v in (r.he or {}).items() if isinstance(v, dict) and v.get("czi")}
    for v in he.values():
        if not Path(v["czi"]).expanduser().is_file() or (v["annotations"] and not Path(v["annotations"]).expanduser().is_file()):
            raise HTTPException(400, f"H&E file not found: {v['czi'] if not Path(v['czi']).expanduser().is_file() else v['annotations']}")

    hreg = _he_reg_args(r.he_reg_folder, r.he_reg_index)

    def prepare():
        rep = batch.run_batch(runs, out, _log, spec, r.merge, r.transcripts, he, r.tx_table, r.workers, r.ifimg, hreg)
        keep = rep["merged"] or rep["first"]
        if keep is None:
            raise importer.ImportError_("no run could be converted; see the log above")
        return keep
    open_project(str(out), False, prepare)
    return {"ok": True}


def _osascript_pick(kind, exts):
    """Native file / folder dialog through macOS' osascript (for when the page is opened in a normal browser, which cannot give a local path).
    Returns the chosen path, "" when cancelled."""
    import subprocess, sys
    if sys.platform != "darwin":
        raise HTTPException(501, "no native file dialog available here: type the path instead")
    if kind == "dir":
        script = 'POSIX path of (choose folder with prompt "Choose a folder")'
    else:
        ty = (" of type {" + ", ".join('"%s"' % re.sub(r"[^A-Za-z0-9]", "", e) for e in exts) + "}") if exts else ""
        script = f'POSIX path of (choose file with prompt "Choose a file"{ty})'
    r = subprocess.run(["osascript", "-e", "tell current application to activate", "-e", script], capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        return ""                                          # cancelled (-128) or unavailable
    return r.stdout.strip().rstrip("/") if kind == "dir" else r.stdout.strip()


@app.get("/api/pick", response_class=ORJSONResponse)
def pick(kind: str = "file", types: str = ""):
    if kind not in ("file", "dir"):
        raise HTTPException(400, "kind must be file or dir")
    return {"path": _osascript_pick(kind, [t for t in types.split(",") if t][:12])}


def _qc_marker():
    f = Path(PROJ) / xenium.MARKER if PROJ else None
    try:
        return json.loads(f.read_text()) if f is not None and f.exists() else None
    except Exception:
        return None


@app.get("/qc-summary")
def qc_summary_page():
    """All-tissue QC summary: one strip plot per metric with a 1st-percentile line."""
    if S is None:
        return Response("no project is open", status_code=409)
    return Response(qc.cohort_html(S, TX, S.meta.get("name", "") if isinstance(getattr(S, "meta", None), dict) else ""), media_type="text/html")


@app.get("/qc/{sample}")
def qc_page(sample: str):
    """One-page QC report of a tissue (HTML; print it to PDF from the browser)."""
    if S is None:
        return Response("no project is open", status_code=409)
    sl(sample)
    return Response(qc.report_html(S, sample, TX, marker=_qc_marker()), media_type="text/html")


@app.get("/api/qc/{sample}", response_class=ORJSONResponse)
def qc_metrics(sample: str):
    sl(sample)
    return qc.metrics(S, sample, TX)


class QCX:
    state, done, total, error, folder = "idle", 0, 0, None, None


class QCExportReq(BaseModel):
    folder: str
    samples: list = []                  # empty = every tissue


@app.get("/api/qc-folder", response_class=ORJSONResponse)
def qc_default_folder():
    """Default folder for the exported QC reports: <cache folder>/qc_reports (the WEBBSCOPE_CACHE the app was started with)."""
    return {"folder": str((importer.cache_root() / "qc_reports").expanduser().resolve())}


@app.post("/api/qc/export", response_class=ORJSONResponse)
def qc_export(r: QCExportReq):
    """Write <folder>/<tissue>.html for every (or the listed) tissue plus qc_summary.csv, in the background."""
    if QCX.state == "running":
        raise HTTPException(409, "an export is already running")
    if not r.folder.strip():
        raise HTTPException(400, "give a folder")
    names = r.samples or list(S.samples)
    bad = [n for n in names if n not in S.slices]
    if bad:
        raise HTTPException(404, f"unknown tissues: {bad[:3]}")
    out = Path(r.folder).expanduser()
    QCX.state, QCX.done, QCX.total, QCX.error, QCX.folder = "running", 0, len(names), None, str(out)
    snap = (S, TX, _qc_marker())

    def work():
        st, tx, marker = snap
        try:
            out.mkdir(parents=True, exist_ok=True)
            (out / "qc_summary.csv").write_text(qc.summary_csv(st, tx))
            (out / "qc_summary.html").write_text(qc.cohort_html(st, tx), encoding="utf-8")
            for n in names:
                safe = re.sub(r"[^A-Za-z0-9._-]+", "_", n).strip("_") or "tissue"
                (out / f"{safe}.html").write_text(qc.report_html(st, n, tx, marker=marker, nav=False), encoding="utf-8")
                QCX.done += 1
            QCX.state = "done"
        except Exception as e:
            traceback.print_exc()
            QCX.state, QCX.error = "error", str(e)
    threading.Thread(target=work, daemon=True).start()
    return {"ok": True, "total": len(names)}


@app.get("/api/qc/export/status", response_class=ORJSONResponse)
def qc_export_status():
    return {"state": QCX.state, "done": QCX.done, "total": QCX.total, "error": QCX.error, "folder": QCX.folder}


@app.get("/api/profile", response_class=ORJSONResponse)
def profile():
    """Time and memory of the last project open (log lines with timestamps, named steps) and per-endpoint request timings since then."""
    return PROF.report()


@app.post("/api/profile/reset", response_class=ORJSONResponse)
def profile_reset():
    """Clear the request timings (start a fresh measurement, e.g. before trying one window)."""
    with PROF.lock:
        PROF.requests = {}
    return {"ok": True}


@app.get("/profile")
def profile_page():
    return FileResponse(WEB / "profile.html")


@app.get("/welcome")
def welcome():
    return FileResponse(WEB / "welcome.html")


@app.get("/xenium-batch")
def xenium_batch_page():
    return FileResponse(WEB / "batch.html")


@app.get("/xenium")
def xenium_page():
    return FileResponse(WEB / "xenium.html")


@app.get("/")
def index():
    return FileResponse(WEB / ("index.html" if S is not None else "welcome.html"))


app.mount("/static", StaticFiles(directory=WEB), name="static")
