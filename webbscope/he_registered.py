"""Registered H&E images: OME-TIFF (one per tissue, already warped into the Xenium frame) -> the viewer's he/ pyramids.

Each .ome.tif holds one tissue as a single tiled 8-bit RGB image (deflate) with the pixel size and the position of its top-left corner in the Xenium
frame (micrometres) in the OME header. An index.json maps tissues to images: {sample_id: [{key, tile, w, h, px, x0, y0, tsize, nz}]}, where `tile` is the
image's file name without "_warped_he.ome.tif". This module matches the images in a folder to the samples by that name, tiles them into
he/<sample>__<tile>/<level>/<x>_<y>.jpg (level 0 = full resolution, every further level halves it; 512 px tiles, as in he_slide.py) and merges the entries
into he/index.json of the project, which the viewer reads ("H&E" toggle next to the cells).

The pixel size and position come from the image itself; index.json is only used for the name -> sample mapping, and a difference between the two is reported.
Only samples that are in the project's cell table are written unless --all is given. Needs Pillow only: tiles are read straight from the TIFF and
inflated, so a 2 GB image is never held in memory (one row of 1024 px tiles at a time). Command line:
    python -m webbscope.he_registered IMAGE_FOLDER --index index.json --project PROJECT_FOLDER [--all] [--force] [--workers 4]
"""
import argparse, json, math, re, shutil, time, zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from .importer import ImportError_
from .he_slide import TSIZE, _levels

SUFFIX = "_warped_he.ome.tif"
QUALITY = 88
BACKGROUND = (233, 228, 232)       # the warped images are black where the warp has no data; shown as this slide-background colour (the value in the project's earlier tiles)
_NS = {"o": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}


def _tile_name(p):
    n = Path(p).name
    return re.sub(r"(_warped_he)?(\.ome)?\.tiff?$", "", n, flags=re.I)      # _warped_he.ome.tif / .ome.tiff, .ome.tif, .tif ...


def discover(folder):
    """{tile name: path} of the .ome.tif / .ome.tiff images in `folder` and its sub-folders (several images with one name are an error)."""
    out = {}
    for p in sorted(Path(folder).rglob("*")):
        if p.is_file() and re.search(r"\.ome\.tiff?$", p.name, flags=re.I) and not p.name.startswith("."):
            t = _tile_name(p)
            if t in out:
                raise ImportError_(f"two images are named {t}: {out[t]} and {p}")
            out[t] = p
    return out


def read_header(path):
    """Size, pixel size and top-left position (um) of an OME-TIFF from its first image directory, plus what the tile reader needs. Nothing is decoded."""
    from PIL import Image
    import xml.etree.ElementTree as ET
    Image.MAX_IMAGE_PIXELS = None
    try:
        im = Image.open(path)
    except Exception as e:
        raise ImportError_(f"cannot read {Path(path).name} as a TIFF file: {e}")
    t = im.tag_v2
    if t.get(322) is None or t.get(259) != 8 or t.get(258) != (8, 8, 8) or t.get(284, 1) != 1 or t.get(317, 1) != 1:
        raise ImportError_(f"{Path(path).name}: only tiled, 8-bit RGB, deflate-compressed TIFFs without a predictor are read "
                           f"(found compression {t.get(259)}, bits {t.get(258)}, tiled {t.get(322) is not None}, planar {t.get(284)}, predictor {t.get(317)})")
    try:
        px = ET.fromstring(t[270]).find(".//o:Pixels", _NS)
        pl = px.find("o:Plane", _NS)
        unit = px.get("PhysicalSizeXUnit", "µm")
        if "m" not in unit.lower() or px.get("PhysicalSizeX") is None:
            raise ValueError(f"physical size unit {unit!r}")
        ps, x0, y0 = float(px.get("PhysicalSizeX")), float(pl.get("PositionX")), float(pl.get("PositionY"))
    except Exception as e:
        raise ImportError_(f"{Path(path).name}: the OME header has no usable pixel size / position ({e})")
    return {"w": im.size[0], "h": im.size[1], "px": ps, "x0": x0, "y0": y0, "tw": t[322], "th": t[323], "off": list(t[324]), "cnt": list(t[325])}


def _row(path, hd, r):
    """Pixels of tile row r of the TIFF as a (rows, w, 3) uint8 array; pixels with no data (stored as black, or in tiles that are not stored) are BACKGROUND."""
    tw, th, w, h = hd["tw"], hd["th"], hd["w"], hd["h"]
    nx = math.ceil(w / tw)
    rows = min(th, h - r * th)
    out = np.empty((rows, w, 3), np.uint8)
    out[:] = hd.get("fill", BACKGROUND)
    dec = hd.get("decode")                                      # a scan (he_slide.py) brings its own decoder; registered images are deflate
    with open(path, "rb") as f:
        for c in range(nx):
            i = r * nx + c
            if not hd["cnt"][i]:
                continue
            f.seek(hd["off"][i])
            raw = f.read(hd["cnt"][i])
            a = dec(raw, hd) if dec else np.frombuffer(zlib.decompress(raw), np.uint8).reshape(th, tw, 3)
            cw = min(tw, w - c * tw)
            out[:, c * tw:c * tw + cw] = a[:rows, :cw]
    if hd.get("black_is_empty", True):
        out[~out.any(axis=2)] = BACKGROUND                     # registered images: pure black (0, 0, 0) = no data
    return out


def _save(a, d, x, y):
    from PIL import Image
    d.mkdir(parents=True, exist_ok=True)
    Image.fromarray(np.ascontiguousarray(a)).save(d / f"{x}_{y}.jpg", quality=QUALITY)


def _level0(path, hd, base, r):
    a = _row(path, hd, r)                                      # one row of source tiles = th / 512 rows of output tiles
    n = 0
    for j in range(0, a.shape[0], TSIZE):
        ty = (r * hd["th"] + j) // TSIZE
        for i in range(0, a.shape[1], TSIZE):
            _save(a[j:j + TSIZE, i:i + TSIZE], base / "0", i // TSIZE, ty)
            n += 1
    return n


def _band(path, hd, base, u):
    """Output tile row u (rows u*512 ...) assembled from the source tile rows / strips that overlap it: for sources whose tiles are not a whole number of 512 px high."""
    w, h, th = hd["w"], hd["h"], hd["th"]
    y0, y1 = u * TSIZE, min(h, (u + 1) * TSIZE)
    parts = []
    for r in range(y0 // th, (y1 - 1) // th + 1):
        a = _row(path, hd, r)
        ys = r * th
        parts.append(a[max(y0, ys) - ys:min(y1, ys + a.shape[0]) - ys])
    band = np.concatenate(parts, 0)
    n = 0
    for i in range(0, w, TSIZE):
        _save(band[:, i:i + TSIZE], base / "0", i // TSIZE, u)
        n += 1
    return n


def _up(base, z, w, h, ty):
    """Row ty of level z: each tile is the 2 x 2 tiles below it, halved."""
    from PIL import Image
    lw, lh = max(1, w >> z), max(1, h >> z)
    pw, ph = max(1, w >> (z - 1)), max(1, h >> (z - 1))
    n = 0
    for tx in range(math.ceil(lw / TSIZE)):
        tw, th = min(TSIZE, lw - tx * TSIZE), min(TSIZE, lh - ty * TSIZE)
        canvas = Image.new("RGB", (2 * TSIZE, 2 * TSIZE), (255, 255, 255))
        for dy in range(2):
            for dx in range(2):
                if (2 * tx + dx) * TSIZE < pw and (2 * ty + dy) * TSIZE < ph:
                    with Image.open(base / str(z - 1) / f"{2 * tx + dx}_{2 * ty + dy}.jpg") as c:
                        canvas.paste(c, (dx * TSIZE, dy * TSIZE))
        sw, sh = min(2 * TSIZE, pw - 2 * tx * TSIZE), min(2 * TSIZE, ph - 2 * ty * TSIZE)
        t = canvas.crop((0, 0, sw, sh)).resize((tw, th), Image.BOX)
        (base / str(z)).mkdir(parents=True, exist_ok=True)
        t.save(base / str(z) / f"{tx}_{ty}.jpg", quality=QUALITY)
        n += 1
    return n


def build_one(path, hd, base, log, workers=4):
    """Write the pyramid of one image into the folder `base`; returns the number of tiles."""
    shutil.rmtree(base, ignore_errors=True)
    w, h = hd["w"], hd["h"]
    nz = _levels(w, h)
    whole = hd["th"] % TSIZE == 0                                       # source tile rows are whole multiples of the output tile: one source row at a time
    nrow = math.ceil(h / hd["th"]) if whole else math.ceil(h / TSIZE)
    tiles = 0
    with ThreadPoolExecutor(workers) as ex:
        for k, n in enumerate(ex.map(lambda r: (_level0 if whole else _band)(path, hd, base, r), range(nrow)), 1):
            tiles += n
            if k % 10 == 0:
                log(f"    level 0: {k}/{nrow} tile rows")
        for z in range(1, nz):
            for n in ex.map(lambda ty: _up(base, z, w, h, ty), range(math.ceil(max(1, h >> z) / TSIZE))):
                tiles += n
    return nz, tiles


def project_samples(project):
    """Sample ids in the project's cell table (the column named by the sample role of project.json)."""
    import pyarrow.parquet as pq
    project = Path(project)
    role = json.loads((project / "project.json").read_text()).get("roles", {}).get("sample")
    if not role:
        raise ImportError_("project.json has no sample role")
    return {str(s) for s in pq.read_table(project / "cells.parquet", columns=[role]).column(0).unique().to_pylist()}


def match(images, index):
    """Returns ({sample: [(tile name, path, index entry)]}, tiles named in the index with no image, images the index does not list)."""
    by_tile = {e["tile"]: (s, e) for s, es in index.items() for e in es}
    got, unlisted = {}, []
    for t, p in images.items():
        if t in by_tile:
            got.setdefault(by_tile[t][0], []).append((t, p, by_tile[t][1]))
        else:
            unlisted.append(t)
    return got, sorted(set(by_tile) - set(images)), sorted(unlisted)


def run(folder, index_path, project, log=print, every_sample=False, force=False, workers=4, out_root=None, samples=None):
    """Tile the registered images of `folder` for the samples of `project`. index_path may be empty: an image is then matched to the sample whose id is its file name.
    Tiles and he/index.json go to `out_root` (default: the project folder); `samples` (a set) replaces the sample ids read from the project's cell table."""
    t0 = time.time()
    project, folder = Path(project), Path(folder)
    out_root = Path(out_root) if out_root else project
    images = discover(folder)
    have = None if every_sample else (set(samples) if samples is not None else project_samples(project))
    if index_path:
        index = json.loads(Path(index_path).read_text())
        what = f"listed in {Path(index_path).name}"
    else:                                                      # no index: the image file name (without _warped_he.ome.tif / .ome.tif) is the sample id
        index = {t: [{"tile": t}] for t in images if have is not None and t in have}
        what = "named like a sample"
    got, no_image, unlisted = match(images, index)
    log(f"registered H&E: {len(images)} image(s) in {folder}; {sum(len(v) for v in got.values())} of them are {what}")
    if unlisted:
        log(f"  not in the index, so not matched to any sample ({len(unlisted)}): {', '.join(unlisted)}")
    if no_image:
        log(f"  in the index but no image in the folder ({len(no_image)}): {', '.join(no_image)}")
    out_f = out_root / "he" / "index.json"
    out = json.loads(out_f.read_text()) if out_f.exists() else {}
    todo = {s: v for s, v in got.items() if have is None or s in have}
    if have is not None and got:
        log(f"  samples in the project: {len(have)}; images for them: {sum(len(v) for v in todo.values())}"
            + (f" (images of {len(got) - len(todo)} other sample(s) skipped; --all writes them too)" if len(got) > len(todo) else ""))
    if not index_path and have is not None and not got:
        log(f"  no image is named like a sample of the project (samples look like: {', '.join(sorted(have)[:3])}; images look like: {', '.join(sorted(images)[:3])}); give an index.json to map them")
    written = []
    for s, items in todo.items():
        ents = []
        for t, p, ie in items:
            hd = read_header(p)
            key = f"{s}__{t}"
            ent = {"key": key, "tile": t, "w": hd["w"], "h": hd["h"], "px": hd["px"], "x0": hd["x0"], "y0": hd["y0"], "tsize": TSIZE, "nz": _levels(hd["w"], hd["h"])}
            diff = [f"{k} {ie.get(k)} (index) vs {ent[k]} (image)" for k in ("w", "h", "px", "x0", "y0") if abs(float(ie.get(k, math.nan)) - ent[k]) > 1e-6]
            if diff:
                log(f"  WARNING {t}: index.json and the image disagree: {'; '.join(diff)}. The image's own values are used.")
            base = out_root / "he" / key
            done = base / "done.json"
            if not force and done.exists() and json.loads(done.read_text()) == ent:
                log(f"  {key}: tiles from an earlier run are current, not rebuilt")
            else:
                log(f"  {key}: {hd['w']} x {hd['h']} px at {hd['px']} um/px, {ent['nz']} levels")
                nz, n = build_one(p, hd, base, log, workers)
                done.write_text(json.dumps(ent))
                log(f"    {n:,} tiles")
            ents.append(ent)
            written.append(key)
        out[s] = ents                                          # this sample's entries are replaced as a whole
    if written:
        (out_root / "he").mkdir(parents=True, exist_ok=True)
        out_f.write_text(json.dumps(out))
    log(f"registered H&E: {len(written)} image(s) attached in {time.time() - t0:.0f} s" if written else "registered H&E: nothing to attach")
    return {"written": written, "unlisted": unlisted, "no_image": no_image}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(prog="python -m webbscope.he_registered", description="Add registered H&E OME-TIFFs to a project folder, matched to samples by an index.json")
    ap.add_argument("folder"); ap.add_argument("--index", required=True); ap.add_argument("--project", required=True)
    ap.add_argument("--all", action="store_true", help="also write images of samples that are not in the project's cell table")
    ap.add_argument("--force", action="store_true"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    try:
        run(a.folder, a.index, a.project, every_sample=a.all, force=a.force, workers=a.workers)
    except ImportError_ as e:
        raise SystemExit(f"error: {e}")
