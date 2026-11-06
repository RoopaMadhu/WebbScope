"""H&E slide images: Zeiss CZI or OME-TIFF -> tiled pyramids (the viewer's he/ format) and ImageScope annotations -> GeoJSON in micrometres.

A slide scan (CZI) holds one or more scenes (here: the tissue sections of one slide). Each scene becomes its own pyramid in the project's
he/<key>/<level>/<x>_<y>.jpg (level 0 = full resolution, every further level halves it; 512 px tiles, the format the viewer already reads),
described by he_slide.json. The slide frame is the pixel grid of the scan's bounding box (origin = its top-left corner, y down), in micrometres.
Annotations (ImageScope .annotations XML: named polygons, some flagged negative) are expressed in that same frame and written to
he_slide_annotations.geojson with the properties the region importer already uses (annotation_name, region_id, negative_roa, line_color_hex).

An OME-TIFF scan (.ome.tif / .ome.tiff; one whole-slide RGB image, one scene) is read the same way: tiled or in strips, 8-bit RGB, uncompressed, deflate
(no predictor) or JPEG; the pixel size comes from the OME header, the frame is its own pixel grid, and annotations are taken to be in its pixel coordinates.
Only the full-resolution image is read (the file's own reduced levels are not used).

This does not register the H&E to the Xenium coordinates; it only reads and shows the image. Needs the pylibCZIrw (CZI) and Pillow packages (reads only 8-bit RGB
brightfield scans). Command line:  python -m webbscope.he_slide scan.czi|scan.ome.tif [--annotations scan.annotations] --project PROJECT_FOLDER
"""
import argparse, json, math, re, shutil, threading, time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from .importer import ImportError_

TSIZE = 512
BLOCK = 4096                       # output pixels read from the CZI at once (8 x 8 tiles)
INDEX, ANN = "he_slide.json", "he_slide_annotations.geojson"


def _czi_module():
    try:
        from pylibCZIrw import czi
    except ImportError:
        raise ImportError_("Reading .czi files needs the pylibCZIrw package (pip install pylibCZIrw).")
    return czi


def _safe(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "scene"


def czi_info(path):
    """Scenes (rectangles in the scan's pixel coordinates), pixel size (um) and a few facts from the metadata."""
    path = Path(path)
    if not path.is_file():
        raise ImportError_(f"not found: {path}")
    czi = _czi_module()
    try:
        with czi.open_czi(str(path)) as d:
            bb = d.total_bounding_rectangle
            scenes = {int(k): v for k, v in (d.scenes_bounding_rectangle or {}).items()}
            ptypes = {int(k): v for k, v in d.pixel_types.items()}
            meta = d.raw_metadata
    except Exception as e:
        raise ImportError_(f"cannot read {path.name} as a CZI file: {e}")
    if not scenes:
        scenes = {0: bb}
    if set(ptypes.values()) - {"Bgr24"}:
        raise ImportError_(f"only 8-bit RGB (Bgr24) brightfield scans are supported; this file has {sorted(set(ptypes.values()))}")
    r = ET.fromstring(meta)
    px = None
    for e in r.iterfind(".//Scaling/Items/Distance"):
        if e.get("Id") == "X" and e.findtext("Value"):
            px = float(e.findtext("Value")) * 1e6
    if not px or not 0.01 < px < 20:
        raise ImportError_("the CZI metadata has no usable pixel size")
    names = {int(e.get("Index")): e.get("Name") for e in r.iterfind(".//Dimensions/S/Scenes/Scene") if e.get("Index") is not None}
    return {"file": path.name, "px_um": px, "bbox": [bb.x, bb.y, bb.w, bb.h],
            "scenes": [{"index": k, "name": names.get(k) or f"scene{k}", "x": v.x, "y": v.y, "w": v.w, "h": v.h} for k, v in sorted(scenes.items())],
            "objective": r.findtext(".//ObjectiveName"), "magnification": r.findtext(".//NominalMagnification"),
            "acquired": r.findtext(".//AcquisitionDateAndTime")}


OME_EXT = (".ome.tif", ".ome.tiff")
_UNIT_UM = {"µm": 1.0, "um": 1.0, "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6}


def is_ome(path):
    return str(path).lower().endswith(OME_EXT)


def _ome_decode(raw, hd):
    """One tile (or strip) of an RGB OME-TIFF -> (th, tw, 3) uint8."""
    import io, zlib
    from PIL import Image
    th, tw, comp = hd["th"], hd["tw"], hd["comp"]
    if comp == 1:
        a = np.frombuffer(raw, np.uint8)
    elif comp in (8, 32946):
        a = np.frombuffer(zlib.decompress(raw), np.uint8)
    elif comp == 7:                                                    # JPEG: the tables shared by all tiles (tag 347) go in front of each tile's data
        tb = hd.get("jpeg_tables")
        data = (tb[:-2] + raw[2:]) if tb else raw
        im = Image.open(io.BytesIO(data)).convert("RGB")
        a = np.asarray(im, np.uint8)
        if a.shape[0] < th or a.shape[1] < tw:                         # a strip / edge tile may be stored smaller than the nominal size
            p = np.zeros((th, tw, 3), np.uint8)
            p[:a.shape[0], :a.shape[1]] = a[:th, :tw]
            return p
        return a[:th, :tw]
    else:
        raise ImportError_(f"TIFF compression {comp} is not supported")
    if a.size < th * tw * 3:                                           # the last strip is stored short
        a = np.concatenate([a, np.zeros(th * tw * 3 - a.size, np.uint8)])
    return a[:th * tw * 3].reshape(th, tw, 3)


def ome_info(path):
    """Header of an OME-TIFF scan, in the shape czi_info returns (one scene = the whole image), plus what the tile reader needs (key "hd")."""
    from . import if_image
    path = Path(path)
    if not path.is_file():
        raise ImportError_(f"not found: {path}")
    try:
        t = if_image._Tiff(path)
    except Exception as e:
        raise ImportError_(f"cannot read {path.name} as a TIFF file: {e}")
    try:
        tags, _ = t.ifd(t.first)
        g = lambda k, d=None: t.val(tags[k]) if k in tags else d
        w, h = g(256)[0], g(257)[0]
        comp, spp, bits = g(259, (1,))[0], g(277, (1,))[0], tuple(g(258, (8,)))
        if spp != 3 or set(bits) != {8}:
            raise ImportError_(f"{path.name}: only 8-bit RGB H&E scans are supported (found {spp} sample(s) per pixel, {bits} bits)")
        if comp not in (1, 7, 8, 32946):
            raise ImportError_(f"{path.name}: TIFF compression {comp} is not supported; uncompressed, deflate and JPEG are")
        if g(284, (1,))[0] != 1 or g(317, (1,))[0] != 1:
            raise ImportError_(f"{path.name}: planar-separate TIFFs and TIFFs with a predictor are not supported")
        if 322 in tags:
            tw, th, off, cnt = g(322)[0], g(323)[0], list(g(324)), list(g(325))
        elif 273 in tags:                                              # strips: a strip is a tile as wide as the image
            tw, th, off, cnt = w, min(g(278, (h,))[0], h), list(g(273)), list(g(279))
        else:
            raise ImportError_(f"{path.name}: the image has no tiles or strips")
        desc = g(270, "") or ""
        jt = bytes(g(347)) if 347 in tags else None                     # JPEGTables: shared by all JPEG tiles
    finally:
        t.close()
    m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', desc)
    u = re.search(r'PhysicalSizeXUnit="([^"]*)"', desc)
    px = float(m.group(1)) * _UNIT_UM.get(u.group(1) if u else "µm", 1.0) if m else None
    if not px or not 0.01 < px < 20:
        raise ImportError_(f"{path.name}: the OME header has no usable pixel size (PhysicalSizeX {m.group(1) if m else 'missing'})")
    nm = re.search(r'<Image[^>]*?Name="([^"]*)"', desc)
    hd = {"w": w, "h": h, "tw": tw, "th": th, "off": off, "cnt": cnt, "comp": comp, "jpeg_tables": jt, "decode": _ome_decode, "fill": (255, 255, 255), "black_is_empty": False}
    return {"file": path.name, "px_um": px, "bbox": [0, 0, w, h], "scenes": [{"index": 0, "name": (nm.group(1) if nm and nm.group(1) else "scan"), "x": 0, "y": 0, "w": w, "h": h}],
            "objective": None, "magnification": None, "acquired": None, "hd": hd}


def scan_info(path):
    """Scenes, pixel size and a few facts of an H&E scan: a CZI or an OME-TIFF."""
    return ome_info(path) if is_ome(path) else czi_info(path)


def _hex(v):
    v = int(v)
    return "#%02x%02x%02x" % (v & 255, (v >> 8) & 255, (v >> 16) & 255)          # ImageScope LineColor = R + 256 G + 65536 B


def read_annotations(path):
    """ImageScope XML -> [{name, color, visible, regions: [{negative, pts (n, 2) in scan pixels}]}]; non-polygon regions are counted, not read."""
    path = Path(path)
    if not path.is_file():
        raise ImportError_(f"not found: {path}")
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        raise ImportError_(f"{path.name} is not valid XML: {e}")
    if root.tag != "Annotations":
        raise ImportError_(f"{path.name} is not an ImageScope annotations file (root element <{root.tag}>)")
    out, skipped = [], 0
    for a in root.iterfind("Annotation"):
        regs = []
        for reg in a.iterfind("./Regions/Region"):
            if reg.get("Type") != "Polygon":
                skipped += 1
                continue
            v = np.array([[float(q.get("X")), float(q.get("Y"))] for q in reg.iterfind("./Vertices/V")], float)
            if len(v) < 3:
                skipped += 1
                continue
            regs.append({"negative": reg.get("NegativeROA") == "1", "pts": v})
        out.append({"name": (a.get("Name") or "").strip() or "unnamed", "color": int(a.get("LineColor") or 0), "visible": a.get("Visible") != "False", "regions": regs})
    return out, skipped


def annotations_geojson(anns, info):
    """Polygons in micrometres of the slide frame (origin = top-left of the scan's bounding box, y down) + a per-annotation summary."""
    bx, by = info["bbox"][0], info["bbox"][1]
    px = info["px_um"]
    rects = [(s["x"] - bx, s["y"] - by, s["x"] - bx + s["w"], s["y"] - by + s["h"], s["index"]) for s in info["scenes"]]
    feats, summ, n_in, n_all = [], [], 0, 0
    for a in anns:
        k = 0
        for r in a["regions"]:
            p = r["pts"]
            ins = np.zeros(len(p), bool)
            for x0, y0, x1, y1, _ in rects:
                ins |= (p[:, 0] >= x0) & (p[:, 0] <= x1) & (p[:, 1] >= y0) & (p[:, 1] <= y1)
            n_in += int(ins.sum()); n_all += len(p)
            c = p.mean(0)
            sc = next((i for x0, y0, x1, y1, i in rects if x0 <= c[0] <= x1 and y0 <= c[1] <= y1), None)
            ring = (p * px).round(3).tolist()
            ring.append(ring[0])
            feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
                          "properties": {"annotation_name": a["name"], "region_id": f"{a['name']}_{k}", "negative_roa": r["negative"],
                                         "line_color_hex": _hex(a["color"]), "scene": sc}})
            k += 1
        summ.append({"name": a["name"], "color": _hex(a["color"]), "n_regions": len(a["regions"]), "n_negative": sum(r["negative"] for r in a["regions"])})
    return {"type": "FeatureCollection", "features": feats}, summ, {"vertices": n_all, "vertices_inside_scenes": n_in}


def inspect(czi_path, ann_path=None):
    """Quick summary for the import page (header + annotation file only; no pixels are read)."""
    info = scan_info(czi_path)
    info = {k: v for k, v in info.items() if k != "hd"}                       # the tile-reader header is not for the page
    out = {"czi": info, "annotations": None}
    if ann_path:
        anns, skipped = read_annotations(ann_path)
        _, summ, chk = annotations_geojson(anns, info)
        out["annotations"] = {"file": Path(ann_path).name, "names": summ, "n_regions": sum(s["n_regions"] for s in summ), "skipped_regions": skipped, **chk}
    return out


# ---------------------------------------------------------------- pyramid
_tl = threading.local()


def _reader(path):
    if getattr(_tl, "path", None) != path:
        _tl.cm = _czi_module().open_czi(path)
        _tl.d = _tl.cm.__enter__()
        _tl.path = path
    return _tl.d


def _levels(w, h):
    return max(1, math.ceil(math.log2(max(w, h) / TSIZE)) + 1) if max(w, h) > TSIZE else 1


def _job(a):
    from PIL import Image                                                         # imported here so the server starts even without Pillow
    path, out, key, sx, sy, z, tx0, ty0, ow, oh, quality = a
    s = 2 ** z
    arr = _reader(path).read(roi=(sx + tx0 * TSIZE * s, sy + ty0 * TSIZE * s, ow * s, oh * s), zoom=1.0 / s, plane={"C": 0},
                                                  background_pixel=(1.0, 1.0, 1.0))   # unscanned parts of the mosaic: white at every level (the reader's default is black)
    arr = arr[:, :, ::-1]                                                          # BGR -> RGB
    n = 0
    for j in range(0, arr.shape[0], TSIZE):
        for i in range(0, arr.shape[1], TSIZE):
            t = np.ascontiguousarray(arr[j:j + TSIZE, i:i + TSIZE])
            if not t.size:
                continue
            d = Path(out) / "he" / key / str(z)
            d.mkdir(parents=True, exist_ok=True)
            Image.fromarray(t).save(d / f"{tx0 + i // TSIZE}_{ty0 + j // TSIZE}.jpg", quality=quality)
            n += 1
    return n


def build(czi_path, out_dir, ann_path=None, log=print, quality=88, workers=4):
    """Write he/<key>/... pyramids, he_slide.json and (with ann_path) he_slide_annotations.geojson into the project folder out_dir."""
    t0 = time.time()
    out_dir = Path(out_dir)
    info = scan_info(czi_path)
    px, bx, by = info["px_um"], info["bbox"][0], info["bbox"][1]
    if is_ome(czi_path):
        return _build_ome(czi_path, out_dir, info, ann_path, log, workers, t0)
    log(f"H&E: {info['file']}, {len(info['scenes'])} scene(s), pixel size {px:.4f} um ({info['magnification'] or '?'}x)")
    scenes, jobs = [], []
    for s in info["scenes"]:
        key = _safe(f"he_{Path(czi_path).stem}__{s['name']}")
        shutil.rmtree(out_dir / "he" / key, ignore_errors=True)
        nz = _levels(s["w"], s["h"])
        scenes.append({"key": key, "name": s["name"], "scene": s["index"], "w": s["w"], "h": s["h"], "px": px, "x0": round((s["x"] - bx) * px, 3),
                       "y0": round((s["y"] - by) * px, 3), "tsize": TSIZE, "nz": nz})
        for z in range(nz):
            lw, lh = max(1, s["w"] >> z), max(1, s["h"] >> z)
            for ty0 in range(0, math.ceil(lh / TSIZE), BLOCK // TSIZE):
                for tx0 in range(0, math.ceil(lw / TSIZE), BLOCK // TSIZE):
                    jobs.append((str(czi_path), str(out_dir), key, s["x"], s["y"], z, tx0, ty0,
                                 min(BLOCK, lw - tx0 * TSIZE), min(BLOCK, lh - ty0 * TSIZE), quality))
    log(f"H&E: writing {len(jobs)} blocks of tiles ({sum(sc['nz'] for sc in scenes)} pyramid levels in all) with {workers} threads")
    (out_dir / "he").mkdir(parents=True, exist_ok=True)
    done = tiles = 0
    last = time.time()
    with ThreadPoolExecutor(workers) as ex:
        for n in ex.map(_job, jobs):
            done += 1; tiles += n
            if time.time() - last > 5:
                log(f"  H&E tiles: {done}/{len(jobs)} blocks, {tiles:,} tiles, {time.time() - t0:.0f} s")
                last = time.time()
    index = {"source": info["file"], "px": px, "frame": "micrometres from the top-left of the scan's bounding box, y down (not registered to Xenium)",
             "bounds_um": [0, 0, round(info["bbox"][2] * px, 1), round(info["bbox"][3] * px, 1)], "scenes": scenes,
             "objective": info["objective"], "magnification": info["magnification"], "acquired": info["acquired"], "annotations": None}
    if ann_path:
        anns, skipped = read_annotations(ann_path)
        gj, summ, chk = annotations_geojson(anns, info)
        (out_dir / ANN).write_text(json.dumps(gj, separators=(",", ":")))
        index["annotations"] = {"file": Path(ann_path).name, "names": summ, "skipped_regions": skipped, **chk}
        log(f"H&E annotations: {sum(s['n_regions'] for s in summ)} polygons in {len(summ)} annotation names ({chk['vertices_inside_scenes']:,} of {chk['vertices']:,} vertices fall inside the scanned scenes)")
    elif (out_dir / ANN).exists():
        (out_dir / ANN).unlink()
    (out_dir / INDEX).write_text(json.dumps(index, indent=1))
    log(f"H&E: {tiles:,} tiles written in {time.time() - t0:.0f} s")
    return {"tiles": tiles, "scenes": len(scenes), "seconds": round(time.time() - t0, 1)}


def _build_ome(path, out_dir, info, ann_path, log, workers, t0):
    """The OME-TIFF scan as one scene: level 0 is read tile row by tile row, the further levels are halved from it (he_registered.build_one)."""
    from . import he_registered
    px, s, hd = info["px_um"], info["scenes"][0], info["hd"]
    log(f"H&E: {info['file']} (OME-TIFF), {s['w']:,} x {s['h']:,} px, pixel size {px:.4f} um")
    key = _safe(f"he_{Path(path).name.split('.')[0]}__{s['name']}")
    (out_dir / "he").mkdir(parents=True, exist_ok=True)
    nz, tiles = he_registered.build_one(str(path), hd, out_dir / "he" / key, lambda m: log("  " + m.strip()), workers)
    scene = {"key": key, "name": s["name"], "scene": 0, "w": s["w"], "h": s["h"], "px": px, "x0": 0.0, "y0": 0.0, "tsize": TSIZE, "nz": nz}
    index = {"source": info["file"], "px": px, "frame": "micrometres from the top-left of the scan, y down (not registered to Xenium)",
             "bounds_um": [0, 0, round(s["w"] * px, 1), round(s["h"] * px, 1)], "scenes": [scene],
             "objective": None, "magnification": None, "acquired": None, "annotations": None}
    if ann_path:
        anns, skipped = read_annotations(ann_path)
        gj, summ, chk = annotations_geojson(anns, {k: v for k, v in info.items() if k != "hd"})
        (out_dir / ANN).write_text(json.dumps(gj, separators=(",", ":")))
        index["annotations"] = {"file": Path(ann_path).name, "names": summ, "skipped_regions": skipped, **chk}
        log(f"H&E annotations: {sum(x['n_regions'] for x in summ)} polygons in {len(summ)} annotation names ({chk['vertices_inside_scenes']:,} of {chk['vertices']:,} vertices fall inside the image; "
            "annotations are taken to be in the image's pixel coordinates)")
    elif (out_dir / ANN).exists():
        (out_dir / ANN).unlink()
    (out_dir / INDEX).write_text(json.dumps(index, indent=1))
    log(f"H&E: {tiles:,} tiles written in {time.time() - t0:.0f} s")
    return {"tiles": tiles, "scenes": 1, "seconds": round(time.time() - t0, 1)}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m webbscope.he_slide", description="Add an H&E slide scan (CZI or OME-TIFF) and optional ImageScope annotations to a project folder")
    ap.add_argument("czi", help="the scan: .czi, .ome.tif or .ome.tiff"); ap.add_argument("--annotations"); ap.add_argument("--project", required=True, help="project folder to add to (nothing else in it is changed)")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    try:
        build(a.czi, a.project, a.annotations, workers=a.workers)
    except ImportError_ as e:
        raise SystemExit(f"error: {e}")


if __name__ == "__main__":
    main()
