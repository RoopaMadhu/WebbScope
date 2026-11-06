"""Imaging mass cytometry (IMC) slides from a Standard BioTools (Fluidigm) Hyperion .mcd file, shown beside the cells (not aligned to them).

An .mcd file holds acquisitions (ROIs). Each acquisition is a list of pixel records [X, Y, Z, channel 1 ... channel n] of float32, row by row,
1 um per pixel; an acquisition that was stopped early simply has fewer records. The metadata (acquisition names, sizes, positions on the slide,
channel names) is read with `readimc`; the pixel data is read here, in bands of rows, so that a 10 GB acquisition never has to fit in memory.

convert() writes, for every acquisition, a pyramid in the cache (not next to the .mcd, which is never changed):
    <store>/<key>/l<z>.npy   float16 array (channels, y, x): level 0 = the raw counts, level z = the mean of 2 x 2 pixels of level z-1
    <store>/<key>/meta.json  size, number of rows present, levels, default intensity limit per channel
float16 keeps the counts exactly up to 2048 and to 0.05% above (values above 65504, which only an argon or other gas channel reaches, are clipped to 65504); it is a display copy, not for quantification (read the .mcd for that).
Tiles are cut from these arrays on request and returned as JPEG: several channels, each in its own colour, are added up in one tile.
"""
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
from PIL import Image

from . import if_image

TILE = 512
_EXTRA = 3                       # X, Y, Z in front of the channels in every pixel record
_STAMP = 4                       # bump to rebuild stores written by an older version of the converter


class IMCError(Exception):
    pass


def _name_of(label):
    """'169Tm_CD45' -> 'CD45'; a label without a marker part ('127I', '80ArAr') -> None."""
    return label.split("_", 1)[1] if "_" in label else None


def channel_info(label, name, mass):
    nm = _name_of(label)
    if nm is None:
        kind, text = "other", label
    elif re.match(r"DNA", nm, re.I):
        kind, text = "dna", f"DNA ({label.split('_')[0]})"
    elif re.match(r"Cell_Seg", nm, re.I):
        kind, text = "seg", f"{nm.replace('_', ' ')} ({label.split('_')[0]})"
    else:
        kind, text = "marker", f"{nm} ({label.split('_')[0]})"
    return {"label": label, "name": text, "kind": kind, "mass": mass}


def scan(path):
    """Metadata of an .mcd file: {rois: [...], channels: [...], slide_um: [w, h]}. Nothing is converted."""
    try:
        from readimc import MCDFile
    except ImportError as e:                                   # pragma: no cover
        raise IMCError("reading .mcd files needs the 'readimc' package (pip install readimc)") from e
    p = Path(path).expanduser()
    if not p.is_file():
        raise IMCError(f"not a file: {p}")
    size = p.stat().st_size
    rois, chans = [], None
    try:
        with MCDFile(p) as f:
            for s in f.slides:
                for a in s.acquisitions:
                    md = a.metadata
                    try:
                        w, h, nc = int(md["MaxX"]), int(md["MaxY"]), a.num_channels
                        st, en, vb = int(md["DataStartOffset"]), int(md["DataEndOffset"]), int(md["ValueBytes"])
                    except (KeyError, ValueError):
                        continue
                    labels = list(a.channel_labels)
                    cur = [channel_info(lb, nm, ms) for lb, nm, ms in zip(labels, a.channel_names, a.channel_masses)]
                    if chans is None:
                        chans = cur
                    elif [c["label"] for c in cur] != [c["label"] for c in chans]:
                        raise IMCError("the acquisitions of this file have different channel lists; that is not supported yet")
                    n = (en - st) // ((nc + _EXTRA) * vb) if en <= size else 0
                    x0 = min(float(md.get("ROIStartXPosUm", 0)), float(md.get("ROIEndXPosUm", 0)))
                    y0 = min(float(md.get("ROIStartYPosUm", 0)), float(md.get("ROIEndYPosUm", 0)))
                    rois.append({"key": f"a{a.id}", "id": a.id, "name": a.description or f"acquisition {a.id}", "w": w, "h": h, "records": n,
                                 "rows": min(h, n // w) if w else 0, "complete": n == w * h, "start": st, "end": en, "channels": nc,
                                 "x_um": x0, "y_um": y0, "pixel_um": float(getattr(a, "pixel_size_x_um", 1) or 1)})
    except IMCError:
        raise
    except Exception as e:
        raise IMCError(f"{p.name} could not be read as an .mcd file: {e}") from e
    if not rois:
        raise IMCError(f"{p.name} has no acquisitions")
    return {"file": str(p), "size": size, "mtime": p.stat().st_mtime_ns, "rois": rois, "channels": chans}


def store_dir(root, mcd):
    p = Path(mcd).expanduser().resolve()
    k = hashlib.sha1(f"{p}|{p.stat().st_size}|{p.stat().st_mtime_ns}|{_STAMP}".encode()).hexdigest()[:12]
    return Path(root) / "imc" / f"{p.stem}-{k}"


def _f16(a):
    return np.clip(np.nan_to_num(a, nan=0.0, posinf=65504.0), 0, 65504).astype(np.float16)


def _levels(w, h):
    n = 1
    while max(w, h) > TILE * 2 ** (n - 1):
        n += 1
    return n


def _convert_roi(mcd, r, d, nch, log, progress):
    d.mkdir(parents=True, exist_ok=True)
    w, h, n, rec = r["w"], r["h"], r["records"], nch + _EXTRA
    mm = np.memmap(mcd, dtype=np.float32, mode="r", offset=r["start"], shape=(n, rec))
    l0 = np.lib.format.open_memmap(d / "l0.npy", mode="w+", dtype=np.float16, shape=(nch, h, w))
    band = max(1, int(3e7 // (w * rec)))                      # rows per band: about 120 MB of the file at a time
    done_rows = 0
    for y0 in range(0, h, band):
        i0, i1 = y0 * w, min(n, (y0 + band) * w)
        if i0 >= n:
            break
        blk = np.asarray(mm[i0:i1])
        full = (i1 - i0) // w
        if full:
            xs, ys = blk[:full * w, 0], blk[:full * w, 1]
            if np.array_equal(xs, np.tile(np.arange(w, dtype=np.float32), full)) and np.array_equal(ys, np.repeat(np.arange(y0, y0 + full, dtype=np.float32), w)):
                l0[:, y0:y0 + full, :] = _f16(blk[:full * w, _EXTRA:].reshape(full, w, nch).transpose(2, 0, 1))   # the usual row-by-row order
            else:                                              # records not in raster order: place each by its own X, Y
                xi, yi = xs.astype(np.int64), ys.astype(np.int64)
                ok = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h)
                l0[:, yi[ok], xi[ok]] = _f16(blk[:full * w, _EXTRA:][ok].T)
        rest = blk[full * w:]
        if len(rest):                                          # the last, partly acquired row
            xi, yi = rest[:, 0].astype(np.int64), rest[:, 1].astype(np.int64)
            ok = (xi >= 0) & (xi < w) & (yi >= 0) & (yi < h)
            l0[:, yi[ok], xi[ok]] = _f16(rest[:, _EXTRA:][ok].T)
        done_rows = y0 + band
        progress(min(1.0, done_rows / h) * 0.8)
    l0.flush()
    nz = _levels(w, h)
    prev = l0
    for z in range(1, nz):
        ph, pw = prev.shape[1:]
        nh, nw = (ph + 1) // 2, (pw + 1) // 2
        cur = np.lib.format.open_memmap(d / f"l{z}.npy", mode="w+", dtype=np.float16, shape=(nch, nh, nw))
        for c in range(nch):
            a = np.zeros((nh * 2, nw * 2), np.float32)
            a[:ph, :pw] = prev[c]
            cur[c] = ((a[0::2, 0::2] + a[1::2, 0::2] + a[0::2, 1::2] + a[1::2, 1::2]) / 4).astype(np.float16)
        cur.flush()
        prev = cur
        progress(0.8 + 0.2 * z / nz)
    hi, frac = [], []
    for z in range(nz):                                        # default intensity limit per level and channel: the 99.5th percentile of the non-zero pixels (a strided sample).
        lv = l0 if z == 0 else np.load(d / f"l{z}.npy", mmap_mode="r")   # Coarser levels are means of 2 x 2 pixels, so sparse counts get dimmer; each level has its own limit.
        st = max(1, int(np.sqrt(lv.shape[1] * lv.shape[2] / 1.5e6)))
        row = []
        for c in range(nch):
            v = np.asarray(lv[c, ::st, ::st], np.float32).ravel()
            v = v[v > 0]
            if z == 0:
                frac.append(float(v.size) / max(1, (lv.shape[1] // st + 1) * (lv.shape[2] // st + 1)))
            row.append(float(np.percentile(v, 99.5)) if v.size > 50 else (float(v.max()) if v.size else 1.0))
        hi.append(row)
    meta = {"key": r["key"], "id": r["id"], "name": r["name"], "w": w, "h": h, "rows": r["rows"], "complete": r["complete"], "nz": nz, "tile": TILE,
            "channels": nch, "hi": hi, "frac_nonzero": frac, "stamp": _STAMP}
    (d / "meta.json").write_text(json.dumps(meta))
    return meta


def convert(mcd, root, log=print, progress=lambda f, msg="": None, info=None):
    """Convert every acquisition of the .mcd into the cache store; returns the store folder. Acquisitions already converted are kept."""
    info = info or scan(mcd)
    d = store_dir(root, mcd)
    d.mkdir(parents=True, exist_ok=True)
    nch = len(info["channels"])
    tot = sum(r["records"] for r in info["rois"]) or 1
    done = 0
    for r in info["rois"]:
        rd = d / r["key"]
        ok = False
        try:
            ok = json.loads((rd / "meta.json").read_text()).get("stamp") == _STAMP
        except (OSError, ValueError):
            pass
        share = r["records"] / tot
        if ok:
            log(f"IMC {r['name']} ({r['key']}): already converted")
        else:
            note = "" if r["complete"] else f" (incomplete: {r['rows']:,} of {r['h']:,} rows)"
            log(f"IMC {r['name']} ({r['key']}): {r['w']:,} x {r['h']:,} um, {nch} channels{note}")
            _convert_roi(Path(mcd), r, rd, nch, log, lambda f, d0=done, sh=share: progress(d0 + f * sh, r["name"]))
        done += share
        progress(done, r["name"])
    (d / "index.json").write_text(json.dumps({"file": info["file"], "size": info["size"], "mtime": info["mtime"], "channels": info["channels"],
                                              "rois": [{k: v for k, v in r.items() if k not in ("start", "end")} for r in info["rois"]]}))
    return d


def read_store(d):
    try:
        return json.loads((Path(d) / "index.json").read_text())
    except (OSError, ValueError):
        return None


_mm = {}


def _level(d, key, z):
    k = (str(d), key, z)
    a = _mm.get(k)
    if a is None:
        if len(_mm) > 64:
            _mm.clear()
        a = _mm[k] = np.load(Path(d) / key / f"l{z}.npy", mmap_mode="r")
    return a


def _meta(d, key):
    return json.loads((Path(d) / key / "meta.json").read_text())


def _compose(a, chs, cols, his, gamma):
    out = np.zeros(a.shape[1:] + (3,), np.float32)
    for i, (c, col, hi) in enumerate(zip(chs, cols, his)):
        f = np.clip(a[i].astype(np.float32) / max(hi, 1e-3), 0, 1)
        if gamma != 1.0:
            f = f ** gamma
        out += f[..., None] * np.array(if_image._rgb(col), np.float32)
    return np.clip(out, 0, 255).astype(np.uint8)


def _jpeg(rgb, quality=85):
    b = io.BytesIO()
    Image.fromarray(rgb, "RGB").save(b, "JPEG", quality=quality)
    return b.getvalue()


def _his(m, z, chs, bright):
    """Intensity limit of each channel at level z: the level's default divided by the brightness factor (1 = default; 2 = twice as bright)."""
    return [m["hi"][z][c] / max(b, 1e-3) for c, b in zip(chs, bright)]


def tile(d, key, z, tx, ty, chs, cols, bright, gamma=1.0):
    """JPEG of one tile (TILE pixels) of the chosen channels added up, each in its colour; None outside the image."""
    m = _meta(d, key)
    if not 0 <= z < m["nz"] or tx < 0 or ty < 0:
        return None
    arr = _level(d, key, z)
    h, w = arr.shape[1:]
    x0, y0 = tx * TILE, ty * TILE
    if x0 >= w or y0 >= h:
        return None
    sub = np.stack([np.asarray(arr[c, y0:y0 + TILE, x0:x0 + TILE]) for c in chs]) if chs else np.zeros((0, min(TILE, h - y0), min(TILE, w - x0)), np.float16)
    return _jpeg(_compose(sub, chs, cols, _his(m, z, chs, bright), gamma))


def fit(d, key, maxw, maxh, chs, cols, bright, gamma=1.0):
    """JPEG of the whole acquisition, at most maxw x maxh pixels (aspect ratio kept), for the small side views of the tile windows."""
    m = _meta(d, key)
    z = 0
    while z + 1 < m["nz"] and (m["w"] / 2 ** z > 2 * maxw or m["h"] / 2 ** z > 2 * maxh):
        z += 1
    arr = _level(d, key, z)
    sub = np.stack([np.asarray(arr[c]) for c in chs]) if chs else np.zeros((0,) + arr.shape[1:], np.float16)
    rgb = Image.fromarray(_compose(sub, chs, cols, _his(m, z, chs, bright), gamma), "RGB")
    s = min(maxw / rgb.width, maxh / rgb.height, 1.0)
    if s < 1.0:
        rgb = rgb.resize((max(1, round(rgb.width * s)), max(1, round(rgb.height * s))), Image.LANCZOS)
    b = io.BytesIO()
    rgb.save(b, "JPEG", quality=85)
    return b.getvalue()
