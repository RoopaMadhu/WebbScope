"""Xenium morphology (immunofluorescence) images: pyramidal OME-TIFF, 16-bit, tiled, JPEG 2000 compressed (TIFF compression 34712).

The image is already in the Xenium frame (pixel = x um / pixel size, y down), so it needs no registration. It is not copied or converted:
tiles are read straight from the file and rendered on request (8-bit, chosen colour and intensity range), using the file's own pyramid (SubIFDs).
Pillow cannot open these files as images, so the TIFF container is parsed here; the JPEG 2000 tiles themselves are decoded by Pillow.

Read: BigTIFF and classic TIFF, little or big endian, one grey channel per IFD, compression JPEG 2000 / none / deflate.
Checked on one real file set (Xenium 1.3 output: morphology_focus.ome.tif and morphology_mip.ome.tif, 51,148 x 47,696 px, 0.2125 um, DAPI).
The z-stack morphology.ome.tif is not used. Multi-channel files (channels as separate IFDs) are read by IFD order but were not checked on real data.
"""
import io, json, os, re, struct, threading, zlib
from collections import OrderedDict
from pathlib import Path
import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
INDEX = "if.json"
ROLE_FILES = {"morphology_focus.ome.tif": "if_focus", "morphology_mip.ome.tif": "if_mip"}
ROLE_RE = re.compile(r"^morphology_focus/([^/]+\.ome\.tif)$|^(morphology_focus_\d+\.ome\.tif)$")
COLORS = {"white": (255, 255, 255), "blue": (80, 140, 255), "cyan": (0, 220, 255), "green": (60, 230, 90), "red": (255, 70, 70), "magenta": (255, 60, 220),
          "yellow": (255, 230, 60), "orange": (255, 150, 40)}
_T = {1: ("B", 1), 7: ("B", 1), 3: ("H", 2), 4: ("I", 4), 8: ("h", 2), 9: ("i", 4), 16: ("Q", 8), 17: ("q", 8), 18: ("Q", 8)}


class IFError(Exception):
    pass


def role_of(rel):
    """Role of a path inside a Xenium output (for the folder / tar scan), or None."""
    if rel in ROLE_FILES:
        return ROLE_FILES[rel]
    m = ROLE_RE.match(rel)
    return "if_" + Path(m.group(1) or m.group(2)).name.replace(".ome.tif", "") if m else None


def pick(src):
    """The image to use from a {role: Path} of a run: the focus image (else the first morphology_focus_* file, else the MIP). None if there is none."""
    for r in ("if_focus",) + tuple(sorted(k for k in src if k.startswith("if_") and k not in ("if_focus", "if_mip"))) + ("if_mip",):
        if r in src:
            return Path(src[r])
    return None


def pick_all(src):
    """All the images to use as channels from a {role: Path} of a run, in order. Newer Xenium output has one single-channel file per stain in morphology_focus/
    (ch0000_dapi.ome.tif, ch0001_atp1a1_cd45_e-cadherin.ome.tif, ...): all of them. Older output has one morphology_focus.ome.tif; else the MIP. [] if none."""
    if "if_focus" in src:
        return [Path(src["if_focus"])]
    many = sorted(k for k in src if k.startswith("if_") and k not in ("if_focus", "if_mip"))
    if many:
        return [Path(src[k]) for k in many]
    return [Path(src["if_mip"])] if "if_mip" in src else []


_STAIN = {"dapi": "DAPI", "atp1a1": "ATP1A1", "cd45": "CD45", "e-cadherin": "E-cadherin", "18s": "18S", "alphasma": "αSMA", "vimentin": "Vimentin"}


def stain_name(path):
    """Display name from a morphology file name: ch0001_atp1a1_cd45_e-cadherin.ome.tif -> 'ATP1A1 / CD45 / E-cadherin'."""
    stem = re.sub(r"\.ome\.tiff?$", "", Path(path).name, flags=re.I)
    m = re.match(r"^ch\d+_(.+)$", stem)
    if not m:
        return stem.replace("morphology_", "").replace("_", " ")
    return " / ".join(_STAIN.get(t.lower(), t) for t in m.group(1).split("_"))


class _Tiff:
    def __init__(self, path):
        self.path = str(path)
        self.f = open(path, "rb")
        h = self.f.read(16)
        if h[:2] not in (b"II", b"MM"):
            raise IFError(f"{Path(path).name} is not a TIFF file")
        self.e = "<" if h[:2] == b"II" else ">"
        magic = struct.unpack(self.e + "H", h[2:4])[0]
        if magic not in (42, 43):
            raise IFError(f"{Path(path).name} is not a TIFF file")
        self.big = magic == 43
        self.first = struct.unpack(self.e + "Q", h[8:16])[0] if self.big else struct.unpack(self.e + "I", h[4:8])[0]

    def ifd(self, off):
        e, f = self.e, self.f
        f.seek(off)
        n, = struct.unpack(e + ("Q" if self.big else "H"), f.read(8 if self.big else 2))
        es = 20 if self.big else 12
        raw = f.read(n * es)
        nxt, = struct.unpack(e + ("Q" if self.big else "I"), f.read(8 if self.big else 4))
        tags = {}
        for i in range(n):
            b = raw[i * es:(i + 1) * es]
            tag, typ, cnt = struct.unpack(e + ("HHQ" if self.big else "HHI"), b[:12 if self.big else 8])
            tags[tag] = (typ, cnt, b[12:] if self.big else b[8:])
        return tags, nxt

    def val(self, t):
        typ, cnt, vb = t
        if typ == 2:
            n = cnt
            raw = vb[:n] if n <= (8 if self.big else 4) else self._at(vb, n)
            return raw.decode("utf8", "replace").rstrip("\x00")
        fmt, sz = _T.get(typ, (None, 0))
        if not sz:
            raise IFError(f"unsupported TIFF field type {typ}")
        raw = vb[:cnt * sz] if cnt * sz <= (8 if self.big else 4) else self._at(vb, cnt * sz)
        return struct.unpack(self.e + fmt * cnt, raw)

    def _at(self, vb, n):
        o, = struct.unpack(self.e + ("Q" if self.big else "I"), vb[:8 if self.big else 4])
        self.f.seek(o)
        return self.f.read(n)

    def close(self):
        self.f.close()


class Level:
    def __init__(self, tags, tif):
        g = lambda k, d=None: tif.val(tags[k]) if k in tags else d
        self.w, self.h = g(256)[0], g(257)[0]
        self.bits, self.comp = g(258, (16,))[0], g(259, (1,))[0]
        self.tw, self.th = (g(322, (0,))[0], g(323, (0,))[0])
        if not self.tw or 324 not in tags:
            raise IFError("the image is stored in strips, not tiles; only tiled morphology images are supported")
        if g(277, (1,))[0] != 1:
            raise IFError("the image has several samples per pixel (RGB); only single-channel images are supported")
        self.offs, self.cnts = g(324), g(325)
        self.nx, self.ny = -(-self.w // self.tw), -(-self.h // self.th)


class Image16:
    """One channel of one file: the levels of its pyramid (level 0 = full resolution)."""
    def __init__(self, path, ifd_index=0):
        self.path = str(path)
        t = _Tiff(path)
        try:
            off = t.first
            for _ in range(ifd_index):
                _, off = t.ifd(off)
                if not off:
                    raise IFError(f"{Path(path).name} has no image number {ifd_index}")
            tags, _ = t.ifd(off)
            self.levels = [Level(tags, t)]
            for so in (t.val(tags[330]) if 330 in tags else ()):
                st, _ = t.ifd(so)
                self.levels.append(Level(st, t))
            self.ome = t.val(tags[270]) if 270 in tags else ""
        finally:
            t.close()
        L0 = self.levels[0]
        if L0.comp not in (1, 8, 32946, 34712):
            raise IFError(f"{Path(path).name}: TIFF compression {L0.comp} is not supported (JPEG 2000, deflate and none are)")
        self.w, self.h, self.nz, self.tsize = L0.w, L0.h, len(self.levels), L0.tw
        if any(lv.tw != L0.tw or lv.th != L0.th for lv in self.levels) or L0.tw != L0.th:
            raise IFError("the pyramid levels do not share one square tile size")

    def raw(self, z, tx, ty):
        """The tile as a 2-D array of the file's integer type (zeros outside the image), or None when the tile is not in the file."""
        lv = self.levels[z]
        if not (0 <= tx < lv.nx and 0 <= ty < lv.ny):
            return None
        k = ty * lv.nx + tx
        with open(self.path, "rb") as f:
            f.seek(lv.offs[k])
            b = f.read(lv.cnts[k])
        dt = np.uint16 if lv.bits == 16 else np.uint8
        if lv.comp == 34712:
            return np.asarray(Image.open(io.BytesIO(b)))
        if lv.comp in (8, 32946):
            b = zlib.decompress(b)
        return np.frombuffer(b, dtype=np.dtype(dt).newbyteorder("<")).reshape(lv.th, lv.tw)


def info(path):
    """What the viewer needs of one morphology file: size, pixel size (um), tile size, number of levels, channel names."""
    im = Image16(path)
    px = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', im.ome)
    names = re.findall(r'<Channel[^>]*?Name="([^"]*)"', im.ome)
    return {"w": im.w, "h": im.h, "px": float(px.group(1)) if px else None, "tsize": im.tsize, "nz": im.nz, "names": names, "bits": im.levels[0].bits}


def default_range(path, ifd_index=0):
    """(lo, hi): a sensible intensity window from a coarse level (the coarsest one that is still at least 6000 px across: blurring at coarser levels lowers the peaks, finer ones cost more tiles), non-zero pixels only."""
    im = Image16(path, ifd_index)
    z = max((i for i, lv in enumerate(im.levels) if max(lv.w, lv.h) >= 6000), default=len(im.levels) - 1)
    lv = im.levels[z]
    vals = []
    for ty in range(lv.ny):
        for tx in range(lv.nx):
            a = im.raw(z, tx, ty)
            if a is not None and a.size:
                vals.append(a[a > 0].ravel())
    v = np.concatenate(vals) if vals else np.zeros(0)
    if not v.size:
        return 0, 1
    lo, hi = float(np.percentile(v, 1)), float(np.percentile(v, 99.5))
    return int(lo), int(max(hi, lo + 1))


def n_ifds(path, cap=64):
    """Number of top-level image directories (planes / channels) the file really holds."""
    t = _Tiff(path)
    try:
        n, off = 0, t.first
        while off and n < cap:
            _, off = t.ifd(off)
            n += 1
        return n
    finally:
        t.close()


def describe(path, rel_to=None, index=0):
    """Entry of if.json for one file: file (absolute, plus the path relative to rel_to when that is possible), geometry and channels with their default range.
    The channels are the image directories the file holds. The OME header can list more stains than the file contains (Xenium 2.x: every morphology_focus_000N.ome.tif
    lists all four stains but holds one); such a file is one channel, named by `index`, its position among the files."""
    path = Path(path).resolve()
    meta = info(path)
    if not meta["px"]:
        raise IFError(f"{path.name}: the pixel size is not in the file's OME metadata")
    names = meta["names"] or [f"Channel {index + 1}"]                 # no stain in the metadata: just the channel number
    n_dir = n_ifds(path)
    if n_dir < len(names):
        if n_dir == 1:       # a header that lists several stains is shared by all the files (this file is the one at `index`); a header with one name names this file
            names = [names[index]] if index < len(names) else [f"Channel {index + 1}"]
        else:
            names = names[:n_dir]
    elif n_dir == 1 and len(names) == 1 and index > 0 and not meta["names"]:
        names = [f"Channel {index + 1}"]
    chans = []
    for i, n in enumerate(names):
        lo, hi = default_range(path, i)
        chans.append({"name": f"{n} ({path.name.replace('.ome.tif', '').replace('morphology_', '')})", "ome_name": n, "ifd": i, "lo": lo, "hi": hi})
    ent = {"file": str(path), **{k: meta[k] for k in ("w", "h", "px", "tsize", "nz")}, "channels": chans}
    if rel_to is not None:
        try:
            ent["rel"] = os.path.relpath(path, Path(rel_to).resolve())
        except ValueError:
            pass
    return ent


def describe_many(paths, rel_to=None):
    """if.json entry for several single-channel files of the same image (one per stain): every file becomes a channel, named after its stain, with its own file."""
    ents = []
    for i, p in enumerate(paths):
        e = describe(p, rel_to, i)
        if ents and (e["w"], e["h"], e["nz"]) != (ents[0]["w"], ents[0]["h"], ents[0]["nz"]):
            raise IFError(f"{Path(p).name} is {e['w']} x {e['h']} px, but {Path(ents[0]['file']).name} is {ents[0]['w']} x {ents[0]['h']}: the channel files are not the same image")
        ents.append(e)
    if not ents:
        raise IFError("no IF image files")
    out = dict(ents[0])
    out["channels"] = []
    for e, src in zip(ents, paths):
        for c in e["channels"]:
            c = dict(c)
            c["file"], c["rel"] = e["file"], e.get("rel")
            if len(ents) > 1 or not c.get("name"):
                by_file = re.match(r"^ch\d+_", Path(src).name) is not None                 # the stain is in the file name (ch0001_atp1a1_cd45_e-cadherin): use it; else the header's name for this file
                one = stain_name(src) if (by_file or not c.get("ome_name")) else c["ome_name"]
                c["name"] = one if len(e["channels"]) == 1 else f"{c['name']} ({stain_name(src)})"      # named from the path as given (a symlink's target may be called anything)
            out["channels"].append(c)
    return out


def write_index(folder, images, samples):
    """<folder>/if.json: images = {key: describe()}, samples = {sample: key}."""
    Path(folder, INDEX).write_text(json.dumps({"images": images, "samples": samples}, indent=1))


def read_index(folder):
    f = Path(folder) / INDEX
    if not f.is_file():
        return None
    d = json.loads(f.read_text())
    for ent in d["images"].values():          # the file where it was recorded, else the same place relative to the project folder
        if not Path(ent["file"]).is_file() and ent.get("rel") and (Path(folder) / ent["rel"]).is_file():
            ent["file"] = str((Path(folder) / ent["rel"]).resolve())
        for c in ent["channels"]:             # channels that come from their own file
            if c.get("file") and not Path(c["file"]).is_file() and c.get("rel") and (Path(folder) / c["rel"]).is_file():
                c["file"] = str((Path(folder) / c["rel"]).resolve())
    return d


def sample_list(idx, sample):
    """Entries for /api/if: one per channel of the sample's image, in the shape of the registered-H&E list (key, w, h, px, x0, y0, tsize, nz)."""
    k = idx["samples"].get(sample)
    ent = idx["images"].get(k) if k else None
    if not ent or not all(Path(c.get("file") or ent["file"]).is_file() for c in ent["channels"]):
        return []
    return [{"key": f"if:{k}:{i}", "name": c["name"], "w": ent["w"], "h": ent["h"], "px": ent["px"], "x0": 0, "y0": 0, "tsize": ent["tsize"], "nz": ent["nz"],
             "lo": c["lo"], "hi": c["hi"]} for i, c in enumerate(ent["channels"])]


# ---------------------------------------------------------------- tiles on request
_IMG, _RAW, _BLACK = {}, OrderedDict(), {}
_LOCK = threading.Lock()
RAW_CAP = 300 << 20                                           # decoded tiles kept, bytes


def _image(path, ifd):
    k = (str(path), ifd)
    with _LOCK:
        if k not in _IMG:
            _IMG[k] = Image16(path, ifd)
        return _IMG[k]


def _raw(path, ifd, z, tx, ty):
    k = (str(path), ifd, z, tx, ty)
    with _LOCK:
        a = _RAW.get(k)
        if a is not None:
            _RAW.move_to_end(k)
            return a
    a = _image(path, ifd).raw(z, tx, ty)
    if a is None:
        return None
    with _LOCK:
        _RAW[k] = a
        while sum(v.nbytes for v in _RAW.values()) > RAW_CAP and len(_RAW) > 1:
            _RAW.popitem(last=False)
    return a


def _rgb(c):
    if c in COLORS:
        return COLORS[c]
    if re.fullmatch(r"[0-9a-fA-F]{6}", c or ""):
        return tuple(int(c[i:i + 2], 16) for i in (0, 2, 4))
    return COLORS["white"]


def render(path, ifd, z, tx, ty, color, lo, hi, gamma=1.0, quality=85):
    """JPEG bytes of one tile, or None when it is outside the image. Intensities lo..hi map to black..colour."""
    im = _image(path, ifd)
    if not 0 <= z < im.nz:
        return None
    a = _raw(path, ifd, z, tx, ty)
    if a is None:
        return None
    rgb = _rgb(color)
    f = np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)
    if gamma != 1.0:
        f = f ** gamma
    out = (f[..., None] * np.array(rgb, np.float32)).astype(np.uint8)
    b = io.BytesIO()
    Image.fromarray(out, "RGB").save(b, "JPEG", quality=quality)
    return b.getvalue()
