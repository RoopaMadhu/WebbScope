"""Minimal reader for R .rds files (saveRDS, serialization version 2 or 3, gzip / bzip2 / xz / uncompressed), in pure Python + numpy / scipy / pandas.

Reads what a plain list of tables is made of: lists, numeric / integer / logical / character vectors, matrices, factors, data frames, and the Matrix
package's sparse matrices (dgCMatrix, dgTMatrix, dgRMatrix, dgeMatrix, lgCMatrix, ngCMatrix). Anything else (environments, functions, S4 classes such as
a Seurat object) is rejected with a message that names the class, so that the user can export a plain list instead.

read_rds(path) -> Python objects: a list with names becomes a dict; a matrix a numpy array / scipy sparse matrix plus its dimnames (see `Matrix`);
a data frame a pandas DataFrame (row names as index); a factor a list of strings; a plain vector a numpy array (NA of a numeric vector = NaN,
of an integer / logical vector = NaN after conversion to float) or a list of strings.
"""
import bz2, gzip, io, lzma, struct
import numpy as np
import pandas as pd
import scipy.sparse as sp
from .importer import ImportError_

NA_INT = -2147483648


class RObj:
    """A parsed R object before conversion: kind (R type name), value, attrs ({name: RObj / value})."""
    __slots__ = ("kind", "value", "attrs")

    def __init__(self, kind, value=None, attrs=None):
        self.kind, self.value, self.attrs = kind, value, attrs or {}


class Matrix:
    """A matrix with its row and column names (None when the object had none)."""
    def __init__(self, x, rownames=None, colnames=None):
        self.x, self.rownames, self.colnames = x, rownames, colnames

    @property
    def shape(self):
        return self.x.shape


class _Reader:
    def __init__(self, b):
        self.b, self.p, self.refs = b, 0, []

    def int(self):
        v = struct.unpack_from(">i", self.b, self.p)[0]; self.p += 4; return v

    def take(self, n):
        v = self.b[self.p:self.p + n]
        if len(v) != n:
            raise ImportError_("the .rds file ends early (truncated or damaged)")
        self.p += n
        return v

    def header(self):
        if self.take(2) != b"X\n":
            raise ImportError_("not an R serialization (.rds) file, or saved as ascii: use saveRDS(x, file) without ascii = TRUE")
        ver = self.int(); self.int(); self.int()
        if ver == 3:
            self.take(self.int())                      # native encoding
        elif ver != 2:
            raise ImportError_(f"unsupported R serialization version {ver}")

    def length(self):
        n = self.int()
        if n == -1:                                    # long vector: two ints
            n = (self.int() << 32) + (self.int() & 0xFFFFFFFF)
        return n

    def string(self):                                   # a CHARSXP after its flags have been read
        n = self.int()
        if n == -1:
            return None
        raw = self.take(n)
        try:
            return raw.decode("utf-8")
        except UnicodeDecodeError:
            return raw.decode("latin-1")

    def pairlist_tail(self, flags):
        """LISTSXP / LANGSXP / ... cells: [attr] [tag] car, then cdr. Returns [(tag, value)] and the attributes of the first cell."""
        out, first_attr = [], None
        while True:
            attr = self.item() if flags & 0x200 else None
            if first_attr is None and attr is not None:
                first_attr = attr
            tag = self.item() if flags & 0x400 else None
            car = self.item()
            out.append((tag.value if isinstance(tag, RObj) and tag.kind == "symbol" else None, car))
            flags = self.int()
            t = flags & 0xFF
            if t == 254:                                # NILVALUE: end of list
                return out, first_attr
            if t not in (2, 6, 239, 240):               # a non-list tail (dotted pair): read it and stop
                self.p -= 4
                self.item()
                return out, first_attr

    def attrs_of(self, a):
        d = {}
        if isinstance(a, list):
            for tag, v in a:
                if tag is not None:
                    d[tag] = v
        return d

    def item(self):
        flags = self.int()
        t = flags & 0xFF
        if t == 254 or t == 242 or t == 241 or t == 253 or t == 252 or t == 251 or t == 250:
            return RObj("nil" if t == 254 else "env-constant")
        if t == 255:                                    # reference
            i = flags >> 8
            if i == 0:
                i = self.int()
            return self.refs[i - 1]
        if t == 1:                                      # symbol
            self.int()                                  # the CHARSXP flags
            s = RObj("symbol", self.string())
            self.refs.append(s)
            return s
        if t in (249, 248, 247):                        # namespace / package / persistent: STRSXP-like info
            self.int()
            n = self.int()
            info = []
            for _ in range(n):
                self.int(); info.append(self.string())
            o = RObj("env-constant", info)
            self.refs.append(o)
            return o
        if t == 4:                                      # environment
            o = RObj("environment")
            self.refs.append(o)
            self.int()
            for _ in range(4):
                self.item()
            return o
        if t in (2, 6, 239, 240):                       # pairlist / language
            cells, attr = self.pairlist_tail(flags)
            return RObj("pairlist" if t in (2, 239) else "language", cells, self.attrs_of(attr.value) if attr is not None else {})
        if t == 238:                                    # ALTREP
            return self.altrep()
        if t == 9:                                      # CHARSXP
            return RObj("char", self.string())
        if t in (3, 5, 7, 8, 12, 21, 22, 23):           # closures, promises, builtins, bytecode, ...
            raise ImportError_("the .rds contains an R function / environment object, which cannot be read: save a plain list (saveRDS(list(...), file))")
        n = None
        if t == 25:                                     # S4 object: no data, attributes only
            o = RObj("S4")
        elif t in (10, 13):
            n = self.length(); v = np.frombuffer(self.take(4 * n), ">i4").astype(np.int64); o = RObj("logical" if t == 10 else "integer", v)
        elif t == 14:
            n = self.length(); o = RObj("real", np.frombuffer(self.take(8 * n), ">f8").astype(np.float64))
        elif t == 15:
            n = self.length(); self.take(16 * n); o = RObj("complex", None)
        elif t == 16:
            n = self.length(); o = RObj("character", [self.item().value for _ in range(n)])
        elif t in (19, 20):
            n = self.length(); o = RObj("list", [self.item() for _ in range(n)])
        elif t == 24:
            n = self.length(); o = RObj("raw", bytes(self.take(n)))
        else:
            raise ImportError_(f"the .rds contains an object of R type {t} that cannot be read")
        if flags & 0x200:
            o.attrs = self.attrs_of(self.item().value)
        return o

    def altrep(self):
        info = self.item()
        state = self.item()
        attr = self.item()
        cls = info.value[0][1].value if info.kind == "pairlist" else ""
        if cls in ("compact_intseq", "compact_realseq"):
            n, start, inc = state.value[0], state.value[1], state.value[2]
            v = start + inc * np.arange(int(n))
            o = RObj("integer" if cls == "compact_intseq" else "real", v.astype(np.int64 if cls == "compact_intseq" else np.float64))
        elif cls == "deferred_string":
            arg = _first(state)
            o = RObj("character", [None if (isinstance(x, float) and np.isnan(x)) or x == NA_INT else (str(int(x)) if float(x).is_integer() else repr(float(x))) for x in arg.value])
        elif cls.startswith("wrap_"):
            o = _first(state)
            o = RObj(o.kind, o.value, dict(o.attrs))
        else:
            raise ImportError_(f"the .rds uses an R vector representation ({cls or 'unknown'}) that cannot be read: re-save with saveRDS(x, file, version = 2)")
        if attr.kind == "pairlist":
            o.attrs.update(self.attrs_of(attr.value))
        return o


def _first(state):
    """The first element of an ALTREP state (a pairlist or a list)."""
    v = state.value[0]
    return v[1] if isinstance(v, tuple) else v


def _names(o):
    a = o.attrs.get("names")
    return list(a.value) if a is not None else None


def _dimnames(o):
    dn = o.attrs.get("dimnames") or o.attrs.get("Dimnames")
    if dn is None or dn.kind != "list":
        return None, None
    r = [None if e.kind == "nil" else list(e.value) for e in dn.value]
    return (r + [None, None])[0], (r + [None, None])[1]


def _cls(o):
    c = o.attrs.get("class")
    return list(c.value) if c is not None and c.kind == "character" else []


def _vec(o):
    v = o.value
    if o.kind in ("integer", "logical"):
        v = v.astype(np.float64)
        v[o.value == NA_INT] = np.nan
    return v


class LazyDict(dict):
    """A named R list whose elements are converted when first read, so that an element this program never asks for (and cannot read) does no harm."""
    def __init__(self, names, objs):
        super().__init__((n, o) for n, o in zip(names, objs))
        self._done = {}

    def __getitem__(self, k):
        if k not in self._done:
            self._done[k] = convert(super().__getitem__(k))
        return self._done[k]

    def get(self, k, d=None):
        return self[k] if k in self else d


def convert(o):
    """RObj -> Python (dict / DataFrame / Matrix / numpy array / list of str / None)."""
    k, cls = o.kind, _cls(o)
    if k == "nil":
        return None
    if k == "S4":
        c = cls[0] if cls else "?"
        g = lambda n: o.attrs[n].value if n in o.attrs else None
        dims = tuple(int(i) for i in g("Dim")) if g("Dim") is not None else None
        rn, cn = _dimnames(o)
        if c in ("dgCMatrix", "lgCMatrix", "ngCMatrix"):
            x = g("x") if g("x") is not None else np.ones(len(g("i")))
            m = sp.csc_matrix((np.asarray(x, float), np.asarray(g("i"), np.int64), np.asarray(g("p"), np.int64)), shape=dims)
        elif c in ("dgTMatrix", "lgTMatrix", "ngTMatrix"):
            x = g("x") if g("x") is not None else np.ones(len(g("i")))
            m = sp.coo_matrix((np.asarray(x, float), (np.asarray(g("i"), np.int64), np.asarray(g("j"), np.int64))), shape=dims).tocsc()
        elif c in ("dgRMatrix", "lgRMatrix", "ngRMatrix"):
            x = g("x") if g("x") is not None else np.ones(len(g("j")))
            m = sp.csr_matrix((np.asarray(x, float), np.asarray(g("j"), np.int64), np.asarray(g("p"), np.int64)), shape=dims).tocsc()
        elif c == "dgeMatrix":
            m = np.asarray(g("x"), float).reshape(dims, order="F")
        else:
            raise ImportError_(f"an object of R class {c} cannot be read: save plain matrices / data frames in a list (as.matrix(), as.data.frame())")
        return Matrix(m, rn, cn)
    if k == "list":
        names = _names(o)
        if "data.frame" in cls:
            cols = {n: convert(v) for n, v in zip(names, o.value)}
            cols = {n: (v["values"] if isinstance(v, dict) and "__names__" in v else v) for n, v in cols.items()}
            rn = o.attrs.get("row.names")
            idx = None
            if rn is not None and rn.kind == "character":
                idx = list(rn.value)
            df = pd.DataFrame({n: (v if not isinstance(v, np.ndarray) or v.ndim == 1 else list(v)) for n, v in cols.items()})
            if idx is not None and len(idx) == len(df):
                df.index = idx
            return df
        return LazyDict(names, o.value) if names else [convert(v) for v in o.value]
    if k == "character":
        return list(o.value)
    if k in ("integer", "logical", "real"):
        if "factor" in cls:
            lv = list(o.attrs["levels"].value)
            return [None if i == NA_INT else lv[i - 1] for i in o.value]
        dim = o.attrs.get("dim")
        v = _vec(o)
        if dim is not None and len(dim.value) == 2:
            rn, cn = _dimnames(o)
            return Matrix(v.reshape(tuple(int(i) for i in dim.value), order="F"), rn, cn)
        nm = _names(o)
        return v if nm is None else {"__names__": nm, "values": v}
    if k == "pairlist":
        return {t: convert(v) for t, v in o.value if t}
    raise ImportError_(f"an R object of kind {k} cannot be read")


def read_rds(path):
    """The converted content of an .rds file."""
    try:
        raw = open(path, "rb").read()
    except OSError as e:
        raise ImportError_(f"cannot read {path}: {e}")
    for opener in (gzip.decompress, bz2.decompress, lzma.decompress, lambda b: b):
        try:
            b = opener(raw)
        except Exception:
            continue
        if b[:2] == b"X\n":
            break
    else:
        raise ImportError_(f"{path} is not an .rds file (it does not start with R's serialization header)")
    r = _Reader(b)
    r.header()
    return convert(r.item())
