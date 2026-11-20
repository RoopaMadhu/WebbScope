"""GO and Reactome gene sets (MSigDB, bundled in data/genesets.json.gz; rebuilt with tools/build_genesets.py) offered as signatures.

Every set stays in the app; a project lists only the sets with at least MIN_GENES genes on its panel (the rest can still be looked up by id).
Ids are "go:<MSigDB name>" and "reactome:<MSigDB name>". Gene symbols are matched as the curated signatures are: the panel's symbol, or its
older / newer alternative (signature_defs.ALT_SYMBOLS), and EXCLUDED_GENES are left out.
"""
import gzip, json, threading
from pathlib import Path
from .signature_defs import ALT_SYMBOLS, EXCLUDED_GENES

MIN_GENES = 5                       # a set is listed for a project only with at least this many genes on its panel
DATA = Path(__file__).resolve().parent / "data" / "genesets.json.gz"
PREFIX = {"GO": "go", "Reactome": "reactome"}
ASPECT = {"GOBP": "GO biological process", "GOCC": "GO cellular component", "GOMF": "GO molecular function", "REACTOME": "Reactome pathway"}
_lock, _raw = threading.Lock(), None


def raw():
    """The bundled file, read once: {release, source, license, collections: {GO: [[name, [genes]]], Reactome: [...]}}."""
    global _raw
    with _lock:
        if _raw is None:
            with gzip.open(DATA, "rt") as f:
                _raw = json.load(f)
        return _raw


def pretty(name):
    """GOBP_RESPONSE_TO_TYPE_I_INTERFERON -> "Response to type I interferon" (roman numerals and short all-caps words kept)."""
    words = name.split("_", 1)[1].split("_") if "_" in name else [name]
    keep = {"I", "II", "III", "IV", "V", "MHC", "DNA", "RNA", "TGF", "NF", "TNF", "EGFR", "FGFR", "PI3K", "AKT", "MAPK", "ERK", "JAK", "STAT", "TLR", "IL", "ER", "ECM", "EMT", "GPCR", "ATP", "GTP", "NAD", "CD", "T", "B", "NK", "RHO", "SUMO", "MTOR", "AMPK", "HIF"}
    out = [w if (w in keep or any(c.isdigit() for c in w)) else w.lower() for w in words]
    s = " ".join(out)
    return s[:1].upper() + s[1:]


class Resolved:
    """The gene sets resolved against one project's panel."""

    def __init__(self, panel_genes):
        panel = set(panel_genes)
        r = raw()
        self.release = r.get("release")
        self.by_id, self.listed, self.n_small = {}, {c: [] for c in PREFIX}, {c: 0 for c in PREFIX}
        for coll, sets in r["collections"].items():
            for name, genes in sets:
                on = []
                for g in genes:
                    f = g if g in panel else ALT_SYMBOLS.get(g) if ALT_SYMBOLS.get(g) in panel else None
                    if f and f not in EXCLUDED_GENES:
                        on.append(f)
                on = list(dict.fromkeys(on))
                sid = f"{PREFIX[coll]}:{name}"
                self.by_id[sid] = (coll, name, on, genes)
                if len(on) >= MIN_GENES:
                    self.listed[coll].append([sid, pretty(name), len(on)])
                else:
                    self.n_small[coll] += 1
        for c in self.listed:
            self.listed[c].sort(key=lambda x: x[1].lower())

    def get(self, sid):
        """The full signature record of a set (the same fields as a curated signature), or None."""
        e = self.by_id.get(sid)
        if e is None:
            return None
        coll, name, on, allg = e
        have = set(on) | {ALT_SYMBOLS.get(g) for g in on}
        missing = [g for g in allg if g not in have and g not in EXCLUDED_GENES]
        aspect = ASPECT.get(name.split("_", 1)[0], coll)
        return {"id": sid, "group": coll, "name": pretty(name), "msigdb_name": name, "genes": list(on), "missing": missing, "n_defined": len(on) + len(missing),
                "desc": f"{aspect}: {name} (MSigDB {self.release}).", "note": "" if len(on) >= MIN_GENES else f"Fewer than {MIN_GENES} of its genes are on this panel, so it is not listed for this project.",
                "url": f"https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/{name}.html"}

    def listing(self):
        return {"release": self.release, "min_genes": MIN_GENES, "collections": self.listed, "n_small": self.n_small, "license": raw().get("license")}
