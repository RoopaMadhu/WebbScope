"""Build webbscope/engine/data/genesets.json.gz, the GO and Reactome gene sets the viewer offers as signatures, from MSigDB's human GMT files.

    python tools/build_genesets.py FOLDER_WITH_GMTS [RELEASE]

FOLDER_WITH_GMTS holds the MSigDB human gene-symbol files of one release, e.g. for 2026.2.Hs (from
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/):
    c5.go.bp.v2026.2.Hs.symbols.gmt  c5.go.cc.v2026.2.Hs.symbols.gmt  c5.go.mf.v2026.2.Hs.symbols.gmt  c2.cp.reactome.v2026.2.Hs.symbols.gmt
MSigDB is distributed under CC BY 4.0 (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp); GO and Reactome are its C5 and C2:CP:REACTOME collections.
"""
import gzip, json, re, sys
from pathlib import Path

FILES = {"GO": ["c5.go.bp", "c5.go.cc", "c5.go.mf"], "Reactome": ["c2.cp.reactome"]}


def main(argv=None):
    a = sys.argv[1:] if argv is None else argv
    if not a:
        sys.exit(__doc__)
    d = Path(a[0])
    release = a[1] if len(a) > 1 else None
    out = {"release": None, "source": "MSigDB (Broad Institute), human gene symbols", "license": "CC BY 4.0", "collections": {}}
    for coll, stems in FILES.items():
        sets = []
        for stem in stems:
            hits = sorted(d.glob(f"{stem}.v*.Hs.symbols.gmt"))
            if release:
                hits = [h for h in hits if f".v{release}." in h.name]
            if len(hits) != 1:
                sys.exit(f"need exactly one {stem}.v<release>.Hs.symbols.gmt in {d}, found {[h.name for h in hits]}")
            out["release"] = re.search(r"\.v(.+?)\.symbols\.gmt$", hits[0].name).group(1)
            for line in hits[0].read_text().splitlines():
                p = line.split("\t")
                if len(p) >= 3:
                    sets.append([p[0], sorted(set(g for g in p[2:] if g))])
        out["collections"][coll] = sorted(sets)
        print(f"{coll}: {len(sets)} sets")
    f = Path(__file__).resolve().parent.parent / "webbscope" / "engine" / "data" / "genesets.json.gz"
    f.parent.mkdir(exist_ok=True)
    with gzip.open(f, "wt") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"wrote {f} ({f.stat().st_size / 1e6:.1f} MB), MSigDB {out['release']}")


if __name__ == "__main__":
    main()
