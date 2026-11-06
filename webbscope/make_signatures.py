"""Add (or refresh) signatures.json in an existing project folder: python -m webbscope.make_signatures PROJECT_FOLDER
The curated signatures are resolved against the project's own gene list (features.tsv.gz); the viewer re-imports the project on the next open."""
import gzip, sys
from pathlib import Path
from . import xenium


def main(argv=None):
    a = sys.argv[1:] if argv is None else argv
    if len(a) != 1 or not (Path(a[0]) / "features.tsv.gz").exists():
        sys.exit("usage: python -m webbscope.make_signatures PROJECT_FOLDER   (a folder with features.tsv.gz)")
    with gzip.open(Path(a[0]) / "features.tsv.gz", "rt") as f:
        genes = [l.rstrip("\n").split("\t")[0] for l in f if l.strip()]
    xenium.write_signatures(a[0], genes)


if __name__ == "__main__":
    main()
