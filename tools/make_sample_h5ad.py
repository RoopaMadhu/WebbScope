"""Write the sample project as an AnnData .h5ad: the example of the .h5ad input (open the folder or the file on the open page).

    python tools/make_sample_h5ad.py [--src sample_project] [--out sample_project_h5ad] [--boundaries FILE]

Needs anndata (pip install anndata), which WebbScope itself does not use. Reads only --src; writes --out/sample.h5ad and copies the cell
outlines next to it as cell_boundaries.parquet (--out is created), so the folder shows the cells with their outlines. The outlines are
--boundaries if given, else --src/cell_boundaries.parquet (one row per outline point). The example folder holds instead the GeoParquet
form (one row per cell), exported from the R object the sample was cut from:
    Rscript tools/export_cell_boundaries.R sample_project/coarse_obj_for_claude_kidney_edp_1_2.qs sample_project_h5ad/cell_boundaries.parquet --cells sample_project/cells.parquet
  X                raw counts, cells x genes, sparse CSR (float32, whole numbers)
  obs              every column of cells.parquet except the cell id (the index) and x / y
  obsm["spatial"]  x, y of each cell (µm), as scanpy / squidpy store coordinates
  var              gene symbols as the index
  cell_boundaries.parquet   the cell outlines (a separate file: an .h5ad has no place for them), matched to the cells by cell id
The embeddings, regions and other extras of the sample project are left out.
"""
import argparse, gzip, shutil
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.io import mmread

root = Path(__file__).resolve().parent.parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--src", default=str(root / "sample_project"))
ap.add_argument("--out", default=str(root / "sample_project_h5ad"))
ap.add_argument("--boundaries", default=None, help="cell_boundaries.parquet to copy (default: the one in --src)")
a = ap.parse_args()
import anndata as ad

src, out = Path(a.src), Path(a.out)
lines = lambda f: [l.rstrip("\n").split("\t")[0] for l in gzip.open(f, "rt") if l.strip()]
genes, ids = lines(src / "features.tsv.gz"), lines(src / "barcodes.tsv.gz")
m = mmread(src / "matrix.mtx.gz").tocsr()                       # genes x cells
X = sp.csr_matrix(m.T, dtype=np.float32)                         # cells x genes
cells = pd.read_parquet(src / "cells.parquet").set_index("cell_id", drop=False).loc[ids]
obs = cells.drop(columns=["cell_id", "x", "y"])
obs.index = pd.Index(ids, name="cell_id")
for c in obs.columns:                                            # text columns as categories (smaller, and how scanpy keeps them)
    if obs[c].dtype == object or pd.api.types.is_string_dtype(obs[c].dtype):
        obs[c] = obs[c].astype(str).astype("category")
adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
adata.obsm["spatial"] = cells[["x", "y"]].to_numpy(np.float64)
adata.uns["source"] = "WebbScope sample_project (proof-of-concept data): raw counts and cell metadata"
out.mkdir(parents=True, exist_ok=True)
f = out / "sample.h5ad"
adata.write_h5ad(f, compression="gzip")
bnd = Path(a.boundaries) if a.boundaries else src / "cell_boundaries.parquet"
if bnd.resolve() != (out / "cell_boundaries.parquet").resolve():
    shutil.copyfile(bnd, out / "cell_boundaries.parquet")
print(f"wrote {f} ({f.stat().st_size / 1e6:.1f} MB): {adata.n_obs:,} cells x {adata.n_vars:,} genes, {X.nnz:,} non-zero counts, {obs.shape[1]} metadata columns")
