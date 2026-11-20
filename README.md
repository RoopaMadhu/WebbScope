# WebbScope

A browser-based viewer for spatial-transcriptomics data. The inputs can be a counts matrix and metadata or raw output (.tar.gz or folders) from Xenium SpaceRanger.

## What it does

- Visualizes cells in space colored by cell type, cell state, gene expression, gene signatures, and regions.
- Can view H&E and IF overlay or side-by-side.
- Calculates spatially smoothed expression of a gene or gene signature.
- Determines the neighbors of a cell type or region of interest and calculates the cell type or cell state composition of the neighborhood. *No stats*
- Visualize the neighborhoods
- In the Gallery view, visualize all the tissues at a glance.

## Install


### A. macOS app

Download WebbScope-v0.1.0-mac.zip: The app is unsigned, so macOS throws a security warning.


### B. Linux or HPC cluster (server only)


Download WebbScope-v0.1.0-src.tar.gz to the HPC.

```bash
# on the cluster (a compute node is best)
tar xzf WebbScope-v0.1.0-src.tar.gz
cd WebbScope
bash build_linux.sh      # needs Python >=3.10
```

Run it inside a job:

```bash
srun --pty -c 8 --mem 70G bash
WEBBSCOPE_CACHE=/path/to/webbscope_cache \
  dist/WebbScope/WebbScope --port 8765 /path/to/project
```

It prints the `ssh -N -L …` command to run on your computer. Run it, then open `http://localhost:8765`.

Linux has no **Browse…** dialog, so type the project path.

### C. From source (any OS with Python 3.10+)

```bash
python -m venv .venv
.venv/bin/pip install -r requirements-linux.txt     # on macOS also: pip install pywebview
.venv/bin/python run_app.py /path/to/project_folder --browser
```

## Inputs

The main input is **one** of the two rows of the first table. *Batch import* says whether the same input can be given on the **Batch import** page (a folder of many runs). Everything can be given for a single run: raw output on the single-run import page, a project folder on the open-project page.

### Main input

| Input | Formats | Batch import |
|---|---|---|
| **1. Raw Xenium output** (cells, boundaries, counts, clusters) | An extracted folder, a `.tar.gz` / `.tar` bundle or a `.zip` bundle. Files read: `cell_feature_matrix/` (the `matrix.mtx.gz` form; `.h5` is not read), `cells.parquet`, `cell_boundaries.parquet`, `experiment.xenium` and the 10x analysis files. The Atera `.zarr.zip` layout is not supported. | **Yes.** A folder holding many runs; all folders or all bundles, not a mix. |
| **2. Project folder** (a matrix + cell table, for data that did not come from Xenium output) | A folder with `features.tsv.gz` + `barcodes.tsv.gz`, the cell metadata (`cells_metadata.parquet` or `cells.parquet`) and the matrix: `counts.mtx` or `matrix.mtx` (raw counts) and/or `counts_logcpx.mtx` (log-normalised values); each also as `.mtx.gz` or a sparse `.npz`. **Or** one AnnData `.h5ad` file: give its folder to also read the folder's cell outlines and other extras, or the file itself (**Browse .h5ad…**) to read that file alone. Optional: `project.json`, `cell_boundaries.parquet` (without it cells are drawn as dots), embeddings, regions, signatures. See [Project folder format](#project-folder-format). | **No.** |

### Auxiliary information

Optional additions to the main input.

| Input | Formats | Batch import |
|---|---|---|
| Transcripts | `transcripts.parquet` inside the raw Xenium output. | **Yes** |
| IF (morphology) image | The run's `morphology_focus` files (`.ome.tif`, one per channel or one file), read from the folder, the bundle or a path you give; not copied. | **Yes**, read from each run |
| Registered H&E | A folder of `.ome.tif` / `.ome.tiff`, one per tissue, already in the Xenium frame; optional `index.json` mapping image names to sample ids. | **Yes.** Images are matched to `<run>_T<i>` by file name. |
| H&E slide scan (not registered) | `.czi`, `.ome.tif` or `.ome.tiff`; optional ImageScope `.annotations` files. | **Yes.** A folder of scans, matched to the runs by the numbers in the names. |
| IMC slide (imaging mass cytometry) | An `.mcd` file of the slide (all its ROIs and channels). Converted once into the cache; the file is never changed. See [IMC slide](#imc-slide-mcd). | **No.** Single-run import page, open-project page or the Images menu. |
| ROI object (ROI heterogeneity) | A `region_heterogeneity` object: `.rds` (or `.qs2`, converted through R). | **No.** |
| Pathologist regions, signatures | `regions.geojson` and `signatures.json` in a project folder. | **No.** Project folder only. |

## Quick start

1. Start the app (see above). Either click **Browse…** and choose a project folder, or import a Xenium output from the links on the first page.
2. For a project folder, the app asks which columns of the cell table are the **sample id** and the **donor id**. It suggests the columns whose names contain "sample" and "donor"; correct them if needed and press **Open** again. A folder without `project.json`, or an `.h5ad`, gets a longer form (see [Project folder format](#project-folder-format)). Projects this app wrote itself are opened without the question.
3. Pick **Colour by** in the left panel, type a gene, open **Slide**, **Gallery** or **QC report** from the toolbar.

### Import your own Xenium data

Pages reached from the open-project page:

| Page | Use it for |
|---|---|
| **Import a Xenium output bundle…** (`/xenium`) | One run: a `.tar.gz` or `.zip` bundle, an extracted folder, or the files separately. Shows a tissue preview where you can rename, merge, split or switch off pieces. |
| **Batch import** (`/xenium-batch`) | A folder of many runs (folders, `.tar.gz` or `.zip` bundles, not a mix). Each run is converted with automatic tissue detection and all tissues are merged into `<out>/_merged`, which then opens. Optionally a folder of H&E scans, matched to the runs by name (see below). |
| **Open existing merged project** | Reopen an earlier `_merged` folder without merging again. |

Tick **Include the transcripts** to enable the transcript overlay and the transcript-based QC panels (quality values, assigned share, negative-control density).

Tick **Write `transcript_tissue.parquet`** to get, in each run's folder, a table to join to that run's raw `transcripts.parquet`. It has one row per transcript that belongs to a tissue:

| Column | Meaning |
|---|---|
| `row` | 0-based row number in the raw `transcripts.parquet` (always present) |
| `transcript_id` | the raw file's `transcript_id`, when the file has that column |
| `tissue` | tissue within the run (`T0`, …) |
| `sample` | `<run folder>_<tissue>`, the sample name in the merged project |
| `assignment` | `cell` or `nearest_cell` |


Imports write a project folder (default `~/.webbscope/projects/<run>`). The viewer cache goes to `~/.webbscope/cache`. Set `WEBBSCOPE_CACHE` to move it. The input data is never modified.

## Using the viewer

| Toolbar | What it does |
|---|---|
| **Slide** | All tissues of one slide in labelled boxes. Click a box to open that tissue. The first selector entry, **All slides**, scrolls through every slide on a common scale. Save PNG gives one picture, or `slides.zip` for all slides. |
| **Gallery** | Every tissue as a tile; *Save PNG → Everything* writes numbered pages in a zip. |
| **QC report** | One-page HTML report per tissue (prints to one A4 page). *Export all tissues…* writes every report plus `qc_summary.csv` and `qc_summary.html`. How every plot and metric is calculated: [docs/QC.md](docs/QC.md). |
| **Figure…** | Publication PNG or PDF of the current view or the whole tissue, at 300 or 600 dpi.
| **Profile** | Timings and memory of the last project load and of each request. |

In the Mac app window, **QC report**, **Profile** and the validation-report link open in your default browser (the window itself cannot open a new tab); in a browser they open in a new tab.

### Gene signatures

Every signature picker lists the SLE ISG-amplification list, then **GO** and **Reactome**: click either to open its gene sets, or type in the search box to search all of them. The GO and Reactome sets come from [MSigDB](https://www.gsea-msigdb.org/) release 2026.2 (human gene symbols; GO biological process, cellular component and molecular function, and Reactome pathways), bundled in `webbscope/engine/data/genesets.json.gz` and rebuilt with `tools/build_genesets.py`. Only the sets with at least 5 genes on the project's panel are listed. MSigDB is distributed under [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) (Liberzon et al. 2011, *Bioinformatics*; Gene Ontology Consortium 2023, *Genetics*; Milacic et al. 2024, *Nucleic Acids Res*). The older hand-curated signatures are hidden while they are reviewed.



## H&E next to the cells

There are two separate H&E features, with two kinds of input. **Raw, unregistered scans** (one whole-slide image, `.czi` or `.ome.tif` / `.ome.tiff`) are given in the H&E field of the import pages and open in the *H&E slide* window. **Registered images** (one OME-TIFF per tissue, already in the Xenium frame, `.ome.tif` or `.ome.tiff`) have their own fields on the open-project page, the single-run import page and the batch page: a folder of images and an optional `index.json` that maps image names to sample ids. They show under or beside the cells (the *Images* button). They can also be added with `python -m webbscope.he_registered`, or as a ready-made `he/` folder in the project. Neither is produced by the Xenium import unless you add the slide scan.

## IMC slide (.mcd)

An imaging mass cytometry file (for example the IMC run on the same slide after Xenium) can be shown beside the cells. It is **not aligned** to the Xenium coordinates, so, like the H&E slide scan, it only goes in the pane beside the cells (with its own pan and zoom, rotate and flip) and shows the **raw counts**.

1. Give the `.mcd` in the **IMC slide** field of the open-project page or of the single-run Xenium import page (not the batch page), or in the viewer in the **Images** menu (*Add IMC file*). It is converted once, in the background, into the cache; the file is never changed.
2. In the **Images** menu tick **IMC** under *Beside*. Each ROI (acquisition) of the file is a separate image, so pick the ROI that belongs to the tissue in the list next to the orientation buttons (the choice is remembered). Incomplete acquisitions are listed and marked.
3. **IMC channels** opens a list of all channels with a search box. Tick one or several (up to 16): each has its own colour and brightness and the ticked channels are added up as light.

The Gallery, Regions, Glomeruli and similar tile windows show the chosen ROI beside each tissue; the Slide window does not offer IMC. Details, limits and what was checked: [docs/NOTES.md](docs/NOTES.md#imc-slide-mcd).

## Command line

```
WebbScope [project] [--port N] [--host ADDR] [--browser] [--headless] [--window] [--force] [--version]
```

- `--headless` serves only, with no window and no browser. This is the default on Linux.
- `--host` defaults to `127.0.0.1`, so the app is reachable only from that machine or through an SSH tunnel.
- `--force` redoes the import even if the cache is current.

## Project folder format

A project is a folder with `features.tsv.gz`, `barcodes.tsv.gz`, the cell metadata and a matrix, or with one AnnData `.h5ad` file. The names are WebbScope's own; 10x's names are accepted too, so a `cell_feature_matrix` and `cells.parquet` taken from a Xenium output work as they are.

| File | Holds | Accepted names |
|---|---|---|
| raw counts | whole numbers | `counts.mtx` or `matrix.mtx` (10x's name) |
| log-normalised values | log1p of normalised counts | `counts_logcpx.mtx` |
| cell metadata | one row per cell: id, x, y, sample and any other columns | `cells_metadata.parquet` or `cells.parquet` (Xenium's name) |

Give the raw counts, the log-normalised values or both (see below). Each matrix can also be `.mtx.gz`, or `.npz` written by `scipy.sparse.save_npz` (sparse; compressed by default, the smallest and fastest to read; a dense array is refused); genes x cells or cells x genes, with the genes and cells of `features.tsv.gz` / `barcodes.tsv.gz`. Only one file of each kind may be in the folder (for example not both `counts.mtx` and `matrix.mtx`). `project.json`, `cell_boundaries.parquet`, embeddings, regions, signatures, transcripts and H&E tiles are optional. The full table is in [docs/NOTES.md](docs/NOTES.md#project-folder-format).

**Without `project.json`** (and for an `.h5ad`) the open page asks which column is which: cell id, x, y and sample (needed), donor, cell type, cell state, three region levels and total counts per cell (optional), suggested from the column names. With no cell type every cell is "unlabelled"; with no cell state the cell type is used. *Save these choices as project.json in the folder* (ticked by default) writes them into the folder, so they travel with the data to another computer; untick it (or if the folder is read-only) and they are kept in WebbScope's cache only. Either way the next opening goes straight in; **Change columns…** under the path asks again.

Example: `sample_project_h5ad/` is the sample project as an `.h5ad` (`sample.h5ad`: raw counts in `X`, the cell metadata in `obs`, coordinates in `obsm["spatial"]`) with its `cell_boundaries.parquet` next to it as GeoParquet (an `.h5ad` has no place for outlines; they are matched to the cells by cell id).

**Cell outlines** (`cell_boundaries.parquet`, optional) can be in either layout:
- one row per outline point: `cell_id`, `vertex_x`, `vertex_y` (as Xenium writes it; `sample_project/` uses it);
- one row per cell: a cell-id column and the outline in a `geometry`, `shape` or `geom` column, as GeoParquet (WKB; what geopandas and `sf` + `arrow` write) or as WKT text (`POLYGON ((x y, ...))`). The cell-id column is the one whose values match the project's cells. Of a multipolygon the largest part is drawn, holes are left out, and an empty or non-polygon outline leaves the cell as a dot; the log counts each.

`tools/export_cell_boundaries.R` writes the per-cell layout from an R object (`.qs`, `.qs2` or `.rds`) whose table (default the list's `metadata` element) has a cell id and an `sf` geometry column: `Rscript tools/export_cell_boundaries.R OBJECT cell_boundaries.parquet [--element metadata] [--id cell_id] [--cells cells.parquet] [--wkt]`. Rebuild it with `python tools/make_sample_h5ad.py` (needs `pip install anndata`).

**Folder or file.** Given the folder, the `.h5ad` is read together with the folder's other files (`cell_boundaries.parquet`, `regions.geojson`, `project.json` …), and the choices can be saved there as `project.json`. Given the `.h5ad` file itself, only that file is read (the cells are drawn as dots, since an `.h5ad` holds no outlines) and the choices are kept in WebbScope's cache. A folder with several `.h5ad` files can only be opened file by file.

**`.h5ad`**: read with `h5py` (anndata itself is not needed; files from anndata 0.7 or later). Choose the matrix with the raw counts and/or the one with log-normalised values (`X`, a layer such as `layers/counts`, or `raw/X`; suggested from the values; when `raw/X` has more genes than `X`, the genes of `X` are used) and the gene names (the gene table's index, or a column such as `gene_symbols` when the index holds Ensembl ids). The cell table is `obs` plus its index; coordinates come from `obs` columns or an `obsm` array named like `spatial` (offered as `spatial_x` / `spatial_y`). `obsm["X_umap"]` and `obsm["X_pca"]` (or their `_harmony` versions) fill the Embeddings window.

**Raw counts, log-normalised values, or both.**
- *Raw counts only*: WebbScope computes log1p(counts / total counts of the cell × scale). The scale comes from `project.json`, else it can be typed when opening; empty = the median total counts per cell. (Xenium imports write 10,000 into their `project.json`.)
- *Log-normalised values only*: stored and shown as they are, never recalculated. There is no "raw counts" layer, single-gene violins have no raw-count histogram, and the QC report's gene plots sum the log-normalised values; transcripts per cell appear only if a total-counts column is chosen.
- *Both*: the log-normalised values are shown as they are; the counts give the "raw counts" layer, the raw-count histograms and the QC numbers. When the two have the same non-zero entries only the values are stored twice (about a third more cache).
- Scaled or z-scored values (with negative numbers) are refused.

## Build from source

| Target | Command |
|---|---|
| macOS app | `./build_app.sh` (needs a `.venv` with `pyinstaller` and `pywebview`) → `dist/WebbScope.app` |
| Linux binary | `bash build_linux.sh` on a Linux machine → `dist/WebbScope-v0.1.0-linux-x86_64.tar.gz` |
| Source tarball | `./pack_source.sh` → `WebbScope-v0.1.0-src.tar.gz` |


More detail, including what was checked and how, is in [docs/NOTES.md](docs/NOTES.md).

## Layout of this repository

```
webbscope/        Python server (importer, QC, transcripts) and the web front end (webbscope/web)
run_app.py        entry point
build_app.sh      macOS build
build_linux.sh    Linux build
pack_source.sh    source tarball for the cluster
docs/QC.md        how each QC plot and metric is calculated
docs/NOTES.md     detailed notes, formats and test record
```
