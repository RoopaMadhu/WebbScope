# WebbScope

A browser-based viewer for spatial-transcriptomics data, built for 10x **Xenium** data. The inputs can be a counts matrix and metadata or raw Xenium output.

## What it does

- Visualizes cells in space colored by cell type, cell state, gene expression, gene signatures, and regions.
- Can view H&E and IF overlay or side-by-side.
- Calculates spatially smoothed expression of a gene or gene signature.
- Determines the neighbors of a cell type or region of interest and calculates the cell type or cell state composition of the neighborhood. *No stats*
- Visualize the neighborhoods
- In the Gallery view, visualize all the tissues at a glance.

## Install

Pick one. All three run the same code.

### A. macOS app

1. Download `WebbScope-mac.zip` and unzip it.
2. Right-click `WebbScope.app` → **Open**. The app is unsigned, so macOS asks once.
3. Click **Browse…** and choose a project folder.

To build the app yourself, see [Build from source](#build-from-source).

### B. Linux or HPC cluster (server only)

A cluster has no screen, so WebbScope serves the page and you view it in a browser on your own computer through an SSH tunnel. PyInstaller cannot cross-compile, so build the binary on the cluster:

```bash
# on your computer
./pack_source.sh                      # makes webbscope-src.tar.gz
scp webbscope-src.tar.gz <cluster>:~/

# on the cluster (a compute node is best)
tar xzf webbscope-src.tar.gz
module load python                    # Python 3.10 or newer; module names vary
bash build_linux.sh                   # needs internet or a pip mirror once
```

This produces `dist/WebbScope/` and `dist/WebbScope-linux-x86_64.tar.gz`. The build needs `binutils` (`objdump`). Build on the oldest machine you will run on.

Run it inside a job:

```bash
srun --pty -c 8 --mem 32G bash
WEBBSCOPE_CACHE=/path/to/scratch/webbscope_cache \
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
| **2. Project folder** (counts + cell table, for data that did not come from Xenium output) | A folder with `project.json`, `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` (or `matrix.npz`), `cells.parquet` and `cell_boundaries.parquet`; optional embeddings, regions, signatures. See [Project folder format](#project-folder-format). | **No.** |

### Auxiliary information

Optional additions to the main input.

| Input | Formats | Batch import |
|---|---|---|
| Transcripts | `transcripts.parquet` inside the raw Xenium output. | **Yes** |
| IF (morphology) image | The run's `morphology_focus` files (`.ome.tif`, one per channel or one file), read from the folder, the bundle or a path you give; not copied. | **Yes**, read from each run |
| Registered H&E | A folder of `.ome.tif` / `.ome.tiff`, one per tissue, already in the Xenium frame; optional `index.json` mapping image names to sample ids. | **Yes.** Images are matched to `<run>_tissue_N` by file name. |
| H&E slide scan (not registered) | `.czi`, `.ome.tif` or `.ome.tiff`; optional ImageScope `.annotations` files. | **Yes.** A folder of scans, matched to the runs by the numbers in the names. |
| ROI object (ROI heterogeneity) | A `region_heterogeneity` object: `.rds` (or `.qs2`, converted through R). | **No.** |
| Pathologist regions, signatures | `regions.geojson` and `signatures.json` in a project folder. | **No.** Project folder only. |

## Quick start

1. Start the app (see above). Either click **Browse…** and choose a project folder, or import a Xenium output from the links on the first page.
2. For a project folder, the app asks which columns of the cell table are the **sample id** and the **donor id**. It suggests the columns whose names contain "sample" and "donor"; correct them if needed and press **Open** again. Projects this app wrote itself are opened without the question.
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
| `tissue` | tissue within the run (`tissue_1`, …) |
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



## H&E next to the cells

There are two separate H&E features, with two kinds of input. **Raw, unregistered scans** (one whole-slide image, `.czi` or `.ome.tif` / `.ome.tiff`) are given in the H&E field of the import pages and open in the *H&E slide* window. **Registered images** (one OME-TIFF per tissue, already in the Xenium frame, `.ome.tif` or `.ome.tiff`) have their own fields on the open-project page, the single-run import page and the batch page: a folder of images and an optional `index.json` that maps image names to sample ids. They show under or beside the cells (the *Images* button). They can also be added with `python -m webbscope.he_registered`, or as a ready-made `he/` folder in the project. Neither is produced by the Xenium import unless you add the slide scan.

## Command line

```
WebbScope [project] [--port N] [--host ADDR] [--browser] [--headless] [--window] [--force] [--version]
```

- `--headless` serves only, with no window and no browser. This is the default on Linux.
- `--host` defaults to `127.0.0.1`, so the app is reachable only from that machine or through an SSH tunnel.
- `--force` redoes the import even if the cache is current.

## Project folder format

A project is a folder with `project.json`, a count matrix (`matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz`, or `matrix.npz`), `cells.parquet` and `cell_boundaries.parquet`. Embeddings, regions, signatures, transcripts and H&E tiles are optional. The full table is in [docs/NOTES.md](docs/NOTES.md#project-folder-format).

## Build from source

| Target | Command |
|---|---|
| macOS app | `./build_app.sh` (needs a `.venv` with `pyinstaller` and `pywebview`) → `dist/WebbScope.app` |
| Linux binary | `bash build_linux.sh` on a Linux machine → `dist/WebbScope-linux-x86_64.tar.gz` |
| Source tarball | `./pack_source.sh` → `webbscope-src.tar.gz` |


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
