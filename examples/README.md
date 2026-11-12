# Examples: public 10x Xenium human lung

Everything here was made from the public **Xenium Preview Human Non-diseased Lung With Add-on FFPE** run (10x Genomics),
imported with WebbScope's single-run import (295,817 cells, 392 genes, 1 tissue; the 10x graph clusters stand in for cell types).
It is the screen layout of the viewer and a sample of what can be done with it, not a biological analysis.

## Open the HTML files (no app needed)

| File | What it is |
|---|---|
| [`lung/qc_report_tissue_1.html`](lung/qc_report_tissue_1.html) | The one-page QC report (toolbar → **QC report**). Self-contained; open it in any browser; prints to one A4 page. |
| [`lung/qc_summary.html`](lung/qc_summary.html) | The cohort page of *Export all tissues…*. With a single tissue it is short; it compares tissues when there are several. |

## Screenshots of the viewer (`lung/*.png`)

| File | Shows |
|---|---|
| `viewer_default.png` | Default layout: left panel (donor, colour by, display, filters), toolbar on top, cells in the middle, current view, legend and composition on the right. |
| `viewer_gene_EPCAM.png` | **Colour by → Gene expression**, with the layer, smoothing, z-score and clip controls. |
| `viewer_dual_EPCAM_CD68.png` | **Two genes (dual colour)**: EPCAM red, CD68 green, both yellow. |
| `viewer_zoom_cell_boundaries.png` | Zoomed in (scale bar 20 µm): the segmentation polygons are drawn; hover shows the cell. |
| `window_gallery.png` | **Gallery**: every tissue as a tile. |
| `window_embeddings.png` | **Embeddings**: the 10x UMAP in a floating window (not batch-corrected). |
| `window_neighbors.png` | **Neighbors**: cell-type composition within 25 µm of the cells. |
| `window_violins.png` | **Violins**: genes per cell type. Sparse here because there is one tissue. |

## What is not shown, and why

- **H&E, IF, IMC, pathologist regions, Slide window, multi-tissue Gallery and Violins**: this run has no H&E or regions and only one tissue.
- **Signatures**: the built-in signatures are mostly kidney genes. On this lung panel most have 0 genes, so the window is empty and was left out.
- Known problem found while making these: with the default signature chip, the **Violins** window shows
  `Error: supplied range of [nan, nan] is not finite` when that signature has 0 genes on the panel (removing the chip fixes it).

## Recreate

```bash
.venv/bin/python run_app.py ~/.webbscope/projects/<imported lung project> --browser
```
The QC pages come from `/qc/tissue_1` and `/qc-summary` on the running server.
The raw run is 17 GB and is not copied here (`Xenium_Preview_*/` is git-ignored).
