# QC report: what each plot and metric means

WebbScope makes a one-page QC report for each tissue (toolbar → **QC report**, or `/qc/<tissue>`), and an all-tissue summary (`/qc-summary`). This page explains how every number and plot is calculated. The code is `webbscope/qc.py`.

**Read this first**

- Everything is computed from what the project stores: the count matrix, the cell table, the cell boundaries and, if you imported them, the transcript coordinates. Nothing is read from the original Xenium output at report time.
- Flags compare a tissue **with the other tissues of the same project**. They are not a pass or fail against a 10x specification.
- "Transcripts" in this report mean **panel-gene** transcripts unless a section says otherwise. Control probes are not included in the gene counts.
- Cells that were left out at import (switched-off pieces, or cells below the minimum transcript count) are not in the project, so they are not in the report either.
- Some panels need transcripts or columns that a project may not have. The panel then says so instead of drawing. Projects imported before a feature existed need a re-import for that feature (see the last section).

## Contents

1. [Header and cards](#1-header-and-cards)
2. [Map of cells](#2-map-of-cells)
3. [Histograms: transcripts, genes, cell area](#3-histograms)
4. [Segmentation and nuclei](#4-segmentation-and-nuclei)
5. [Transcripts vs cell area](#5-transcripts-vs-cell-area)
6. [Genes ranked by counts](#6-genes-ranked-by-counts)
7. [Detection vs expression](#7-detection-vs-expression)
8. [Negative controls (cell counts)](#8-negative-controls-per-cell-counts)
9. [Transcript quality (QV)](#9-transcript-quality-qv)
10. [Transcript assignment](#10-transcript-assignment)
11. [Transcript density map](#11-transcript-density-map)
12. [Negative-control density map](#12-negative-control-density-map)
13. [This tissue against the other tissues, and the flags](#13-this-tissue-against-the-other-tissues-and-the-flags)
14. [Footer](#14-footer)
15. [All-tissue summary](#15-all-tissue-summary-and-qc_summarycsv)
16. [What needs a re-import](#16-what-needs-a-re-import)

## 1. Header and cards

The header shows the project name, run, slide, panel and software version when the import recorded them, and the number of panel genes.

| Card | How it is calculated |
|---|---|
| **cells** | Cells of the tissue in the project. |
| **area (approx.)** | The tissue is divided into 25 µm × 25 µm squares. The area is the number of squares that hold at least one cell, times 625 µm², in mm². This is a rough area: holes smaller than 25 µm count as tissue, and empty squares inside the tissue are not counted. |
| **cells per mm²** | Cells ÷ the area above. |
| **median transcripts / cell** | Median over cells of the **panel-gene count per cell**: the sum of the cell's counts over all genes in the matrix. Control probes and codewords are not in the matrix columns that are summed. |
| **median genes / cell** | Median over cells of the number of genes with a count of at least 1. |
| **transcripts per mm²** (last card) | Stored transcripts of the tissue ÷ the area. Only with transcripts imported. Without them this card shows the **share of cells with fewer than 20 transcripts** instead. |

## 2. Map of cells

Cells are drawn from the centroid positions, y pointing down like the Xenium images, on a grid of about 760 pixels along the longer side of the tissue.

- Colour = the mean of log10(transcripts per cell) of the cells in that pixel (a count of 0 is treated as 1).
- The colour range runs from the 1st to the 99th percentile of the pixel values, so a few extreme pixels do not wash out the map.
- The colour bar has five labels at equal steps of log10, so read it as a log scale.
- If a project has no transcript-count column the map uses genes per cell instead and says so on the bar.
- When there are few pixels per cell, each occupied pixel is widened by one pixel so isolated cells stay visible.

## 3. Histograms

| Plot | What is drawn |
|---|---|
| **Transcripts per cell** | Histogram of the panel-gene count per cell, on a log10 axis (cells with 0 are left out of this plot). |
| **Genes detected per cell** | Histogram of the number of genes with a count of at least 1, linear axis. |
| **Cell area** | Histogram of the cell-area column from the segmentation (µm²). Needs a column named `cell_area`, `area_um` or `area`; otherwise the panel says it is missing. |

For all three: 40 bins between the 0.2th and the 99.8th percentile. Values outside that range are placed in the end bins, so a few extreme cells do not stretch the axis. The dashed red line marks the **median**, with its value.

## 4. Segmentation and nuclei

Horizontal bars, each a percentage of the tissue's cells:

- one bar for each **segmentation method** recorded for the cells (the six most common), taken from the `segmentation_method` column, with the words "Segmented by" removed;
- **cells with a nucleus**: the share of cells whose `nucleus_area` is greater than 0.

Panels without those columns say so.

## 5. Transcripts vs cell area

A 2-D density plot, cell area (x) against transcripts per cell (y), both on log10 axes.

- 36 × 36 bins, with axes cut at the 0.5th and 99.5th percentile.
- Bin colour = log(1 + cells in the bin), scaled to the fullest bin. Empty bins are blank.
- Cells with a zero in either value are left out, and the plot needs at least 10 cells.

A healthy tissue shows transcripts rising with area. A cloud of large cells with few transcripts can point to merged or empty cells.

## 6. Genes ranked by counts

For each gene, the total count over the tissue's cells. Genes are sorted from most to least counts. The plot shows rank (x) against total counts (y, log10). Genes with 0 counts cannot be drawn on a log axis, so the plot states how many of the genes have no counts.

## 7. Detection vs expression

One dot per gene:

- **x** = mean counts per cell (total counts ÷ number of cells), log10 axis;
- **y** = percentage of the tissue's cells in which the gene has a count of at least 1.

Genes with a mean of 0 are left out. Genes follow a rising curve: the more a gene is expressed on average, the more cells detect it. Dots well below the curve are genes whose counts are concentrated in few cells (cell-type-specific genes or bursty expression); dots on the curve are broadly detected genes. Very low-expressed genes sit at the left where detection is close to the mean count, because most cells carry 0 or 1 counts.

## 8. Negative controls (per-cell counts)

This box appears when the import recorded non-gene features. It uses the **count matrix**, per cell, and is separate from the transcript-position map in section 12.

The table has one row per feature type found in the cell table: negative control probes, negative control codewords, unassigned codewords and deprecated codewords. For each type:

| Column | Meaning |
|---|---|
| **n** | Number of features of that type in the run. |
| **mean/cell** | Mean count per cell of that type. |
| **% ≥ 1** | Share of the tissue's cells with at least one count of that type. |

The small bar chart **Neg. controls per cell** is the share of cells with 0, 1, 2, or 3 or more negative-control counts, probes plus codewords together.

**Background level** = (mean negative-control counts per cell ÷ number of negative-control features) ÷ (mean panel-gene counts per cell ÷ number of panel genes) × 100.
In words: the average counts that one control feature collects, as a percentage of what an average gene collects. A value near 0 is clean. It is an estimate of background, because a control feature has no real transcript to detect, so counts on it are false positives.

## 9. Transcript quality (QV)

Needs transcripts imported with their quality values. The plot is a histogram of the QV of the **panel-gene transcripts assigned to a cell of this tissue**, from QV 0 to 40, one bar per integer (QV is rounded down, and values above 40 go in the last bar). The dashed line is QV 20 (Q20). Bars at 20 and above are blue and below are grey.

The note gives **% at Q20 or better**: transcripts with QV ≥ 20 ÷ all transcripts in the histogram.

This histogram counts transcripts of **any** quality. The transcripts actually stored and used in every other transcript-based number are only those with QV ≥ 20.

## 10. Transcript assignment

Among panel-gene transcripts with QV ≥ 20 that belong to this tissue, the percentage assigned to a cell.

- A transcript belongs to the tissue when it is assigned to one of its cells, or it is **unassigned**, in which case it goes to the tissue of its nearest cell at any distance. (Projects imported before this changed only kept unassigned ones within 25 µm; the report says how many are within 25 µm and how many farther.)
- **Tissue-level percentage** = assigned ÷ (assigned + unassigned).
- The second, smaller number is the same ratio for the **whole run**: all QV ≥ 20 panel-gene transcripts in the file, including those outside every tissue.
- The top card **transcripts in this tissue (all qualities)** counts the rows of `transcripts.parquet` assigned to this tissue's cells, whatever their quality or feature (controls included). It is recorded at import (re-import older projects to get it). The total of the whole run's file is shown in the viewer's *Current view* as **Transcripts in this slide**.

The tissue-level number can differ from the slide-wide figure, because each unassigned transcript is counted in the tissue of its nearest cell. Compare tissues with each other using the same number, not with the run figure.

## 11. Transcript density map

Needs transcripts. It shows **all stored transcripts** of the tissue: panel genes with QV ≥ 20, assigned to a cell of the tissue, or unassigned and given to the tissue of the nearest cell at any distance (older imports kept only those within 25 µm). Control probes and low-quality transcripts are not in this map.

- Transcripts are counted in **20 µm × 20 µm bins** and divided by the bin area, giving transcripts per mm².
- The bin size is the same for every tissue, so maps can be compared between tissues.
- Colour = log(1 + density), from 0 up to the **99th percentile of the occupied bins**. The colour bar has five labels. Empty bins are grey.
- A warning appears if fewer transcripts could be read than are recorded for the tissue.

Without transcripts, this slot shows the **largest 10x clusters** instead (top eight, as % of the tissue's cells).

## Expression vs quality (QV)

Needs transcripts imported by a version that records the average QV per gene (older projects show a note asking for a re-import). One dot per gene.

- **y** = the gene's total counts in the tissue (log scale), the same totals as *Genes ranked by counts*.
- **x** = the average QV of all of that gene's transcripts assigned to a cell of this tissue, at any quality (QV 0 to 40), recorded at import.
- The **dashed red line** is at QV 20, the quality threshold the stored transcripts use. A gene left of it has transcripts that are on average below Q20.
- Genes without counts or without transcripts are not drawn (the plot says how many are shown). Hover a dot to see the gene's name, average QV and total counts.

## 12. Negative-control density map

Needs transcripts imported by a version that stores control positions. It is drawn like section 11, but only for **negative-control transcripts**.

- A transcript is a negative control when its `codeword_category` is `genomic_control_probe`, `deprecated_codeword`, `negative_control_probe`, `negative_control_codeword` or `unassigned_codeword`. `predesigned_gene` and `custom_gene` are the positive (gene) probes. If the file has no `codeword_category` column, the feature-name prefixes are used instead (`GenomicControlProbe`, `DeprecatedCodeword`, `NegControlProbe`, `NegControlCodeword`, `UnassignedCodeword`).
- The same filters as the gene transcripts apply: QV ≥ 20, and assigned to a cell of the tissue or unassigned within 25 µm of one.
- **Colour = the number of negative-control transcripts in a bin** (1, 2, 3 …, linear; the last label means that many or more, up to 8+). A **grey** bin is a tissue bin (a bin that holds gene transcripts) with none; the pale background is outside the tissue. The unit is per bin, not per mm², so a single transcript is not shown as thousands per mm². The per-mm² figure for the whole tissue is in the note.
- The note gives the count per class and the overall density (count ÷ area, per mm²).

Negative controls should be rare and spread evenly. A patch or an edge with a much higher density than the rest points to a local problem such as debris, a bubble or an edge effect.

This map counts **transcript positions** (section 12). The table in section 8 counts per cell from the **matrix**. They are related but not the same numbers: the matrix has four feature types and the transcript file has five classes.

## 13. This tissue against the other tissues, and the flags

Two parts show the same comparison.

**Strip plot.** One row per metric. Grey dots are the other tissues, the red dot is this tissue, and the numbers at the ends are the smallest and largest values in the project. A row is drawn only if every tissue in the project has that metric.

**Flags table.** For each metric the report gives the value, the **percentile** and a **flag**:

- **Percentile** is the share of tissues with a lower value, with ties counted as half, so a tie reads 50. It is shown whenever the project has more than one tissue.
- **Robust z-score** = (value − median) ÷ (1.4826 × median absolute deviation), using all tissues of the project, including this one. The 1.4826 factor makes the spread comparable to a standard deviation for normally distributed data.
- **Flag** is `high` when z > 3, `low` when z < −3, `ok` otherwise, and `n/a` when the project has fewer than **6** tissues.

The metrics that can be flagged:

| Metric | Meaning | Worse when |
|---|---|---|
| median transcripts / cell | section 1 | low |
| median genes / cell | section 1 | low |
| cells per mm² | section 1 | low |
| cells with < 20 transcripts | share of cells (0 to 1) with fewer than 20 panel-gene transcripts | high |
| median cell area (µm²) | median of the cell-area column | low (over-segmentation). Unusually large cells are not flagged. |
| background vs genes (%) | the background level of section 8 | high |
| assigned to cells (%) | section 10, tissue level | low |
| at Q20 or better (%) | section 9 | low |

The flag says **high** or **low**, not good or bad: whether it is bad depends on the metric (for example, a high background value is bad and a high median transcripts value is not). A metric missing from the project (no cell area, no transcripts) is not in the table.

Below the table the report lists **cell integrity**: how many cells have no boundary or a degenerate one (fewer than 3 vertices), and the minimum transcript count that was applied at import.

## 14. Footer

The footer lists which blocks are missing and why, the **top 10 genes** by total counts in the tissue with their counts, how the area was estimated, and when the report was made.

## 15. All-tissue summary and `qc_summary.csv`

`/qc-summary` (button in the report) and `qc_summary.html` (written by **Export all tissues…**) show one strip plot per metric with **every tissue as a dot**.

- A panel is drawn only when at least 3 tissues have a value.
- The metrics are the eight in section 13 plus three more: **transcripts per mm²** (section 1), **neg. controls per mm²** (negative-control transcripts of section 12 ÷ area, worse when high) and **transcripts stored** (Q20+ panel-gene transcripts of the tissue). These three are not part of the per-tissue flags.
- The axis is log10 when the largest value is more than 30 times the smallest. The three transcript-based metrics always use a log axis.
- The dashed red line is the **1st percentile on the bad side**: the lowest 1% for metrics where higher is better, the highest 1% for metrics where lower is better. Tissues beyond the line are red and named, and listed in a table under the plot.
- With fewer than 100 tissues the line can only catch the single most extreme tissue, since 1% of fewer than 100 is less than one tissue. Read it as "the most extreme tissue", not as a threshold.
- Each panel has a plain-language description and the direction that counts as worse. A final table lists the tissues that fall beyond the line for any metric, with how many metrics and which.
- Names are shortened by removing the words that all tissue names share at the beginning. The full names show on screen and when hovering.

**Export all tissues…** also writes `qc_summary.csv`: one row per tissue with the columns `n_cells`, `area_mm2`, `density`, `transcripts_total`, `tx_density`, `nc_density`, `median_tx`, `median_genes`, `low_frac`, `median_area`, `no_boundary_frac`, `neg_ratio`, `assigned_frac`, `q20_frac`, and a `flag_…` column (`ok`, `high`, `low` or `n/a`) for each of the eight flagged metrics. A value that does not exist for a tissue is left empty.

## 16. What needs a re-import

Some panels depend on what was stored at import time. Projects imported before a feature was added show a note, not a chart.

| Panel | Needs |
|---|---|
| Section 8 and the background level | non-gene feature counts stored as cell columns when the run was imported |
| Sections 9, 10, 11, 12, the transcript cards | transcripts ticked at import |
| Section 12 only | positions of negative controls, which are stored only by imports made after that feature was added |
| Cell area, segmentation, nuclei | the matching columns in the cell table |

To add one, import the run again (with **Include the transcripts** for the transcript panels).

## Limits

- The area is an approximation from 25 µm squares. The cells-per-mm² and per-mm² transcript numbers inherit that error, mostly in small or fragmented tissues.
- The flags use the project's own tissues as the reference. In a project with few tissues or tissues that are all alike or all bad, a flag says little.
- The 1st-percentile line in the summary is a ranking tool, not a pass/fail limit.
- These panels were checked on one real run and on synthetic runs. They have not been checked on large real cohorts.
