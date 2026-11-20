# Export the cell outlines of an R object as cell_boundaries.parquet for WebbScope: one row per cell, the cell id and its outline
# (GeoParquet: the outline as WKB in a "geometry" column with the "geo" metadata; or, with --wkt, as WKT text).
#
#   Rscript tools/export_cell_boundaries.R OBJECT OUT.parquet [--element metadata] [--id cell_id] [--cells CELLS] [--wkt]
#
# OBJECT    .qs / .qs2 (read with qs2::qs_read), .rds, or a data frame saved in either; a list is read at --element (default metadata).
#           The table needs a cell-id column (--id, default cell_id) and an sf geometry column named geometry, shape or geom.
# --cells   keep only these cells: a .parquet with a cell_id column (e.g. the project's cells.parquet) or a text file, one id per line.
# Needs the R packages sf and arrow (and qs2 for .qs / .qs2). The coordinates are written as they are (WebbScope expects µm).
# Cells whose outline is empty are written with an empty polygon; WebbScope draws them as dots.
suppressMessages({library(sf); library(arrow)})
a <- commandArgs(TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
pos <- a[!startsWith(a, "--") & !(seq_along(a) - 1) %in% which(a %in% c("--element", "--id", "--cells"))]
if (length(pos) < 2) stop("usage: Rscript tools/export_cell_boundaries.R OBJECT OUT.parquet [--element metadata] [--id cell_id] [--cells CELLS] [--wkt]")
src <- pos[1]; out <- pos[2]
element <- opt("--element", "metadata"); idcol <- opt("--id", "cell_id"); cells <- opt("--cells"); wkt <- "--wkt" %in% a

x <- if (grepl("\\.qs2?$", src, ignore.case = TRUE)) qs2::qs_read(src, validate_checksum = FALSE) else readRDS(src)
m <- if (is.data.frame(x)) x else x[[element]]
if (is.null(m) || !is.data.frame(m)) stop("no table at element '", element, "' (elements: ", paste(names(x), collapse = ", "), ")")
gcol <- intersect(c("geometry", "shape", "geom"), names(m))[1]
if (is.na(gcol)) stop("the table has no geometry, shape or geom column")
if (!idcol %in% names(m)) stop("the table has no column '", idcol, "' (use --id)")
geom <- m[[gcol]]
if (!inherits(geom, "sfc")) geom <- sf::st_as_sfc(geom)          # e.g. WKT text
ids <- as.character(m[[idcol]])
if (anyDuplicated(ids)) stop("the cell ids in '", idcol, "' are not unique")
if (!is.null(cells)) {
  keep <- if (grepl("\\.parquet$", cells)) as.character(arrow::read_parquet(cells, col_select = "cell_id")$cell_id) else readLines(cells)
  miss <- sum(!keep %in% ids)
  if (miss) message(miss, " of the ", length(keep), " cells in ", cells, " are not in the table")
  sel <- match(keep[keep %in% ids], ids); ids <- ids[sel]; geom <- geom[sel]
}
types <- table(as.character(sf::st_geometry_type(geom)))
if (wkt) {
  tab <- arrow::arrow_table(cell_id = ids, geometry = sf::st_as_text(geom, digits = 15))   # full precision (the default rounds to 7 digits)
} else {
  wkb <- lapply(sf::st_as_binary(geom), as.raw)
  tab <- arrow::arrow_table(cell_id = ids, geometry = arrow::Array$create(wkb, type = arrow::binary()))
  gt <- paste0('"', sub("^MULTIPOLYGON$", "MultiPolygon", sub("^POLYGON$", "Polygon", names(types))), '"', collapse = ",")
  tab$metadata$geo <- sprintf('{"version":"1.0.0","primary_column":"geometry","columns":{"geometry":{"encoding":"WKB","geometry_types":[%s],"crs":null}}}', gt)
}
arrow::write_parquet(tab, out)
cat(sprintf("wrote %s: %d cells (%s), %d with an empty outline, as %s\n", out, length(ids), paste(names(types), types, sep = " ", collapse = ", "),
            sum(sf::st_is_empty(geom)), if (wkt) "WKT text" else "GeoParquet (WKB)"))
