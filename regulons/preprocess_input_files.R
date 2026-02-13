suppressPackageStartupMessages({
  library(argparse)
  library(qs)
  library(arrow)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(cli)
})

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
parser <- ArgumentParser()

required <- parser$add_argument_group("Required arguments")
required$add_argument("--sce_file",     help = "Path to *_sce.qs", required = TRUE)
required$add_argument("--feather_file", help = "Path to pySCENIC ranking database (.feather)", required = TRUE)
required$add_argument("--celltype",     help = "Cell type label used for output filename", required = TRUE)
required$add_argument("--outdir",       help = "Output directory", required = TRUE)

args <- parser$parse_args()

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------
cli::cli_h2("Loading inputs")

sce <- qs::qread(args$sce_file)
if (!"gene" %in% colnames(rowData(sce))) {
  stop("rowData(sce)$gene not found. Please ensure rowData contains a 'gene' column.")
}
rownames(sce) <- rowData(sce)$gene

ranking_db <- arrow::read_feather(args$feather_file)
db_genes <- colnames(ranking_db)
if (is.null(db_genes) || length(db_genes) == 0) {
  stop("No gene columns found in the feather ranking database: ", args$feather_file)
}

# ------------------------------------------------------------------------------
# Filter genes expressed in >= 10% of cells (non-zero logcounts)
# ------------------------------------------------------------------------------
cli::cli_h2("Filtering genes")

if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
  stop("Assay 'logcounts' not found in SCE.")
}

logc <- SummarizedExperiment::assay(sce, "logcounts")
if (!inherits(logc, "Matrix")) logc <- as.matrix(logc)

min_cells_pc <- 0.10
min_cells <- ceiling(ncol(logc) * min_cells_pc)

keep <- Matrix::rowSums(logc != 0) >= min_cells
sce_f <- sce[keep, , drop = FALSE]

cli::cli_text(
  "Kept {.val {nrow(sce_f)}} / {.val {nrow(sce)}} genes expressed in >= {.val {min_cells}} cells ({.val {min_cells_pc * 100}}%)."
)

# ------------------------------------------------------------------------------
# Keep only overlap with pySCENIC ranking DB
# ------------------------------------------------------------------------------
cli::cli_h2("Intersecting with ranking database")

overlap <- intersect(rownames(sce_f), db_genes)
if (length(overlap) == 0) {
  stop("No overlap between SCE genes and ranking DB genes.")
}

sce_f <- sce_f[overlap, , drop = FALSE]
cli::cli_text("Overlap genes retained: {.val {nrow(sce_f)}}")

# ------------------------------------------------------------------------------
# Export logcounts matrix to TSV (genes x cells)
# ------------------------------------------------------------------------------
cli::cli_h2("Writing output")

out_path <- file.path(args$outdir, sprintf("%s.tsv", args$celltype))

mat <- SummarizedExperiment::assay(sce_f, "logcounts")
mat <- as.matrix(mat)

# Write with gene names as first column (more robust than rownames in TSV)
out_df <- data.frame(gene = rownames(mat), mat, check.names = FALSE)

write.table(
  out_df,
  file = out_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cli::cli_alert_success("Wrote: {out_path}")
