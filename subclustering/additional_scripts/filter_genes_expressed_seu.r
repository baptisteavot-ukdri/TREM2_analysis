suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(cli)
})

filter_seu_genes <- function(seu,
                             min_counts = 1L,
                             min_cells_pc = 0.1,
                             assay = "RNA",
                             slot = "counts",
                             verbose = TRUE) {
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(is.numeric(min_counts), length(min_counts) == 1, min_counts >= 0)
  stopifnot(is.numeric(min_cells_pc), length(min_cells_pc) == 1, min_cells_pc >= 0, min_cells_pc <= 1)
  
  n_genes_before <- nrow(seu)
  n_cells <- ncol(seu)
  
  min_cells <- max(1L, floor(n_cells * min_cells_pc))
  
  counts <- Seurat::GetAssayData(seu, assay = assay, slot = slot)
  
  # Keep genes expressed (>= min_counts) in at least min_cells cells
  keep_genes <- Matrix::rowSums(counts >= min_counts) >= min_cells
  
  seu_f <- seu[keep_genes, , drop = FALSE]
  n_genes_after <- nrow(seu_f)
  
  if (isTRUE(verbose)) {
    cli::cli_text(
      "Selected {n_genes_after} / {n_genes_before} genes with >= {min_counts} count(s) in >= {min_cells} cells ({round(min_cells_pc * 100, 2)}%)."
    )
  }
  
  seu_f
}
