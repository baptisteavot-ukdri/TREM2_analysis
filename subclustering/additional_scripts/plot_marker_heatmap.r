suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

scale_data <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s  <- stats::sd(x, na.rm = TRUE)
  
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  
  (x - mu) / s
}

plot_marker_heatmap <- function(seu,
                                markers,
                                cluster_col,
                                assay = "RNA",
                                slot = "data",
                                cell_id_col = "barcode",
                                show_row_names = FALSE,
                                show_column_dend = TRUE,
                                cluster_columns = TRUE,
                                column_names_rot = 45,
                                column_names_side = "top",
                                row_names_side = "left",
                                ht_width = 4,
                                fontsize = 8) {
  
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(is.character(cluster_col), length(cluster_col) == 1)
  
  meta <- as.data.frame(seu@meta.data)
  if (!cluster_col %in% colnames(meta)) stop("cluster_col not found in seu@meta.data: ", cluster_col)
  
  # Determine cell IDs to use for grouping
  if (!cell_id_col %in% colnames(meta)) {
    meta[[cell_id_col]] <- colnames(seu)
  }
  if (!all(meta[[cell_id_col]] %in% colnames(seu))) {
    stop("cell_id_col ('", cell_id_col, "') contains IDs not present in colnames(seu).")
  }
  
  if (!"gene" %in% colnames(markers)) stop("markers must contain a 'gene' column.")
  genes <- unique(markers$gene)
  genes <- genes[genes %in% rownames(seu)]
  if (length(genes) == 0) stop("None of markers$gene are present in rownames(seu).")
  
  expr_mat <- as.matrix(GetAssayData(seu, assay = assay, slot = slot))
  
  # Mean expression per cluster (genes x clusters)
  cell_groups <- split(meta[[cell_id_col]], meta[[cluster_col]])
  mean_expr <- vapply(
    cell_groups,
    FUN = function(cells) {
      Matrix::rowMeans(expr_mat[, as.character(cells), drop = FALSE])
    },
    FUN.VALUE = numeric(nrow(expr_mat))
  )
  rownames(mean_expr) <- rownames(expr_mat)
  
  mat <- mean_expr[genes, , drop = FALSE]
  mat <- t(apply(mat, 1, scale_data)) # genes x clusters (z-score per gene)
  
  # Column clustering (optional)
  col_clust <- if (isTRUE(cluster_columns)) stats::hclust(stats::dist(t(mat))) else FALSE
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Z-score",
    col = circlize::colorRamp2(
      breaks = c(min(mat, na.rm = TRUE), mean(mat, na.rm = TRUE), max(mat, na.rm = TRUE)),
      colors = c("navy", "white", "firebrick")
    ),
    column_names_side = column_names_side,
    column_names_rot = column_names_rot,
    cluster_rows = FALSE,
    cluster_columns = col_clust,
    show_row_dend = FALSE,
    show_column_dend = show_column_dend,
    show_row_names = show_row_names,
    row_names_side = row_names_side,
    width = grid::unit(ht_width, "cm"),
    row_names_gp = grid::gpar(fontsize = fontsize)
  )
  
  ht
}
