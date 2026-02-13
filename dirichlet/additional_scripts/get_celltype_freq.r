suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(cli)
})

tally_cells <- function(metadata,
                        unique_id_var = "manifest",
                        celltype_var = "cluster_celltype") {
  unique_id_sym <- rlang::ensym(unique_id_var)
  celltype_sym  <- rlang::ensym(celltype_var)
  
  mat <- metadata %>%
    dplyr::select(!!unique_id_sym, !!celltype_sym) %>%
    dplyr::mutate(
      !!unique_id_sym := as.character(!!unique_id_sym),
      !!celltype_sym  := as.character(!!celltype_sym)
    ) %>%
    dplyr::count(!!unique_id_sym, !!celltype_sym, name = "n") %>%
    tidyr::pivot_wider(
      names_from  = !!celltype_sym,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::arrange(!!unique_id_sym) %>%
    tibble::column_to_rownames(var = rlang::as_string(unique_id_sym)) %>%
    as.matrix()
  
  mat
}

get_celltype_freq <- function(metadata,
                              unique_id_var = "manifest",
                              celltype_var = "cluster_celltype",
                              verbose = TRUE) {
  counts_mat <- tally_cells(metadata, unique_id_var, celltype_var)
  prop_counts_mat <- prop.table(counts_mat, margin = 1)
  
  if (isTRUE(verbose)) {
    cli::cli_alert_info(
      "Cell frequencies calculated across {.val {ncol(counts_mat)}} cell types."
    )
  }
  
  list(
    counts_mat = counts_mat,
    prop_counts_mat = prop_counts_mat
  )
}
