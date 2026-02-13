suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

process_dirichlet_fit <- function(fit, col_var) {
  stopifnot(!missing(fit), !missing(col_var), is.character(col_var), length(col_var) == 1)
  
  u <- summary(fit)
  
  coef_mat <- as.data.frame(u$coef.mat) %>%
    tibble::rownames_to_column("term_raw")
  
  n_celltypes <- length(u$varnames)
  n_rows <- nrow(coef_mat)
  
  if (n_celltypes == 0 || n_rows %% n_celltypes != 0) {
    stop("Unexpected coefficient matrix layout: cannot map coefficients to cell types.")
  }
  
  rows_per_celltype <- n_rows / n_celltypes
  
  coef_mat %>%
    mutate(
      celltype = rep(u$varnames, each = rows_per_celltype),
      term = str_replace_all(term_raw, fixed("get(col_var)"), col_var),
      model = str_replace_all(as.character(fit$formula)[3], fixed("get(col_var)"), col_var)
    ) %>%
    rename(pval = `Pr(>|z|)`) %>%
    mutate(
      label = case_when(
        is.na(pval)    ~ "",
        pval <= 0.001  ~ "***",
        pval <= 0.01   ~ "**",
        pval <= 0.05   ~ "*",
        TRUE           ~ ""
      )
    ) %>%
    select(celltype, term, Estimate, StdErr, z.value, pval, label, model)
}
