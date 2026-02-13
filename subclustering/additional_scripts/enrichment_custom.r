suppressPackageStartupMessages({
  library(dplyr)
  library(cli)
})

enrichment_custom <- function(genes,
                              reference,
                              genesets,
                              adj = "fdr",
                              verbose = FALSE) {
  stopifnot(is.character(genes), is.character(reference))
  stopifnot(is.list(genesets), length(genesets) > 0)
  
  genes <- unique(na.omit(genes))
  reference <- unique(na.omit(reference))
  
  if (length(genes) == 0) stop("`genes` is empty after filtering.")
  if (length(reference) == 0) stop("`reference` is empty after filtering.")
  
  # Use only genes present in the reference universe
  genes <- intersect(genes, reference)
  
  if (length(genes) == 0) {
    return(tibble(
      TermID = names(genesets),
      genes = 0L, all = 0L, geneset_size = lengths(genesets),
      pct_overlap = NA_real_, enrichment_ratio = NA_real_,
      pval = NA_real_, padj = NA_real_, overlap_genes = ""
    ))
  }
  
  # Split reference into "genes" vs "background" (excluding genes to avoid double counting)
  reference_bg <- setdiff(reference, genes)
  
  ra <- length(genes) / length(reference)  # expected fraction in a random draw
  
  term_ids <- names(genesets)
  if (is.null(term_ids)) term_ids <- paste0("Term", seq_along(genesets))
  
  res <- lapply(seq_along(genesets), function(i) {
    term <- unique(na.omit(as.character(genesets[[i]])))
    term <- intersect(term, reference)
    
    if (isTRUE(verbose)) {
      cli::cli_text("Processing term {i}/{length(genesets)}: {.emph {term_ids[[i]]}}")
    }
    
    geneset_size <- length(term)
    expect <- geneset_size * ra
    
    GinSet <- sum(genes %in% term)
    RinSet <- sum(reference_bg %in% term)
    
    GninSet <- length(genes) - GinSet
    RninSet <- length(reference_bg) - RinSet
    
    # Handle edge cases
    enrichment_ratio <- if (expect > 0) round(GinSet / expect, 2) else NA_real_
    overlap <- intersect(genes, term)
    overlap_genes <- paste(overlap, collapse = ",")
    
    fmat <- matrix(
      c(GinSet, RinSet,
        GninSet, RninSet),
      nrow = 2, byrow = FALSE,
      dimnames = list(c("genes", "reference"), c("inSet", "ninSet"))
    )
    
    fish <- fisher.test(fmat, alternative = "greater")
    pval <- as.numeric(fish$p.value)
    
    inSet_total <- RinSet + GinSet
    pct_overlap <- if (inSet_total > 0) round((GinSet / inSet_total) * 100, 2) else NA_real_
    
    tibble(
      TermID = term_ids[[i]],
      genes = GinSet,
      all = inSet_total,
      geneset_size = geneset_size,
      pct_overlap = pct_overlap,
      enrichment_ratio = enrichment_ratio,
      pval = pval,
      overlap_genes = overlap_genes
    )
  })
  
  out <- bind_rows(res) %>%
    arrange(pval) %>%
    mutate(
      padj = p.adjust(pval, method = adj),
      pval = as.numeric(format(pval, format = "e", digits = 2)),
      padj = as.numeric(format(padj, format = "e", digits = 2))
    ) %>%
    select(-overlap_genes, everything(), overlap_genes)
  
  out
}
