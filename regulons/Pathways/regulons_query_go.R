# ---- packages ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(stringr)
})

csv_path <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_enriched_scflow/pySCENIC/out/wholeMicro/Micro.regulons.filtered.csv"

lines <- readLines(csv_path, warn = FALSE)

# Find the two header lines
i_names <- which(str_detect(lines, "^,,AUC,NES,"))          # line with AUC/NES/...
i_head  <- which(str_detect(lines, "^TF,MotifID,"))         # line with TF,MotifID,...

if (length(i_names) == 0 || length(i_head) == 0) {
  stop("Couldn't find the expected header lines: ',,AUC,NES,...' and/or 'TF,MotifID,...'")
}
i_names <- i_names[1]
i_head  <- i_head[1]

# Parse header fields
fields_names <- str_split(lines[i_names], ",", simplify = TRUE)
fields_head  <- str_split(lines[i_head],  ",", simplify = TRUE)

# Build final column names: first 2 from TF/MotifID, rest from AUC..RankAtMax
col_names <- c(fields_head[1, 1:2], fields_names[1, 3:ncol(fields_names)])

# Read data starting AFTER the TF,MotifID header line
reg <- read_csv(
  csv_path,
  skip = i_head,              # skip through the header line itself
  col_names = col_names,
  show_col_types = FALSE
)

colnames(reg)


# sanity check expected columns exist
stopifnot(all(c("TF", "TargetGenes") %in% colnames(reg)))

# ---- 2) parse TargetGenes -> vector of gene symbols (safe regex) ----
# TargetGenes looks like: "[('GENE1', 1.23), ('GENE2', 4.56), ...]"
extract_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unique(str_match_all(x, "'([^']+)'")[[1]][, 2]) %>% na.omit()
}

reg2 <- reg %>%
  mutate(target_vec = map(TargetGenes, extract_genes))

# ---- 3) merge rows with same TF (union of all targets) ----
tf_targets <- reg2 %>%
  group_by(TF) %>%
  summarise(
    targets = list(sort(unique(unlist(target_vec)))),
    n_targets = length(targets[[1]]),
    .groups = "drop"
  )

# optional: a “flat” two-column table (comma-separated targets)
tf_targets_flat <- tf_targets %>%
  mutate(targets_csv = map_chr(targets, ~ paste(.x, collapse = ","))) %>%
  dplyr::select(TF, n_targets, targets_csv)

# write it if you want
write_csv(tf_targets_flat, "/rds/general/user/bavot/home/TREM2/Micro/regulons/TF_targets_merged.csv")

# ---- 4) GO enrichment per TF for BP, MF, CC ----
background_symbols <- sort(unique(unlist(tf_targets$targets)))

sym2entrez <- function(sym) {
  out <- bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  unique(out$ENTREZID)
}

bg_entrez <- sym2entrez(background_symbols)

run_enrichGO <- function(sym_vec, universe_entrez, ont) {
  gene_entrez <- sym2entrez(sym_vec)
  if (length(gene_entrez) < 5) return(NULL)
  
  enrichGO(
    gene          = gene_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,          # "BP", "MF", "CC"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    minGSSize     = 10,
    readable      = TRUE
  )
}

# run BP/MF/CC for each TF and collect results
onts <- c("BP", "MF", "CC")

go_results_all <- pmap_dfr(
  list(tf_targets$TF, tf_targets$targets),
  function(tf, targets_sym) {
    map_dfr(onts, function(ont) {
      ego <- run_enrichGO(targets_sym, bg_entrez, ont)
      if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(tibble())
      
      as.data.frame(ego) %>%
        as_tibble() %>%
        mutate(TF = tf, Ontology = ont) %>%
        arrange(p.adjust) %>%
        slice_head(n = 10) %>%   # top 10 per TF per ontology
        dplyr::select(TF, Ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, geneID)
    })
  }
)

write_csv(go_results_all, "/rds/general/user/bavot/home/TREM2/Micro/regulons/TF_GO_enrichment_BP_MF_CC_top10.csv")
