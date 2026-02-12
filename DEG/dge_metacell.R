suppressPackageStartupMessages({
  library(argparse)
  library(cli)
  library(dplyr)
  library(stringr)
  library(qs)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(edgeR)
  library(limma)
})

source("/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/pseudobulk_dge/generate_pseudobulk_limma_trend_deg_gazestani.r")

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
parser <- ArgumentParser()

required <- parser$add_argument_group("Required arguments")
optional <- parser$add_argument_group("Optional arguments")

required$add_argument("--sce", help = "Path to SCE .qs object", required = TRUE)
required$add_argument("--dependent_var", help = "Dependent variable (colData)", required = TRUE)
required$add_argument("--ref_class", help = "Reference class for dependent_var (or 'NULL')", required = TRUE)
required$add_argument("--confounding_vars", help = "Comma-separated confounders (colData)", required = TRUE)
optional$add_argument("--stratification_var", help = "Stratification variable (colData) or 'NULL'", default = "NULL")
required$add_argument("--output_dir", help = "Output directory", required = TRUE)
required$add_argument("--mod_name", help = "Model name (used in folder naming)", required = TRUE)

args <- parser$parse_args()

args$confounding_vars <- strsplit(args$confounding_vars, ",", fixed = TRUE)[[1]] %>%
  str_trim() %>%
  discard(~ .x == "")

celltype <- gsub("_sce\\.qs$", "", basename(args$sce))
mod_name <- args$mod_name

# ------------------------------------------------------------------------------
# Read input
# ------------------------------------------------------------------------------
cli::cli_text(c(
  "Reading {.strong {celltype}} from ",
  cli::col_green("{args$sce}")
))

sce <- qs::qread(args$sce)
rownames(sce) <- rowData(sce)$gene
sce$manifest <- droplevels(sce$manifest)

# Validate key columns early
stopifnot("manifest" %in% colnames(colData(sce)))

missing_cols <- setdiff(c(args$dependent_var, args$confounding_vars), colnames(colData(sce)))
if (length(missing_cols) > 0) {
  stop("Missing colData columns in SCE: ", paste(missing_cols, collapse = ", "))
}

# Stratification config
stratify <- isTRUE(args$stratification_var != "NULL") &&
  args$stratification_var %in% colnames(colData(sce))

if (stratify) {
  outdir <- file.path(
    args$output_dir,
    sprintf("de_%s_%s_%s", args$stratification_var, args$dependent_var, mod_name),
    celltype
  )
} else {
  outdir <- file.path(
    args$output_dir,
    sprintf("de_%s_%s", args$dependent_var, mod_name),
    celltype
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$output_dir, "dge_metacell_sce"), recursive = TRUE, showWarnings = FALSE)

cli::cli_text(c("Creating ", cli::col_green("{outdir}")))

# ------------------------------------------------------------------------------
# Pseudobulk generation
# ------------------------------------------------------------------------------
cli::cli_h2("Generating pseudobulk")

mean_pseudocell_size <- 30
min_pseudocell_size  <- 10

pseudocell_data_all <- .sconline.PseudobulkGeneration(
  argList = NULL,
  n_clusters = NULL,
  parsing.col.names = c("manifest"),
  use.sconline.cluster4parsing = FALSE,
  cluster_obj = NULL,
  pseudocell.size = mean_pseudocell_size,
  inputExpData = sce,
  min_size_limit = min_pseudocell_size,
  inputPhenoData = as.data.frame(colData(sce)),
  inputEmbedding = NULL,
  tol_level = 0.9,
  use.sconline.embeddings = FALSE,
  nPCs = 30,
  ncores = 1,
  rand_pseudobulk_mod = TRUE,
  organism = "Human"
)

cli::cli_text(c(
  "Pseudocells: {.val {ncol(pseudocell_data_all)}} | Individuals: {.val {length(unique(pseudocell_data_all$manifest))}}"
))
cli::cli_text("Pseudocells per individual:")
print(table(pseudocell_data_all$manifest))

# Optional releveling if those columns exist
if ("NeuropathologicalDiagnosis" %in% colnames(colData(pseudocell_data_all))) {
  pseudocell_data_all$NeuropathologicalDiagnosis <- relevel(
    factor(pseudocell_data_all$NeuropathologicalDiagnosis),
    ref = "Control"
  )
}
if ("BraakGroup" %in% colnames(colData(pseudocell_data_all))) {
  pseudocell_data_all$BraakGroup <- relevel(
    factor(pseudocell_data_all$BraakGroup),
    ref = "Control"
  )
}

# Background genes (filter)
tmp_nonzero <- apply(counts(pseudocell_data_all), 1, function(x) sum(x > 0))
tmp_cpm_sum <- rowSums(edgeR::cpm(as.matrix(counts(pseudocell_data_all))))

keep <- tmp_cpm_sum > max(0.1 * ncol(pseudocell_data_all), min(20, ncol(pseudocell_data_all) / 3))
keep <- keep & tmp_nonzero > max(0.05 * ncol(pseudocell_data_all), min(20, ncol(pseudocell_data_all) / 2))

bkg_genes <- rownames(pseudocell_data_all)[keep]

qs::qsave(
  pseudocell_data_all,
  file = file.path(args$output_dir, "dge_metacell_sce", sprintf("%s_metasce.qs", celltype))
)

# ------------------------------------------------------------------------------
# DE helper functions
# ------------------------------------------------------------------------------
default_vars <- c("pseudocell_size_scale", "QC_MT.pct", "nUMI_scaled")

scale_numeric <- function(x) as.numeric(scale(as.numeric(x)))

prepare_qc_covariates <- function(pb) {
  pb$pseudocell_size_scale <- scale_numeric(pb$pseudocell_size)
  
  pb$nUMI <- colSums(counts(pb))
  pb$nUMI_scaled <- log2(pb$nUMI)
  
  pb$nGene <- apply(counts(pb), 2, function(x) sum(x > 0))
  pb$nGene_scaled <- log2(pb$nGene)
  
  if ("mean_pc_mito" %in% colnames(colData(pb))) {
    pb$QC_MT.pct <- log2(pb$mean_pc_mito + 1)
  } else {
    pb$QC_MT.pct <- NA_real_
  }
  
  if (all(c("nGene", "nUMI") %in% colnames(colData(pb)))) {
    pc <- prcomp(as.matrix(as.data.frame(colData(pb))[, c("nGene", "nUMI")]), scale. = TRUE)
    pb$nGeneUMI <- pc$x[, 1]
  } else {
    pb$nGeneUMI <- NA_real_
  }
  
  pb
}

encode_dependent <- function(pb, dep_var, ref_class) {
  dep <- pb[[dep_var]]
  
  if (is.character(dep)) dep <- factor(dep)
  if (is.factor(dep)) {
    dep <- droplevels(dep)
    
    if (!is.null(ref_class) && !is.na(ref_class) && ref_class != "NULL" && ref_class %in% levels(dep)) {
      dep <- relevel(dep, ref = ref_class)
    }
    
    if (nlevels(dep) != 2) {
      stop("dependent_var as factor must have exactly 2 levels (or provide a numeric dependent variable).")
    }
    
    dep_num <- as.numeric(dep == levels(dep)[2])
    pb[[dep_var]] <- scale_numeric(dep_num)
    return(pb)
  }
  
  pb[[dep_var]] <- scale_numeric(dep)
  pb
}

drop_bad_covariates <- function(pb, covars) {
  exclude <- character(0)
  
  for (v in covars) {
    pb <- pb[, !is.na(pb[[v]])]
    
    if (is.character(pb[[v]]) || is.factor(pb[[v]])) {
      pb[[v]] <- droplevels(factor(pb[[v]]))
      if (nlevels(pb[[v]]) < 2) {
        exclude <- c(exclude, v)
      }
    }
  }
  
  list(pb = pb, exclude = unique(exclude))
}

run_limma_trend <- function(pb, dep_var, covariates, bkg_genes, random_effect = "manifest") {
  fit <- .sconline.fitLimmaFn(
    inputExpData = pb,
    covariates = covariates,
    randomEffect = random_effect,
    DEmethod = "Trend",
    prior.count = 1,
    bkg_genes = bkg_genes
  )
  
  fit2 <- limma::eBayes(fit$fit, trend = TRUE, robust = TRUE)
  
  res <- limma::topTable(
    fit2,
    number = nrow(fit2),
    adjust.method = "BH",
    coef = dep_var,
    confint = TRUE
  )
  
  res
}

format_results <- function(res, dep_var, covariates) {
  model_str <- paste0("~ ", paste(covariates, collapse = " + "))
  
  res %>%
    mutate(
      gene = rownames(.),
      contrast = dep_var,
      model = model_str
    ) %>%
    rename(pval = P.Value, padj = adj.P.Val) %>%
    select(gene, logFC, pval, padj, AveExpr, CI.L, CI.R, t, B, contrast, model)
}

run_de_block <- function(pseudocell_data, subset_label) {
  
  pb <- prepare_qc_covariates(pseudocell_data)
  
  # encode dependent variable (scaled numeric or 0/1)
  pb <- encode_dependent(pb, args$dependent_var, args$ref_class)
  
  # drop samples with NA in dependent/confounders, and drop confounders with <2 levels
  covars_to_check <- c(args$dependent_var, args$confounding_vars)
  tmp <- drop_bad_covariates(pb, covars_to_check)
  pb <- tmp$pb
  exclude_vars <- tmp$exclude
  
  n_samples <- length(unique(pb$manifest))
  if (n_samples < 5) {
    cli::cli_alert_warning("Not enough samples for DE ({n_samples} unique manifests). Skipping.")
    note_file <- if (stratify) {
      file.path(outdir, sprintf("%s_%s_in_%s.txt", celltype, args$dependent_var, subset_label))
    } else {
      file.path(outdir, sprintf("%s_%s.txt", celltype, args$dependent_var))
    }
    writeLines(paste0("Only ", n_samples, " unique samples"), con = note_file, useBytes = FALSE)
    return(invisible(NULL))
  }
  
  confounders <- setdiff(args$confounding_vars, exclude_vars)
  if (length(exclude_vars) > 0) {
    cli::cli_alert_info("Dropping confounders with <2 levels: {.emph {exclude_vars}}")
  }
  
  covariates <- c(args$dependent_var, default_vars, confounders)
  
  cli::cli_text(c(
    "Covariates: ",
    cli::col_blue("{paste(covariates, collapse = ', ')}")
  ))
  
  res <- run_limma_trend(
    pb = pb,
    dep_var = args$dependent_var,
    covariates = covariates,
    bkg_genes = bkg_genes,
    random_effect = "manifest"
  )
  
  res_fmt <- format_results(res, dep_var = args$dependent_var, covariates = covariates)
  
  out_file <- if (stratify) {
    file.path(outdir, sprintf("%s_%s_in_%s.tsv", celltype, args$dependent_var, subset_label))
  } else {
    file.path(outdir, sprintf("%s_%s.tsv", celltype, args$dependent_var))
  }
  
  write.table(res_fmt, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cli::cli_alert_success("Wrote: {out_file}")
  
  invisible(res_fmt)
}

# ------------------------------------------------------------------------------
# Run: stratified or not
# ------------------------------------------------------------------------------
if (stratify) {
  lvls <- unique(sce[[args$stratification_var]])
  for (lvl in lvls) {
    cli::cli_h2("Subsetting {.strong {lvl}} from {.emph {args$stratification_var}}")
    pb_sub <- pseudocell_data_all[, pseudocell_data_all[[args$stratification_var]] %in% lvl]
    run_de_block(pb_sub, subset_label = as.character(lvl))
  }
} else {
  cli::cli_h2("Running without stratification")
  run_de_block(pseudocell_data_all, subset_label = "all")
}
