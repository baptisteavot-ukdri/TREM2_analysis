suppressPackageStartupMessages({
  library(qs)
  library(vroom)
  library(tidyverse)
  library(Seurat)
  library(SingleCellExperiment)
  library(stringr)
  library(future)
  library(BiocParallel)
  library(variancePartition)
  library(scCustomize)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "config.R"))

paths <- list(
  micro_sce     = cfg$sce_purified,
  micro_seu     = cfg$seu_round2,
  auc_csv       = cfg$auc_csv,
  out_base      = cfg$pyscenic_out_dir,
  out_markers   = file.path(cfg$pyscenic_out_dir, "out/wilcox_subclusters_Micro.tsv"),
  out_seu_auc   = file.path(cfg$pyscenic_out_dir, "diffmod_subcelltypes/seu_with_auc_Total4G8Density.qs"),
  out_pb_R47H   = file.path(cfg$pyscenic_out_dir, "TREM2/out/wilcox_subclusters_Micro_pb_R47H.tsv"),
  out_dream_dir = file.path(cfg$pyscenic_out_dir, "diffmod_subcelltypes/Total4G8Density"),
  out_umap_pdf  = file.path(cfg$pyscenic_out_dir, "diffmod_subcelltypes/UMAP/modulescore_Micro_UMAP.pdf")
)

dir.create(paths$out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$out_markers), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$out_seu_auc), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$out_pb_R47H), recursive = TRUE, showWarnings = FALSE)
dir.create(paths$out_dream_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$out_umap_pdf), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
stop_if_missing <- function(x, cols, what = "object") {
  missing <- setdiff(cols, colnames(x))
  if (length(missing) > 0) stop(what, " is missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

read_auc_matrix <- function(auc_csv) {
  auc <- vroom::vroom(auc_csv, show_col_types = FALSE) %>%
    tibble::column_to_rownames("Regulon") %>%
    as.data.frame(check.names = FALSE)
  
  as.matrix(auc)
}

add_auc_assay <- function(seu, auc_mat, assay_name = "regulons_auc") {
  overlap <- intersect(colnames(seu), colnames(auc_mat))
  if (length(overlap) < 50) {
    stop("Very small overlap between Seurat cells and AUC columns (n=", length(overlap), "). Check inputs.")
  }
  
  auc_mat <- auc_mat[, overlap, drop = FALSE]
  seu <- seu[, overlap]
  
  assay_obj <- CreateAssay5Object(data = auc_mat)
  seu[[assay_name]] <- assay_obj
  seu
}

make_manifest_subcluster <- function(seu, subcluster_col = "subclusters_label", manifest_col = "manifest") {
  stop_if_missing(seu@meta.data, c(subcluster_col, manifest_col), "Seurat meta.data")
  seu$manifest_subcluster <- paste0(seu[[manifest_col]][, 1], "_", seu[[subcluster_col]][, 1])
  seu
}

split_manifest_subcluster <- function(pb, col = "orig.ident") {
  x <- str_split_fixed(pb@meta.data[[col]], "_", 2)
  pb@meta.data$manifest <- x[, 1]
  pb@meta.data$subcluster_label <- x[, 2]
  pb
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# -------------------------------------------------------------------------
# Load objects
# -------------------------------------------------------------------------
micro_sce <- qs::qread(paths$micro_sce)  # used for colData / covariates
micro_seu <- qs::qread(paths$micro_seu)

# Ensure we have the same cells between SCE and Seurat (for metadata transfer)
overlap_cells <- intersect(colnames(micro_seu), colnames(micro_sce))
if (length(overlap_cells) == 0) stop("No overlap between Seurat and Micro_sce cells.")
micro_seu <- micro_seu[, overlap_cells]
micro_sce <- micro_sce[, overlap_cells]

# -------------------------------------------------------------------------
# Read AUC matrix and add to Seurat as an assay
# -------------------------------------------------------------------------
auc_mat <- read_auc_matrix(paths$auc_csv)
micro_seu <- add_auc_assay(micro_seu, auc_mat, assay_name = "regulons_auc")

# -------------------------------------------------------------------------
# Marker detection on regulon AUC (per subcluster)
# -------------------------------------------------------------------------
stop_if_missing(micro_seu@meta.data, c("subclusters_label"), "Seurat meta.data")

markers <- FindAllMarkers(
  object   = micro_seu,
  assay    = "regulons_auc",
  group.by = "subclusters_label"
)

write_tsv(markers, paths$out_markers)

# Save Seurat with AUC
qs::qsave(micro_seu, paths$out_seu_auc)

# -------------------------------------------------------------------------
# Pseudobulk (R47H example) + markers on pseudobulked AUC
# -------------------------------------------------------------------------
micro_seu <- qs::qread(paths$out_seu_auc)

seu_R47H <- subset(micro_seu, subset = TREM2Variant == "R47H")

seu_R47H <- make_manifest_subcluster(seu_R47H, subcluster_col = "subclusters_label", manifest_col = "manifest")

pb_auc <- PseudobulkExpression(
  object = seu_R47H,
  assay = "regulons_auc",
  method = "average",
  group.by = "manifest_subcluster",
  return.seurat = TRUE
)

pb_auc <- split_manifest_subcluster(pb_auc, col = "orig.ident")
Idents(pb_auc) <- "subcluster_label"

markers_pb <- FindAllMarkers(pb_auc, group.by = "subcluster_label")
write_tsv(markers_pb, paths$out_pb_R47H)

# -------------------------------------------------------------------------
# Add Total4G8Density metadata from Micro_sce to Seurat (if needed)
# -------------------------------------------------------------------------
if ("Total4G8Density" %in% colnames(colData(micro_sce))) {
  micro_seu <- AddMetaData(micro_seu, metadata = colData(micro_sce)$Total4G8Density, col.name = "Total4G8Density")
}

# -------------------------------------------------------------------------
# dream: association between Total4G8Density and regulon AUC
#   - within each subcluster and stratification layer
# -------------------------------------------------------------------------
strats <- c("TREM2Variant", "APOEgroup", "CD33Group")
celltypes <- c("ProliferativeMicro", "DAM", "DAM_HLA", "CRM", "PVM")

# Parallel backend for dream
param <- SnowParam(workers = max(1, future::availableCores() - 1), type = "SOCK", progressbar = TRUE)
register(param)

make_formula <- function(strat) {
  if (strat == "APOEgroup") {
    ~ Total4G8Density + (1 | manifest) + total_features_by_counts + Sex + Age + PostMortemInterval + BrainRegion + TREM2Variant + CD33Group
  } else if (strat == "TREM2Variant") {
    ~ Total4G8Density + (1 | manifest) + total_features_by_counts + Sex + Age + PostMortemInterval + BrainRegion + APOEgroup + CD33Group
  } else if (strat == "CD33Group") {
    ~ Total4G8Density + (1 | manifest) + total_features_by_counts + Sex + Age + PostMortemInterval + BrainRegion + TREM2Variant + APOEgroup
  } else {
    stop("Unknown stratification var: ", strat)
  }
}

layers_for_strat <- function(strat) {
  if (strat == "APOEgroup") return(c("APOE4-neg", "APOE4-pos"))
  if (strat == "CD33Group") return(c("CD33var", "CV"))
  if (strat == "TREM2Variant") return(c("R62H", "CV", "R47H"))
  stop("Unknown strat: ", strat)
}

for (strat in strats) {
  form <- make_formula(strat)
  layers <- layers_for_strat(strat)
  
  out_dir_strat <- file.path(paths$out_dream_dir, strat)
  dir.create(out_dir_strat, recursive = TRUE, showWarnings = FALSE)
  
  for (layer in layers) {
    for (subcelltype in celltypes) {
      
      # Subset
      seu_subset <- subset(
        micro_seu,
        subset = subclusters_label == subcelltype & .data[[strat]] == layer
      )
      
      # Skip empty / tiny subsets
      if (ncol(seu_subset) < 50) {
        cli::cli_alert_warning("Skipping {subcelltype} / {strat}={layer} (n={ncol(seu_subset)})")
        next
      }
      
      meta <- seu_subset@meta.data %>%
        mutate(
          Age = as.vector(scale(Age)),
          total_features_by_counts = as.vector(scale(total_features_by_counts)),
          PostMortemInterval = as.vector(scale(PostMortemInterval)),
          NeuropathologicalDiagnosis = factor(NeuropathologicalDiagnosis, levels = c("Control", "AD"))
        )
      
      if (!("Total4G8Density" %in% colnames(meta))) {
        cli::cli_alert_warning("Missing Total4G8Density; skipping {subcelltype} / {strat}={layer}")
        next
      }
      
      # Expression matrix for dream: genesets x samples
      auc_expr <- as.matrix(seu_subset@assays$regulons_auc$data)
      
      fitmm <- try(dream(auc_expr, form, meta), silent = TRUE)
      if (inherits(fitmm, "try-error")) {
        cli::cli_alert_warning("dream failed for {subcelltype} / {strat}={layer}")
        next
      }
      
      res <- variancePartition::topTable(fitmm, coef = 2, number = Inf) %>%
        mutate(
          geneset = rownames(.),
          form = paste0(as.character(form), collapse = ""),
          subcluster_label = subcelltype,
          stratification_var = strat,
          stratification_level = layer
        )
      
      out_file <- file.path(out_dir_strat, sprintf("res_dream_%s_%s.tsv", subcelltype, layer))
      write.table(res, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      cli::cli_alert_success("Wrote {out_file}")
    }
  }
}

# -------------------------------------------------------------------------
# Module scores + UMAP feature plots
# -------------------------------------------------------------------------
micro_seu <- qs::qread(paths$out_seu_auc)
DefaultAssay(micro_seu) <- "regulons_auc"

regulons_ms <- paste0(rownames(micro_seu), "_ModuleScore")

micro_seu <- AddModuleScore(
  object = micro_seu,
  features = as.list(rownames(micro_seu)),   # one feature set per regulon
  layer = "data",
  assay = "regulons_auc",
  nbin = 10,
  ctrl = 4,
  name = "Regulon"
)

# Normalise the column names created by AddModuleScore
colnames(micro_seu@meta.data) <- gsub("_ModuleScore[0-9]+", "_ModuleScore", colnames(micro_seu@meta.data))

# Example plot
scCustomize::FeaturePlot_scCustom(micro_seu, features = "CEBPZ(+)_ModuleScore", na_cutoff = 0.1)

pdf(paths$out_umap_pdf, width = 7, height = 7)
for (feat in regulons_ms) {
  p <- scCustomize::FeaturePlot_scCustom(micro_seu, features = feat, na_cutoff = 0.05)
  ggsave(filename = file.path(dirname(paths$out_umap_pdf), paste0(feat, ".png")),
         plot = p, height = 10, width = 12, dpi = 300)
}
dev.off()
