# =============================================================================
# TREM2_analysis — Central Configuration
# =============================================================================
#
# Before running any script, set the repository root:
#
#   export TREM2_ANALYSIS_ROOT=/path/to/TREM2_analysis   # bash / PBS / SLURM
#
# Optionally override individual base paths:
#
#   export TREM2_PROJECT_RDS=/rds/general/project/ukdrmultiomicsproject/live
#   export TREM2_USER_RDS=/rds/general/user/bavot/home
#   export TREM2_MNT_DATA=/mnt/data/bavot/TREM2
#
# Then source this file from any script with:
#   source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "config.R"))
# =============================================================================

cfg <- list()

# ── Helper: prefer env-var, fall back to default ──────────────────────────────
.env_or <- function(var, default) {
  v <- Sys.getenv(var, unset = "")
  if (nchar(v) > 0) v else default
}

# ── HPC base paths ─────────────────────────────────────────────────────────────
cfg$project_rds <- .env_or(
  "TREM2_PROJECT_RDS",
  "/rds/general/project/ukdrmultiomicsproject/live"
)
cfg$user_rds <- .env_or(
  "TREM2_USER_RDS",
  "/rds/general/user/bavot/home"
)
cfg$mnt_data <- .env_or(
  "TREM2_MNT_DATA",
  "/mnt/data/bavot/TREM2"
)

# ── Input data ────────────────────────────────────────────────────────────────
cfg$sce_original <- file.path(
  cfg$project_rds,
  "MAP_analysis/TREM2_enriched_scflow/dge/split_sce/celltype_sce/celltype_sce/original/Micro_sce.qs"
)
cfg$sce_purified <- file.path(
  cfg$project_rds,
  "MAP_analysis/TREM2_enriched_scflow/dge/split_sce/celltype_sce/celltype_sce/purified_with_Total4G8Density/Micro_sce.qs"
)
cfg$seu_micro  <- file.path(cfg$mnt_data, "SEU/Micro_seu.qs")
cfg$seu_round2 <- file.path(
  cfg$project_rds,
  "MAP_analysis/TREM2_enriched_scflow/subclustering/Micro/subclusters_round2/v2/seu.qs"
)

# ── pySCENIC ──────────────────────────────────────────────────────────────────
cfg$pyscenic_resources <- file.path(
  cfg$project_rds,
  "MAP_pipelines/snRNAseq/pySCENIC/resources"
)
cfg$auc_csv      <- file.path(cfg$user_rds, "pySCENIC/out/Micro_no_hvg.auc_mtx.csv")
cfg$regulons_csv <- file.path(
  cfg$project_rds,
  "MAP_analysis/TREM2_enriched_scflow/pySCENIC/out/wholeMicro/Micro.regulons.filtered.csv"
)

# ── MAP pipeline helper scripts ───────────────────────────────────────────────
cfg$pseudobulk_script   <- file.path(
  cfg$project_rds,
  "MAP_pipelines/snRNAseq/pseudobulk_dge/generate_pseudobulk_limma_trend_deg_gazestani.r"
)
cfg$heatmap_script      <- file.path(
  cfg$project_rds,
  "MAP_pipelines/snRNAseq/additional_scripts/plot_marker_heatmap.r"
)
cfg$celltype_map_script <- file.path(
  cfg$project_rds,
  "MAP_pipelines/snRNAseq/additional_scripts/map_celltype_seu.r"
)
cfg$ctd_dir <- file.path(
  cfg$project_rds,
  "MAP_pipelines/snRNAseq/assets/ctd"
)

# ── Output directories ────────────────────────────────────────────────────────
cfg$subclustering_dir  <- file.path(cfg$user_rds, "scflow/subclustering")
cfg$pyscenic_out_dir   <- file.path(cfg$user_rds, "pySCENIC")
cfg$regulons_out_dir   <- file.path(cfg$user_rds, "TREM2/Micro/regulons")
cfg$dirichlet_work_dir <- file.path(cfg$user_rds, "TREM2/Micro/dirichlet")

# ── Regulon visualisation inputs ──────────────────────────────────────────────
cfg$astro_markers_dir <- file.path(
  cfg$mnt_data, "regulon_analysis/Astro_subcelltype_specific_markers"
)
cfg$micro_markers_dir <- file.path(
  cfg$mnt_data, "regulon_analysis/Micro_subcelltype_specific_markers"
)

# ── Analysis parameters ───────────────────────────────────────────────────────
cfg$params <- list(
  # Subclustering
  nfeatures           = 2000L,
  pca_dims_harmony    = 1:30,
  umap_dims           = 1:30,
  neighbors_dims      = 1:20,
  umap_epochs         = 200L,
  cluster_resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5),

  # Pseudobulk / DEG
  mean_pseudocell_size = 30L,
  min_pseudocell_size  = 10L,

  # pySCENIC preprocessing
  min_pct_cells_expressed = 0.10,

  # Differential modularity
  min_cells_per_subset = 50L,

  # Dirichlet
  trem2_variants = c("CV", "R47H", "R62H"),
  facet_vars     = c("TREM2Variant", "BrainRegion", "APOEgroup", "CD33Group"),

  # Visualisation thresholds
  p_adj_thresh = 0.05,
  lfc_thresh   = 0.5
)

invisible(cfg)
