suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(tidyverse)
  library(harmony)
  library(ggplot2)
  library(ComplexHeatmap)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
work_dir <- "/rds/general/user/bavot/home/scflow/subclustering/"
setwd(work_dir)

celltype <- "Micro"
clustering_round <- "round1"

sce_path <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_enriched_scflow/dge/split_sce/celltype_sce/celltype_sce/original/Micro_sce.qs"

outdir <- file.path(work_dir, paste0("subclustering_", celltype, "_full_cohort"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sub_dir <- file.path(outdir, sprintf("subclusters_%s", clustering_round))
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

# Harmony / clustering settings
nfeatures <- 2000
pca_dims_harmony <- 1:30
umap_dims <- 1:30
neighbors_dims <- 1:20
umap_epochs <- 200
cluster_resolution <- 0.5
cluster_col <- paste0("RNA_snn_res.", cluster_resolution)

# ------------------------------------------------------------------------------
# Load SCE and build Seurat object
# ------------------------------------------------------------------------------
sce <- qs::qread(sce_path)
rownames(sce) <- rowData(sce)$gene
sce$manifest <- droplevels(sce$manifest)

seu <- CreateSeuratObject(
  counts = SingleCellExperiment::counts(sce),
  meta.data = as.data.frame(SingleCellExperiment::colData(sce))
)

# Add gene-level metadata to the RNA assay
seu[["RNA"]] <- AddMetaData(
  object = seu[["RNA"]],
  metadata = as.data.frame(SummarizedExperiment::rowData(sce))
)

# Remove any pre-existing reduction that could cause confusion
seu@reductions$umap_by_individual <- NULL

# ------------------------------------------------------------------------------
# Preprocessing: normalisation, HVGs, scaling, PCA
# ------------------------------------------------------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu))

pdf(file.path(outdir, "elbow_plot.pdf"))
print(ElbowPlot(seu))
dev.off()

# ------------------------------------------------------------------------------
# Integration + clustering
# ------------------------------------------------------------------------------
seu <- RunHarmony(seu, group.by.vars = "manifest")
seu <- RunUMAP(seu, reduction = "harmony", dims = umap_dims, n.epochs = umap_epochs)
seu <- FindNeighbors(seu, reduction = "harmony", dims = neighbors_dims)
seu <- FindClusters(seu, resolution = cluster_resolution)

Idents(seu) <- cluster_col

# ------------------------------------------------------------------------------
# QC / overview plots
# ------------------------------------------------------------------------------
pdf(file.path(outdir, "UMAP.pdf"))
print(
  DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.01) +
    ggtitle("Clustering round 1")
)

for (feat in c("pct4G8PositiveArea", "pc_mito", "total_features_by_counts")) {
  if (feat %in% colnames(seu@meta.data)) {
    print(FeaturePlot(seu, reduction = "umap", features = feat))
  }
}

for (grp in c("NeuropathologicalDiagnosis", "TREM2Variant")) {
  if (grp %in% colnames(seu@meta.data)) {
    print(DimPlot(seu, reduction = "umap", group.by = grp, pt.size = 0.01))
  }
}
dev.off()

# Save intermediate object
qs::qsave(seu, file = file.path(sub_dir, "seu_round1.qs"))

# ------------------------------------------------------------------------------
# Differential markers (per cluster)
# ------------------------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
Idents(seu) <- cluster_col

markers <- FindAllMarkers(
  object = seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
) %>%
  select(gene, everything())

markers_path <- file.path(sub_dir, sprintf("markers_%s.tsv", cluster_col))
write.table(markers, file = markers_path, sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# Marker heatmap
# ------------------------------------------------------------------------------
heatmap_fun <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/plot_marker_heatmap.r"
source(heatmap_fun)

# pick up to 10 significant markers per cluster, then de-duplicate genes
markers_for_ht <- read.delim(markers_path) %>%
  group_by(cluster) %>%
  filter(p_val_adj <= 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = TRUE) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)

ht <- plot_marker_heatmap(
  seu = seu,
  marker = markers_for_ht,
  cluster_col = cluster_col,
  show_row_names = TRUE
)

png(file.path(sub_dir, sprintf("ht_%s.png", cluster_col)),
    height = 20, width = 10, units = "in", res = 300)
draw(ht)
dev.off()

# UMAP labelled by clusters
p_cluster <- DimPlot(seu, reduction = "umap", group.by = cluster_col, label = TRUE, pt.size = 0.01)
ggsave(filename = file.path(sub_dir, sprintf("dimplot_%s.png", cluster_col)),
       plot = p_cluster, height = 7, width = 7, units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# EWCE-based cell type annotation
# ------------------------------------------------------------------------------
annot_fun <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/map_celltype_seu.r"
source(annot_fun)

seu <- map_celltypes_seu(
  seu,
  ctd_dir = "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/assets/ctd",
  clusters_colname = cluster_col
)

Idents(seu) <- seu$subcluster_celltype

p_ewce <- DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.01) +
  ggtitle("UMAP subcluster cell type (EWCE)")

ggsave(filename = file.path(sub_dir, "umap_subcluster_celltype_ewce.png"),
       plot = p_ewce, height = 5, width = 7, units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# Manual relabel (optional, keeps your original mapping)
# ------------------------------------------------------------------------------
if ("seurat_clusters" %in% colnames(seu@meta.data)) {
  seu@meta.data <- seu@meta.data %>%
    mutate(
      subcluster_celltype = case_when(
        seurat_clusters %in% c(13) ~ "Peri",
        seurat_clusters %in% c(3, 4, 7, 10) ~ "Micro",
        seurat_clusters %in% c(11, 12) ~ "Endo",
        TRUE ~ as.character(subcluster_celltype)
      )
    )
}

p_relab <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "subcluster_celltype",
  label = TRUE,
  pt.size = 0.01
)

ggsave(filename = file.path(sub_dir, "umap_subcluster_celltype.png"),
       plot = p_relab, height = 5, width = 7, units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# Save final object
# ------------------------------------------------------------------------------
qs::qsave(seu, file = file.path(sub_dir, "seu_relabelled.qs"))
