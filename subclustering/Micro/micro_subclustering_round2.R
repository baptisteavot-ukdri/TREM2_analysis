suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(tidyverse)
  library(harmony)
  library(ggplot2)
  library(ComplexHeatmap)
  library(scCustomize)
  library(paletteer)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
work_dir <- "/rds/general/user/bavot/home/scflow/subclustering/"
setwd(work_dir)

celltype <- "Micro"
outdir <- file.path(work_dir, paste0("subclustering_", celltype, "_full_cohort"))

clustering_round <- "round2"
sub_dir <- file.path(outdir, sprintf("subclusters_%s", clustering_round))

dirs <- list(
  sub_dir = sub_dir,
  dimplot = file.path(sub_dir, "dimplot"),
  featureplot = file.path(sub_dir, "featureplot"),
  marker_genes = file.path(sub_dir, "marker_genes"),
  heatmap = file.path(sub_dir, "heatmap"),
  marker_enrichment = file.path(sub_dir, "marker_enrichment")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Inputs
seu_round1_path <- file.path(outdir, "subclusters_round1", "seu_relabelled.qs")

geneset_paths <- list(
  Gazestani_2024 = "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/manuscript_geneset/Gazestani_2024.qs",
  Mancuso_2024   = "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/manuscript_geneset/Mancuso_2024.qs"
)

# External helpers
plot_marker_heatmap_path <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/plot_marker_heatmap.r"
filter_genes_path        <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/filter_genes_expressed_seu.r"
enrichment_path          <- "/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/enrichment_custom.r"

source(plot_marker_heatmap_path)
source(filter_genes_path)
source(enrichment_path)

# Analysis settings
umap_dims <- 1:20
umap_epochs <- 200
res_grid <- seq(0.1, 0.5, by = 0.1)
min_marker_fdr <- 0.05
min_marker_logfc <- 0.25
markers_per_cluster <- 10

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) read.delim(path, check.names = FALSE)

select_heatmap_markers <- function(markers, n = 10, fdr = 0.05) {
  markers %>%
    group_by(cluster) %>%
    filter(p_val_adj <= fdr) %>%
    slice_max(order_by = avg_log2FC, n = n, with_ties = TRUE) %>%
    ungroup() %>%
    distinct(gene, .keep_all = TRUE)
}

plot_and_save <- function(p, path, h = 5, w = 7, dpi = 300) {
  ggsave(filename = path, plot = p, height = h, width = w, units = "in", dpi = dpi)
}

# ------------------------------------------------------------------------------
# Load round 1 object and prepare for round 2 clustering
# ------------------------------------------------------------------------------
seu <- qs::qread(seu_round1_path)

# Drop a specific cluster (kept from your original script)
if ("seurat_clusters" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = seurat_clusters != 16)
}

p_init <- DimPlot(seu, group.by = "subcluster_celltype", reduction = "umap", label = TRUE)
plot_and_save(p_init, file.path(dirs$dimplot, "initial_dimplot.png"), h = 5, w = 7)

# Keep Micro only
seu <- subset(seu, subset = subcluster_celltype %in% celltype)
seu$manifest <- droplevels(seu$manifest)

p_micro <- DimPlot(seu, group.by = "subcluster_celltype", reduction = "umap", label = TRUE)
plot_and_save(p_micro, file.path(dirs$dimplot, "micro_only.png"), h = 5, w = 7)

# ------------------------------------------------------------------------------
# Filter low-nuclei samples by manifest
# ------------------------------------------------------------------------------
sample_pc <- table(seu$manifest)
sample_threshold <- ceiling(quantile(sample_pc, probs = 0.05))  # bottom 5%
keep_manifests <- names(sample_pc)[sample_pc > sample_threshold]

seu <- subset(seu, subset = manifest %in% keep_manifests)
seu$manifest <- droplevels(seu$manifest)

# Save distribution plot
png(file.path(dirs$dimplot, "nuclei_per_manifest_hist.png"), width = 1600, height = 1200, res = 200)
hist(sample_pc, breaks = 100, main = "Nuclei per manifest (pre-filter)", xlab = "Nuclei")
dev.off()

# ------------------------------------------------------------------------------
# Recompute UMAP from Harmony reduction (already present from round 1)
# ------------------------------------------------------------------------------
seu <- RunUMAP(seu, reduction = "harmony", dims = umap_dims, n.epochs = umap_epochs)

p_harmony <- DimPlot(seu, group.by = "subcluster_celltype", reduction = "umap", label = TRUE)
plot_and_save(p_harmony, file.path(dirs$dimplot, "harmony_reduction.png"), h = 5, w = 7)

# Example marker feature
p_p2ry12 <- FeaturePlot(seu, features = "P2RY12", reduction = "umap")
plot_and_save(p_p2ry12, file.path(dirs$featureplot, "P2RY12.png"), h = 5, w = 7)

# ------------------------------------------------------------------------------
# Neighbours + clustering grid search
# ------------------------------------------------------------------------------
seu <- FindNeighbors(seu, reduction = "harmony", dims = umap_dims)

nclust_dt <- map_dfr(res_grid, function(r) {
  seu_tmp <- FindClusters(seu, resolution = r)
  cl_col <- paste0("RNA_snn_res.", r)
  tibble(resolution = r, number_cluster = n_distinct(seu_tmp[[cl_col]][, 1]))
})

write_tsv(nclust_dt, file.path(sub_dir, "number_of_cluster_by_resolution.tsv"))

p_nclust <- ggplot(nclust_dt, aes(x = resolution, y = number_cluster)) +
  geom_point() +
  theme_bw()

plot_and_save(p_nclust, file.path(sub_dir, "number_of_cluster_by_resolution.png"), h = 5, w = 7)

# Compute clusters at all resolutions on the same object (Seurat adds columns)
for (r in res_grid) {
  seu <- FindClusters(seu, resolution = r)
}

qs::qsave(seu, file.path(sub_dir, "seu.qs"))

# ------------------------------------------------------------------------------
# Marker genes and heatmaps per resolution
# ------------------------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
clustering_cols <- paste0("RNA_snn_res.", res_grid)

for (cl_col in clustering_cols) {
  Idents(seu) <- cl_col
  
  markers <- FindAllMarkers(
    object = seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  ) %>%
    select(gene, everything())
  
  markers_path <- file.path(dirs$marker_genes, sprintf("markers_%s.tsv", cl_col))
  write_tsv(markers, markers_path)
  
  # Dimplot
  p <- DimPlot(seu, reduction = "umap", group.by = cl_col, label = TRUE, pt.size = 0.1)
  plot_and_save(p, file.path(dirs$dimplot, sprintf("%s.png", cl_col)), h = 5, w = 7)
  
  # Heatmap from top markers per cluster
  markers_ht <- select_heatmap_markers(read_tsv(markers_path), n = markers_per_cluster, fdr = min_marker_fdr)
  
  ht <- plot_marker_heatmap(
    seu = seu,
    marker = markers_ht,
    cluster_col = cl_col,
    show_row_names = TRUE,
    fontsize = 6
  )
  
  png(file.path(dirs$heatmap, sprintf("ht_%s.png", cl_col)),
      height = 20, width = 10, units = "in", res = 300)
  draw(ht)
  dev.off()
}

# Reload clean object for downstream steps
seu <- qs::qread(file.path(sub_dir, "seu.qs"))

# ------------------------------------------------------------------------------
# Gene set module scores (Gazestani + Mancuso)
# ------------------------------------------------------------------------------
geneset_l <- lapply(geneset_paths, qs::qread)

for (nm in names(geneset_l)) {
  seu <- AddModuleScore(seu, features = geneset_l[nm], name = nm)
}

# Example feature plot
p_dam1 <- FeaturePlot_scCustom(seu, features = "DAM1", reduction = "umap")
plot_and_save(p_dam1, file.path(dirs$featureplot, "DAM1.png"), h = 5, w = 7)

# Plot each geneset module score (AddModuleScore appends "1")
for (score_col in paste0(names(geneset_l), "1")) {
  file_stub <- str_replace(score_col, "1$", "")
  p <- FeaturePlot_scCustom(seu, features = score_col, reduction = "umap") + ggtitle(file_stub)
  plot_and_save(p, file.path(dirs$featureplot, sprintf("%s.png", file_stub)), h = 5, w = 7)
}

# ------------------------------------------------------------------------------
# Marker enrichment against expressed genes (per resolution)
# ------------------------------------------------------------------------------
seu_f <- filter_seu_genes(seu, min_counts = 1, min_cells_pc = 0.05)
expressed_genes <- rownames(seu_f) %>% as.character()
rm(seu_f)

for (cl_col in clustering_cols) {
  
  markers_path <- file.path(dirs$marker_genes, sprintf("markers_%s.tsv", cl_col))
  markers <- read_tsv(markers_path) %>%
    filter(p_val_adj <= min_marker_fdr, avg_log2FC >= min_marker_logfc) %>%
    arrange(desc(avg_log2FC))
  
  markers_l <- split(markers$gene, markers$cluster)
  
  enr <- map_dfr(names(markers_l), function(k) {
    res <- enrichment_custom(
      genes = markers_l[[k]],
      reference = expressed_genes,
      genesets = geneset_l
    )
    res$cluster <- k
    res
  })
  
  out_tsv <- file.path(dirs$marker_enrichment, sprintf("marker_enrichment_%s.tsv", cl_col))
  write_tsv(enr, out_tsv)
  
  # Plot
  plot_df <- enr %>%
    filter(padj <= 0.05, enrichment_ratio >= 1.5) %>%
    mutate(
      geneset_name = case_when(
        TermID %in% names(geneset_l$Gazestani_2024) ~ "Gazestani et al. 2024",
        TermID %in% names(geneset_l$Mancuso_2024)   ~ "Mancuso et al. 2023",
        TRUE ~ "Other"
      )
    )
  
  p <- ggplot(plot_df, aes(x = cluster, y = TermID)) +
    geom_point(aes(size = enrichment_ratio, fill = padj), shape = 21, color = "black") +
    labs(x = "Cluster (this study)", y = NULL, title = cl_col) +
    facet_wrap(~ geneset_name, scales = "free_y") +
    scale_fill_gradient(
      low = "gold", high = "navy", name = "FDR",
      guide = guide_colorbar(reverse = TRUE),
      limits = c(0, 0.05)
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", size = 14),
      axis.text  = element_text(colour = "black", size = 12),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(dirs$marker_enrichment, sprintf("marker_enrichment_%s.pdf", cl_col)),
         plot = p, width = 10, height = 4, units = "in", dpi = 300, device = "pdf")
}

# ------------------------------------------------------------------------------
# Split UMAP by TREM2Variant (optional panel)
# ------------------------------------------------------------------------------
if ("TREM2Variant" %in% colnames(seu@meta.data)) {
  p_split <- DimPlot(seu, reduction = "umap", split.by = "TREM2Variant")
  plot_and_save(p_split, file.path(dirs$dimplot, "post_enrichment.png"), h = 5, w = 7)
}

# ------------------------------------------------------------------------------
# Manual labelling (based on a chosen resolution)
# ------------------------------------------------------------------------------
label_resolution <- 0.5
label_col <- paste0("RNA_snn_res.", label_resolution)

stopifnot(label_col %in% colnames(seu@meta.data))

seu$subclusters_manual <- as.character(seu[[label_col]][, 1])
seu$subclusters_manual <- case_when(
  seu$subclusters_manual %in% c("0", "1", "3") ~ "Micro1",
  seu$subclusters_manual %in% c("4", "5")      ~ "Micro2",
  seu$subclusters_manual %in% c("6")           ~ "Micro3",
  seu$subclusters_manual %in% c("2")           ~ "Micro4",
  seu$subclusters_manual %in% c("7")           ~ "PVM",
  seu$subclusters_manual %in% c("8")           ~ "ProliferativeMicro",
  TRUE ~ seu$subclusters_manual
)

seu$subclusters_label <- as.character(seu[[label_col]][, 1])
seu$subclusters_label <- case_when(
  seu$subclusters_label %in% c("0", "1", "3") ~ "HM",
  seu$subclusters_label %in% c("4", "5")      ~ "DAM_HLA",
  seu$subclusters_label %in% c("6")           ~ "CRM",
  seu$subclusters_label %in% c("2")           ~ "DAM",
  seu$subclusters_label %in% c("7")           ~ "PVM",
  seu$subclusters_label %in% c("8")           ~ "ProliferativeMicro",
  TRUE ~ seu$subclusters_label
)

palette_choice <- paletteer::paletteer_d("ggsci::lanonc_lancet")

p_manual <- DimPlot(seu, reduction = "umap", group.by = "subclusters_manual", pt.size = 0.1) +
  scale_color_manual(values = palette_choice[c(3, 1, 2, 4, 5, 6)], aesthetics = c("colour", "fill"), name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    legend.text = element_text(color = "black", size = 12)
  )

plot_and_save(p_manual, file.path(sub_dir, "umap_subclusters_manual.png"), h = 5, w = 8)

p_label <- DimPlot(seu, reduction = "umap", group.by = "subclusters_label", label = TRUE, pt.size = 0.1) +
  scale_color_manual(values = palette_choice[c(2, 4, 1, 3, 5, 6)], aesthetics = c("colour", "fill"), name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 14),
    panel.grid = element_blank(),
    legend.text = element_text(color = "black", size = 12)
  )

plot_and_save(p_label, file.path(sub_dir, "umap_subclusters_label.png"), h = 5, w = 8)

# Save final object
qs::qsave(seu, file.path(sub_dir, "seu.qs"))

# ------------------------------------------------------------------------------
# Markers + heatmap for manual labels
# ------------------------------------------------------------------------------
Idents(seu) <- "subclusters_manual"
DefaultAssay(seu) <- "RNA"

markers_manual <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
) %>%
  select(gene, everything())

write_tsv(markers_manual, file.path(sub_dir, "markers_subclusters_manual.tsv"))

markers_ht <- markers_manual %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = TRUE) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE)

ht_manual <- plot_marker_heatmap(
  seu = seu,
  marker = markers_ht,
  cluster_col = "subclusters_manual",
  show_row_names = TRUE,
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 45,
  ht_width = 3
)

png(file.path(sub_dir, "ht_subclusters_manual.png"),
    height = 7, width = 4, units = "in", res = 300)
draw(ht_manual)
dev.off()

# ------------------------------------------------------------------------------
# Marker enrichment for manual labels
# ------------------------------------------------------------------------------
markers_l <- markers_manual %>%
  filter(p_val_adj <= min_marker_fdr, avg_log2FC >= min_marker_logfc) %>%
  arrange(desc(avg_log2FC)) %>%
  split(.$cluster) %>%
  map(~ .x$gene)

enr_manual <- map_dfr(names(markers_l), function(k) {
  res <- enrichment_custom(genes = markers_l[[k]], reference = expressed_genes, genesets = geneset_l)
  res$cluster <- k
  res
})

write_tsv(enr_manual, file.path(sub_dir, "marker_enrichment_subclusters_manual.tsv"))

plot_df <- enr_manual %>%
  filter(padj <= 0.05, enrichment_ratio >= 1.5) %>%
  mutate(
    geneset_name = case_when(
      TermID %in% names(geneset_l$Gazestani_2024) ~ "Gazestani et al. 2024",
      TermID %in% names(geneset_l$Mancuso_2024)   ~ "Mancuso et al. 2023",
      TRUE ~ "Other"
    )
  )

p <- ggplot(plot_df, aes(x = cluster, y = TermID)) +
  geom_point(aes(size = enrichment_ratio, fill = padj), shape = 21, color = "black") +
  labs(x = "Cluster (this study)", y = NULL, title = "subclusters_manual") +
  facet_wrap(~ geneset_name, scales = "free_y") +
  scale_fill_gradient(
    low = "gold", high = "navy", name = "FDR",
    guide = guide_colorbar(reverse = TRUE),
    limits = c(0, 0.05)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(file.path(sub_dir, "marker_enrichment_subclusters_manual.pdf"),
       plot = p, width = 12, height = 4, units = "in", dpi = 300, device = "pdf")
