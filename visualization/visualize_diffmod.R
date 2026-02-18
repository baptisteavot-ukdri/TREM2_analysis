###############################################################################
#   Astro (CV | R47H) stacked on top of Micro (CV | R47H)
#
# Features (same as your Astro plot):
# - Horizontal panels per cell type: CV | R47H
# - Y-axis = UNION of genes across CV+R47H (per cell type), shown only on LEFT panel
# - X-axis kept even if missing clusters after filtering (drop = FALSE)
# - ONE shared legend (collected) across ALL panels
# - PNG output only
###############################################################################

###############################################################################
# Publication-ready plots:
#  1) Astro:  CV | R47H | R62H   (fixed Astro1..Astro5 x-order)
#  2) Micro:  CV | R47H          (cluster order inferred from data)
#
# Features (same as the “initial Astro” style you requested):
# - Horizontal panels per variant
# - Y-axis = UNION of genes across variants (per cell type), shown only on LEFT
#   and aligned across panels (genes kept even if missing in a given variant)
# - X-axis categories kept even if empty (drop = FALSE)
# - ONE shared legend per plot (collected), global ranges across that plot
# - PNG only
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scico)
  library(patchwork)
  library(scales)
  library(stringr)
})

# -----------------------------
# INPUT PATHS
# -----------------------------
source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "config.R"))

paths_astro <- c(
  CV   = file.path(cfg$astro_markers_dir, "wilcox_subclusters_Astro_pb_CV.tsv"),
  R47H = file.path(cfg$astro_markers_dir, "wilcox_subclusters_Astro_pb_R47H.tsv"),
  R62H = file.path(cfg$astro_markers_dir, "wilcox_subclusters_Astro_pb_R62H.tsv")
)

paths_micro <- c(
  CV   = file.path(cfg$micro_markers_dir, "wilcox_subclusters_Micro_pb_CV.tsv"),
  R47H = file.path(cfg$micro_markers_dir, "wilcox_subclusters_Micro_pb_R47H.tsv")
)

out_png_astro <- file.path(cfg$astro_markers_dir, "astro_one_vs_all_CV_R47H_R62H_pubready_horizontal.png")
out_png_micro <- file.path(cfg$micro_markers_dir, "micro_one_vs_all_CV_R47H_pubready_horizontal.png")

# -----------------------------
# FILTERING THRESHOLDS
# -----------------------------
p_adj_thresh <- 0.05
lfc_thresh   <- 0.5

# Fixed Astro cluster order
astro_cluster_levels <- c("Astro1","Astro2","Astro3","Astro4","Astro5")

# -----------------------------
# IO helpers
# -----------------------------
read_tsv_variant <- function(path, variant){
  read.delim(path, sep = "\t", check.names = FALSE) |>
    as_tibble() |>
    mutate(
      variant      = variant,
      neglog10padj = -log10(p_val_adj)
    )
}

load_from_paths <- function(paths_named){
  purrr::imap_dfr(paths_named, read_tsv_variant)
}

# -----------------------------
# Plot builder: N horizontal panels for a cell type
# -----------------------------
build_variant_panel_plot <- function(df_raw,
                                     variants,
                                     title,
                                     cluster_levels = NULL,
                                     out_png,
                                     width_in = 28,
                                     height_in = 14) {
  
  # Keep only desired variants
  df_raw <- df_raw |>
    mutate(variant = factor(variant, levels = variants)) |>
    filter(!is.na(variant))
  
  # Filter significant rows
  df_sig <- df_raw |>
    filter(p_val_adj < p_adj_thresh, abs(avg_log2FC) > lfc_thresh)
  
  if (nrow(df_sig) == 0) {
    stop(title, ": no rows pass thresholds (p_val_adj < 0.05 & |avg_log2FC| > 0.5).")
  }
  
  # UNION genes across all included variants (post-filter) + stable order
  gene_order <- df_sig |>
    group_by(gene) |>
    summarise(score = mean(abs(avg_log2FC), na.rm = TRUE), .groups = "drop") |>
    arrange(score) |>
    pull(gene)
  
  # Cluster levels:
  # - if provided: enforce fixed order
  # - else infer and natural-sort (handles "Micro1..MicroN" nicely)
  if (is.null(cluster_levels)) {
    cluster_levels <- df_sig |>
      distinct(cluster) |>
      pull(cluster) |>
      as.character() |>
      str_sort(numeric = TRUE)
  }
  
  df_sig <- df_sig |>
    mutate(
      cluster = factor(cluster, levels = cluster_levels),
      gene    = factor(gene, levels = gene_order)
    )
  
  # Complete grid so every panel shares the same Y and X (empty slots allowed)
  df_plot <- df_sig |>
    select(variant, cluster, gene, avg_log2FC, p_val_adj, neglog10padj) |>
    tidyr::complete(
      variant = factor(variants, levels = variants),
      cluster = factor(cluster_levels, levels = cluster_levels),
      gene    = factor(gene_order, levels = gene_order)
    )
  
  # GLOBAL legend limits across the whole plot
  fill_lim <- range(df_sig$avg_log2FC, na.rm = TRUE)
  size_lim <- range(df_sig$neglog10padj, na.rm = TRUE)
  size_cap <- min(50, size_lim[2])
  size_breaks <- c(1, 2, 5, 10, 20, 50)
  size_breaks <- size_breaks[size_breaks <= size_cap]
  
  make_panel <- function(v, show_y = FALSE) {
    
    d <- df_plot |> filter(variant == v)
    
    p <- ggplot(d, aes(x = cluster, y = gene)) +
      geom_point(
        data = d |> filter(!is.na(avg_log2FC), !is.na(neglog10padj)),
        aes(size = pmin(neglog10padj, size_cap), fill = avg_log2FC),
        shape  = 21,
        colour = "black",
        stroke = 0.35
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      scale_fill_scico(
        palette   = "vik",
        direction = 1,
        midpoint  = 0,
        limits    = fill_lim,
        oob       = squish,
        name      = "avg_log2FC"
      ) +
      scale_size(
        limits = c(0, size_cap),
        range  = c(2.0, 10.0),
        breaks = size_breaks,
        name   = expression(-log[10]("adj p"))
      ) +
      guides(
        fill = guide_colorbar(
          title.position = "top",
          barwidth  = unit(18, "cm"),  # <-- make gradient longer
          barheight = unit(0.7, "cm")  # <-- slightly taller
        ),
        size = guide_legend(
          title.position = "top",
          nrow = 1
        )
      ) +
      theme_bw(base_size = 18) +
      theme(
        legend.position = "bottom",
        legend.box      = "vertical",          # stack fill + size nicely
        legend.spacing.y = unit(0.3, "cm"),     # space between legends
        legend.margin    = margin(t = 6, r = 6, b = 6, l = 6),
        legend.title     = element_text(size = 14),
        legend.text      = element_text(size = 12)
      )
    
    
    if (!show_y) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    p
  }
  
  panels <- purrr::map(seq_along(variants), function(i){
    v <- variants[i]
    make_panel(v, show_y = (i == 1))
  })
  
  p_combined <-
    patchwork::wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  p_combined <- p_combined +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
    )
  
  ggsave(out_png, p_combined, width = width_in, height = height_in, units = "in", dpi = 350, bg = "white")
  p_combined
}

# -----------------------------
# ASTRO PLOT (CV | R47H | R62H)
# -----------------------------
astro_raw <- load_from_paths(paths_astro)

p_astro <- build_variant_panel_plot(
  df_raw         = astro_raw,
  variants       = c("CV","R47H","R62H"),
  title          = "Wilcoxon one-vs-all (Astro): CV, R47H, R62H",
  cluster_levels = astro_cluster_levels,
  out_png        = out_png_astro,
  width_in       = 28,
  height_in      = 14
)

# -----------------------------
# MICRO PLOT (CV | R47H)
# -----------------------------
micro_raw <- load_from_paths(paths_micro)

p_micro <- build_variant_panel_plot(
  df_raw         = micro_raw,
  variants       = c("CV","R47H"),
  title          = "Wilcoxon one-vs-all (Micro): CV, R47H",
  cluster_levels = NULL,   # inferred from data
  out_png        = out_png_micro,
  width_in       = 22,
  height_in      = 14
)

# Print last plot (R prints p_micro by default if it's last)
p_micro
