suppressPackageStartupMessages({
  library(tidyverse)
  library(DirichletReg)
  library(Seurat)
  library(qs)
  library(paletteer)
  library(ggpubr)
  library(stringr)
})

# -----------------------------
# Paths + inputs
# -----------------------------
source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "config.R"))

work_dir <- cfg$dirichlet_work_dir
setwd(work_dir)

source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "dirichlet/additional_scripts/get_celltype_freq.r"))
source(file.path(Sys.getenv("TREM2_ANALYSIS_ROOT"), "dirichlet/additional_scripts/process_dirichlet_fit.r"))

celltype_label <- "Micro"
seu_path <- cfg$seu_micro

outdir <- file.path(sprintf("subclustering_%s_full_cohort", celltype_label), "dirichlet")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "table"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "plot"),  recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helpers
# -----------------------------
make_sig_labels <- function(p) {
  padj <- p.adjust(p, method = "fdr")
  tibble(
    padj = padj,
    label = case_when(
      is.na(padj)      ~ " ",
      padj <= 0.001    ~ "***",
      padj <= 0.01     ~ "**",
      padj <= 0.05     ~ "*",
      TRUE             ~ " "
    )
  )
}

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

plot_props_with_pvals <- function(
    sn_ct_prop,
    fit_res,
    col_var,
    facet_var = NULL,
    out_png,
    title_text = NULL,
    ncol_facet = 1,
    height = 6,
    width = 8,
    dpi = 300
) {
  ctype <- sum(vapply(sn_ct_prop, is.numeric, logical(1))) - 0  # numeric columns are celltype proportions
  stopifnot(ctype > 1)
  
  sn_ct_prop_long <- sn_ct_prop %>%
    pivot_longer(cols = 1:ctype, names_to = "celltype", values_to = "ct_prop")
  
  # y-position per celltype (+ facet if needed)
  if (is.null(facet_var)) {
    y_dt <- sn_ct_prop_long %>%
      group_by(celltype) %>%
      summarise(y.position = max(ct_prop * 100, na.rm = TRUE) + 3, .groups = "drop")
  } else {
    y_dt <- sn_ct_prop_long %>%
      group_by(celltype, .data[[facet_var]]) %>%
      summarise(y.position = max(ct_prop * 100, na.rm = TRUE) + 3, .groups = "drop") %>%
      rename(!!facet_var := .data[[facet_var]])
  }
  
  # Build annotation table for ggpubr
  ann_text <- fit_res %>%
    filter(str_detect(term, fixed(col_var))) %>%
    mutate(.row = row_number()) %>%
    bind_cols(make_sig_labels(.$pval)) %>%
    select(-.row) %>%
    mutate(
      group1 = levels(factor(sn_ct_prop_long[[col_var]]))[1],
      group2 = str_replace(term, fixed(col_var), ""),
      x = as.numeric(factor(celltype)),
      xmin = if_else(label == " ", NA_real_, x - 0.2),
      xmax = x + 0.2
    ) %>%
    left_join(y_dt, by = c("celltype", if (!is.null(facet_var)) facet_var))
  
  # Palette
  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
  dodge <- position_dodge(width = 0.8)
  
  p <- ggplot(sn_ct_prop_long, aes(x = celltype, y = ct_prop * 100, fill = .data[[col_var]])) +
    geom_bar(stat = "summary", fun = "mean", position = dodge, width = 0.75) +
    geom_errorbar(stat = "summary", position = dodge, width = 0.25) +
    geom_point(position = dodge, size = 0.2, alpha = 0.7) +
    scale_color_manual(values = palette_choice[c(2, 1)], aesthetics = c("colour", "fill"), name = NULL) +
    labs(x = NULL, y = "% Cell type", title = title_text %||% fit_res$model[1]) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title  = element_text(size = 14, color = "black"),
      legend.text = element_text(size = 12, color = "black"),
      strip.text  = element_text(size = 12, color = "black"),
      plot.title  = element_text(size = 11, color = "black")
    ) +
    ggpubr::stat_pvalue_manual(
      ann_text,
      label = "label",
      tip.length = 0.01,
      hide.ns = TRUE
    )
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(vars(.data[[facet_var]]), ncol = ncol_facet, scales = "free_x")
  }
  
  ggsave(filename = out_png, plot = p, height = height, width = width, units = "in", dpi = dpi)
  invisible(p)
}

fit_dirichlet <- function(sn_ct_prop, ctype, col_var, confounders) {
  dt <- DirichletReg::DR_data(sn_ct_prop[, 1:ctype])
  form <- as.formula(sprintf("dt ~ %s", paste(c(col_var, confounders), collapse = " + ")))
  fit <- DirichletReg::DirichReg(formula = form, data = sn_ct_prop)
  fit
}

# -----------------------------
# Load Seurat + filter
# -----------------------------
seu <- qs::qread(seu_path)

# Exclude samples without CD33 values
seu <- subset(seu, subset = !is.na(CD33Group))
seu$manifest <- droplevels(seu$manifest)

# Optional: exclude low-nuclei samples (was previously commented out)
sample_pc <- table(seu$manifest)
sample_outlier <- unname(quantile(sample_pc, probs = 0.01))  # bottom 1%
keep_manifests <- names(sample_pc)[sample_pc > sample_outlier]
# seu <- subset(seu, subset = manifest %in% keep_manifests)
# seu$manifest <- droplevels(seu$manifest)

# -----------------------------
# Cell-type proportions per manifest
# -----------------------------
res <- get_celltype_freq(
  metadata = seu@meta.data,
  unique_id_var = "manifest",
  celltype_var  = "subclusters_label"
)

sn_ct_prop <- as.data.frame(res$prop_counts_mat)
sn_ct_prop <- sn_ct_prop[, order(colnames(sn_ct_prop)), drop = FALSE]
sn_ct_prop$manifest <- rownames(sn_ct_prop)

ctype <- length(unique(seu$subclusters_label))

metadata <- seu@meta.data %>%
  select(
    sample, SampleID, StudyID, Age, Sex, APOE, TREM2Variant,
    BrainRegion, Braak, NeuropathologicalDiagnosis,
    PostMortemInterval, BrainBankNetworkID, BrainBankNetworkIDFormatted,
    CD33Group, APOEgroup,
    pctAT8PositiveArea, pctPHF1PositiveArea, pct4G8PositiveArea,
    manifest
  ) %>%
  distinct() %>%
  mutate(
    NeuropathologicalDiagnosis = factor(NeuropathologicalDiagnosis, levels = c("Control", "AD")),
    Age = as.vector(scale(Age)),
    PostMortemInterval = as.vector(scale(PostMortemInterval))
  )

sn_ct_prop <- sn_ct_prop %>%
  select(1:(ctype + 1)) %>%
  left_join(metadata, by = "manifest")

# -----------------------------
# 1) Full-cohort model: Diagnosis effect
# -----------------------------
col_var <- "NeuropathologicalDiagnosis"
confounders <- c("Sex", "Age", "PostMortemInterval", "BrainRegion", "TREM2Variant", "APOEgroup", "CD33Group")

dirichlet_fit <- fit_dirichlet(sn_ct_prop = sn_ct_prop, ctype = ctype, col_var = col_var, confounders = confounders)
fit_res <- process_dirichlet_fit(dirichlet_fit, col_var = col_var)

write_tsv(
  fit_res,
  file.path(outdir, "table", sprintf("%s_dirichlet_pval.tsv", col_var))
)

plot_props_with_pvals(
  sn_ct_prop = sn_ct_prop,
  fit_res = fit_res,
  col_var = col_var,
  facet_var = NULL,
  out_png = file.path(outdir, "plot", sprintf("%s.png", col_var)),
  title_text = str_wrap(fit_res$model[1], width = 60),
  height = 5,
  width = 8
)

# -----------------------------
# 2) Stratified analyses: Diagnosis within each level of facet_var
# -----------------------------
facet_vars <- c("TREM2Variant", "BrainRegion", "APOEgroup", "CD33Group")

for (facet_var in facet_vars) {
  
  message("Stratified by: ", facet_var)
  
  confounders_base <- c("Sex", "Age", "PostMortemInterval", "TREM2Variant", "BrainRegion", "APOEgroup", "CD33Group")
  confounders_strat <- setdiff(confounders_base, facet_var)
  
  model_formula <- as.formula(sprintf(
    "count ~ %s + %s",
    col_var,
    paste(confounders_strat, collapse = " + ")
  ))
  
  fit_list <- vector("list", length = length(unique(sn_ct_prop[[facet_var]])))
  names(fit_list) <- as.character(unique(sn_ct_prop[[facet_var]]))
  
  for (lvl in names(fit_list)) {
    dt <- sn_ct_prop %>%
      filter(.data[[facet_var]] %in% lvl) %>%
      mutate(!!facet_var := droplevels(.data[[facet_var]]))
    
    count <- DirichletReg::DR_data(dt[, 1:ctype])
    
    fit <- DirichletReg::DirichReg(formula = model_formula, data = dt)
    fit_tbl <- process_dirichlet_fit(fit, col_var = col_var)
    fit_tbl[[facet_var]] <- lvl
    
    fit_list[[lvl]] <- fit_tbl
  }
  
  fit_res_strat <- bind_rows(fit_list)
  
  write_tsv(
    fit_res_strat,
    file.path(outdir, "table", sprintf("%s_dirichlet_pval_stratified_by_%s.tsv", col_var, facet_var))
  )
  
  plot_props_with_pvals(
    sn_ct_prop = sn_ct_prop,
    fit_res = fit_res_strat,
    col_var = col_var,
    facet_var = facet_var,
    out_png = file.path(outdir, "plot", sprintf("%s_stratified_by_%s.png", col_var, facet_var)),
    title_text = str_wrap(fit_res_strat$model[1], width = 60),
    ncol_facet = 1,
    height = 10,
    width = 8
  )
  
  prop_means <- sn_ct_prop %>%
    pivot_longer(cols = 1:ctype, names_to = "celltype", values_to = "ct_prop") %>%
    group_by(celltype, .data[[facet_var]], .data[[col_var]]) %>%
    summarise(Mean = round(mean(ct_prop * 100, na.rm = TRUE)), .groups = "drop")
  
  write_tsv(
    prop_means,
    file.path(outdir, "table", sprintf("%s_mean_prop_stratified_by_%s.tsv", col_var, facet_var))
  )
}

# -----------------------------
# 3) Stratified by TREM2Variant: only selected variants
# -----------------------------
facet_var <- "TREM2Variant"
variants_keep <- c("CV", "R47H", "R62H")

sn_ct_prop_sub <- sn_ct_prop %>%
  filter(.data[[facet_var]] %in% variants_keep) %>%
  mutate(!!facet_var := droplevels(.data[[facet_var]]))

confounders_tv <- c("Sex", "Age", "PostMortemInterval", "BrainRegion", "APOEgroup", "CD33Group")
model_formula <- as.formula(sprintf("count ~ %s + %s", col_var, paste(confounders_tv, collapse = " + ")))

fit_list <- vector("list", length = length(variants_keep))
names(fit_list) <- variants_keep

for (lvl in variants_keep) {
  dt <- sn_ct_prop_sub %>%
    filter(.data[[facet_var]] %in% lvl) %>%
    mutate(!!facet_var := droplevels(.data[[facet_var]]))
  
  count <- DirichletReg::DR_data(dt[, 1:ctype])
  fit <- DirichletReg::DirichReg(formula = model_formula, data = dt)
  
  fit_tbl <- process_dirichlet_fit(fit, col_var = col_var)
  fit_tbl[[facet_var]] <- lvl
  fit_list[[lvl]] <- fit_tbl
}

fit_res_tv <- bind_rows(fit_list)

write_tsv(
  fit_res_tv,
  file.path(outdir, "table", sprintf("%s_dirichlet_pval_stratified_by_%s.tsv", col_var, facet_var))
)

plot_props_with_pvals(
  sn_ct_prop = sn_ct_prop_sub,
  fit_res = fit_res_tv,
  col_var = col_var,
  facet_var = facet_var,
  out_png = file.path(outdir, "plot", sprintf("%s_stratified_by_%s.png", col_var, facet_var)),
  title_text = str_wrap(fit_res_tv$model[1], width = 60),
  ncol_facet = 1,
  height = 10,
  width = 8
)

prop_means_tv <- sn_ct_prop_sub %>%
  pivot_longer(cols = 1:ctype, names_to = "celltype", values_to = "ct_prop") %>%
  group_by(celltype, .data[[facet_var]], .data[[col_var]]) %>%
  summarise(Mean = round(mean(ct_prop * 100, na.rm = TRUE)), .groups = "drop")

write_tsv(
  prop_means_tv,
  file.path(outdir, "table", sprintf("%s_mean_prop_stratified_by_%s.tsv", col_var, facet_var))
)

# -----------------------------
# 4) Compare TREM2Variant within AD only
# -----------------------------
facet_var <- "NeuropathologicalDiagnosis"
col_var2 <- "TREM2Variant"

sn_ct_prop_ad <- sn_ct_prop %>%
  filter(.data[[facet_var]] %in% "AD") %>%
  mutate(!!facet_var := droplevels(.data[[facet_var]]))

confounders_ad <- c("Sex", "Age", "PostMortemInterval", "BrainRegion", "APOEgroup", "CD33Group")
dirichlet_fit_ad <- fit_dirichlet(
  sn_ct_prop = sn_ct_prop_ad,
  ctype = ctype,
  col_var = col_var2,
  confounders = confounders_ad
)

fit_res_ad <- process_dirichlet_fit(dirichlet_fit_ad, col_var = col_var2)
fit_res_ad[[facet_var]] <- "AD"

write_tsv(
  fit_res_ad,
  file.path(outdir, "table", sprintf("%s_dirichlet_pval_stratified_by_%s.tsv", col_var2, facet_var))
)

# Plot within AD
plot_props_with_pvals(
  sn_ct_prop = sn_ct_prop_ad,
  fit_res = fit_res_ad,
  col_var = col_var2,
  facet_var = facet_var,
  out_png = file.path(outdir, "plot", sprintf("%s_stratified_by_%s.png", col_var2, facet_var)),
  title_text = str_wrap(fit_res_ad$model[1], width = 60),
  ncol_facet = 1,
  height = 5,
  width = 8
)

# -----------------------------
# 5) Convenience plot: stratify by CD33Group only (no p-value annotation)
# -----------------------------
facet_var <- "CD33Group"
col_var <- "NeuropathologicalDiagnosis"

ctype_long <- sn_ct_prop %>%
  pivot_longer(cols = 1:ctype, names_to = "celltype", values_to = "ct_prop")

palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")
dodge <- position_dodge(width = 0.8)

p_cd33 <- ggplot(ctype_long, aes(x = celltype, y = ct_prop * 100, fill = .data[[col_var]])) +
  geom_bar(stat = "summary", fun = "mean", position = dodge, width = 0.75) +
  geom_errorbar(stat = "summary", position = dodge, width = 0.25) +
  geom_point(position = dodge, size = 0.2, alpha = 0.7) +
  scale_color_manual(values = palette_choice[c(2, 1)], aesthetics = c("colour", "fill"), name = NULL) +
  labs(x = NULL, y = "% Cell type", title = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    strip.text  = element_text(size = 12, color = "black")
  ) +
  facet_wrap(vars(.data[[facet_var]]), scales = "free_x", ncol = 1)

ggsave(
  filename = file.path(outdir, "plot", sprintf("%s_stratified_by_%s.png", col_var, facet_var)),
  plot = p_cd33,
  height = 15, width = 8, units = "in", dpi = 300
)
