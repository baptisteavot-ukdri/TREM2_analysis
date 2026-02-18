# TREM2_analysis

Bioinformatics pipeline for investigating the effects of **TREM2** and **ABCA7** genetic variants on microglial and astrocytic transcriptional regulation in post-mortem brain tissue from Alzheimer's disease patients.

**Institute:** UK Dementia Research Institute at Imperial College London
**PI:** Prof. Paul Matthews (OBE, MD, DPhil, FMedSci)
**Data:** snRNA-seq and snATAC-seq — MAP (Multimodal Alzheimer's disease Profiling) cohort

---

## Overview

Four complementary analytical modules:

| Module | Description |
|--------|-------------|
| **Subclustering** | Harmony-integrated microglial subclustering (2 rounds) |
| **DEG** | Pseudobulk differential expression with limma-trend |
| **Regulons** | pySCENIC gene regulatory network inference + differential modularity |
| **Dirichlet** | Cell type proportion modelling across conditions |

---

## Workflow

```
snRNA-seq SCE objects (post-QC, per cell type)
         │
         ├─── SUBCLUSTERING ──────────────────────────────────────────────────┐
         │    Round 1: Harmony integration, UMAP, Seurat clustering (res=0.5) │
         │    Round 2: Micro-only, resolution grid (0.1–0.5), manual labels   │
         │    Output: seu.qs with subclusters_label                           │
         │                                                                    │
         ├─── DEG ────────────────────────────────────────────────────────────┤
         │    Pseudobulk (mean size=30), limma-trend, random effect on sample │
         │    Contrasts: AD vs Control; TREM2 variants; stratified analyses   │
         │    Output: {celltype}_{contrast}.tsv (logFC, padj, CI)             │
         │                                                                    │
         ├─── REGULONS ───────────────────────────────────────────────────────┤
         │    pySCENIC: GRN → motif annotation (CTX) → AUCell scoring         │
         │    Differential modularity (dream): regulon AUC ~ pathology metric │
         │    GO enrichment on regulon target genes                            │
         │    Output: AUC matrix, dream results, GO enrichment tables          │
         │                                                                    │
         └─── DIRICHLET ──────────────────────────────────────────────────────┘
              Cell type proportion ~ Diagnosis + covariates (DirichletReg)
              Stratified by: TREM2Variant, BrainRegion, APOEgroup, CD33Group
              Output: proportion tables, p-value TSVs, bar plots
```

---

## Repository Structure

```
TREM2_analysis/
├── config.R                           ← Central path & parameter configuration
├── README.md
│
├── DEG/
│   └── dge_metacell.R                 ← CLI: pseudobulk + limma-trend DEG
│
├── dirichlet/
│   ├── dirichlet_micro.R              ← Dirichlet regression on cell proportions
│   └── additional_scripts/
│       ├── get_celltype_freq.r        ← Helper: tally cell type frequencies
│       └── process_dirichlet_fit.r    ← Helper: extract Dirichlet coefficients
│
├── regulons/
│   ├── pySCENIC.sh                    ← PBS job: GRN → CTX → AUCell
│   ├── preprocess_input_files.R       ← CLI: prepare expression matrix for pySCENIC
│   ├── Pathways/
│   │   └── regulons_query_go.R        ← GO enrichment on regulon target genes
│   ├── diffmod/
│   │   └── diffmod_micro.R            ← dream: regulon AUC vs pathology metrics
│   └── resources/
│       ├── allTFs_hg38.txt
│       ├── motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
│       └── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
│
├── subclustering/
│   ├── Micro/
│   │   ├── micro_subclustering_round1.R
│   │   └── micro_subclustering_round2.R
│   └── additional_scripts/
│       ├── enrichment_custom.r
│       ├── filter_genes_expressed_seu.r
│       └── plot_marker_heatmap.r
│
└── visualization/
    └── visualize_diffmod.R            ← Publication-ready regulon heatmaps
```

---

## Setup

### 1. Set environment variables

All scripts load paths from `config.R`, which reads from environment variables. Set these before running:

```bash
# Required: path to this repository
export TREM2_ANALYSIS_ROOT=/path/to/TREM2_analysis

# Optional: override default base paths
export TREM2_PROJECT_RDS=/rds/general/project/ukdrmultiomicsproject/live
export TREM2_USER_RDS=/rds/general/user/bavot/home
export TREM2_MNT_DATA=/mnt/data/bavot/TREM2
```

Add these to your `~/.bashrc`, PBS job header (`#PBS -v`), or SLURM script (`#SBATCH --export`).

### 2. R environment

**R >= 4.2** required. Install key packages:

```r
# Bioconductor
BiocManager::install(c(
  "SingleCellExperiment", "SummarizedExperiment",
  "edgeR", "limma", "variancePartition",
  "clusterProfiler", "org.Hs.eg.db"
))

# CRAN
install.packages(c(
  "Seurat", "harmony", "qs", "argparse", "DirichletReg",
  "tidyverse", "patchwork", "ComplexHeatmap", "scico",
  "vroom", "scCustomize", "paletteer", "ggpubr", "BiocParallel"
))
```

### 3. pySCENIC environment

```bash
conda create -n pyscenic python=3.8
conda activate pyscenic
pip install pyscenic==0.12.0
```

Or use the Singularity container referenced in `regulons/pySCENIC.sh` (`aertslab/pyscenic:0.12.0`).

### 4. External MAP pipeline scripts

Several scripts depend on internal MAP pipeline functions sourced via `config.R`:

| Function | Source script | Used in |
|----------|---------------|---------|
| `.sconline.PseudobulkGeneration()` | `generate_pseudobulk_limma_trend_deg_gazestani.r` | `DEG/dge_metacell.R` |
| `.sconline.fitLimmaFn()` | same as above | `DEG/dge_metacell.R` |
| `plot_marker_heatmap()` | `plot_marker_heatmap.r` | `subclustering/Micro/micro_subclustering_round1.R` |
| `map_celltypes_seu()` | `map_celltype_seu.r` | `subclustering/Micro/micro_subclustering_round1.R` |

---

## Running the Pipeline

### Step 1 — Subclustering

```bash
Rscript subclustering/Micro/micro_subclustering_round1.R
Rscript subclustering/Micro/micro_subclustering_round2.R
```

Round 1 produces a coarsely clustered Seurat object with EWCE-based cell type labels.
Round 2 refines microglia into 6 subclusters: `HM`, `DAM_HLA`, `CRM`, `DAM`, `PVM`, `ProliferativeMicro`.

### Step 2 — Differential Expression (DEG)

The DEG script is CLI-driven via `argparse`:

```bash
Rscript DEG/dge_metacell.R \
  --sce /path/to/Micro_sce.qs \
  --dependent_var NeuropathologicalDiagnosis \
  --ref_class Control \
  --confounding_vars "Sex,Age,PostMortemInterval,BrainRegion,TREM2Variant,APOEgroup,CD33Group" \
  --output_dir /path/to/output \
  --mod_name full_cohort
```

Add `--stratification_var TREM2Variant` to run within each variant group separately.

### Step 3 — pySCENIC (regulon inference)

```bash
# Preprocess input expression matrix
Rscript regulons/preprocess_input_files.R \
  --sce_file /path/to/Micro_sce.qs \
  --feather_file regulons/resources/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --celltype Micro \
  --outdir /path/to/pyscenic_input

# Submit pySCENIC job array (PBS)
qsub -v INPUT_DIR=/path/to/pyscenic_input,OUTDIR=/path/to/pyscenic_out,N=8 \
  regulons/pySCENIC.sh
```

### Step 4 — GO Enrichment on Regulons

```bash
Rscript regulons/Pathways/regulons_query_go.R
```

### Step 5 — Differential Modularity

```bash
Rscript regulons/diffmod/diffmod_micro.R
```

Associates regulon AUC scores with pathology metrics (e.g. `Total4G8Density`) using `dream`, stratified by `TREM2Variant`, `APOEgroup`, and `CD33Group`.

### Step 6 — Dirichlet Regression

```bash
Rscript dirichlet/dirichlet_micro.R
```

### Step 7 — Visualization

```bash
Rscript visualization/visualize_diffmod.R
```

---

## Key Outputs

| Analysis | Output | Location |
|----------|--------|----------|
| Subclustering R1 | `seu_relabelled.qs`, marker TSV, plots | `subclustering_Micro_full_cohort/subclusters_round1/` |
| Subclustering R2 | `seu.qs`, marker TSVs, enrichment plots | `subclustering_Micro_full_cohort/subclusters_round2/` |
| DEG | `{celltype}_{contrast}.tsv` | `de_{contrast}_{model}/` |
| pySCENIC | `Micro.regulons.filtered.csv`, `Micro.auc_mtx.csv` | `pySCENIC/out/` |
| GO enrichment | `TF_GO_enrichment_BP_MF_CC_top10.csv` | `TREM2/Micro/regulons/` |
| Differential modularity | `res_dream_{subcluster}_{layer}.tsv` | `diffmod_subcelltypes/Total4G8Density/` |
| Dirichlet | `*_dirichlet_pval.tsv`, `*.png` | `subclustering_Micro_full_cohort/dirichlet/` |
| Visualization | `*_pubready_horizontal.png` | regulon marker directories |

---

## Metadata Columns

| Column | Type | Description |
|--------|------|-------------|
| `manifest` | string | Unique sample/donor ID |
| `NeuropathologicalDiagnosis` | factor | `Control` / `AD` |
| `TREM2Variant` | factor | `CV`, `R47H`, `R62H`, others |
| `APOEgroup` | factor | `APOE4-neg` / `APOE4-pos` |
| `CD33Group` | factor | `CV` / `CD33var` |
| `BrainRegion` | factor | Brain region of dissection |
| `Age`, `PostMortemInterval` | numeric | Donor demographics |
| `Sex` | factor | `M` / `F` |
| `Total4G8Density` | numeric | Amyloid pathology metric (4G8 immunoreactivity) |
| `subclusters_label` | factor | Microglial subcluster identity |

---

## Microglial Subclusters

| Label | Identity |
|-------|----------|
| `HM` | Homeostatic Microglia |
| `DAM` | Disease-Associated Microglia |
| `DAM_HLA` | DAM with high HLA expression |
| `CRM` | Cycling/Reactive Microglia |
| `PVM` | Perivascular Macrophages |
| `ProliferativeMicro` | Proliferative Microglia |

---

## Citation

> Avot B. et al. (in preparation). Effects of TREM2 variants on glial transcriptional regulation in Alzheimer's disease. UK Dementia Research Institute at Imperial College London.

---

## Contact

Baptiste Avot · [LinkedIn](https://www.linkedin.com/in/baptiste-avot) · [Google Scholar](https://scholar.google.com/citations?user=h878R24AAAAJ&hl=en) · [GitHub](https://github.com/baptisteavot-ukdri)
