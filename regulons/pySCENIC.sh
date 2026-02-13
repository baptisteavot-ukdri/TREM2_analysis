#!/usr/bin/env bash
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -N pySCENIC
#PBS -J 1-5

set -euo pipefail

cd "${PBS_O_WORKDIR:?PBS_O_WORKDIR is not set}"

# ------------------------------------------------------------------------------
# User-provided env vars (required)
# ------------------------------------------------------------------------------
: "${INPUT_DIR:?Set INPUT_DIR to the directory containing SCE .qs files}"
: "${LOGDIR:?Set LOGDIR to a writable log directory}"
: "${TABLE_DIR:?Set TABLE_DIR to a writable table output directory}"
: "${OUTDIR:?Set OUTDIR to a writable pySCENIC output directory}"
: "${N:?Set N to number of threads}"

input_dir="${INPUT_DIR%/}"
log_dir="${LOGDIR%/}"
table_dir="${TABLE_DIR%/}"
out_dir="${OUTDIR%/}"
n_threads="${N}"

resources_dir="/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/pySCENIC/resources"
preprocess_r="/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/pySCENIC/R/preprocess_input_files.R"

ranking_feather="${resources_dir}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
tfs_list="${resources_dir}/allTFs_hg38.txt"
motif_tbl="${resources_dir}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# ------------------------------------------------------------------------------
# Pick file for this array index
# ------------------------------------------------------------------------------
mapfile -t files < <(find "${input_dir}" -maxdepth 1 -type f -name "*.qs" -printf "%f\n" | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "No .qs files found in: ${input_dir}" >&2
  exit 1
  fi
  
  idx=$((PBS_ARRAY_INDEX - 1))
  if (( idx < 0 || idx >= ${#files[@]} )); then
    echo "PBS_ARRAY_INDEX=${PBS_ARRAY_INDEX} out of range (1..${#files[@]})" >&2
    exit 1
    fi
    
    input_sce="${files[$idx]}"
    input_path="${input_dir}/${input_sce}"
    
    subcluster="$(basename "${input_sce}")"
    subcluster="${subcluster%%_*}"  # take string before first underscore
    
    # ------------------------------------------------------------------------------
    # Create output dirs
    # ------------------------------------------------------------------------------
    mkdir -p "${log_dir}" "${table_dir}" "${out_dir}"
    
    echo "Job started: $(date)"
    echo "Input SCE:   ${input_path}"
    echo "Subcluster:  ${subcluster}"
    echo "Threads:     ${n_threads}"
    
    # ------------------------------------------------------------------------------
    # Preprocess (R) inside conda env
    # ------------------------------------------------------------------------------
    eval "$("${HOME}/anaconda3/bin/conda" shell.bash hook)"
    conda activate preprocess_pySCENIC
    
    Rscript "${preprocess_r}" \
    --sce_file "${input_path}" \
    --feather_file "${ranking_feather}" \
    --celltype "${subcluster}" \
    --outdir "${table_dir}" \
    &> "${log_dir}/${subcluster}.preprocessing.out"
    
    conda deactivate
    
    expr_tsv="${table_dir}/${subcluster}.tsv"
    adj_tsv="${out_dir}/${subcluster}.adjacencies.tsv"
    reg_csv="${out_dir}/${subcluster}.regulons.csv"
    auc_csv="${out_dir}/${subcluster}.auc_mtx.csv"
    
    # ------------------------------------------------------------------------------
    # pySCENIC via Singularity (docker image)
    # ------------------------------------------------------------------------------
    img="docker://aertslab/pyscenic:0.12.0"
    
    singularity run "${img}" pyscenic grn \
    --num_workers "${n_threads}" \
    "${expr_tsv}" \
    "${tfs_list}" \
    -o "${adj_tsv}" \
    &> "${log_dir}/${subcluster}.grn.out"
    
    singularity run "${img}" pyscenic ctx \
    "${adj_tsv}" \
    "${ranking_feather}" \
    --annotations_fname "${motif_tbl}" \
    --expression_mtx_fname "${expr_tsv}" \
    --mode "custom_multiprocessing" \
    --output "${reg_csv}" \
    --num_workers "${n_threads}" \
    &> "${log_dir}/${subcluster}.ctx.out"
    
    singularity run "${img}" pyscenic aucell \
    "${expr_tsv}" \
    "${reg_csv}" \
    -o "${auc_csv}" \
    --num_workers "${n_threads}" \
    &> "${log_dir}/${subcluster}.aucell.out"
    
    echo "Job ended: $(date)"
    