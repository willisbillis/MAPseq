#!/bin/bash
#
# run_mapseq_aggr.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs the MAPseq aggregation pipeline.

# NOTICE: At this point, the user has run through all steps in at least two
#       RNA and ATAC runs and is ready to aggregate the runs for a single 
#       project. Additionally, the user has set and checked all variables
#       in the project_config.txt and populated the sample hashtag
#       sheets with the correct hashtags for demultiplexing the pools.
################################################################################
set -euo pipefail

# Import all the global variables for this project
source ./project_config.txt
# Export global variables for R scripts
export PROJECT_PATH=$PROJECT_PATH
export PROJECT_NAME=$PROJECT_NAME
export GEX_NAMING_ID=$GEX_NAMING_ID
export GEX_FEAT_NAMING_ID=$GEX_FEAT_NAMING_ID
export ATAC_NAMING_ID=$ATAC_NAMING_ID
export ASAP_NAMING_ID=$ASAP_NAMING_ID
################################################################################
AGGR_DIR="$PROJECT_PATH/$PROJECT_NAME"
RNA_DIR="$AGGR_DIR/pipeline/RNA.FB.VDJ"
ATAC_DIR="$AGGR_DIR/pipeline/ATAC.ASAP"

run_step() {
    local step_name=$1
    local working_dir=$2
    shift 2
    echo "$(date) Running ${step_name}..."
    if (cd "$working_dir" && "$@"); then
        return 0
    else
        local status=$?
        echo "[ERROR] ${step_name} failed with status ${status}" >&2
        exit "${status:-1}"
    fi
}

has_manifest_samples() {
    awk -F, 'NR > 1 { gsub(/\r/, "", $1); if ($1 != "") found = 1 }
        END { exit !found }' "$1"
}

RNA_MANIFEST="$RNA_DIR/hashtag_ref_rna.csv"
ATAC_MANIFEST="$ATAC_DIR/hashtag_ref_atac.csv"
[[ -r $RNA_MANIFEST ]] || { echo "[ERROR] Cannot read $RNA_MANIFEST" >&2; exit 1; }
[[ -r $ATAC_MANIFEST ]] || { echo "[ERROR] Cannot read $ATAC_MANIFEST" >&2; exit 1; }

METRICS_PATH="$PROJECT_PATH/reports/sample_metrics.csv"
LOCK_PATH="$METRICS_PATH.lock"
REPORTS_DIR=$(dirname "$METRICS_PATH")
if [[ -L $REPORTS_DIR ]]; then
    echo "[ERROR] Metrics reports directory must not be a symlink: $REPORTS_DIR" >&2
    exit 1
fi
mkdir -p "$REPORTS_DIR"
lock_started=$SECONDS
until mkdir "$LOCK_PATH" 2>/dev/null; do
    if (( SECONDS - lock_started >= 300 )); then
        echo "[ERROR] Timed out waiting for metrics lock: $LOCK_PATH" >&2
        exit 1
    fi
    sleep 0.1
done
release_metrics_lock() { rmdir "$LOCK_PATH"; }
trap release_metrics_lock EXIT

if [[ -L $METRICS_PATH ]]; then
    echo "[ERROR] Metrics output must not be a symlink: $METRICS_PATH" >&2
    exit 1
fi
if [[ -e $METRICS_PATH ]]; then
    mv -f -- "$METRICS_PATH" "$METRICS_PATH.stale"
fi

# Phase-2 QC outputs describe the previous aggregation and must not survive reset.
for analysis_dir in "$AGGR_DIR/analysis/RNA.FB.VDJ" "$AGGR_DIR/analysis/ATAC.ASAP"; do
    [[ -d $analysis_dir ]] || continue
    find "$analysis_dir" -type f \( -name 'qc_rna_*.RDS' -o -name 'qc_atac_*.RDS' -o -name 'QC.rna_sampleID_filtering.stats.csv' -o -name 'QC.atac_sampleID_filtering.stats.csv' -o -name 'scatter_rna_*.png' -o -name 'scatter_atac_*.png' -o -name 'clustree_rna_*.png' -o -name 'clustree_atac_*.png' -o -name 'umap_rna_*.pdf' -o -name 'umap_atac_*.pdf' -o -name 'pca_variance_explained.csv' -o -name 'correlation_rna_filtered_cells_vs_depth.*' -o -name 'correlation_atac_filtered_cells_vs_depth.*' -o -name 'DEG_RNA_clusters.csv' -o -name 'DAR_ATAC_clusters.csv' -o -name 'lsi_variance_explained.csv' \) -delete
done

rm -f "$RNA_DIR/aggr.csv" "$RNA_DIR/summary_metrics_rna.csv"
rm -f "$RNA_DIR"/demux_metrics.*.csv
rm -f "$ATAC_DIR/ATAC/aggr.csv" "$ATAC_DIR/ATAC/summary_metrics_atac.csv"
rm -f "$ATAC_DIR"/demux_metrics.*.csv "$ATAC_DIR"/library_stats.*.csv
rm -f "$AGGR_DIR"/analysis/RNA.FB.VDJ/QC.*_sampleID_filtering.stats.csv
rm -f "$AGGR_DIR"/analysis/ATAC.ASAP/QC.*_sampleID_filtering.stats.csv
rm -f "$AGGR_DIR/analysis/RNA.FB.VDJ/data/raw_rna.hto.adt_${PROJECT_NAME}.RDS"
rm -f "$AGGR_DIR/analysis/ATAC.ASAP/data/raw_atac.hto_${PROJECT_NAME}.RDS"

if has_manifest_samples "$RNA_MANIFEST"; then
    run_step "RNA Cell Ranger aggregation" "$RNA_DIR" ./run_cellranger.aggr_rna.sh
    run_step "RNA souporcell demultiplexing" "$RNA_DIR" \
        conda run -n souporcell --live-stream ./run_souporcell.demux_rna.sh
    run_step "RNA Seurat demultiplexing" "$RNA_DIR" Rscript run_seurat.demux_rna.R
fi

if has_manifest_samples "$ATAC_MANIFEST"; then
    run_step "ATAC Cell Ranger aggregation" "$ATAC_DIR" ./run_cellranger.aggr_atac.sh
    run_step "ATAC souporcell demultiplexing" "$ATAC_DIR" \
        conda run -n souporcell --live-stream ./run_souporcell.demux_atac.sh
    run_step "ATAC Seurat demultiplexing" "$ATAC_DIR" Rscript run_seurat.demux_atac.R
fi

echo "Compiling available mapping, demultiplexing, and QC metrics..."
python3 "$AGGR_DIR/phase2_scripts/compile_pipeline_metrics.py" "$PROJECT_PATH" \
    --project-dir "$AGGR_DIR" --project-id "$PROJECT_NAME" --lock-held
rm -f -- "$METRICS_PATH.stale"