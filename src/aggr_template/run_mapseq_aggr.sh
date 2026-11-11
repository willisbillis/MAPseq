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
    --project-dir "$AGGR_DIR" --project-id "$PROJECT_NAME"