#!/bin/bash
#
# run_cellranger.aggr_atac.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger-atac aggr pipeline formatted for the ATAC + ASAP
# side of the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has run through all steps in at least two
#     ATAC + ASAP runs and is ready to aggregate the runs for a single project.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
SAMPLES_ARRAY=($(ls -d $PROJECT_PATH/*/pipeline/ATAC.ASAP/ATAC/*/))
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/ATAC
OUTPUT_FILE=$OUTPUT_DIR/cellranger_atac_aggr.log
AGGR_CSV=$OUTPUT_DIR/aggr.csv
METRICS_CSV=$OUTPUT_DIR/summary_metrics_atac.csv
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

CR_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
echo "$(date) Running Cell Ranger ATAC version $CR_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE

## Create aggregation csv
printf '%s\n' library_id fragments cells | paste -sd ',' >> $AGGR_CSV
for sample_path in "${SAMPLES_ARRAY[@]}"; do
  if [[ $(basename $sample_path) != reports ]]; then
    sample_name=$(basename $sample_path)
    fragments_file=${sample_path}outs/fragments.tsv.gz
    cells_file=${sample_path}outs/singlecell.csv

    if [[ -f $fragments_file && -f $cells_file ]]; then
      printf '%s\n' $sample_name $fragments_file $cells_file | paste -sd ',' >> $AGGR_CSV
    else
      echo "Warning: Missing files for $sample_name. Skipping." >> $OUTPUT_FILE
    fi
    summary_metrics_file=${sample_path}outs/summary.csv
    # append summary metrics to the metrics CSV
    if [ -f "$summary_metrics_file" ]; then
      if [ ! -f $METRICS_CSV ]; then
        printf '%s\n' $(head -n 1 $summary_metrics_file) >> $METRICS_CSV
      fi
      tail -n +2 $summary_metrics_file >> $METRICS_CSV
    else
      echo "Warning: $summary_metrics_file does not exist, skipping $sample_name" >> $OUTPUT_FILE
    fi
  fi
done

## Run aggregation pipeline
echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE
cellranger-atac aggr --id=${PROJECT_NAME}_aggr \
  --csv $AGGR_CSV --reference $ATAC_REF_PATH --normalize none \
  --localcores=$NCPU --localmem=$MEM
