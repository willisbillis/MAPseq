#!/bin/bash
#
# run_cellranger.aggr_rna.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger aggr pipeline formatted for the RNA + HTO + ADT + VDJ
# side of the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has run through all steps in at least two
#     RNA + HTO + ADT (+ VDJ) runs and is ready to aggregate the runs for a
#     single project.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
SAMPLES_ARRAY=($(ls -d $PROJECT_PATH/*/pipeline/RNA.FB.VDJ/*/))
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ
OUTPUT_FILE=$OUTPUT_DIR/cellranger_rna_aggr.log
AGGR_CSV=$OUTPUT_DIR/aggr.csv
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

## Create aggregation csv
printf '%s\n' sample_id molecule_h5 | paste -sd ',' >> $AGGR_CSV
for sample_path in "${SAMPLES_ARRAY[@]}"; do
  if [ $(basename $sample_path) != reports ]; then
    if [ -d ${sample_path}outs/per_sample_outs ]; then
      multi_samples_array=($(ls -d ${sample_path}outs/per_sample_outs/*/))
      for multi_path in "${multi_samples_array[@]}"; do
        sample_name=$(basename $multi_path)
        count_molecule_file=${multi_path}count/sample_molecule_info.h5
        printf '%s\n' $sample_name $count_molecule_file | paste -sd ',' >> $AGGR_CSV
      done
    else
      sample_name=$(basename $sample_path)
      count_molecule_file=${sample_path}outs/molecule_info.h5
      printf '%s\n' $sample_name $count_molecule_file | paste -sd ',' >> $AGGR_CSV
    fi
  fi
done

## Run aggregation pipeline
echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE
cellranger aggr --id=${PROJECT_NAME}_aggr \
	--csv $AGGR_CSV  --normalize none \
	--localcores=$NCPU --localmem=$MEM
