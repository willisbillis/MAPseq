#!/bin/bash
#
# run_cellranger_RNA.FB_multi.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger multi pipeline formatted for the RNA + HTO + ADT + BCR side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the fastqs and
#       has formatted a feature reference table for the 
#       HTO/ADT demultiplexing.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
PIPELINE_NAME=${PROJECT_NAME}_ATAC
OUTPUT_DIR=$PROJECT_DIR/pipeline/ATAC.ASAP/ATAC
OUTPUT_FILE=$OUTPUT_DIR/cellranger_atac_mapping.log
################################################################################
mkdir -p $OUTPUT_DIR

FASTQ_PATH=$PROJECT_PATH/data/LRA001_ATAC/outs/fastq_path
################################################################################
mkdir -p $OUTPUT_DIR


CR_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')

echo "$(date) Running Cell Ranger ATAC version $CR_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE
cd $OUTPUT_DIR

for sample in $(ls -d $SAMPLES_PATH)
do
  echo "$(date) Running sample $(basename $sample)..." >> $OUTPUT_FILE

  cellranger-atac count --id=$sample --sample=$sample --reference=$REF_PATH --fastqs=$FASTQ_PATH --localcores=$NCPU --localmem=$MEM
done

printf '%s\n' library_id fragments cells | paste -sd ',' >> $AGGR_CSV
for sample_path in $(ls -d $OUTPUT_DIR/*/)
do
  sample_name=$(basename $sample_path)
  fragments_file=$(realpath $sample_path)/outs/fragments.tsv.gz
  cells_file=$(realpath $sample_path)/outs/singlecell.csv
  printf '%s\n' $sample_name $fragments_file $cells_file | paste -sd ',' >> $AGGR_CSV
done

echo "$(date) Running sample aggregation..." >> $OUTPUT_FILE

cellranger-atac aggr --id=${PROJECT_NAME}_aggr --csv $AGGR_CSV --reference $REF_PATH --localcores=$NCPU --localmem=$MEM
