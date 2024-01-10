#!/bin/bash
#
# run_cellranger_ATAC.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger count pipeline formatted for the ATAC + ASAP side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the ATAC and ASAP 
#       fastqs and has formatted an ASAP feature reference
#       table for the HTO demultiplexing.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/$ATAC_FLOWCELL_ID
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ATAC
OUTPUT_FILE=$OUTPUT_DIR/cellranger_atac_mapping.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
sample_names=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv | uniq)

CR_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
echo "$(date) Running Cell Ranger ATAC version $CR_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE
ATAC_REF_version=$(grep -Po '(?<=refdata-cellranger-arc-)+' $ATAC_REF_PATH)
echo "$(date) Using epigenome reference $ATAC_REF_version located at $ATAC_REF_PATH" >> $OUTPUT_FILE

for sample in $(grep '*${ATAC_NAMING_ID}*' $sample_names); do
  echo "$(date) Running sample ${sample}..." >> $OUTPUT_FILE

  cellranger-atac count --id=$sample --sample=$sample \
    --reference=$ATAC_REF_PATH --fastqs=$FASTQ_PATH \
    --localcores=$NCPU --localmem=$MEM
done