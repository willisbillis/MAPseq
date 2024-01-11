#!/bin/bash
#
# run_mkfastq.sh - written by MEW (https://github.com/willisbillis) Jan 2023
# This script runs the mkfastq pipeline formatted for the LRA MAPseq project started in 2023.

# NOTICE: At this point, the user has the bcl files from
#       the sequencer. There should be a corresponding 
#       sample sheet with 10x indices to demultiplex the bcls.
################################################################################
# Import all the global variables for this project
source ../project_config.txt

# Set all the local variables for this pipeline
OUTPUT_DIR=$PROJECT_PATH/data
OUTPUT_FILE=$OUTPUT_DIR/cellranger_mkfastq.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

## ATAC.ASAP mkfastq demultiplexing
CR_atac_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
echo "$(date) Running Cell Ranger ATAC version $CR_atac_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE

cellranger-atac mkfastq --id=${PROJECT_NAME}_ATAC --run=$ATAC_DIR \
    --csv=${PROJECT_NAME}.ATAC.sampleManifest.csv \
    --delete-undetermined \
    --localcores=$NCPU --localmem=$MEM

# Save Flowcell ID to project config for finding fastqs laters
#     NB: This checks to see if this variable is already set.
#     See https://stackoverflow.com/a/13864829
if [[ -z ${ATAC_FLOWCELL_ID+x} ]]; then
    ATAC_FC_PATH=$(ls -d $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    echo "ATAC_FLOWCELL_ID=$(basename $ATAC_FC_PATH)" >> ../project_config.txt
fi

## RNA.FB.BCR mkfastq demultiplexing
CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

cellranger mkfastq --id=${PROJECT_NAME}_RNA --run=$RNA_DIR \
    --csv=${PROJECT_NAME}.RNA.sampleManifest.csv \
    --delete-undetermined \
    --localcores=$NCPU --localmem=$MEM

# Save Flowcell ID to project config for finding fastqs later
#     NB: This checks to see if this variable is already set.
#     See https://stackoverflow.com/a/13864829
if [[ -z ${RNA_FLOWCELL_ID+x} ]]; then
    RNA_FC_ID=$(ls -d $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    echo "RNA_FLOWCELL_ID=$(basename $RNA_FC_PATH)" >> ../project_config.txt
fi