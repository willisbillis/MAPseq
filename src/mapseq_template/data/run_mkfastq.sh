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
mkdir -p $OUTPUT_DIR/reports

## ATAC.ASAP mkfastq demultiplexing
if [[ $(wc -l ${PROJECT_NAME}.ATAC.sampleManifest.csv) -gt 1 ]]; then
    CR_atac_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
    echo "$(date) Running Cell Ranger ATAC version $CR_atac_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE

    cellranger-atac mkfastq --id=${PROJECT_NAME}_ATAC --run=$ATAC_DIR \
        --csv=${PROJECT_NAME}.ATAC.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM

    # Save Flowcell ID to project config for finding fastqs laters
    #     NB: This checks to see if this variable is already set.
    #     See https://stackoverflow.com/a/13864829
    ATAC_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_ATAC/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    if [[ -z ${ATAC_FLOWCELL_ID+x} ]]; then
        echo "ATAC_FLOWCELL_ID=$(basename $ATAC_FC_PATH)" >> ../project_config.txt
    fi

    mkfastq_report_dir=$OUTPUT_DIR/${PROJECT_NAME}_ATAC/outs/fastq_path/Reports/html/$(basename $ATAC_FC_PATH)/all/all/all
    cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_ATAC.html
    cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_ATAC.html
fi

## RNA.FB.BCR mkfastq demultiplexing
if [[ $(wc -l ${PROJECT_NAME}.RNA.sampleManifest.csv) -gt 1 ]]; then
    CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
    echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

    cellranger mkfastq --id=${PROJECT_NAME}_RNA --run=$RNA_DIR \
        --csv=${PROJECT_NAME}.RNA.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM

    # Save Flowcell ID to project config for finding fastqs later
    #     NB: This checks to see if this variable is already set.
    #     See https://stackoverflow.com/a/13864829
    RNA_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_RNA/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    if [[ -z ${RNA_FLOWCELL_ID+x} ]]; then
        echo "RNA_FLOWCELL_ID=$(basename $RNA_FC_PATH)" >> ../project_config.txt
    fi

    mkfastq_report_dir=$OUTPUT_DIR/${PROJECT_NAME}_ATAC/outs/fastq_path/Reports/html/$(basename $RNA_FC_PATH)/all/all/all
    cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_RNA.html
    cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_RNA.html
fi