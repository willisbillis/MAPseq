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
if [ $(wc -l < ${PROJECT_NAME}.ATAC.sampleManifest.csv) -gt 1 ]; then
    CR_atac_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
    echo "$(date) Running Cell Ranger ATAC version $CR_atac_version using binary $(which cellranger-atac)" >> $OUTPUT_FILE

    cellranger-atac mkfastq --id=${PROJECT_NAME}_ATAC --run=$ATAC_DIR \
        --csv=${PROJECT_NAME}.ATAC.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM


    ATAC_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_ATAC/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    # Standardize where fastqs live between given FQs and non-demuxed FQs
    NEW_FQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/$ATAC_FLOWCELL_ID/* $NEW_FQ_PATH
    # rearrange mkfastq outputs
    mkdir -p $NEW_FQ_PATH/mkfastq_outputs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/* $NEW_FQ_PATH/mkfastq_outputs
    rm -r $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path

    mkfastq_report_dir=$NEW_FQ_PATH/mkfastq_outputs/Reports/html/$(basename $ATAC_FC_PATH)/all/all/all
    cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_ATAC.html
    cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_ATAC.html
fi

## RNA.FB.BCR mkfastq demultiplexing
if [ $(wc -l < ${PROJECT_NAME}.RNA.sampleManifest.csv) -gt 1 ]; then
    CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
    echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" >> $OUTPUT_FILE

    cellranger mkfastq --id=${PROJECT_NAME}_RNA --run=$RNA_DIR \
        --csv=${PROJECT_NAME}.RNA.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM

    RNA_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_RNA/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    # Standardize where fastqs live between given FQs and non-demuxed FQs
    NEW_FQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/$RNA_FLOWCELL_ID/* $NEW_FQ_PATH
    # rearrange mkfastq outputs
    mkdir -p $NEW_FQ_PATH/mkfastq_outputs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/* $NEW_FQ_PATH/mkfastq_outputs
    rm -r $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path

    mkfastq_report_dir=$NEW_FQ_PATH/mkfastq_outputs/Reports/html/$(basename $RNA_FC_PATH)/all/all/all
    cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_RNA.html
    cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_RNA.html
fi