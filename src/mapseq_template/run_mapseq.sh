#!/bin/bash
#
# run_mapseq.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs the entire MAPseq pipeline.

# NOTICE: At this point, the user has set and checked all variables
#       in the project_config.txt and populated the sample sheets
#       in the data directory with the correct indexes for demultiplexing
#       the bcl files. Additionally, the feature references in the 
#       RNA and ATAC pipeline directories have been checked to ensure
#       they are accurate for the project.
################################################################################
# Import all the global variables for this project
source ./project_config.txt
################################################################################
# TODO: add unit tests here

rna_fqs=$(ls $RNA_DIR/*fastq*)
atac_fqs=%(ls $ATAC_DIR/*fastq*)

if [ $(wc -c <<< $rna_fqs) > 0 | $(wc -c <<< $atac_fqs) > 0 ]; then
    # demultiplex any fastqs available on RNA.FB.VDJ or ATAC.ASAP side
    cd $PROJECT_PATH/data && $PROJECT_PATH/data/run_mkfastq.sh
else
    # mv the fastqs from the data directory
    if [ $(wc -c <<< $rna_fqs) > 0  ]; then
        RNA_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
        mkdir -p $RNA_FASTQ_PATH
        cp $RNA_DIR/*fastq* $RNA_FASTQ_PATH
    fi
    if [ $(wc -c <<< $atac_fqs) > 0  ]; then
        ATAC_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
        mkdir -p $ATAC_FASTQ_PATH
        cp $ATAC_DIR/*fastq* $ATAC_FASTQ_PATH
    fi
fi

# check for any fastqs from RNA.FB.VDJ
rna_fqs=$(ls $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/*fastq*)
if [ $(wc -c <<< $rna_fqs) > 0 ]; then
    cd $PROJECT_PATH/pipeline/RNA.FB.VDJ && $PROJECT_PATH/pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh &
fi

# check for any fastqs from ATAC.ASAP
atac_fqs=%(ls $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/*fastq*)
if [ $(wc -c <<< $atac_fqs) > 0 ]; then
    cd $PROJECT_PATH/pipeline/ATAC.ASAP
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_asap_to_kite.sh && \
        $PROJECT_PATH/pipeline/ATAC.ASAP/run_kite.sh &
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_cellranger_ATAC.sh &
fi

wait

# compile reports from each step
mkdir -p $PROJECT_PATH/reports
cp $PROJECT_PATH/data/reports/* $PROJECT_PATH/reports
cp $PROJECT_PATH/pipeline/RNA.FB.VDJ/reports/* $PROJECT_PATH/reports
cp $PROJECT_PATH/pipeline/ATAC.ASAP/reports/* $PROJECT_PATH/reports