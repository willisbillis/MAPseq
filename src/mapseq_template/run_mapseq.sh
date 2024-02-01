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

# demultiplex any fastqs available on RNA.FB.VDJ or ATAC.ASAP side
$PROJECT_PATH/data/run_mkfastq.sh

# check for any fastqs from ATAC.ASAP
if [ -d $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs ]; then
    $PROJECT_PATH/pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh &
fi

# check for any fastqs from RNA.FB.VDJ
if [ -d $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs ]; then
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_asap_to_kite.sh && \
        $PROJECT_PATH/pipeline/ATAC.ASAP/run_kite.sh &
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_cellranger_ATAC.sh &
fi
