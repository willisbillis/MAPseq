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
source ../project_config.txt
################################################################################
# TODO: add unit tests here

$PROJECT_PATH/data/run_mkfastq.sh

$PROJECT_PATH/pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh &
$PROJECT_PATH/pipeline/ATAC.ASAP/run_asap_to_kite.sh && \
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_kite.sh &
$PROJECT_PATH/pipeline/ATAC.ASAP/run_cellranger_ATAC.sh &
