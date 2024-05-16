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
echo "Performing MS pipeline input checks..."
cd $PROJECT_PATH/pipeline && $PROJECT_PATH/pipeline/preflight_checks.sh

# create directory to compile reports from each step
mkdir -p $PROJECT_PATH/reports

# demultiplex any fastqs available on RNA.FB.VDJ or ATAC.ASAP side
cd $PROJECT_PATH/data && $PROJECT_PATH/data/run_genfastq.sh

# check for any fastqs from RNA.FB.VDJ
rna_fqs=$(ls $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/*fastq.gz 2>/dev/null)
if [[ $(wc -c <<< $rna_fqs) -gt 1 ]]; then
    cd $PROJECT_PATH/pipeline/RNA.FB.VDJ && $PROJECT_PATH/pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh &&
        cp $PROJECT_PATH/pipeline/RNA.FB.VDJ/reports/* $PROJECT_PATH/reports
fi

# check for any fastqs from ATAC.ASAP
atac_fqs=$(ls $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/*fastq.gz 2>/dev/null)
if [[ $(wc -c <<< $atac_fqs) -gt 1 ]]; then
    cd $PROJECT_PATH/pipeline/ATAC.ASAP
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_asap_to_kite.sh && \
        $PROJECT_PATH/pipeline/ATAC.ASAP/run_kite.sh
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_cellranger_ATAC.sh &&
        cp $PROJECT_PATH/pipeline/ATAC.ASAP/ATAC/reports/* $PROJECT_PATH/reports
fi