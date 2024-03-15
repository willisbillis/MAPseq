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

rna_fqs=$(ls $RNA_DIR/*fastq.gz 2> /dev/null)
atac_fqs=$(ls $ATAC_DIR/*fastq.gz 2> /dev/null)

if [[ $(wc -c <<< $rna_fqs) == 1 && $(wc -c <<< $atac_fqs) == 1 ]]; then
    # demultiplex any fastqs available on RNA.FB.VDJ or ATAC.ASAP side
    cd $PROJECT_PATH/data && $PROJECT_PATH/data/run_mkfastq.sh
fi

# mv the fastqs from the data directory
RNA_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
if [[ $(wc -c <<< $rna_fqs) -gt 1 && $(wc -c <<< $RNA_FASTQ_PATH/*gz) == 1 ]]; then
    mkdir -p $RNA_FASTQ_PATH
    rna_samples=$(ls $RNA_DIR/*fastq.gz)
    for sample in "${rna_samples[@]}"; do
        ln -s $(realpath $sample) $RNA_FASTQ_PATH
    done
fi
ATAC_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
if [[ $(wc -c <<< $atac_fqs) -gt 1 && $(wc -c <<< $ATAC_FASTQ_PATH/*gz) == 1 ]]; then
    mkdir -p $ATAC_FASTQ_PATH
    atac_samples=($(ls $ATAC_DIR/*fastq.gz))
    for sample in "${atac_samples[@]}"; do
        ln -s $(realpath $sample) $ATAC_FASTQ_PATH
    done
fi

# check for any fastqs from RNA.FB.VDJ
rna_fqs=$(ls $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/*fastq.gz 2> /dev/null)
if [[ $(wc -c <<< $rna_fqs) -gt 1 ]]; then
    cd $PROJECT_PATH/pipeline/RNA.FB.VDJ && $PROJECT_PATH/pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh &
fi

# check for any fastqs from ATAC.ASAP
atac_fqs=$(ls $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/*fastq.gz 2> /dev/nulls)
if [[ $(wc -c <<< $atac_fqs) -gt 1 ]]; then
    cd $PROJECT_PATH/pipeline/ATAC.ASAP
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_asap_to_kite.sh && \
        $PROJECT_PATH/pipeline/ATAC.ASAP/run_kite.sh &
    $PROJECT_PATH/pipeline/ATAC.ASAP/run_cellranger_ATAC.sh &
fi

wait

# compile reports from each step
mkdir -p $PROJECT_PATH/reports
if [[ -d $PROJECT_PATH/data/reports ]]; then
    cp $PROJECT_PATH/data/reports/* $PROJECT_PATH/reports
fi
if [[ -d $PROJECT_PATH/pipeline/RNA.FB.VDJ/reports ]]; then
    cp $PROJECT_PATH/pipeline/RNA.FB.VDJ/reports/* $PROJECT_PATH/reports
fi
if [[ -d $PROJECT_PATH/pipeline/ATAC.ASAP/reports ]]; then
    cp $PROJECT_PATH/pipeline/ATAC.ASAP/reports/* $PROJECT_PATH/reports
fi