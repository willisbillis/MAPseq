#!/bin/bash
#
# run_cellranger_ATAC.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the asap_to_kite_v2.py help script written by
# Caleb Lareau (https://github.com/caleblareau) to prepare ASAP fastq
# files for the kite (https://github.com/pachterlab/kite) pipeline.
# Formatted for the ATAC + ASAP side of the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the ATAC and ASAP fastqs.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ASAP
OUTPUT_FILE=$OUTPUT_DIR/asap_to_kite.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

python_version=$(python --version | grep -Po '(?<=Python )[^;]+')
echo "$(date) Running asap_to_kite_v2.py using python version $python_version and binary $(which python)" >> $OUTPUT_FILE

for sample in $(grep '*${ASAP_NAMING_ID}*' $sample_names); do
    python asap_to_kite_v2.py -f $FASTQ_PATH \
        -s $sample -o $sample -j TotalSeqB \
        -c $NCPU
done
