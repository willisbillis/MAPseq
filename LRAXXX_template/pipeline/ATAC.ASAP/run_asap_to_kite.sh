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
TOOL_PATH=$PROJECT_PATH/pipeline/ATAC.ASAP/tools
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ASAP
OUTPUT_FILE=$OUTPUT_DIR/asap_to_kite.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv)
sample_names=$(printf -- '%s ' "${sample_name_col[@]}" | grep -v Sample | uniq)
asap_samples=($(printf -- '%s ' "${sample_names[@]}" | grep .*${ASAP_NAMING_ID}.*))

python_version=$(python --version | grep -Po '(?<=Python )[^;]+')
echo "$(date) Running asap_to_kite_v2.py using python version $python_version and binary $(which python)" >> $OUTPUT_FILE

for sample in "${asap_samples[@]}"; do
    python $TOOL_PATH/asap_to_kite_v2.py -f $FASTQ_PATH \
        -s $sample -o $OUTPUT_DIR/$sample -j TotalSeqB \
        -c $NCPU
done
