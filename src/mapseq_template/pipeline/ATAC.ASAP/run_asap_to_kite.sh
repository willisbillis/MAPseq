#!/bin/bash
# run_cellranger_ATAC.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the asap_to_kite_v2.py help script written by
# Caleb Lareau (https://github.com/caleblareau) to prepare ASAP fastq
# files for the kite (https://github.com/pachterlab/kite) pipeline.
#
# NOTICE: At this point, the user has generated the ATAC and ASAP fastqs.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
TOOL_PATH=$PROJECT_PATH/pipeline/ATAC.ASAP/tools
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ASAP
OUTPUT_FILE=$OUTPUT_DIR/asap_to_kite.log
################################################################################
# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Change to the output directory
cd $OUTPUT_DIR

# Extract sample names from the sample manifest
sample_name_col=()
while IFS=',' read -ra array; do
  sample_name_col+=("${array[1]}")
done < $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv

# Filter for ASAP samples and remove duplicates
asap_samples=()
for sample in "${sample_name_col[@]}"; do
  if [[ $sample =~ .*$ASAP_NAMING_ID.* ]]; then
    asap_samples+=("$sample")
  fi
done

# Get Python version
python_version=$(python --version | grep -Po '(?<=Python )[^;]+')

# Log information about the script execution
echo "[INFO] $(date) Running asap_to_kite_v2.py using python version $python_version and binary $(which python)" &>> $OUTPUT_FILE

# Loop through each ASAP sample
for sample in "${asap_samples[@]}"; do
    # Run the asap_to_kite_v2.py script
    python $TOOL_PATH/asap_to_kite_v2.py -f $FASTQ_PATH \
        -s $sample -o $OUTPUT_DIR/$sample -j TotalSeqB \
        -c $NCPU &>> $OUTPUT_FILE
done
