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

FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
ASAP_FASTQ_PATH=$FASTQ_PATH/ASAP
ATAC_FASTQ_PATH=$FASTQ_PATH/ATAC
TOOL_PATH=$PROJECT_PATH/pipeline/ATAC.ASAP/tools
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ASAP
OUTPUT_FILE=$OUTPUT_DIR/asap_to_kite.log
################################################################################

# Create output and FASTQ subdirectories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $ASAP_FASTQ_PATH
mkdir -p $ATAC_FASTQ_PATH

# Move ASAP and ATAC FASTQs into their respective subfolders if not already separated
for fq in $FASTQ_PATH/*.fastq.gz; do
  fname=$(basename "$fq")
  # Move ASAP FASTQs if the filename contains the ASAP_NAMING_ID substring
  if [[ "$fname" == *"$ASAP_NAMING_ID"* ]]; then
    if [[ ! -f "$ASAP_FASTQ_PATH/$fname" ]]; then
      mv "$fq" "$ASAP_FASTQ_PATH/"
    fi
  # Move ATAC FASTQs if the filename contains the ATAC_NAMING_ID substring
  elif [[ "$fname" == *"$ATAC_NAMING_ID"* ]]; then
    if [[ ! -f "$ATAC_FASTQ_PATH/$fname" ]]; then
      mv "$fq" "$ATAC_FASTQ_PATH/"
    fi
  fi
done

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

# Check for correct ASAP FASTQ files (R1, R2, R3) in the separated ASAP folder
for sample in "${asap_samples[@]}"; do
    r1_file=$(ls $ASAP_FASTQ_PATH/${sample}*R1*.fastq.gz 2>/dev/null | head -n1)
    r2_file=$(ls $ASAP_FASTQ_PATH/${sample}*R2*.fastq.gz 2>/dev/null | head -n1)
    r3_file=$(ls $ASAP_FASTQ_PATH/${sample}*R3*.fastq.gz 2>/dev/null | head -n1)
    if [[ -z "$r3_file" ]]; then
        echo "[ERROR] ASAP sample '$sample' is missing an R3 FASTQ file. If you sequenced on a NextSeq 2000 or similar, you must pull ATAC I2 as ASAP R2 and ASAP R2 as ASAP R3. See the documentation for details." | tee -a $OUTPUT_FILE >&2
        exit 1
    fi
    # Run the asap_to_kite_v2.py script
    python $TOOL_PATH/asap_to_kite_v2.py -f $ASAP_FASTQ_PATH \
        -s $sample -o $OUTPUT_DIR/$sample -j TotalSeqB \
        -c $NCPU &>> $OUTPUT_FILE
done
