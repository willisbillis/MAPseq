# run_souporcell.demux_atac.sh
# written by MEW (https://github.com/willisbillis) Sep 2024
# This script runs the souporcell demultiplexing software to assign metadata
# information to cells sequenced with ATAC+ASAPseq.
#
# NOTICE: At this point, the user has run through all steps in at least two
#     ATAC + ASAP runs.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
SAMPLES_ARRAY=($(ls -d $PROJECT_PATH/*/pipeline/ATAC.ASAP/ATAC/*/))
OUTPUT_DIR=$PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/ATAC_demuxing
OUTPUT_FILE=$OUTPUT_DIR/souporcell_atac.log
SOUPORCELL_PATH=/home/Apps/github/souporcell
FASTA=/home/Apps/genomes/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa
HTO_REF=$PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/hashtag_ref_atac.csv
################################################################################
# Preflight checks
if [ ! -d $SOUPORCELL_PATH ]; then
  echo "ERROR: souporcell repo not found at $SOUPORCELL_PATH"
  echo "Follow instructions at \
    https://github.com/wheaton5/souporcell?tab=readme-ov-file#hard-install to \
    complete installation."
  exit 1
fi

#  code to check whether the fasta file exists
if [ ! -f $FASTA ]; then
  echo "ERROR: FASTA file not found at $FASTA"
  echo "Please provide a valid path to the FASTA file."
  exit 1
fi

#  code to check whether the HTO reference file exists
if [ ! -f $HTO_REF ]; then
  echo "ERROR: HTO_REF file not found at $HTO_REF"
  echo "Please provide a valid path to the HTO_REF file."
  exit 1
fi

#  code to check whether a log file has been generated already
if [ -f $OUTPUT_FILE ]; then
  echo "Removing old log file generated at $OUTPUT_FILE"
  rm $OUTPUT_FILE
fi

# Create output directory
mkdir -p $OUTPUT_DIR
################################################################################
for sample_path in "${SAMPLES_ARRAY[@]}"; do
  if [[ $(basename $sample_path) != reports ]]; then
    sample_name=$(basename $sample_path)

    # Count lines in HTO_REF matching the sample_name
    N=$(awk -F',' -v sn="$sample_name" '$1 == sn {count++} END {print count}' $HTO_REF) 

    # Check if N is 0 or unset after attempting to count
    if [ -z "$N" ] || [ "$N" -eq 0 ]; then
      echo "ERROR: No matching HTOs found for sample '$sample_name' in '$HTO_REF'."
      echo "Please check your HTO_REF file and ensure it contains entries for this sample."
      exit 1 
    fi

    BAM=${sample_path}outs/possorted_bam.bam
    BARCODES=${sample_path}outs/filtered_peak_bc_matrix/barcodes.tsv

    #  code to check whether the bam file exists
    if [ ! -f $BAM ]; then
      echo "ERROR: BAM file not found at $BAM"
      echo "Please provide a valid path to the BAM file."
      exit 1
    fi

    #  code to check whether the barcodes file exists
    if [ ! -f $BARCODES ]; then
      echo "ERROR: Barcodes file not found at $BARCODES"
      echo "Please provide a valid path to the barcodes file."
      exit 1
    fi

    echo "Demultiplexing $N samples in pool $sample_name..."
    $SOUPORCELL_PATH/souporcell_pipeline.py \
      -i $BAM \
      -b $BARCODES \
      -f $FASTA \
      -t $NCPU \
      -o $OUTPUT_DIR/$sample_name \
      -k $N \
      --no_umi True >> $OUTPUT_FILE 2>&1
  fi
done