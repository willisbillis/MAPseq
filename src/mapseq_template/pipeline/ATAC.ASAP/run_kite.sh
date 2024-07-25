# run_kite.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the kite pipeline formatted for the LRA MAPseq project 
# started in 2023.

# NOTICE: At this point, the user has generated the ATAC and ASAP fastqs and
#	    the asap_to_kite_v2.py generation of kite formatted fastqs should have 
#     been run.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
TOOL_PATH=$PROJECT_PATH/pipeline/ATAC.ASAP/tools
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ASAP
OUTPUT_FILE=$OUTPUT_DIR/kite_asap_mapping.log
################################################################################
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

echo "$(date) Using HTO feature reference located at $ASAP_FEAT_REF_PATH" &>> $OUTPUT_FILE
python_version=$(python --version | grep -Po '(?<=Python )[^;]+')
echo "$(date) Running kite using python version $python_version and binary $(which python)" &>> $OUTPUT_FILE
kallisto_version=$(kallisto version | grep -Po '(?<=version )[^;]+')
echo "$(date) Running kallisto using version $kallisto_version and binary $(which kallisto)" &>> $OUTPUT_FILE
bustools_version=$(bustools version | grep -Po '(?<=version )[^;]+')
echo "$(date) Running bustools using version $bustools_version and binary $(which bustools)" &>> $OUTPUT_FILE

for sample in "${asap_samples[@]}"; do
  WD=$OUTPUT_DIR/$sample
  mkdir -p $WD
  cd $WD
  # cut just the HTO names and sequences from the full table
  awk -v OFS="\t" -F"," '{print $5,$2}' $ASAP_FEAT_REF_PATH > FeatureBarcodes.tsv
  # Remove header from TSV
  sed '1d' < FeatureBarcodes.tsv > FeatureBarcodes_noheader.tsv
  rm FeatureBarcodes.tsv

  kb ref -i mismatch.idx -f1 mismatch.fa -g t2g.txt --kallisto /usr/local/bin/kallisto --workflow kite FeatureBarcodes_noheader.tsv &>> $OUTPUT_FILE
  kb count --verbose --cellranger --kallisto /usr/local/bin/kallisto --workflow kite:10xFB -i mismatch.idx -g t2g.txt -x 10XV3 ../${sample}*.gz &>> $OUTPUT_FILE
done
