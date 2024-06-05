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
cd $OUTPUT_DIR
sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv)
sample_names=$(printf -- '%s ' "${sample_name_col[@]}" | grep -v Sample | uniq)
asap_samples=($(printf -- '%s ' "${sample_names[@]}" | grep .*${ASAP_NAMING_ID}.*))

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
  barcodes_csv=FeatureBarcodes.csv
  # cut just the HTO names and sequences from the full table
  cut -d, -f2,5 < $ASAP_FEAT_REF_PATH > $barcodes_csv
  # generate the mismatch FASTA and t2g files (for following commands, see tutorial at https://github.com/pachterlab/kite)
  python $TOOL_PATH/kite/featuremap/featuremap.py $barcodes_csv --header &>> $OUTPUT_FILE
  # build kallisto index with mismatch fasta and a k-mer length equal to the length of the Feature Barcodes (of the HTO)
  kallisto index -i FeaturesMismatch.idx -k 15 FeaturesMismatch.fa &>> $OUTPUT_FILE
  # pseudoalign the reads
  kallisto bus -i FeaturesMismatch.idx -o ./ -x 10xv3 -t $NCPU $OUTPUT_DIR/${sample}*fastq.gz &>> $OUTPUT_FILE
  # run bustools (note we are NOT running the whitelist filtering command from the tutorial,
  # we are trusting the barcodes we have are good 10x barcodes. Removed because this filtered too many cells in the past)
  bustools sort -t $NCPU -o output_sorted.bus output.bus &>> $OUTPUT_FILE
  mkdir -p featurecounts/
  # generate the counts matrix using bustools
  bustools count -o featurecounts/featurecounts --genecounts -g FeaturesMismatch.t2g -e matrix.ec -t transcripts.txt output_sorted.bus &>> $OUTPUT_FILE
done
