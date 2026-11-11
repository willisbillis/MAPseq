# run_cellranger_ATAC.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger count pipeline formatted for the ATAC + ASAP side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the ATAC and ASAP 
#       fastqs and has formatted an ASAP feature reference
#       table for the HTO demultiplexing.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
OUTPUT_DIR=$PROJECT_PATH/pipeline/ATAC.ASAP/ATAC
OUTPUT_FILE=$OUTPUT_DIR/cellranger_atac_mapping.log
################################################################################
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/reports
rm -f "$OUTPUT_DIR"/reports/mapping.report_*.html \
  "$OUTPUT_DIR"/reports/mapping.metrics_*.csv

sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv)
sample_names=$(printf -- '%s ' "${sample_name_col[@]}" | grep -v Sample | uniq)
atac_samples=($(printf -- '%s ' "${sample_names[@]}" | grep .*${ATAC_NAMING_ID}.*))

CR_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
echo "$(date) Running Cell Ranger ATAC version $CR_version using binary $(which cellranger-atac)" &>> $OUTPUT_FILE
ATAC_REF_version=$(echo $ATAC_REF_PATH | grep -Po '(?<=refdata-cellranger-arc-)[^;]+')
echo "$(date) Using epigenome reference $ATAC_REF_version located at $ATAC_REF_PATH" &>> $OUTPUT_FILE

for sample in "${atac_samples[@]}"; do
  echo "$(date) Running sample ${sample}..." &>> $OUTPUT_FILE
  report_html="$OUTPUT_DIR/reports/mapping.report_${sample}.html"
  report_metrics="$OUTPUT_DIR/reports/mapping.metrics_${sample}.csv"
  rm -f "$report_html" "$report_metrics"

  if cellranger-atac count --id="$sample" \
    --sample=$sample \
    --reference=$ATAC_REF_PATH --fastqs=$FASTQ_PATH \
    --localcores=$NCPU --localmem=$MEM &>> $OUTPUT_FILE; then
    :
  else
    status=$?
    echo "$(date) ERROR: Cell Ranger ATAC failed for ${sample} (status ${status})" &>> "$OUTPUT_FILE"
    exit "$status"
  fi

  cp "$OUTPUT_DIR/$sample/outs/web_summary.html" "$report_html" || exit 1

  metrics_csv=$(find "$OUTPUT_DIR/$sample/outs" -type f \
    -name metrics_summary.csv -print -quit)
  if [[ -z $metrics_csv ]]; then
    metrics_csv=$(find "$OUTPUT_DIR/$sample/outs" -type f \
      -name summary.csv -print -quit)
  fi
  if [[ -n $metrics_csv ]]; then
    cp -- "$metrics_csv" "$report_metrics" || exit 1
  else
    echo "$(date) WARNING: No Cell Ranger ATAC metrics CSV found for ${sample}" &>> "$OUTPUT_FILE"
  fi
done