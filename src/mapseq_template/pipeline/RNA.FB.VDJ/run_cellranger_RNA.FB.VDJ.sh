# run_cellranger_RNA.FB.VDJ.sh - written by MEW (https://github.com/willisbillis) Jan 2024
# This script runs the cellranger multi pipeline formatted for the RNA + HTO + ADT (+ VDJ) side of
# the LRA MAPseq project started in 2023.
#
# NOTICE: At this point, the user has generated the fastqs and
#       has formatted a feature reference table if HTO/ADT data are present.
################################################################################
# Import all the global variables for this project
source ../../project_config.txt

# Set all the local variables for this pipeline
FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
OUTPUT_DIR=$PROJECT_PATH/pipeline/RNA.FB.VDJ
OUTPUT_FILE=$OUTPUT_DIR/cellranger_rna.fb.vdj_mapping.log
################################################################################
mkdir -p "$OUTPUT_DIR/reports" || exit 1
cd "$OUTPUT_DIR" || exit 1

if [[ -z $GEX_NAMING_ID || -z $VDJ_NAMING_ID || -z $GEX_FEAT_NAMING_ID ]]; then
    echo "[ERROR] GEX_NAMING_ID, VDJ_NAMING_ID and GEX_FEAT_NAMING_ID must not be empty" >&2
    exit 1
fi

manifest=$PROJECT_PATH/data/${PROJECT_NAME}.RNA.sampleManifest.csv
if [[ ! -r $manifest ]]; then
    echo "[ERROR] Cannot read sample manifest '$manifest'" >&2
    exit 1
fi
declare -a rna_samples=() vdj_samples=() feature_samples=()
declare -A seen_samples=()
while IFS=, read -r lane sample index || [[ -n $lane$sample$index ]]; do
    lane=${lane%$'\r'}
    sample=${sample%$'\r'}
    index=${index%$'\r'}
    [[ $lane == Lane && $sample == Sample ]] && continue
    [[ -z $lane$sample$index ]] && continue
    if [[ -z $sample ]]; then
        echo "[ERROR] Empty sample ID in $manifest" >&2
        exit 1
    fi
    [[ -n ${seen_samples[$sample]+present} ]] && continue
    seen_samples["$sample"]=1
    matches=0
    [[ $sample == *"$GEX_NAMING_ID"* ]] && ((matches+=1))
    [[ $sample == *"$VDJ_NAMING_ID"* ]] && ((matches+=1))
    [[ $sample == *"$GEX_FEAT_NAMING_ID"* ]] && ((matches+=1))
    if (( matches > 1 )); then
        echo "[ERROR] Sample '$sample' matches multiple RNA, VDJ or feature naming IDs" >&2
        exit 1
    fi
    if [[ $sample == *"$GEX_NAMING_ID"* ]]; then
        rna_samples+=("$sample")
    fi
    if [[ $sample == *"$VDJ_NAMING_ID"* ]]; then
        vdj_samples+=("$sample")
    fi
    if [[ $sample == *"$GEX_FEAT_NAMING_ID"* ]]; then
        feature_samples+=("$sample")
    fi
done < "$manifest"
if (( ${#rna_samples[@]} > 1 || ${#vdj_samples[@]} > 1 || ${#feature_samples[@]} > 1 )); then
    echo "[ERROR] Expected at most one RNA, one VDJ and one feature sample in $manifest" >&2
    exit 1
fi
if (( ${#rna_samples[@]} == 0 )); then
    echo "[WARN] No RNA samples matching '$GEX_NAMING_ID' in $manifest" >&2
fi

CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
echo "$(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" &>> $OUTPUT_FILE
GEX_REF_version=$(echo $GEX_REF_PATH | grep -Po '(?<=refdata-gex-)[^;]+')
echo "$(date) Using transcriptome reference $GEX_REF_version located at $GEX_REF_PATH" &>> $OUTPUT_FILE
if (( ${#feature_samples[@]} > 0 )); then
    if [[ ! -f $GEX_FEAT_REF_PATH ]]; then
        echo "[ERROR] Feature reference '$GEX_FEAT_REF_PATH' not found" >&2
        exit 1
    fi
    echo "$(date) Using HTO/ADT feature reference located at $GEX_FEAT_REF_PATH" &>> "$OUTPUT_FILE"
fi
if [ ${#vdj_samples[@]} != 0 ]; then
    VDJ_REF_version=$(echo $VDJ_REF_PATH | grep -Po '(?<=refdata-cellranger-vdj-)[^;]+')
    echo "$(date) Using vdj reference $VDJ_REF_version located at $VDJ_REF_PATH" &>> $OUTPUT_FILE
fi

for sample in "${rna_samples[@]}"; do
    echo "$(date) Running sample ${sample}..." &>> $OUTPUT_FILE
    SAMPLE_CONFIG_CSV=$OUTPUT_DIR/${sample}_config.csv
    vdj_sample=${vdj_samples[0]-}
    feature_sample=${feature_samples[0]-}
    {
        echo "[gene-expression]"
        printf 'reference,%s\ncreate-bam,true\n' "$GEX_REF_PATH"
        if [[ -n $feature_sample ]]; then
            echo "[feature]"
            printf 'reference,%s\n' "$GEX_FEAT_REF_PATH"
        fi
        if [[ -n $vdj_sample ]]; then
            echo "[vdj]"
            printf 'reference,%s\n' "$VDJ_REF_PATH"
        fi
        echo "[libraries]"
        echo "fastq_id,fastqs,feature_types"
        printf '%s,%s,Gene Expression\n' "$sample" "$FASTQ_PATH"
        if [[ -n $feature_sample ]]; then
            printf '%s,%s,Antibody Capture\n' "$feature_sample" "$FASTQ_PATH"
        fi
        if [[ -n $vdj_sample ]]; then
            printf '%s,%s,VDJ\n' "$vdj_sample" "$FASTQ_PATH"
        fi
    } > "$SAMPLE_CONFIG_CSV" || exit 1

    # Run the Cell Ranger multi command for the sample
    cellranger multi --id $sample \
    --csv $SAMPLE_CONFIG_CSV \
    --localcores $NCPU --localmem $MEM &>> $OUTPUT_FILE

    cp $sample/outs/per_sample_outs/$sample/web_summary.html $OUTPUT_DIR/reports/mapping.report_${sample}.html

    metrics_csv=$(find "$sample/outs" -type f \
        -name metrics_summary.csv -print -quit)
    if [[ -z $metrics_csv ]]; then
        metrics_csv=$(find "$sample/outs" -type f \
            -name summary.csv -print -quit)
    fi
    if [[ -n $metrics_csv ]]; then
        cp -- "$metrics_csv" "$OUTPUT_DIR/reports/mapping.metrics_${sample}.csv"
    else
        echo "$(date) WARNING: No Cell Ranger metrics CSV found for ${sample}" &>> "$OUTPUT_FILE"
    fi
done