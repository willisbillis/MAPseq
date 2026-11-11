# run_mkfastq.sh - written by MEW (https://github.com/willisbillis) Jan 2023
# This script runs the mkfastq pipeline formatted for the LRA MAPseq project started in 2023.

# NOTICE: At this point, the user has the bcl files from
#       the sequencer. There should be a corresponding 
#       sample sheet with 10x indices to demultiplex the bcls.
################################################################################
# Import all the global variables for this project
source ../project_config.txt

# Set all the local variables for this pipeline
OUTPUT_DIR=$PROJECT_PATH/data
OUTPUT_FILE=$OUTPUT_DIR/cellranger_mkfastq.log
################################################################################
cd "$OUTPUT_DIR" || exit 1
mkdir -p "$OUTPUT_DIR/reports" || exit 1

declare -a rna_dirs=() atac_dirs=() rna_fqs=() atac_fqs=()

collect_fastqs() {
    local configured_dirs=$1 dirs_name=$2 fastqs_name=$3 dir
    local -n dirs=$dirs_name fastqs=$fastqs_name
    local IFS=,
    local -a configured_dir_list
    read -r -a configured_dir_list <<< "$configured_dirs"
    dirs=()
    fastqs=()
    for dir in "${configured_dir_list[@]}"; do
        dir=${dir#"${dir%%[![:space:]]*}"}
        dir=${dir%"${dir##*[![:space:]]}"}
        [[ -z $dir ]] && continue
        dirs+=("$dir")
        while IFS= read -r -d '' fastq; do
            fastqs+=("$fastq")
        done < <(find "$dir" \( -type f -o -type l \) -name '*.fastq.gz' -print0 2>/dev/null)
    done
}

collect_fastqs "$RNA_DIR" rna_dirs rna_fqs
collect_fastqs "$ATAC_DIR" atac_dirs atac_fqs

copy_manifest_fastqs() {
    local modality=$1 destination=$2 manifest=$3
    local lane sample index fastq filename suffix matched missing=0 naming_id
    shift 3
    local -a naming_ids=()
    while [[ $1 != -- ]]; do
        naming_ids+=("$1")
        shift
    done
    shift
    local -a fastqs=("$@")

    mkdir -p "$destination" || return 1
    while IFS=, read -r lane sample index || [[ -n $lane$sample$index ]]; do
        lane=${lane%$'\r'}
        sample=${sample%$'\r'}
        index=${index%$'\r'}
        [[ $lane == Lane && $sample == Sample ]] && continue
        [[ -z $lane$sample$index ]] && continue
        if [[ -z $sample ]]; then
            echo "[ERROR] Empty sample ID in $manifest" >&2
            return 1
        fi

        local naming_match=0
        for naming_id in "${naming_ids[@]}"; do
            if [[ -n $naming_id && $sample == *"$naming_id"* ]]; then
                naming_match=1
                break
            fi
        done
        (( naming_match )) || continue

        matched=0
        for fastq in "${fastqs[@]}"; do
            filename=${fastq##*/}
            [[ $filename == "$sample"* ]] || continue
            suffix=${filename#"$sample"}
            if [[ $suffix =~ ^(_S[0-9]+)?_L[0-9]{3}_[IR][0-9]_[0-9]{3}\.fastq\.gz$ ]]; then
                cp -- "$fastq" "$destination/" || return 1
                matched=1
            fi
        done
        if (( ! matched )); then
            echo "[ERROR] No FASTQs found for $modality sample '$sample'" >&2
            missing=1
        fi
    done < "$manifest"
    (( ! missing ))
}

## ATAC.ASAP fastq generation
if [[ $(wc -l < ${PROJECT_NAME}.ATAC.sampleManifest.csv) -gt 1 ]] && [[ ${#atac_fqs[@]} -eq 0 ]]; then
    if (( ${#atac_dirs[@]} != 1 )); then
        echo "[ERROR] Cell Ranger mkfastq requires one ATAC BCL run directory." >&2
        exit 1
    fi
    # Need to demultiplex fastqs from bcl files, run cellranger mkfastq
    CR_atac_version=$(cellranger-atac --version | grep -Po '(?<=cellranger-atac-)[^;]+')
    echo "[INFO] $(date) Running Cell Ranger ATAC version $CR_atac_version using binary $(which cellranger-atac)" &>> $OUTPUT_FILE

    cellranger-atac mkfastq --id=${PROJECT_NAME}_ATAC --run=${atac_dirs[0]} \
        --csv=${PROJECT_NAME}.ATAC.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM &>> $OUTPUT_FILE

    ATAC_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_ATAC/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    ATAC_FLOWCELL_ID=$(basename $ATAC_FC_PATH)
    # Standardize where fastqs live between received FQs and non-demuxed FQs
    NEW_FQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
    mkdir -p $NEW_FQ_PATH
    mv $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/$ATAC_FLOWCELL_ID/*/*.gz $NEW_FQ_PATH
    mv $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/*.gz $NEW_FQ_PATH
    # rearrange mkfastq outputs
    mkdir -p $NEW_FQ_PATH/mkfastq_outputs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path/* $NEW_FQ_PATH/mkfastq_outputs
    rm -r $PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs/fastq_path

    # MEW (07/19/2024)
    # Commented out - not working anymore, maybe mkfastq doesn't generate htmls anymore?
    #mkfastq_report_dir=$NEW_FQ_PATH/mkfastq_outputs/Reports/html/$ATAC_FLOWCELL_ID/all/all/all
    #cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_ATAC.html
    #cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_ATAC.html

    #cp $OUTPUT_DIR/reports/* $PROJECT_PATH/reports
elif [[ $(wc -l < ${PROJECT_NAME}.ATAC.sampleManifest.csv) -gt 1 ]]; then
    # Copy only FASTQs for samples in the ATAC manifest.
    ATAC_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_ATAC/outs
    copy_manifest_fastqs ATAC "$ATAC_FASTQ_PATH" \
        "${PROJECT_NAME}.ATAC.sampleManifest.csv" \
        "$ATAC_NAMING_ID" "$ASAP_NAMING_ID" -- "${atac_fqs[@]}" || exit 1
fi

## RNA.FB.BCR fastq generation
if [ $(wc -l < ${PROJECT_NAME}.RNA.sampleManifest.csv) -gt 1 ] && [ ${#rna_fqs[@]} -eq 0 ]; then
    if (( ${#rna_dirs[@]} != 1 )); then
        echo "[ERROR] Cell Ranger mkfastq requires one RNA BCL run directory." >&2
        exit 1
    fi
    # Need to demultiplex fastqs from bcl files, run cellranger mkfastq
    CR_version=$(cellranger --version | grep -Po '(?<=cellranger-)[^;]+')
    echo "[INFO] $(date) Running Cell Ranger version $CR_version using binary $(which cellranger)" &>> $OUTPUT_FILE

    cellranger mkfastq --id=${PROJECT_NAME}_RNA --run=${rna_dirs[0]} \
        --csv=${PROJECT_NAME}.RNA.sampleManifest.csv \
        --delete-undetermined \
        --localcores=$NCPU --localmem=$MEM &>> $OUTPUT_FILE

    RNA_FC_PATH=$(ls -d $OUTPUT_DIR/${PROJECT_NAME}_RNA/outs/fastq_path/*/ | grep -v "Reports\|Stats")
    RNA_FLOWCELL_ID=$(basename $RNA_FC_PATH)
    # Standardize where fastqs live between received FQs and non-demuxed FQs
    NEW_FQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
    mkdir -p $NEW_FQ_PATH
    mv $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/$RNA_FLOWCELL_ID/*.gz $NEW_FQ_PATH
    # rearrange mkfastq outputs
    mkdir -p $NEW_FQ_PATH/mkfastq_outputs
    mv $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path/* $NEW_FQ_PATH/mkfastq_outputs
    rm -r $PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs/fastq_path

    # MEW (07/19/2024)
    # Commented out - not working anymore, maybe mkfastq doesn't generate htmls anymore?
    #mkfastq_report_dir=$NEW_FQ_PATH/mkfastq_outputs/Reports/html/$RNA_FLOWCELL_ID/all/all/all
    #cp $mkfastq_report_dir/lane.html $OUTPUT_DIR/reports/lane.stats_${PROJECT_NAME}_RNA.html
    #cp $mkfastq_report_dir/laneBarcode.html $OUTPUT_DIR/reports/laneBarcode.stats_${PROJECT_NAME}_RNA.html

    #cp $OUTPUT_DIR/reports/* $PROJECT_PATH/reports
elif [[ $(wc -l < ${PROJECT_NAME}.RNA.sampleManifest.csv) -gt 1 ]]; then
    # Copy only FASTQs for samples in the RNA manifest.
    RNA_FASTQ_PATH=$PROJECT_PATH/data/${PROJECT_NAME}_RNA/outs
    copy_manifest_fastqs RNA "$RNA_FASTQ_PATH" \
        "${PROJECT_NAME}.RNA.sampleManifest.csv" \
        "$GEX_NAMING_ID" "$GEX_FEAT_NAMING_ID" "$VDJ_NAMING_ID" -- "${rna_fqs[@]}" || exit 1
fi
