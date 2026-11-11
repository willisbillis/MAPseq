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
    local configured_dirs=$1 dirs_name=$2 fastqs_name=$3 maxdepth=${4:-} dir root fastq resolved
    local -n dirs=$dirs_name fastqs=$fastqs_name
    local IFS=,
    local -a configured_dir_list find_args
    local -A resolved_fastqs=()
    local fastq_list
    read -r -a configured_dir_list <<< "$configured_dirs"
    dirs=()
    fastqs=()
    for dir in "${configured_dir_list[@]}"; do
        dir=${dir#"${dir%%[![:space:]]*}"}
        dir=${dir%"${dir##*[![:space:]]}"}
        [[ -z $dir ]] && continue
        dirs+=("$dir")
        root=$(realpath -e -- "$dir") || {
            echo "[ERROR] Cannot resolve FASTQ source directory '$dir'." >&2
            return 1
        }
        find_args=(-L "$dir")
        [[ -n $maxdepth ]] && find_args+=(-maxdepth "$maxdepth")
        find_args+=(-type f -name '*.fastq.gz' -print0)
        fastq_list=$(mktemp) || return 1
        if ! find "${find_args[@]}" > "$fastq_list"; then
            echo "[ERROR] Failed to fully traverse FASTQ source directory '$dir'." >&2
            rm -f -- "$fastq_list"
            return 1
        fi
        while IFS= read -r -d '' fastq; do
            [[ -L $fastq ]] && continue
            resolved=$(realpath -e -- "$fastq") || {
                echo "[ERROR] Cannot resolve FASTQ '$fastq'." >&2
                rm -f -- "$fastq_list"
                return 1
            }
            case "$resolved" in
                "$root"/*)
                    if [[ -z ${resolved_fastqs[$resolved]+present} ]]; then
                        resolved_fastqs[$resolved]=1
                        fastqs+=("$fastq")
                    fi
                    ;;
            esac
        done < "$fastq_list"
        rm -f -- "$fastq_list"
    done
}

collect_fastqs "$RNA_DIR" rna_dirs rna_fqs || exit 1
collect_fastqs "$ATAC_DIR" atac_dirs atac_fqs || exit 1

copy_manifest_fastqs() {
    local modality=$1 destination=$2 manifest=$3
    local lane sample index fastq filename suffix matched naming_id
    shift 3
    local -a naming_ids=()
    while [[ $1 != -- ]]; do
        naming_ids+=("$1")
        shift
    done
    shift
    local -a fastqs=("$@")
    local missing=0
    local -A basename_sources=() selected_sources=()
    local -a selected_fastqs=()

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
        if (( ! naming_match )); then
            echo "[ERROR] Manifest sample '$sample' does not match any configured $modality naming ID." >&2
            missing=1
            continue
        fi

        matched=0
        for fastq in "${fastqs[@]}"; do
            filename=${fastq##*/}
            [[ $filename == "$sample"* ]] || continue
            suffix=${filename#"$sample"}
            if [[ $suffix =~ ^(_S[0-9]+)?_L[0-9]{3}_[IR][0-9]_[0-9]{3}\.fastq\.gz$ ]]; then
                if [[ -n ${basename_sources[$filename]+present} && ${basename_sources[$filename]} != "$fastq" ]]; then
                    echo "[ERROR] Duplicate $modality FASTQ destination basename '$filename': '${basename_sources[$filename]}' and '$fastq'." >&2
                    return 1
                fi
                basename_sources[$filename]=$fastq
                if [[ -z ${selected_sources[$fastq]+present} ]]; then
                    selected_sources[$fastq]=1
                    selected_fastqs+=("$fastq")
                fi
                matched=1
            fi
        done
        if (( ! matched )); then
            echo "[ERROR] No FASTQs found for $modality sample '$sample'" >&2
            missing=1
        fi
    done < "$manifest"
    (( ! missing )) || return 1
    mkdir -p "$destination" || return 1
    for fastq in "${selected_fastqs[@]}"; do
        cp -- "$fastq" "$destination/" || return 1
    done
}

filter_bcl_fastqs() {
    local modality=$1 destination=$2 manifest=$3
    shift 3
    local stage_dir
    local -a bcl_dirs=() bcl_fastqs=()

    if ! collect_fastqs "$destination" bcl_dirs bcl_fastqs 1; then
        find "$destination" -maxdepth 1 -type f -name '*.fastq.gz' -delete || \
            echo "[ERROR] Failed to remove unfiltered $modality BCL FASTQs." >&2
        return 1
    fi
    stage_dir=$(mktemp -d "$destination/.filtered-fastqs.XXXXXX") || return 1
    if ! copy_manifest_fastqs "$modality" "$stage_dir" "$manifest" "$@" -- "${bcl_fastqs[@]}"; then
        rm -f -- "$stage_dir"/*.fastq.gz
        rmdir -- "$stage_dir"
        find "$destination" -maxdepth 1 -type f -name '*.fastq.gz' -delete || \
            echo "[ERROR] Failed to remove unfiltered $modality BCL FASTQs." >&2
        return 1
    fi
    if ! find "$destination" -maxdepth 1 -type f -name '*.fastq.gz' -delete; then
        echo "[ERROR] Failed to remove unfiltered $modality BCL FASTQs." >&2
        rm -f -- "$stage_dir"/*.fastq.gz
        rmdir -- "$stage_dir"
        return 1
    fi
    if ! mv -- "$stage_dir"/*.fastq.gz "$destination/"; then
        echo "[ERROR] Failed to stage filtered $modality BCL FASTQs." >&2
        rm -f -- "$stage_dir"/*.fastq.gz
        rmdir -- "$stage_dir"
        return 1
    fi
    rmdir -- "$stage_dir"
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
    filter_bcl_fastqs ATAC "$NEW_FQ_PATH" "${PROJECT_NAME}.ATAC.sampleManifest.csv" \
        "$ATAC_NAMING_ID" "$ASAP_NAMING_ID" || exit 1

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
    filter_bcl_fastqs RNA "$NEW_FQ_PATH" "${PROJECT_NAME}.RNA.sampleManifest.csv" \
        "$GEX_NAMING_ID" "$GEX_FEAT_NAMING_ID" "$VDJ_NAMING_ID" || exit 1

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
