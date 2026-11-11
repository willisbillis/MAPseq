# preflight_checks.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs unit tests to ensure the MAPseq pipeline is ready to be run.
################################################################################
# FUNCTIONS
function check_dir_var() {
    if [ ! -d "$1" ]; then
        echo "[ERROR] Directory $1 not found! Please check your project config file."
        return 1
    elif [ ! -r "$1" ] || [ ! -x "$1" ]; then
        echo "[ERROR] Directory $1 is not readable or searchable. Please check permissions." >&2
        return 1
    else
        return 0
    fi
}

function check_dir_list_var() {
    local var_name=$1 configured_dirs=$2 dir found=0
    local IFS=,
    local -a dirs
    read -r -a dirs <<< "$configured_dirs"
    for dir in "${dirs[@]}"; do
        dir=${dir#"${dir%%[![:space:]]*}"}
        dir=${dir%"${dir##*[![:space:]]}"}
        [[ -z $dir ]] && continue
        found=1
        check_dir_var "$dir" || return 1
    done
    if (( ! found )); then
        echo "[ERROR] $var_name must contain at least one directory." >&2
        return 1
    fi
}

function check_path_var() {
    if [ ! -f "$1" ]; then
        echo "[ERROR] File $1 not found! Please check your project config file."
        return 1
    else
        return 0
    fi
}
################################################################################
# Get global variables
source ../project_config.txt
# check that the input parameters are valid
declare -a dir_vars=("$DATA_DOWNLOADS_DIR" "$GEX_REF_PATH" \
    "$VDJ_REF_PATH" "$ATAC_REF_PATH")
declare -a path_vars=("$ASAP_FEAT_REF_PATH")
rna_manifest=$PROJECT_PATH/data/${PROJECT_NAME}.RNA.sampleManifest.csv
if [[ ! -r $rna_manifest || -z $GEX_FEAT_NAMING_ID ]]; then
    echo "[ERROR] RNA manifest must be readable and GEX_FEAT_NAMING_ID must not be empty."
    exit 1
fi
while IFS=, read -r lane sample index || [[ -n $lane$sample$index ]]; do
    [[ $lane == Lane && $sample == Sample ]] && continue
    if [[ $sample == *"$GEX_FEAT_NAMING_ID"* ]]; then
        path_vars+=("$GEX_FEAT_REF_PATH")
        break
    fi
done < "$rna_manifest"

for dir in "${dir_vars[@]}"; do
    check_dir_var "$dir" || exit 1
done
check_dir_list_var ATAC_DIR "$ATAC_DIR" || exit 1
check_dir_list_var RNA_DIR "$RNA_DIR" || exit 1

for path in "${path_vars[@]}"; do
    check_path_var "$path" || exit 1
done

# check to make sure manifests aren't both empty
manifest_path=$PROJECT_PATH/data/
rna_mf=$manifest_path/${PROJECT_NAME}.RNA.sampleManifest.csv
atac_mf=$manifest_path/${PROJECT_NAME}.ATAC.sampleManifest.csv
if [ $(wc -l < $rna_mf) -eq 1 ] && [ $(wc -l < $atac_mf) -eq 1 ]; then
   echo "[ERROR] Both manifests are empty! Make sure to populate these for every run."
   exit 1
fi

# check that there aren't existing runs
declare -a existing_projects=($(find "$PROJECT_PATH/pipeline/" -mindepth 2 -maxdepth 2 -type d ! -name "tools" ! -name "reports" -print))

if [ ${#existing_projects[@]} -gt 0 ]; then
  echo "[ERROR] Found existing mapping projects:"
  printf '%s\n' "${existing_projects[@]}"
  echo "Run 'clean_ms_tree $(basename $PROJECT_PATH)' to remove and try running the pipeline again."
  exit 1
fi

echo "[INFO] Pre-flight checks complete! Running MAPseq pipeline..."