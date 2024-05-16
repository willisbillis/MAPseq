#!/bin/bash
#
# preflight_checks.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs unit tests to ensure the MAPseq pipeline is ready to be run.
################################################################################
# FUNCTIONS
function check_dir_var() {
    if [ ! -d $1 ]; then
        echo "Directory $1 not found! Please check your project config file."
        return 1
    else
        return 0
    fi
}

function check_path_var() {
    if [ ! -f $1 ]; then
        echo "File $1 not found! Please check your project config file."
        return 1
    else
        return 0
    fi
}
################################################################################
# Get global variables
source ../project_config.txt
# check that the input parameters are valid
declare -a dir_vars=($DATA_DOWNLOADS_DIR $ATAC_DIR $RNA_DIR $GEX_REF_PATH \
    $VDJ_REF_PATH $ATAC_REF_PATH)
declare -a path_vars=($GEX_FEAT_REF_PATH $ASAP_FEAT_REF_PATH)

for dir in "${dir_vars[@]}"; do
    check_dir_var $dir || exit 1
done

for path in "${path_vars[@]}"; do
    check_path_var $path || exit 1
done


# check that there aren't existing runs
declare -a existing_projects=($(find "$PROJECT_PATH/pipeline/" -mindepth 2 -maxdepth 2 -type d ! -name "tools" ! -name "reports" -print))

if [ ${#existing_projects[@]} -gt 0 ]; then
  echo "Found existing mapping projects:"
  printf '%s\n' "${existing_projects[@]}"
  echo "Run 'clean_ms_tree $(basename $PROJECT_PATH)' to remove and try running the pipeline again."
  exit 1
fi

echo "Pre-flight checks complete! Running MAPseq pipeline..."