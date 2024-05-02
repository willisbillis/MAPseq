#!/bin/bash
#
# unit_tests.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs unit tests to ensure the MAPseq pipeline is ready to be run.
################################################################################
# FUNCTIONS
function check_dir_var() {
    if ![ -d $1 ]; then
        echo "Directory $1 not found! Please check your project config file."
        return 1
    else;
        return 0
    fi
}

function check_path_var() {
    if ![ -f $1 ]; then
        echo "File $1 not found! Please check your project config file."
        return 1
    else;
        return 0
    fi
}
################################################################################
# check that the input parameters are valid
dir_vars=($DATA_DOWNLOADS_DIR $ATAC_DIR $RNA_DIR $GEX_REF_PATH \
    $VDJ_REF_PATH $ATAC_REF_PATH)
path_vars=($GEX_FEAT_REF_PATH $ASAP_FEAT_REF_PATH)

for dir in "${dir_vars[@]}"; do
    check_dir_var $dir
done

for path in "${path_vars[@]}"; do
    check_path_var $path
done


# check that there aren't existing runs
existing_projects=$(find "$PROJECT_PATH/pipeline/" -mindepth 2 -maxdepth 2 -type d ! -name "tools" ! -name "reports" -printf '%f ')
if [ ${#existing_projects[@]} > 0 ]; then
  echo -e "Found existing mapping projects: $existing_projects \nRun 'clean_ms_tree $(basename $PROJECT_PATH)' to remove before running pipeline again."
fi

echo "Pre-flight checks complete! Running MAPseq pipeline..."