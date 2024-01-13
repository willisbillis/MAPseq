#!/bin/bash

function create_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 argument required!"
        return 1
    fi
    if [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    fi

    ln -s $MAPSEQ_REPO_PATH/LRAXXX_template $1
    echo "Lane,Sample,Index" > $1/data/$1.ATAC.sampleManifest.csv
    echo "Lane,Sample,Index" > $1/data/$1.RNA.sampleManifest.csv
    cp $MAPSEQ_REPO_PATH/src/references/HTOB_feature_ref.csv $1/pipeline/ATAC.ASAP/HTOB_feature_ref.csv
    cp $MAPSEQ_REPO_PATH/src/references/TSC_feature_ref.csv $1/pipeline/RNA.FB.BCR/TSC_feature_ref.csv

    cat $MAPSEQ_REPO_PATH/src/project_config_header.txt >> $1/project_config.txt
    echo "PROJECT_NAME=$1" >> $1/project_config.txt
    echo "PROJECT_PATH=$PWD/$1" >> $1/project_config.txt
    cat $MAPSEQ_REPO_PATH/src/project_config_body.txt >> $1/project_config.txt
}

function create_aggr_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 arguments required!"
        return 1
    fi
    if [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    fi

    ln -s $MAPSEQ_REPO_PATH/LRA_all_template $1
}