#!/bin/bash

function create_ms_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 argument required!"
        return 1
    elif [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    elif [ -d $1 ]; then
        echo "Run directory $PWD/$1 already exists!"
        return 1
    fi

    mkdir -p $1/data
    mkdir -p $1/pipeline/RNA.FB.VDJ
    mkdir -p $1/pipeline/ATAC.ASAP
    ln -s $MAPSEQ_REPO_PATH/src/mapseq_template/run_mapseq.sh $1/run_mapseq.sh
    ln -s $MAPSEQ_REPO_PATH/src/mapseq_template/data/* $1/data
    ln -s $MAPSEQ_REPO_PATH/src/mapseq_template/pipeline/RNA.FB.VDJ/* $1/pipeline/RNA.FB.VDJ
    ln -s $MAPSEQ_REPO_PATH/src/mapseq_template/pipeline/ATAC.ASAP/* $1/pipeline/ATAC.ASAP
    echo "Lane,Sample,Index" > $1/data/$1.ATAC.sampleManifest.csv
    echo "Lane,Sample,Index" > $1/data/$1.RNA.sampleManifest.csv
    cp $MAPSEQ_REPO_PATH/src/references/HTOB_feature_ref.csv $1/pipeline/ATAC.ASAP/HTOB_feature_ref.csv
    cp $MAPSEQ_REPO_PATH/src/references/TSC_feature_ref.csv $1/pipeline/RNA.FB.VDJ/TSC_feature_ref.csv

    cat $MAPSEQ_REPO_PATH/src/project_config_header.txt >> $1/project_config.txt
    echo -e "\n" >> $1/project_config.txt
    echo "PROJECT_NAME=$1" >> $1/project_config.txt
    echo "PROJECT_PATH=$PWD/$1" >> $1/project_config.txt
    cat $MAPSEQ_REPO_PATH/src/project_config_body.txt >> $1/project_config.txt
    echo -e "\n" >> $1/project_config.txt
}

function create_ms_aggr_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 argument required!"
        return 1
    elif [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    elif [ -d $1 ]; then
        echo "Run directory $PWD/$1 already exists!"
        return 1
    fi

    mkdir -p $1/pipeline/RNA.FB.VDJ
    mkdir -p $1/pipeline/ATAC.ASAP
    ln -s $MAPSEQ_REPO_PATH/src/mapseq_template/run_mapseq_aggr.sh $1/run_mapseq_aggr.sh
    printf '%s\n' library_id hashtag patient_id | paste -sd ',' >> $1/pipeline/RNA.FB.VDJ/hashtag_ref_rna.csv
    ln -s $MAPSEQ_REPO_PATH/src/aggr_template/pipeline/RNA.FB.VDJ/* $1/pipeline/RNA.FB.VDJ
    printf '%s\n' library_id hashtag patient_id | paste -sd ',' >> $1/pipeline/ATAC.ASAP/hashtag_ref_atac.csv
    ln -s $MAPSEQ_REPO_PATH/src/aggr_template/pipeline/ATAC.ASAP/* $1/pipeline/ATAC.ASAP

    cat $MAPSEQ_REPO_PATH/src/project_config_header.txt >> $1/project_config.txt
    echo -e "\n" >> $1/project_config.txt
    echo "PROJECT_NAME=$1" >> $1/project_config.txt
    echo "PROJECT_PATH=$PWD" >> $1/project_config.txt
    tail -n 9 $MAPSEQ_REPO_PATH/src/project_config_body.txt >> $1/project_config.txt
}