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
    cp $MAPSEQ_REPO_PATH/src/mapseq_template/run_mapseq.sh $1/run_mapseq.sh
    cp -r $MAPSEQ_REPO_PATH/src/mapseq_template/data/* $1/data
    cp $MAPSEQ_REPO_PATH/src/mapseq_template/pipeline/preflight_checks.sh $1/pipeline
    cp -r $MAPSEQ_REPO_PATH/src/mapseq_template/pipeline/RNA.FB.VDJ/* $1/pipeline/RNA.FB.VDJ
    cp -r $MAPSEQ_REPO_PATH/src/mapseq_template/pipeline/ATAC.ASAP/* $1/pipeline/ATAC.ASAP
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
    cp $MAPSEQ_REPO_PATH/src/aggr_template/run_mapseq_aggr.sh $1/run_mapseq_aggr.sh
    echo "library_id,hashtag,patient_id" > $1/pipeline/RNA.FB.VDJ/hashtag_ref_rna.csv
    cp $MAPSEQ_REPO_PATH/src/aggr_template/pipeline/RNA.FB.VDJ/* $1/pipeline/RNA.FB.VDJ
    echo "library_id,hashtag,patient_id" > $1/pipeline/ATAC.ASAP/hashtag_ref_atac.csv
    cp $MAPSEQ_REPO_PATH/src/aggr_template/pipeline/ATAC.ASAP/* $1/pipeline/ATAC.ASAP

    mkdir -p $1/phase2_scripts
    cp -r $MAPSEQ_REPO_PATH/src/phase2_scripts/* $1/phase2_scripts

    cat $MAPSEQ_REPO_PATH/src/project_config_header.txt >> $1/project_config.txt
    echo -e "\n" >> $1/project_config.txt
    echo "PROJECT_NAME=$1" >> $1/project_config.txt
    echo "PROJECT_PATH=$PWD" >> $1/project_config.txt
    echo "GEX_NAMING_ID=GEX" >> $1/project_config.txt
    echo "GEX_FEAT_NAMING_ID=CSP" >> $1/project_config.txt
    echo "ATAC_NAMING_ID=scATAC" >> $1/project_config.txt
    echo "ASAP_NAMING_ID=ASAP" >> $1/project_config.txt
    tail -n 9 $MAPSEQ_REPO_PATH/src/project_config_body.txt >> $1/project_config.txt
}

function clean_ms_tree() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 argument required!"
        return 1
    elif [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    elif [ ! -d $1 ]; then
        echo "Run directory $PWD/$1 does not exist!"
        return 1
    fi

    directory_to_process=$PWD/$1
    source $directory_to_process/project_config.txt
    preserve_list=(
        "run_mapseq.sh"
        "project_config.txt"
        "data"
        "data/run_genfastq.sh"
        "data/$PROJECT_NAME.RNA.sampleManifest.csv"
        "data/$PROJECT_NAME.ATAC.sampleManifest.csv"
        "pipeline"
        "pipeline/preflight_checks.sh"
        "pipeline/ATAC.ASAP"
        $(find "$directory_to_process/pipeline/ATAC.ASAP/tools/*")
        "pipeline/ATAC.ASAP/run_asap_to_kite.sh"
        "pipeline/ATAC.ASAP/run_cellranger_ATAC.sh"
        "pipeline/ATAC.ASAP/run_kite.sh"
        "pipeline/RNA.FB.VDJ"
        "pipeline/RNA.FB.VDJ/run_cellranger_RNA.FB.VDJ.sh"
    )
    preserve_list=(
        $directory_to_process
        "${preserve_list[@]/#/$directory_to_process\/}"
        $GEX_FEAT_REF_PATH
        $ASAP_FEAT_REF_PATH
    )

    find "$directory_to_process" $(printf "! -wholename %s " "${preserve_list[@]}") -delete
}

function update_ms_tree() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 argument required!"
        return 1
    elif [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    elif [ ! -d $1 ]; then
        echo "Run directory $PWD/$1 does not exist!"
        return 1
    fi

    create_ms_run temp_MS_run
    declare -a update_list=($(find "temp_MS_run" ! -name "project_config.txt" ! -name "*sampleManifest.csv" ! -name "pipeline/ATAC.ASAP/tools/*" -print))

    for full_path in "${update_list[@]}"; do
        relative_path=$(echo $full_path | sed -e 's/.*temp_MS_run\///g')
        if [ -d "$full_path" ]; then
            mkdir -p $1/$relative_path
        elif [ ! -f "$1/$relative_path" ]; then
            echo "Creating new file at $1/$relative_path"
            cp $full_path $1/$relative_path
        fi
    done
    clean_ms_tree $1
    rm -r temp_MS_run
}
