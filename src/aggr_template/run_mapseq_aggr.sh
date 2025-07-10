#!/bin/bash
#
# run_mapseq_aggr.sh - written by MEW (https://github.com/willisbillis) Feb 2023
# This script runs the MAPseq aggregation pipeline.

# NOTICE: At this point, the user has run through all steps in at least two
#       RNA and ATAC runs and is ready to aggregate the runs for a single 
#       project. Additionally, the user has set and checked all variables
#       in the project_config.txt and populated the sample hashtag
#       sheets with the correct hashtags for demultiplexing the pools.
################################################################################
# Import all the global variables for this project
source ./project_config.txt
# Export global variables for R scripts
export PROJECT_PATH=$PROJECT_PATH
export PROJECT_NAME=$PROJECT_NAME
export GEX_NAMING_ID=$GEX_NAMING_ID
export GEX_FEAT_NAMING_ID=$GEX_FEAT_NAMING_ID
export ATAC_NAMING_ID=$ATAC_NAMING_ID
export ASAP_NAMING_ID=$ASAP_NAMING_ID
################################################################################
# TODO: add unit tests here

# check to see if RNA hashtag reference has been populated to run RNA aggr
rna_sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.RNA.sampleManifest.csv)
rna_sample_names=$(printf -- '%s ' "${rna_sample_name_col[@]}" | grep -v library_id | uniq)
# check to see if ATAC hashtag reference has been populated to run ATAC aggr
atac_sample_name_col=$(cut -d, -f2 $PROJECT_PATH/data/${PROJECT_NAME}.ATAC.sampleManifest.csv)
atac_sample_names=$(printf -- '%s ' "${atac_sample_name_col[@]}" | grep -v library_id | uniq)

if [ ${#rna_sample_names[@]} != 0 ]; then
    # run cellranger aggr on the RNA+FB samples
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ && ./run_cellranger.aggr_rna.sh
    # run souporcell genotype demultiplexing for RNA samples in souporcell conda environment
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ && conda run -n souporcell --live-stream ./run_souporcell.demux_rna.sh
    # load RNA cellranger matrices into Seurat and generate a demultiplexed raw Seurat object
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ && Rscript run_seurat.demux_rna.R
fi

if [ ${#atac_sample_names[@]} != 0 ]; then
    # run cellranger aggr on the ATAC samples
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP && ./run_cellranger.aggr_atac.sh
    # run souporcell genotype demultiplexing for ATAC samples in souporcell conda environment
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP && conda run -n souporcell --live-stream ./run_souporcell.demux_atac.sh
    # load ATAC cellranger/kite matrices into Seurat and generate a demultiplexed raw Seurat object
    cd $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP && Rscript run_seurat.demux_atac.R
fi