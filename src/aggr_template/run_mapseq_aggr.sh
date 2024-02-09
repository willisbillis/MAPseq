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
################################################################################
# TODO: add unit tests here

# run cellranger aggr on the RNA+FB samples
cd $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ && $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ/run_cellranger.aggr_rna.sh &
# run cellranger aggr on the ATAC samples
cd $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP && $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/run_cellranger.aggr_atac.sh &
wait

# load RNA cellranger matrices into Seurat and generate a demultiplexed raw Seurat object
cd $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ && Rscript $PROJECT_PATH/$PROJECT_NAME/pipeline/RNA.FB.VDJ/run_seurat.demux_rna.R &
# load ATAC cellranger/kite matrices into Seurat and generate a demultiplexed raw Seurat object
cd $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP && Rscript $PROJECT_PATH/$PROJECT_NAME/pipeline/ATAC.ASAP/run_seurat.demux_atac.R &
