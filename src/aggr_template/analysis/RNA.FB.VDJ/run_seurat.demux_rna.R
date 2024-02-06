#!/bin/bash
#
# run_seurat.demux_rna.R - written by MEW (https://github.com/willisbillis) Feb 2024
# This script creates a demultiplexed raw data Seurat object for downstream
# QC and analysis.
#
# NOTICE: At this point, the user has aggregated multiple runs using cellranger
#       aggr into a single counts matrix.
################################################################################
# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

p_load(Seurat)

set.seed(1234)                # set seed for reproducibility
## Library descriptions ##
# Seurat: functions for single cell data
################################################################################
# Import all the global variables for this project
system(command = "source", args = c("../../project_config.txt"))
PROJECT_PATH = Sys.getenv("PROJECT_PATH")
PROJECT_NAME = Sys.getenv("PROJECT_NAME")
GEX_NAMING_ID = Sys.getenv("GEX_NAMING_ID")
GEX_FEAT_NAMING_ID = Sys.getenv("GEX_FEAT_NAMING_ID")

# Set all the local variables for this pipeline
HTO_DEMUX_CSV = paste0(PROJECT_PATH, "/analysis/RNA.FB.VDJ/hashtag_demux_ref.csv")
OUTS_DIR = paste0(PROJECT_PATH,,"/",PROJECT_NAME"/pipeline/RNA.FB.VDJ/",PROJECT_NAME,"_aggr/outs")
OUTPUT_DIR = paste0(PROJECT_PATH, "/analysis/RNA.FB.VDJ")
################################################################################

aggr_df = read.csv(paste0(OUTS_DIR, "/aggregation_csv.csv"))
barcodes = read.csv(paste0(OUTS_DIR, "/filtered_peak_bc_matrix/barcodes.tsv"), header=F)
bcs$library_id = aggr_df$library_id[as.numeric(gsub(".*\\-", "", bcs$V1))]

for (idx in seq_along(asap_samples_list)) {
  asap_sample_id = asap_samples_list[idx]
  atac_sample_id = gsub("ASAP","scATAC",asap_sample_id)
  features_path = paste0(PROJECT_DIR,"/pipeline/ATAC.ASAP/ASAP/",asap_sample_id,"/featurecounts")
  hto <- import_kite_counts(features_path)
  cells = bcs$V1[bcs$library_id == atac_sample_id]
  sample_suffix = match(atac_sample_id, aggr_df$library_id)
  colnames(hto) = paste0(colnames(hto), "-", sample_suffix)
  cmat <- hto[,colnames(hto) %in% cells]
  print(paste0(dim(cmat)[2]," overlapping cells. (",dim(hto)[2]," in HTO, ",length(cells)," in scATAC)"))
  
  sample_summary = data.frame(sample=gsub("ASAP_","",asap_sample_id), HTO_cells=dim(hto)[2], scATAC_cells=length(cells), overlap=dim(cmat)[2], overlap_pct=100*dim(cmat)[2]/length(cells))
  if (!exists("master_ht")) {
    master_ht = cmat
  } else {
    master_ht = cbind(master_ht, cmat)
  }
  if (!exists("summary_table")) {
    summary_table = sample_summary
  } else {
    summary_table = rbind(summary_table, sample_summary)
  }
}

write.csv(summary_table, paste0("summary_table.",PROJECT_NAME,".csv"), quote=F, row.names=F)

data_dir = paste0(PROJECT_DIR,"/analysis/ATAC.ASAP/data/")
dir.create(data_dir, recursive = T, showWarnings = F)

hashtag_obj_list = list()

for (idx in seq_along(asap_samples_list)) {
  asap_sample_id = asap_samples_list[idx]
  atac_sample_id = gsub("ASAP","scATAC",asap_sample_id)
  cells = bcs$V1[bcs$library_id == atac_sample_id]
  htos = hto_demux_df$hashtag[hto_demux_df$library_id == atac_sample_id]

  library_ht_adt = master_ht[!grepl("HTO",rownames(master_ht)),colnames(master_ht) %in% cells]
  library_ht_hto = master_ht[htos,colnames(master_ht) %in% cells]
  hashtag <- CreateSeuratObject(counts = library_ht_adt, assay = "ADT")
  hashtag[["HTO"]] = CreateAssayObject(counts = library_ht_hto)
  hashtag <- NormalizeData(hashtag, assay = "HTO", normalization.method = "CLR")
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99)
  hashtag$library_id = atac_sample_id

  saveRDS(hashtag, paste0(data_dir,"hto.adt_",asap_sample_id,".RDS"))

  hashtag_obj_list[[idx]] = hashtag
}

merged_hashtag = merge(hashtag_obj_list[[1]], c(hashtag_obj_list[2:idx]))
merged_hashtag = JoinLayers(merged_hashtag)

saveRDS(merged_hashtag, paste0(data_dir,"hto.adt_",PROJECT_NAME,".RDS"))
