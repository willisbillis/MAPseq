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
PROJECT_PATH = Sys.getenv("PROJECT_PATH")
PROJECT_NAME = Sys.getenv("PROJECT_NAME")
GEX_NAMING_ID = Sys.getenv("GEX_NAMING_ID")
GEX_FEAT_NAMING_ID = Sys.getenv("GEX_FEAT_NAMING_ID")

# Set all the local variables for this pipeline
HTO_DEMUX_CSV = paste0(PROJECT_PATH,"/",PROJECT_NAME,"/pipeline/RNA.FB.VDJ/hashtag_ref_rna.csv")
OUTS_DIR = paste0(PROJECT_PATH,"/",PROJECT_NAME,"/pipeline/RNA.FB.VDJ/",PROJECT_NAME,"_aggr/outs")
OUTPUT_DIR = paste0(PROJECT_PATH,"/",PROJECT_NAME,"/analysis/RNA.FB.VDJ/data")
################################################################################
dir.create(OUTPUT_DIR, showWarnings = F, recursive = T)

sc.data = Read10X(data.dir=paste0(OUTS_DIR, "/count/filtered_feature_bc_matrix/"))
sc_total = CreateSeuratObject(counts=sc.data[["Gene Expression"]],
                              assay="RNA",
                              project=PROJECT_NAME)
aggr_df = read.csv(paste0(OUTS_DIR, "/aggregation.csv"))
new_sample_names = factor(aggr_df$sample_id, levels = aggr_df$sample_id, ordered = TRUE)
sc_total$library_id = new_sample_names[as.integer(sub("*.-","",names(sc_total$nCount_RNA)))]

adt.data = sc.data[["Antibody Capture"]]
sc_total[["HTO"]] = CreateAssayObject(counts = adt.data[grepl("HT", rownames(adt.data)),])
sc_total[["ADT"]] = CreateAssayObject(counts = adt.data[!grepl("HT", rownames(adt.data)),])

data_dir = paste0(OUTPUT_DIR,"/data/")
dir.create(data_dir, recursive = T, showWarnings = F)
hto_reference = read.csv(HTO_DEMUX_CSV)

sub_obj_list = list()

for (idx in seq_len(nrow(aggr_df))) {
  rna_library_id = aggr_df[idx, "sample_id"]
  run_id = basename(gsub("\\/pipeline.*","",aggr_df[idx, "fragments"]))

  hto_reference_sub = hto_reference[hto_reference$library_id == rna_library_id,]
  # ensure input HTOs match Seurat's replacement of underscores with dashes
  htos = gsub("_","-",hto_reference_sub$hashtag)
  sub = subset(sc_total, library_id == rna_library_id)
  DefaultAssay(sub) = "HTO"
  print(rownames(sub))
  print(htos)
  print(rowSums(sub))
  sub = subset(sub, features = htos)
  print(sub)
  sub <- NormalizeData(sub, assay = "HTO", normalization.method = "CLR")
  sub <- HTODemux(sub, assay = "HTO", positive.quantile = 0.99)

  sub$patient_id = hto_reference_sub$patient_id[match(htos, sub$HTO_maxID)]
  sub$library_id = rna_library_id
  sub$run_id = run_id

  sub_obj_list[[idx]] = sub
}

merged = merge(sub_obj_list[[1]], c(sub_obj_list[2:idx]))
merged = JoinLayers(merged)
DefaultAssay(merged) = "RNA"
merged = JoinLayers(merged)

saveRDS(merged, paste0(data_dir,"raw_rna.hto.adt_",PROJECT_NAME,".RDS"))
