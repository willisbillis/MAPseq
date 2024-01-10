library(data.table)
library(Seurat)
library(Matrix)
library(Signac)
library(data.table)
library(ggplot2)
library(stringr)

import_kite_counts <- function(data_path){
  mtx <- fread(paste0(data_path, "/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]], dims=c(dim[[1]],dim[[2]]))
  rownames(matx) <- fread(paste0(data_path, "/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- fread(paste0(data_path, "/featurecounts.genes.txt"), header = FALSE)[[1]]
  return(t(matx))
}

PROJECT_DIR = "/home/boss_lab/Projects/Scharer_sc/LRA.MAPseq/LRA_all"
PROJECT_NAME = "LRA_all"
asap_samples_list = c("LRA001_ASAP","LRA002_ASAP_LTR","LRA002_ASAP_SF","LRA003_ASAP_NRF","LRA003_ASAP_BEL")
aggr_df = read.csv(paste0(PROJECT_DIR,"/pipeline/ATAC.ASAP/ATAC/LRA_all_aggr/outs/aggregation_csv.csv"))
bcs = read.csv(paste0(PROJECT_DIR,"/pipeline/ATAC.ASAP/ATAC/LRA_all_aggr/outs/filtered_peak_bc_matrix/barcodes.tsv"), header=F)
bcs$library_id = aggr_df$library_id[as.numeric(gsub(".*\\-","", bcs$V1))]
hto_demux_df = read.csv("path/to/demux_ref_table.csv")

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
