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

p_load(Seurat, hdf5r)
p_load_gh("samuel-marsh/scCustomize")

set.seed(1234) # set seed for reproducibility
## Library descriptions ##
# Seurat: functions for single cell data
################################################################################
# Import all the global variables for this project
PROJECT_PATH <- Sys.getenv("PROJECT_PATH")
PROJECT_NAME <- Sys.getenv("PROJECT_NAME")
GEX_NAMING_ID <- Sys.getenv("GEX_NAMING_ID")
GEX_FEAT_NAMING_ID <- Sys.getenv("GEX_FEAT_NAMING_ID")

# Set all the local variables for this pipeline
HTO_DEMUX_PATH <- paste0(PROJECT_PATH, "/", PROJECT_NAME, "/pipeline/RNA.FB.VDJ/hashtag_ref_rna.csv")
OUTS_DIR <- paste0(PROJECT_PATH, "/", PROJECT_NAME, "/pipeline/RNA.FB.VDJ/", PROJECT_NAME, "_aggr/outs")
OUTPUT_DIR <- paste0(PROJECT_PATH, "/", PROJECT_NAME, "/analysis/RNA.FB.VDJ")
################################################################################
dir.create(OUTPUT_DIR, showWarnings = F, recursive = T)

sc.data <- Read_CellBender_h5_Mat(paste0(OUTS_DIR, "/count/cellbender_feature_bc_matrix.h5"))
sc.data = list("Gene Expression" = sc.data[!(grepl("^anti-", rownames(sc.data)) & grepl("^HTC", rownames(sc.data))), ],
               "Antibody Capture" = sc.data[grepl("^anti-", rownames(sc.data)), ],
               "Hashtag" = sc.data[grepl("^HTC", rownames(sc.data)), ])
sc_total <- CreateSeuratObject(
  counts = sc.data$`Gene Expression`,
  assay = "RNA",
  project = PROJECT_NAME
)
aggr_df <- read.csv(paste0(OUTS_DIR, "/aggregation.csv"))
new_sample_names <- factor(aggr_df$sample_id, levels = aggr_df$sample_id, ordered = TRUE)
sc_total$library_id <- new_sample_names[as.integer(gsub(".*-", "", colnames(sc_total)))]

adt.data <- sc.data$`Antibody Capture`
hto.data = sc.data$Hashtag

sc_total[["HTO"]] <- CreateAssay5Object(counts = hto.data)
sc_total[["ADT"]] <- CreateAssay5Object(counts = adt.data)

data_dir <- paste0(OUTPUT_DIR, "/data/")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
hto_reference <- read.csv(HTO_DEMUX_PATH)
hto_reference = hto_reference[hto_reference$library_id %in% aggr_df$sample_id, ]

sub_obj_list <- list()

for (idx in seq_len(nrow(aggr_df))) {
  rna_library_id <- aggr_df[idx, "sample_id"]
  run_id <- basename(gsub("\\/pipeline.*", "", aggr_df[idx, "molecule_h5"]))

  print(paste("Demultiplexing", rna_library_id))

  hto_ref_sub <- hto_reference[hto_reference$library_id == rna_library_id, ]
  # ensure input HTOs match Seurat's replacement of underscores with dashes
  hto_ref_sub$hashtag <- gsub("_", "-", hto_ref_sub$hashtag)
  sc_sub <- sc_total[, colnames(sc_total)[sc_total$library_id ==
                                            rna_library_id]]
  DefaultAssay(sc_sub) <- "HTO"

  hto_counts <- sc_sub@assays$HTO@layers$counts
  rownames(hto_counts) = rownames(sc_sub)
  colnames(hto_counts) = colnames(sc_sub) 
  hto_counts <- hto_counts[hto_ref_sub$hashtag, ]
  # check to see if there are more than 0 non-zero hashtag counts cells
  if (sum(colSums(hto_counts) > 0) > 0) {
    hashtag <- CreateSeuratObject(counts = hto_counts, assay = "HTO")
    hto_count_sums = rowSums(hashtag@assays$HTO@layers$counts)
    names(hto_count_sums) = rownames(hashtag)

    # check for failed hashtags (< 3 HTO count per cell per HTO on average)
    if (sum(hto_count_sums < (3 * ncol(hashtag) / nrow(hto_ref_sub))) > 0) {
      print(paste("[WARNING] Hashtag staining failed for the",
                  "following hashtags! Excluding from final object."))
      failed_htos = names(hto_count_sums[hto_count_sums <
                                           3 * (ncol(hashtag) /
                                                  nrow(hto_ref_sub))])
      print("patient_id:")
      print(hto_ref_sub$patient_id[match(failed_htos,
                                        hto_ref_sub$hashtag)])
      print(failed_htos)
      hto_ref_sub = hto_ref_sub[!(hto_ref_sub$hashtag %in% failed_htos), ]
    }

    if (nrow(hto_ref_sub) > 1) {
      hashtag = subset(hashtag, features = hto_ref_sub$hashtag)
      hashtag <- NormalizeData(hashtag, assay = "HTO",
                               normalization.method = "CLR",
                               verbose = FALSE)
      cells_keep = colnames(hashtag)[colSums(hashtag) > 0 &
                                       hashtag$nCount_HTO >
                                         summary(hashtag$nCount_HTO)[2]]
      cells_keep = cells_keep[!is.na(cells_keep)]
      hashtag = subset(hashtag, cells = cells_keep)

      if (length(cells_keep) > nrow(hashtag)) {
        hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = 0.99,
                            verbose = FALSE)

        hashtag$patient_id <- hto_ref_sub$patient_id[match(hashtag$hash.ID,
                                                           hto_ref_sub$hashtag)]
        hashtag$library_id <- rna_library_id
        hashtag$run_id <- run_id

        sub_obj_list[[idx]] <- hashtag
      }
    }
  }
}

merged_hashtag <- merge(sub_obj_list[[1]], c(sub_obj_list[2:idx]))
merged_hashtag <- JoinLayers(merged_hashtag)

sc_total[["HTO"]] <- CreateAssay5Object(counts = merged_hashtag[["HTO"]]$counts, data = merged_hashtag[["HTO"]]$data)
sc_total <- AddMetaData(sc_total, merged_hashtag@meta.data)
DefaultAssay(sc_total) <- "RNA"

saveRDS(sc_total, paste0(data_dir, "raw_rna.hto.adt_", PROJECT_NAME, ".RDS"))

# Save the R session environment information
capture.output(sessionInfo(),
               file=paste0(OUTPUT_DIR, "/",
                           PROJECT_NAME,
                           ".Rsession.Info.",
                           gsub("\\D", "", Sys.time()), ".txt"))