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
# make sure Seurat v5 is used
options(Seurat.object.assay.version = "v5")
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
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

#h5_path = paste0(OUTS_DIR, "/count/cellbender_feature_bc_matrix_filtered.h5")
#sc_data <- Read_CellBender_h5_Mat(h5_path)
#sc_data = list("Gene Expression" =
#                 sc_data[!(grepl("^anti-", rownames(sc_data)) |
#                           grepl("^HTC", rownames(sc_data))), ],
#               "Antibody Capture" =
#                sc_data[grepl("^anti-", rownames(sc_data)) |
#                         grepl("^HTC", rownames(sc_data)), ])

sc_data <- Read10X_h5(paste0(OUTS_DIR, "/count/filtered_feature_bc_matrix.h5"))
sc_total <- CreateSeuratObject(
  counts = sc_data$`Gene Expression`,
  assay = "RNA",
  project = PROJECT_NAME
)
aggr_df <- read.csv(paste0(OUTS_DIR, "/aggregation.csv"))
new_sample_names <- factor(aggr_df$sample_id,
                           levels = aggr_df$sample_id,
                           ordered = TRUE)
sc_total$library_id <- new_sample_names[as.integer(gsub(".*-", "",
                                                        colnames(sc_total)))]

adt_data <- sc_data$`Antibody Capture`
hto_data = adt_data[grepl("^HTC", rownames(adt_data)), ]
adt_data = adt_data[!grepl("^HTC", rownames(adt_data)), ]

sc_total[["HTO"]] <- CreateAssay5Object(counts = hto_data)
sc_total[["ADT"]] <- CreateAssay5Object(counts = adt_data)

data_dir <- paste0(OUTPUT_DIR, "/data/")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
hto_reference <- read.csv(HTO_DEMUX_PATH)
hto_reference = hto_reference[hto_reference$library_id %in%
                                intersect(hto_reference$library_id,
                                          aggr_df$sample_id), ]

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
  # Check if there are enough cells to normalize data and demultiplex
  if (ncol(sc_sub) > nrow(sc_sub)) {
    sc_sub = NormalizeData(sc_sub, normalization.method = "CLR",
                           verbose = FALSE)
    hashtag = MULTIseqDemux(sc_sub, autoThresh = TRUE, verbose = TRUE)
    successful_htos = unique(hashtag$MULTI_ID[!(hashtag$MULTI_ID %in%
                                                  c("Doublet", "Negative"))])
    failed_htos = hto_ref_sub$hashtag[!(hto_ref_sub$hashtag %in%
                                          successful_htos)]
    if (length(failed_htos) > 0) {
      print(paste("[WARNING] Demultiplexing failed for the",
                  "following hashtags! Excluding from final object."))
      print("patient_id:")
      print(hto_ref_sub$patient_id[match(failed_htos,
                                         hto_ref_sub$hashtag)])
      print(failed_htos)
    }
    hashtag$patient_id <- hto_ref_sub$patient_id[match(hashtag$MULTI_ID,
                                                       hto_ref_sub$hashtag)]
    hashtag$library_id <- rna_library_id
    hashtag$run_id <- run_id
    if (ncol(hto_ref_sub) > 3) {
      for (metadata_col in colnames(hto_ref_sub)[4:ncol(hto_ref_sub)]) {
        hashtag@meta.data[[metadata_col]] =
          hto_ref_sub[[metadata_col]][match(hashtag$MULTI_ID,
                                            hto_ref_sub$hashtag)]
      }
    }

    rna_sub = subset(sc_total, library_id == {{rna_library_id}})
    rna_sub[["HTO"]] = CreateAssay5Object(counts = hashtag[["HTO"]]$counts,
                                          data = hashtag[["HTO"]]$data)
    rna_sub <- AddMetaData(rna_sub, hashtag@meta.data)

    sub_obj_list[[idx]] <- rna_sub
  } else {
    print("[WARNING] Pool failed. Too few cells to demultiplex.")
  }
}

sub_obj_list = sub_obj_list[lengths(sub_obj_list) != 0]
sc_total <- merge(sub_obj_list[[1]], c(sub_obj_list[2:length(sub_obj_list)]))
sc_total <- JoinLayers(sc_total, assay = "HTO")
sc_total[["HTO"]] = subset(sc_total[["HTO"]],
                           features = 
                             rownames(sc_total)[rownames(sc_total) %in%
                                                (hto_reference$hashtag)])
sc_total <- JoinLayers(sc_total, assay = "RNA")
sc_total <- JoinLayers(sc_total, assay = "ADT")
DefaultAssay(sc_total) <- "RNA"

saveRDS(sc_total, paste0(data_dir, "raw_rna.hto.adt_", PROJECT_NAME, ".RDS"))

# Save the R session environment information
capture.output(sessionInfo(),
               file = paste0(OUTPUT_DIR, "/",
                             PROJECT_NAME,
                             ".Rsession.Info.",
                             gsub("\\D", "", Sys.time()), ".txt"))