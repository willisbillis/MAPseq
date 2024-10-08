# run_seurat.demux_atac.R
# written by MEW (https://github.com/willisbillis) Feb 2024
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

p_load(Seurat, Signac, data.table, Matrix)

set.seed(1234)                # set seed for reproducibility

## Library descriptions ##
# Seurat: functions for single cell data
# Signac: create chromatin assay
# data.table: faster reading/writing of csv tables
# Matrix: processing ASAP counts matrices

################################################################################
# Import all the global variables for this project
PROJECT_PATH = Sys.getenv("PROJECT_PATH")[1]
PROJECT_NAME = Sys.getenv("PROJECT_NAME")[1]
ATAC_NAMING_ID = Sys.getenv("ATAC_NAMING_ID")[1]
ASAP_NAMING_ID = Sys.getenv("ASAP_NAMING_ID")[1]

# Set all the local variables for this pipeline
HTO_DEMUX_PATH = paste0(PROJECT_PATH, "/", PROJECT_NAME,
                        "/pipeline/ATAC.ASAP/hashtag_ref_atac.csv")
OUTS_DIR = paste0(PROJECT_PATH, "/", PROJECT_NAME, "/pipeline/ATAC.ASAP/ATAC/",
                  PROJECT_NAME, "_aggr/outs")
OUTPUT_DIR = paste0(PROJECT_PATH, "/", PROJECT_NAME, "/analysis/ATAC.ASAP")
################################################################################
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

aggr_df = read.csv(paste0(OUTS_DIR, "/aggregation_csv.csv"))
barcodes = read.csv(paste0(OUTS_DIR, "/filtered_peak_bc_matrix/barcodes.tsv"),
                    header = FALSE)
barcodes$library_id = aggr_df$library_id[as.numeric(gsub(".*\\-", "",
                                                         barcodes$V1))]

metadata_df = data.frame(run_id = character(),
                         asap_id = character(),
                         atac_id = character(),
                         HTO_cells = numeric(),
                         ATAC_cells = numeric(),
                         overlap = numeric(),
                         overlap_pct = numeric())

asap_library_paths = Sys.glob(paste0(PROJECT_PATH,
                                     "/*/pipeline/ATAC.ASAP/ASAP/*/"))

for (idx in seq_along(asap_library_paths)) {
  asap_path = asap_library_paths[idx]
  asap_id = basename(asap_path)
  atac_id = gsub(ASAP_NAMING_ID, ATAC_NAMING_ID, asap_id)
  run_id = basename(gsub("\\/pipeline.*", "", asap_path))
  metadata_df[idx, c("run_id", "asap_id", "atac_id")] = c(run_id,
                                                          asap_id,
                                                          atac_id)
}

for (idx in seq_len(nrow(metadata_df))) {
  asap_lib_id = metadata_df[idx, "asap_id"]
  atac_lib_id = metadata_df[idx, "atac_id"]
  run_id = metadata_df[idx, "run_id"]
  features_path = paste0(PROJECT_PATH, "/", run_id,
                         "/pipeline/ATAC.ASAP/ASAP/", asap_lib_id,
                         "/counts_unfiltered/cellranger")
  hto <- Read10X(data.dir = features_path, strip.suffix = TRUE)
  cells = barcodes$V1[barcodes$library_id == atac_lib_id]
  library_suffix = match(atac_lib_id, aggr_df$library_id)
  colnames(hto) = paste0(colnames(hto), "-", library_suffix)
  cmat <- hto[, colnames(hto) %in% cells]

  metadata_df[idx, c("HTO_cells",
                     "ATAC_cells",
                     "overlap",
                     "overlap_pct")] = c(dim(hto)[2],
                                         length(cells),
                                         dim(cmat)[2],
                                         100 * dim(cmat)[2] / length(cells))

  if (!exists("master_ht")) {
    master_ht = cmat
  } else {
    master_ht = cbind(master_ht, cmat)
  }
}

write.csv(metadata_df, paste0("library_stats.", PROJECT_NAME, ".csv"),
          quote = FALSE, row.names = FALSE)

data_dir = paste0(OUTPUT_DIR, "/data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
hto_reference = read.csv(HTO_DEMUX_PATH)

#### LOAD ATAC DATA ####
peak_matrix_file = paste0(OUTS_DIR, "/filtered_peak_bc_matrix.h5")
metadata_file = paste0(OUTS_DIR, "/singlecell.csv")
frag_file = paste0(OUTS_DIR, "/fragments.tsv.gz")
counts <- Read10X_h5(filename = peak_matrix_file)
metadata <- read.csv(file = metadata_file, header = TRUE, row.names = 1)

atac_obj_list = list()

for (idx in seq_len(nrow(metadata_df))) {
  run_id = metadata_df[idx, "run_id"]
  asap_lib_id = metadata_df[idx, "asap_id"]
  atac_lib_id = metadata_df[idx, "atac_id"]

  cells = barcodes$V1[barcodes$library_id == atac_lib_id]
  hto_ref_sub = hto_reference[hto_reference$library_id == atac_lib_id, ]
  # ensure input HTOs match Seurat's replacement of underscores with dashes
  hto_ref_sub$hashtag = gsub("_", "-", hto_ref_sub$hashtag)

  if (length(hto_ref_sub$hashtag) > 1) {
    print(paste("Demultiplexing", asap_lib_id))

    library_ht_hto = master_ht[hto_ref_sub$hashtag,
                               colnames(master_ht) %in% cells]
    hashtag <- CreateSeuratObject(counts = library_ht_hto, assay = "HTO")

    hto_count_sums = rowSums(hashtag@assays$HTO@layers$counts)
    names(hto_count_sums) = rownames(hashtag)

    # Check if there are enough cells to normalize data and demultiplex
    if (ncol(hashtag) > nrow(hashtag)) {
      hashtag <- NormalizeData(hashtag, assay = "HTO",
                               normalization.method = "CLR",
                               verbose = FALSE)
      if (ncol(hashtag) < 25000) {
        hashtag = HTODemux(hashtag, kfunc = "kmeans")
        hashtag$MULTI_ID = hashtag$hash.ID
        hashtag$MULTI_classification = hashtag$HTO_classification
        print(table(hashtag$MULTI_ID))
        print(summary(t(hashtag@assays$HTO@layers$data)))
      } else {
        hashtag = MULTIseqDemux(hashtag, autoThresh = TRUE, verbose = TRUE)
      }
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
    } else {
      print("[WARNING] Pool failed. Too few cells to demultiplex.")
    }
  } else {
    print(paste0("Only one hashtag in ASAP sample", asap_lib_id,
                 ". No demultiplexing performed."))
    library_ht_hto = master_ht[, cells]
    hashtag <- CreateSeuratObject(counts = library_ht_hto, assay = "HTO")
    hashtag <- NormalizeData(hashtag, assay = "HTO",
                             normalization.method = "CLR",
                             verbose = FALSE)
    hashtag$patient_id = hto_ref_sub$patient_id
  }

  hashtag$atac_id = atac_lib_id
  hashtag$asap_id = asap_lib_id
  hashtag$run_id = run_id
  if (ncol(hto_ref_sub) > 3) {
    for (metadata_col in colnames(hto_ref_sub)[4:ncol(hto_ref_sub)]) {
      hashtag@meta.data[[metadata_col]] =
        hto_ref_sub[[metadata_col]][match(hashtag$MULTI_ID,
                                          hto_ref_sub$hashtag)]
    }
  }

  atac_cts_sub = CreateChromatinAssay(counts = counts[, cells],
                                      sep = c(":", "-"),
                                      fragments = frag_file)
  metadata_sub = metadata[cells, ]
  atac_sub = CreateSeuratObject(counts = atac_cts_sub,
                                assay = "ATAC",
                                meta.data = metadata_sub,
                                project = PROJECT_NAME)

  hashtag = hashtag[, intersect(colnames(atac_sub), colnames(hashtag))]
  atac_sub[["HTO"]] = CreateAssay5Object(counts = hashtag[["HTO"]]$counts,
                                         data = hashtag[["HTO"]]$data)
  atac_sub = AddMetaData(atac_sub, hashtag@meta.data)

  atac_obj_list[[idx]] = atac_sub
}

sc_total = merge(atac_obj_list[[1]], c(atac_obj_list[2:idx]))
sc_total = JoinLayers(sc_total, assay = "HTO")
sc_total = sc_total[, !is.na(sc_total$nCount_HTO)]

DefaultAssay(sc_total) = "ATAC"

saveRDS(sc_total, paste0(data_dir, "/raw_atac.hto_", PROJECT_NAME, ".RDS"))

# Save the R session environment information
capture.output(sessionInfo(),
               file = paste0(OUTPUT_DIR, "/",
                             PROJECT_NAME,
                             ".Rsession.Info.",
                             gsub("\\D", "", Sys.time()), ".txt"))