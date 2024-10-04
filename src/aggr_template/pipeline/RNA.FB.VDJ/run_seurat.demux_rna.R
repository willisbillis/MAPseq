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
  m <- GetAssayData(sc_sub, assay = "HTO", layer = "counts")
  m <- m + 1
  new_assay <- CreateAssayObject(counts = m[hto_ref_sub$hashtag, ])
  sc_sub[["HTO"]] = NULL
  sc_sub[["HTO"]] <- new_assay
  DefaultAssay(sc_sub) <- "HTO"
  # Check if there are enough cells to normalize data and demultiplex
  if (ncol(sc_sub) > nrow(sc_sub)) {
    sc_sub = NormalizeData(sc_sub, normalization.method = "CLR",
                           verbose = FALSE)
    hashtag = HTODemux(sc_sub, verbose = TRUE)
    successful_htos = unique(hashtag$HTO_maxID[!(hashtag$HTO_maxID %in%
                                                   c("Doublet", "Negative"))])
    failed_htos = hto_ref_sub$hashtag[!(hto_ref_sub$hashtag %in%
                                          successful_htos)]
    if (length(failed_htos) > 0) {
      print(paste("[WARNING] Demultiplexing failed for the",
                  "following hashtags!"))
      print("patient_id:")
      print(hto_ref_sub$patient_id[match(failed_htos,
                                         hto_ref_sub$hashtag)])
      print(failed_htos)
    }
    hashtag$patient_id <- hto_ref_sub$patient_id[match(hashtag$hash.ID,
                                                         hto_ref_sub$hashtag)]
    if (ncol(hto_ref_sub) > 3) {
        for (metadata_col in colnames(hto_ref_sub)[4:ncol(hto_ref_sub)]) {
          hashtag@meta.data[[metadata_col]] =
            hto_ref_sub[[metadata_col]][match(hashtag$hash.ID,
                                              hto_ref_sub$hashtag)]
        }
      }

    souporcell_clusters = paste0(PROJECT_PATH, "/", PROJECT_NAME,
                                 "/pipeline/RNA.FB.VDJ/RNA_demuxing/",
                                 rna_library_id, "/clusters.tsv")
    if (file.exists(souporcell_clusters)) {
      # Load data
      clusters_data <- read.table(souporcell_clusters, header = TRUE, sep = "\t")
      hashtag$barcode = gsub("\\-.*", "-1", colnames(hashtag))

      # Create combined data frame
      combined_data <- data.frame(
        barcode = clusters_data$barcode[clusters_data$status == "singlet"],
        cluster = clusters_data$assignment[clusters_data$status == "singlet"]
      )

      # Copy metadata columns from hto_ref_sub to combined_data
      combined_data$HTO_maxID = hashtag$HTO_maxID[match(combined_data$barcode, hashtag$barcode)]
      extra_metadata = c()
      for (metadata_col in colnames(hto_ref_sub)[3:ncol(hto_ref_sub)]) {
        combined_data[[metadata_col]] <- hto_ref_sub[[metadata_col]][match(combined_data$HTO_maxID, hto_ref_sub$hashtag)]
        
        if ("unique_sample_id" %in% colnames(combined_data)) {
          combined_data$unique_sample_id = paste0(combined_data$unique_sample_id, "-",
                                                  combined_data[[metadata_col]])
        } else {
          combined_data$unique_sample_id = combined_data[[metadata_col]]
        }
        extra_metadata = c(extra_metadata, metadata_col)
      }

      # Calculate proportions for each unique_sample_id within each cluster
      proportions_df <- combined_data %>%
        filter(!is.na(patient_id)) %>%
        group_by(cluster) %>%
        mutate(total_cluster_cells = n()) %>%
        group_by(unique_sample_id, cluster) %>%
        summarize(
          cluster_proportion = n() / first(total_cluster_cells),
          total_sample_cells = n(),
          .groups = "drop"
        ) %>%
        group_by(unique_sample_id) %>%
        mutate(sample_proportion = total_sample_cells / sum(total_sample_cells)) %>%
        ungroup()

      # Initialize cluster mapping
      cluster_mapping <- data.frame(cluster = integer(), unique_sample_id = character())

      # Rank Choice Voting (Criteria 1 & 2)
      unassigned_clusters <- unique(proportions_df$cluster)
      while(length(unassigned_clusters) > 0) {
        # Rank by criteria 1 then criteria 2
        ranked_df <- proportions_df %>%
          filter(cluster %in% unassigned_clusters) %>%
          arrange(cluster, desc(cluster_proportion), desc(sample_proportion))

        # Iterate through ranked samples to handle potential pre-assignments
        for (i in 1:nrow(ranked_df)) {
          current_row <- ranked_df[i, ]
          current_cluster <- current_row$cluster
          current_sample <- current_row$unique_sample_id

          # Check if the current sample is already assigned to a cluster
          if (!(current_sample %in% cluster_mapping$unique_sample_id)) {
            # Assign the sample to the cluster
            cluster_mapping <- rbind(cluster_mapping, 
                                      data.frame(cluster = current_cluster, 
                                                unique_sample_id = current_sample))

            # Remove assigned cluster and sample from further consideration
            unassigned_clusters <- setdiff(unassigned_clusters, current_cluster)
            proportions_df <- proportions_df %>%
              filter(!(cluster == current_cluster & unique_sample_id == current_sample))

            # Break the inner loop as we've assigned the cluster
            break
          }

          # Check if all samples for this cluster have been assigned
          if (all(ranked_df$unique_sample_id[ranked_df$cluster == current_cluster] %in% cluster_mapping$unique_sample_id)) {
            # If all samples are assigned, remove the cluster from unassigned_clusters
            unassigned_clusters <- setdiff(unassigned_clusters, current_cluster)
            break # Break the inner loop and move to the next cluster
          }
        }
      }

      # Process of Elimination (Criteria 3)
      unassigned_clusters <- setdiff(unique(combined_data$cluster), cluster_mapping$cluster)
      unassigned_samples <- setdiff(unique(combined_data$unique_sample_id), cluster_mapping$unique_sample_id)

      if (length(unassigned_samples) == 1) {
        remaining_sample <- unassigned_samples[1]

        # Find the best cluster for the remaining sample based on cell count
        best_cluster <- combined_data %>%
          filter(unique_sample_id == remaining_sample, cluster %in% unassigned_clusters) %>%
          group_by(cluster) %>%
          summarize(cell_count = n()) %>%
          filter(cell_count > 10) %>%
          arrange(desc(cell_count)) %>%
          slice_head(n = 1) %>%
          pull(cluster)

        # Assign if a suitable cluster is found
        if (length(best_cluster) > 0) {
          cluster_mapping <- rbind(cluster_mapping, data.frame(cluster = best_cluster, unique_sample_id = remaining_sample))
        }
      } else if (length(unassigned_samples) > 1) {
        print(paste0("[WARNING] Unable to assign all clusters for ",
                     "samples using Process of Elimination."))
        print(paste0("[WARNING] Unable to assign cluster for ",
                     "sample ", paste(unassigned_samples), "."))
      }

      # expand out cluster mapping unique_sample_id into original metadata
      cluster_mapping <- separate(cluster_mapping, col = "unique_sample_id",
                                  into = extra_metadata,
                                  sep = "-",
                                  remove = TRUE)

      # Add cluster and status information to hashtag object
      hashtag$genotype_cluster = combined_data$cluster[match(hashtag$barcode, combined_data$barcode)]
      hashtag$genotype_status = clusters_data$status[match(hashtag$barcode, clusters_data$barcode)]

      # Update hashtag metadata
      for (metadata_col in extra_metadata) {
        for (geno_cl in unique(cluster_mapping$cluster)) {
          mask = (hashtag$genotype_cluster == geno_cl) & (is.na(hashtag@meta.data[[metadata_col]]))
          hashtag@meta.data[[metadata_col]][mask] = cluster_mapping[[metadata_col]][cluster_mapping$cluster == geno_cl]
        }
      }

      # Remove extraneous metadata column from hashtag object
      hashtag$barcode = NULL
    }

    hashtag$library_id <- rna_library_id
    hashtag$run_id <- run_id
    if (ncol(hto_ref_sub) > 3) {
      for (metadata_col in colnames(hto_ref_sub)[4:ncol(hto_ref_sub)]) {
        hashtag@meta.data[[metadata_col]] =
          hto_ref_sub[[metadata_col]][match(hashtag$HTO_maxID,
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