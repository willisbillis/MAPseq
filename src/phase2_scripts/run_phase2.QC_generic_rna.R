# Generic MAPseq RNA phase-2 QC and clustering workflow.
# Review the configuration and QC plots before accepting filtering thresholds.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(Seurat, ggplot2, clustree, future, parallel)
options(Seurat.object.assay.version = "v5")
options(future.rng.onMisuse = "ignore")
set.seed(1234)

###############################################################################
#### USER CONFIGURATION ####
###############################################################################
PROJECT_ROOT <- "/path/to/project_parent"
PROJECT_NAME <- "aggregation_project"
MAPSEQ_REPO_PATH <- Sys.getenv("MAPSEQ_REPO_PATH", unset = "/path/to/MAPseq")
GENOME <- "hg38"                         # hg38 or mm10
SAMPLE_ID_COLUMNS <- c("library_id", "patient_id")
MAX_CORES <- 8
MAX_MEM_GB <- 32
MAX_PCS <- 50
VARIANCE_TARGET <- 0.80
CLUSTER_RESOLUTIONS <- c(0.25, 0.5, 1)

# These are starting points, not universal cutoffs. Inspect the joint QC plots
# and adjust for the biology, tissue, and sequencing protocol in this dataset.
MAX_PCT_MT <- 10
MIN_FEATURES <- 200
MAX_FEATURES <- 10000

###############################################################################
#### INPUTS AND OUTPUTS ####
###############################################################################
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)
PROJECT_DIR <- file.path(PROJECT_ROOT, PROJECT_NAME,
                         "analysis", "RNA.FB.VDJ")
DATA_DIR <- file.path(PROJECT_DIR, "data")
RAW_SEURAT_PATH <- file.path(DATA_DIR,
                             paste0("raw_rna.hto.adt_", PROJECT_NAME, ".RDS"))
METRICS_PATH <- file.path(PROJECT_ROOT, "reports", "sample_metrics.csv")
HELPER_PATH <- file.path(MAPSEQ_REPO_PATH, "src", "phase2_scripts",
                         "phase2_metrics_helpers.R")
if (!file.exists(HELPER_PATH)) {
  stop("Metrics helper not found at ", HELPER_PATH,
       "; set MAPSEQ_REPO_PATH to the local MAPseq repository")
}
source(HELPER_PATH)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(PROJECT_DIR)

if (MAX_MEM_GB > 0) {
  options(future.globals.maxSize = (MAX_MEM_GB / MAX_CORES) * 1024^3)
}
future::plan("multicore", workers = MAX_CORES)

if (!file.exists(RAW_SEURAT_PATH)) {
  stop("Raw Seurat object not found: ", RAW_SEURAT_PATH)
}
if (!GENOME %in% c("hg38", "mm10")) {
  stop("GENOME must be 'hg38' or 'mm10'")
}

sc_total <- readRDS(RAW_SEURAT_PATH)
if (!"RNA" %in% names(sc_total@assays)) {
  stop("The input Seurat object does not contain an RNA assay")
}
DefaultAssay(sc_total) <- "RNA"
metadata <- sc_total@meta.data
if (!all(SAMPLE_ID_COLUMNS %in% colnames(metadata))) {
  stop("Update SAMPLE_ID_COLUMNS; at least one configured column is missing")
}

mt_pattern <- if (GENOME == "hg38") "^MT-" else "^mt-"
sc_total$percent.mt <- PercentageFeatureSet(sc_total, pattern = mt_pattern)
sc_total$percent.Ig <- PercentageFeatureSet(
  sc_total,
  pattern = if (GENOME == "hg38") "^IGH|^IGK|^IGL" else "^Igh|^Igk|^Igl"
)

classification_column <- intersect(
  c("HTO_classification.global", "HTO_classification"), colnames(metadata)
)[1]
if (is.na(classification_column)) {
  doublet_status <- negative_status <- rep(FALSE, ncol(sc_total))
  message("No HTO classification metadata found; HTO doublet/negative filtering skipped")
} else {
  classification <- as.character(metadata[[classification_column]])
  doublet_status <- classification == "Doublet"
  negative_status <- classification == "Negative"
  doublet_status[is.na(doublet_status)] <- FALSE
  negative_status[is.na(negative_status)] <- FALSE
}
sc_total$doublet_status <- doublet_status
sc_total$negative_status <- negative_status
sc_total$sample_id <- make_sample_ids(sc_total@meta.data, SAMPLE_ID_COLUMNS)

###############################################################################
#### INSPECT QC AND FILTER ####
###############################################################################
qc_plot <- DensityScatter(sc_total, "nCount_RNA", "nFeature_RNA",
                          quantiles = TRUE, log_x = TRUE, log_y = TRUE)
ggsave("scatter_rna_counts_vs_features.png", qc_plot, width = 7, height = 6,
       dpi = 300)
qc_plot <- DensityScatter(sc_total, "nFeature_RNA", "percent.mt",
                          quantiles = TRUE, log_x = TRUE)
qc_plot <- qc_plot +
  geom_hline(yintercept = MAX_PCT_MT, linetype = "dashed") +
  geom_vline(xintercept = c(MIN_FEATURES, MAX_FEATURES), linetype = "dashed")
ggsave("scatter_rna_features_vs_mito.png", qc_plot, width = 7, height = 6,
       dpi = 300)

metadata <- sc_total@meta.data
keep_cells <- metadata$nFeature_RNA > MIN_FEATURES &
  metadata$nFeature_RNA < MAX_FEATURES &
  metadata$percent.mt < MAX_PCT_MT &
  !metadata$doublet_status & !metadata$negative_status
keep_cells[is.na(keep_cells)] <- FALSE
sc <- subset(sc_total, cells = rownames(metadata)[keep_cells])
if (ncol(sc) < 3) {
  stop("Fewer than three cells remain after QC; review the filtering thresholds")
}

sample_ids <- unique(metadata$sample_id)
qc_stats <- do.call(rbind, lapply(sample_ids, function(sample) {
  all_cells <- metadata[metadata$sample_id == sample, , drop = FALSE]
  retained <- all_cells[keep_cells[metadata$sample_id == sample], , drop = FALSE]
  run_id <- if ("run_id" %in% colnames(all_cells)) {
    unique(as.character(all_cells$run_id[!is.na(all_cells$run_id)]))
  } else {
    character()
  }
  data.frame(
    sample_id = sample,
    run_id = if (length(run_id) == 1) run_id else "",
    unfiltered_cells = nrow(all_cells),
    filtered_cells = nrow(retained),
    doublet_cells = sum(all_cells$doublet_status),
    negative_cells = sum(all_cells$negative_status),
    mean_depth_before = mean(all_cells$nCount_RNA, na.rm = TRUE),
    mean_depth_after = if (nrow(retained)) mean(retained$nCount_RNA, na.rm = TRUE) else NA_real_,
    mean_features_after = if (nrow(retained)) mean(retained$nFeature_RNA, na.rm = TRUE) else NA_real_,
    median_percent_mt_after = if (nrow(retained)) median(retained$percent.mt, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
stats_path <- "QC.rna_sampleID_filtering.stats.csv"
write.csv(qc_stats, stats_path, quote = FALSE, row.names = FALSE, na = "")

###############################################################################
#### NORMALIZE, SELECT FEATURES, PCA, AND UMAP ####
###############################################################################
DefaultAssay(sc) <- "RNA"
sc <- NormalizeData(sc, normalization.method = "LogNormalize", verbose = FALSE)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000,
                           verbose = FALSE)
features <- VariableFeatures(sc)
sc <- ScaleData(sc, features = features, verbose = FALSE)
max_pcs_available <- min(MAX_PCS, length(features) - 1, ncol(sc) - 1)
if (max_pcs_available < 2) {
  stop("Not enough cells or variable features to run PCA")
}
sc <- RunPCA(sc, features = features, npcs = max_pcs_available,
             verbose = FALSE)

pca_sd <- sc[["pca"]]@stdev
scaled_data <- GetAssayData(sc, assay = "RNA", layer = "scale.data")
row_means <- Matrix::rowMeans(scaled_data)
total_variance <- sum(Matrix::rowSums(scaled_data^2) -
                        ncol(scaled_data) * row_means^2) /
  (ncol(scaled_data) - 1)
pca_variance <- pca_sd^2 / total_variance
cumulative_variance <- cumsum(pca_variance)
n_dims_keep <- which(cumulative_variance >= VARIANCE_TARGET)[1]
if (is.na(n_dims_keep)) {
  n_dims_keep <- length(pca_sd)
  warning("The variance target was not reached; using all computed PCs")
}
write.csv(data.frame(
  pc = seq_along(pca_sd), variance_explained = pca_variance,
  cumulative_variance = cumulative_variance
), "pca_variance_explained.csv", row.names = FALSE)
message("Using ", n_dims_keep, " PCs (", round(100 * cumulative_variance[n_dims_keep], 1),
        "% of scaled-feature variance)")

sc <- FindNeighbors(sc, dims = seq_len(n_dims_keep), verbose = FALSE)
sc <- RunUMAP(sc, dims = seq_len(n_dims_keep), reduction.name = "umap.rna",
              return.model = TRUE, verbose = FALSE)
sc <- FindClusters(sc, resolution = CLUSTER_RESOLUTIONS,
                   algorithm = 4, verbose = FALSE)
cluster_columns <- grep("^RNA_snn_res\\.", colnames(sc@meta.data), value = TRUE)
if (length(cluster_columns) > 0) {
  p <- clustree(sc, prefix = "RNA_snn_res.")
  ggsave("clustree_rna_resolutions.png", p, width = 8, height = 6, dpi = 300)
  selected_cluster_column <- cluster_columns[which.min(vapply(
    cluster_columns, function(column) {
      abs(as.numeric(sub(".*res\\.", "", column)) - 0.5)
    }, numeric(1)
  ))]
  sc$seurat_clusters <- factor(sc@meta.data[[selected_cluster_column]])
  Idents(sc) <- "seurat_clusters"
  markers <- FindAllMarkers(sc, assay = "RNA", min.pct = 0.1,
                            logfc.threshold = 0.25, verbose = FALSE)
  write.csv(markers, "DEG_RNA_clusters.csv", row.names = FALSE)
}

p <- DimPlot(sc, reduction = "umap.rna", group.by = "seurat_clusters",
             label = TRUE)
ggsave("umap_rna_clusters.pdf", p, width = 8, height = 6)
p <- DimPlot(sc, reduction = "umap.rna", group.by = "sample_id")
ggsave("umap_rna_by_sample.pdf", p, width = 8, height = 6)
p <- FeaturePlot(sc, features = c("nCount_RNA", "percent.mt"),
                 reduction = "umap.rna")
ggsave("umap_rna_qc_metrics.pdf", p, width = 9, height = 5)

###############################################################################
#### BACK UP AND UPDATE SAMPLE METRICS ####
###############################################################################
append_qc_metrics(
  METRICS_PATH, qc_stats, "RNA", PROJECT_NAME,
  file.path(PROJECT_NAME, "analysis", "RNA.FB.VDJ", stats_path),
  units = c(unfiltered_cells = "cells", filtered_cells = "cells",
            doublet_cells = "cells", negative_cells = "cells",
            mean_depth_before = "counts_per_cell",
            mean_depth_after = "counts_per_cell",
            mean_features_after = "features_per_cell",
            median_percent_mt_after = "percent"),
  project_root = PROJECT_ROOT
)
plot_depth_cell_correlation(qc_stats, "RNA", PROJECT_DIR)

saveRDS(sc_total, RAW_SEURAT_PATH)
saveRDS(sc, file.path(DATA_DIR, paste0("qc_rna_", PROJECT_NAME, ".RDS")))
capture.output(sessionInfo(), file = paste0(PROJECT_NAME, ".Rsession.Info.",
                                            format(Sys.time(), "%Y%m%d%H%M%S"), ".txt"))