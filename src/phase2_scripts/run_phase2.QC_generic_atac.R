# Generic MAPseq ATAC phase-2 QC and clustering workflow.
# Review the configuration and QC plots before accepting filtering thresholds.
# Signac LSI is retained for MAPseq compatibility; inspect depth on the UMAP,
# since LSI components can remain correlated with fragment depth.

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(Seurat, Signac, ggplot2, clustree, future, parallel)
options(Seurat.object.assay.version = "v5")
options(future.rng.onMisuse = "ignore")
set.seed(1234)

###############################################################################
#### USER CONFIGURATION ####
###############################################################################
PROJECT_ROOT <- "/path/to/project_parent"
PROJECT_NAME <- "aggregation_project"
MAPSEQ_REPO_PATH <- Sys.getenv("MAPSEQ_REPO_PATH", unset = "/path/to/MAPseq")
SAMPLE_ID_COLUMNS <- c("atac_id", "patient_id")
MAX_CORES <- 8
MAX_MEM_GB <- 32
MAX_LSI_DIMS <- 50
VARIANCE_TARGET <- 0.80
CLUSTER_RESOLUTIONS <- c(0.25, 0.5, 1)

# These are Signac-specific starting points, not universal cutoffs. TSS and
# fragment scores depend on protocol and method; inspect joint QC plots first.
MIN_PEAK_FRAGMENTS <- 300
MIN_PCT_READS_IN_PEAKS <- 60
MAX_BLACKLIST_RATIO <- 1
MAX_NUCLEOSOME_SIGNAL <- 1
MIN_TSS_ENRICHMENT <- 4

###############################################################################
#### INPUTS AND OUTPUTS ####
###############################################################################
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)
PROJECT_DIR <- file.path(PROJECT_ROOT, PROJECT_NAME,
                         "analysis", "ATAC.ASAP")
DATA_DIR <- file.path(PROJECT_DIR, "data")
RAW_SEURAT_PATH <- file.path(DATA_DIR,
                             paste0("raw_atac.hto_", PROJECT_NAME, ".RDS"))
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

sc_total <- readRDS(RAW_SEURAT_PATH)
if (!"ATAC" %in% names(sc_total@assays)) {
  stop("The input Seurat object does not contain an ATAC assay")
}
DefaultAssay(sc_total) <- "ATAC"
metadata <- sc_total@meta.data
required_metrics <- c("peak_region_fragments", "passed_filters",
                      "blacklist_region_fragments")
if (!all(required_metrics %in% colnames(metadata))) {
  stop("The ATAC metadata is missing one or more fragment QC columns: ",
       paste(setdiff(required_metrics, colnames(metadata)), collapse = ", "))
}
if (!all(SAMPLE_ID_COLUMNS %in% colnames(metadata))) {
  stop("Update SAMPLE_ID_COLUMNS; at least one configured column is missing")
}

sc_total <- NucleosomeSignal(sc_total, verbose = FALSE)
sc_total <- TSSEnrichment(sc_total, verbose = FALSE)
sc_total$pct_reads_in_peaks <- ifelse(
  sc_total$passed_filters > 0,
  100 * sc_total$peak_region_fragments / sc_total$passed_filters,
  NA_real_
)
sc_total$blacklist_ratio <- ifelse(
  sc_total$peak_region_fragments > 0,
  sc_total$blacklist_region_fragments / sc_total$peak_region_fragments,
  NA_real_
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
if ("patient_id" %in% colnames(sc_total@meta.data)) {
  sc_total$patient_id[doublet_status] <- "Doublet"
  sc_total$patient_id[negative_status] <- "Negative"
}
sc_total$sample_id <- make_sample_ids(sc_total@meta.data, SAMPLE_ID_COLUMNS)

###############################################################################
#### INSPECT QC AND FILTER ####
###############################################################################
qc_plot <- DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                          quantiles = TRUE, log_x = TRUE)
qc_plot <- qc_plot +
  geom_hline(yintercept = MIN_TSS_ENRICHMENT, linetype = "dashed") +
  geom_vline(xintercept = MIN_PEAK_FRAGMENTS, linetype = "dashed")
ggsave("scatter_atac_fragments_vs_tss.png", qc_plot, width = 7, height = 6,
       dpi = 300)
qc_plot <- DensityScatter(sc_total, "peak_region_fragments",
                          "pct_reads_in_peaks", quantiles = TRUE,
                          log_x = TRUE)
ggsave("scatter_atac_fragments_vs_frip.png", qc_plot, width = 7, height = 6,
       dpi = 300)
qc_plot <- DensityScatter(sc_total, "peak_region_fragments",
                          "nucleosome_signal", quantiles = TRUE,
                          log_x = TRUE)
ggsave("scatter_atac_fragments_vs_nucleosome.png", qc_plot,
       width = 7, height = 6, dpi = 300)
qc_plot <- DensityScatter(sc_total, "peak_region_fragments", "blacklist_ratio",
                          quantiles = TRUE, log_x = TRUE)
ggsave("scatter_atac_fragments_vs_blacklist.png", qc_plot,
       width = 7, height = 6, dpi = 300)

metadata <- sc_total@meta.data
keep_cells <- metadata$peak_region_fragments > MIN_PEAK_FRAGMENTS &
  metadata$pct_reads_in_peaks > MIN_PCT_READS_IN_PEAKS &
  metadata$blacklist_ratio < MAX_BLACKLIST_RATIO &
  metadata$nucleosome_signal < MAX_NUCLEOSOME_SIGNAL &
  metadata$TSS.enrichment > MIN_TSS_ENRICHMENT &
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
    mean_depth_before = mean(all_cells$peak_region_fragments, na.rm = TRUE),
    mean_depth_after = if (nrow(retained)) mean(retained$peak_region_fragments, na.rm = TRUE) else NA_real_,
    mean_tss_after = if (nrow(retained)) mean(retained$TSS.enrichment, na.rm = TRUE) else NA_real_,
    mean_pct_reads_in_peaks_after = if (nrow(retained)) mean(retained$pct_reads_in_peaks, na.rm = TRUE) else NA_real_,
    mean_nucleosome_signal_after = if (nrow(retained)) mean(retained$nucleosome_signal, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
stats_path <- "QC.atac_sampleID_filtering.stats.csv"
write.csv(qc_stats, stats_path, quote = FALSE, row.names = FALSE, na = "")

###############################################################################
#### TF-IDF, LSI, UMAP, AND CLUSTERING ####
###############################################################################
DefaultAssay(sc) <- "ATAC"
sc <- RunTFIDF(sc, min.cells = 1)
sc <- FindTopFeatures(sc, min.cutoff = "q0")
max_lsi_available <- min(MAX_LSI_DIMS, nrow(sc[["ATAC"]]) - 1, ncol(sc) - 1)
if (max_lsi_available < 3) {
  stop("Not enough cells or peaks to run LSI")
}
sc <- RunSVD(sc, n = max_lsi_available, reduction.name = "lsi",
             verbose = FALSE)

lsi_sd <- sc[["lsi"]]@stdev
lsi_components <- seq.int(2, length(lsi_sd))
lsi_variance <- lsi_sd[lsi_components]^2 /
  sum(lsi_sd[lsi_components]^2)
cumulative_variance <- cumsum(lsi_variance)
n_dims_keep <- which(cumulative_variance >= VARIANCE_TARGET)[1]
if (is.na(n_dims_keep)) {
  n_dims_keep <- length(lsi_components)
  warning("The variance target was not reached; using all computed non-first LSI dimensions")
}
selected_dims <- lsi_components[seq_len(n_dims_keep)]
write.csv(data.frame(
  lsi_component = lsi_components, variance_explained = lsi_variance,
  cumulative_variance = cumulative_variance
), "lsi_variance_explained.csv", row.names = FALSE)
message("Using LSI components ", min(selected_dims), "-", max(selected_dims),
        " (", round(100 * cumulative_variance[n_dims_keep], 1),
        "% of computed non-first-component variance)")

sc <- FindNeighbors(sc, reduction = "lsi", dims = selected_dims,
                    graph.name = c("ATAC_nn", "ATAC_snn"), verbose = FALSE)
sc <- RunUMAP(sc, reduction = "lsi", dims = selected_dims,
              reduction.name = "umap.atac", return.model = TRUE,
              verbose = FALSE)
sc <- FindClusters(sc, graph.name = "ATAC_snn",
                   resolution = CLUSTER_RESOLUTIONS, algorithm = 4,
                   verbose = FALSE)
cluster_columns <- grep("^ATAC_snn_res\\.", colnames(sc@meta.data), value = TRUE)
if (length(cluster_columns) > 0) {
  p <- clustree(sc, prefix = "ATAC_snn_res.")
  ggsave("clustree_atac_resolutions.png", p, width = 8, height = 6, dpi = 300)
  selected_cluster_column <- cluster_columns[which.min(vapply(
    cluster_columns, function(column) {
      abs(as.numeric(sub(".*res\\.", "", column)) - 0.5)
    }, numeric(1)
  ))]
  sc$seurat_clusters <- factor(sc@meta.data[[selected_cluster_column]])
  Idents(sc) <- "seurat_clusters"
  markers <- FindAllMarkers(sc, assay = "ATAC", min.pct = 0.1,
                            logfc.threshold = 0.25, verbose = FALSE)
  write.csv(markers, "DAR_ATAC_clusters.csv", row.names = FALSE)
}

p <- DimPlot(sc, reduction = "umap.atac", group.by = "seurat_clusters",
             label = TRUE)
ggsave("umap_atac_clusters.pdf", p, width = 8, height = 6)
p <- DimPlot(sc, reduction = "umap.atac", group.by = "sample_id")
ggsave("umap_atac_by_sample.pdf", p, width = 8, height = 6)
p <- FeaturePlot(sc, features = c("peak_region_fragments", "TSS.enrichment"),
                 reduction = "umap.atac")
ggsave("umap_atac_qc_metrics.pdf", p, width = 9, height = 5)

###############################################################################
#### BACK UP AND UPDATE SAMPLE METRICS ####
###############################################################################
project_id <- PROJECT_NAME
existing_metrics <- read.csv(METRICS_PATH, stringsAsFactors = FALSE,
                             check.names = FALSE)
if ("project_id" %in% colnames(existing_metrics)) {
  project_ids <- unique(existing_metrics$project_id[!is.na(existing_metrics$project_id) &
                                                        existing_metrics$project_id != ""])
  if (length(project_ids) == 1) project_id <- project_ids
}
append_qc_metrics(
  METRICS_PATH, qc_stats, "ATAC", project_id,
  file.path(PROJECT_NAME, "analysis", "ATAC.ASAP", stats_path),
  units = c(unfiltered_cells = "cells", filtered_cells = "cells",
            doublet_cells = "cells", negative_cells = "cells",
            mean_depth_before = "fragments_per_cell",
            mean_depth_after = "fragments_per_cell",
            mean_tss_after = "tss_enrichment",
            mean_pct_reads_in_peaks_after = "percent",
            mean_nucleosome_signal_after = "score")
)
plot_depth_cell_correlation(qc_stats, "ATAC", PROJECT_DIR)

saveRDS(sc_total, RAW_SEURAT_PATH)
saveRDS(sc, file.path(DATA_DIR, paste0("qc_atac_", PROJECT_NAME, ".RDS")))
capture.output(sessionInfo(), file = paste0(PROJECT_NAME, ".Rsession.Info.",
                                            format(Sys.time(), "%Y%m%d%H%M%S"), ".txt"))