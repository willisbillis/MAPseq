# run_phase2.QC_LRA_atac.R
# created by M Elliott Williams (https://github.com/willisbillis) Apr 2024

# Install required packages using the package manager 'pacman'
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
p_load(Seurat, Signac, GenomeInfoDb, AnnotationHub, biovizBase, ggplot2,
       clustree, dplyr, future, parallel, reticulate, harmony, tidyverse,
       Azimuth, scDblFinder, BiocParallel, aggregation)
p_load_gh("SGDDNB/ShinyCell")
p_load_gh("cellgeni/sceasy")

# Set python path to ensure reticulate packages can be used
python_path = system("which python", intern = TRUE)
use_python(python_path)
# silence random number warning
options(future.rng.onMisuse = "ignore")
set.seed(1234)                # set seed for reproducibility

## Library descriptions ##
# Seurat: functions for single cell data
# Signac: DensityScatter function for QC, ATAC processing
# GenomeInfoDb: database for genomic annotations
# AnnotationHub: database for annotations
# biovizBase: necessary for plots
# ggplot2: functions for plotting
# clustree: plotting clusters vs resolution
# dplyr: pipe command '%>%'
# future: multiprocessing limits
# parallel: multiprocessing limits
# reticulate: set which python to use
# harmony: integration method
# tidyverse: separate function
# Azimuth: cell type annotation
# scDblFinder: doublet detection
# BiocParallel: parallelization for scDblFinder
# aggregation: fisher test for scDblFinder
# ShinyCell: Interact with your data

###############################################################################
#### SET RESOURCE LIMITS ####
###############################################################################
max_cores = 32
max_mem = 128
if (max_cores == -1) {
  max_cores = detectCores()
}
if (max_mem != -1) {
  options("future.globals.maxSize" = (max_mem / max_cores) * 1024^3)
}

plan("multicore", workers = max_cores)
###############################################################################
#### OPTIONS ####
###############################################################################
# REPLACE, must be the same as used in MAPseq pipeline
PROJECT_NAME = "LRA_all"
# REPLACE, path to ATAC.ASAP analysis dir from MAPseq pipeline
PROJECT_DIR = paste0("/home/Projects/Scharer_sc/LRA.MAPseq",
                     "/LRA_all/analysis/ATAC.ASAP")
RAW_SEURAT_PATH = paste0(PROJECT_DIR,
                         "/data/raw_atac.hto_", PROJECT_NAME, ".RDS")
HTO_DEMUX_PATH = paste0(PROJECT_DIR,
                        "/../../pipeline/ATAC.ASAP/hashtag_ref_atac.csv")

GENOME = "GRCh38"                     # REPLACE (GRCh38 or GRCm39)
OUTPUT_FIG_WIDTH =  8                 # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8                 # inches, height of output figures
###############################################################################
#### LOAD DATA ####
###############################################################################
setwd(PROJECT_DIR)
sc_total = readRDS(RAW_SEURAT_PATH)
hto_reference = read.csv(HTO_DEMUX_PATH)
###############################################################################
#### PLOT DEMULTIPLEXING RESULTS ####
###############################################################################
DefaultAssay(sc_total) = "HTO"
ncol = ceiling(nrow(sc_total[["HTO"]]) / 3)

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "hash.ID",
            pt.size = 0)
ggsave(paste0("vln_called_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_maxID",
            pt.size = 0)
ggsave(paste0("vln_max_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_classification.global",
            pt.size = 0)
ggsave(paste0("vln_classification_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol*0.5))

###############################################################################
#### CALCULATE QC METRICS (HTO) ####
###############################################################################
sc_total$combo_id = paste0(sc_total$HTO_maxID, "_", sc_total$HTO_secondID)
for (hto1 in unique(sc_total$HTO_maxID)) {
  for (hto2 in unique(sc_total$HTO_secondID)) {
    if (sum(order(c(hto1, hto2)) == c(1, 2)) != 2) {
      combo_id = paste0(hto2, "_", hto1)
      sc_total$combo_id[sc_total$HTO_secondID == hto2 &
                          sc_total$HTO_maxID == hto1] = combo_id
    }
  }
}
margin_stats = aggregate(sc_total$HTO_margin,
                         by = list(sc_total$combo_id,
                                   sc_total$atac_id),
                         FUN = mean)

margin_stats = separate(margin_stats, Group.1,
                        into = c("HT_1st", "HT_2nd"), sep = "_")
colnames(margin_stats) = c("HT_1st", "HT_2nd",
                           "library_id", "hto_separation")
margin_stats = margin_stats[complete.cases(margin_stats), ]
margin_stats = margin_stats[margin_stats$HT_1st != margin_stats$HT_2nd, ]
# sort and grab top pairs and worst pairs
margin_stats = margin_stats[order(margin_stats$hto_separation), ]
best_htos = c(margin_stats$HT_1st[1],
              margin_stats$HT_2nd[1],
              "Doublet")
margin_stats = margin_stats[order(-margin_stats$hto_separation), ]
worst_htos = c(margin_stats$HT_1st[1],
               margin_stats$HT_2nd[1],
               "Doublet")

write.csv(margin_stats,
          paste0("HTB.combos_", PROJECT_NAME, "_metrics.csv"),
          quote = FALSE, row.names = FALSE)

p = ggplot(margin_stats, aes(hto_separation, library_id,
                             group = c(HT_1st, HT_2nd))) +
  geom_boxplot() + labs(title = "HTB Demultiplexing Margins") +
  xlab("Average Margin Between Cells") +
  ylab("HTO Combination") +
  xlim(c(0, max(margin_stats$hto_separation))) +
  theme_linedraw()
ggsave(paste0("boxplot_HTB_hto.separation", PROJECT_NAME, ".png"), p,
       width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

Idents(sc_total) = "hash.ID"

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in%
                                                          best_htos],
                   feature1 = best_htos[1], feature2 = best_htos[2])
ggsave(paste0("scatter_best.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

p = FeatureScatter(sc_total, cells = colnames(sc_total)[sc_total$hash.ID %in%
                                                          worst_htos],
                   feature1 = worst_htos[1], feature2 = worst_htos[2])
ggsave(paste0("scatter_worst.hto.separation_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT, width = OUTPUT_FIG_WIDTH)

sc_total = subset(sc_total, HTO_classification.global != "Doublet")
###############################################################################
#### ATTACH LATEST GENE ANNOTATIONS TO ATAC DATA ####
###############################################################################
DefaultAssay(sc_total) = "ATAC"
# Add the gene annotations to the peaks in the assay
ah <- AnnotationHub()
if (GENOME %in% c("GRCm39")) {
  ahDb <- query(ah, pattern = c("Mus Musculus", "EnsDb"))
} else if (GENOME %in% c("GRCh38")) {
  ahDb <- query(ah, pattern = c("Homo Sapiens", "EnsDb"))
}  else {
  print(paste("[ERROR] Requested genome,",
              GENOME,
              ",not available or not implemented. Exiting."))
  quit()
}

ahDb = ahDb[ahDb$genome == GENOME] # Filter to specific genome
ahDb = ahDb[ahDb$rdatadateadded == ahDb@date] # Find latest version
GENOME_DATA = ahDb[[names(ahDb@.db_uid)]]

# extract gene annotations from EnsDb
annotations = GetGRangesFromEnsDb(ensdb = GENOME_DATA)

# change to UCSC style
annotations <- renameSeqlevels(annotations,
                               mapSeqlevels(seqlevels(annotations), "UCSC"))

# add the gene information to the object
Annotation(sc_total) = annotations
###############################################################################
#### CALCULATE QC METRICS (ATAC) ####
###############################################################################
# compute nucleosome signal score per cell
sc_total = NucleosomeSignal(object=sc_total, verbose = FALSE)

# compute TSS enrichment score per cell
sc_total = TSSEnrichment(object=sc_total, verbose = FALSE)

# add blacklist ratio and fraction of reads in peaks
sc_total$pct_reads_in_peaks = sc_total$peak_region_fragments /
  sc_total$passed_filters * 100
sc_total$blacklist_ratio = sc_total$blacklist_region_fragments /
  sc_total$peak_region_fragments

# Doublet Detection
sc_v3 = sc_total
sc_v3[["HTO"]] = NULL
sc_v3[["ATAC"]] = as(sc_v3[["ATAC"]], Class = "ChromatinAssay")
sce <- scDblFinder(as.SingleCellExperiment(sc_v3), artificialDoublets = 1,
                   aggregateFeatures = TRUE, samples = "library_id",
                   nfeatures = 25, processing = "normFeatures",
                   BPPARAM = MulticoreParam(max_cores))

to_exclude <- GRanges(c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"),
                      IRanges(1L, width = 10^8))
res <- amulet(Fragments(sc_total)@path, regionsToExclude = to_exclude)
res$scDblFinder.p <- 1 - colData(sce)[row.names(res), "scDblFinder.score"]
res$combined <- apply(res[, c("scDblFinder.p", "p.value")], 1,
  FUN = function(x) {
    x[x < 0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
    suppressWarnings(aggregation::fisher(x))
  }
)
sc_total$scDblFinder.score <- res$combined

p = DensityScatter(sc_total, "peak_region_fragments", "scDblFinder.score",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.dbl_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.TSSe_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "pct_reads_in_peaks",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.prip_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "nucleosome_signal",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.nsig_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

p = DensityScatter(sc_total, "peak_region_fragments", "blacklist_ratio",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
ggsave("scatter_peakfrags.v.blacklist_alldata.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### ATAC QC CUTOFFS ####
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MIN_PEAK_FRAGMENTS = 750    # REPLACE, minimum peak fragments per cell
MIN_PCT_RiP = 65            # REPLACE, minimum percent reads in peaks per cell
MAX_BLACKLIST_RATIO = 1.0   # REPLACE, maximum blacklist ratio per cell
MAX_NUCLEOSOME_SIG = 1      # REPLACE, maximum nucleosome signal per cell
MIN_TSS = 4                 # REPLACE, minimum TSS enrichment score per cell
DBL_LIMIT = 0.25            # REPLACE, minimum Doublet score accepted

p = DensityScatter(sc_total, "peak_region_fragments", "TSS.enrichment",
                   quantiles = TRUE, log_x = TRUE, log_y = FALSE)
p = p +
  geom_hline(yintercept=MIN_TSS, linetype = "dashed") +
  geom_vline(xintercept=MIN_PEAK_FRAGMENTS,
             linetype = "dashed")
ggsave("scatter_peakfrags.v.TSSe_filtered.png",
       p, width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### QUANTIFY QC FILTERING ####
###############################################################################
# adjust metadata to accomodate Seurat's AggregateExpression
sc_total$atac_id = gsub("_", "-", sc_total$atac_id)
sc_total$patient_id = gsub("_", "-", sc_total$patient_id)
hto_reference$library_id = gsub("_", "-", hto_reference$library_id)
hto_reference$patient_id = gsub("_", "-", hto_reference$patient_id)
# pair hto reference with seurat object
hto_reference$match_id = paste(hto_reference$library_id,
                               hto_reference$patient_id,
                               hto_reference$hashtag,
                               sep = "-")
sc_total$match_id = paste(sc_total$atac_id,
                          sc_total$patient_id,
                          sc_total$hash.ID,
                          sep = "-")
# add metadata from hto reference to seurat object
for (col_id in names(hto_reference)[4:ncol(hto_reference)]) {
  sc_id = sc_total$match_id
  hto_id = hto_reference$match_id
  sc_total@meta.data[[col_id]] = hto_reference[[col_id]][match(sc_id, hto_id)]
}
# create new column for unique sample ID - adjust as needed for each dataset
hto_reference$sample_id = paste(hto_reference$library_id,
                                hto_reference$patient_id,
                                sep = "-")
sc_total$sample_id = paste(sc_total$atac_id,
                           sc_total$patient_id,
                           sep = "-")
neg_cells_mask = sc_total$HTO_classification.global == "Negative"
sc_total$sample_id[neg_cells_mask] = "Negative"
stats = data.frame(match_id = hto_reference$match_id)

if (ncol(hto_reference) > 3) {
  stats = merge(stats, hto_reference[, 3:ncol(hto_reference)], by="match_id")
} else {
  stats$sample_id = hto_reference$sample_id[match(stats$match_id,
                                                  hto_reference$match_id)]
}
stats$match_id = NULL
sc_total$match_id = NULL
neg_df = data.frame(sample_id = "Negative")
neg_df[names(stats)[names(stats) != "sample_id"]] = NA
stats = rbind(stats, neg_df)
sample_id_counts = as.data.frame(table(sc_total$sample_id))
stats$Unfiltered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                     sample_id_counts$Var1)]
stats[is.na(stats)] = 0
cell_ct = AggregateExpression(sc_total,
                              group.by = "sample_id")$ATAC %>%
  colSums
stats$Unfiltered_Avg_Accessibility = round(cell_ct[match(stats$sample_id,
                                                         names(cell_ct))] /
                                             stats$Unfiltered_Cells)
cell_ct = AggregateExpression(sc_total,
                              group.by = "sample_id")$HTO %>%
  colSums
stats$Unfiltered_Avg_HTO = round(cell_ct[match(stats$sample_id,
                                               names(cell_ct))] /
                                   stats$Unfiltered_Cells)

sc <- subset(sc_total,
             peak_region_fragments > MIN_PEAK_FRAGMENTS &
               pct_reads_in_peaks > MIN_PCT_RiP &
               blacklist_ratio < MAX_BLACKLIST_RATIO &
               nucleosome_signal < MAX_NUCLEOSOME_SIG &
               TSS.enrichment > MIN_TSS &
               scDblFinder.score > DBL_LIMIT)

sample_id_counts = as.data.frame(table(sc$sample_id))
stats$Filtered_Cells = sample_id_counts$Freq[match(stats$sample_id,
                                                   sample_id_counts$Var1)]
stats[is.na(stats)] = 0

cell_ct = AggregateExpression(sc,
                              group.by = "sample_id")$ATAC %>%
  colSums
stats$Filtered_Avg_Accessibility = round(cell_ct[match(stats$sample_id,
                                                       names(cell_ct))] /
                                           stats$Filtered_Cells)
cell_ct = AggregateExpression(sc,
                              group.by = "sample_id")$HTO %>%
  colSums
stats$Filtered_Avg_HTO = round(cell_ct[match(stats$sample_id,
                                             names(cell_ct))] /
                                 stats$Filtered_Cells)
stats[is.na(stats)] = 0

write.csv(stats, "QC.atac_sampleID_filtering.stats.csv",
          quote = FALSE, row.names = FALSE)
###############################################################################
# SAVE RAW SEURAT OBJECT
###############################################################################
saveRDS(sc_total,
        paste0(PROJECT_DIR, "/data/raw_atac.hto_", PROJECT_NAME, ".RDS"))
###############################################################################
#### BATCH CORRECTION (OPTIONAL) ####
###############################################################################
sc_na = sc[, colnames(sc)[is.na(sc$endotype)]]
sc <- sc[, colnames(sc)[!is.na(sc$endotype)]]
sc_na$label_exists = FALSE
sc$label_exists = TRUE

if (FALSE) {
  DefaultAssay(sc) = "ATAC"
  batch_column = "endotype"
  n_dims = 30

  sc_na = sc[, colnames(sc)[is.na(sc[[batch_column]])]]
  sc <- sc[, colnames(sc)[!is.na(sc[[batch_column]])]]
  sc_na$label_exists = FALSE
  sc$label_exists = TRUE

  sc <- RunTFIDF(sc, min.cells = 1)
  sc <- FindTopFeatures(sc, min.cutoff = "q0")
  sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
               verbose = FALSE)

  sc <- FindNeighbors(sc, dims = 2:n_dims, reduction = "atac.lsi",
                      verbose = FALSE)
  sc <- RunUMAP(sc, dims = 2:n_dims, reduction = "atac.lsi",
                reduction.name = "umap.atac.unintegrated",
                return.model = TRUE, verbose = FALSE)
  # visualize by batch annotations
  p = DimPlot(sc, reduction = "umap.atac.unintegrated", group.by = batch_column)
  ggsave(paste0("umap_atac.unintegrated_", batch_column, ".pdf"), p,
         width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
  ggsave(paste0("umap_atac.unintegrated_", batch_column, ".png"), p,
         width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

  sc <- RunHarmony(
    object = sc,
    group.by.vars = batch_column,
    reduction.use = "atac.lsi",
    assay.use = "ATAC",
    project.dim = FALSE
  )
  sc <- FindNeighbors(sc, reduction = "harmony",
                      dims = 2:n_dims, verbose = FALSE)
  sc <- RunUMAP(sc, reduction = "harmony",
                dims = 2:n_dims, reduction.name = "umap.atac",
                return.model = TRUE, verbose = FALSE)

  p <- DimPlot(sc, reduction = "umap.atac", group.by = batch_column)
  ggsave(paste0("umap_atac.integrated_", batch_column, "_known.labels.pdf"), p,
         width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
  ggsave(paste0("umap_atac.integrated_", batch_column, "_known.labels.png"), p,
         width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
}
###############################################################################
#### NON-BATCH CORRECTED DIMENSIONAL REDUCTION ####
###############################################################################
DefaultAssay(sc) = "ATAC"
batch_column = "endotype"
n_dims = 30

sc <- RunTFIDF(sc, min.cells = 1)
sc <- FindTopFeatures(sc, min.cutoff = "q0")
sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
             verbose = FALSE)

sc <- FindNeighbors(sc, dims = 2:n_dims, reduction = "atac.lsi",
                    verbose = FALSE)
sc <- RunUMAP(sc, dims = 2:n_dims, reduction = "atac.lsi",
              reduction.name = "umap.atac.unintegrated",
              return.model = TRUE, verbose = FALSE)
# visualize by batch annotations
p = DimPlot(sc, reduction = "umap.atac.unintegrated", group.by = batch_column)
ggsave(paste0("umap_atac.unintegrated_", batch_column, ".pdf"), p,
       width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
ggsave(paste0("umap_atac.unintegrated_", batch_column, ".png"), p,
       width = OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### LABEL TRANSFER ON NEGATIVE CELLS (OPTIONAL) ####
###############################################################################
if (FALSE) {
  if (!exists("sc_na")) {
    batch_column = "endotype"
    n_dims = 30
    sc_na = sc[, colnames(sc)[is.na(sc[[batch_column]])]]
    sc <- sc[, colnames(sc)[!is.na(sc[[batch_column]])]]
    sc_na$label_exists = FALSE
    sc$label_exists = TRUE
  }

  # find transfer anchors
  DefaultAssay(sc) = "ATAC"
  sc_na = RunTFIDF(sc_na, min.cells = 1)
  sc_na = FindTopFeatures(sc_na, min.cutoff = "q0")
  sc_na <- RunSVD(sc_na, n = n_dims, reduction.name = "atac.lsi",
                  verbose = FALSE)

  transfer_anchors <- FindTransferAnchors(
    reference = sc,
    query = sc_na,
    reference.reduction = "atac.lsi",
    reduction = "lsiproject",
    dims = 2:n_dims
  )

  # map query onto the reference dataset
  sc_na <- MapQuery(
    anchorset = transfer_anchors,
    reference = sc,
    query = sc_na,
    refdata = sc[[batch_column]],
    reference.reduction = "atac.lsi",
    new.reduction.name = "ref.lsi",
    reduction.model = "umap.atac.unintegrated"
  )
  sc_na[[batch_column]] = sc_na$predicted.id
  sc <- merge(sc, sc_na)
  sc = JoinLayers(sc, assay = "HTO")

  # Run dimensionality reduction and clustering on the new
  #      fully annotated dataset
  sc <- RunTFIDF(sc, min.cells = 1)
  sc <- FindTopFeatures(sc, min.cutoff = "q0")
  sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
               verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 2:n_dims, reduction = "atac.lsi",
                      verbose = FALSE)
  sc <- FindClusters(sc, resolution = 2, algorithm = 4,
                     cluster.name = "unintegrated_atac.clusters",
                     verbose = FALSE)
  sc <- RunUMAP(sc, dims = 2:n_dims, reduction = "atac.lsi",
                reduction.name = "umap.atac.unintegrated",
                verbose = FALSE)

  # visualize by batch annotations
  p = DimPlot(sc, reduction = "umap.atac.unintegrated", group.by = batch_column,
              split.by = "label_exists")
  ggsave(paste0("umap_atac.unintegrated_", batch_column, "_label.split.pdf"), p,
         width = 2 * OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
  ggsave(paste0("umap_atac.unintegrated_", batch_column, "_label.split.png"), p,
         width = 2 * OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)

  sc <- RunHarmony(
    object = sc,
    group.by.vars = batch_column,
    reduction.use = "atac.lsi",
    assay.use = "ATAC",
    project.dim = FALSE
  )
  sc <- FindNeighbors(sc, reduction = "harmony",
                      dims = 2:n_dims, verbose = FALSE)
  sc <- RunUMAP(sc, reduction = "harmony",
                dims = 2:n_dims, reduction.name = "umap.atac",
                return.model = TRUE, verbose = FALSE)

  p <- DimPlot(sc, reduction = "umap.atac", group.by = batch_column,
               split.by = "label_exists")
  ggsave(paste0("umap_atac.integrated_", batch_column, ".pdf"), p,
         width = 2 * OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
  ggsave(paste0("umap_atac.integrated_", batch_column, ".png"), p,
         width = 2 * OUTPUT_FIG_WIDTH, height = OUTPUT_FIG_HEIGHT)
}
###############################################################################
#### CLUSTERING AND ANNOTATION ####
###############################################################################
# Annotate PBMC cell types using Azimuth's PBMC reference
sc_v3 = sc
sc_v3[["HTO"]] = NULL
sc_v3[["ATAC"]] = as(sc_v3[["ATAC"]], Class = "Assay")
# REPLACE AZIMUTH REFERENCE WITH APPROPRIATE DATASET
sc_v3 <- RunAzimuth(sc_v3, query.modality = "ATAC", reference = "pbmcref")

sc$predicted.celltype.l1 = sc_v3$predicted.celltype.l1
sc$predicted.celltype.l2 = sc_v3$predicted.celltype.l2
sc$predicted.celltype.l3 = sc_v3$predicted.celltype.l3
sc$predicted.celltype.l1.score = sc_v3$predicted.celltype.l1.score
sc$predicted.celltype.l2.score = sc_v3$predicted.celltype.l2.score
sc$predicted.celltype.l3.score = sc_v3$predicted.celltype.l3.score

sc@reductions$ref.umap = sc_v3@reductions$ref.umap

graph = "ATAC_snn"
for (res in c(1, 0.5, 0.25, 0.1, 0.05)) {
  sc <- FindClusters(sc, resolution = res, graph.name = graph,
                     algorithm = 4, verbose = FALSE)
}

p = clustree(sc, prefix = paste0(graph, "_res."))
ggsave("clustree_clusters_atac.png", p,
       width=OUTPUT_FIG_WIDTH, height=OUTPUT_FIG_HEIGHT)

sc$seurat_clusters = sc[[paste0(graph, "_res.", 0.25)]]
Idents(sc) = "seurat_clusters"
sc$seurat_clusters = factor(sc$seurat_clusters)

# DAR testing between clusters (one vs all)
all_markers = FindAllMarkers(sc, verbose = FALSE, assay = "ATAC")
all_markers = all_markers[all_markers$p_val_adj < 0.05, ]
closest_feats = ClosestFeature(sc, regions=rownames(all_markers))
all_markers$gene = closest_feats$gene_name[match(rownames(all_markers),
                                                 closest_feats$query_region)]
write.csv(all_markers, paste0("DAR_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)
###############################################################################
# save Seurat object
saveRDS(sc, paste0(PROJECT_DIR,"/data/qc_atac.hto_", PROJECT_NAME, ".RDS"))
# Save the R session environment information
capture.output(sessionInfo(),
               file=paste0(PROJECT_DIR, "/",
                           PROJECT_NAME,
                           ".Rsession.Info.",
                           gsub("\\D", "", Sys.time()), ".txt"))
###############################################################################
###############################################################################
# create ShinyCell app with data - MUST pre-authenticate using shinyapps.io
#      token with rsconnect
if (FALSE) {
  sc <- AddMetaData(sc, t(LayerData(sc, assay = "HTO")))
  DefaultAssay(sc) = "ATAC"
  gene_activities = GeneActivity(sc)
  sc[["pseudoRNA"]] <- CreateAssayObject(counts = gene_activities)
  sc <- NormalizeData(
    object = sc,
    assay = "pseudoRNA",
    normalization.method = "LogNormalize",
    scale.factor = median(sc$nCount_pseudoRNA)
  )
  DefaultAssay(sc) = "pseudoRNA"
  sc <- FindVariableFeatures(sc)
  sc_conf = createConfig(sc)
  makeShinyApp(sc, sc_conf, gene.mapping = TRUE,
               shiny.title = paste0(PROJECT_NAME, " ATAC (pseudoRNA) + HTO"),
               shiny.dir = paste0("shiny_", PROJECT_NAME, "_atac"),
               gex.assay="pseudoRNA")
  rsconnect::deployApp(paste0("shiny_", PROJECT_NAME, "_atac"))
}