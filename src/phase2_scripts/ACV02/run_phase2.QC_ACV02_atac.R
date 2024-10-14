# run_phase2.QC_ACV02_atac.R
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
PROJECT_NAME = "ACV02_all"
# REPLACE, path to ATAC.ASAP analysis dir from MAPseq pipeline
PROJECT_DIR = paste0("/home/Projects/Scharer_sc/ACV02",
                     "/ACV02_all/analysis/ATAC.ASAP")
RAW_SEURAT_PATH = paste0(PROJECT_DIR,
                         "/data/raw_atac.hto_", PROJECT_NAME, ".RDS")

GENOME = "GRCh38"                     # REPLACE (GRCh38 or GRCm39)
OUTPUT_FIG_WIDTH =  8                 # inches, width of output figures
OUTPUT_FIG_HEIGHT = 8                 # inches, height of output figures
###############################################################################
#### LOAD DATA ####
###############################################################################
setwd(PROJECT_DIR)
sc_total = readRDS(RAW_SEURAT_PATH)
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
       width = OUTPUT_FIG_WIDTH * floor(ncol * 0.5))

p = VlnPlot(sc_total,
            features = rownames(sc_total[["HTO"]]),
            ncol = ncol,
            group.by = "HTO_classification",
            pt.size = 0)
ggsave(paste0("vln_classification_", PROJECT_NAME, ".png"),
       p, height = OUTPUT_FIG_HEIGHT,
       width = OUTPUT_FIG_WIDTH * floor(ncol * 0.5))
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
# Label doublets from object
if ("genotype_status" %in% colnames(sc_total)) {
  sc_total$doublet_status = (sc_total$genotype_status == "doublet") &
    (sc_total$HTO_classification.global == "Doublet")
} else {
  sc_total$doublet_status = (sc_total$HTO_classification.global == "Doublet")
}
sc_total$negative_status = is.na(sc_total$patient_id)

# compute nucleosome signal score per cell
sc_total = NucleosomeSignal(object = sc_total, verbose = FALSE)
# compute TSS enrichment score per cell
sc_total = TSSEnrichment(object = sc_total, verbose = FALSE)
# add blacklist ratio and fraction of reads in peaks
sc_total$pct_reads_in_peaks = sc_total$peak_region_fragments /
  sc_total$passed_filters * 100
sc_total$blacklist_ratio = sc_total$blacklist_region_fragments /
  sc_total$peak_region_fragments

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
# SAVE RAW SEURAT OBJECT
###############################################################################
saveRDS(sc_total,
        paste0(PROJECT_DIR, "/data/raw_atac.hto_", PROJECT_NAME, ".RDS"))
###############################################################################
#### ATAC QC CUTOFFS ####
###############################################################################
# PAUSE, view scatter figures above and determine appropriate cutoffs below
MIN_PEAK_FRAGMENTS = 500   # REPLACE, minimum peak fragments per cell
MIN_PCT_RiP = 65            # REPLACE, minimum percent reads in peaks per cell
MAX_BLACKLIST_RATIO = 1.0   # REPLACE, maximum blacklist ratio per cell
MAX_NUCLEOSOME_SIG = 1      # REPLACE, maximum nucleosome signal per cell
MIN_TSS = 4                 # REPLACE, minimum TSS enrichment score per cell

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
sc_total$patient_id[sc_total$doublet_status] = "Doublet"
sc_total$patient_id[sc_total$negative_status] = "Negative"
sc_total$patient_id[is.na(sc_total$patient_id)] = "Negative"
# create new column for unique sample ID - adjust as needed for each dataset
sc_total$sample_id = paste(sc_total$atac_id,
                           sc_total$patient_id,
                           sc_total$treatment,
                           sc_total$visit,
                           sep = "-")
stats = data.frame(sample_id = unique(sc_total$sample_id))
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
               doublet_status == FALSE)

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

# Create Filtered Feature Scatter of ATAC and HTO fragments/reads
cells_use = colnames(sc[, !is.na(sc$patient_id)])
p = FeatureScatter(sc, "nCount_ATAC", "nCount_HTO", group.by = "run_id",
                   split.by = "HTO_maxID", cells = cells_use, ncol = ncol) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = run_id),
               sc@meta.data[!is.na(sc$patient_id), ]) +
  labs(title = paste(PROJECT_NAME, "ATAC vs HTO Read Depths"))
ggsave("featscatter_nctATAC.v.nctHTO_HTO_maxid_filtereddata.png", p,
       width = OUTPUT_FIG_WIDTH * 2, height = OUTPUT_FIG_HEIGHT)

p = FeatureScatter(sc, "nCount_ATAC", "nCount_HTO", group.by = "run_id",
                   split.by = "hash.ID", cells = cells_use, ncol = ncol) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  stat_ellipse(aes(group = run_id),
               sc@meta.data[!is.na(sc$patient_id), ]) +
  labs(title = paste(PROJECT_NAME, "ATAC vs HTO Read Depths"))
ggsave("featscatter_nctATAC.v.nctHTO_calledHT_filtereddata.png", p,
       width = OUTPUT_FIG_WIDTH * 2, height = OUTPUT_FIG_HEIGHT)
###############################################################################
#### BATCH CORRECTION (OPTIONAL) ####
###############################################################################
sc_na = sc[, colnames(sc)[is.na(sc$endotype)]]
sc <- sc[, colnames(sc)[!is.na(sc$endotype)]]
sc_na$label_exists = FALSE
sc$label_exists = TRUE

if (FALSE) {
  DefaultAssay(sc) = "ATAC"
  batch_column = "run_id"
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
n_dims = 30

sc <- RunTFIDF(sc, min.cells = 1)
sc <- FindTopFeatures(sc, min.cutoff = "q0")
sc <- RunSVD(sc, n = n_dims, reduction.name = "atac.lsi",
             verbose = FALSE)

## LOOK AT THIS PLOT AND SET VARIABLE ##
p = ElbowPlot(sc, ndims = n_dims, reduction = "atac.lsi")
ggsave("elbow_plot.png", p, width = OUTPUT_FIG_WIDTH,
       height = OUTPUT_FIG_HEIGHT)
n_dims_keep = 10
########################################

sc <- FindNeighbors(sc, dims = 2:n_dims_keep, reduction = "atac.lsi",
                    verbose = FALSE)
sc <- RunUMAP(sc, dims = 2:n_dims_keep, reduction = "atac.lsi",
              reduction.name = "umap.atac",
              return.model = TRUE, verbose = FALSE)
###############################################################################
#### CLUSTERING AND ANNOTATION ####
###############################################################################
# Annotate PBMC cell types using Azimuth's PBMC reference
# REPLACE AZIMUTH REFERENCE WITH APPROPRIATE DATASET

if (!file.exists("annotation_reference/ext.Rds")) {
  dir.create(file.path("annotation_reference"))
  download.file("https://zenodo.org/records/7770374/files/ext.Rds",
                destfile = "annotation_reference/ext.Rds") # PBMC reference
}
DefaultAssay(sc) = "ATAC"
new_frags = CreateFragmentObject(Fragments(sc)[[1]]@path)
Fragments(sc) = NULL
Fragments(sc) = new_frags
sc <- RunAzimuth(sc, query.modality = "ATAC",
                 reference = file.path("annotation_reference"))

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
closest_feats = ClosestFeature(sc, regions = rownames(all_markers))
all_markers$gene = closest_feats$gene_name[match(rownames(all_markers),
                                                 closest_feats$query_region)]
write.csv(all_markers, paste0("DAR_", graph, ".clusters.res0.25.csv"),
          row.names = FALSE, quote = FALSE)
###############################################################################
# save Seurat object
saveRDS(sc, paste0(PROJECT_DIR, "/data/qc_atac.hto_", PROJECT_NAME, ".RDS"))
# Save the R session environment information
capture.output(sessionInfo(),
               file = paste0(PROJECT_DIR, "/",
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
###############################################################################
# Save h5ad for CellxGene use (https://github.com/chanzuckerberg/cellxgene)
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
  sceasy::convertFormat(sc, from = "seurat", to = "anndata",
                        outFile = paste0("qc_atac.pseudorna_", PROJECT_NAME, ".h5ad"))
}